"""Run the passive device-scale wick model for a named case.

    python analysis/run_passive_device.py --case parameters/cases/vaartstra2019_fig4.yaml --preflight \
        --output data/processed/device_runs/vaartstra_preflight
    python analysis/run_passive_device.py --case parameters/cases/vaartstra2019_fig4.yaml --search-dryout \
        --output data/processed/device_runs/vaartstra_prediction
    python analysis/run_passive_device.py --case tests/fixtures/device/manufactured_case.yaml --search-dryout \
        --output data/processed/device_runs/verification_fixture_run

A literature case never falls back to fixtures or to the legacy analytical model: without a qualified
physical cell-table dataset it stops with MISSING_VALIDATED_INPUTS (exit 3) and writes a report
naming the missing items. Exit codes: 0 admissible / limit located, 1 invalid case, 2 usage,
3 MISSING_VALIDATED_INPUTS, 4 OUTSIDE_MODEL_REGIME, 5 TABLE_COVERAGE_LIMIT, 6 NUMERICAL_FAILURE,
7 prescribed forward load at or above the located limit.
Each invocation writes a new run directory <output>/<run_id>/ (never overwriting) and plots to
plots/device/<run_id>/.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

from _common import REPO

from nanofab.device.case import CaseError, build_problem, load_case, preflight
from nanofab.device.diagnostics import energy_balance, mass_balance, regime_diagnostics
from nanofab.device.model import EXIT_CODES, EXIT_CONFIG_ERROR, EXIT_LOAD_ABOVE_LIMIT, CaseRole, SolverSettings, SolverStatus
from nanofab.device.output import new_run_dir, run_header, write_json, write_solution
from nanofab.device.solver import search_dryout, solve_prescribed_load
from nanofab.device.tables import sha256_file


def _rel(p: Path) -> str:
    p = Path(p).resolve()
    return p.relative_to(REPO).as_posix() if p.is_relative_to(REPO) else str(p)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, type=Path, help="case manifest (parameters/cases/*.yaml)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true", help="report readiness only; never solves")
    mode.add_argument("--search-dryout", action="store_true", help="locate the limiting load")
    mode.add_argument("--load-W-cm2", type=float, help="forward solve at a prescribed mean input flux")
    ap.add_argument("--output", required=True, type=Path, help="parent directory for the run directory")
    ap.add_argument("--Nx", type=int, help="override the case's axial cell count")
    ap.add_argument("--Nz", type=int, help="override the case's substrate cell count")
    ap.add_argument("--rtol", type=float, help="override the coupled-iteration tolerance")
    ap.add_argument("--plot-root", type=Path, default=None, help="plot parent directory (default plots/device)")
    return ap.parse_args(argv)


def settings_for(case, args) -> SolverSettings:
    num = case.raw.get("numerics") or {}
    return SolverSettings(Nx=args.Nx or int(num.get("Nx", 200)), Nz=args.Nz or int(num.get("Nz", 32)),
                          rtol=args.rtol or float(num.get("rtol", 1e-10)))


def comparison_block(case, q_W_m2, status) -> list[dict]:
    """Literature comparison, kept separate from verification and never fed back into inputs."""
    out = []
    for t in case.raw.get("comparison_targets") or []:
        if t.get("parameter") is None:
            out.append({"quantity": t["quantity"], "status": "NO DATA", "note": t.get("note")})
            continue
        ref = case.ps.meta(t["parameter"])
        rec = {"quantity": t["quantity"], "reference_SI": ref.value, "reference_unit": ref.si_unit,
               "reference_provenance": ref.provenance, "reference_source": ref.source,
               "use": "COMPARISON_ONLY", "threshold_rel": t.get("investigation_threshold_rel")}
        if q_W_m2 is None or status not in (SolverStatus.PHYSICAL_CAPILLARY_LIMIT, SolverStatus.MENISCUS_GEOMETRY_LIMIT):
            rec["status"] = "NOT COMPARED (no located physical limit)"
        else:
            rel = (q_W_m2 - ref.value) / ref.value
            rec.update({"model_SI": q_W_m2, "relative_difference": rel,
                        "status": "INVESTIGATE" if abs(rel) > t.get("investigation_threshold_rel", 0.05) else "WITHIN THRESHOLD"})
        out.append(rec)
    return out


def main(argv=None) -> int:
    args = parse_args(argv)
    argv_list = sys.argv[:] if argv is None else ["analysis/run_passive_device.py", *argv]
    try:
        case = load_case(args.case)
    except CaseError as exc:
        print(f"invalid case: {exc}", file=sys.stderr)
        return EXIT_CONFIG_ERROR
    verification_only = case.role == CaseRole.VERIFICATION_FIXTURE
    mode = "preflight" if args.preflight else ("search_dryout" if args.search_dryout else "forward_load")
    run_id, run_dir, plot_dir = new_run_dir(args.output, mode, args.plot_root)
    header = run_header(argv_list, f"{mode} / case_role={case.role.value}", verification_only)
    pf = preflight(case)
    dataset = pf.pop("_dataset")
    manifest = {
        "run_id": run_id, **header,
        "case": {"path": _rel(case.path), "sha256": sha256_file(case.path), "case_id": case.case_id,
                 "case_role": case.role.value, "raw": case.raw},
        "parameters": case.parameter_snapshot(),
        "liquid_properties": case.liquid.as_dict(),
        "boundary_conditions": {
            "x0": f"reservoir edge: Pc = {case.raw['topology']['feed_suction_Pa']} Pa, insulated substrate end",
            "xL": "symmetry: Q = 0, insulated substrate end",
            "substrate_bottom": "prescribed q_in''(x)", "substrate_top": "Robin q_e'' = h_e(Pc)(T_e - T_v) at the pillar-base plane",
            "liquid_inlet_temperature": "T_v (V2019)",
        },
        "preflight": pf,
        "plots_dir": _rel(plot_dir),
    }
    print(f"run {run_id}: case {case.case_id} ({case.role.value}); preflight {pf['status']}")
    if verification_only:
        print("VERIFICATION ONLY: manufactured inputs; results are not physical predictions")
    if not pf["ready_for_prediction"]:
        status = SolverStatus.MISSING_VALIDATED_INPUTS
        manifest["status"] = status.value
        manifest["exit_code"] = EXIT_CODES[status]
        write_json(run_dir / "missing_inputs.json", {k: pf[k] for k in ("missing_inputs", "blocking_source_gaps",
                                                                        "provisional_inputs", "data_inventory")})
        write_json(run_dir / "run_manifest.json", manifest)
        print(f"{status.value}: nothing was solved. Missing:")
        for m in pf["missing_inputs"]:
            print(f"  - {m['input']}: {m['reason']}" + ("  [non-blocking: " + m["blocks"] + "]" if "blocks" in m else ""))
        for g in pf["blocking_source_gaps"]:
            print(f"  - source gap: {g['gap']}")
        print(f"report: {_rel(run_dir / 'missing_inputs.json')}")
        return EXIT_CODES[status]
    if args.preflight:
        manifest["status"] = pf["status"]
        manifest["exit_code"] = 0
        write_json(run_dir / "run_manifest.json", manifest)
        print(f"preflight passed; report: {_rel(run_dir / 'run_manifest.json')}")
        return 0

    settings = settings_for(case, args)
    problem = build_problem(case, dataset, settings.endpoint_snap_rel)
    manifest["settings"] = settings.as_dict()
    manifest["dataset"] = {"report": dataset.report.as_dict(),
                           "file_hashes": {k: v.get("sha256") for k, v in dataset.manifest["tables"].items()},
                           "closure_limits": asdict(problem.closure.limits),
                           "closure_range_Pa": [problem.closure.Pc_lo, problem.closure.Pc_hi]}
    manifest["capillary_limit"] = {"Pc_Pa": problem.Pc_physical_limit, "criterion": problem.physical_limit_criterion}
    if args.search_dryout:
        res = search_dryout(problem, settings, regime_check=lambda s: regime_diagnostics(problem, s))
        status = res.status
        manifest["result"] = res.summary()
        sols = {"endpoint": res.endpoint}
        q_report = res.upper_bracket if res.upper_bracket is not None else res.lower_operating_point
        manifest["literature_comparison"] = comparison_block(case, q_report, status) if not verification_only else \
            "not applicable (VERIFICATION ONLY)"
    else:
        sol = solve_prescribed_load(problem, args.load_W_cm2 * 1e4, settings)
        status = sol.status
        manifest["result"] = {"status": status.value, "reason": sol.reason, "load_W_m2": sol.load,
                              "iterations": sol.iterations}
        sols = {"forward": sol}
        if sol.converged:
            manifest["result"]["regime"] = regime_diagnostics(problem, sol)
    for tag, sol in sols.items():
        if sol is None or not sol.converged:
            continue
        manifest.setdefault("balances", {})[tag] = {"mass": mass_balance(problem, sol),
                                                    "energy": energy_balance(problem, sol)}
        manifest.setdefault("files", {})[tag] = write_solution(run_dir, plot_dir, problem, sol, tag, verification_only)
        manifest.setdefault("convergence", {})[tag] = {"iterations": sol.iterations,
                                                       "final": sol.history[-1] if sol.history else None}
    exit_code = EXIT_CODES[status]
    if not args.search_dryout and not sols["forward"].converged and exit_code == 0:
        exit_code = EXIT_LOAD_ABOVE_LIMIT
    manifest["status"] = status.value
    manifest["exit_code"] = exit_code
    write_json(run_dir / "run_manifest.json", manifest)
    r = manifest["result"]
    print(f"status: {status.value}")
    if not args.search_dryout:
        print(f"reason: {r['reason']}")
    if args.search_dryout:
        if r.get("upper_bracket_W_m2"):
            print(f"limit bracket: [{r['lower_operating_point_W_m2']/1e4:.6g}, {r['upper_bracket_W_m2']/1e4:.6g}] W/cm^2"
                  f" (relative width {r['relative_interval']:.2e})")
        elif r.get("lower_operating_point_W_m2"):
            print(f"certified lower operating point: {r['lower_operating_point_W_m2']/1e4:.6g} W/cm^2")
        print(f"reason: {r['reason']}")
    print(f"run directory: {_rel(run_dir)}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
