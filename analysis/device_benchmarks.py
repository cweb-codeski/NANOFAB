"""Numerical verification benchmarks for the passive device model (predeclared gates).

    python analysis/device_benchmarks.py --suite all --output data/processed/device_runs/verification
    python analysis/device_benchmarks.py --suite thermal --output data/processed/device_runs/verification

Suites: analytic, accounting, zero_load, hydraulics, thermal, coupled, statuses, reference_scales.
All inputs are ANALYTICAL_EXAMPLE or MANUFACTURED_FIXTURE (VERIFICATION ONLY). Passing gates means
"numerically verified for these cases"; it is not a literature comparison or a validation.
Writes <output>/<run_id>/{run_manifest.json, gates.csv, <suite>.json} and plots to
plots/device/<run_id>/. Exit 0 if every gate passes, 1 otherwise.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

import _common  # noqa: F401  (puts src/ on sys.path)

from nanofab.device import verification as V
from nanofab.device.output import new_run_dir, rel_path, run_header, write_json


def plot_refinement(results: dict, plot_dir: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    files = []
    if "thermal" in results and "hydraulics" in results:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        th = pd.DataFrame(results["thermal"]["details"]["levels"])
        ax[0].loglog(th.Nx, th.err_T_e_K, "o-", label="T_e (top face)")
        ax[0].loglog(th.Nx, th.err_T_source_K, "s-", label="T_source (bottom face)")
        ax[0].loglog(th.Nx, th.err_T_e_K.iloc[0] * (th.Nx.iloc[0] / th.Nx) ** 2, "k--", lw=0.8, label="slope 2")
        ax[0].set(xlabel="Nx (Nz = Nx/6.25)", ylabel="max error [K]", title="manufactured x-z conduction")
        ax[0].legend()
        hy = pd.DataFrame(results["hydraulics"]["details"]["levels"])
        ax[1].loglog(hy.Nx, hy.err_continuous_q_Pa, "o-", label="Pc vs ODE reference, continuous q(x)")
        ax[1].loglog(hy.Nx, hy.err_continuous_q_Pa.iloc[0] * (hy.Nx.iloc[0] / hy.Nx) ** 2, "k--", lw=0.8, label="slope 2")
        ax[1].set(xlabel="Nx", ylabel="max |Pc error| [Pa]", title="manufactured-table hydraulics")
        ax[1].legend()
        fig.suptitle("VERIFICATION ONLY (manufactured inputs)")
        fig.tight_layout()
        p = plot_dir / "refinement_orders.png"
        fig.savefig(p, dpi=120)
        plt.close(fig)
        files.append(rel_path(p))
    return files


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suite", default="all", choices=["all", *V.SUITES])
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--plot-root", type=Path, default=None, help="plot parent directory (default plots/device)")
    args = ap.parse_args(argv)
    argv_list = sys.argv[:] if argv is None else ["analysis/device_benchmarks.py", *argv]
    run_id, run_dir, plot_dir = new_run_dir(args.output, f"bench_{args.suite}", args.plot_root)
    names = list(V.SUITES) if args.suite == "all" else [args.suite]
    wd = run_dir / "fixture_variants"
    results, rows = {}, []
    for name in names:
        t0 = time.time()
        try:
            r = V.SUITES[name](wd)
        except Exception as exc:  # noqa: BLE001 - a crash is a failed suite, recorded as such
            r = {"suite": name, "passed": False, "gates": [V.gate("suite raised", repr(exc), "no exception", "equal")],
                 "details": {}}
        r["seconds"] = time.time() - t0
        results[name] = r
        write_json(run_dir / f"{name}.json", r)
        for g in r["gates"]:
            rows.append({"suite": name, **{k: (str(v) if isinstance(v, (tuple, list)) else v) for k, v in g.items()}})
        print(f"[{'PASS' if r['passed'] else 'FAIL'}] {name} ({r['seconds']:.1f} s, {len(r['gates'])} gates)")
        for g in r["gates"]:
            if not g["passed"]:
                print(f"    FAILED: {g['gate']}: {g['value']} (limit {g['limit']})")
    pd.DataFrame(rows).to_csv(run_dir / "gates.csv", index=False)
    plots = plot_refinement(results, plot_dir)
    passed = all(r["passed"] for r in results.values())
    write_json(run_dir / "run_manifest.json", {
        "run_id": run_id, **run_header(argv_list, f"benchmarks:{args.suite}", True),
        "gates_predeclared": V.GATES, "refinement_levels": V.REFINEMENT,
        "suites": {k: {"passed": v["passed"], "n_gates": len(v["gates"]), "seconds": v["seconds"]} for k, v in results.items()},
        "all_passed": passed, "plots": plots,
        "claim": "numerically verified for these cases (ANALYTICAL_EXAMPLE / MANUFACTURED_FIXTURE inputs only)",
    })
    print(f"{'ALL GATES PASSED' if passed else 'GATE FAILURES'}; {len(rows)} gates; run directory: "
          f"{rel_path(run_dir)}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
