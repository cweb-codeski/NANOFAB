"""Named device cases: loading, parameter resolution, preflight and problem construction.

A case manifest (``parameters/cases/*.yaml``) references entries in ``parameters/*.yaml`` by name, so
every scalar keeps its provenance and source location. The loader rejects design targets as inputs
and refuses comparison targets that are not marked COMPARISON_ONLY. A literature case never falls
back to fixtures or to the legacy analytical model: missing qualified tables produce
MISSING_VALIDATED_INPUTS with a report naming exactly what is missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from ..params import ParameterSet, load_parameters
from ..provenance import REPO_ROOT, load_sources
from .model import CaseRole, DeviceGeometry, LiquidProperties, SolverStatus
from .solver import CoupledProblem
from .tables import CellTableDataset, load_dataset, match_case, sha256_file

CASE_SCHEMA = "nanofab.device_case/1"
DESIGN_FILE = "design_targets.yaml"
REQUIRED_SYMBOLS = ("d", "l", "h", "L", "W", "t_s", "theta_rec", "T_v", "fluid", "k_s")
DATA_INVENTORY = REPO_ROOT / "config" / "data_inventory.yaml"


class CaseError(ValueError):
    pass


def liquid_properties(fluid_name: str, T: float) -> LiquidProperties:
    """Saturated properties at T from CoolProp (IAPWS formulations for water)."""
    import CoolProp
    from CoolProp.CoolProp import PropsSI
    from ..fluids import UNIVERSAL_GAS_CONSTANT, Fluid
    f = Fluid(fluid_name)
    s = f.saturated(T)
    cp_v = PropsSI("C", "T", T, "Q", 1, fluid_name)
    cv_v = PropsSI("O", "T", T, "Q", 1, fluid_name)
    return LiquidProperties(T_v=T, p_v=s.p_sat, rho=s.rho_l, mu=s.mu_l, h_fg=s.h_lv, sigma=s.surface_tension,
                            cp=s.cp_l, k=s.k_l, rho_v=s.rho_v, gamma_v=cp_v / cv_v,
                            R_specific=UNIVERSAL_GAS_CONSTANT / f.molar_mass,
                            source=f"IAPWS via CoolProp {CoolProp.__version__} at saturation, T = {T} K")


def eq3_capillary_limit(sigma: float, d: float, l: float, theta: float) -> float:
    """V2019 Eq. (3): Pc_limit = sigma pi d cos(theta_rec) / (l^2 - pi d^2/4) [Pa], positive suction."""
    return sigma * math.pi * d * math.cos(theta) / (l**2 - math.pi * d**2 / 4)


@dataclass
class LoadedCase:
    path: Path
    raw: dict
    role: CaseRole
    params: dict
    ps: ParameterSet
    liquid: LiquidProperties
    geometry: DeviceGeometry
    Pc_limit_eq3: float
    notes: list[str] = field(default_factory=list)

    @property
    def case_id(self) -> str:
        return self.raw["case_id"]

    def parameter_snapshot(self) -> list[dict]:
        out = []
        for sym, pname in self.raw["parameters"].items():
            p = self.ps.meta(pname)
            out.append({"symbol": sym, "parameter": pname, "value_SI": p.value, "si_unit": p.si_unit,
                        "raw_value": p.raw_value, "raw_unit": p.raw_unit, "provenance": p.provenance,
                        "source": p.source, "source_location": p.source_location, "file": p.file,
                        "notes": p.notes})
        return out

    def case_requirements(self) -> dict:
        req = dict((self.raw.get("cell_tables") or {}).get("requirements") or {})
        P = self.params
        out = {"d": P["d"], "l": P["l"], "h": P["h"], "theta_rec_deg": math.degrees(P["theta_rec"]),
               "fluid": P["fluid"], "T_v": P["T_v"], "rho": self.liquid.rho, "mu": self.liquid.mu,
               "sigma": self.liquid.sigma, "h_fg": self.liquid.h_fg, "k_s": P["k_s"],
               "require_verified_accommodation": bool(req.get("require_verified_accommodation")),
               "kinetic_model": req.get("kinetic_model")}
        if "cell_substrate" in P:
            out["cell_substrate"] = P["cell_substrate"]
        if (self.raw.get("capillary_limit") or {}).get("equation") == "V2019_Eq3":
            out["Pc_limit"] = self.Pc_limit_eq3
        return out


def load_case(path: str | Path, parameter_dir: str | Path | None = None) -> LoadedCase:
    path = Path(path)
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if raw.get("schema_version") != CASE_SCHEMA:
        raise CaseError(f"{path}: not a device case manifest (schema_version {raw.get('schema_version')!r})")
    try:
        role = CaseRole(raw["case_role"])
    except (KeyError, ValueError) as exc:
        raise CaseError(f"{path}: invalid case_role") from exc
    ps = load_parameters(parameter_dir)
    pmap = raw.get("parameters") or {}
    missing = [s for s in REQUIRED_SYMBOLS if s not in pmap]
    if missing:
        raise CaseError(f"{path}: missing parameter references {missing}")
    params = {}
    for sym, pname in pmap.items():
        if pname not in ps:
            raise CaseError(f"{path}: parameter {pname!r} (for {sym}) not found in parameters/*.yaml")
        meta = ps.meta(pname)
        if meta.file == DESIGN_FILE:
            raise CaseError(f"{path}: {pname} is a design target; design targets are never case inputs")
        params[sym] = meta.value
    for t in raw.get("comparison_targets") or []:
        if t.get("use") != "COMPARISON_ONLY":
            raise CaseError(f"{path}: comparison target {t.get('quantity')} must be marked COMPARISON_ONLY")
    if raw.get("source"):
        if raw["source"] not in load_sources():
            raise CaseError(f"{path}: source {raw['source']!r} not registered")
    liquid = liquid_properties(params["fluid"], params["T_v"])
    geom = DeviceGeometry(L=params["L"], W=params["W"], h=params["h"], t_s=params["t_s"], k_s=params["k_s"],
                          d=params["d"], l=params["l"])
    Pc_lim = eq3_capillary_limit(liquid.sigma, params["d"], params["l"], params["theta_rec"])
    return LoadedCase(path, raw, role, params, ps, liquid, geom, Pc_lim)


def _resolve(p: str | None) -> Path | None:
    if not p:
        return None
    q = Path(p)
    return q if q.is_absolute() else REPO_ROOT / q


def data_inventory() -> list[dict]:
    if not DATA_INVENTORY.is_file():
        return []
    inv = yaml.safe_load(DATA_INVENTORY.read_text(encoding="utf-8")) or {}
    out = []
    for key, e in inv.items():
        p = REPO_ROOT / e["path"]
        rec = {"id": key, **e, "present": p.exists()}
        if p.is_file():
            rec["sha256"] = sha256_file(p)
        out.append(rec)
    return out


def preflight(case: LoadedCase) -> dict:
    """Machine-readable readiness report. Never runs a solve; never substitutes inputs."""
    tables = case.raw.get("cell_tables") or {}
    ds_path = _resolve(tables.get("dataset"))
    missing, blocking_gaps, provisional, available = [], [], [], []
    for sym, pname in case.raw["parameters"].items():
        m = case.ps.meta(pname)
        available.append({"symbol": sym, "parameter": pname, "value_SI": m.value, "unit": m.si_unit,
                          "provenance": m.provenance, "source": m.source, "source_location": m.source_location})
    for sym, text in (case.raw.get("provisional_inputs") or {}).items():
        provisional.append({"symbol": sym, "status": "PROVISIONAL (non-blocking limitation)", "note": " ".join(text.split())})
    provisional.append({"symbol": "fluid_properties", "status": "AVAILABLE (runtime CoolProp)",
                        "note": "V2019 does not identify its property database; a remaining limitation."})
    dataset: CellTableDataset | None = None
    report = None
    if ds_path is None:
        for name in tables.get("required") or []:
            missing.append({"input": f"cell table: {name}", "reason": "no dataset configured in the case "
                            "(cell_tables.dataset is null)", "needed": _NEEDED[name]})
    else:
        dataset = load_dataset(ds_path, repo_root=REPO_ROOT)
        match_case(dataset, case.case_requirements())
        report = dataset.report
        ok, reasons = report.admitted_for(case.role)
        if not ok:
            missing.append({"input": f"cell table dataset {ds_path}", "reason": "; ".join(reasons),
                            "details": report.as_dict()})
    if case.role == CaseRole.LITERATURE_REPLICATION:
        req = (tables.get("requirements") or {})
        if req.get("require_verified_accommodation"):
            src = load_sources().get("meland_ytrehus2003", {})
            if not src.get("content_verified"):
                blocking_gaps.append({
                    "gap": "meland_ytrehus2003 accommodation treatment",
                    "effect": "a self-generated thermal table cannot be qualified as V2019 replication until "
                              "the nonunity-accommodation treatment is checked against Meland & Ytrehus (2003), "
                              "Phys. Fluids 15, 1348-1350, DOI 10.1063/1.1564097",
                })
        for t in case.raw.get("comparison_targets") or []:
            if t.get("file") is None and t.get("parameter") is None:
                missing.append({"input": f"comparison data: {t['quantity']}", "reason": "not available",
                                "needed": "digitized curves with documented procedure and uncertainty",
                                "blocks": "literature profile comparison only (not the solver)"})
    blocking = [m for m in missing if "blocks" not in m] + blocking_gaps
    status = SolverStatus.MISSING_VALIDATED_INPUTS if blocking else SolverStatus.ADMISSIBLE
    return {
        "case_id": case.case_id,
        "case_role": case.role.value,
        "status": status.value,
        "ready_for_prediction": not blocking,
        "missing_inputs": missing,
        "blocking_source_gaps": blocking_gaps,
        "provisional_inputs": provisional,
        "available_inputs": available,
        "capillary_limit_V2019_Eq3_Pa": case.Pc_limit_eq3,
        "liquid_properties": case.liquid.as_dict(),
        "dataset_report": report.as_dict() if report else None,
        "data_inventory": data_inventory() if case.role == CaseRole.LITERATURE_REPLICATION else [],
        "fallbacks_used": "none (fixtures and the legacy analytical model are never substituted)",
        "_dataset": dataset,
    }


_NEEDED = {
    "hydraulic": ("K_bulk(Pc) or tau(Pc) from three-dimensional laminar-flow cell calculations on the pinned "
                  "menisci, Pc = 0 to the physical endpoint, with source Q_cell, Delta_p_cell, ell, b, h, mu, "
                  ">= 3 meshes, a pressure-gradient sweep covering device gradients (<= 1% K variation) and a "
                  "Stokes vs Navier-Stokes comparison"),
    "meniscus": ("pinned-meniscus states from V2019 Eqs. (1)-(2), Pc = 0 to the physical endpoint, with "
                 "minimum clearance, interface area, contact-angle/limit measure, admissibility and mesh "
                 "convergence; A_liquid_section if intrinsic velocities are wanted"),
    "thermal": ("h_e(Pc) referenced to the pillar-base plane from conjugate cell conduction with a "
                "source-verified V2019 interfacial boundary (linearized BTE, a = 0.052 with the verified "
                "accommodation treatment), Q_in, Q_evap, T_base, T_v, footprint area, energy balance <= 1e-3, "
                ">= 3 meshes, and loads covering the device operating range (around 50 and 82 W/cm^2)"),
}


def build_problem(case: LoadedCase, dataset: CellTableDataset, snap_rel: float = 1e-9) -> CoupledProblem:
    closure = dataset.to_closure(snap_rel)
    cap = case.raw.get("capillary_limit") or {}
    if cap.get("equation") == "V2019_Eq3":
        Pc_lim, crit = case.Pc_limit_eq3, "V2019 Eq. (3)"
    elif cap.get("equation") == "table_endpoint":
        Pc_lim, crit = closure.limits.Pc_physical_limit, closure.limits.physical_limit_criterion
        if Pc_lim is None:
            Pc_lim = float(dataset.manifest["physical_endpoint"]["Pc_Pa"])
            crit += " (declared, not demonstrated within coverage)"
    else:
        raise CaseError(f"unknown capillary_limit.equation {cap.get('equation')!r}")
    shape = (case.raw.get("heat_load") or {}).get("shape", "uniform")
    if shape != "uniform":
        raise CaseError("only uniform heat-load shapes are defined in case files at present")
    return CoupledProblem(case.geometry, case.liquid, closure, Pc_feed=float(case.raw["topology"]["feed_suction_Pa"]),
                          Pc_physical_limit=Pc_lim, physical_limit_criterion=crit)
