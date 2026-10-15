"""Strict cell-property table datasets: schema, validation, qualification and closure construction.

The contract is documented in docs/cell_table_contract.md. Three questions are kept separate:

1. Schema validity: the files are well formed, hashed, unit-labelled, sign/area/reference
   conventions are the ones this model uses, and the values are internally consistent.
2. Numerical qualification: the dataset is a physical cell calculation with mesh, gradient and
   energy-balance evidence meeting the predeclared targets and a demonstrated physical endpoint.
3. Case match: geometry, fluid, vapor, properties and endpoint agree with the case being run.

A manifest flag such as ``validated: true`` is ignored; status is computed from evidence.
Manufactured fixtures and analytical examples are never admitted to literature cases.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from ..provenance import load_sources, validate_provenance
from .model import CaseRole, DataRole, VERIFICATION_ONLY_ROLES
from .transport import ClosureLimits, TabulatedClosure

SCHEMA_VERSION = "nanofab.cell_tables/1"
MANIFEST_NAME = "manifest.yaml"

PRESSURE_CONVENTION = "Pc = p_v - p_liquid [Pa], positive suction"
HYDRAULIC_NORMALIZATION = "bulk: Q = (W h K_bulk / mu) dPc/dx, A_ref = W h"
THERMAL_REFERENCE_PLANE = "pillar_base_plane"
AREA_BASIS = "footprint"
VAPOR_COMPOSITION = "pure_saturated_vapor"

COLUMN_UNITS = {
    "Pc_Pa": "Pa", "K_bulk_m2": "m2", "tau_m3": "m3", "K_pore_m2": "m2",
    "Q_cell_m3_s": "m3/s", "Delta_p_cell_Pa": "Pa", "ell_m": "m", "b_m": "m", "h_m": "m", "mu_Pa_s": "Pa*s",
    "h_e_W_m2K": "W/m2/K", "Q_in_W": "W", "Q_evap_W": "W", "A_footprint_m2": "m2", "T_base_K": "K",
    "T_v_K": "K", "q_in_W_m2": "W/m2",
    "state_id": "-", "branch_id": "-", "min_clearance_m": "m", "interface_area_m2": "m2",
    "contact_angle_measure_deg": "deg", "admissible": "-", "A_liquid_section_m2": "m2", "section_location": "-",
}
STRING_COLUMNS = {"state_id", "branch_id", "section_location"}
REQUIRED_COLUMNS = {
    "hydraulic": {"Pc_Pa"},
    "thermal": {"Pc_Pa", "h_e_W_m2K", "Q_in_W", "Q_evap_W", "A_footprint_m2", "T_base_K", "T_v_K", "q_in_W_m2"},
    "meniscus": {"Pc_Pa", "state_id", "branch_id", "min_clearance_m", "interface_area_m2",
                 "contact_angle_measure_deg", "admissible"},
}
COMSOL_HYDRAULIC_COLUMNS = {"Q_cell_m3_s", "Delta_p_cell_Pa", "ell_m", "b_m", "h_m", "mu_Pa_s"}

REQUIRED_KEYS = [
    "schema_version", "dataset_id", "data_role", "provenance", "source_ids", "source_locations", "geometry",
    "fluid", "solid", "vapor", "generation", "conventions", "cell", "boundary_conditions", "kinetic_model",
    "pressure_interval_Pa", "physical_endpoint", "geometry_limit", "mesh_evidence", "gradient_study",
    "supported_ranges", "energy_balance_tolerance", "tables", "limitations",
]

# Predeclared numerical qualification targets (docs/cell_table_contract.md).
QUALIFICATION_TARGETS = {
    "mesh_levels_min": 3,
    "mesh_finest_pair_rel_change_max": {"hydraulic": 1e-3, "thermal": 1e-3, "meniscus": 1e-3},
    "gradient_K_variation_max": 0.01,
    "thermal_energy_balance_max": 1e-3,
    "consistency_rel": 1e-6,
}
CASE_MATCH_REL = {"geometry": 1e-9, "properties": 1e-3, "temperature_K": 0.01, "conductivity": 1e-3,
                  "endpoint": 1e-3}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_hash(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


@dataclass
class ValidationReport:
    dataset_id: str | None
    data_role: str | None
    path: str
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    qualification_failures: list[str] = field(default_factory=list)
    case_mismatches: list[str] = field(default_factory=list)
    numerically_qualified: bool = False
    case_matched: bool | None = None
    dataset_hash: str | None = None

    @property
    def schema_valid(self) -> bool:
        return not self.errors

    @property
    def verification_only(self) -> bool:
        return self.data_role in {r.value for r in VERIFICATION_ONLY_ROLES}

    def admitted_for(self, case_role: CaseRole) -> tuple[bool, list[str]]:
        reasons = []
        if not self.schema_valid:
            reasons.append("schema invalid")
        if case_role == CaseRole.LITERATURE_REPLICATION:
            if self.data_role != DataRole.PHYSICAL_CELL_CALCULATION.value:
                reasons.append(f"data_role {self.data_role} is never admitted to a literature case")
            if not self.numerically_qualified:
                reasons.append("not numerically qualified")
        if self.case_matched is False:
            reasons.append("does not match the case")
        if self.case_matched is None:
            reasons.append("case match not evaluated")
        return (not reasons), reasons

    def as_dict(self) -> dict:
        return {
            "dataset_id": self.dataset_id, "data_role": self.data_role, "path": self.path,
            "dataset_hash": self.dataset_hash, "schema_valid": self.schema_valid,
            "verification_only": self.verification_only, "numerically_qualified": self.numerically_qualified,
            "case_matched": self.case_matched, "errors": self.errors, "warnings": self.warnings,
            "qualification_failures": self.qualification_failures, "case_mismatches": self.case_mismatches,
        }


@dataclass
class CellTableDataset:
    path: Path
    manifest: dict
    frames: dict[str, pd.DataFrame]
    report: ValidationReport

    @property
    def dataset_id(self) -> str:
        return self.manifest["dataset_id"]

    def to_closure(self, snap_rel: float = 1e-9, knot_stride: int = 1) -> TabulatedClosure:
        """Piecewise-linear closure. ``knot_stride > 1`` drops interior knots for sampling studies."""
        if not self.report.schema_valid:
            raise ValueError(f"dataset {self.path} is schema-invalid: {self.report.errors}")
        m = self.manifest
        hyd = self.frames["hydraulic"]
        Pc_h = hyd["Pc_Pa"].to_numpy(float)
        tau = bulk_transmissivity(hyd, m)
        th = self.frames["thermal"]
        lo, hi = m["pressure_interval_Pa"]
        snap = snap_rel * max(hi - lo, 1.0)
        men = self.frames["meniscus"]
        adm = men[men["admissible"].astype(int) == 1]
        Pc_area = area = None
        if "A_liquid_section_m2" in men.columns:
            b = float(m["cell"]["represented_width_b_m"])
            Pc_area = adm["Pc_Pa"].to_numpy(float)
            area = adm["A_liquid_section_m2"].to_numpy(float) / b
        pe, gl, sr = m["physical_endpoint"], m["geometry_limit"], m["supported_ranges"]
        limits = ClosureLimits(
            Pc_physical_limit=float(pe["Pc_Pa"]) if pe.get("demonstrated") else None,
            physical_limit_criterion=str(pe.get("criterion", "")),
            physical_limit_demonstrated=bool(pe.get("demonstrated")),
            Pc_geometry_limit=None if gl.get("Pc_Pa") is None else float(gl["Pc_Pa"]),
            geometry_limit_criterion=str(gl.get("criterion", "")),
            geometry_limit_demonstrated=bool(gl.get("demonstrated")),
            supported_gradient_Pa_m=tuple(float(v) for v in sr["pressure_gradient_Pa_m"]),
            supported_heat_flux_W_m2=tuple(float(v) for v in sr["heat_flux_W_m2"]),
        )
        # meniscus admissibility bounds the usable pressure interval
        Pc_men_hi = float(adm["Pc_Pa"].max())
        keep_h = Pc_h <= Pc_men_hi + snap
        keep_t = th["Pc_Pa"].to_numpy(float) <= Pc_men_hi + snap
        Pc_t, he = th["Pc_Pa"].to_numpy(float)[keep_t], th["h_e_W_m2K"].to_numpy(float)[keep_t]
        Pc_h, tau = Pc_h[keep_h], tau[keep_h]
        if knot_stride > 1:
            Pc_h, tau = _stride(Pc_h, tau, knot_stride)
            Pc_t, he = _stride(Pc_t, he, knot_stride)
            if area is not None:
                Pc_area, area = _stride(Pc_area, area, knot_stride)
        return TabulatedClosure(Pc_h, tau, Pc_t, he, limits, snap=snap,
                                Pc_area=Pc_area, area_per_width=area,
                                name=self.dataset_id + (f"[stride {knot_stride}]" if knot_stride > 1 else ""))


def _stride(x: np.ndarray, y: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Every k-th knot, always keeping both endpoints (table-sampling sensitivity only)."""
    idx = np.unique(np.r_[np.arange(0, x.size, k), x.size - 1])
    return x[idx], y[idx]


def bulk_transmissivity(hyd: pd.DataFrame, manifest: dict) -> np.ndarray:
    """tau = h K_bulk [m^3], converting a pore-normalized K only through an explicit stated area ratio."""
    h = float(manifest["geometry"]["h_m"])
    spec = manifest["tables"]["hydraulic"]
    if "tau_m3" in hyd.columns:
        return hyd["tau_m3"].to_numpy(float)
    if "K_bulk_m2" in hyd.columns:
        return h * hyd["K_bulk_m2"].to_numpy(float)
    if "K_pore_m2" in hyd.columns:
        ratio = float(spec["pore_to_bulk_area_ratio"])
        return h * ratio * hyd["K_pore_m2"].to_numpy(float)
    raise ValueError("hydraulic table has no K_bulk_m2, tau_m3 or K_pore_m2 column")


def _num(x) -> bool:
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _pos(x) -> bool:
    return _num(x) and float(x) > 0


def _rel(a, b) -> float:
    a, b = float(a), float(b)
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def load_dataset(path: str | Path, sources: dict | None = None, repo_root: Path | None = None) -> CellTableDataset:
    """Load and schema-validate a dataset directory. Never raises for content problems; see report."""
    path = Path(path)
    report = ValidationReport(None, None, str(path))
    frames: dict[str, pd.DataFrame] = {}
    mpath = path / MANIFEST_NAME
    if not mpath.is_file():
        report.errors.append(f"missing {MANIFEST_NAME}")
        return CellTableDataset(path, {}, frames, report)
    try:
        m = yaml.safe_load(mpath.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        report.errors.append(f"manifest YAML error: {exc}")
        return CellTableDataset(path, {}, frames, report)
    if not isinstance(m, dict):
        report.errors.append("manifest must be a mapping")
        return CellTableDataset(path, {}, frames, report)
    report.dataset_id = m.get("dataset_id")
    report.data_role = m.get("data_role")
    E = report.errors

    for k in REQUIRED_KEYS:
        if k not in m:
            E.append(f"missing manifest key '{k}'")
    if E:
        return CellTableDataset(path, m, frames, report)
    if "validated" in m:
        report.warnings.append("manifest 'validated' flag ignored; qualification is computed from evidence")
    if m["schema_version"] != SCHEMA_VERSION:
        E.append(f"schema_version {m['schema_version']!r} != {SCHEMA_VERSION!r}")
    try:
        role = DataRole(m["data_role"])
    except ValueError:
        E.append(f"unknown data_role {m['data_role']!r}")
        role = None
    try:
        validate_provenance(str(m["provenance"]))
    except ValueError as exc:
        E.append(str(exc))
    sources = sources if sources is not None else load_sources()
    sids = m["source_ids"] or []
    if not isinstance(sids, list):
        E.append("source_ids must be a list")
        sids = []
    for sid in sids:
        if sid not in sources:
            E.append(f"source id {sid!r} not registered in config/sources.yaml")
    if role in (DataRole.PHYSICAL_CELL_CALCULATION, DataRole.LITERATURE_DIGITIZATION) and not sids:
        E.append("physical/literature data need at least one registered source id")
    if role in VERIFICATION_ONLY_ROLES and not m.get("generating_expressions"):
        E.append("fixtures and analytical examples must record their generating_expressions")

    conv = m["conventions"] or {}
    for key, expected in (("pressure", PRESSURE_CONVENTION), ("hydraulic_normalization", HYDRAULIC_NORMALIZATION),
                          ("thermal_reference_plane", THERMAL_REFERENCE_PLANE), ("area_basis", AREA_BASIS)):
        if conv.get(key) != expected:
            E.append(f"conventions.{key} = {conv.get(key)!r}; expected {expected!r} "
                     "(no automatic reconciliation of sign, area or reference-plane conventions)")
    if not conv.get("curvature_normal"):
        E.append("conventions.curvature_normal must state the Eq. (2) normal convention")

    geo = m["geometry"] or {}
    for k in ("d_m", "l_m", "h_m"):
        if not _pos(geo.get(k)):
            E.append(f"geometry.{k} must be positive and finite")
    if _pos(geo.get("d_m")) and _pos(geo.get("l_m")) and float(geo["d_m"]) >= float(geo["l_m"]):
        E.append("geometry: pillars overlap (d >= l)")
    for k in ("geometry_id", "surface_id", "array"):
        if not geo.get(k):
            E.append(f"geometry.{k} required")
    if not _num(geo.get("receding_contact_angle_deg")):
        E.append("geometry.receding_contact_angle_deg required")
    fl = m["fluid"] or {}
    for k in ("rho_kg_m3", "mu_Pa_s", "sigma_N_m", "h_fg_J_kg", "k_W_mK", "cp_J_kgK", "T_ref_K"):
        if not _pos(fl.get(k)):
            E.append(f"fluid.{k} must be positive and finite")
    for k in ("name", "property_source", "property_version"):
        if not fl.get(k):
            E.append(f"fluid.{k} required")
    so = m["solid"] or {}
    if not _pos(so.get("k_W_mK")) or not so.get("name") or not so.get("k_status"):
        E.append("solid.name, solid.k_W_mK (>0) and solid.k_status are required")
    va = m["vapor"] or {}
    if va.get("composition") != VAPOR_COMPOSITION:
        E.append(f"vapor.composition must be {VAPOR_COMPOSITION!r}")
    for k in ("T_v_K", "p_v_Pa"):
        if not _pos(va.get(k)):
            E.append(f"vapor.{k} must be positive")
    gen = m["generation"] or {}
    for k in ("method", "software", "software_version"):
        if not gen.get(k):
            E.append(f"generation.{k} required")
    cell = m["cell"] or {}
    for k in ("flow_length_ell_m", "represented_width_b_m", "footprint_area_m2", "thermal_cell_substrate_m"):
        if not _pos(cell.get(k)):
            E.append(f"cell.{k} must be positive")
    if not _pos(cell.get("symmetry_factor")):
        E.append("cell.symmetry_factor must be positive")
    bc = m["boundary_conditions"] or {}
    for k in ("hydraulic", "thermal", "meniscus"):
        if not bc.get(k):
            E.append(f"boundary_conditions.{k} required")
    km = m["kinetic_model"] or {}
    for k in ("name", "accommodation_treatment"):
        if not km.get(k):
            E.append(f"kinetic_model.{k} required")
    if not isinstance(km.get("accommodation_treatment_verified"), bool):
        E.append("kinetic_model.accommodation_treatment_verified must be true/false")
    pi = m["pressure_interval_Pa"]
    if not (isinstance(pi, list) and len(pi) == 2 and _num(pi[0]) and _num(pi[1]) and float(pi[1]) > float(pi[0])):
        E.append("pressure_interval_Pa must be [lo, hi] with hi > lo")
        pi = None
    pe = m["physical_endpoint"] or {}
    if not isinstance(pe.get("demonstrated"), bool) or not pe.get("criterion"):
        E.append("physical_endpoint needs criterion and demonstrated (true/false)")
    if pe.get("demonstrated") and not (_num(pe.get("Pc_Pa")) and pe.get("evidence")):
        E.append("a demonstrated physical_endpoint needs Pc_Pa and evidence")
    gl = m["geometry_limit"] or {}
    if not isinstance(gl.get("demonstrated"), bool):
        E.append("geometry_limit.demonstrated must be true/false")
    sr = m["supported_ranges"] or {}
    for k in ("pressure_gradient_Pa_m", "heat_flux_W_m2", "T_v_K"):
        v = sr.get(k)
        if not (isinstance(v, list) and len(v) == 2 and all(_num(x) for x in v) and float(v[1]) >= float(v[0])):
            E.append(f"supported_ranges.{k} must be [lo, hi]")
    if not _pos(m["energy_balance_tolerance"]):
        E.append("energy_balance_tolerance must be positive")
    if not isinstance(m["limitations"], list):
        E.append("limitations must be a list (may be empty only with justification)")

    tables = m["tables"] or {}
    for name in ("hydraulic", "thermal", "meniscus"):
        spec = tables.get(name)
        if not spec:
            E.append(f"tables.{name} missing")
            continue
        df = _load_table(path, name, spec, E)
        if df is not None:
            frames[name] = df
    if E:
        return CellTableDataset(path, m, frames, report)

    _check_values(m, frames, pi, report)
    if not E:
        report.dataset_hash = canonical_hash(m)
        _qualify(m, path, report, repo_root)
    return CellTableDataset(path, m, frames, report)


def _load_table(path: Path, name: str, spec: dict, E: list) -> pd.DataFrame | None:
    f = spec.get("file")
    if not f:
        E.append(f"tables.{name}.file required")
        return None
    fp = path / f
    if not fp.is_file():
        E.append(f"tables.{name}: file {f} not found")
        return None
    digest = sha256_file(fp)
    if spec.get("sha256") != digest:
        E.append(f"tables.{name}: sha256 mismatch (manifest {spec.get('sha256')}, file {digest})")
    cols = spec.get("columns") or {}
    if not isinstance(cols, dict) or not cols:
        E.append(f"tables.{name}.columns must map column -> unit")
        return None
    for c, u in cols.items():
        if c not in COLUMN_UNITS:
            E.append(f"tables.{name}: unknown column {c!r}")
        elif u != COLUMN_UNITS[c]:
            E.append(f"tables.{name}: column {c!r} declared unit {u!r}, required {COLUMN_UNITS[c]!r}")
    missing = REQUIRED_COLUMNS[name] - set(cols)
    if missing:
        E.append(f"tables.{name}: missing required columns {sorted(missing)}")
    if name == "hydraulic" and not ({"K_bulk_m2", "tau_m3", "K_pore_m2"} & set(cols)):
        E.append("tables.hydraulic needs K_bulk_m2, tau_m3 or (with a stated area ratio) K_pore_m2")
    try:
        df = pd.read_csv(fp, comment="#")
    except Exception as exc:  # noqa: BLE001 - report any parse failure
        E.append(f"tables.{name}: CSV parse error {exc}")
        return None
    extra = [c for c in df.columns if c not in cols]
    absent = [c for c in cols if c not in df.columns]
    shared_ok = spec.get("shared_file", False)
    if extra and not shared_ok:
        E.append(f"tables.{name}: undeclared columns {extra}")
    if absent:
        E.append(f"tables.{name}: declared columns missing from file {absent}")
        return None
    df = df[list(cols)]
    for c in cols:
        if c in STRING_COLUMNS:
            df[c] = df[c].astype(str)
            continue
        vals = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        if not np.all(np.isfinite(vals)):
            E.append(f"tables.{name}: nonfinite or nonnumeric values in {c}")
        df[c] = vals
    return df


def _check_values(m: dict, frames: dict, pi, report: ValidationReport) -> None:
    E, W = report.errors, report.warnings
    geo, fl, va, cell = m["geometry"], m["fluid"], m["vapor"], m["cell"]
    tol = QUALIFICATION_TARGETS["consistency_rel"]
    for name, df in frames.items():
        Pc = df["Pc_Pa"].to_numpy(float)
        if np.any(Pc < 0):
            E.append(f"{name}: negative Pc rows; Pc = p_v - p_liquid must be >= 0 for this passive case "
                     "(check the pressure sign convention)")
        if np.any(np.diff(Pc) <= 0):
            E.append(f"{name}: Pc must be strictly increasing (duplicate or unsorted states)")
        if pi is not None and (Pc.min() > float(pi[0]) + 1e-9 or Pc.max() < float(pi[1]) - 1e-9):
            E.append(f"{name}: rows span [{Pc.min()}, {Pc.max()}] Pa but pressure_interval_Pa is {pi}")
    hyd = frames["hydraulic"]
    h = float(geo["h_m"])
    for c in ("K_bulk_m2", "tau_m3", "K_pore_m2"):
        if c in hyd.columns and np.any(hyd[c].to_numpy(float) <= 0):
            E.append(f"hydraulic: {c} must be positive")
    if "K_pore_m2" in hyd.columns and not _pos(m["tables"]["hydraulic"].get("pore_to_bulk_area_ratio")):
        E.append("hydraulic: K_pore_m2 requires tables.hydraulic.pore_to_bulk_area_ratio (stated area)")
    if "tau_m3" in hyd.columns and "K_bulk_m2" in hyd.columns:
        if np.max(np.abs(hyd["tau_m3"] - h * hyd["K_bulk_m2"]) / hyd["tau_m3"]) > tol:
            E.append("hydraulic: tau_m3 != h K_bulk_m2")
    present = COMSOL_HYDRAULIC_COLUMNS & set(hyd.columns)
    if present and present != COMSOL_HYDRAULIC_COLUMNS:
        E.append(f"hydraulic: incomplete COMSOL source columns, missing {sorted(COMSOL_HYDRAULIC_COLUMNS - present)}")
    elif present:
        K_rec = hyd["mu_Pa_s"] * hyd["Q_cell_m3_s"] * hyd["ell_m"] / (hyd["b_m"] * hyd["h_m"] * hyd["Delta_p_cell_Pa"])
        if np.any(hyd["Delta_p_cell_Pa"] <= 0) or np.any(hyd["Q_cell_m3_s"] <= 0):
            E.append("hydraulic: Delta_p_cell_Pa and Q_cell_m3_s must be positive (flow in +x for +Delta p)")
        if "K_bulk_m2" in hyd.columns and np.max(np.abs(K_rec - hyd["K_bulk_m2"]) / hyd["K_bulk_m2"]) > tol:
            E.append("hydraulic: K_bulk_m2 != mu Q ell / (b h Delta_p) (area or length normalization mismatch)")
        if np.max(np.abs(hyd["h_m"] - h)) > 1e-9 * h:
            E.append("hydraulic: h_m differs from geometry.h_m")
        if np.max(np.abs(hyd["ell_m"] - float(cell["flow_length_ell_m"]))) > 1e-12 or \
                np.max(np.abs(hyd["b_m"] - float(cell["represented_width_b_m"]))) > 1e-12:
            E.append("hydraulic: ell_m/b_m differ from the cell definition")
        if np.max(np.abs(hyd["mu_Pa_s"] - float(fl["mu_Pa_s"])) / float(fl["mu_Pa_s"])) > CASE_MATCH_REL["properties"]:
            E.append("hydraulic: mu_Pa_s differs from fluid.mu_Pa_s")
    th = frames["thermal"]
    if np.any(th["h_e_W_m2K"] <= 0):
        E.append("thermal: h_e must be positive")
    if np.any(th["Q_in_W"] <= 0):
        E.append("thermal: zero or negative cell load; h_e at zero load is 0/0 and needs a documented "
                 "small-load limit, not a 0/0 row")
    if np.any(th["T_base_K"] <= th["T_v_K"]):
        E.append("thermal: T_base_K must exceed T_v_K (evaporating branch)")
    if np.max(np.abs(th["T_v_K"] - float(va["T_v_K"]))) > CASE_MATCH_REL["temperature_K"]:
        E.append("thermal: T_v_K column differs from vapor.T_v_K")
    if np.max(np.abs(th["A_footprint_m2"] - float(cell["footprint_area_m2"])) / float(cell["footprint_area_m2"])) > tol:
        E.append("thermal: A_footprint_m2 differs from cell.footprint_area_m2 (area basis mismatch)")
    if not E:
        he_rec = th["Q_in_W"] / (th["A_footprint_m2"] * (th["T_base_K"] - th["T_v_K"]))
        if np.max(np.abs(he_rec - th["h_e_W_m2K"]) / th["h_e_W_m2K"]) > tol:
            E.append("thermal: h_e != Q_in / (A_footprint (T_base - T_v))")
        q_rec = th["Q_in_W"] / th["A_footprint_m2"]
        if np.max(np.abs(q_rec - th["q_in_W_m2"]) / th["q_in_W_m2"]) > tol:
            E.append("thermal: q_in_W_m2 != Q_in / A_footprint")
        eb = np.abs(th["Q_in_W"] - th["Q_evap_W"]) / th["Q_in_W"]
        if np.max(eb) > float(m["energy_balance_tolerance"]):
            E.append(f"thermal: energy balance {np.max(eb):.3g} exceeds declared tolerance")
    men = frames["meniscus"]
    if men["state_id"].duplicated().any():
        E.append("meniscus: duplicate state_id")
    if not set(np.unique(men["admissible"])) <= {0.0, 1.0}:
        E.append("meniscus: admissible must be 0 or 1")
    adm = men[men["admissible"] == 1]
    if len(adm) < 2:
        E.append("meniscus: fewer than two admissible states")
    elif np.any(adm["min_clearance_m"] <= 0):
        E.append("meniscus: admissible states must have positive liquid clearance")
    if np.any(men["interface_area_m2"] <= 0):
        E.append("meniscus: interface_area_m2 must be positive")
    if "A_liquid_section_m2" in men.columns and "section_location" not in men.columns:
        E.append("meniscus: A_liquid_section_m2 needs section_location")
    if len(adm) >= 1 and adm["Pc_Pa"].iloc[0] > float(m["pressure_interval_Pa"][0]) + 1e-9:
        E.append("meniscus: first admissible state above the start of pressure_interval_Pa")
    if "A_liquid_section_m2" not in men.columns:
        W.append("no A_liquid_section_m2: intrinsic section velocity diagnostics are unavailable")


def _qualify(m: dict, path: Path, report: ValidationReport, repo_root: Path | None) -> None:
    Q = report.qualification_failures
    role = m["data_role"]
    if role != DataRole.PHYSICAL_CELL_CALCULATION.value:
        Q.append(f"data_role {role} is not a physical cell calculation")
    if m["provenance"] not in ("COMSOL_EXTRACTED", "MODEL_CALCULATED"):
        Q.append(f"provenance {m['provenance']} is not a computed cell result")
    gen = m["generation"]
    if not gen.get("model_file") or not gen.get("export_files"):
        Q.append("generation.model_file and export_files are required for qualification")
    for ef in gen.get("export_files") or []:
        p = Path(ef.get("path", ""))
        cands = [p] if p.is_absolute() else [path / p] + ([repo_root / p] if repo_root else [])
        hit = next((c for c in cands if c.is_file()), None)
        if hit is None:
            Q.append(f"raw export {p} not available for hash verification")
        elif sha256_file(hit) != ef.get("sha256"):
            Q.append(f"raw export {p} hash mismatch")
    if not gen.get("run_log"):
        Q.append("generation.run_log (actual solver log) required")
    me = m["mesh_evidence"] or {}
    for name in ("hydraulic", "thermal", "meniscus"):
        ev = me.get(name) or {}
        levels = ev.get("levels") or []
        if len(levels) < QUALIFICATION_TARGETS["mesh_levels_min"]:
            Q.append(f"mesh_evidence.{name}: fewer than {QUALIFICATION_TARGETS['mesh_levels_min']} mesh levels")
        rc = ev.get("finest_pair_rel_change")
        lim = QUALIFICATION_TARGETS["mesh_finest_pair_rel_change_max"][name]
        if not _num(rc) or float(rc) > lim:
            Q.append(f"mesh_evidence.{name}: finest-pair change {rc} not <= {lim}")
        if ev.get("systematic") is not True:
            Q.append(f"mesh_evidence.{name}: systematic convergence not demonstrated")
    gs = m["gradient_study"] or {}
    var = gs.get("max_rel_K_variation")
    if not _num(var) or float(var) > QUALIFICATION_TARGETS["gradient_K_variation_max"]:
        Q.append(f"gradient_study: K variation {var} not <= {QUALIFICATION_TARGETS['gradient_K_variation_max']} "
                 "(a one-dimensional K(Pc) table is not qualified; Q(Pc, gradient) data would be needed)")
    if len(gs.get("Pc_values_Pa") or []) < 3:
        Q.append("gradient_study: needs representative Pc values spanning the interval")
    if gs.get("stokes_vs_navier_stokes_rel_diff") is None:
        Q.append("gradient_study: Stokes vs Navier-Stokes comparison missing")
    if float(m["energy_balance_tolerance"]) > QUALIFICATION_TARGETS["thermal_energy_balance_max"]:
        Q.append("energy_balance_tolerance looser than the qualification target")
    if not m["physical_endpoint"].get("demonstrated"):
        Q.append("physical endpoint not demonstrated (final CSV row is not evidence of a physical endpoint)")
    report.numerically_qualified = not Q


def match_case(ds: CellTableDataset, case_req: dict) -> None:
    """Compare the dataset with a case's requirements; fills report.case_matched/case_mismatches."""
    r = ds.report
    if not r.schema_valid:
        r.case_matched = False
        r.case_mismatches.append("schema invalid")
        return
    m, X = ds.manifest, r.case_mismatches
    g = m["geometry"]
    for key, val in (("d_m", case_req["d"]), ("l_m", case_req["l"]), ("h_m", case_req["h"])):
        if _rel(g[key], val) > CASE_MATCH_REL["geometry"]:
            X.append(f"geometry {key} {g[key]} != case {val} (mixed geometries are not interpolated)")
    if abs(float(g["receding_contact_angle_deg"]) - case_req["theta_rec_deg"]) > 1e-9:
        X.append("receding contact angle differs from the case")
    fl, va = m["fluid"], m["vapor"]
    if str(fl["name"]).lower() != str(case_req["fluid"]).lower():
        X.append(f"fluid {fl['name']} != case {case_req['fluid']}")
    if abs(float(va["T_v_K"]) - case_req["T_v"]) > CASE_MATCH_REL["temperature_K"]:
        X.append(f"vapor T_v {va['T_v_K']} K != case {case_req['T_v']} K")
    if abs(float(fl["T_ref_K"]) - case_req["T_v"]) > CASE_MATCH_REL["temperature_K"]:
        X.append("fluid property reference temperature differs from the case T_v")
    for key, val in (("rho_kg_m3", case_req["rho"]), ("mu_Pa_s", case_req["mu"]), ("sigma_N_m", case_req["sigma"]),
                     ("h_fg_J_kg", case_req["h_fg"])):
        if _rel(fl[key], val) > CASE_MATCH_REL["properties"]:
            X.append(f"fluid {key} {fl[key]} differs from case {val} by more than "
                     f"{CASE_MATCH_REL['properties']:.0e} (property-source mismatch)")
    if "k_s" in case_req and _rel(m["solid"]["k_W_mK"], case_req["k_s"]) > CASE_MATCH_REL["conductivity"]:
        X.append(f"cell solid conductivity {m['solid']['k_W_mK']} != case substrate {case_req['k_s']}")
    if "cell_substrate" in case_req and _rel(m["cell"]["thermal_cell_substrate_m"], case_req["cell_substrate"]) > 1e-9:
        X.append("thermal cell substrate thickness differs from the case")
    if case_req.get("require_verified_accommodation") and not m["kinetic_model"].get("accommodation_treatment_verified"):
        X.append("accommodation treatment is not source-verified (V2019 replication requires it)")
    if case_req.get("kinetic_model") and m["kinetic_model"]["name"] != case_req["kinetic_model"]:
        X.append(f"kinetic model {m['kinetic_model']['name']} != required {case_req['kinetic_model']}")
    pe = m["physical_endpoint"]
    Pc_lim = case_req.get("Pc_limit")
    if Pc_lim is not None and pe.get("demonstrated") and _rel(pe["Pc_Pa"], Pc_lim) > CASE_MATCH_REL["endpoint"]:
        X.append(f"table physical endpoint {pe['Pc_Pa']} Pa inconsistent with case limit {Pc_lim:.6g} Pa")
    r.case_matched = not X


def write_dataset(out_dir: Path, manifest: dict, frames: dict[str, pd.DataFrame]) -> Path:
    """Write tables and a manifest with hashes. Refuses to overwrite an existing dataset."""
    out_dir = Path(out_dir)
    if (out_dir / MANIFEST_NAME).exists():
        raise FileExistsError(f"{out_dir} already holds a dataset; choose a new DATASET_ID")
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(json.dumps(manifest, default=float))
    for name, df in frames.items():
        fname = f"{name}.csv"
        df.to_csv(out_dir / fname, index=False, float_format="%.12g")
        manifest["tables"][name]["file"] = fname
        manifest["tables"][name]["sha256"] = sha256_file(out_dir / fname)
        manifest["tables"][name]["columns"] = {c: COLUMN_UNITS[c] for c in df.columns}
    (out_dir / MANIFEST_NAME).write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return out_dir / MANIFEST_NAME
