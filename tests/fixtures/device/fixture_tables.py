"""Deterministic MANUFACTURED_FIXTURE cell tables for verification. VERIFICATION ONLY.

Nothing here is a physical cell calculation. The functions below are chosen to be smooth, positive
and to vary over the pressure interval by amounts comparable to those V2019 describes
qualitatively (h_e rising with suction, K falling), so that the coupled solver is exercised. They
are not fits to any paper, and must never be used as inputs to a literature case.

Run ``python tests/fixtures/device/fixture_tables.py`` to regenerate the on-disk fixtures.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

from nanofab.device.tables import (AREA_BASIS, HYDRAULIC_NORMALIZATION, PRESSURE_CONVENTION,  # noqa: E402
                                   SCHEMA_VERSION, THERMAL_REFERENCE_PLANE, VAPOR_COMPOSITION, write_dataset)

P_STAR = 2200.0  # Pa, manufactured endpoint
K0 = 3.0e-11  # m^2
HE0 = 6.0e4  # W/(m^2 K)
CELL_LOAD = 50e4  # W/m^2, manufactured cell operating flux
D, LP, H = 10e-6, 30e-6, 25e-6
T_V = 373.15

GENERATING_EXPRESSIONS = {
    "K_bulk_m2": "3.0e-11 * (1 - 0.2 * (Pc/2200)^2)",
    "h_e_W_m2K": "6.0e4 * (1 + 0.3 * Pc/2200)",
    "min_clearance_m": "h * (1 - Pc/(2*2200))  [geometry-limit variant: h * (1 - Pc/2000)]",
    "interface_area_m2": "(l/2)^2 - pi d^2/16 times (1 + 0.1 (Pc/2200)^2)",
    "contact_angle_measure_deg": "90 - 75 * Pc/2200",
    "A_liquid_section_m2": "0.8 * b * h * (1 - 0.15 * Pc/2200)",
    "Q_cell_m3_s": "K_bulk b h Delta_p / (mu ell), Delta_p = 1 Pa",
    "T_base_K": "T_v + q_in / h_e, q_in = 5e5 W/m^2; Q_evap = Q_in",
}


def K_fn(Pc):
    return K0 * (1 - 0.2 * (np.asarray(Pc, float) / P_STAR) ** 2)


def he_fn(Pc):
    return HE0 * (1 + 0.3 * np.asarray(Pc, float) / P_STAR)


def water_props(T=T_V) -> dict:
    from nanofab.fluids import Fluid
    import CoolProp
    w = Fluid("Water")
    s = w.saturated(T)
    return {"name": "Water", "property_source": "IAPWS via CoolProp", "property_version": CoolProp.__version__,
            "T_ref_K": T, "rho_kg_m3": s.rho_l, "mu_Pa_s": s.mu_l, "sigma_N_m": s.surface_tension,
            "h_fg_J_kg": s.h_lv, "k_W_mK": s.k_l, "cp_J_kgK": s.cp_l, "p_sat_Pa": s.p_sat}


def build(pc_max: float = P_STAR, step: float = 100.0, geometry_limit: float | None = None,
          gradient_hi: float = 1e8, heat_flux_range=(0.0, 1e7), props: dict | None = None,
          endpoint_demonstrated: bool | None = None, with_area: bool = True):
    props = props or water_props()
    Pc = np.arange(0.0, pc_max + 0.5 * step, step)
    ell = b = LP / 2
    A_fp = b * ell
    mu = props["mu_Pa_s"]
    K = K_fn(Pc)
    hyd = pd.DataFrame({"Pc_Pa": Pc, "K_bulk_m2": K, "Q_cell_m3_s": K * b * H * 1.0 / (mu * ell),
                        "Delta_p_cell_Pa": 1.0, "ell_m": ell, "b_m": b, "h_m": H, "mu_Pa_s": mu})
    he = he_fn(Pc)
    Q_in = CELL_LOAD * A_fp
    th = pd.DataFrame({"Pc_Pa": Pc, "h_e_W_m2K": Q_in / (A_fp * (CELL_LOAD / he)), "Q_in_W": Q_in, "Q_evap_W": Q_in,
                       "A_footprint_m2": A_fp, "T_base_K": T_V + CELL_LOAD / he, "T_v_K": T_V,
                       "q_in_W_m2": CELL_LOAD})
    if geometry_limit is None:
        clear = H * (1 - Pc / (2 * P_STAR))
        adm = np.ones_like(Pc)
    else:
        clear = H * (1 - Pc / 2000.0)
        adm = (Pc <= geometry_limit + 1e-9).astype(float)
    men = pd.DataFrame({"Pc_Pa": Pc, "state_id": [f"mf{i:03d}" for i in range(Pc.size)],
                        "branch_id": "pinned_top_edge", "min_clearance_m": clear,
                        "interface_area_m2": (A_fp - math.pi * D**2 / 16) * (1 + 0.1 * (Pc / P_STAR) ** 2),
                        "contact_angle_measure_deg": 90 - 75 * Pc / P_STAR, "admissible": adm})
    if with_area:
        men["A_liquid_section_m2"] = 0.8 * b * H * (1 - 0.15 * Pc / P_STAR)
        men["section_location"] = "x = l/2 mid-gap plane (manufactured)"
    demonstrated = (pc_max >= P_STAR) if endpoint_demonstrated is None else endpoint_demonstrated
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": "fixture_manufactured",
        "data_role": "MANUFACTURED_FIXTURE",
        "provenance": "MODEL_CALCULATED",
        "verification_only": "VERIFICATION ONLY - manufactured functions, not a physical calculation",
        "generating_expressions": GENERATING_EXPRESSIONS,
        "source_ids": [],
        "source_locations": [],
        "geometry": {"geometry_id": "d10_l30_h25_square", "array": "square", "d_m": D, "l_m": LP, "h_m": H,
                     "surface_id": "fixture_surface", "receding_contact_angle_deg": 15.0},
        "fluid": props,
        "solid": {"name": "silicon", "k_W_mK": 148.0, "k_status": "fixture value (matches the provisional case k_s)"},
        "vapor": {"composition": VAPOR_COMPOSITION, "T_v_K": T_V, "p_v_Pa": props["p_sat_Pa"]},
        "generation": {"method": "manufactured analytic functions", "software": "tests/fixtures/device/fixture_tables.py",
                       "software_version": "1", "model_file": None, "export_files": [], "run_log": None},
        "conventions": {"pressure": PRESSURE_CONVENTION, "hydraulic_normalization": HYDRAULIC_NORMALIZATION,
                        "thermal_reference_plane": THERMAL_REFERENCE_PLANE, "area_basis": AREA_BASIS,
                        "curvature_normal": "V2019 Eq. (2): H > 0 for a convex liquid surface; Pc = -2 sigma H"},
        "cell": {"flow_length_ell_m": ell, "represented_width_b_m": b, "symmetry_factor": 4.0,
                 "footprint_area_m2": A_fp, "thermal_cell_substrate_m": 1e-6},
        "boundary_conditions": {"hydraulic": "manufactured", "thermal": "manufactured", "meniscus": "manufactured"},
        "kinetic_model": {"name": "manufactured", "accommodation_coefficient": None,
                          "accommodation_treatment": "none (manufactured h_e)",
                          "accommodation_treatment_verified": False},
        "pressure_interval_Pa": [0.0, float(Pc[-1])],
        "physical_endpoint": {"Pc_Pa": P_STAR, "criterion": "MANUFACTURED declared endpoint",
                              "demonstrated": bool(demonstrated),
                              "evidence": "fixture definition" if demonstrated else None},
        "geometry_limit": ({"Pc_Pa": geometry_limit, "criterion": "MANUFACTURED clearance threshold",
                            "demonstrated": True} if geometry_limit is not None
                           else {"Pc_Pa": None, "criterion": "none reached", "demonstrated": False}),
        "mesh_evidence": {"hydraulic": {}, "thermal": {}, "meniscus": {}},
        "gradient_study": {"Pc_values_Pa": [], "gradient_range_Pa_m": None, "max_rel_K_variation": None,
                           "stokes_vs_navier_stokes_rel_diff": None},
        "supported_ranges": {"pressure_gradient_Pa_m": [0.0, gradient_hi], "heat_flux_W_m2": list(heat_flux_range),
                             "T_v_K": [T_V, T_V]},
        "energy_balance_tolerance": 1e-9,
        "limitations": ["VERIFICATION ONLY: manufactured functions"],
        "tables": {"hydraulic": {}, "thermal": {}, "meniscus": {}},
    }
    return manifest, {"hydraulic": hyd, "thermal": th, "meniscus": men}


def write_fixture(out_dir, **kw) -> Path:
    manifest, frames = build(**kw)
    return write_dataset(Path(out_dir), manifest, frames).parent


def write_raw_export_fixture(out_dir: Path) -> Path:
    """COMSOL-format raw exports plus a raw-export manifest (MANUFACTURED_FIXTURE)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest, frames = build(step=200.0)
    hdr = ["% Model,fixture_manufactured (no .mph: MANUFACTURED FIXTURE)",
           "% Version,MANUFACTURED FIXTURE - not a COMSOL run",
           "% Date,deterministic"]
    hyd, th, men = frames["hydraulic"], frames["thermal"], frames["meniscus"]
    (out_dir / "flow_raw.csv").write_text("\n".join(hdr + ['% Pc (Pa),"Q, flow rate (m^3/s)",dp (Pa)'] + [
        f"{p:.12g},{q:.12g},{d:.12g}" for p, q, d in zip(hyd.Pc_Pa, hyd.Q_cell_m3_s, hyd.Delta_p_cell_Pa)]) + "\n")
    (out_dir / "thermal_raw.csv").write_text("\n".join(hdr + ["% Pc (Pa),Q_in (W),Q_evap (W),T_base (K)"] + [
        f"{p:.12g},{a:.12g},{b:.12g},{t:.12g}" for p, a, b, t in zip(th.Pc_Pa, th.Q_in_W, th.Q_evap_W, th.T_base_K)])
        + "\n")
    (out_dir / "meniscus_raw.csv").write_text("\n".join(hdr + [
        "% Pc (Pa),zclear (m),A_int (m^2),theta_meas (deg)"] + [
        f"{p:.12g},{c:.12g},{a:.12g},{t:.12g}" for p, c, a, t in
        zip(men.Pc_Pa, men.min_clearance_m, men.interface_area_m2, men.contact_angle_measure_deg)]) + "\n")
    meta = {k: v for k, v in manifest.items() if k not in ("schema_version", "dataset_id", "tables")}
    raw = {
        "schema_version": "nanofab.raw_export_manifest/1",
        "dataset_id": "fixture_raw_import",
        "metadata": meta,
        "exports": {
            "hydraulic": {"file": "flow_raw.csv",
                          "columns": {"Pc_Pa": "Pc", "Q_cell_m3_s": "Q, flow rate", "Delta_p_cell_Pa": "dp"},
                          "constants": {"ell_m": LP / 2, "b_m": LP / 2, "h_m": H, "mu_Pa_s": meta["fluid"]["mu_Pa_s"]}},
            "thermal": {"file": "thermal_raw.csv",
                        "columns": {"Pc_Pa": "Pc", "Q_in_W": "Q_in", "Q_evap_W": "Q_evap", "T_base_K": "T_base"},
                        "constants": {"A_footprint_m2": (LP / 2) ** 2, "T_v_K": T_V}},
            "meniscus": {"file": "meniscus_raw.csv", "state_id_prefix": "raw",
                         "columns": {"Pc_Pa": "Pc", "min_clearance_m": "zclear", "interface_area_m2": "A_int",
                                     "contact_angle_measure_deg": "theta_meas"},
                         "constants": {"branch_id": "pinned_top_edge", "admissible": 1.0}},
        },
    }
    p = out_dir / "export_manifest.yaml"
    p.write_text("# VERIFICATION ONLY: manufactured raw exports for the import CLI demonstration.\n"
                 + yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return p


if __name__ == "__main__":
    import shutil
    target = HERE / "manufactured_tables"
    if target.exists():
        shutil.rmtree(target)
    print("wrote", write_fixture(target))
    print("wrote", write_raw_export_fixture(HERE / "raw_export_fixture"))
