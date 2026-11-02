"""Generate validation/VALIDATION_REPORT.md for the whole-wick model.

    python validation/run_validation.py

Verification (does the code solve its equations?) and validation (do the equations describe the
device?) are reported separately. Passing verification checks does not make the model validated.
"""

from __future__ import annotations

import csv
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nanofab import dryout as dryout1d  # noqa: E402
from nanofab.capillary import pillar_max_capillary_pressure  # noqa: E402
from nanofab.params import load_parameters  # noqa: E402
from nanofab.permeability import pillar_wick_permeability  # noqa: E402
from nanofab.whole_wick import cos_theta_V, find_dryout, load_case_set, solve_coupled  # noqa: E402
from nanofab.whole_wick.electrowetting import voltage_for_full_wetting  # noqa: E402
from nanofab.whole_wick.fv import face_conductances, net_face_inflow  # noqa: E402
from nanofab.whole_wick.report import markdown_table, write_csv  # noqa: E402

OUT = ROOT / "validation"
RECIPE = "cells.property_recipes.analytical_flat_placeholder_he"
CURVATURE_K = {
    "method": "analytical_curvature",
    "reduction_at_pc_max": {"value": 0.37, "unit": "-", "provenance": "ASSUMED", "source": "vaartstra2019",
                            "source_location": "calibrated to Sec. 3.2 82 W/cm2"},
}


def grid_study(case: str, heat_map: str, sizes: list[tuple[int, int]], extra: dict | None = None) -> list[dict]:
    rows = []
    for nx, ny in sizes:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cs = load_case_set(overrides={**(extra or {}), "wick.grid.nx.value": nx, "wick.grid.ny.value": ny})
            m = cs.build(case, heat_map)
            t0 = time.perf_counter()
            d = find_dryout(m)
            n = solve_coupled(m, 1.0)
        rows.append({
            "nx": nx, "ny": ny, "dx_um": m.grid.dx * 1e6,
            "dx_over_pitch": m.grid.dx / 30e-6,
            "dryout_W_cm2": d.mean_heat_flux / 1e4,
            "T_max_C": n.thermal.T_max - 273.15,
            "Pc_max_ratio_nominal": float(n.Pc.max() / n.props.Pc_max.max()),
            "seconds": time.perf_counter() - t0,
            "W2_warning": any("pillar pitches" in str(w.message) for w in caught),
        })
    return rows


def observed_order(rows: list[dict], key: str) -> tuple[float, float, float]:
    """Observed order p, Richardson estimate and GCI (fine grid, Fs = 1.25) from the 3 finest grids."""
    f1, f2, f3 = (rows[-1][key], rows[-2][key], rows[-3][key])
    r = rows[-1]["nx"] / rows[-2]["nx"]
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0 or e32 / e21 <= 0:
        return math.nan, f1, math.nan
    p = math.log(e32 / e21) / math.log(r)
    if p < 0.5:  # not in the asymptotic range (e.g. the governing cell moves between grids)
        return math.nan, f1, math.nan
    f_ext = f1 + (f1 - f2) / (r**p - 1)
    gci = 1.25 * abs((f1 - f2) / f1) / (r**p - 1)
    return p, f_ext, gci


def solver_sensitivity() -> list[dict]:
    rows = []
    for label, extra in (("flat K", {}), ("K(Pc), a = 0.37", {f"{RECIPE}.K": CURVATURE_K})):
        for tol in (1e-3, 1e-5, 1e-8, 1e-11):
            for relax in (0.6, 1.0):
                cs = load_case_set(overrides={**extra, "wick.solver.tol.value": tol, "wick.solver.relax.value": relax})
                m = cs.build("passive_uniform", "gaussian_hotspot")
                n = solve_coupled(m, 1.0)
                d = find_dryout(m)
                rows.append({"K model": label, "tol": tol, "relax": relax, "iterations": n.iterations,
                             "converged": n.converged, "T_max_C": n.thermal.T_max - 273.15,
                             "max_Pc_Pa": float(n.Pc.max()), "dryout_W_cm2": d.mean_heat_flux / 1e4})
    return rows


def conservation_table(cs) -> list[dict]:
    rows = []
    for hm in cs.heat_map_names:
        for case in cs.case_names:
            n = solve_coupled(cs.build(case, hm), 1.0)
            rows.append({"heat_map": hm, "case": case,
                         "mass_rel": n.darcy.mass.relative_residual,
                         "mass_cell_max_kg_s": n.darcy.mass.max_cell_residual,
                         "energy_rel": n.thermal.energy.relative_residual,
                         "Q_in_W": n.thermal.energy.heat_input, "Q_evap_W": n.thermal.energy.evaporative,
                         "inlet_kg_s": n.darcy.mass.inlet, "evap_kg_s": n.darcy.mass.evaporation})
    return rows


def h_scaling() -> list[dict]:
    rows = []
    for h in (10, 25, 50):
        geo = {}
        for key, val in (("pillar_diameter", 0.4 * h), ("pillar_pitch", 1.2 * h), ("pillar_height", h)):
            geo[f"cells.cell_types.si_d10_l30_h25.{key}"] = {
                "value": val, "unit": "um", "provenance": "PRIMARY_SOURCE", "source": "vaartstra2019",
                "source_location": "Sec. 3.3, Fig. 5 (d/h = 0.4, l/h = 1.2)"}
        cs = load_case_set(overrides=geo)
        d = find_dryout(cs.build("passive_uniform", "uniform_50"))
        rows.append({"h_um": h, "dryout_W_cm2": d.mean_heat_flux / 1e4})
    hs = np.log([r["h_um"] for r in rows])
    qs = np.log([r["dryout_W_cm2"] for r in rows])
    slope = float(np.polyfit(hs, qs, 1)[0])
    for r in rows:
        r["fitted_exponent"] = slope
    return rows


def lateral_conduction_ratio(result) -> float:
    """max |net lateral conduction| / local heat input. Analogous to vaartstra2019 A.1, not identical."""
    m = result.model
    G = m.substrate_k * m.substrate_t + result.props.k_eff * m.per_cell("h")
    Cx, Cy = face_conductances(G, m.grid, m.options.face_average)
    lat = net_face_inflow(result.thermal.T, Cx, Cy) / m.grid.cell_area
    return float(np.max(np.abs(lat) / np.maximum(result.q_in, 1e-30)))


def stage1_rows() -> list[dict]:
    with open(ROOT / "data" / "processed" / "permeability_stage1_flow.csv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    ps = load_parameters()
    cs = load_case_set()
    md = [
        "# Whole-wick model: validation report",
        "",
        "**Status: NOT VALIDATED.** The code passes its verification checks: conservation, analytical",
        "limits and grid convergence. Validation against independent device-scale data has not been",
        "done. h_e(Pc) is an ASSUMED placeholder, and K(Pc) is still the flat-meniscus analytical",
        "model. Electrowetting uses Lippmann-Young with an ASSUMED ε_r and saturation angle, and no",
        "electrowetting measurement on this wick exists. The explicit-array COMSOL comparison has not",
        "been done.",
        "",
        "Generated by `python validation/run_validation.py`. Raw tables are in `validation/data/`.",
        "",
        "| # | Item | Kind | Status |",
        "| --- | --- | --- | --- |",
        "| 1 | Grid refinement | verification | done (below) |",
        "| 2 | Nonlinear solver sensitivity (steady model, so no time step) | verification | done |",
        "| 3 | Conservation errors | verification | done |",
        "| 3b | Electrowetting implementation (V = 0 limit, dryout ∝ cos θ_V) | verification | done |",
        "| 4 | Unit-cell COMSOL outputs | input check | Stage 1 K only; Stages 3-5 pending |",
        "| 5 | vaartstra2019 modeling trends | partial validation | mixed; see the per-trend labels |",
        "| 6 | Homogenized vs explicit small-array COMSOL | validation of homogenization | **pending**, protocol below |",
        "",
    ]

    # 1. grid refinement
    K = pillar_wick_permeability(ps["pillar_diameter"], ps["pillar_pitch"], ps["pillar_height"])
    pcm = pillar_max_capillary_pressure(cs.fluid.surface_tension, ps["pillar_diameter"], ps["pillar_pitch"],
                                        ps["receding_contact_angle"])
    q1d = dryout1d.dryout_heat_flux(K, ps["pillar_height"], ps["wick_feed_length"], pcm, cs.fluid.rho_l,
                                    cs.fluid.mu_l, cs.fluid.h_lv) / 1e4
    sizes_1d = [(10, 2), (20, 2), (40, 2), (80, 2), (160, 2), (320, 2)]
    g_const = grid_study("passive_uniform", "uniform_50", sizes_1d, {f"{RECIPE}.h_e.rise_at_pc_max.value": 0.0})
    write_csv(g_const, OUT / "data" / "grid_uniform_constant_he.csv")
    p_c, ext_c, gci_c = observed_order(g_const, "dryout_W_cm2")
    g_uni = grid_study("passive_uniform", "uniform_50", sizes_1d)
    g_hot = grid_study("passive_uniform", "gaussian_hotspot", [(10, 10), (20, 20), (40, 40), (80, 80), (160, 160)])
    write_csv(g_uni, OUT / "data" / "grid_uniform.csv")
    write_csv(g_hot, OUT / "data" / "grid_hotspot.csv")
    cols = [("nx", "nx", "d"), ("ny", "ny", "d"), ("dx_um", "dx [µm]", ".1f"), ("dx_over_pitch", "dx / pitch", ".1f"),
            ("dryout_W_cm2", "dryout [W/cm²]", ".4f"), ("T_max_C", "T_max at 50 W/cm² [°C]", ".4f"),
            ("W2_warning", "below 3 pitches (W2)", ""), ("seconds", "time [s]", ".1f")]
    p_u, ext_u, gci_u = observed_order(g_uni, "dryout_W_cm2")
    p_h, ext_h, gci_h = observed_order(g_hot, "dryout_W_cm2")
    pt_h, extT_h, gciT_h = observed_order(g_hot, "T_max_C")
    def order_text(p, ext, gci, unit, ref=None):
        if math.isnan(p):
            return "not in the asymptotic range (non-monotone), so no order or GCI is reported"
        s = f"observed order {p:.2f}, Richardson extrapolation {ext:.4f} {unit}, GCI(fine) {100 * gci:.4f}%"
        if ref is not None:
            s += f", {100 * (ext / ref - 1):+.4f}% vs exact 1D"
        return s

    hot_change = 100 * max(abs(r["dryout_W_cm2"] / g_hot[-1]["dryout_W_cm2"] - 1) for r in g_hot[-3:])
    md += [
        "## 1. Grid refinement",
        "",
        "**(a) Exact-solution check.** Uniform wick, uniform 50 W/cm², h_e held constant (rise = 0),",
        "so evaporation is uniform and the 1D analytical solution of `nanofab.dryout` is exact:",
        f"**{q1d:.4f} W/cm²** (MODEL_CALCULATED, same K and Pc_max).",
        "",
        markdown_table(g_const, cols),
        "",
        f"Largest deviation from the exact value on any grid: "
        f"{100 * max(abs(r['dryout_W_cm2'] / q1d - 1) for r in g_const):.1e}%. For this quadratic pressure",
        "profile, the cell-centered scheme's half-cell boundary error and its peak-sampling error",
        "cancel, so the match is exact on every grid. The check confirms assembly, source terms and",
        "boundary conditions; the convergence *rate* is shown by (b) and (c).",
        "",
        "**(b) Baseline inputs.** The same wick with the placeholder h_e(Pc). It converges to a",
        "slightly lower value than (a). h_e rises toward the center, so lateral conduction moves",
        "evaporation inward and raises the pressure drop. That is a model effect, not a grid error.",
        "",
        markdown_table(g_uni, cols),
        "",
        "Dryout: " + order_text(p_u, ext_u, gci_u, "W/cm²") + ".",
        "",
        "**(c) Gaussian hotspot (2D).**",
        "",
        markdown_table(g_hot, cols),
        "",
        "Dryout: " + order_text(p_h, ext_h, gci_h, "W/cm²")
        + f"; the largest change over the three finest grids is {hot_change:.3f}%. The first-dry cell moves"
        " between grids, so the onset value converges non-monotonically.",
        "T_max: " + order_text(pt_h, extT_h, gciT_h, "°C") + ".",
        "",
        "Caveat (W2): below about 3 pitches per cell the homogenized properties are not",
        "representative. Refinement past that point converges the continuum equations, not the",
        "physics of the pillar array, so production runs should keep dx ≥ 3 pitches (≥ 90 µm here).",
        "The discretization error at nx = 40 is already below the model-form uncertainty from K(Pc)",
        "and h_e.",
        "",
    ]

    # 2. solver sensitivity
    sens = solver_sensitivity()
    write_csv(sens, OUT / "data" / "solver_sensitivity.csv")
    md += [
        "## 2. Nonlinear solver sensitivity",
        "",
        "The model is steady, so there is no time step. The nonlinearity is the Picard coupling",
        "through K(Pc) and h_e(Pc). Case: passive uniform wick, Gaussian hotspot map. The dryout",
        "search uses Brent's method with rtol = 1e-5 on the load scale.",
        "",
        markdown_table(sens, [("K model", "K model", ""), ("tol", "tol", ".0e"), ("relax", "relax", ".1f"),
                              ("iterations", "iterations", "d"), ("converged", "converged", ""),
                              ("T_max_C", "T_max [°C]", ".6f"), ("max_Pc_Pa", "max Pc [Pa]", ".3f"),
                              ("dryout_W_cm2", "dryout [W/cm²]", ".4f")]),
        "",
    ]

    # 3. conservation
    cons = conservation_table(cs)
    write_csv(cons, OUT / "data" / "conservation.csv")
    md += [
        "## 3. Conservation errors (nominal load, every case and heat map)",
        "",
        "The mass residual is recomputed from face fluxes after the solve, so it does not rely on",
        "the linear-solver residual.",
        "",
        markdown_table(cons, [("heat_map", "heat map", ""), ("case", "case", ""), ("mass_rel", "mass (rel.)", ".1e"),
                              ("mass_cell_max_kg_s", "worst cell [kg/s]", ".1e"), ("energy_rel", "energy (rel.)", ".1e"),
                              ("Q_in_W", "Q_in [W]", ".2f"), ("Q_evap_W", "Q_evap [W]", ".2f")]),
        "",
        f"Worst relative mass imbalance {max(r['mass_rel'] for r in cons):.1e}; worst energy imbalance "
        f"{max(r['energy_rel'] for r in cons):.1e}. Both are at machine precision.",
        "",
    ]

    # 3b. electrowetting verification
    gamma = cs.gamma.value
    base_ct = next(c for c in cs.cell_types if c.name == "si_d10_l30_h25")
    pc_eq3 = pillar_max_capillary_pressure(gamma, base_ct.d, base_ct.l, base_ct.theta_0)
    const_he = {f"{RECIPE}.h_e.rise_at_pc_max.value": 0.0, "wick.grid.ny.value": 2}
    ew_rows = []
    for V in (0.0, 5.0, 10.0, 20.0):
        cs_v = load_case_set(overrides={**const_he, "cases.voltage_maps.global.V.value": V})
        q0 = find_dryout(cs_v.build("passive_uniform", "uniform_50")).mean_heat_flux
        mv = cs_v.build("global_ew", "uniform_50")
        qv = find_dryout(mv).mean_heat_flux
        c, _ = cos_theta_V(base_ct.theta_0, V, base_ct.coating, gamma)
        expected = float(c) / math.cos(base_ct.theta_0)
        ew_rows.append({"V": V, "theta_V_deg": math.degrees(math.acos(float(c))), "expected_ratio": expected,
                        "model_ratio": qv / q0, "rel_error": qv / q0 / expected - 1})
    write_csv(ew_rows, OUT / "data" / "electrowetting_verification.csv")
    md += [
        "## 3b. Electrowetting implementation",
        "",
        f"At V = 0, Pc_max(0) = 2γ cos θ_0 / r_eff = {base_ct.pc_max0(gamma):.4f} Pa. vaartstra2019 Eq. 3 gives "
        f"{pc_eq3:.4f} Pa, so r_eff = {base_ct.r_eff * 1e6:.3f} µm reproduces Eq. 3 exactly.",
        "",
        "With constant h_e and flat K, the 1D dryout load is proportional to Pc_max. Globally actuated",
        "dryout should therefore scale by exactly cos θ_V / cos θ_0:",
        "",
        markdown_table(ew_rows, [("V", "V [V]", ".1f"), ("theta_V_deg", "θ_V [°]", ".3f"),
                                 ("expected_ratio", "cos θ_V / cos θ_0", ".6f"),
                                 ("model_ratio", "dryout(V) / dryout(0)", ".6f"), ("rel_error", "rel. error", ".1e")]),
        "",
        f"Upper bound for any voltage: cos θ_V ≤ 1, so the gain is ≤ 1/cos θ_0 − 1 = "
        f"{100 * (1 / math.cos(base_ct.theta_0) - 1):.2f}% at θ_0 = {math.degrees(base_ct.theta_0):.0f}°. "
        f"The 1 µm SiO₂ coating (ε_r ASSUMED) reaches that bound at "
        f"{voltage_for_full_wetting(base_ct.theta_0, base_ct.coating, gamma):.2f} V; higher voltage adds nothing.",
        "",
    ]

    # 4. unit-cell COMSOL
    s1 = stage1_rows()
    finest = min(s1, key=lambda r: float(r["hmax"]))
    K_c = float(finest["K_comsol_m2"])
    cs_c = load_case_set(overrides={"cells.cell_types.si_d10_l30_h25.property_recipe": "comsol_stage1_K"})
    q_c = find_dryout(cs_c.build("passive_uniform", "uniform_50")).mean_heat_flux / 1e4
    q_a = g_uni[2]["dryout_W_cm2"]
    md += [
        "## 4. Comparison with unit-cell COMSOL outputs",
        "",
        "| Quantity | COMSOL (COMSOL_EXTRACTED) | Analytical (MODEL_CALCULATED) | Difference |",
        "| --- | --- | --- | --- |",
        f"| K, d10/l30/h25, flat top, Stage 1 (hmax = {float(finest['hmax']) * 1e6:.1f} µm) | {K_c:.4e} m² | {K:.4e} m² | {100 * (K_c / K - 1):+.2f}% |",
        f"| whole-wick dryout, uniform 50 W/cm² map, nx = 40 | {q_c:.2f} W/cm² | {q_a:.2f} W/cm² | {100 * (q_c / q_a - 1):+.2f}% |",
        "",
        "Stage 1 mesh study (`data/processed/permeability_stage1_flow.csv`): K changes by "
        + ", ".join(f"{float(r['diff_pct']):+.2f}%" for r in s1)
        + " vs analytical as hmax goes from 2 to 0.5 µm. The ~1.3% gap is model-form error of the",
        "Sangani-Acrivos + Brinkman formula, which is not yet confirmed independently.",
        "",
        "**Pending.** K(Pc) (Stage 3/5), h_e(Pc) (Stage 4/5), and the Stage 3 depinning pressure to",
        "check Eq. 3's Pc_max. Until these exist, the whole-wick model's h_e and its K(Pc) are not",
        "connected to COMSOL. `unit_cell/README.md` has the hand-off steps.",
        "",
    ]

    # 5. V2019 trends
    cs_k = load_case_set(overrides={f"{RECIPE}.K": CURVATURE_K})
    d_k = find_dryout(cs_k.build("passive_uniform", "uniform_50"))
    hs = h_scaling()
    write_csv(hs, OUT / "data" / "h_scaling.csv")
    r82 = solve_coupled(cs.build("passive_uniform", "uniform_50"), 82.0 / 50.0)
    T = r82.thermal.T
    dT_inlet_center = float(T[:, 0].mean() - T[:, T.shape[1] // 2].mean())
    he = r82.h_e_used
    he_var = float((he.max() - he.min()) / he.min())
    x_first = find_dryout(cs.build("passive_uniform", "uniform_50")).location[0] * 1e3
    lat_u = lateral_conduction_ratio(solve_coupled(cs.build("passive_uniform", "uniform_50"), 1.0))
    lat_h = lateral_conduction_ratio(solve_coupled(cs.build("passive_uniform", "v2019_center_hotspot"), 1.0))
    md += [
        "## 5. Comparison with Vaartstra, Lu & Wang (2019) modeling trends",
        "",
        "Labels: **independent**: the trend is not built into the inputs, so agreement is evidence.",
        "**structural**: it follows from the analytical K and Pc_max scaling, so agreement only checks",
        "consistency. **input-driven**: the trend was put in through an ASSUMED placeholder, so",
        "agreement is circular. **calibrated**: a parameter was fitted to this target.",
        "",
        "| Trend (vaartstra2019, PRIMARY_SOURCE) | Paper | This model | Label |",
        "| --- | --- | --- | --- |",
        f"| Dryout heat flux, baseline device (Sec. 3.2) | 82 W/cm² | {q_a:.1f} W/cm² (flat K), {100 * (q_a / 82 - 1):+.1f}% | independent; disagrees as expected for flat K (A2) |",
        f"| Same, with K(Pc) = K_flat(1 - 0.37 Pc/Pc_max) | 82 W/cm² | {d_k.mean_heat_flux / 1e4:.1f} W/cm² | **calibrated** (not evidence) |",
        f"| Independent model of same device (zhou_adera2026) | 84 W/cm² | {q_a:.1f} W/cm² | independent cross-check; same flat-K caveat |",
        f"| q_dry ∝ h² at d/h = 0.4, l/h = 1.2 (Sec. 3.3, Fig. 5) | exponent 2 | exponent {hs[0]['fitted_exponent']:.3f} | structural |",
        "| h_e decreases with h (Fig. 5) | decreasing | not representable (placeholder h_e has no geometry dependence) | pending COMSOL h_e |",
        f"| h_e variation across device at 82 W/cm² (Sec. 1, 3.2) | ~30% | {100 * he_var:.0f}% | **input-driven** (W18) |",
        f"| T_e at inlet minus center at 82 W/cm² (Fig. 4a) | ~3 K | {dT_inlet_center:.2f} K | **input-driven** (follows from the h_e placeholder) |",
        f"| First dryout at the point farthest from the feed (Sec. 3.2) | device center | x = {x_first:.2f} mm of 10 | structural |",
        f"| Lateral heat flow small vs local load (A.1: 0.19% uniform, 1% hotspot) | 0.19% / 1% | {100 * lat_u:.2f}% / {100 * lat_h:.1f}% | metric only analogous (net lateral conduction / local q''_in) |",
        "",
        "The lateral-heat-flow row uses a different metric from the paper. The paper compares",
        "gradients inside one unit cell; this model reports net in-plane conduction per cell. The",
        "large hotspot value comes from the patch edges, where 100 µm of Si spreads heat over about",
        "sqrt(G/h_e) ≈ 0.5 mm, the same width as the hotspot. The 2D model resolves this spreading.",
        "A single unit cell with symmetry walls (assumption A12) would miss it near hotspot edges.",
        "",
        "Fig. 3 (h_e vs Δp for pitch and height sweeps), Fig. 6 (Adera, Zhu data) and Fig. A2",
        "(hotspot superheat) are planned digitizations (`data/literature/vaartstra2019_targets.csv`)",
        "and are not compared yet.",
        "",
        "Conclusion: the only independent comparison available today (dryout) disagrees by "
        f"{100 * (q_a / 82 - 1):+.0f}%. That is the known flat-meniscus bias, not a validated result.",
        "",
    ]

    md += [
        "## 6. Homogenized model vs explicit small-array COMSOL (pending)",
        "",
        "Nothing has been run yet. This is the only test of the homogenization itself (W2, W19), so it",
        "should come before trusting any graded-wick result. Proposed protocol:",
        "",
        "1. Explicit COMSOL strip: 2 mm × 0.3 mm (~65 × 10 pillars at l = 30 µm), real pillars, Laminar",
        "   Flow with a prescribed uniform evaporative sink on the meniscus. Start with a flat top (Stage 1",
        "   style), then add the Stage 3 meniscus. Feed at x = 0, symmetry on the long sides.",
        "2. Variant B: the same strip with a pitch step l = 30 → 22 µm at x = 1 mm (a graded transition).",
        "3. Homogenized model: identical domain, feed and sink, nx ≥ 20.",
        "4. Compare Pc(x) along the strip and the outlet-to-feed flow. Acceptance (proposed): |ΔPc| /",
        "   Pc_max < 5% more than 3 pitches from the transition. Report the width of the zone where it",
        "   exceeds 5%.",
        "5. If the transition error is large, refine the averaging rule (D5) or add a transition cell",
        "   type, and re-examine whether 2D is adequate (D6).",
        "",
        "## What would count as validation",
        "",
        "- Dryout heat flux and superheat vs heat flux against device experiments (adera2016 A1/A2,",
        "  zhu2016 B1/B2 from vaartstra2019 Fig. 6), with h_e and K(Pc) from COMSOL, not placeholders.",
        "- The explicit-array comparison of section 6 for a graded transition.",
        "- For electrowetting cases: measured receding angle vs voltage on the coated pillar surface",
        "  (θ_0, saturation, hysteresis), and a measured dryout shift under actuation.",
    ]
    (OUT / "VALIDATION_REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"wrote {OUT / 'VALIDATION_REPORT.md'}")


if __name__ == "__main__":
    main()
