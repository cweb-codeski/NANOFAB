"""Compare the analytical models against values reported by Vaartstra, Lu & Wang (IJHMT 2019).

Every printed number carries a provenance label (docs/SOURCE_AUTHORITY.md):
    PRIMARY_SOURCE    stated in a held paper (section given)
    MODEL_CALCULATED  produced by this repository's analytical models
    PROJECT_DERIVED   inferred by us from a paper's stated values
    ASSUMED           a choice not supported by a held source

The analytical models are simpler than the paper's (flat-meniscus permeability, uniform
evaporation, no unit-cell conduction), so they are NOT expected to match exactly.

Usage:
    python analysis/vaartstra_replication.py
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from _common import PROCESSED_DIR
from nanofab import capillary, dryout, evaporation, permeability
from nanofab.estimates import interface_area_fraction
from nanofab.fluids import Fluid
from nanofab.params import load_parameters
from nanofab.provenance import ASSUMED, MODEL_CALCULATED, PRIMARY_SOURCE, PROJECT_DERIVED

W_CM2 = 1e4

# Published values used below (PRIMARY_SOURCE). Parameters come from parameters/*.yaml.
V2019_HE_LOWER_BOUND = 60e3  # W/m^2K, vaartstra2019 Sec. 3.2: "h_e > 60 kW/m^2K"
ZHOU_ADERA_DRYOUT_D10_L30 = 84e4  # W/m^2, zhou_adera2026 abstract / Sec. 3.2 (d=10, l=30, h=25 um, 1 atm)
ADERA_CHAMBER_PRESSURE = 3000.0  # Pa, vaartstra2019 Sec. 3.4
ZHU_AMBIENT_T = 298.15  # K, ASSUMED: not stated in vaartstra2019 (air ambient, room temperature)


def label(provenance: str, text: str) -> None:
    print(f"   [{provenance}] {text}")


def vaartstra_eq3(gamma: float, d: float, l: float, theta_rec: float) -> float:
    """Paper Eq. (3): dp_dryout = 4 gamma cos(theta) / (d (4/pi (l/d)^2 - 1))."""
    return 4 * gamma * math.cos(theta_rec) / (d * (4 / math.pi * (l / d) ** 2 - 1))


def analytical_dryout(fluid: Fluid, T: float, d: float, l: float, h: float, L: float, theta: float) -> dict:
    s = fluid.saturated(T)
    Pc = capillary.pillar_max_capillary_pressure(s.surface_tension, d, l, theta)
    K = permeability.pillar_wick_permeability(d, l, h)
    q = dryout.dryout_heat_flux(K, h, L, Pc, s.rho_l, s.mu_l, s.h_lv)
    return {"Pc_dryout_Pa": Pc, "K_flat_m2": K, "q_dry_W_cm2": q / W_CM2}


def main() -> None:
    ps = load_parameters()
    fluid = Fluid(ps["working_fluid"])
    T = ps["vapor_temperature"]
    theta = ps["receding_contact_angle"]
    d, l, h, L = ps["pillar_diameter"], ps["pillar_pitch"], ps["pillar_height"], ps["wick_feed_length"]
    a = ps["accommodation_coefficient"]
    q_ref = ps["reference_dryout_heat_flux"]
    s = fluid.saturated(T)

    print("=" * 78)
    print("Vaartstra, Lu & Wang (2019) vs analytical models")
    print(f"inputs (all PRIMARY_SOURCE, vaartstra2019): water at {T-273.15:.0f} C, d/l/h = "
          f"{d*1e6:.0f}/{l*1e6:.0f}/{h*1e6:.0f} um, L = {L*1e3:.0f} mm, theta_rec = {math.degrees(theta):.0f} deg, "
          f"a = {a}")
    print("water properties: AUTHORITATIVE_DATABASE (IAPWS via CoolProp); vaartstra2019 does not state its source")
    print("=" * 78)

    # 1: baseline dryout
    base = analytical_dryout(fluid, T, d, l, h, L, theta)
    eq3 = vaartstra_eq3(s.surface_tension, d, l, theta)
    print("\n1) Dryout capillary pressure")
    label(PRIMARY_SOURCE, "equation: vaartstra2019 Eq. (3) (the paper gives no numeric value)")
    label(MODEL_CALCULATED, f"Eq. (3) evaluated: {eq3:.1f} Pa; nanofab.capillary: {base['Pc_dryout_Pa']:.1f} Pa")

    print("\n2) Dryout heat flux, baseline device")
    rel = base["q_dry_W_cm2"] * W_CM2 / q_ref - 1
    label(MODEL_CALCULATED, f"q''_dry = {base['q_dry_W_cm2']:.1f} W/cm^2 "
          "(1D Darcy, flat-meniscus K from Sangani-Acrivos + Brinkman, uniform evaporation)")
    label(PRIMARY_SOURCE, f"vaartstra2019 target = {q_ref/W_CM2:.0f} W/cm^2 (Sec. 3.2)")
    label(PRIMARY_SOURCE, f"independent cross-check: zhou_adera2026 = {ZHOU_ADERA_DRYOUT_D10_L30/W_CM2:.0f} W/cm^2 "
          "(same geometry, 1 atm; different model)")
    print(f"   relative error vs vaartstra2019: {100*rel:+.1f}%")
    K_implied = base["K_flat_m2"] / (1 + rel)
    label(MODEL_CALCULATED, f"K_flat = {base['K_flat_m2']:.3e} m^2; K that would make our 1D model give "
          f"82 W/cm^2 = {K_implied:.3e} m^2 ({100*(1-K_implied/base['K_flat_m2']):.0f}% lower)")
    print("   Interpretation (PROJECT_DERIVED hypothesis, not a paper statement): vaartstra2019 Sec. 2.2")
    print("   states that curvature lowers K, so part of the gap is the curved meniscus. The rest may be")
    print("   error in our analytical K model (coefficients not yet verified) or the paper's non-uniform")
    print("   evaporation. COMSOL Stages 3 and 5 will separate these.")

    # 3: height scaling at fixed ratios
    print("\n3) Height scaling at d/h = 0.4, l/h = 1.2, L = 5 mm")
    label(PRIMARY_SOURCE, "vaartstra2019 Sec. 3.3 / Fig. 5: q''_dry increases monotonically, 'scaling with h^2'")
    hs = np.array([5, 10, 20, 30, 40, 50]) * 1e-6
    rows = [analytical_dryout(fluid, T, 0.4 * hh, 1.2 * hh, hh, L, theta) for hh in hs]
    qd = np.array([r["q_dry_W_cm2"] for r in rows])
    slope = np.polyfit(np.log(hs), np.log(qd), 1)[0]
    for hh, r in zip(hs, rows):
        label(MODEL_CALCULATED, f"h = {hh*1e6:4.0f} um: q''_dry = {r['q_dry_W_cm2']:8.1f} W/cm^2")
    label(MODEL_CALCULATED, f"fitted exponent = {slope:.2f} (exactly 2 by geometric similarity in our model)")

    # 4: kinetic share of the unit-cell resistance
    print("\n4) Kinetic interfacial resistance vs the paper's total h_e")
    frac = interface_area_fraction(d, l)
    h_i_schrage = evaporation.interfacial_htc(T, a, fluid)
    dT = 0.1
    h_i_bte = evaporation.sone_bte_heat_flux(T + dT, T, a, fluid) / dT
    R_int = 1 / (h_i_bte * frac)
    label(MODEL_CALCULATED, f"h_i per interface area: linearized BTE {h_i_bte:.3e}, Schrage {h_i_schrage:.3e} W/m^2K "
          "(BTE form from vaartstra2019 Eqs. 5-7; accommodation correction form unverified)")
    label(MODEL_CALCULATED, f"R''_int per footprint = {R_int*1e6:.2f} mm^2K/W (ASSUMED flat interface area)")
    label(PRIMARY_SOURCE, f"vaartstra2019 Sec. 3.2: h_e > {V2019_HE_LOWER_BOUND/1e3:.0f} kW/m^2K "
          f"-> total 1/h_e < {1e6/V2019_HE_LOWER_BOUND:.1f} mm^2K/W")
    label(MODEL_CALCULATED, f"kinetic share of total unit-cell resistance >= ~{100*R_int*V2019_HE_LOWER_BOUND:.0f}% "
          "(conditional LOWER bound: R_int/R_total = R_int h_e > R_int x 60000 W/m^2K; valid only if our "
          "R''_int is right; not stated by the paper)")
    print("   Caveats: R''_int uses an ASSUMED flat interface area and an UNVERIFIED nonunity-accommodation")
    print("   correction (Meland & Ytrehus 2003, source gap). The fitted a = 0.052 is a single-paper value.")

    # 5: validation geometries
    print("\n5) Validation geometries (geometry: PRIMARY_SOURCE, vaartstra2019 Fig. 6 caption)")
    T_adera = fluid.T_sat(ADERA_CHAMBER_PRESSURE)
    label(PRIMARY_SOURCE, f"Adera chamber pressure {ADERA_CHAMBER_PRESSURE/1e3:.0f} kPa (Sec. 3.4)")
    label(MODEL_CALCULATED, f"T_sat(3 kPa) = {T_adera:.2f} K (IAPWS via CoolProp)")
    label(ASSUMED, f"Zhu air-ambient cases evaluated at {ZHU_AMBIENT_T:.2f} K (not stated in the paper)")
    cases = [
        ("A1 Adera", 5e-6, 12e-6, 82e-6, T_adera, MODEL_CALCULATED, "3 kPa vapor, fed all sides (2D flow)"),
        ("A2 Adera", 12e-6, 20e-6, 90e-6, T_adera, MODEL_CALCULATED, "3 kPa vapor, fed all sides (2D flow)"),
        ("B1 Zhu", 7e-6, 20e-6, 20e-6, ZHU_AMBIENT_T, ASSUMED, "air ambient, fed from one side"),
        ("B2 Zhu", 6e-6, 50e-6, 19e-6, ZHU_AMBIENT_T, ASSUMED, "air ambient, fed from one side"),
    ]
    out = []
    for name, dd, ll, hh, T_case, T_prov, note in cases:
        for L_case in (2e-3, 5e-3, 10e-3):
            r = analytical_dryout(fluid, T_case, dd, ll, hh, L_case, theta)
            out.append({
                "case": name, "d_um": dd * 1e6, "l_um": ll * 1e6, "h_um": hh * 1e6,
                "geometry_provenance": PRIMARY_SOURCE,
                "T_K": T_case, "T_provenance": T_prov,
                "L_mm": L_case * 1e3, "L_provenance": ASSUMED,
                **r, "result_provenance": MODEL_CALCULATED, "notes": note,
            })
    df = pd.DataFrame(out)
    show = df[["case", "d_um", "l_um", "h_um", "T_K", "L_mm", "Pc_dryout_Pa", "K_flat_m2", "q_dry_W_cm2"]]
    print("   [MODEL_CALCULATED] results (wicking lengths L are ASSUMED evaluation points):")
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("   Zhu cases: compare against Fig. 6(b) once digitized (digitized values will be PROJECT_DERIVED).")
    print("   Adera cases: a 1D L is only indicative (their wicks were fed from all sides).")
    print(f"\n   Note: {PROJECT_DERIVED} = inferred by us from a paper; see config/model_provenance.yaml.")

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    dest = PROCESSED_DIR / "vaartstra2019_analytical_comparison.csv"
    df.to_csv(dest, index=False)
    print(f"\nWrote {dest} (all result columns MODEL_CALCULATED)")


if __name__ == "__main__":
    main()
