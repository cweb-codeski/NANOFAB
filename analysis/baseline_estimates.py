"""Print the analytical sanity-check table for the baseline micropillar wick.

Usage:
    python analysis/baseline_estimates.py
    python analysis/baseline_estimates.py --psat-table     # also write p_sat(T) for COMSOL

What it reports (all from parameters/*.yaml + CoolProp):
    - capillary pressure, permeability, capillary-limited dryout heat flux
    - dimensionless numbers that justify COMSOL modeling choices (Stokes flow, no gravity, ...)
    - area-specific thermal-resistance breakdown (substrate / wick / interface)
    - sensitivity of the interfacial resistance to the accommodation coefficient
    - the capillary-pressure vs permeability tradeoff over pillar pitch
    - the expected flow rate for the COMSOL Stage 1 quarter-cell permeability model
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from nanofab.estimates import (  # noqa: E402
    accommodation_sensitivity,
    compute_baseline,
    pitch_tradeoff,
    unit_cell_flow_check,
)
from nanofab.evaporation import comsol_bte_expression, comsol_schrage_expression  # noqa: E402
from nanofab.fluids import Fluid  # noqa: E402
from nanofab.params import load_parameters  # noqa: E402
from nanofab.provenance import ASSUMED  # noqa: E402

W_CM2 = 1e4  # W/m^2 per W/cm^2


def flag(ok: bool) -> str:
    return "ok" if ok else "CHECK"


def write_psat_table(fluid: Fluid, T_center: float, path: Path) -> None:
    T = np.arange(T_center - 20, T_center + 80.0001, 0.5)
    p = [fluid.p_sat(t) for t in T]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"% p_sat of {fluid.name} from CoolProp. Columns: T (K), p_sat (Pa)\n")
        for t, pv in zip(T, p):
            fh.write(f"{t:.2f}\t{pv:.6g}\n")
    print(f"\nWrote COMSOL interpolation table: {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", type=Path, default=None, help="parameter directory (default: parameters/)")
    ap.add_argument("--psat-table", action="store_true", help="write data/processed/psat_<fluid>.txt")
    args = ap.parse_args()

    ps = load_parameters(args.params)
    fluid = Fluid(ps["working_fluid"])
    b = compute_baseline(ps, fluid)
    props = fluid.saturated(b.T_v)

    print("=" * 78)
    print("NANOFAB baseline analytical estimates")
    print("=" * 78)

    print("Input provenance (parameters/*.yaml):")
    for category in sorted({p.provenance for p in ps.parameters.values()}):
        names = [p.name for p in ps.by_provenance(category)]
        print(f"  {category:24s} {len(names):2d}: {', '.join(names)}")
    if ps.by_provenance(ASSUMED):
        print("  -> ASSUMED inputs must be swept, not trusted as single values.")
    print("Every number printed below is MODEL_CALCULATED (this repo's analytical models) unless it")
    print("is explicitly labeled otherwise. None of them are literature values.")

    print("\n-- Geometry and state --")
    print(f"pillar d / l / h          : {b.d*1e6:.1f} / {b.l*1e6:.1f} / {b.h*1e6:.1f} um")
    print(f"solid fraction phi        : {b.phi:.3f}")
    print(f"fluid / vapor T / p_v     : {fluid.name} / {b.T_v-273.15:.1f} C / {b.p_v:.0f} Pa")
    print(f"gamma, mu, rho_l, h_lv    : {props.surface_tension:.4f} N/m, {props.mu_l:.3e} Pa s, "
          f"{props.rho_l:.1f} kg/m3, {props.h_lv/1e6:.3f} MJ/kg")

    print("\n-- Capillary pumping and dryout --")
    print(f"P_c,max (V2019 Eq. 3)     : {b.Pc_max/1e3:.2f} kPa")
    print(f"K_s (Sangani-Acrivos)     : {b.K_s:.3e} m^2")
    print(f"K_eff (Brinkman, height)  : {b.K_eff:.3e} m^2   (K_eff/K_s = {b.K_eff/b.K_s:.3f})")
    print(f"feed length L             : {ps['wick_feed_length']*1e3:.2f} mm")
    print(f"q''_dry (1D capillary)    : {b.q_dry/W_CM2:.1f} W/cm^2")
    print(f"target q'' [DESIGN_TARGET]: {b.q_target/W_CM2:.1f} W/cm^2   "
          f"-> margin q''_dry/q'' = {b.q_dry/b.q_target:.2f}  [{flag(b.q_dry > b.q_target)}]")
    print(f"wick dP at target q''     : {b.dP_at_target/1e3:.2f} kPa  (vs P_c,max {b.Pc_max/1e3:.2f} kPa)")

    print("\n-- Dimensionless checks (at target q'', at the feed where velocity is largest) --")
    print(f"superficial / pore vel.   : {b.U_feed*1e3:.3f} / {b.u_pore*1e3:.3f} mm/s")
    print(f"Re_h (V2019 def., U h)    : {b.Re_h:.2e}   laminar but not creeping -> use COMSOL Laminar Flow, as V2019"
          f" does [{flag(b.Re_h < 100)}]")
    print("                            (PRIMARY_SOURCE: Re_h ~ 10^1, Vaartstra thesis p. 23, 51)")
    print(f"Re (pore vel., gap scale) : {b.Re_pore:.2e}   Creeping Flow exact only if << 1 [{flag(b.Re_pore < 0.1)}]")
    print(f"Ca                        : {b.Ca:.2e}   static meniscus OK if << 1 [{flag(b.Ca < 1e-2)}]")
    print(f"Pe_h (vertical, V2019)    : {b.Pe_vertical:.2e}   convection toward interface negligible "
          f"[{flag(b.Pe_vertical < 0.1)}]")
    print(f"Ja (dT = 10 K, V2019)     : {b.Ja:.2e}   liquid enters at T_sat, lateral convection negligible "
          f"[{flag(b.Ja < 0.05)}]")
    print(f"Pe (lateral, pore vel.)   : {b.Pe_pore:.2e}   (large values are harmless while Ja << 1: the liquid"
          " carries almost no lateral temperature gradient)")
    print(f"Bo (pitch)                : {b.Bo_pitch:.2e}   gravity negligible at cell scale [{flag(b.Bo_pitch < 0.1)}]")
    print(f"rho g L / P_c,max         : {b.gravity_over_Pc_vertical:.2f}   (only matters if wick is vertical)")
    kn_note = (
        "continuum vapor OK" if b.Kn_pitch < 0.01
        else "vapor not continuum at cell scale -> use kinetic Schrage BC, don't mesh the vapor"
    )
    print(f"vapor mean free path      : {b.mean_free_path*1e6:.2f} um, Kn(pitch) = {b.Kn_pitch:.2f}  ({kn_note})")

    print("\n-- Thermal resistance breakdown (area-specific, per footprint) --")
    rows = [
        ("substrate conduction (t/k_Si)", b.R_substrate),
        ("wick conduction, parallel bound", b.R_wick_parallel),
        ("wick conduction, series bound", b.R_wick_series),
        (f"interface kinetic (s={b.accommodation:g})", b.R_interface),
    ]
    for name, R in rows:
        print(f"{name:34s}: {R*1e6:9.3f} mm^2K/W -> dT = {R*b.q_target:7.2f} K at target q''")
    print(f"{'total (optimistic, parallel wick)':34s}: {b.R_total_optimistic*1e6:9.3f} mm^2K/W -> "
          f"dT = {b.dT_at_target:7.2f} K")
    print(f"h_evap upper bound (wick+interface): {b.h_evap_upper_bound/1e4:.1f} W/cm^2K")
    print("  The true wick resistance lies between the bounds; resolving it is the job of the COMSOL")
    print("  unit cell (Stage 2-5), because heat reaches the interface mostly through the pillars and")
    print("  the thin liquid film near the contact line.")

    print("\n-- Accommodation-coefficient sensitivity (interfacial resistance) --")
    print(f"{'s':>6} {'h_i [W/m^2K]':>14} {'R_int [mm^2K/W]':>17} {'dT_int [K]':>11}")
    for s, h_i, R in accommodation_sensitivity(b.T_v, fluid, b.d, b.l):
        print(f"{s:6.2f} {h_i:14.3e} {R*1e6:17.3f} {R*b.q_target:11.2f}")

    print("\n-- Capillary pressure vs permeability tradeoff (vary pitch, fixed d, h, L) --")
    print(f"{'l [um]':>8} {'P_c [kPa]':>10} {'K [m^2]':>11} {'q_dry [W/cm2]':>14}")
    lo, hi = ps.meta("pillar_pitch").sweep or (1.5 * b.d, 10 * b.d)
    for l, Pc, K, qd in pitch_tradeoff(ps, fluid, np.linspace(max(lo, 1.2 * b.d), hi, 8)):
        print(f"{l*1e6:8.1f} {Pc/1e3:10.2f} {K:11.3e} {qd/W_CM2:14.1f}")

    dP_cell = 1.0
    Q = unit_cell_flow_check(b.K_eff, b.l, b.h, props.mu_l, dP_cell)
    print("\n-- COMSOL Stage 1 expectation (quarter cell, flat slip top) --")
    print(f"apply dP = {dP_cell:.0f} Pa across length l/2 = {b.l/2*1e6:.1f} um -> expect outlet flow "
          f"Q ~ {Q:.3e} m^3/s (from K_eff = {b.K_eff:.3e} m^2)")
    print(f"Re in that model ~ {props.rho_l * (Q/((b.l/2)*b.h)) / (1-b.phi) * (b.l-b.d) / props.mu_l:.2e}")

    print("\n-- COMSOL Stage 4 evaporative heat flux on the meniscus (define as variable q_evap) --")
    print("linearized BTE (V2019 Eqs. 5-7, paper-exact):")
    print("  q_evap = " + comsol_bte_expression(b.accommodation, fluid.molar_mass, b.T_v, props.h_lv))
    print("Schrage (alternative; within a few % of BTE at small accommodation):")
    print("  q_evap = " + comsol_schrage_expression(b.accommodation, fluid.molar_mass, b.T_v, props.h_lv)[1:])
    print("Boundary condition: Heat Flux, General inward heat flux q0 = -q_evap")
    print("(requires an Interpolation function psat(T); create it with --psat-table)")

    if args.psat_table:
        write_psat_table(fluid, b.T_v, REPO / "data" / "processed" / f"psat_{fluid.name.lower()}.txt")


if __name__ == "__main__":
    main()
