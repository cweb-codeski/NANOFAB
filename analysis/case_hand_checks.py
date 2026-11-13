"""Hand-calculation equivalents for every COMSOL stage, side by side for the reference and proposed cells.

All numbers are MODEL_CALCULATED estimates for planning and for recognizing a wrong COMSOL result.
They are not predictions and must never be entered as COMSOL results. Inputs come from
parameters/unit_cell/cases.yaml (geometry, angle, T_v, q0), parameters/*.yaml (k_Si, accommodation)
and CoolProp (IAPWS water).

Sections
  A  geometry and capillary limit (V2019 Eq. 3), meniscus sag estimate, array extent
  B  Stage 1 expectations: K (shear-free and no-slip flat top), quarter-cell Q at dp = 1 Pa, Re
  C  capillary-limit scaling K*Pc_max and q''_dry*L^2 (feed length is UNRESOLVED for the proposed wick)
  D  Stage 2 expectations: Q_in, footprint, conduction bounds, interfacial resistance by closure
  E  Stage 2A expectations: F, 2D field-treatment limits, Debye length / gap screening

Usage:
    python analysis/case_hand_checks.py
    python analysis/case_hand_checks.py --feed-lengths-mm 2.5 5 10 --zeta-mV -50
"""

from __future__ import annotations

import argparse
import math

from _common import REPO  # noqa: F401  (sets sys.path)
from nanofab.capillary import pillar_max_capillary_pressure, spherical_cap_depth
from nanofab.conduction import parallel_conductivity, series_conductivity
from nanofab.dryout import dryout_heat_flux
from nanofab.electroosmosis import physics
from nanofab.electroosmosis.formation_factor import formation_factor_maxwell_2d
from nanofab.electroosmosis.planform2d import field_treatment
from nanofab.evaporation import interfacial_htc, sone_bte_mass_flux
from nanofab.fluids import Fluid
from nanofab.geometry import porosity, solid_fraction
from nanofab.geometry_cases import (
    UnresolvedInputError,
    array_extent,
    load_geometry_cases,
    load_screening_inputs,
    start_reference_interpretations,
)
from nanofab.params import load_parameters
from nanofab.permeability import pillar_wick_permeability, sangani_acrivos_permeability

CONCENTRATIONS_M = (1e-7, 1e-6, 1e-5, 1e-4, 1e-3)  # mol/L; 1e-7 ~ intrinsic ions of pure water at 25 C


def bte_htc(T: float, a: float, water: Fluid) -> float:
    """Slope of the reduced-BTE flux (ASSUMED Meland-Ytrehus form) at T, numerically [W/m2/K]."""
    dT = 1e-3
    return sone_bte_mass_flux(T + dT, T, a, water) * water.h_lv(T) / dT


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--feed-lengths-mm", type=float, nargs="+", default=[2.5, 5.0, 10.0],
                    help="ILLUSTRATIVE feed lengths for the proposed wick (its feed length is UNRESOLVED)")
    ap.add_argument("--zeta-mV", type=float, default=-50.0, help="ILLUSTRATIVE zeta potential for screening")
    ap.add_argument("--E-V-per-m", type=float, default=1e4, help="ILLUSTRATIVE mean field for screening")
    args = ap.parse_args()

    cases = load_geometry_cases()
    ps = load_parameters()
    scr = load_screening_inputs()
    water = Fluid(ps["working_fluid"])
    k_si, k_si_sweep = ps["silicon_thermal_conductivity"], ps.meta("silicon_thermal_conductivity").sweep
    acc = ps["accommodation_coefficient"]
    T_room = scr["room_temperature"].value

    print("All values MODEL_CALCULATED from the inputs listed in the docstring. Not COMSOL results.\n")
    for name, c in cases.items():
        d, l, hp = c.d, c.l, c.hp
        Tv = c.si("vapor_temperature")
        sat, room = water.saturated(Tv), water.saturated(T_room)
        theta = c.si("receding_contact_angle")
        phi, eps = solid_fraction(d, l), porosity(d, l)
        print("=" * 100)
        print(f"{name}  [{c.role}]  d/l/hp = {d*1e6:g}/{l*1e6:g}/{hp*1e6:g} um, "
              f"theta_rec = {math.degrees(theta):g} deg ({c.meta('receding_contact_angle').provenance})")
        print("=" * 100)

        # A. geometry and capillary limit
        pc_v = pillar_max_capillary_pressure(sat.surface_tension, d, l, theta)
        pc_r = pillar_max_capillary_pressure(room.surface_tension, d, l, theta)
        half_span = l / math.sqrt(2) - d / 2
        sag = spherical_cap_depth(2 * sat.surface_tension / pc_v, half_span)
        print("A  solid fraction phi = {:.4f}, porosity = {:.4f}, gap l-d = {:.1f} um".format(phi, eps, (l - d) * 1e6))
        print(f"   Pc_max (Eq. 3) = {pc_v:.0f} Pa at T_v = {Tv - 273.15:.0f} C; {pc_r:.0f} Pa at {T_room - 273.15:.0f} C")
        print(f"   meniscus sag at Pc_max, spherical-cap estimate over the diagonal gap = {sag*1e6:.2f} um "
              f"({sag / hp:.1%} of hp). Floor contact before Eq. 3 is not expected; Stage 3 must confirm.")
        if c.role == "PROPOSED_DESIGN":
            ext = array_extent(int(c.si("array_columns")), int(c.si("array_rows")), l, d)
            print(f"   array {int(c.si('array_columns'))} x {int(c.si('array_rows'))}: center span "
                  f"{ext.center_span_x*1e3:.3f} x {ext.center_span_y*1e3:.3f} mm, outer pillar-edge span "
                  f"{ext.edge_span_x*1e3:.3f} x {ext.edge_span_y*1e3:.3f} mm")
            for r in start_reference_interpretations(c):
                print(f"     start reference read as {r['reading']}: x {r['x_min']*1e3:+.3f}..{r['x_max']*1e3:+.3f} mm, "
                      f"y {r['y_min']*1e3:+.3f}..{r['y_max']*1e3:+.3f} mm; offset from coupon center "
                      f"({r['offset_from_coupon_center_x']*1e3:+.2f}, {r['offset_from_coupon_center_y']*1e3:+.2f}) mm")
            print("     -> neither reading centers the array on the coupon; CAD scale/units are UNRESOLVED.")

        # B. Stage 1
        K_s = sangani_acrivos_permeability(d, l)
        K_sf = pillar_wick_permeability(d, l, hp, top="shear_free")
        K_ns = pillar_wick_permeability(d, l, hp, top="no_slip")
        mu = sat.mu_l
        Q1 = K_sf * hp * 1.0 / mu
        G1 = 1.0 / (l / 2)
        U1 = K_sf * G1 / mu
        print(f"B  K_s (infinite array) = {K_s:.3e}; K flat shear-free = {K_sf:.3e}; K flat no-slip top = {K_ns:.3e} m^2 "
              f"(no-slip/shear-free = {K_ns / K_sf:.3f})")
        print(f"   Stage 1 quarter cell, dp = 1 Pa over l/2, mu = {mu:.4e} Pa s: Q ~ {Q1:.3e} m^3/s "
              f"(analytical; COMSOL is expected within a few %)")
        print(f"   Stage 1 test gradient {G1:.2e} Pa/m: U = {U1:.2e} m/s, Re_h = {sat.rho_l*U1*hp/mu:.2f}, "
              f"Re_K = {sat.rho_l*U1*math.sqrt(K_sf)/mu:.3f}")

        # C. capillary-limit scaling
        try:
            Ls = [c.si("feed_length")]
            tag = f"{c.meta('feed_length').provenance}"
        except UnresolvedInputError:
            Ls = [x * 1e-3 for x in args.feed_lengths_mm]
            tag = "ILLUSTRATIVE (feed length UNRESOLVED)"
        qL2 = dryout_heat_flux(K_sf, hp, 1.0, pc_v, sat.rho_l, mu, sat.h_lv)
        qL2_ns = dryout_heat_flux(K_ns, hp, 1.0, pc_v, sat.rho_l, mu, sat.h_lv)
        print(f"C  K*Pc_max = {K_sf*pc_v:.3e} N (shear-free); q''_dry*L^2 = {qL2:.3e} W (shear-free), "
              f"{qL2_ns:.3e} W (no-slip top) [flat K, uniform flux, one feed edge, closed far end]")
        for L in Ls:
            Gdev = 2 * pc_v / L
            Udev = K_sf * Gdev / mu
            print(f"   L = {L*1e3:g} mm [{tag}]: q''_dry = {qL2/L**2/1e4:.0f} W/cm^2 (no-slip top {qL2_ns/L**2/1e4:.0f}); "
                  f"feed gradient at the limit {Gdev:.2e} Pa/m -> Re_h = {sat.rho_l*Udev*hp/mu:.1f}, "
                  f"Re_K = {sat.rho_l*Udev*math.sqrt(K_sf)/mu:.2f}")

        # D. Stage 2
        q0 = c.si("unit_cell_heat_flux")
        A_fp = (l / 2) ** 2
        A_lv = A_fp - math.pi * d**2 / 16
        t_sub = c.si("unit_cell_substrate_thickness")
        k_l = sat.k_l
        R_sub = t_sub / k_si
        R_par = hp / parallel_conductivity(phi, k_si, k_l)
        R_ser = hp / series_conductivity(phi, k_si, k_l)
        print(f"D  A_fp = (l/2)^2 = {A_fp:.4e} m^2 (Surface Integration of 1 over the T_base selection must equal this); "
              f"flat liquid interface area = {A_lv:.4e} m^2; Q_in = q0*A_fp = {q0*A_fp:.4e} W")
        print(f"   pillar volume in the quarter cell = pi d^2 hp / 16 = {math.pi*d**2*hp/16*1e18:.1f} um^3 "
              "(Geometry > Measure; a larger value means the cylinder was not trimmed to the cell)")
        print(f"   R'' substrate {t_sub*1e6:g} um = {R_sub*1e6:.4f}; wick parallel bound = {R_par*1e6:.2f}; "
              f"series bound = {R_ser*1e6:.1f} mm^2K/W (k_Si = {k_si:g}, ASSUMED; sweep {k_si_sweep})")
        for label, hi in (("reduced BTE, ASSUMED Meland-Ytrehus form, a = %.3g" % acc, bte_htc(Tv, acc, water)),
                          ("Schrage, a = %.3g" % acc, interfacial_htc(Tv, acc, water)),
                          ("Schrage, a = 1", interfacial_htc(Tv, 1.0, water))):
            R_int = A_fp / (hi * A_lv)
            print(f"   h_i [{label}] = {hi:.3e} W/m^2K -> footprint R''_int = {R_int*1e6:.3f} mm^2K/W; "
                  f"T_base - T_v between {q0*(R_sub+R_par+R_int):.2f} and {q0*(R_sub+R_ser+R_int):.1f} K")

        # E. Stage 2A
        ft = field_treatment(d, l, n_modes=24)
        F2 = 1.0 / ft.inv_F
        print(f"E  F (Maxwell 2D) = {formation_factor_maxwell_2d(phi):.4f}; F (2D collocation) = {F2:.4f}. "
              "Flat top + insulating floor/top -> field is z-invariant, so the 3D F equals the 2D F.")
        print(f"   Exact Case B check: U_sup/(mu_EO E0) = 1/F = {ft.inv_F:.4f}. "
              f"2D tall-pillar limit Q_uniform/Q_solved = {ft.q_uniform/ft.q_solved:.3f} (floor adds flux in 3D)")
        a_half = (l - d) / 2
        zeta = args.zeta_mV * 1e-3
        mob = physics.electroosmotic_mobility(scr["water_relative_permittivity_298K"].value, zeta, room.mu_l)
        U_eo = mob * args.E_V_per_m / F2
        print(f"   ILLUSTRATIVE at {T_room-273.15:.0f} C, zeta = {args.zeta_mV:g} mV, E = {args.E_V_per_m:.0e} V/m: "
              f"mu_EO = {mob:.2e} m^2/Vs, U_EO = {U_eo:.2e} m/s, equivalent gradient mu U_EO/K = "
              f"{room.mu_l*U_eo/K_sf:.2e} Pa/m")
        print(f"   Debye length at {T_room-273.15:.0f} C (eps_r {scr['water_relative_permittivity_298K'].value:g}, ASSUMED) "
              f"vs half-gap {a_half*1e6:.0f} um:")
        for cM in CONCENTRATIONS_M:
            lam = physics.debye_length(scr["water_relative_permittivity_298K"].value, T_room, cM * 1e3)
            Du = physics.dukhin_number_lower_bound(lam, zeta, T_room, a_half)
            corr = physics.thin_edl_flux_factor(lam, zeta, T_room, a_half)
            print(f"     c = {cM:.0e} M: lambda_D = {lam*1e9:7.1f} nm, lambda_D/a = {lam/a_half:.1e}, "
                  f"finite-EDL flux factor = {corr:.3f}, Du lower bound = {Du:.2e}")
        print()

    print("Interpretation limits: analytical K uses dilute-array Sangani-Acrivos + Brinkman (unverified coefficients);")
    print("the sag estimate assumes a spherical cap; Debye-length screening ignores CO2, buffer ions and temperature;")
    print("zeta and E are illustrative. Use these only to recognize order-of-magnitude errors in COMSOL output.")


if __name__ == "__main__":
    main()
