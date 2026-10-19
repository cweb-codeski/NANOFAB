"""Passive vs electroosmotically assisted micropillar wick: benchmarks, sweeps, falsification, plots.

Case A: passive wick. Case B: the SAME wick object plus an in-plane field (mu_EO, E). K, F, Pc_max,
geometry and liquid are shared; only the forcing differs (enforced by wick1d.assert_identical...).

Usage:
    python analysis/eo_assisted_wick.py                 # everything -> results/electroosmosis/
    python analysis/eo_assisted_wick.py --quick         # coarser formation-factor grids

Every number written is MODEL_CALCULATED from inputs in parameters/electroosmosis/eo_inputs.yaml
(provenance in results/electroosmosis/parameter_provenance.csv). Electrokinetic values come from
pujahari2022 (a simulation paper), not from measurements on our device.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import replace
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from _common import PROCESSED_DIR, REPO  # noqa: F401  (sets sys.path)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nanofab.capillary import pillar_max_capillary_pressure  # noqa: E402
from nanofab.electroosmosis import physics as ek  # noqa: E402
from nanofab.electroosmosis import wick1d as w1  # noqa: E402
from nanofab.electroosmosis.formation_factor import formation_factor_maxwell_2d, formation_factor_numerical  # noqa: E402
from nanofab.electroosmosis.inputs import (  # noqa: E402
    baseline_liquid,
    baseline_wick,
    load_eo_inputs,
    scenario_mobilities,
)
from nanofab.geometry import solid_fraction  # noqa: E402
from nanofab.params import load_parameters  # noqa: E402
from nanofab.permeability import pillar_wick_permeability  # noqa: E402

OUT = REPO / "results" / "electroosmosis"
W_CM2 = 1e4
SC = "literature_scenarios.pujahari2022"
COLORS = {"A": "0.15", -0.1: "#1f77b4", -0.2: "#d62728"}
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3, "savefig.dpi": 200,
                     "figure.dpi": 100, "axes.titlesize": 10.5, "legend.fontsize": 8.5})

FOOT_BASE = ("Geometry d/l/h = 10/30/25 um, L = 5 mm (vaartstra2019, PRIMARY_SOURCE); K = 2.97e-11 m2 "
             "(COMSOL Stage 1, flat top); water at 100 C (IAPWS via CoolProp).")
FOOT_EK = ("zeta = -100 / -200 mV (1e-4 / 1e-6 M KCl), eps_r = 55.4, E <= 1e5 V/m: pujahari2022 "
           "(simulation inputs, PRIMARY_SOURCE for 'stated in paper', NOT measured on our device).")


def footer(fig, *lines, y=0.005):
    fig.text(0.01, y, "\n".join(lines), fontsize=6.3, color="0.35", va="bottom", ha="left")


def save(fig, name, bottom=0.17):
    fig.subplots_adjust(bottom=bottom)
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)


def write_csv(name, rows):
    pd.DataFrame(rows).to_csv(OUT / f"{name}.csv", index=False)


# ---------------------------------------------------------------------------------------------
# 1. passive hydraulic benchmark
# ---------------------------------------------------------------------------------------------

def k_mesh_extrapolation(df: pd.DataFrame) -> dict:
    """K(h) = K_inf - C h^p through the three finest meshes (non-uniform refinement ratio)."""
    d = df.sort_values("hmax")
    h = d["hmax"].to_numpy()[:3]
    K = d["K_comsol_m2"].to_numpy()[:3]  # fine -> coarse: h[0] finest
    target = (K[2] - K[1]) / (K[1] - K[0])

    def f(p):
        return (h[2] ** p - h[1] ** p) / (h[1] ** p - h[0] ** p) - target

    p = brentq(f, 0.2, 6.0)
    C = (K[0] - K[1]) / (h[1] ** p - h[0] ** p)
    K_inf = K[0] + C * h[0] ** p
    return {"order": p, "K_inf": K_inf, "K_finest": K[0], "rel_extrapolation": (K_inf - K[0]) / K[0]}


def passive_benchmark(inp, liquid, wick, quick: bool) -> dict:
    df = pd.read_csv(PROCESSED_DIR / "permeability_stage1_flow.csv")
    mesh = k_mesh_extrapolation(df)
    K_an = pillar_wick_permeability(wick.d, wick.l, wick.h)
    levels = (64, 128, 256) if quick else (128, 256, 512)
    Fst = formation_factor_numerical(wick.d, wick.l, levels)
    A = w1.Case("A", wick, liquid, w1.Forcing())
    q_dry_comsol = w1.passive_dryout_heat_flux(A)
    q_dry_an = w1.passive_dryout_heat_flux(replace(A, wick=replace(wick, K=K_an)))
    q_dry_fv = w1.dryout_heat_flux_fv(A, 400)
    ps = load_parameters()
    out = {
        "K_comsol_finest_m2": mesh["K_finest"], "K_comsol_mesh_extrapolated_m2": mesh["K_inf"],
        "K_mesh_observed_order": mesh["order"], "K_mesh_extrapolation_rel": mesh["rel_extrapolation"],
        "K_analytical_m2": K_an, "K_comsol_vs_analytical_rel": (wick.K - K_an) / K_an,
        "F_maxwell": Fst.maxwell_2d, "F_numerical_extrapolated": Fst.F_extrapolated,
        "F_numerical_levels": list(Fst.n), "F_numerical_values": list(Fst.F), "F_observed_order": Fst.observed_order,
        "F_numerical_vs_maxwell_rel": Fst.rel_diff_vs_maxwell,
        "Pc_max_Pa": wick.Pc_max, "phi": wick.phi,
        "q_dry_passive_comsolK_Wcm2": q_dry_comsol / W_CM2, "q_dry_passive_fv_Wcm2": q_dry_fv / W_CM2,
        "q_dry_passive_analyticalK_Wcm2": q_dry_an / W_CM2,
        "q_dry_vaartstra2019_Wcm2 [PRIMARY_SOURCE]": ps["reference_dryout_heat_flux"] / W_CM2,
        "U_cap_m_s": wick.K * wick.Pc_max / (liquid.mu * wick.L),
    }

    # Figure: benchmark (K mesh, Q vs grad p linear, Pc profile vs 2D-equivalent FV)
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.9))
    a = ax[0]
    d = df.sort_values("hmax")
    a.plot(d["hmax"] * 1e6, d["K_comsol_m2"] * 1e11, "o-", label="COMSOL Stage 1")
    a.axhline(mesh["K_inf"] * 1e11, ls="--", c="k", lw=1, label=f"mesh-extrapolated ({mesh['order']:.2f} order)")
    a.axhline(K_an * 1e11, ls=":", c="tab:green", label="Sangani-Acrivos + Brinkman")
    a.set_xlabel("max element size hmax [um]")
    a.set_ylabel("intrinsic permeability K [1e-11 m$^2$]")
    a.set_title("(a) K mesh convergence (geometric property)")
    a.legend()
    a = ax[1]
    G = np.linspace(0, wick.Pc_max / wick.L, 50)
    r = w1.forced_flow(A, G * wick.L)
    a.plot(G / 1e5, r["Q_total"] * 1e9, c=COLORS["A"], lw=2, label="Case A: Darcy, Q = K A dP/(mu L)")
    a.plot(G[::7] / 1e5, wick.K * wick.area * G[::7] / liquid.mu * 1e9, "o", mfc="none", c="tab:orange",
           label="hand check K A G / mu")
    a.set_xlabel("imposed pressure gradient -dp/dx [1e5 Pa/m]")
    a.set_ylabel("volumetric flow Q [mm$^3$/s]")
    a.set_title("(b) passive forced flow: linear, Q(0) = 0")
    a.legend()
    a = ax[2]
    for q, ls in ((50e4, "-"), (q_dry_comsol, "--")):
        s = w1.solve_fv(A, q, 200)
        xx = np.linspace(0, wick.L, 200)
        a.plot(xx * 1e3, w1.capillary_pressure_profile(A, q, xx), ls, c="0.2", label=f"closed form, {q/W_CM2:.0f} W/cm$^2$")
        a.plot(s.x[::10] * 1e3, s.Pc[::10], "o", ms=3.5, c="tab:purple", mfc="none")
    a.axhline(wick.Pc_max, c="tab:red", lw=1, label=f"Pc,max = {wick.Pc_max:.0f} Pa (V2019 Eq. 3)")
    a.plot([], [], "o", c="tab:purple", mfc="none", ms=3.5, label="independent FV solve")
    a.set_xlabel("distance from feed x [mm]")
    a.set_ylabel("capillary pressure p$_v$ - p$_l$ [Pa]")
    a.set_title(f"(c) passive evaporating wick: dryout {q_dry_comsol/W_CM2:.1f} W/cm$^2$")
    a.legend(loc="upper left")
    fig.tight_layout()
    footer(fig, FOOT_BASE,
           f"K: finest-mesh change 0.10%, extrapolated +{100*mesh['rel_extrapolation']:.2f}%, analytical "
           f"{100*(K_an-wick.K)/wick.K:+.1f}%. Dryout is the uniform-flux flat-meniscus 1D capillary limit "
           f"(MODEL_CALCULATED); vaartstra2019 reports 82 W/cm$^2$ for this geometry (PRIMARY_SOURCE).")
    save(fig, "fig1_passive_benchmark", bottom=0.2)
    return out


# ---------------------------------------------------------------------------------------------
# 2-5. transport comparisons
# ---------------------------------------------------------------------------------------------

def transport_figures(inp, liquid, wick, scen) -> list[dict]:
    A = w1.Case("A", wick, liquid, w1.Forcing())
    Gmax = wick.Pc_max / wick.L
    G = np.linspace(0, Gmax, 120)
    E_list = list(inp.si(f"{SC}.electric_field"))
    rows = []

    # Fig 2: Q vs pressure gradient, passive vs assisted
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for a, s in zip(ax, scen):
        rA = w1.forced_flow(A, G * wick.L)
        lo = rA["Q_total"] * 0.63
        a.fill_between(G / 1e5, lo * 1e9, rA["Q_total"] * 1e9, color="0.85", label="Case A, K band 0.63-1.0 K (ASSUMED)")
        a.plot(G / 1e5, rA["Q_total"] * 1e9, c=COLORS["A"], lw=2.2, label="Case A: passive")
        for E, ls in zip(E_list, (":", "--", "-")):
            B = w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], E))
            w1.assert_identical_except_forcing(A, B)
            rB = w1.forced_flow(B, G * wick.L)
            a.plot(G / 1e5, rB["Q_total"] * 1e9, ls, c=COLORS[round(s["zeta"], 3)], lw=1.8,
                   label=f"Case B: E = {E:.0e} V/m")
            for g, qp, qe, qt in zip(G[::12], rB["Q_pressure"][::12], rB["Q_eo"][::12], rB["Q_total"][::12]):
                rows.append({"zeta_V": s["zeta"], "E_V_per_m": E, "minus_dpdx_Pa_per_m": g, "Q_pressure_m3_s": qp,
                             "Q_eo_m3_s": qe, "Q_total_m3_s": qt, "Q_total_over_Q_passive": qt / qp if qp else np.inf})
        a.set_title(f"{s['label']}")
        a.set_xlabel("imposed pressure gradient -dp/dx [1e5 Pa/m]  (right end = Pc,max / L)")
        a.legend(loc="upper left")
    ax[0].set_ylabel("Q_total [mm$^3$/s]  (W x h = 10 mm x 25 um)")
    fig.suptitle("Forced flow through the SAME wick: Q_total = K A (-dp/dx)/mu + (mu_EO E / F) A", y=0.99)
    fig.tight_layout()
    footer(fig, FOOT_BASE, FOOT_EK,
           "Pressure-driven part identical in A and B (K unchanged). EO part is the zero-gradient intercept. F = 1.191 (Maxwell 2D, verified numerically).")
    save(fig, "fig2_Q_vs_pressure_gradient", bottom=0.2)

    # Fig 3: EO contribution and velocities vs E (linear axes, literature range only)
    E = np.linspace(0, 1e5, 101)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    U_c = 2 * wick.K * wick.Pc_max / (liquid.mu * wick.L)
    U50 = float(w1.feed_velocity(A, 50e4))
    for s in scen:
        c = COLORS[round(s["zeta"], 3)]
        v = [w1.eo_velocities(w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], e))) for e in E]
        ax[0].plot(E / 1e3, [x["hs_slip"] * 1e3 for x in v], ":", c=c, label=f"HS slip mu_EO E, {s['label']}")
        ax[0].plot(E / 1e3, [x["intrinsic_pore_average"] * 1e3 for x in v], "--", c=c, label="pore-average U_EO/(1-phi)")
        ax[0].plot(E / 1e3, [x["superficial"] * 1e3 for x in v], "-", c=c, lw=2, label="superficial U_EO = mu_EO E/F")
        for f, ls in ((0.0, "-"), (0.5, "--"), (1.0, ":")):
            Qp = float(w1.pressure_superficial_velocity(A, f * wick.Pc_max / wick.L)) * wick.area
            ax[1].plot(E / 1e3, (Qp + np.array([x["superficial"] for x in v]) * wick.area) * 1e9, ls, c=c,
                       lw=2 if f == 0 else 1.3, label=f"{s['zeta']*1e3:.0f} mV, -dp/dx = {f:g} Pc,max/L")
    ax[0].axhline(U_c * 1e3, c="k", lw=1)
    ax[0].text(2, U_c * 1e3 * 1.02, "superficial feed velocity at passive dryout", fontsize=7.5)
    ax[0].axhline(U50 * 1e3, c="k", lw=1, ls="-.")
    ax[0].text(2, U50 * 1e3 * 1.02, "feed velocity needed at 50 W/cm$^2$", fontsize=7.5)
    ax[0].set_xlabel("applied field E [kV/m]  (pujahari2022 range 5-100 kV/m)")
    ax[0].set_ylabel("electroosmotic liquid velocity [mm/s]")
    ax[0].set_title("(a) three distinct EO velocities (do not confuse)")
    ax[0].set_ylim(0, U_c * 1e3 * 1.15)
    ax[0].legend(fontsize=6.6, ncol=2, loc="upper left", bbox_to_anchor=(0.0, 0.85))
    ax[1].set_xlabel("applied field E [kV/m]   (over 5 mm: 100 kV/m = 500 V)")
    ax[1].set_ylabel("Q_total [mm$^3$/s]  (solid: Q_EO alone, zero gradient)")
    ax[1].set_title("(b) Q_total vs E at fixed -dp/dx (E = 0 is Case A)")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    footer(fig, FOOT_BASE, FOOT_EK, "Superficial = per total cross-section (Darcy). HS slip = velocity outside a thin EDL in an open channel.")
    save(fig, "fig3_eo_contribution_vs_E", bottom=0.2)

    # Fig 4: enhancement ratio vs E at fixed gradients + apparent permeability
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    fracs = (0.1, 0.5, 1.0)
    for s in scen:
        c = COLORS[round(s["zeta"], 3)]
        for f, ls in zip(fracs, (":", "--", "-")):
            ratio = [1 + w1.eo_superficial_velocity(w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], e)))
                     / float(w1.pressure_superficial_velocity(A, f * Gmax)) for e in E]
            ax[0].plot(E / 1e3, ratio, ls, c=c, label=f"{s['zeta']*1e3:.0f} mV, -dp/dx = {f:g} Pc,max/L")
        Gs = np.logspace(3, math.log10(Gmax), 100)
        for e, ls in ((1e4, "--"), (1e5, "-")):
            B = w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], e))
            ax[1].loglog(Gs, w1.apparent_permeability(B, Gs) / wick.K, ls, c=c, label=f"{s['zeta']*1e3:.0f} mV, E = {e:.0e} V/m")
    ax[0].set_xlabel("applied field E [kV/m]")
    ax[0].set_ylabel("Q_total / Q_passive  (same wick, same dP)")
    ax[0].set_title("(a) transport enhancement ratio = 1 + U_EO / U_pressure")
    ax[0].legend(fontsize=7)
    ax[1].axhline(1, c="k", lw=1.5, label="intrinsic K (both cases)")
    ax[1].set_xlabel("imposed pressure gradient -dp/dx [Pa/m]")
    ax[1].set_ylabel("apparent K_app / K = mu U_total / (K (-dp/dx))")
    ax[1].set_title("(b) 'effective permeability' is NOT a property: diverges as dp/dx -> 0")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    footer(fig, FOOT_BASE, FOOT_EK, "K_app is a derived apparent response; K itself is the COMSOL geometric value in both cases.")
    save(fig, "fig4_enhancement_ratio", bottom=0.2)
    return rows


# ---------------------------------------------------------------------------------------------
# 6. sensitivity, falsification and regime maps
# ---------------------------------------------------------------------------------------------

def geometry_wick(d, l, h, L, W, theta, gamma):
    return w1.Wick.analytical(d, l, h, L, W, theta, gamma)


def falsification(inp, liquid, wick, scen) -> dict:
    theta0 = inp.si("wick.receding_contact_angle")
    mob = {round(s["zeta"], 3): s["mobility"] for s in scen}
    U_cap0 = wick.K * wick.Pc_max / (liquid.mu * wick.L)
    lit_E = inp.si(f"{SC}.electric_field")
    Emin, Emax = float(min(lit_E)), float(max(lit_E))

    # universal regime map: Pi = U_EO / U_cap
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    Ucap = np.logspace(-5, 0, 200)
    Ueo = np.logspace(-6, 0, 200)
    X, Y = np.meshgrid(Ucap, Ueo)
    Pi = Y / X
    cs = ax.contourf(X, Y, np.log10(Pi), levels=np.linspace(-4, 4, 33), cmap="RdBu_r", alpha=0.75, extend="both")
    cb = fig.colorbar(cs, ax=ax)
    cb.set_label("log10 Pi = log10(U_EO / U_cap)")
    for val, lab in ((0.05, "Pi = 0.05: dryout gain ~5% (negligible below)"), (1.0, "Pi = 1: EO pressure = Pc,max (dominant above)")):
        ax.loglog(Ucap, val * Ucap, "k", lw=1.6 if val == 1 else 1.0, ls="-" if val == 1 else "--")
        ax.text(Ucap[25], val * Ucap[25] * 1.4, lab, fontsize=7.5, rotation=33)
    # literature EO band
    for z, m in mob.items():
        lo, hi = m * Emin / wick.F, m * Emax / wick.F
        ax.fill_between([1e-5, 1], lo, hi, color=COLORS[z], alpha=0.12)
        ax.text(1.2e-5, hi * 0.8, f"pujahari2022 zeta = {z*1e3:.0f} mV, E = 5-100 kV/m", fontsize=7, color=COLORS[z])
    # geometry sweeps (analytical K, Maxwell F)
    gamma = liquid.surface_tension
    sweeps = {
        "pitch 15-60 um": [geometry_wick(wick.d, l, wick.h, wick.L, wick.W, theta0, gamma) for l in np.linspace(15e-6, 60e-6, 10)],
        "height 5-50 um": [geometry_wick(wick.d, wick.l, h, wick.L, wick.W, theta0, gamma) for h in np.linspace(5e-6, 50e-6, 10)],
        "feed length 1-10 mm": [replace(wick, L=L) for L in np.linspace(1e-3, 10e-3, 10)],
    }
    marks = {"pitch 15-60 um": "s", "height 5-50 um": "^", "feed length 1-10 mm": "D"}
    sweep_rows = []
    for name, wl in sweeps.items():
        uc = [w.K * w.Pc_max / (liquid.mu * w.L) for w in wl]
        y_hi = [mob[-0.2] * Emax / w.F for w in wl]
        y_lo = [mob[-0.1] * Emin / w.F for w in wl]
        ax.plot(uc, y_hi, marks[name] + "-", ms=4.5, lw=0.8, c="0.15", label=f"{name}: strongest literature EO (-200 mV, 1e5 V/m)")
        ax.plot(uc, y_lo, marks[name] + "-", ms=4.5, lw=0.8, c="0.15", mfc="white", label=f"{name}: weakest (-100 mV, 5e3 V/m)")
        for w, u in zip(wl, uc):
            row = {"sweep": name, "d_um": w.d * 1e6, "l_um": w.l * 1e6, "h_um": w.h * 1e6, "L_mm": w.L * 1e3,
                   "K_m2": w.K, "Pc_max_Pa": w.Pc_max, "F": w.F, "U_cap_m_s": u}
            for z, m in mob.items():
                for E in (Emin, 1e4, Emax):
                    P = m * E / w.F / u
                    row[f"Pi_z{abs(z*1e3):.0f}_E{E:.0e}"] = P
                    row[f"G_z{abs(z*1e3):.0f}_E{E:.0e}"] = float(w1.dryout_gain(P))
                # field for 5 % gain and for Pi = 1
                row[f"E_for_5pct_gain_z{abs(z*1e3):.0f}_V_per_m"] = brentq(lambda P: float(w1.dryout_gain(P)) - 1.05, 0, 1) * u * w.F / m
                row[f"E_for_Pi1_z{abs(z*1e3):.0f}_V_per_m"] = u * w.F / m
            sweep_rows.append(row)
    ax.plot([U_cap0, U_cap0], [mob[-0.1] * Emin / wick.F, mob[-0.2] * Emax / wick.F], "*", ms=14, c="gold", mec="k",
            label="baseline wick (d/l/h = 10/30/25 um, L = 5 mm, COMSOL K)")
    ax.axvline(U_cap0, c="gold", lw=1)
    ax.set_xlim(1e-5, 1)
    ax.set_ylim(1e-6, 1)
    ax.set_xlabel("capillary velocity scale U_cap = K Pc,max / (mu L) [m/s]   (wick + fluid)")
    ax.set_ylabel("superficial EO velocity U_EO = mu_EO E / F [m/s]   (electrokinetics)")
    ax.set_title("Falsification map: dryout gain G = (1 + Pi + sqrt(1 + 2 Pi))/2 depends only on Pi = U_EO/U_cap")
    ax.legend(loc="lower right", fontsize=6.3)
    fig.tight_layout()
    footer(fig, "Geometry sweeps use analytical K (Sangani-Acrivos + Brinkman; -1.3% vs COMSOL at baseline only) and Maxwell F; "
           "Pc,max from V2019 Eq. 3 with theta_r = 15 deg.", FOOT_EK)
    save(fig, "fig5_regime_map", bottom=0.13)
    write_csv("geometry_regimes", sweep_rows)

    # field required for 5 % gain / Pi = 1 vs geometry (both zetas)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.0), sharey=True)
    for a, (name, key, scale, xl) in zip(ax, (("pitch 15-60 um", "l_um", 1, "pitch l [um]"),
                                              ("height 5-50 um", "h_um", 1, "pillar height h [um]"),
                                              ("feed length 1-10 mm", "L_mm", 1, "feed length L [mm]"))):
        sub = [r for r in sweep_rows if r["sweep"] == name]
        x = [r[key] for r in sub]
        for z in mob:
            zz = abs(z * 1e3)
            a.semilogy(x, [r[f"E_for_5pct_gain_z{zz:.0f}_V_per_m"] for r in sub], "--", c=COLORS[z], label=f"5% gain, {z*1e3:.0f} mV")
            a.semilogy(x, [r[f"E_for_Pi1_z{zz:.0f}_V_per_m"] for r in sub], "-", c=COLORS[z], label=f"Pi = 1, {z*1e3:.0f} mV")
        a.axhspan(Emin, Emax, color="0.85", alpha=0.6, label="pujahari2022 field range")
        a.set_xlabel(xl)
        a.set_title(name)
    ax[0].set_ylabel("required field E [V/m]")
    ax[0].legend(fontsize=7)
    fig.suptitle("Field needed for EO to matter (5% dryout gain) and to dominate (Pi = 1), same wick in A and B", y=0.99)
    fig.tight_layout()
    footer(fig, "Each point is a different wick; within each point Case A and Case B share it. Analytical K/Maxwell F off-baseline.", FOOT_EK)
    save(fig, "fig6_required_field_vs_geometry", bottom=0.2)

    # one-at-a-time sensitivity (tornado) at E = Emax, zeta = -200 mV, and interval bounds
    base_m = mob[-0.2]
    def G_of(w, m=base_m, E=Emax, lq=liquid):
        return float(w1.dryout_gain(w1.eo_pressure_ratio(w1.Case("B", w, lq, w1.Forcing(m, E)))))
    def qdry(w, m=base_m, E=Emax):
        return w1.assisted_dryout_heat_flux(w1.Case("B", w, liquid, w1.Forcing(m, E))) / W_CM2
    G0 = G_of(wick)
    gam = liquid.surface_tension
    grid = np.linspace(0, 1, 25)

    def span(values):
        return min(values), max(values)

    items = [
        ("zeta: -100 / -200 mV (pujahari2022 pair)", *span([G_of(wick, m=m) for m in mob.values()])),
        ("K x 0.63-1.0 (curvature band, ASSUMED)", *span([G_of(replace(wick, K=f * wick.K)) for f in 0.63 + 0.37 * grid])),
        ("F: Maxwell vs numerical (verified)", *span([G_of(replace(wick, F=1.19131)), G0])),
        ("theta_r 0-40 deg (materials.yaml sweep)", *span([G_of(replace(wick, Pc_max=pillar_max_capillary_pressure(
            gam, wick.d, wick.l, math.radians(40 * g)))) for g in grid])),
        ("feed length 1-10 mm (fixed E)", *span([G_of(replace(wick, L=(1 + 9 * g) * 1e-3)) for g in grid])),
        ("pitch 15-60 um (analytical K)", *span([G_of(geometry_wick(wick.d, (15 + 45 * g) * 1e-6, wick.h, wick.L, wick.W, theta0, gam))
                                                 for g in grid])),
        ("height 5-50 um (analytical K)", *span([G_of(geometry_wick(wick.d, wick.l, (5 + 45 * g) * 1e-6, wick.L, wick.W, theta0, gam))
                                                 for g in grid])),
        ("E 5e3-1e5 V/m (pujahari2022 range)", *span([G_of(wick, E=Emin + (Emax - Emin) * g) for g in grid])),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    for i, (lab, a_, b_) in enumerate(items):
        lo, hi = min(a_, b_), max(a_, b_)
        ax.barh(i, hi - lo, left=lo, color="#d62728", alpha=0.6)
        ax.text(hi + 0.02, i, f"{lo:.2f} - {hi:.2f}", va="center", fontsize=7.5)
    ax.axvline(G0, c="k", lw=1)
    ax.axvline(1.0, c="0.5", lw=1, ls="--")
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([x[0] for x in items], fontsize=8)
    ax.set_xlabel("dryout gain G = q''_dry(EO) / q''_dry(passive)")
    ax.set_title(f"One-at-a-time sensitivity at zeta = -200 mV, E = 1e5 V/m (G = {G0:.2f})\nsolid: baseline; dashed: no benefit (G = 1)")
    fig.tight_layout()
    footer(fig, FOOT_BASE, "Ranges are bounding intervals, not probability distributions. K and F changes apply to BOTH cases.")
    save(fig, "fig7_sensitivity_tornado", bottom=0.17)

    # interval bound over property uncertainties at fixed baseline geometry
    corners = []
    for m in mob.values():
        for kf in (0.63, 1.0):
            for th in (0.0, math.radians(40)):
                w = replace(wick, K=kf * wick.K, Pc_max=pillar_max_capillary_pressure(liquid.surface_tension, wick.d, wick.l, th))
                corners.append({"mobility": m, "K_factor": kf, "theta_deg": math.degrees(th), "G_at_1e5": G_of(w, m=m),
                                "q_dry_B_Wcm2": qdry(w, m=m),
                                "q_dry_A_Wcm2": w1.passive_dryout_heat_flux(w1.Case("A", w, liquid, w1.Forcing())) / W_CM2})
    write_csv("sensitivity_corners", corners)
    write_csv("sensitivity_tornado", [{"item": a, "G_low_end": b, "G_high_end": c} for a, b, c in items])
    return {"G_baseline_z200_E1e5": G0, "tornado": items, "corner_G_range": [min(c["G_at_1e5"] for c in corners),
                                                                             max(c["G_at_1e5"] for c in corners)],
            "U_cap_baseline": U_cap0}


def voltage_scale(inp, liquid, wick, scen) -> dict:
    """Pi = U_EO/U_cap = V eps|zeta| / (F K Pc_max) with V = E L: the needed VOLTAGE is L-independent."""
    T_eps = liquid.relative_permittivity * ek.EPSILON_0
    theta0 = inp.si("wick.receding_contact_angle")
    gamma = liquid.surface_tension
    Pi5 = brentq(lambda P: float(w1.dryout_gain(P)) - 1.05, 0, 1)
    V_dec = 1.23  # textbook standard water-decomposition voltage at 25 C; reference line only (unsourced)

    def V_cap(w, zeta):
        return w.F * w.K * w.Pc_max / (T_eps * abs(zeta))

    out = {"Pi_for_5pct_gain": Pi5}
    for s in scen:
        out[f"V_cap_baseline_z{abs(s['zeta']*1e3):.0f}_V"] = V_cap(wick, s["zeta"])
        out[f"V_for_5pct_baseline_z{abs(s['zeta']*1e3):.0f}_V"] = Pi5 * V_cap(wick, s["zeta"])
        # check: identical to E_5% * L from the field-based route
        E5 = Pi5 * (wick.K * wick.Pc_max / (liquid.mu * wick.L)) * wick.F / s["mobility"]
        assert abs(E5 * wick.L / (Pi5 * V_cap(wick, s["zeta"])) - 1) < 1e-12

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3))
    pitches = np.linspace(15e-6, 60e-6, 30)
    heights = np.linspace(5e-6, 50e-6, 30)
    rows = []
    for s in scen:
        c = COLORS[round(s["zeta"], 3)]
        Vp = [V_cap(geometry_wick(wick.d, l, wick.h, wick.L, wick.W, theta0, gamma), s["zeta"]) for l in pitches]
        Vh = [V_cap(geometry_wick(wick.d, wick.l, h, wick.L, wick.W, theta0, gamma), s["zeta"]) for h in heights]
        for a, x, V in ((ax[0], pitches, Vp), (ax[1], heights, Vh)):
            a.semilogy(x * 1e6, Pi5 * np.array(V), "--", c=c, label=f"5% dryout gain, {s['zeta']*1e3:.0f} mV")
            a.semilogy(x * 1e6, V, "-", c=c, label=f"Pi = 1 (EO pressure = Pc,max), {s['zeta']*1e3:.0f} mV")
    for a, xl in ((ax[0], "pitch l [um] (d = 10, h = 25 um)"), (ax[1], "pillar height h [um] (d = 10, l = 30 um)")):
        a.axhline(V_dec, c="k", ls=":", label="~1.23 V water decomposition (textbook, unsourced)")
        a.set_xlabel(xl)
        a.set_ylabel("voltage across the wick V = E L [V]  (independent of L)")
        a.legend(fontsize=6.8)
    ax[0].set_title("(a) voltage needed vs pitch")
    ax[1].set_title("(b) voltage needed vs height")

    # geometric similarity: scale d, l, h by s at fixed L = 5 mm
    a = ax[2]
    sc = np.logspace(-2, 0, 60)
    qA = []
    for k, s in enumerate(scen):
        V, Gat, edl = [], [], []
        for f in sc:
            w = geometry_wick(wick.d * f, wick.l * f, wick.h * f, wick.L, wick.W, theta0, gamma)
            V.append(V_cap(w, s["zeta"]))
            Gat.append(float(w1.dryout_gain(V_dec / V[-1])))
            edl.append(s["debye_length"] / ((w.l - w.d) / 2))
            if k == 0:
                qA.append(w1.passive_dryout_heat_flux(w1.Case("A", w, liquid, w1.Forcing())) / W_CM2)
            rows.append({"scale": f, "zeta_V": s["zeta"], "l_um": w.l * 1e6, "V_cap_V": V[-1], "G_at_1p23V": Gat[-1],
                         "debye_over_halfgap": edl[-1], "q_dry_passive_Wcm2": w1.passive_dryout_heat_flux(
                             w1.Case("A", w, liquid, w1.Forcing())) / W_CM2})
        c = COLORS[round(s["zeta"], 3)]
        edl = np.array(edl)
        a.loglog(sc, Pi5 * np.array(V), "-", c=c, label=f"V for 5% gain, {s['zeta']*1e3:.0f} mV")
        bad = edl > 0.1
        if bad.any():
            a.loglog(sc[bad], Pi5 * np.array(V)[bad], "x", ms=3, c=c, label=f"lambda_D > 0.1 half-gap (thin EDL invalid)")
    a.axhline(V_dec, c="k", ls=":")
    a.set_xlabel("geometric scale factor s (d, l, h = s x 10/30/25 um; L = 5 mm)")
    a.set_ylabel("voltage for 5% dryout gain [V]")
    a2 = a.twinx()
    a2.loglog(sc, qA, c="0.4", lw=2, alpha=0.6)
    a2.set_ylabel("passive dryout q''_dry,A [W/cm$^2$] (grey)", color="0.35")
    a.set_title("(c) V_cap = F K Pc,max/(eps|zeta|) ~ s; passive q''_dry ~ s$^2$")
    a.legend(fontsize=6.5, loc="upper left")
    fig.suptitle("Falsification: EO gain is set by the VOLTAGE across the wick, Pi = V eps|zeta| / (F K Pc,max)", y=0.99)
    fig.tight_layout()
    footer(fig, "Analytical K (Sangani-Acrivos + Brinkman), Maxwell F, V2019 Eq. 3 Pc,max (theta_r = 15 deg), IAPWS water at 100 C. "
           "Electrodes assumed at the feed and at the closed end (V = E L).", FOOT_EK)
    save(fig, "fig11_voltage_scale", bottom=0.17)
    write_csv("voltage_scale_similarity", rows)
    return out


# ---------------------------------------------------------------------------------------------
# 7-8. dryout: passive and preliminary assisted
# ---------------------------------------------------------------------------------------------

def dryout_figures(inp, liquid, wick, scen) -> list[dict]:
    A = w1.Case("A", wick, liquid, w1.Forcing())
    qA = w1.passive_dryout_heat_flux(A)
    E = np.linspace(0, 1e5, 101)
    rows = []
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a = ax[0]
    a.axhline(qA / W_CM2, c=COLORS["A"], lw=2, label=f"Case A passive: {qA/W_CM2:.1f} W/cm$^2$")
    a.axhline(0.63 * qA / W_CM2, c=COLORS["A"], lw=1, ls="--", label="Case A passive if K x 0.63 (ASSUMED band)")
    for s in scen:
        c = COLORS[round(s["zeta"], 3)]
        qd, qf, qd_lo = [], [], []
        for e in E:
            B = w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], e))
            qd.append(w1.assisted_dryout_heat_flux(B) / W_CM2)
            qf.append(w1.flooding_heat_flux(B) / W_CM2)
            Blo = w1.Case("B", replace(wick, K=0.63 * wick.K), liquid, w1.Forcing(s["mobility"], e))
            qd_lo.append(w1.assisted_dryout_heat_flux(Blo) / W_CM2)
            rows.append({"zeta_V": s["zeta"], "E_V_per_m": e, "Pi": w1.eo_pressure_ratio(B),
                         "G": float(w1.dryout_gain(w1.eo_pressure_ratio(B))), "q_dry_B_Wcm2": qd[-1],
                         "q_flood_B_Wcm2": qf[-1], "q_dry_A_Wcm2": qA / W_CM2,
                         "q_dry_B_K063_Wcm2": qd_lo[-1],
                         "q_dry_A_K063_Wcm2": w1.passive_dryout_heat_flux(replace(A, wick=Blo.wick)) / W_CM2})
        a.plot(E / 1e3, qd, c=c, lw=2, label=f"Case B dryout, {s['zeta']*1e3:.0f} mV")
        a.fill_between(E / 1e3, qf, qd, color=c, alpha=0.12)
        a.plot(E / 1e3, qf, ":", c=c, label=f"Case B flooding limit Pc(L) = 0, {s['zeta']*1e3:.0f} mV")
        a.plot(E / 1e3, qd_lo, "--", c=c, lw=1, label="Case B dryout if K x 0.63 (ASSUMED band)")
    a.set_xlabel("applied field E [kV/m]")
    a.set_ylabel("heat flux q'' [W/cm$^2$]")
    a.set_title("(a) PRELIMINARY: operating window between flooding and dryout")
    a.legend(fontsize=6.8, loc="upper left")
    a = ax[1]
    s = scen[1]
    xx = np.linspace(0, wick.L, 300)
    for q, ls in ((80e4, "-"), (40e4, "--")):
        a.plot(xx * 1e3, w1.capillary_pressure_profile(A, q, xx) / wick.Pc_max, ls, c=COLORS["A"], label=f"A, {q/W_CM2:.0f} W/cm$^2$")
        for e in (1e4, 1e5):
            B = w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], e))
            a.plot(xx * 1e3, w1.capillary_pressure_profile(B, q, xx) / wick.Pc_max, ls, c=COLORS[-0.2],
                   alpha=0.45 if e == 1e4 else 1.0, label=f"B, {q/W_CM2:.0f} W/cm$^2$, E = {e:.0e}")
    a.axhline(1, c="k", lw=1)
    a.axhline(0, c="0.4", lw=0.8)
    a.text(3.2, 1.02, "dryout: Pc = Pc,max", fontsize=7.5)
    a.text(0.1, -0.08, "below 0: liquid pressure > vapor pressure (convex meniscus / flooding)", fontsize=7)
    a.set_xlabel("distance from feed x [mm]")
    a.set_ylabel("Pc(x) / Pc,max")
    a.set_title("(b) mechanism: EO lowers the required Pc (-200 mV)")
    a.legend(fontsize=6.8, loc="upper left")
    fig.tight_layout()
    footer(fig, FOOT_BASE, FOOT_EK,
           "Valid only under: thin EDL, uniform zeta on pillars AND floor, uniform E, flat meniscus K, no salt accumulation, no electrolysis gas, "
           "no field-effect on zeta (see validity diagnostics). NOT a device prediction.")
    save(fig, "fig8_preliminary_dryout", bottom=0.2)
    return rows


# ---------------------------------------------------------------------------------------------
# cross-check against pujahari2022's reported gains
# ---------------------------------------------------------------------------------------------

def pujahari_cross_check(inp, quick: bool) -> list[dict]:
    b = f"{SC}.property_set"
    eps_r = inp.si(f"{SC}.relative_permittivity")
    liq = w1.Liquid(mu=inp.si(f"{b}.mu"), rho=inp.si(f"{b}.rho"), h_fg=inp.si(f"{b}.h_fg"),
                    surface_tension=inp.si(f"{b}.surface_tension"), relative_permittivity=eps_r)
    cc = f"{SC}.cross_check"
    d, L, th = inp.si(f"{cc}.pillar_diameter"), inp.si(f"{cc}.feed_length"), inp.si(f"{cc}.receding_contact_angle")
    E = inp.si(f"{cc}.field_for_reported_gain")
    zetas = list(inp.si(f"{SC}.zeta_potential"))
    rows = []
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    levels = (64, 128, 256) if quick else (128, 256, 512)
    for a, key in zip(ax, ("tall", "short")):
        h = inp.si(f"{cc}.{key}.pillar_height")
        lo_l, hi_l = inp.si(f"{cc}.{key}.pitch_range")
        pitches = np.linspace(lo_l, hi_l, 8)
        for z in zetas:
            m = ek.electroosmotic_mobility(eps_r, z, liq.mu)
            rep = inp.si(f"{cc}.{key}.reported_gain_zeta{abs(z*1e3):.0f}")
            Gm, Gn = [], []
            for l in pitches:
                Fnum = _F_numerical(float(d), float(l), levels)
                wm = w1.Wick.analytical(d, l, h, L, 1e-3, th, liq.surface_tension)
                wn = replace(wm, F=Fnum)
                Gm.append(float(w1.dryout_gain(w1.eo_pressure_ratio(w1.Case("B", wm, liq, w1.Forcing(m, E))))))
                Gn.append(float(w1.dryout_gain(w1.eo_pressure_ratio(w1.Case("B", wn, liq, w1.Forcing(m, E))))))
                lam = inp.si(f"{SC}.debye_length")[zetas.index(z)]
                rows.append({"set": key, "h_um": h * 1e6, "l_um": l * 1e6, "zeta_V": z, "phi": solid_fraction(d, l),
                             "K_analytical_m2": wm.K, "Pc_max_Pa": wm.Pc_max, "F_maxwell": wm.F, "F_numerical": Fnum,
                             "Pi": w1.eo_pressure_ratio(w1.Case("B", wn, liq, w1.Forcing(m, E))), "G_model": Gn[-1],
                             "reported_lo": rep[0], "reported_hi": rep[1],
                             "debye_over_halfgap": lam / ((l - d) / 2), "debye_over_h": lam / h,
                             "meniscus_sag_over_h_at_Pcmax": _sag(liq.surface_tension, wm.Pc_max, d, l) / h})
            c = COLORS[round(z, 3)]
            a.plot(pitches * 1e6, Gn, "o-", c=c, label=f"this model, {z*1e3:.0f} mV (numerical F)")
            a.plot(pitches * 1e6, Gm, ":", c=c, lw=1, label="this model, Maxwell F")
            a.axhspan(rep[0], rep[1], color=c, alpha=0.18, label=f"pujahari2022 stated range, {z*1e3:.0f} mV")
        a.set_title(f"d = 5 um, h = {h*1e6:.0f} um, L = 5 mm, E = 1e5 V/m")
        a.set_xlabel("pitch l [um]")
        a.set_ylabel("dryout gain G")
        a.legend(fontsize=7)
    fig.suptitle("Cross-check (model vs model, NOT validation): reduced-order HS + Darcy vs pujahari2022 Poisson-Boltzmann + NS", y=0.99)
    fig.tight_layout()
    footer(fig, "pujahari2022 property set (Table 1) and stated gain ranges (Sec. 3.4); their ranges span the pitch range and are not tied to individual pitches.",
           "This model: flat meniscus, analytical K (Sangani-Acrivos + Brinkman), uniform evaporation. Their model: curved menisci, full EDL body force.")
    save(fig, "fig9_pujahari_cross_check", bottom=0.17)
    write_csv("pujahari_cross_check", rows)
    return rows


@lru_cache(maxsize=None)
def _F_numerical(d, l, levels):
    return formation_factor_numerical(d, l, levels).F_extrapolated


def _sag(gamma, Pc, d, l):
    """Spherical-cap sag over the half-diagonal gap at Pc (rough meniscus-depth indicator)."""
    R = 2 * gamma / Pc
    half = (l * math.sqrt(2) - d) / 2
    return R - math.sqrt(max(R**2 - half**2, 0.0))


# ---------------------------------------------------------------------------------------------
# validity diagnostics
# ---------------------------------------------------------------------------------------------

def diagnostics(inp, liquid, wick, scen) -> list[dict]:
    T = inp.si("operating.vapor_temperature")
    sigma = inp.si(f"{SC}.electrical_conductivity")
    spacing = inp.si(f"{SC}.electrode_spacing")
    eps_ox = inp.si("exploration_ranges.oxide_relative_permittivity")
    t_ox_lo, t_ox_hi = inp.si("exploration_ranges.oxide_thickness")
    halfgap = (wick.l - wick.d) / 2
    rows = []
    for s in scen:
        lam, z, c = s["debye_length"], s["zeta"], s["concentration"]
        sig_native = ek.grahame_surface_charge(liquid.relative_permittivity, T, c, z)
        for E in inp.si(f"{SC}.electric_field"):
            B = w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], E))
            q = w1.assisted_dryout_heat_flux(B)
            U0 = float(w1.feed_velocity(B, q))
            gradp0 = liquid.mu / wick.K * (U0 - w1.eo_superficial_velocity(B))
            V = E * spacing
            re = w1.reynolds_numbers(B, U0)
            rows.append({
                "zeta_mV": z * 1e3, "conc_M": c / 1e3, "E_V_per_m": E, "voltage_over_spacing_V": V,
                "q_dry_B_Wcm2": q / W_CM2, "U_feed_at_dryout_m_s": U0,
                "Re_K_feed": float(re["Re_K"]), "Re_gap_feed": float(re["Re_gap"]),
                "debye_over_halfgap": lam / halfgap, "debye_over_h": lam / wick.h,
                "thin_edl_flux_factor_halfgap": ek.thin_edl_flux_factor(lam, z, T, halfgap),
                "dukhin_lower_bound_halfgap": ek.dukhin_number_lower_bound(lam, z, T, halfgap),
                "joule_over_qdry": float(ek.joule_heat_flux(sigma, E, wick.h, wick.F)) / q,
                "sigma_for_1pct_joule_uS_cm": 0.01 * q * wick.F / (E**2 * wick.h) / 1e-4,
                "streaming_over_conduction_at_feed": float(ek.streaming_to_conduction_ratio(s["mobility"], abs(gradp0), sigma, E)),
                "electroviscous_beta": ek.electroviscous_coefficient(s["mobility"], liquid.mu, wick.K, sigma, wick.F),
                "native_charge_C_m2": sig_native,
                "field_effect_ratio_tox_hi_Vhalf": ek.field_effect_charge_ratio(eps_ox, t_ox_hi, V / 2, sig_native),
                "field_effect_ratio_tox_lo_Vhalf": ek.field_effect_charge_ratio(eps_ox, t_ox_lo, V / 2, sig_native),
                "oxide_field_tox_hi_V_per_um": V / 2 / t_ox_hi / 1e6,
                "salt_Peclet_D_1e-7": U0 * wick.L / 1e-7, "salt_Peclet_D_1e-9": U0 * wick.L / 1e-9,
                # sealed system, zero net ion flux: c U = D dc/dx -> c(L)/c(0) = exp(U0 L / 2D) (migration neglected)
                "log10_closed_end_enrichment_D_1e-7": U0 * wick.L / (2 * 1e-7) / math.log(10),
            })
    write_csv("validity_diagnostics", rows)

    fig, ax = plt.subplots(2, 3, figsize=(12.5, 6.4))
    E = np.logspace(3, 5, 60)
    for s in scen:
        c = COLORS[round(s["zeta"], 3)]
        lam, z, conc = s["debye_length"], s["zeta"], s["concentration"]
        sig_native = ek.grahame_surface_charge(liquid.relative_permittivity, T, conc, z)
        qd = np.array([w1.assisted_dryout_heat_flux(w1.Case("B", wick, liquid, w1.Forcing(s["mobility"], e))) for e in E])
        U0 = qd * wick.L / (liquid.rho * liquid.h_fg * wick.h)
        lab = f"{z*1e3:.0f} mV"
        ax[0, 0].loglog(E, liquid.rho * U0 * math.sqrt(wick.K) / liquid.mu, c=c, label=f"Re_K at feed, {lab}")
        ax[0, 1].semilogx(E, [ek.dukhin_number_lower_bound(lam, z, T, halfgap)] * len(E), c=c, label=f"Du lower bound, {lab}")
        ax[0, 1].semilogx(E, [1 - ek.thin_edl_flux_factor(lam, z, T, halfgap)] * len(E), "--", c=c, label=f"finite-EDL flux error, {lab}")
        ax[0, 2].loglog(E, ek.joule_heat_flux(sigma, E, wick.h, wick.F) / qd, c=c, label=f"Joule / q''_dry (0.79 uS/cm), {lab}")
        gp = liquid.mu / wick.K * np.abs(U0 - s["mobility"] * E / wick.F)
        ax[1, 0].loglog(E, ek.streaming_to_conduction_ratio(s["mobility"], gp, sigma, E), c=c, label=f"i_stream/i_cond at feed, {lab}")
        for t, ls in ((t_ox_hi, "-"), (t_ox_lo, "--")):
            ax[1, 1].loglog(E, [ek.field_effect_charge_ratio(eps_ox, t, e * spacing / 2, sig_native) for e in E], ls, c=c,
                            label=f"t_ox = {t*1e6:.1f} um, {lab}")
    ax[0, 0].axhline(1, c="k", lw=1)
    ax[0, 0].set_title("inertia: Darcy-Forchheimer onset Re_K ~ 1-10")
    ax[0, 1].axhline(0.1, c="k", lw=1)
    ax[0, 1].set_title("thin-EDL / surface conduction (need << 1)")
    ax[0, 2].set_title("Joule heating (negligible)")
    ax[1, 0].axhline(0.1, c="k", lw=1)
    ax[1, 0].set_title("streaming-current coupling (need << 1)")
    ax[1, 1].axhline(1, c="k", lw=1)
    ax[1, 1].set_title("field-effect charge / native charge (need << 1)")
    a = ax[1, 2]
    a.loglog(E, E * spacing, c="k", label="V = E x 5 mm")
    a.axhline(1.23, c="tab:red", ls="--", label="~1.2 V water decomposition (textbook, unsourced)")
    a.set_title("applied voltage vs electrolysis")
    for a in ax.ravel():
        a.set_xlabel("E [V/m]")
        a.legend(fontsize=6.5)
    fig.suptitle("Validity diagnostics of the reduced-order EO model at the predicted Case B dryout point", y=0.995)
    fig.tight_layout()
    footer(fig, FOOT_EK, "Oxide eps_r 3.9 ASSUMED; t_ox swept 0.1-2 um (REQUIRES_SOURCE); substrate assumed floating at mid potential (V/2 across oxide). "
           "Du uses the m = 0 Bikerman lower bound (cited from memory).")
    save(fig, "fig10_validity_diagnostics", bottom=0.1)
    return rows


# ---------------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    inp = load_eo_inputs()
    liquid = baseline_liquid(inp)
    wick = baseline_wick(inp, liquid)
    scen = scenario_mobilities(inp, liquid)
    with open(OUT / "parameter_provenance.csv", "w", newline="", encoding="utf-8") as fh:
        rows = inp.provenance_rows()
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)

    summary = {"wick": wick.__dict__, "liquid": liquid.__dict__,
               "scenarios": [{k: float(v) if not isinstance(v, str) else v for k, v in s.items()} for s in scen]}
    print("[1] passive benchmark ...")
    summary["passive_benchmark"] = passive_benchmark(inp, liquid, wick, args.quick)
    print("[2-5] transport comparisons ...")
    write_csv("forced_flow_comparison", transport_figures(inp, liquid, wick, scen))
    print("[6] sensitivity and falsification ...")
    fals = falsification(inp, liquid, wick, scen)
    summary["falsification"] = {k: v for k, v in fals.items() if k != "tornado"}
    summary["tornado"] = [{"item": a, "lo": float(min(b, c)), "hi": float(max(b, c))} for a, b, c in fals["tornado"]]
    summary["voltage_scale"] = voltage_scale(inp, liquid, wick, scen)
    print("[7-8] dryout ...")
    dr = dryout_figures(inp, liquid, wick, scen)
    write_csv("dryout_vs_field", dr)
    summary["dryout_at_literature_fields"] = [r for r in dr if r["E_V_per_m"] in (0.0, 5e3, 1e4, 1e5)]
    print("[x] pujahari2022 cross-check ...")
    cc = pujahari_cross_check(inp, args.quick)
    summary["pujahari_cross_check"] = pd.DataFrame(cc).groupby(["set", "zeta_V"]).agg(
        G_min=("G_model", "min"), G_max=("G_model", "max"), rep_lo=("reported_lo", "first"),
        rep_hi=("reported_hi", "first"), F_rel_max=("F_numerical", "max")).reset_index().to_dict("records")
    print("[v] validity diagnostics ...")
    summary["diagnostics"] = diagnostics(inp, liquid, wick, scen)
    summary["requires_source"] = [m.__dict__ for m in inp.missing.values()]
    with open(OUT / "summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print(json.dumps({k: summary[k] for k in ("passive_benchmark", "falsification", "pujahari_cross_check")}, indent=1, default=float))
    print(f"\nWrote figures and tables to {OUT}")


if __name__ == "__main__":
    main()
