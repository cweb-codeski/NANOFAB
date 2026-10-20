"""Independent test of the field-treatment hypothesis for the pujahari2022 EO-flux discrepancy.

Case A: uniform imposed field E0 x. Case B: field solved around insulating pillars (same mean field).
Same HS slip, zeta, mu, geometry. 2D planform unit cell (src/nanofab/electroosmosis/planform2d.py):
pillar height does not enter; height dependence is shown only as a labelled heuristic bridge between
two exact limits (h -> infinity: 2D result; h -> 0: ratio -> 1).

Usage:
    python analysis/eo_field_treatment_test.py      -> results/electroosmosis/field_treatment/
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from _common import REPO  # noqa: F401  (sets sys.path)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nanofab.electroosmosis import planform2d as pf  # noqa: E402
from nanofab.electroosmosis.formation_factor import formation_factor_maxwell_2d, formation_factor_numerical  # noqa: E402
from nanofab.geometry import solid_fraction  # noqa: E402
from nanofab.permeability import sangani_acrivos_permeability  # noqa: E402

OUT = REPO / "results" / "electroosmosis" / "field_treatment"
REPRO = REPO / "results" / "electroosmosis" / "pujahari_reproduction.csv"
N_MODES = 24
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3, "savefig.dpi": 200})

GEOMETRIES = {
    "pujahari_h25": {"d": 5e-6, "h": 25e-6, "pitches": [10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45]},
    "pujahari_h5": {"d": 5e-6, "h": 5e-6, "pitches": [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]},
    "nanofab_baseline": {"d": 10e-6, "h": 25e-6, "pitches": [30]},
    "user_mask_r10_h40_l50": {"d": 20e-6, "h": 40e-6, "pitches": [50]},
}


def bridge(ratio_2d, q_uniform, d, l, h):
    """Heuristic height bridge for Q_uniform/Q_solved (labelled; NOT a 3D solution).

    Case B is exactly z-invariant (similitude), so only Case A changes with h. The floor slip of
    Case A reaches a depth fraction f = tanh(h/sqrt(K_s))/(h/sqrt(K_s)) (Brinkman). Two closures:
      exact-endpoint:  ratio = R2D + f (1 - R2D)                         (-> 1 as h -> 0, exact)
      Brinkman shear:  ratio = F [q_A + f (1 - phi - q_A)]               (uniform floor shear)
    """
    bh = h / math.sqrt(sangani_acrivos_permeability(d, l))
    f = math.tanh(bh) / bh
    phi = solid_fraction(d, l)
    F = 1.0 / q_solved_cache[(d, l)]
    return f, ratio_2d + f * (1 - ratio_2d), F * (q_uniform + f * (1 - phi - q_uniform))


q_solved_cache: dict = {}


def run_geometries(repro: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in GEOMETRIES.items():
        for l_um in g["pitches"]:
            d, l, h = g["d"], l_um * 1e-6, g["h"]
            r = pf.field_treatment(d, l, N_MODES)
            q_solved_cache[(d, l)] = r.q_solved
            f, b_exact, b_brink = bridge(r.uniform_over_solved, r.q_uniform, d, l, h)
            R_inf = np.nan
            if name.startswith("pujahari"):
                h_um = g["h"] * 1e6
                sel = repro[(repro.h_um == h_um) & (np.isclose(repro.l_um, l_um))]
                if not sel.empty:
                    R_inf = float(sel.R_inferred.mean())
            rows.append({
                "geometry": name, "d_um": d * 1e6, "l_um": l_um, "h_um": h * 1e6, "d_over_l": d / l,
                "phi": solid_fraction(d, l), "F_solved": 1 / r.inv_F, "F_maxwell": formation_factor_maxwell_2d(solid_fraction(d, l)),
                "q_solved": r.q_solved, "q_uniform": r.q_uniform,
                "R_field_solved_over_uniform": r.R_field, "uniform_over_solved_2D": r.uniform_over_solved,
                "floor_fraction_f": f, "uniform_over_solved_bridge_exact": b_exact,
                "uniform_over_solved_bridge_brinkman": b_brink,
                "R_inferred_pujahari": R_inf,
                "diff_2D_minus_inferred": r.uniform_over_solved - R_inf,
                "max_bc_residual": max(max(r.flow_solved.residuals.values()), max(r.flow_uniform.residuals.values()),
                                       r.field.residual_equipotential, r.field.residual_mirror),
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "R_field_by_geometry.csv", index=False)
    return df


def verification() -> dict:
    """Checks requested before interpreting R (baseline cell and the densest/sparsest Pujahari cells)."""
    out = {"convergence": [], "checks": {}}
    for d, l in ((10e-6, 30e-6), (5e-6, 8e-6), (5e-6, 45e-6), (20e-6, 50e-6)):
        for N in (8, 12, 16, 24, 32):
            r = pf.field_treatment(d, l, N)
            out["convergence"].append({
                "d_um": d * 1e6, "l_um": l * 1e6, "n_modes": N, "collocation_points_per_edge": 4 * N,
                "inv_F": r.inv_F, "q_solved": r.q_solved, "q_uniform": r.q_uniform,
                "uniform_over_solved": r.uniform_over_solved,
                "max_residual_flow": max(max(r.flow_solved.residuals.values()), max(r.flow_uniform.residuals.values())),
                "max_residual_field": max(r.field.residual_equipotential, r.field.residual_mirror)})
        # collocation-density independence at fixed N
        for m in (48, 96, 192):
            fs = pf.solve_field(d / l, N_MODES, m)
            fa = pf.solve_flow(d / l, N_MODES, slip_modes=[-1.0], m_points=m)
            out["convergence"].append({"d_um": d * 1e6, "l_um": l * 1e6, "n_modes": N_MODES,
                                       "collocation_points_per_edge": m, "inv_F": fs.inverse_formation_factor,
                                       "q_uniform": fa.q})
    pd.DataFrame(out["convergence"]).to_csv(OUT / "convergence.csv", index=False)

    d, l = 10e-6, 30e-6
    r = pf.field_treatment(d, l, N_MODES)
    a = d / l
    sections = [1.0, 0.7, 0.5 * a, 0.0]
    currents = [r.field.current(x0) for x0 in sections]
    fluxes_b = [r.flow_solved.section_flux(x0) for x0 in sections]
    fluxes_a = [r.flow_uniform.section_flux(x0) for x0 in sections]
    rev = pf.field_treatment(d, l, N_MODES, E_sign=-1.0)
    zero = pf.solve_flow(a, N_MODES, slip_modes=np.zeros(N_MODES))
    fv = formation_factor_numerical(d, l, (128, 256, 512)).F_extrapolated
    s = np.linspace(-1, 1, 401)
    phi_right = r.field.phi(np.ones_like(s), s)
    phi_left = r.field.phi(-np.ones_like(s), s)
    out["checks"] = {
        "current_conservation_rel_spread": float(np.ptp(currents) / np.mean(currents)),
        "currents_at_sections": currents,
        "F_vs_independent_FV_rel": float((1 / r.inv_F - fv) / fv),
        "mean_potential_drop_over_cell_E0_units": float(np.mean(phi_left - phi_right)),  # expect 2 (= E0 l / (E0 l/2))
        "periodicity_phi_plus_E0x_left_minus_right_max": float(np.max(np.abs((phi_left - 1) - (phi_right + 1)))),
        "mass_conservation_rel_spread_B": float(np.ptp(fluxes_b) / abs(np.mean(fluxes_b))),
        "mass_conservation_rel_spread_A": float(np.ptp(fluxes_a) / abs(np.mean(fluxes_a))),
        "similitude_q_solved_minus_invF": float(r.q_solved - r.inv_F),
        "reversal_B_sum": float(r.q_solved + rev.q_solved), "reversal_A_sum": float(r.q_uniform + rev.q_uniform),
        "zero_field_q": float(zero.q),
        "uniform_slip_residual": r.flow_uniform.residuals["wall_slip"],
        "solved_slip_residual": r.flow_solved.residuals["wall_slip"],
        "K2D_over_sangani_acrivos_baseline": float(pf.permeability_2d(d, l, N_MODES) / sangani_acrivos_permeability(d, l)),
        "vanishing_obstruction": [{"d_over_l": x, "inv_F": pf.field_treatment(x * 30e-6, 30e-6, N_MODES).inv_F,
                                   "uniform_over_solved": pf.field_treatment(x * 30e-6, 30e-6, N_MODES).uniform_over_solved}
                                  for x in (0.01, 0.03, 0.1)],
    }
    return out


def figures(df: pd.DataFrame):
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for a, name, title in ((ax[0], "pujahari_h25", "d = 5 um, h = 25 um"), (ax[1], "pujahari_h5", "d = 5 um, h = 5 um")):
        g = df[df.geometry == name]
        a.plot(g.l_um, g.uniform_over_solved_2D, "k-", lw=2, label="2D planform (exact for h -> infinity)")
        a.fill_between(g.l_um, g.uniform_over_solved_bridge_exact, g.uniform_over_solved_bridge_brinkman, color="tab:green",
                       alpha=0.25, label="height bridge, heuristic band (two closures)")
        a.plot(g.l_um, g.R_inferred_pujahari, "o", c="tab:red", label="inferred from pujahari2022 Fig. 6")
        a.axhline(1.0, c="0.5", ls=":", lw=1)
        a.set_ylim(0.3, 1.15)
        a.set_xlabel("pitch l [um]")
        a.set_ylabel("Q_EO(uniform field) / Q_EO(solved field) = 1 / R_field")
        a.set_title(title)
        sec = a.secondary_yaxis("right", functions=(lambda v: 1 / np.clip(v, 1e-3, None), lambda v: 1 / np.clip(v, 1e-3, None)))
        sec.set_ylabel("R_field = Q_solved / Q_uniform")
        a.legend(fontsize=7.5, loc="upper left")
    for name, mk, lab in (("nanofab_baseline", "*", "NANOFAB baseline d/l/h = 10/30/25"),
                          ("user_mask_r10_h40_l50", "D", "mask r = 10, h = 40, l = 50 um")):
        g = df[df.geometry == name].iloc[0]
        ax[0].plot([g.l_um], [g.uniform_over_solved_2D], mk, ms=11, c="tab:blue", mec="k", label=f"{lab} (2D)")
        ax[0].errorbar([g.l_um], [(g.uniform_over_solved_bridge_exact + g.uniform_over_solved_bridge_brinkman) / 2],
                       yerr=[[abs(g.uniform_over_solved_bridge_exact - g.uniform_over_solved_bridge_brinkman) / 2]],
                       fmt="none", ecolor="tab:blue", capsize=4)
    ax[0].legend(fontsize=7, loc="upper left")
    fig.suptitle("Field-treatment test: uniform imposed field vs field solved around insulating pillars (same HS slip, zeta, mean E)", y=1.0)
    fig.tight_layout()
    fig.text(0.01, 0.0, "2D planform: z-invariant, no floor, no meniscus; height enters only through the labelled heuristic bridge. "
             "Tall-pillar values do not depend on h.", fontsize=7, color="0.35")
    fig.subplots_adjust(bottom=0.14)
    fig.savefig(OUT / "fig_R_vs_pitch.png")
    plt.close(fig)

    # contours for the baseline cell
    d, l = 10e-6, 30e-6
    r = pf.field_treatment(d, l, N_MODES)
    a = d / l
    x = np.linspace(-1, 1, 321)
    X, Y = np.meshgrid(x, x)
    solid = X**2 + Y**2 < a**2
    phi = np.where(solid, np.nan, r.field.phi(X, Y))
    ex, ey = r.field.grad(X, Y)
    Emag = np.where(solid, np.nan, np.hypot(ex, ey))
    uB, vB = r.flow_solved.velocity(X, Y)
    uA, vA = r.flow_uniform.velocity(X, Y)
    uB, uA = np.where(solid, np.nan, uB), np.where(solid, np.nan, uA)
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
    c0 = ax[0].contourf(X, Y, Emag, 30, cmap="viridis")
    ax[0].contour(X, Y, phi, 21, colors="w", linewidths=0.6)
    fig.colorbar(c0, ax=ax[0], label="|E| / E0 (solved field)")
    ax[0].set_title("(a) solved |E|/E0, equipotentials", fontsize=9)
    th = np.linspace(0, np.pi, 181)
    e_modes = r.field.wall_tangential_field_modes()
    Et = sum(e_modes[k] * np.sin(n * th) for k, n in enumerate(r.field.n))
    ax[1].plot(np.degrees(th), -Et, label="solved field: -E_theta/E0 on pillar wall")
    ax[1].plot(np.degrees(th), np.sin(th), "--", label="uniform field: sin(theta)")
    ax[1].plot(np.degrees(th), 2 * np.sin(th), ":", c="0.4", label="isolated cylinder: 2 sin(theta)")
    ax[1].set_xlabel("theta [deg] (0 = downstream)")
    ax[1].set_ylabel("tangential field magnitude / E0")
    ax[1].set_title("(b) wall-tangential field = HS slip/(mu_EO E0)", fontsize=9)
    ax[1].legend(fontsize=7.5)
    lev = np.linspace(0, max(np.nanmax(uB), np.nanmax(uA)), 31)
    for k, (u, lab) in enumerate(((uA, f"(c) A, uniform field: u/(mu_EO E0)\nsuperficial = {r.q_uniform:.3f}"),
                                  (uB, f"(d) B, solved field: u/(mu_EO E0)\nsuperficial = {r.q_solved:.3f} = 1/F"))):
        c = ax[2 + k].contourf(X, Y, u, lev, cmap="magma", extend="both")
        fig.colorbar(c, ax=ax[2 + k])
        ax[2 + k].set_title(lab, fontsize=9)
    for k in (0, 2, 3):
        ax[k].set_aspect("equal")
        ax[k].add_patch(plt.Circle((0, 0), a, color="0.6"))
        ax[k].set_xlabel("x / (l/2)")
        ax[k].set_ylabel("y / (l/2)")
    fig.suptitle("NANOFAB baseline cell d/l = 1/3 (2D planform, tall-pillar limit)", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "fig_baseline_cell_fields.png")
    plt.close(fig)


def write_tables(df: pd.DataFrame, ver: dict):
    lines = ["| geometry | d/l/h [um] | 1/R_field = Q_uniform/Q_solved (2D) | height bridge (heuristic) | R inferred (pujahari2022) | 2D minus inferred |",
             "| --- | --- | --- | --- | --- | --- |"]
    for _, g in df.iterrows():
        inf = "-" if not np.isfinite(g.R_inferred_pujahari) else f"{g.R_inferred_pujahari:.3f}"
        dif = "-" if not np.isfinite(g.diff_2D_minus_inferred) else f"{g.diff_2D_minus_inferred:+.3f}"
        lo, hi = sorted([g.uniform_over_solved_bridge_exact, g.uniform_over_solved_bridge_brinkman])
        lines.append(f"| {g.geometry} | {g.d_um:g}/{g.l_um:g}/{g.h_um:g} | {g.uniform_over_solved_2D:.3f} (R_field {g.R_field_solved_over_uniform:.2f}) "
                     f"| {lo:.2f}-{hi:.2f} | {inf} | {dif} |")
    conv = pd.DataFrame(ver["convergence"])
    cl = ["| d/l | modes N | points/edge | 1/F | Q_uniform/(mu E0 A) | Q_solved/(mu E0 A) | ratio | max BC residual (flow) |",
          "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for _, c in conv.iterrows():
        def fmt(v, f="{:.7f}"):
            return f.format(v) if pd.notna(v) else "-"
        cl.append(f"| {c.d_um / c.l_um:.3f} | {c.n_modes:g} | {c.collocation_points_per_edge:g} | {fmt(c.inv_F)} | {fmt(c.q_uniform)} | "
                  f"{fmt(c.get('q_solved'))} | {fmt(c.get('uniform_over_solved'), '{:.6f}')} | {fmt(c.get('max_residual_flow'), '{:.1e}')} |")
    (OUT / "tables.md").write_text("## Comparison\n\n" + "\n".join(lines) + "\n\n## Convergence\n\n" + "\n".join(cl) + "\n",
                                   encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    repro = pd.read_csv(REPRO)
    df = run_geometries(repro)
    ver = verification()
    figures(df)
    write_tables(df, ver)
    pd.set_option("display.width", 220)
    print(df[["geometry", "l_um", "d_over_l", "uniform_over_solved_2D", "R_field_solved_over_uniform", "floor_fraction_f",
              "uniform_over_solved_bridge_exact", "uniform_over_solved_bridge_brinkman", "R_inferred_pujahari",
              "max_bc_residual"]].round(4).to_string(index=False))
    pj = df[df.geometry.str.startswith("pujahari")].dropna(subset=["R_inferred_pujahari"])
    lo = pj.uniform_over_solved_2D
    hi = pj[["uniform_over_solved_bridge_exact", "uniform_over_solved_bridge_brinkman"]].max(axis=1)
    inside = (pj.R_inferred_pujahari >= lo) & (pj.R_inferred_pujahari <= hi)
    ver["checks"]["inferred_inside_2D_to_bridge_band"] = f"{int(inside.sum())} of {len(pj)}"
    ver["checks"]["outside_band"] = pj.loc[~inside, ["geometry", "l_um", "uniform_over_solved_2D", "R_inferred_pujahari"]].round(3).to_dict("records")
    for k, v in ver["checks"].items():
        print(f"{k}: {v}")
    pd.Series({k: str(v) for k, v in ver["checks"].items()}).to_csv(OUT / "verification_checks.csv", header=["value"])
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
