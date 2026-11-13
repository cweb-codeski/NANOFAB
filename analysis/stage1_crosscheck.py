"""Compare two independent Stage 1 (flat-interface permeability) builds of the reference cell.

Inputs (reference geometry d/l/hp = 10/30/25 um, quarter cell, dp = 1 Pa):
    data/processed/permeability_stage1_flow.csv            Chase, permeability1.mph (stored derived result)
    data/reported/stage1_student_reference_d10_l30_h25.csv  undergraduate, REPORTED_UNVERIFIED

Reported quantities:
    relative change U_i = (K_i - K_{i-1}) / K_{i-1}  (coarser value in the denominator, as reported
                                                     by the student and as in V2019 Eq. A1)
    builder difference at matching hmax            (K_student - K_chase) / K_chase
    observed order p and Richardson estimate K_inf from the three finest levels of each build

Free tetrahedral meshes controlled only by hmax are not geometrically similar refinements, so p and
K_inf are indicative, not formal error bounds (COMSOL KB 1261 discusses refinement strategies).

Usage:
    python analysis/stage1_crosscheck.py
"""

from __future__ import annotations

import math

import pandas as pd
from scipy.optimize import brentq

from _common import REPO
from nanofab.permeability import pillar_wick_permeability

CHASE = REPO / "data" / "processed" / "permeability_stage1_flow.csv"
STUDENT = REPO / "data" / "reported" / "stage1_student_reference_d10_l30_h25.csv"
D, L, HP = 10e-6, 30e-6, 25e-6


def load() -> tuple[pd.Series, pd.Series, pd.Series]:
    chase = pd.read_csv(CHASE)
    chase_K = pd.Series(chase["K_comsol_m2"].to_numpy(), index=(chase["hmax"] * 1e6).round(6))
    stu = pd.read_csv(STUDENT, comment="#")
    stu_K = pd.Series(stu["K_m2"].to_numpy(), index=stu["hmax_um"].astype(float))
    stu_Q = pd.Series(stu["Q_m3_s"].to_numpy(), index=stu["hmax_um"].astype(float))
    return chase_K.sort_index(ascending=False), stu_K.sort_index(ascending=False), stu_Q.sort_index(ascending=False)


def relative_changes(K: pd.Series) -> list[tuple[float, float]]:
    """(hmax, (K_i - K_{i-1})/K_{i-1}) for successive levels, coarse to fine."""
    vals = list(K.items())
    return [(h1, (k1 - k0) / k0) for (_, k0), (h1, k1) in zip(vals, vals[1:])]


def observed_order(h: tuple[float, float, float], K: tuple[float, float, float]) -> tuple[float, float]:
    """(p, K_inf) for K(h) = K_inf + C h^p through three points (h coarse -> fine)."""
    (h1, h2, h3), (k1, k2, k3) = h, K
    target = (k2 - k1) / (k3 - k2)

    def f(p: float) -> float:
        return (h2**p - h1**p) / (h3**p - h2**p) - target

    p = brentq(f, 0.05, 8.0)
    C = (k3 - k2) / (h3**p - h2**p)
    return p, k3 - C * h3**p


def main() -> None:
    chase, stu, stu_Q = load()
    K_an = pillar_wick_permeability(D, L, HP)

    print("Stage 1 reference cell (d/l/hp = 10/30/25 um), flat shear-free top, dp = 1 Pa, quarter cell\n")
    for name, K in (("Chase (permeability1.mph, stored)", chase), ("Student (REPORTED_UNVERIFIED)", stu)):
        print(f"{name}")
        print("  hmax_um   K_m2          rel. change vs coarser   K/K_analytic - 1")
        changes = dict(relative_changes(K))
        for h, k in K.items():
            rc = f"{changes[h]:+.4%}" if h in changes else "   -"
            print(f"  {h:<8g}  {k:.6e}   {rc:>22}   {k / K_an - 1:+.3%}")
        hs = tuple(K.index[-3:])
        ks = tuple(K.to_numpy()[-3:])
        try:
            p, kinf = observed_order(hs, ks)
            print(f"  three finest levels: observed order p = {p:.2f}, Richardson K_inf = {kinf:.5e} "
                  f"(finest is {ks[-1] / kinf - 1:+.3%} from K_inf)")
        except ValueError:
            print("  three finest levels are not monotone; no order estimate")
        print()

    ratio = stu / stu_Q
    print(f"Student K/Q = {ratio.min():.5g} to {ratio.max():.5g} 1/m (spread {ratio.max() / ratio.min() - 1:.1e}); "
          f"mu/(hp dp) = {2.816e-4 / (HP * 1.0):.5g} for mu = 2.816e-4 Pa s, dp = 1 Pa (quarter-cell normalization)")

    print("\nBuilder difference at matching hmax, (K_student - K_chase)/K_chase:")
    for h in sorted(set(chase.index) & set(stu.index), reverse=True):
        print(f"  hmax = {h:g} um: {stu[h] / chase[h] - 1:+.4%}")
    print(f"Finest available: student {stu.index[-1]:g} um vs Chase {chase.index[-1]:g} um: "
          f"{stu.iloc[-1] / chase.iloc[-1] - 1:+.4%}")
    print(f"\nAnalytical flat-top K (Sangani-Acrivos + Brinkman, MODEL_CALCULATED) = {K_an:.4e} m^2")
    half_ulp = 0.5e-4 * 1e-11  # student K reported as d.dddd e-11
    print("Note: student values are rounded to 5 significant figures, so each relative change carries "
          f"a rounding uncertainty of about +/-{half_ulp / stu.iloc[-1] * math.sqrt(2):.1e} (absolute).")


if __name__ == "__main__":
    main()
