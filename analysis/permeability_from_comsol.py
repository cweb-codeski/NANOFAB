"""Permeability from a COMSOL Stage 1 (Creeping Flow, flat slip top) unit-cell export.

Equation (Darcy, superficial velocity over the full cross-section including pillars):
    K = mu * U * L / dP,   U = Q / A

Expected export (Results > Tables > Export, CSV or TXT), one row per solution:
    - a column with unit m^3/s : volumetric flow Q through the outlet face
                                 (Surface Integration of u over the outlet boundary)
    - a column with unit Pa    : applied pressure drop dP = p_in - p_out (global parameter)
    - optional columns d, l, h : geometry parameters if you ran a parametric sweep
    - optional column 'hmax' or 'mesh' : mesh-size parameter for a convergence study

Cell convention (--cell):
    quarter (default): domain l/2 long (flow direction) x l/2 wide x h tall -> L = l/2, A = (l/2) h
    full:              domain l x l x h                                   -> L = l,   A = l h

Usage:
    python analysis/permeability_from_comsol.py data/comsol_exports/stage1_flow.csv
    python analysis/permeability_from_comsol.py export.txt --mu 7.97e-4 --cell full
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import PROCESSED_DIR, find_column, geometry_columns
from nanofab.comsol_io import mesh_convergence_report, read_comsol_table
from nanofab.fluids import Fluid
from nanofab.params import load_parameters
from nanofab.permeability import permeability_from_flow, pillar_wick_permeability


def single_column_with_unit(table, unit: str, override: str | None) -> str:
    if override:
        return find_column(table, override)
    cols = table.columns_with_unit(unit)
    if len(cols) != 1:
        raise KeyError(f"Need exactly one column with unit {unit}, found {cols}; pass it explicitly")
    return cols[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", type=Path)
    ap.add_argument("--q-col", help="name of the outlet flow column (default: the m^3/s column)")
    ap.add_argument("--dp-col", help="name of the pressure-drop column (default: the Pa column)")
    ap.add_argument("--dp", type=float, help="pressure drop [Pa] if it is not a column")
    ap.add_argument("--mu", type=float, help="viscosity [Pa s] used in COMSOL (default: CoolProp at T_v)")
    ap.add_argument("--cell", choices=["quarter", "full"], default="quarter")
    ap.add_argument("--tolerance", type=float, default=0.02, help="mesh-convergence tolerance")
    args = ap.parse_args()

    ps = load_parameters()
    table = read_comsol_table(args.export)
    df = table.data

    q = df[single_column_with_unit(table, "m^3/s", args.q_col)].to_numpy(dtype=float)
    if args.dp is not None:
        dp = np.full_like(q, args.dp)
    else:
        dp = df[single_column_with_unit(table, "Pa", args.dp_col)].to_numpy(dtype=float)

    mu = args.mu if args.mu is not None else Fluid(ps["working_fluid"]).mu_l(ps["vapor_temperature"])
    geo = geometry_columns(table, ps)
    d, l, h = (np.broadcast_to(geo[k], q.shape).astype(float) for k in ("d", "l", "h"))

    if args.cell == "quarter":
        length, area = l / 2, (l / 2) * h
    else:
        length, area = l, l * h

    K_comsol = np.array([permeability_from_flow(qi, ai, li, dpi, mu) for qi, ai, li, dpi in zip(q, area, length, dp)])
    K_analytic = np.array([pillar_wick_permeability(di, pi, hi) for di, pi, hi in zip(d, l, h)])

    out = pd.DataFrame(
        {
            "d_um": d * 1e6,
            "l_um": l * 1e6,
            "h_um": h * 1e6,
            "dp_Pa": dp,
            "Q_m3_s": q,
            "K_comsol_m2": K_comsol,
            "K_analytic_m2": K_analytic,
            "diff_pct": 100 * (K_comsol - K_analytic) / K_analytic,
        }
    )
    mesh_col = find_column(table, "hmax", required=False) or find_column(table, "mesh", required=False)
    if mesh_col:
        out.insert(0, mesh_col, df[mesh_col].to_numpy())

    pd.set_option("display.width", 140)
    print(f"mu used: {mu:.4e} Pa s  (must match the viscosity set in COMSOL)")
    print(f"cell: {args.cell}\n")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\nAnalytical model: Sangani-Acrivos + Brinkman height correction (flat top).")
    print("Differences of 10-30% are plausible (analytical expansion accuracy, mesh); >50% -> check setup.")

    if mesh_col:
        # coarse -> fine: larger hmax is coarser
        ordered = out.sort_values(mesh_col, ascending=("hmax" not in mesh_col.lower()))
        values = {f"{v:g}": k for v, k in zip(ordered[mesh_col], ordered["K_comsol_m2"])}
        print("\nMesh convergence on K:")
        print(mesh_convergence_report(values, args.tolerance))

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    dest = PROCESSED_DIR / f"permeability_{args.export.stem}.csv"
    out.to_csv(dest, index=False)
    print(f"\nWrote {dest}")


if __name__ == "__main__":
    main()
