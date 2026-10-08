"""Effective evaporative heat-transfer coefficient from a COMSOL Stage 4/5 unit-cell export.

Definitions (footprint-based, so results feed a device-scale model directly):
    q''      = Q_in / A_fp                       A_fp = footprint of the modeled cell
    h_evap   = q'' / (T_base - T_v)              T_base = mean temperature at the pillar base
                                                 (top of substrate), T_v = vapor temperature
    R''_wick = (T_base - T_lv) / q''             conduction through pillars + liquid + film
    R''_int  = (T_lv - T_v) / q''                kinetic interfacial resistance
    energy balance error = |Q_in - Q_evap| / Q_in   (should be < 1%)

This matches Vaartstra et al. (2019) Eq. (4), h_e = q''_in / (T_e - T_sat), with q''_in per
footprint area and T_e the average temperature of the evaporator surface (T_base here; see the
COMSOL guide, Stage 2, for which surface that is).

Expected export columns (names matched case-insensitively; override with flags):
    Q_in   (W)  heat entering through the heated bottom boundary (Surface Integration of the
                normal conductive heat flux, or q0 * area)
    Q_evap (W)  heat leaving through the meniscus (Surface Integration of the Schrage heat flux)
    T_base (K)  Surface Average of T on the pillar-base / substrate-top boundary
    T_lv   (K)  Surface Average of T on the meniscus
    optional: T_max (K), and sweep columns such as Pc (Pa), d / l / h (m), accommodation (1)

Usage:
    python analysis/heat_transfer_from_comsol.py data/comsol_exports/stage4_thermal.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import PROCESSED_DIR, find_column, geometry_columns
from nanofab.comsol_io import read_comsol_table
from nanofab.params import load_parameters


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", type=Path)
    ap.add_argument("--qin-col", default="Q_in")
    ap.add_argument("--qevap-col", default="Q_evap")
    ap.add_argument("--base-col", default="T_base")
    ap.add_argument("--lv-col", default="T_lv")
    ap.add_argument("--T-vapor", type=float, help="vapor temperature [K] (default: parameters)")
    ap.add_argument("--cell", choices=["quarter", "full"], default="quarter")
    ap.add_argument("--energy-tol", type=float, default=0.01)
    args = ap.parse_args()

    ps = load_parameters()
    table = read_comsol_table(args.export)
    df = table.data
    T_v = args.T_vapor if args.T_vapor is not None else ps["vapor_temperature"]

    qin_col = find_column(table, args.qin_col)
    base_col = find_column(table, args.base_col)
    lv_col = find_column(table, args.lv_col)
    qevap_col = find_column(table, args.qevap_col, required=False)

    Q_in = df[qin_col].to_numpy(dtype=float)
    T_base = df[base_col].to_numpy(dtype=float)
    T_lv = df[lv_col].to_numpy(dtype=float)
    Q_evap = df[qevap_col].to_numpy(dtype=float) if qevap_col else None

    l = np.broadcast_to(geometry_columns(table, ps)["l"], Q_in.shape).astype(float)
    A_fp = (l / 2) ** 2 if args.cell == "quarter" else l**2

    q = Q_in / A_fp
    used = {qin_col, base_col, lv_col, qevap_col}
    out = pd.DataFrame({c: df[c] for c in df.columns if c not in used})
    out["q_W_cm2"] = q / 1e4
    out["T_base_minus_Tv_K"] = T_base - T_v
    out["h_evap_W_m2K"] = q / (T_base - T_v)
    out["R_wick_mm2K_W"] = (T_base - T_lv) / q * 1e6
    out["R_int_mm2K_W"] = (T_lv - T_v) / q * 1e6
    out["kinetic_fraction"] = (T_lv - T_v) / (T_base - T_v)

    if Q_evap is not None:
        err = np.abs(Q_in - Q_evap) / np.abs(Q_in)
        out["energy_balance_err"] = err
        bad = err > args.energy_tol
        if bad.any():
            print(f"WARNING: {bad.sum()} row(s) violate the energy balance by more than {args.energy_tol:.0%}."
                  " Refine the mesh near the meniscus/contact line or check boundary assignments.")

    if (T_base <= T_v).any():
        print("WARNING: some rows have T_base <= T_v; h_evap is undefined there (check T_v and BC signs).")

    pd.set_option("display.width", 160)
    print(f"T_v = {T_v:.2f} K, cell = {args.cell}, footprint = {'(l/2)^2' if args.cell == 'quarter' else 'l^2'}")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\nkinetic_fraction near 1 -> interface-limited (accommodation coefficient dominates);"
          " near 0 -> conduction-limited (geometry/material dominates).")

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    dest = PROCESSED_DIR / f"h_evap_{args.export.stem}.csv"
    out.to_csv(dest, index=False)
    print(f"\nWrote {dest} (use as h_evap lookup table vs the sweep columns)")


if __name__ == "__main__":
    main()
