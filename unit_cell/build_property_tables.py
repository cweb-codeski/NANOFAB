"""Convert a COMSOL Stage 5 sweep export into a whole-wick lookup table (long-format CSV).

The unit-cell model is the property generator; this script is the hand-off. Output columns:

    Pc_Pa,<property>     (SI; one row per sweep point)

Tables have a Pc axis only. Electrowetting changes Pc_max through the contact angle (closed form in
whole_wick/electrowetting.py), not K or h_e, so no voltage or field column is accepted.
For coated cell types the Pc sweep must reach 2*gamma/r_eff (the Pc_max at cos(theta_V) = 1).

with ``#`` header lines carrying provenance (COMSOL_EXTRACTED) and the source export/model.
Point a cell type at it in ``parameters/whole_wick/cell_types.yaml``::

    h_e: {method: table, file: parameters/whole_wick/tables/si_d10_l30_h25_h_e.csv, column: h_e,
          provenance: COMSOL_EXTRACTED, source: null,
          source_location: "comsol/unit_cell/stage5_sweep.mph, export data/comsol_exports/stage5.csv"}

Usage:
    python unit_cell/build_property_tables.py EXPORT.csv --property h_e \
        --pc-col "dp" --value-col "h_e" --out parameters/whole_wick/tables/X_h_e.csv

The export must already contain the property itself (for K, compute it in COMSOL as
mu*U*L/dp, or post-process the flow rate with analysis/permeability_from_comsol.py first).
Pc must be the capillary pressure p_v - p_l >= 0 (docs/whole_wick_assumptions.md, D7). If the
COMSOL sweep variable is the liquid gauge pressure, pass --negate-pc.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanofab.comsol_io import read_comsol_table  # noqa: E402

UNIT_FACTORS = {
    "Pc": {"Pa": 1.0, "kPa": 1e3},
    "K": {"m^2": 1.0, "m²": 1.0, "um^2": 1e-12, "µm^2": 1e-12},
    "h_e": {"W/(m^2*K)": 1.0, "W/(m²·K)": 1.0, "W/m^2/K": 1.0, "kW/(m^2*K)": 1e3},
}


def _factor(kind: str, unit: str, column: str) -> float:
    try:
        return UNIT_FACTORS[kind][unit]
    except KeyError:
        raise SystemExit(f"column '{column}' has unit '{unit}'; expected one of {list(UNIT_FACTORS[kind])}") from None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", type=Path)
    ap.add_argument("--property", choices=["K", "h_e"], required=True)
    ap.add_argument("--pc-col", required=True, help="keyword(s) matching the Pc column, space-separated")
    ap.add_argument("--value-col", required=True)
    ap.add_argument("--negate-pc", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    t = read_comsol_table(args.export)
    pc_c = t.column(*args.pc_col.split())
    v_c = t.column(*args.value_col.split())
    pc = t.data[pc_c].to_numpy(float) * _factor("Pc", t.units[pc_c], pc_c)
    if args.negate_pc:
        pc = -pc
    if np.any(pc < 0):
        raise SystemExit("negative Pc after conversion; check the sign convention (use --negate-pc?)")
    val = t.data[v_c].to_numpy(float) * _factor(args.property, t.units[v_c], v_c)
    cols = {"Pc_Pa": pc, args.property: val}
    if np.unique(pc).size != pc.size:
        raise SystemExit("repeated Pc values: the export must be a single Pc sweep (no voltage or field axis)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    model = t.metadata.get("Model", "unknown model")
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        fh.write(f"# provenance: COMSOL_EXTRACTED\n# source export: {args.export.as_posix()}\n# COMSOL model: {model}\n")
        fh.write(f"# columns SI; Pc = p_v - p_l >= 0\n{','.join(cols)}\n")
        for row in zip(*cols.values()):
            fh.write(",".join(f"{v:.10g}" for v in row) + "\n")
    print(f"wrote {args.out} ({len(val)} rows, Pc {pc.min():.4g}..{pc.max():.4g} Pa)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
