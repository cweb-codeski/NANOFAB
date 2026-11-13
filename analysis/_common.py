"""Shared helpers for the analysis scripts (path setup, column lookup, geometry resolution)."""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

from nanofab.comsol_io import ComsolTable  # noqa: E402
from nanofab.params import ParameterSet  # noqa: E402

PROCESSED_DIR = REPO / "data" / "processed"


def find_column(table: ComsolTable, name: str, required: bool = True) -> str | None:
    """Column whose name equals ``name`` (case-insensitive), else the unique one containing it."""
    cols = list(table.data.columns)
    exact = [c for c in cols if c.lower() == name.lower()]
    if exact:
        return exact[0]
    partial = [c for c in cols if name.lower() in c.lower()]
    if len(partial) == 1:
        return partial[0]
    if required:
        raise KeyError(f"No unique column for '{name}'. Columns: {cols} (units {table.units})")
    return None


def geometry_columns(table: ComsolTable, ps: ParameterSet, case: str | None = None) -> dict[str, "object"]:
    """Per-row pillar d, l, h [m]: from sweep columns d / l / hp (or h) if present, else a named case.

    ``case`` is a key of parameters/unit_cell/cases.yaml (e.g. ``proposed_d20_l50_h40``). Without it,
    missing columns fall back to the reference geometry in parameters/geometry.yaml and a warning is
    printed, because footprint and cross-section areas would silently be wrong for another geometry.
    COMSOL reserves ``h`` (mesh element size), so the COMSOL guide names the pillar height ``hp``.
    """
    from nanofab.geometry_cases import load_geometry_cases

    fallback = None
    if case is not None:
        c = load_geometry_cases()[case]
        fallback = {"d": c.d, "l": c.l, "h": c.hp}
    out = {}
    used_default = []
    aliases = {"d": ("d",), "l": ("l",), "h": ("hp", "h")}
    for key, pname in (("d", "pillar_diameter"), ("l", "pillar_pitch"), ("h", "pillar_height")):
        col = next((c for c in table.data.columns if c.lower() in aliases[key]), None)
        if col is not None:
            unit = table.units.get(col, "")
            scale = {"m": 1.0, "um": 1e-6, "µm": 1e-6, "mm": 1e-3, "": 1.0}.get(unit)
            if scale is None:
                raise ValueError(f"Column {col} has unsupported unit '{unit}'")
            out[key] = table.data[col].to_numpy(dtype=float) * scale
        elif fallback is not None:
            out[key] = fallback[key]
        else:
            out[key] = ps[pname]
            used_default.append(key)
    if used_default:
        print(f"WARNING: no {'/'.join(used_default)} column in the export and no --case given; using the "
              "REFERENCE geometry from parameters/geometry.yaml. Pass --case for any other geometry.")
    return out
