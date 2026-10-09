"""Read tables and data exported from COMSOL (CSV or TXT).

COMSOL exports start with comment lines beginning with ``%``. Metadata lines look like
``% Model:   stage1_flow.mph`` and the last comment line before the numbers holds the column
headers, e.g. ``% d (m),l (m),Surface Integration 1 (m^3/s)`` (CSV) or columns separated by two or
more spaces (TXT). Units in parentheses are split off into ``ComsolTable.units``.

Both "Results > Tables > Export" and "Results > Export > Data" files are handled.
"""

from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

_UNIT_RE = re.compile(r"^(.*?)\s*\(([^()]*)\)\s*$")


@dataclass
class ComsolTable:
    data: pd.DataFrame
    units: dict[str, str] = field(default_factory=dict)  # column name -> unit string ('' if none)
    metadata: dict[str, str] = field(default_factory=dict)
    path: Path | None = None

    def column(self, *keywords: str, unit: str | None = None) -> str:
        """Find the single column whose name contains all ``keywords`` (case-insensitive).

        Optionally require an exact unit (e.g. ``unit="Pa"``). Raises if zero or several match.
        """
        kws = [k.lower() for k in keywords]
        matches = [
            c
            for c in self.data.columns
            if all(k in c.lower() for k in kws) and (unit is None or self.units.get(c) == unit)
        ]
        if len(matches) != 1:
            raise KeyError(
                f"Expected exactly one column matching {keywords} unit={unit!r}, found {matches}. "
                f"Available: {list(self.data.columns)} with units {self.units}"
            )
        return matches[0]

    def columns_with_unit(self, unit: str) -> list[str]:
        return [c for c, u in self.units.items() if u == unit]


def _split_header(header: str, delimiter: str | None) -> list[str]:
    if delimiter == ",":
        # COMSOL 6.x quotes names that contain commas, e.g. "Q, flow rate (m^3/s)".
        parts = next(csv.reader([header]))
    else:
        parts = re.split(r"\s{2,}|\t", header)
    return [p.strip() for p in parts if p.strip()]


def _split_name_unit(raw: str) -> tuple[str, str]:
    m = _UNIT_RE.match(raw)
    if m:
        return m.group(1).strip(), m.group(2).strip()
    return raw.strip(), ""


def read_comsol_table(path: str | Path) -> ComsolTable:
    """Parse a COMSOL CSV/TXT export into a DataFrame with units and metadata."""
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = text.splitlines()

    comment_lines: list[str] = []
    data_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("%"):
            if data_lines:
                raise ValueError(f"{path}: comment line after data starts; multiple tables in one file?")
            comment_lines.append(stripped.lstrip("%").strip())
        else:
            data_lines.append(stripped)
    if not data_lines:
        raise ValueError(f"{path}: no numeric data found")

    delimiter = "," if "," in data_lines[0] else None

    metadata: dict[str, str] = {}
    header = None
    for c in comment_lines[:-1]:
        # "% Model:  x.mph" (TXT, older CSV) or "% Model,x.mph" (COMSOL 6.x CSV)
        key, sep, value = c.partition(":")
        if not (sep and re.fullmatch(r"[A-Za-z ]+", key.strip())):
            key, sep, value = c.partition(",")
        if sep and re.fullmatch(r"[A-Za-z ]+", key.strip()) and len(key.strip()) < 30:
            metadata[key.strip()] = value.strip().strip('"')
    if comment_lines:
        header = comment_lines[-1]

    n_cols = len(data_lines[0].split(",")) if delimiter == "," else len(data_lines[0].split())
    raw_names = _split_header(header, delimiter) if header else []
    if len(raw_names) != n_cols:
        raw_names = [f"col{i}" for i in range(n_cols)]

    names, units = [], {}
    for raw in raw_names:
        name, unit = _split_name_unit(raw)
        base, k = name, 2
        while name in units:
            name = f"{base}_{k}"
            k += 1
        names.append(name)
        units[name] = unit

    df = pd.read_csv(
        io.StringIO("\n".join(data_lines)),
        sep="," if delimiter == "," else r"\s+",
        header=None,
        names=names,
        engine="python",
    )
    return ComsolTable(data=df, units=units, metadata=metadata, path=path)


def relative_change(fine: float, coarser: float) -> float:
    """|X_fine - X_coarser| / |X_fine| for mesh-convergence tables."""
    if fine == 0:
        raise ZeroDivisionError("fine-mesh value is zero; use an absolute tolerance instead")
    return abs(fine - coarser) / abs(fine)


def mesh_convergence_report(values: dict[str, float], tolerance: float = 0.02) -> str:
    """Text report for an ordered mapping {mesh_label: metric}, coarse -> fine.

    Each line shows the change relative to the next finer mesh and whether it is within tolerance.
    """
    labels = list(values)
    if len(labels) < 2:
        raise ValueError("Need at least two meshes")
    fine_label = labels[-1]
    lines = [f"{'mesh':>10}  {'value':>14}  {'rel. change vs finer':>22}"]
    for i, label in enumerate(labels):
        if i < len(labels) - 1:
            rc = relative_change(values[labels[i + 1]], values[label])
            flag = "OK" if rc < tolerance else "refine"
            lines.append(f"{label:>10}  {values[label]:14.6g}  {rc:21.2%} {flag}")
        else:
            lines.append(f"{label:>10}  {values[label]:14.6g}  {'(finest)':>22}")
    rc_last = relative_change(values[fine_label], values[labels[-2]])
    verdict = "converged" if rc_last < tolerance else "NOT converged"
    lines.append(f"Finest pair change {rc_last:.2%} -> {verdict} at tolerance {100 * tolerance:g}%")
    return "\n".join(lines)
