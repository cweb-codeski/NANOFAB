"""CSV / markdown / figure output for whole-wick runs."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from ..provenance import report_category
from .dryout import margins
from .model import CoupledResult


def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def write_columns(cols: dict[str, np.ndarray], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(cols)
    n = len(next(iter(cols.values())))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for i in range(n):
            w.writerow([cols[k][i] for k in keys])


def _fmt(v, spec: str) -> str:
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (float, np.floating)):
        if math.isnan(v):
            return "-"
        return format(v, spec)
    return str(v)


def markdown_table(rows: list[dict], columns: list[tuple[str, str, str]]) -> str:
    """columns: (key, header, format spec)."""
    head = "| " + " | ".join(h for _, h, _ in columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = ["| " + " | ".join(_fmt(r.get(k, ""), f) for k, _, f in columns) + " |" for r in rows]
    return "\n".join([head, sep, *body])


def provenance_rows(registry: dict) -> list[dict]:
    rows = []
    for path, p in sorted(registry.items()):
        rows.append(
            {
                "parameter": path,
                "value": p.raw_value if not isinstance(p.raw_value, list) else str(p.raw_value),
                "unit": p.raw_unit,
                "provenance": p.provenance,
                "category": report_category(p.provenance),
                "source": p.source,
                "source_location": p.source_location,
            }
        )
    return rows


def plot_result(result: CoupledResult, path: Path, title: str = "") -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = result.model
    g = m.grid
    ext = [0, g.lx * 1e3, 0, g.ly * 1e3]
    mg = margins(result)
    panels = [
        ("heat input q''_in [W/cm²]", result.q_in / 1e4, "inferno"),
        ("cell type", m.type_map.astype(float), "tab10"),
        ("applied voltage V [V]", m.V, "viridis"),
        ("Pc_max(V) [Pa]", result.props.Pc_max, "viridis"),
        ("T - T_sat [K]", result.thermal.T - m.fluid.T_sat, "inferno"),
        ("q''_evap [W/cm²]", result.thermal.q_evap / 1e4, "magma"),
        ("required p_v - p_l [Pa]", result.Pc, "cividis"),
        ("required Δp / Pc_max(V)", mg.pressure_ratio, "RdYlGn_r"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(17, 7), constrained_layout=True)
    for ax, (name, data, cmap) in zip(axes.ravel(), panels):
        kw = {"vmin": 0, "vmax": 1} if name.startswith("required Δp") else {}
        im = ax.imshow(data, origin="lower", extent=ext, cmap=cmap, aspect="auto", **kw)
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        fig.colorbar(im, ax=ax)
        if result.dry.any() and name.startswith("required Δp"):
            ax.contour(g.x * 1e3, g.y * 1e3, result.dry.astype(float), levels=[0.5], colors="k")
    fig.suptitle(title)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)
