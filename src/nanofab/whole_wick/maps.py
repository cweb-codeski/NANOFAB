"""Spatial input maps on the grid: heat flux, cell type (geometry + material), applied voltage.

Every map spec is a mapping with a ``type`` key. Numeric entries are provenance-checked leaves
(see ``config.py``). Heat-flux maps are in W/m^2 per footprint area; voltage maps in V.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from ..params import _to_si
from ..provenance import REPO_ROOT
from .config import si
from .grid import Grid


def _resample_file(path: str, grid: Grid, unit: str) -> np.ndarray:
    """Read a numeric CSV matrix (row 0 = y = 0 edge) and bilinearly resample it onto ``grid``.

    The file's cells are assumed to tile the full wick footprint uniformly.
    """
    full = Path(path) if Path(path).is_absolute() else REPO_ROOT / path
    data = np.loadtxt(full, delimiter=",", comments="#", ndmin=2)
    factor = _to_si(1.0, unit, str(path))[0]
    ny, nx = data.shape
    if (ny, nx) == grid.shape:
        return data * factor
    xs = (np.arange(nx) + 0.5) * grid.lx / nx
    ys = (np.arange(ny) + 0.5) * grid.ly / ny
    if nx == 1 or ny == 1:
        raise ValueError(f"{path}: need at least 2x2 values to resample")
    interp = RegularGridInterpolator((ys, xs), data, bounds_error=False, fill_value=None)
    X, Y = grid.XY
    pts = np.column_stack([np.clip(Y.ravel(), ys[0], ys[-1]), np.clip(X.ravel(), xs[0], xs[-1])])
    return interp(pts).reshape(grid.shape) * factor


def _rect_mask(grid: Grid, rect: dict) -> np.ndarray:
    X, Y = grid.XY
    return (
        (X >= si(rect["x_min"])) & (X <= si(rect["x_max"])) & (Y >= si(rect["y_min"])) & (Y <= si(rect["y_max"]))
    )


def heat_flux_map(spec: dict, grid: Grid) -> np.ndarray:
    """q''_in [W/m^2]. Types: uniform, gaussian, patches, file."""
    kind = spec["type"]
    if kind == "uniform":
        return np.full(grid.shape, si(spec["q"]))
    if kind == "gaussian":
        q = np.full(grid.shape, si(spec["background"]))
        X, Y = grid.XY
        for hs in spec["hotspots"]:
            r2 = (X - si(hs["x"])) ** 2 + (Y - si(hs["y"])) ** 2
            q = q + si(hs["peak_increment"]) * np.exp(-r2 / (2 * si(hs["sigma"]) ** 2))
        return q
    if kind == "patches":
        q = np.full(grid.shape, si(spec["background"]))
        for patch in spec["patches"]:
            q[_rect_mask(grid, patch)] = si(patch["q"])
        return q
    if kind == "file":
        return _resample_file(spec["path"], grid, spec["unit"])
    raise ValueError(f"unknown heat-flux map type '{kind}'")


def cell_type_map(spec: dict, grid: Grid, names: list[str], feed_edges: list[str]) -> np.ndarray:
    """Integer index into ``names`` per cell. Types: uniform, feed_distance_bands, rectangles, file."""
    index = {n: i for i, n in enumerate(names)}

    def idx(name: str) -> int:
        if name not in index:
            raise KeyError(f"cell type '{name}' not defined; known: {names}")
        return index[name]

    kind = spec["type"]
    if kind == "uniform":
        return np.full(grid.shape, idx(spec["cell_type"]), dtype=int)
    if kind == "feed_distance_bands":
        dist = grid.distance_to_edges(feed_edges)
        out = np.full(grid.shape, -1, dtype=int)
        lower = 0.0
        for band in spec["bands"]:
            upper = si(band["max_distance"]) if "max_distance" in band else np.inf
            out[(dist >= lower) & (dist < upper) & (out < 0)] = idx(band["cell_type"])
            lower = upper
        if (out < 0).any():
            raise ValueError("feed_distance_bands leave cells unassigned; give the last band no max_distance")
        return out
    if kind == "rectangles":
        out = np.full(grid.shape, idx(spec["background"]), dtype=int)
        for rect in spec["regions"]:
            out[_rect_mask(grid, rect)] = idx(rect["cell_type"])
        return out
    if kind == "file":
        full = REPO_ROOT / spec["path"]
        rows = [r.split(",") for r in full.read_text(encoding="utf-8").splitlines() if r and not r.startswith("#")]
        data = np.array([[idx(c.strip()) for c in r] for r in rows], dtype=int)
        if data.shape != grid.shape:
            raise ValueError(f"{spec['path']}: shape {data.shape} != grid {grid.shape} (no resampling of type maps)")
        return data
    raise ValueError(f"unknown cell-type map type '{kind}'")


ACTUATION = {"zero": "none", "uniform": "global", "vulnerability_mask": "selective", "rectangles": "custom", "file": "custom"}


def voltage_map(spec: dict | None, grid: Grid, mask: np.ndarray | None = None) -> np.ndarray:
    """Applied electrowetting voltage V [V] per cell.

    Types: zero (passive), uniform (global actuation), vulnerability_mask (selective actuation;
    needs ``mask``), rectangles, file.
    """
    if spec is None or spec["type"] == "zero":
        return np.zeros(grid.shape)
    kind = spec["type"]
    if kind == "uniform":
        return np.full(grid.shape, si(spec["V"]))
    if kind == "rectangles":
        V = np.zeros(grid.shape)
        for rect in spec["regions"]:
            V[_rect_mask(grid, rect)] = si(rect["V"])
        return V
    if kind == "file":
        return _resample_file(spec["path"], grid, spec["unit"])
    if kind == "vulnerability_mask":
        if mask is None:
            raise ValueError("vulnerability_mask voltage map needs a mask from a reference case (cases.py)")
        return np.where(mask, si(spec["V"]), 0.0)
    raise ValueError(f"unknown voltage map type '{kind}'")


def vulnerability_mask(pressure_ratio: np.ndarray, threshold: float, feed_cells: np.ndarray | None = None) -> np.ndarray:
    """Cells whose required Pc / Pc_max(V=0) is at or above ``threshold`` (decision D11)."""
    mask = pressure_ratio >= threshold
    if feed_cells is not None:
        mask &= ~feed_cells
    return mask
