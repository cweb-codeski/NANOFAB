"""Shared finite-volume pieces: face averaging and the 5-point operator on the cell grid."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .grid import EDGES, Grid

AVERAGES = ("harmonic", "arithmetic")


def face_average(a: np.ndarray, b: np.ndarray, kind: str) -> np.ndarray:
    """Face value between two cells (decision D5). Harmonic keeps normal flux continuous."""
    if kind == "harmonic":
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where((a > 0) & (b > 0), 2 * a * b / (a + b), 0.0)
        return out
    if kind == "arithmetic":
        return 0.5 * (a + b)
    raise ValueError(f"face average must be one of {AVERAGES}")


def face_conductances(c: np.ndarray, grid: Grid, kind: str) -> tuple[np.ndarray, np.ndarray]:
    """Interior face conductances for a per-cell coefficient ``c`` (flux = C * difference).

    Returns Cx with shape (ny, nx-1) between [j, i] and [j, i+1], and Cy with shape (ny-1, nx).
    """
    Cx = face_average(c[:, :-1], c[:, 1:], kind) * grid.dy / grid.dx
    Cy = face_average(c[:-1, :], c[1:, :], kind) * grid.dx / grid.dy
    return Cx, Cy


def boundary_conductance(c: np.ndarray, grid: Grid, edge: str) -> np.ndarray:
    """Half-cell conductance from cell centers to one edge (1D array along that edge)."""
    if edge == "left":
        return c[:, 0] * grid.dy / (grid.dx / 2)
    if edge == "right":
        return c[:, -1] * grid.dy / (grid.dx / 2)
    if edge == "bottom":
        return c[0, :] * grid.dx / (grid.dy / 2)
    if edge == "top":
        return c[-1, :] * grid.dx / (grid.dy / 2)
    raise ValueError(f"edge must be one of {EDGES}")


def edge_cells(grid: Grid, edge: str) -> tuple[np.ndarray, np.ndarray]:
    """(j, i) index arrays of the cells adjacent to ``edge``."""
    ny, nx = grid.shape
    if edge == "left":
        return np.arange(ny), np.zeros(ny, int)
    if edge == "right":
        return np.arange(ny), np.full(ny, nx - 1)
    if edge == "bottom":
        return np.zeros(nx, int), np.arange(nx)
    if edge == "top":
        return np.full(nx, ny - 1), np.arange(nx)
    raise ValueError(f"edge must be one of {EDGES}")


def laplacian(Cx: np.ndarray, Cy: np.ndarray, grid: Grid) -> sp.csr_matrix:
    """Matrix of sum_faces C (u_i - u_j) (symmetric, positive semidefinite)."""
    ny, nx = grid.shape
    n = nx * ny
    idx = np.arange(n).reshape(ny, nx)
    rows, cols, vals = [], [], []
    diag = np.zeros(n)
    for C, a, b in ((Cx, idx[:, :-1], idx[:, 1:]), (Cy, idx[:-1, :], idx[1:, :])):
        a, b, C = a.ravel(), b.ravel(), C.ravel()
        rows += [a, b]
        cols += [b, a]
        vals += [-C, -C]
        np.add.at(diag, a, C)
        np.add.at(diag, b, C)
    rows.append(np.arange(n))
    cols.append(np.arange(n))
    vals.append(diag)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def net_face_inflow(u: np.ndarray, Cx: np.ndarray, Cy: np.ndarray) -> np.ndarray:
    """Per-cell sum over interior faces of C (u_neighbor - u_cell), recomputed from face fluxes."""
    out = np.zeros_like(u)
    fx = Cx * (u[:, :-1] - u[:, 1:])  # flow from i to i+1
    fy = Cy * (u[:-1, :] - u[1:, :])
    out[:, :-1] -= fx
    out[:, 1:] += fx
    out[:-1, :] -= fy
    out[1:, :] += fy
    return out


def conductance_sum(Cx: np.ndarray, Cy: np.ndarray) -> np.ndarray:
    """Sum of interior face conductances touching each cell."""
    ny = Cy.shape[0] + 1 if Cy.size else Cx.shape[0]
    nx = Cx.shape[1] + 1 if Cx.size else Cy.shape[1]
    s = np.zeros((ny, nx))
    s[:, :-1] += Cx
    s[:, 1:] += Cx
    s[:-1, :] += Cy
    s[1:, :] += Cy
    return s
