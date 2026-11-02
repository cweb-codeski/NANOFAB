"""Stage C: depth-averaged Darcy flow of the liquid in the wick plane.

Governing equations (2D planform, wick height h; docs/whole_wick_assumptions.md W1, W5):

    u_l = -(K/mu_l) (grad p_l - rho_l g)             superficial velocity [m/s]
    div(rho_l h u_l) = -m''_evap                     m''_evap [kg/m^2/s] per footprint area

which is div(rho_l h u_l) = -m_dot_evap with the volumetric sink m_dot_evap = m''_evap / h.

Pressure convention (D7): p is gauge relative to the uniform vapor pressure p_v, so p = p_l - p_v
and the capillary pressure is Pc = -p >= 0. Feeds are fixed-pressure edges p = -Pc_feed, and
optionally fixed-pressure cells (distributed feed arteries). Every other edge is no-flux.

Finite volume: cell-centered unknowns, face mobility rho K h / mu averaged harmonically (D5),
half-cell conductance to Dirichlet edges, sparse direct solve.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .fv import boundary_conductance, conductance_sum, edge_cells, face_conductances, laplacian
from .grid import Grid


@dataclass
class FeedSpec:
    edges: dict[str, float] = field(default_factory=dict)  # edge -> Pc_feed [Pa]
    cells: np.ndarray | None = None  # bool mask of fixed-pressure feed cells
    cell_pc: float = 0.0  # Pc_feed of the feed cells [Pa]
    max_mass_flow: float | None = None  # optional finite feed capacity [kg/s] (independent supply limit)

    def feed_cell_mask(self, shape) -> np.ndarray:
        return np.zeros(shape, bool) if self.cells is None else self.cells.astype(bool)


@dataclass
class MassReport:
    inlet: float  # kg/s entering through feeds
    outlet: float  # kg/s leaving through feeds (> 0 only if some feed acts as a sink)
    evaporation: float  # kg/s consumed by evaporation
    residual: float  # inlet - outlet - evaporation
    relative_residual: float
    max_cell_residual: float  # kg/s, largest per-cell imbalance recomputed from face fluxes

    def lines(self) -> list[str]:
        return [
            f"inlet flow           {self.inlet:.6e} kg/s",
            f"outlet flow          {self.outlet:.6e} kg/s",
            f"evaporation          {self.evaporation:.6e} kg/s",
            f"residual imbalance   {self.residual:.3e} kg/s (relative {self.relative_residual:.2e})",
            f"max cell imbalance   {self.max_cell_residual:.3e} kg/s",
        ]


@dataclass
class DarcyResult:
    p: np.ndarray  # gauge liquid pressure p_l - p_v [Pa]
    Pc: np.ndarray  # p_v - p_l [Pa]
    ux: np.ndarray  # cell-centered superficial velocity [m/s]
    uy: np.ndarray
    Cx: np.ndarray  # interior face transmissibilities [kg/s/Pa]
    Cy: np.ndarray
    boundary_C: dict[str, np.ndarray]
    conductance_total: np.ndarray  # sum of face + feed-edge transmissibilities per cell
    feed_cells: np.ndarray
    mass: MassReport


def solve_darcy(
    grid: Grid,
    K: np.ndarray,
    h: np.ndarray,
    rho: np.ndarray,
    mu: np.ndarray,
    m_evap: np.ndarray,
    feeds: FeedSpec,
    gravity: tuple[float, float] = (0.0, 0.0),
    average: str = "harmonic",
) -> DarcyResult:
    """Solve for the liquid pressure given the evaporation mass flux m_evap [kg/m^2/s]."""
    if not feeds.edges and (feeds.cells is None or not feeds.cells.any()):
        raise ValueError("Darcy problem needs at least one feed (edge or cell); otherwise it is singular")
    ny, nx = grid.shape
    A = grid.cell_area
    mob = rho * K * h / mu  # kg/(Pa s); face transmissibility = mob_f * dy/dx
    Cx, Cy = face_conductances(mob, grid, average)

    # gravity source: inflow from j to i of C rho_f g.(r_i - r_j)
    gx, gy = gravity
    rho_x = 0.5 * (rho[:, :-1] + rho[:, 1:])
    rho_y = 0.5 * (rho[:-1, :] + rho[1:, :])
    gflow_x = Cx * rho_x * gx * grid.dx  # flow from i to i+1 driven by gravity
    gflow_y = Cy * rho_y * gy * grid.dy
    g_in = np.zeros(grid.shape)
    g_in[:, :-1] -= gflow_x
    g_in[:, 1:] += gflow_x
    g_in[:-1, :] -= gflow_y
    g_in[1:, :] += gflow_y

    M = laplacian(Cx, Cy, grid).tolil()
    rhs = (-m_evap * A + g_in).ravel()
    diag_add = np.zeros(grid.shape)
    boundary_C: dict[str, np.ndarray] = {}
    for edge, pc_feed in feeds.edges.items():
        Cb = boundary_conductance(mob, grid, edge)
        boundary_C[edge] = Cb
        jj, ii = edge_cells(grid, edge)
        np.add.at(diag_add, (jj, ii), Cb)
        rhs_2d = rhs.reshape(grid.shape)
        np.add.at(rhs_2d, (jj, ii), Cb * (-pc_feed))
        rhs = rhs_2d.ravel()
    M.setdiag(M.diagonal() + diag_add.ravel())

    fc = feeds.feed_cell_mask(grid.shape)
    if fc.any():
        rows = np.flatnonzero(fc.ravel())
        for r in rows:
            M.rows[r] = [r]
            M.data[r] = [1.0]
        rhs[rows] = -feeds.cell_pc
    p = spla.spsolve(sp.csr_matrix(M), rhs).reshape(grid.shape)

    # velocities (superficial, cell-centered average of face values)
    fx = Cx * (p[:, :-1] - p[:, 1:]) + gflow_x  # kg/s from i to i+1
    fy = Cy * (p[:-1, :] - p[1:, :]) + gflow_y
    ux = np.zeros(grid.shape)
    uy = np.zeros(grid.shape)
    face_ux = fx / (rho_x * 0.5 * (h[:, :-1] + h[:, 1:]) * grid.dy)
    face_uy = fy / (rho_y * 0.5 * (h[:-1, :] + h[1:, :]) * grid.dx)
    cnt_x = np.zeros(grid.shape)
    cnt_y = np.zeros(grid.shape)
    ux[:, :-1] += face_ux
    ux[:, 1:] += face_ux
    cnt_x[:, :-1] += 1
    cnt_x[:, 1:] += 1
    uy[:-1, :] += face_uy
    uy[1:, :] += face_uy
    cnt_y[:-1, :] += 1
    cnt_y[1:, :] += 1
    ux = np.divide(ux, cnt_x, out=np.zeros_like(ux), where=cnt_x > 0)
    uy = np.divide(uy, cnt_y, out=np.zeros_like(uy), where=cnt_y > 0)

    # mass audit from face fluxes (independent of the linear solve)
    net_in = np.zeros(grid.shape)
    net_in[:, :-1] -= fx
    net_in[:, 1:] += fx
    net_in[:-1, :] -= fy
    net_in[1:, :] += fy
    feed_in = np.zeros(grid.shape)
    for edge, Cb in boundary_C.items():
        jj, ii = edge_cells(grid, edge)
        np.add.at(feed_in, (jj, ii), Cb * (-feeds.edges[edge] - p[jj, ii]))
    evap = m_evap * A
    cell_res = net_in + feed_in - evap
    # fixed-pressure feed cells: the reservoir supplies whatever balances the cell
    reservoir = np.where(fc, -cell_res, 0.0)
    cell_res = np.where(fc, 0.0, cell_res)
    supply = np.concatenate([feed_in[feed_in != 0].ravel(), reservoir[fc].ravel()]) if fc.any() else feed_in.ravel()
    inlet = float(supply[supply > 0].sum())
    outlet = float(-supply[supply < 0].sum()) + 0.0
    evaporation = float(evap.sum())
    residual = inlet - outlet - evaporation
    scale = max(abs(evaporation), inlet, 1e-300)

    cond_total = conductance_sum(Cx, Cy)
    for edge, Cb in boundary_C.items():
        jj, ii = edge_cells(grid, edge)
        np.add.at(cond_total, (jj, ii), Cb)

    return DarcyResult(
        p=p,
        Pc=-p,
        ux=ux,
        uy=uy,
        Cx=Cx,
        Cy=Cy,
        boundary_C=boundary_C,
        conductance_total=cond_total,
        feed_cells=fc,
        mass=MassReport(
            inlet=inlet,
            outlet=outlet,
            evaporation=evaporation,
            residual=residual,
            relative_residual=abs(residual) / scale,
            max_cell_residual=float(np.abs(cell_res).max()),
        ),
    )
