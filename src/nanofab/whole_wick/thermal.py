"""Stage D: in-plane heat conduction with evaporative removal.

    div(G grad T) + q''_in - q''_evap = 0,      q''_evap = h_e (T - T_sat)

T is the pillar-base plane temperature (T_e of vaartstra2019, assumption W6). G [W/K] is the
in-plane sheet conductance of substrate plus wick layer, G = k_sub t_sub + k_eff h. h_e already
contains the through-wick resistance from the unit cell, so no through-plane k_eff term is added.
This is the requested div(k_eff grad T) + q_input - q_evap = 0 written per unit footprint area.

Edge conditions: ``adiabatic`` (default; vaartstra2019 device symmetry) or a fixed temperature.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .fv import boundary_conductance, edge_cells, face_conductances, laplacian
from .grid import Grid


@dataclass
class EnergyReport:
    heat_input: float  # W
    electrical_heat: float  # W; always 0 while electrical power is unresolved (decision D14)
    evaporative: float  # W
    edge_conduction: float  # W leaving through fixed-temperature edges
    residual: float  # input - evaporative - edge
    relative_residual: float

    def lines(self) -> list[str]:
        return [
            f"heat input           {self.heat_input:.6e} W (electrical heating: not modeled, power unresolved)",
            f"evaporative removal  {self.evaporative:.6e} W",
            f"edge conduction out  {self.edge_conduction:.6e} W",
            f"residual             {self.residual:.3e} W (relative {self.relative_residual:.2e})",
        ]


@dataclass
class ThermalResult:
    T: np.ndarray  # K
    q_evap: np.ndarray  # W/m^2
    T_max: float
    T_max_location: tuple[float, float]  # (x, y) m
    energy: EnergyReport


def solve_thermal(
    grid: Grid,
    G: np.ndarray,
    h_e: np.ndarray,
    q_in: np.ndarray,
    T_sat: float,
    edges: dict[str, float] | None = None,
    q_electrical: np.ndarray | None = None,
    average: str = "harmonic",
) -> ThermalResult:
    """Solve for T [K].

    Args:
        G: in-plane sheet conductance per cell [W/K]
        h_e: evaporative coefficient per cell [W/m^2/K] (0 for dry cells)
        q_in: total areal heat input [W/m^2], including any electrical heating
        edges: {edge: T_fixed [K]} for fixed-temperature edges; all others adiabatic
        q_electrical: the electrical part of q_in, only for reporting
    """
    edges = edges or {}
    if not np.any(h_e > 0) and not edges:
        raise ValueError("No heat sink: every cell has h_e = 0 and all edges are adiabatic")
    A = grid.cell_area
    Cx, Cy = face_conductances(G, grid, average)
    M = laplacian(Cx, Cy, grid)
    diag = (h_e * A).ravel()
    rhs = (q_in * A + h_e * A * T_sat).ravel()
    edge_C = {}
    extra = np.zeros(grid.shape)
    rhs2 = rhs.reshape(grid.shape)
    for edge, T_b in edges.items():
        Cb = boundary_conductance(G, grid, edge)
        edge_C[edge] = Cb
        jj, ii = edge_cells(grid, edge)
        np.add.at(extra, (jj, ii), Cb)
        np.add.at(rhs2, (jj, ii), Cb * T_b)
    M = M + sp.diags(diag + extra.ravel())
    T = spla.spsolve(sp.csr_matrix(M), rhs2.ravel()).reshape(grid.shape)

    q_evap = h_e * (T - T_sat)
    Q_in = float((q_in * A).sum())
    Q_el = float((q_electrical * A).sum()) if q_electrical is not None else 0.0
    Q_ev = float((q_evap * A).sum())
    Q_edge = 0.0
    for edge, Cb in edge_C.items():
        jj, ii = edge_cells(grid, edge)
        Q_edge += float((Cb * (T[jj, ii] - edges[edge])).sum())
    res = Q_in - Q_ev - Q_edge
    k = int(np.argmax(T))
    j, i = np.unravel_index(k, grid.shape)
    return ThermalResult(
        T=T,
        q_evap=q_evap,
        T_max=float(T.max()),
        T_max_location=(float(grid.x[i]), float(grid.y[j])),
        energy=EnergyReport(Q_in, Q_el, Q_ev, Q_edge, res, abs(res) / max(abs(Q_in), 1e-300)),
    )
