"""Two-dimensional x-z substrate conduction, cell-centered finite volumes (PROJECT_DERIVED numerics).

Solves div(k_s grad T) = 0 on 0 <= x <= L, -t_s <= z <= 0 (invariant in y), with conductive flux
q_c = -k_s grad T and outward normal n:

* bottom (z = -t_s): q_c . n = -q_in''(x)                (heat enters)
* top    (z = 0):    q_c . n =  q_e''(x) = h_e(x) (T_e(x) - T_v)   (V2019 Eq. 4 / Sec. 2.3)
* x = 0 and x = L:   q_c . n = 0                          (V2019 Sec. 2.3 symmetry/insulation)

The top Robin condition uses the boundary-face temperature T_e, eliminated through the half-cell
conductive resistance dz/(2 k_s): q_e = (T_c - T_v) / (1/h_e + dz/(2 k_s)). T_e is reconstructed
as T_c - q_e dz/(2 k_s), never reported as a cell-centre value. The bottom-face temperature
T_source = T_c + q_in dz/(2 k_s) likewise.

Unknowns are theta = T - T_v, which keeps the linear system well scaled. Fluxes are per unit width
(W cancels); multiply by W for totals.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass(frozen=True)
class SlabMesh:
    L: float
    t_s: float
    Nx: int
    Nz: int

    def __post_init__(self):
        if self.Nx < 1 or self.Nz < 1:
            raise ValueError("Nx and Nz must be >= 1")

    @property
    def x_faces(self):
        return np.linspace(0.0, self.L, self.Nx + 1)

    @property
    def z_faces(self):
        return np.linspace(-self.t_s, 0.0, self.Nz + 1)

    @property
    def dx(self):
        return np.diff(self.x_faces)

    @property
    def dz(self):
        return self.t_s / self.Nz

    @property
    def xc(self):
        f = self.x_faces
        return 0.5 * (f[1:] + f[:-1])

    @property
    def zc(self):
        f = self.z_faces
        return 0.5 * (f[1:] + f[:-1])


@dataclass
class SlabSolution:
    mesh: SlabMesh
    theta: np.ndarray  # (Nx, Nz) cell-centre T - T_v [K]; j = 0 is the bottom row
    q_in: np.ndarray  # (Nx,) W/m^2
    h_e: np.ndarray  # (Nx,) W/(m^2 K)
    q_e: np.ndarray  # (Nx,) W/m^2, evaporative flux leaving the top face
    theta_e: np.ndarray  # (Nx,) top-face superheat T_e - T_v [K]
    theta_source: np.ndarray  # (Nx,) bottom-face superheat T_source - T_v [K]
    linear_residual: float  # ||A theta - b||_inf / max(||b||_inf, floor) [-]
    k_s: float


def solve_slab(mesh: SlabMesh, k_s: float, q_in, h_e) -> SlabSolution:
    """Solve the linear conduction problem for given q_in''(x) and h_e(x) on cells."""
    Nx, Nz, dz = mesh.Nx, mesh.Nz, mesh.dz
    dx = mesh.dx
    q_in = np.asarray(q_in, float)
    h_e = np.asarray(h_e, float)
    if q_in.shape != (Nx,) or h_e.shape != (Nx,):
        raise ValueError("q_in and h_e must have one value per x cell")
    if not (np.all(np.isfinite(q_in)) and np.all(np.isfinite(h_e))):
        raise ValueError("nonfinite q_in or h_e")
    if np.any(h_e <= 0):
        raise ValueError("h_e must be positive on every evaporating cell")
    if k_s <= 0:
        raise ValueError("k_s must be positive")

    n = Nx * Nz
    idx = np.arange(n).reshape(Nx, Nz)
    diag = np.zeros(n)
    b = np.zeros(n)
    rows, cols, vals = [], [], []

    # x-direction interior faces: conductance per unit width, face area dz
    Gx = k_s * dz / (0.5 * (dx[:-1] + dx[1:]))  # (Nx-1,)
    a, c = idx[:-1, :], idx[1:, :]
    G = np.broadcast_to(Gx[:, None], a.shape)
    np.add.at(diag, a.ravel(), G.ravel())
    np.add.at(diag, c.ravel(), G.ravel())
    rows += [a.ravel(), c.ravel()]
    cols += [c.ravel(), a.ravel()]
    vals += [-G.ravel(), -G.ravel()]
    # z-direction interior faces: face area dx[i]
    Gz = k_s * dx / dz
    a, c = idx[:, :-1], idx[:, 1:]
    G = np.broadcast_to(Gz[:, None], a.shape)
    np.add.at(diag, a.ravel(), G.ravel())
    np.add.at(diag, c.ravel(), G.ravel())
    rows += [a.ravel(), c.ravel()]
    cols += [c.ravel(), a.ravel()]
    vals += [-G.ravel(), -G.ravel()]
    G_top = dx / (1.0 / h_e + dz / (2 * k_s))
    b[idx[:, 0]] += q_in * dx
    diag[idx[:, -1]] += G_top

    A = sp.coo_matrix((np.concatenate(vals + [diag]),
                       (np.concatenate(rows + [np.arange(n)]), np.concatenate(cols + [np.arange(n)]))),
                      shape=(n, n)).tocsc()
    theta = spla.spsolve(A, b)
    res = A @ theta - b
    bmax = float(np.max(np.abs(b)))
    lin_res = float(np.max(np.abs(res)) / bmax) if bmax > 0 else float(np.max(np.abs(res)))

    T = theta.reshape(Nx, Nz)
    q_e = T[:, -1] / (1.0 / h_e + dz / (2 * k_s))
    theta_e = T[:, -1] - q_e * dz / (2 * k_s)
    theta_source = T[:, 0] + q_in * dz / (2 * k_s)
    return SlabSolution(mesh, T, q_in, h_e, q_e, theta_e, theta_source, lin_res, k_s)


def reconstruct_energy_balance(sol: SlabSolution) -> dict:
    """Independently rebuild every conductive face flux from the temperature field.

    Returns per-cell residuals (net inflow, W per unit width) and global terms. This does not reuse
    the assembled matrix; it recomputes fluxes from theta, q_in, h_e and the mesh.
    """
    m, k = sol.mesh, sol.k_s
    T, dx, dz = sol.theta, m.dx, m.dz
    Nx, Nz = T.shape
    net = np.zeros_like(T)
    # x faces (insulated ends contribute nothing)
    for i in range(Nx - 1):
        f = -k * (T[i + 1, :] - T[i, :]) / (0.5 * (dx[i] + dx[i + 1])) * dz  # flux in +x
        net[i, :] -= f
        net[i + 1, :] += f
    # z faces
    for j in range(Nz - 1):
        f = -k * (T[:, j + 1] - T[:, j]) / dz * dx  # flux in +z
        net[:, j] -= f
        net[:, j + 1] += f
    net[:, 0] += sol.q_in * dx
    top = (T[:, -1] / (1.0 / sol.h_e + dz / (2 * k))) * dx
    net[:, -1] -= top
    # independent check that the reported T_e satisfies the Robin condition at the face
    robin = sol.h_e * sol.theta_e - sol.q_e
    E_in = float(np.sum(sol.q_in * dx))
    E_out = float(np.sum(top))
    # scale: sum of magnitudes of the balance terms in each cell
    scale_local = np.abs(sol.q_in * dx).max() if np.any(sol.q_in) else np.abs(top).max()
    return {
        "per_cell_net_W_per_m": net,
        "max_local_abs_W_per_m": float(np.max(np.abs(net))),
        "local_scale_W_per_m": float(scale_local),
        "input_W_per_m": E_in,
        "evaporation_W_per_m": E_out,
        "global_abs_W_per_m": E_in - E_out,
        "robin_max_abs_W_m2": float(np.max(np.abs(robin))),
    }
