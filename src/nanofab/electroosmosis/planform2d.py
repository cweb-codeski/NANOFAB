"""2D planform unit cell of a square pillar array: Laplace field and thin-EDL electroosmotic Stokes flow.

Purpose: isolate one mechanism, the electric-field treatment around insulating pillars.
    Case A (uniform):  E = E0 x everywhere in the liquid (pujahari2022-style imposed field)
    Case B (solved):   div(grad phi) = 0, n.grad(phi) = 0 on the pillar, mean gradient -E0
Both cases use the same Helmholtz-Smoluchowski slip on the pillar wall, u_t = mu_EO E_t, the same
geometry and the same mean field. Nondimensional results: superficial velocity in units of mu_EO E0.

LIMITATION (by construction): the planform model is z-invariant. It represents pillars of infinite
height with no floor; pillar height h does NOT enter. Floor slip and the curved meniscus of a real
wick are absent.

Method: boundary collocation (multipole series) in the quarter cell 0 <= x, y <= l/2, lengths in units
of l/2, pillar radius alpha = d/l. Square-array symmetry: x = l/2 is a fore-aft symmetry plane
(equipotential for phi; v = 0 and zero disturbance pressure, i.e. no mean pressure gradient, for the
flow); y = l/2 is a mirror plane (no current; v = 0 and zero shear).

    phi / (E0 l/2) = - sum_k c_k (rho^n + alpha^(2n) rho^-n) cos(n theta),   n = 1, 3, 5, ...
        (insulating pillar condition satisfied term by term)
    psi / (U* l/2) = sum_k [a1 rho^n + a2 (rho/alpha)^-n + a3 rho^(n+2) + a4 g4_n(rho)] sin(n theta)
        g4_n = (rho/alpha)^(2-n) (n >= 3), rho ln rho (n = 1); u = d psi/dy, v = -d psi/dx
        pillar: psi = 0 and -d psi/d rho = u_theta (slip) imposed exactly per mode
    pressure p l/2 / (mu U*): a3 4(n+1) rho^n cos + a4 4(n-1) alpha^-2 (rho/alpha)^-n cos (n >= 3),
        -2 a4 cos(theta)/rho (n = 1); harmonic terms carry no pressure.
The pillar-face conditions are exact; the cell-face conditions are met in least squares and their
residuals are reported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla


def _edge_points(m: int) -> np.ndarray:
    """Chebyshev-Lobatto points on [0, 1] (clustered at the ends)."""
    return 0.5 * (1 - np.cos(np.pi * np.arange(m) / (m - 1)))


def _edges(m: int):
    s = _edge_points(m)
    right = (np.ones(m), s)  # x = 1, y in [0, 1]
    top = (s, np.ones(m))    # y = 1, x in [0, 1]
    return right, top


def _polar(x, y):
    return np.hypot(x, y), np.arctan2(y, x)


# --------------------------------------------------------------------------------------------- Laplace

@dataclass
class FieldSolution:
    alpha: float
    n: np.ndarray
    c: np.ndarray
    residual_equipotential: float
    residual_mirror: float

    def phi(self, x, y):
        r, t = _polar(x, y)
        with np.errstate(divide="ignore", invalid="ignore"):
            H = r[..., None] ** self.n + self.alpha ** (2 * self.n) * r[..., None] ** (-self.n)
            return -np.sum(self.c * H * np.cos(self.n * t[..., None]), axis=-1)

    def grad(self, x, y):
        """(dphi/dx, dphi/dy) in units of E0."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return self._grad(x, y)

    def _grad(self, x, y):
        r, t = _polar(x, y)
        n, a = self.n, self.alpha
        R = r[..., None]
        H = R**n + a ** (2 * n) * R ** (-n)
        Hp = n * R ** (n - 1) - n * a ** (2 * n) * R ** (-n - 1)
        cn, sn = np.cos(n * t[..., None]), np.sin(n * t[..., None])
        dr = -np.sum(self.c * Hp * cn, axis=-1)
        dt = np.sum(self.c * H * n * sn, axis=-1)
        ct, st = np.cos(t), np.sin(t)
        return ct * dr - st / r * dt, st * dr + ct / r * dt

    def wall_tangential_field_modes(self) -> np.ndarray:
        """Sine coefficients e_n of E_theta(alpha, theta)/E0 = sum e_n sin(n theta)."""
        return -2 * self.n * self.c * self.alpha ** (self.n - 1)

    def current(self, x0: float, m: int = 400) -> float:
        """Current through the section x = x0 (liquid part), units sigma E0 (l/2): equals 1/F."""
        y_lo = np.sqrt(max(self.alpha**2 - x0**2, 0.0))
        g, w = np.polynomial.legendre.leggauss(m)
        y = y_lo + (1 - y_lo) * (g + 1) / 2
        ex, _ = self.grad(np.full_like(y, x0), y)
        return float(np.sum(-ex * w) * (1 - y_lo) / 2)

    @property
    def inverse_formation_factor(self) -> float:
        return self.current(1.0)


def solve_field(alpha: float, n_modes: int, m_points: int | None = None) -> FieldSolution:
    if not 0 < alpha < 1:
        raise ValueError("need 0 < d/l < 1")
    n = 2 * np.arange(n_modes) + 1
    m = m_points or 4 * n_modes
    (xr, yr), (xt, yt) = _edges(m)
    rows, rhs = [], []
    # right edge: phi = -1
    r, t = _polar(xr, yr)
    H = r[:, None] ** n + alpha ** (2 * n) * r[:, None] ** (-n)
    rows.append(-H * np.cos(n * t[:, None]))
    rhs.append(-np.ones(m))
    # top edge: dphi/dy = 0
    r, t = _polar(xt, yt)
    R = r[:, None]
    H = R**n + alpha ** (2 * n) * R ** (-n)
    Hp = n * R ** (n - 1) - n * alpha ** (2 * n) * R ** (-n - 1)
    dr = -Hp * np.cos(n * t[:, None])
    dt = H * n * np.sin(n * t[:, None])
    rows.append(np.sin(t)[:, None] * dr + (np.cos(t) / r)[:, None] * dt)
    rhs.append(np.zeros(m))
    A, b = np.vstack(rows), np.concatenate(rhs)
    scale = np.abs(A).max(axis=0)
    c = sla.lstsq(A / scale, b)[0] / scale
    sol = FieldSolution(alpha=alpha, n=n, c=c, residual_equipotential=0.0, residual_mirror=0.0)
    s = np.linspace(0, 1, 1001)
    sol.residual_equipotential = float(np.max(np.abs(sol.phi(np.ones_like(s), s) + 1)))
    sol.residual_mirror = float(np.max(np.abs(sol.grad(s, np.ones_like(s))[1])))
    return sol


# --------------------------------------------------------------------------------------------- Stokes

def _mode_funcs(n: int, alpha: float, rho):
    """Values, d/drho, Laplacian-radial-part and pressure-radial-part of the 4 basis functions of mode n.

    Laplacian of f(rho) sin(n theta) is L(rho) sin(n theta); pressure is P(rho) cos(n theta).
    """
    rho = np.asarray(rho, dtype=float)
    q = rho / alpha
    f = np.empty(rho.shape + (4,))
    fp = np.empty_like(f)
    L = np.zeros_like(f)
    P = np.zeros_like(f)
    f[..., 0], fp[..., 0] = rho**n, n * rho ** (n - 1)
    f[..., 1], fp[..., 1] = q ** (-n), -n * q ** (-n) / rho
    f[..., 2], fp[..., 2] = rho ** (n + 2), (n + 2) * rho ** (n + 1)
    L[..., 2], P[..., 2] = 4 * (n + 1) * rho**n, 4 * (n + 1) * rho**n
    if n == 1:
        f[..., 3], fp[..., 3] = rho * np.log(rho), np.log(rho) + 1
        L[..., 3], P[..., 3] = 2 / rho, -2 / rho
    else:
        f[..., 3], fp[..., 3] = q ** (2 - n), (2 - n) * q ** (2 - n) / rho
        L[..., 3] = 4 * (1 - n) * alpha**-2 * q ** (-n)
        P[..., 3] = 4 * (n - 1) * alpha**-2 * q ** (-n)
    return f, fp, L, P


@dataclass
class FlowSolution:
    alpha: float
    n: np.ndarray
    coef: np.ndarray  # (n_modes, 4)
    q: float  # superficial velocity / U*
    residuals: dict

    def _eval(self, x, y):
        with np.errstate(divide="ignore", invalid="ignore"):  # points inside the pillar are not physical
            return self._eval_raw(x, y)

    def _eval_raw(self, x, y):
        r, t = _polar(np.asarray(x, float), np.asarray(y, float))
        out = {k: np.zeros_like(r) for k in ("psi", "psi_r", "psi_t", "lap", "p")}
        for k, n in enumerate(self.n):
            f, fp, L, P = _mode_funcs(int(n), self.alpha, r)
            c = self.coef[k]
            sn, cn = np.sin(n * t), np.cos(n * t)
            out["psi"] += (f @ c) * sn
            out["psi_r"] += (fp @ c) * sn
            out["psi_t"] += (f @ c) * n * cn
            out["lap"] += (L @ c) * sn
            out["p"] += (P @ c) * cn
        return r, t, out

    def velocity(self, x, y):
        r, t, o = self._eval(x, y)
        with np.errstate(divide="ignore", invalid="ignore"):
            u = np.sin(t) * o["psi_r"] + np.cos(t) / r * o["psi_t"]
            v = -(np.cos(t) * o["psi_r"] - np.sin(t) / r * o["psi_t"])
        return u, v

    def psi(self, x, y):
        return self._eval(x, y)[2]["psi"]

    def pressure(self, x, y):
        return self._eval(x, y)[2]["p"]

    def section_flux(self, x0: float, m: int = 400) -> float:
        """Integral of u over the liquid part of x = x0, 0 <= y <= 1 (equals q by mass conservation)."""
        y_lo = np.sqrt(max(self.alpha**2 - x0**2, 0.0))
        g, w = np.polynomial.legendre.leggauss(m)
        y = y_lo + (1 - y_lo) * (g + 1) / 2
        u, _ = self.velocity(np.full_like(y, x0), y)
        return float(np.sum(u * w) * (1 - y_lo) / 2)


def solve_flow(alpha: float, n_modes: int, slip_modes=None, pressure_drop: float = 0.0,
               m_points: int | None = None) -> FlowSolution:
    """Stokes flow in the quarter cell.

    slip_modes: sine coefficients s_n of the wall tangential velocity u_theta(alpha, theta)/U*
                (length <= n_modes; None = no slip, i.e. the classical no-slip wall).
    pressure_drop: nondimensional pressure at x = 1 is -pressure_drop (0 = no mean gradient).
                   With U* = G (l/2)^2 / mu and pressure_drop = 1, q = K / (l/2)^2.
    """
    n = 2 * np.arange(n_modes) + 1
    s = np.zeros(n_modes)
    if slip_modes is not None:
        sm = np.asarray(slip_modes, float)
        s[: min(len(sm), n_modes)] = sm[:n_modes]
    m = m_points or 4 * n_modes
    nu = 4 * n_modes + 1  # coefficients + q
    # exact pillar conditions, mode by mode: f(alpha) = 0, f'(alpha) = -s_n
    Cg = np.zeros((2 * n_modes, nu))
    dg = np.zeros(2 * n_modes)
    for k, nk in enumerate(n):
        f, fp, _, _ = _mode_funcs(int(nk), alpha, np.array([alpha]))
        Cg[2 * k, 4 * k:4 * k + 4] = f[0]
        Cg[2 * k + 1, 4 * k:4 * k + 4] = fp[0]
        dg[2 * k + 1] = -s[k]
    x_p = sla.lstsq(Cg, dg)[0]
    N = sla.null_space(Cg)

    (xr, yr), (xt, yt) = _edges(m)
    rows, rhs = [], []

    def basis_rows(x, y, kind):
        r, t = _polar(x, y)
        M = np.zeros((x.size, nu))
        for k, nk in enumerate(n):
            f, fp, L, P = _mode_funcs(int(nk), alpha, r)
            sn, cn = np.sin(nk * t)[:, None], np.cos(nk * t)[:, None]
            if kind == "psi":
                M[:, 4 * k:4 * k + 4] = f * sn
            elif kind == "psi_x":
                M[:, 4 * k:4 * k + 4] = np.cos(t)[:, None] * fp * sn - (np.sin(t) / r)[:, None] * f * nk * cn
            elif kind == "lap":
                M[:, 4 * k:4 * k + 4] = L * sn
            elif kind == "p":
                M[:, 4 * k:4 * k + 4] = P * cn
        return M

    # right edge: v = 0 (psi_x = 0) and p = -pressure_drop
    rows.append(basis_rows(xr, yr, "psi_x"))
    rhs.append(np.zeros(m))
    rows.append(basis_rows(xr, yr, "p"))
    rhs.append(np.full(m, -pressure_drop))
    # top edge: psi = q and zero vorticity (Laplacian psi = 0)
    Mt = basis_rows(xt, yt, "psi")
    Mt[:, -1] = -1.0
    rows.append(Mt)
    rhs.append(np.zeros(m))
    rows.append(basis_rows(xt, yt, "lap"))
    rhs.append(np.zeros(m))
    A, b = np.vstack(rows), np.concatenate(rhs)
    AN = A @ N
    scale = np.abs(AN).max(axis=0)
    scale[scale == 0] = 1.0
    z = sla.lstsq(AN / scale, b - A @ x_p)[0] / scale
    x = x_p + N @ z
    sol = FlowSolution(alpha=alpha, n=n, coef=x[:-1].reshape(n_modes, 4), q=float(x[-1]), residuals={})

    sp = np.linspace(0, 1, 801)
    one = np.ones_like(sp)
    u_r, v_r = sol.velocity(one, sp)
    u_t, v_t = sol.velocity(sp, one)
    th = np.linspace(0, np.pi / 2, 401)
    xc, yc = alpha * np.cos(th), alpha * np.sin(th)
    u_c, v_c = sol.velocity(xc * (1 + 1e-12), yc * (1 + 1e-12))
    un = u_c * np.cos(th) + v_c * np.sin(th)
    ut = -u_c * np.sin(th) + v_c * np.cos(th)
    s_target = sum(s[k] * np.sin(nk * th) for k, nk in enumerate(n))
    scale_v = max(abs(sol.q), np.max(np.abs(s)), 1e-30)
    sol.residuals = {
        "right_v": float(np.max(np.abs(v_r))) / scale_v,
        "right_p": float(np.max(np.abs(sol.pressure(one, sp) + pressure_drop))),
        "top_psi": float(np.max(np.abs(sol.psi(sp, one) - sol.q))) / scale_v,
        "top_v": float(np.max(np.abs(v_t))) / scale_v,
        "wall_normal": float(np.max(np.abs(un))) / scale_v,
        "wall_slip": float(np.max(np.abs(ut - s_target))) / scale_v,
    }
    return sol


# --------------------------------------------------------------------------------------------- cases

@dataclass
class FieldTreatmentResult:
    alpha: float
    n_modes: int
    inv_F: float  # 1/F from the solved field (current through x = l/2)
    q_solved: float  # U_EO / (mu_EO E0), Case B
    q_uniform: float  # U_EO / (mu_EO E0), Case A
    field: FieldSolution
    flow_solved: FlowSolution
    flow_uniform: FlowSolution

    @property
    def R_field(self) -> float:
        """User definition: Q_EO(solved field) / Q_EO(uniform field)."""
        return self.q_solved / self.q_uniform

    @property
    def uniform_over_solved(self) -> float:
        """Q_EO(uniform) / Q_EO(solved): the quantity comparable to R_inferred from pujahari2022."""
        return self.q_uniform / self.q_solved


def field_treatment(d: float, l: float, n_modes: int = 24, E_sign: float = 1.0) -> FieldTreatmentResult:
    """Case A (uniform E0 x) and Case B (solved field) on the same 2D cell, same mean field E0."""
    alpha = d / l
    fs = solve_field(alpha, n_modes)
    e_solved = E_sign * fs.wall_tangential_field_modes()  # u_theta / (mu_EO E0) on the wall
    e_uniform = np.zeros(n_modes)
    e_uniform[0] = -E_sign  # E0 x . theta_hat = -E0 sin(theta)
    flow_b = solve_flow(alpha, n_modes, slip_modes=e_solved)
    flow_a = solve_flow(alpha, n_modes, slip_modes=e_uniform)
    return FieldTreatmentResult(alpha=alpha, n_modes=n_modes, inv_F=fs.inverse_formation_factor,
                                q_solved=flow_b.q, q_uniform=flow_a.q, field=fs,
                                flow_solved=flow_b, flow_uniform=flow_a)


def permeability_2d(d: float, l: float, n_modes: int = 24) -> float:
    """Transverse permeability [m^2] of the infinite square array from the same Stokes solver."""
    sol = solve_flow(d / l, n_modes, slip_modes=None, pressure_drop=1.0)
    return sol.q * (l / 2) ** 2
