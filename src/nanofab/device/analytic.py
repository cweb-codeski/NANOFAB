"""Closed-form verification limits for the passive half-device (PROJECT_DERIVED).

These are derived by us from V2019 Eqs. (8)-(9) and the boundary conditions of V2019 Sec. 2.3
(Pc(0) = Pc_feed, Q(L) = 0, insulated substrate ends). They are verification references for the
numerical solver, not V2019 results.

Constant K, constant h_e, uniform q_in'' = q_e'' = q, insulated ends, saturated inlet:

    Q(x)      = W q (L - x) / (rho h_fg)
    p_g(x)    = p_g(0) - mu q (L x - x^2/2) / (rho K h h_fg)
    q_dry     = 2 rho K h h_fg (Pc_limit - Pc_feed) / (mu L^2)
    T_e       = T_v + q / h_e
    T_source  = T_e + q t_s / k_s

Pressure-dependent transmissivity tau(Pc) = h K_bulk(Pc) [m^3], prescribed uniform evaporation:

    q_dry = (2 rho h_fg / L^2) * integral_{Pc_feed}^{Pc_limit} tau(Pc)/mu dPc

General evaporation-moment identity (any q_e''(x) >= 0, Q(L) = 0):

    integral_{Pc(0)}^{Pc(L)} tau/mu dPc = integral_0^L x q_e''(x) / (rho h_fg) dx

The uniform-q closed forms are local identities only because the constant-property thermal
problem is one-dimensional; with variable h_e, q_e'' must come from the conduction solution.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy import integrate


def flow_profile(x, q, L, W, rho, h_fg):
    """Q(x) [m^3/s] for uniform evaporation q [W/m^2]."""
    return W * q * (L - np.asarray(x, float)) / (rho * h_fg)


def pressure_profile_pg(x, q, L, K, h, rho, mu, h_fg, pg0=0.0):
    """p_g(x) = p_liquid - p_v [Pa] for constant K and uniform evaporation."""
    x = np.asarray(x, float)
    return pg0 - mu * q * (L * x - x**2 / 2) / (rho * K * h * h_fg)


def dryout_constant(K, h, L, rho, mu, h_fg, Pc_limit, Pc_feed=0.0):
    """q_dry [W/m^2] for constant K and uniform evaporation."""
    return 2 * rho * K * h * h_fg * (Pc_limit - Pc_feed) / (mu * L**2)


def superheats_constant(q, h_e, t_s, k_s):
    """(T_e - T_v, T_source - T_v) [K] for the uniform constant-property limit."""
    dTe = q / h_e
    return dTe, dTe + q * t_s / k_s


def dryout_from_transmissivity(
    tau: Callable[[float], float], mu: float, rho: float, h_fg: float, L: float,
    Pc_limit: float, Pc_feed: float = 0.0, breakpoints=None,
) -> float:
    """q_dry for uniform evaporation and pressure-dependent tau(Pc), by adaptive quadrature."""
    pts = None
    if breakpoints is not None:
        pts = [p for p in breakpoints if Pc_feed < p < Pc_limit]
    val, _ = integrate.quad(lambda p: tau(p) / mu, Pc_feed, Pc_limit, points=pts, limit=500,
                            epsabs=0.0, epsrel=1e-13)
    return 2 * rho * h_fg * val / L**2


def evaporation_moment(x_faces, q_cells, rho, h_fg) -> float:
    """integral_0^L x q_e''/(rho h_fg) dx for piecewise-constant q_e'' [m^2/s]."""
    x_faces = np.asarray(x_faces, float)
    xc = 0.5 * (x_faces[1:] + x_faces[:-1])
    dx = np.diff(x_faces)
    return float(np.sum(xc * np.asarray(q_cells, float) * dx) / (rho * h_fg))
