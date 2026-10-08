"""1D capillary-limited dryout heat flux of a wick fed from one edge.

Model
-----
A wick of height h and length L is fed with liquid at x = 0 and closed at x = L. Heat flux q''
(per unit footprint area) is uniform and is all removed by evaporation. Mass conservation per unit
width with Darcy flow:

    d/dx (rho U h) = -q'' / h_lv,    U = -(K/mu) dP/dx

Boundary conditions: P(0) = P_feed, dP/dx(L) = 0 (no flow out of the closed end). Solution:

    P(x) = P_feed - (mu q'' / (rho K h h_lv)) * (L x - x^2 / 2)
    dP_total = P(0) - P(L) = mu q'' L^2 / (2 rho K h h_lv)

Dryout occurs when dP_total equals the capillary pressure available (P_c,max minus any adverse
pressure such as gravity or a negative feed pressure):

    q''_dry = 2 rho K h h_lv dP_available / (mu L^2)

Source: standard capillary-limit analysis (e.g. Faghri, Heat Pipe Science and Technology; not
checked). vaartstra2019 Sec. 1 reports that Adera et al., IJHMT 101 (2016) 280 used uniform K
with uniform evaporation, and Zhu et al., Langmuir 32 (2016) 1920 solved uniform evaporation
numerically. The closed form above is our derivation (MODEL_CALCULATED when evaluated).

Assumptions: uniform heat flux and wick properties, no meniscus recession effect on K, 1D flow, all
heat goes to evaporation (no sensible heating of the liquid), steady state. For a wick fed from
both edges, use L = half the wick length.
"""

from __future__ import annotations

import numpy as np


def dryout_heat_flux(
    permeability: float,
    wick_height: float,
    feed_length: float,
    available_pressure: float,
    rho_l: float,
    mu_l: float,
    h_lv: float,
) -> float:
    """Capillary-limited dryout heat flux q''_dry [W/m^2] per unit footprint area."""
    if feed_length <= 0 or wick_height <= 0:
        raise ValueError("feed_length and wick_height must be positive")
    return 2 * rho_l * permeability * wick_height * h_lv * available_pressure / (mu_l * feed_length**2)


def wick_pressure_drop(
    heat_flux: float,
    permeability: float,
    wick_height: float,
    feed_length: float,
    rho_l: float,
    mu_l: float,
    h_lv: float,
) -> float:
    """Liquid pressure drop [Pa] from the feed to the closed end at heat flux q'' [W/m^2]."""
    return mu_l * heat_flux * feed_length**2 / (2 * rho_l * permeability * wick_height * h_lv)


def pressure_profile(
    x: np.ndarray,
    heat_flux: float,
    permeability: float,
    wick_height: float,
    feed_length: float,
    rho_l: float,
    mu_l: float,
    h_lv: float,
    feed_pressure: float = 0.0,
) -> np.ndarray:
    """Liquid pressure P(x) [Pa] along the wick (x = 0 at the feed)."""
    x = np.asarray(x, dtype=float)
    coeff = mu_l * heat_flux / (rho_l * permeability * wick_height * h_lv)
    return feed_pressure - coeff * (feed_length * x - x**2 / 2)


def superficial_velocity_at_feed(heat_flux: float, feed_length: float, wick_height: float, rho_l: float, h_lv: float) -> float:
    """Darcy (superficial) velocity [m/s] entering the wick at x = 0; the maximum along the wick."""
    return heat_flux * feed_length / (rho_l * h_lv * wick_height)


def gravity_head(rho_l: float, gravity: float, elevation: float) -> float:
    """Hydrostatic pressure [Pa] the wick must overcome to lift liquid by ``elevation`` [m]."""
    return rho_l * gravity * elevation
