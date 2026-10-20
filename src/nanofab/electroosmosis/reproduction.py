"""Tools for reproducing pujahari2022 and for the dilute-KCl electrolyte bookkeeping.

Inverse dryout gain
    Pi = G^-1(q''_B / q''_A): the EO pressure ratio a published dryout pair implies, so a published
    curve can be compared with this model one factor at a time.

Implied permeability
    K_implied = q''_A mu L^2 / (2 rho h h_fg Pc_max): the uniform K that makes the 1D capillary limit
    return a published passive dryout flux (absorbs meniscus-curvature effects).

Uniform-imposed-field EO estimate (interpretation of pujahari2022, PROJECT_DERIVED)
    pujahari2022 Sec. 2.1.2 writes the external potential as a 1D Laplace equation and applies "a
    constant electric field ... along each unit cell" in the rho_e E body force. If E is a uniform
    x-field, the slip it drives on an insulating pillar sidewall is mu_EO E sin(theta), whereas the
    field that respects n.E = 0 at an insulating wall gives ~2 mu_EO E sin(theta) (dilute limit). The
    floor still sees the full field. Splitting the flow into floor-driven and pillar-driven Brinkman
    parts (beta = 1/sqrt(K_s), f = tanh(beta h)/(beta h)) gives the ratio to mu_EO E / F:
        R = U_EO,uniform-E / (mu_EO E / F) ~= (1 + phi) f + (1 - f) / 2
    R -> 1/2 for tall dense pillars (pillar walls dominate) and -> 1 + phi as h -> 0 (floor only).
    An estimate to test against the paper, not a replacement model.

Advective salt enrichment (dilute electrolyte, solvent evaporates, ions stay in the liquid)
    With Pe = U L / D >> 1 and no ion loss to the vapor, the salt flux Q(x) c(x) is conserved:
        c(x) / c_feed = Q_0 / Q(x) = 1 / (1 - x/L)     (uniform evaporation, both cases A and B)
    valid where the advective length (L - x) exceeds the diffusion length sqrt(D eps L / U_0).
    Inside that end layer salt accumulates with time (no steady state without a sink).
"""

from __future__ import annotations

import math

from scipy.optimize import brentq

from ..geometry import solid_fraction
from ..permeability import sangani_acrivos_permeability
from .wick1d import dryout_gain


def inverse_dryout_gain(G: float) -> float:
    """Pi such that dryout_gain(Pi) = G (G >= 1)."""
    if G < 1:
        return G - 1.0  # opposing branch, G = 1 + Pi
    if G == 1:
        return 0.0
    hi = 2 * G + 2
    return brentq(lambda P: float(dryout_gain(P)) - G, 0.0, hi)


def implied_permeability(q_dry_passive: float, mu: float, L: float, rho: float, h: float, h_fg: float,
                         Pc_max: float) -> float:
    """Uniform K [m^2] that reproduces a passive dryout flux with the 1D capillary limit."""
    return q_dry_passive * mu * L**2 / (2 * rho * h * h_fg * Pc_max)


def floor_fraction(d: float, l: float, h: float) -> float:
    """f = tanh(beta h)/(beta h), beta = 1/sqrt(K_s): depth-averaged reach of floor-driven slip."""
    bh = h / math.sqrt(sangani_acrivos_permeability(d, l))
    return math.tanh(bh) / bh


def uniform_field_eo_ratio(d: float, l: float, h: float) -> float:
    """R = U_EO under a uniform imposed x-field / (mu_EO E / F) (estimate; see module docstring)."""
    f = floor_fraction(d, l, h)
    return (1 + solid_fraction(d, l)) * f + (1 - f) / 2


def advective_enrichment(x, L: float):
    """c(x)/c_feed = 1/(1 - x/L) for uniform evaporation with Pe >> 1."""
    return 1.0 / (1.0 - x / L)


def enrichment_layer_width(D: float, porosity: float, L: float, U0: float) -> float:
    """Distance from the closed end inside which diffusion matters, sqrt(D eps L / U_0) [m]."""
    return math.sqrt(D * porosity * L / U0)
