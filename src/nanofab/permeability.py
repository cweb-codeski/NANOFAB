"""Permeability K of a micropillar wick for flow parallel to the substrate.

Darcy's law (superficial velocity U, averaged over the full wick cross-section):
    U = -(K / mu) dP/dx

Step 1: in-plane permeability of an infinitely tall square array of cylinders (transverse flow),
Sangani & Acrivos, Int. J. Multiphase Flow 8 (1982) 193, written in the Drummond & Tahir (1984)
form with a = d/2:
    K_s = a^2 / (8 phi) * (-ln(phi) - 1.476 + 2 phi - 1.774 phi^2 + 4.076 phi^3)
Accurate for phi below ~0.4 (dilute to moderately dense arrays).

Step 2: correction for the finite height h (no-slip substrate at z = 0, shear-free flat liquid
surface at z = h), from the Brinkman equation with mu_eff = mu:
    mu d^2u/dz^2 - (mu / K_s) u = dP/dx
    u(z) = -(K_s/mu) dP/dx * [1 - cosh(beta (h - z)) / cosh(beta h)],  beta = 1/sqrt(K_s)
Averaging over height:
    K_eff = K_s * [1 - tanh(beta h) / (beta h)]
Source: Brinkman, Appl. Sci. Res. 1 (1949) 155 (vaartstra2019 ref. [20]); the uniform-K Brinkman
approach for pillar wicks is attributed by vaartstra2019 to Adera et al., IJHMT 101 (2016) 280.
Neither has been checked by us (TODO). The height-averaged closed form above is our derivation.
A citation to "Byon et al., IJHMT 79 (2014)" appeared here previously; it was written from memory
and could not be verified, so it was removed (see docs/SOURCE_AUTHORITY.md, audit).

Assumptions: Stokes flow (Re << 1), flat liquid surface (meniscus curvature ignored, which
overestimates K when the meniscus recedes into the pillars), fully wetted pillars.
"""

from __future__ import annotations

import math
import warnings

from .geometry import check_pillar_geometry, solid_fraction

SANGANI_ACRIVOS_PHI_MAX = 0.4


def sangani_acrivos_permeability(d: float, l: float) -> float:
    """In-plane permeability K_s [m^2] of an infinite square pillar array (no floor, no top)."""
    phi = solid_fraction(d, l)
    if phi > SANGANI_ACRIVOS_PHI_MAX:
        warnings.warn(
            f"Solid fraction phi={phi:.3f} > {SANGANI_ACRIVOS_PHI_MAX}: Sangani-Acrivos dilute "
            "expansion is outside its accuracy range",
            stacklevel=2,
        )
    a = d / 2
    series = -math.log(phi) - 1.476 + 2 * phi - 1.774 * phi**2 + 4.076 * phi**3
    return a**2 / (8 * phi) * series


def brinkman_height_factor(K_s: float, h: float) -> float:
    """Finite-height correction 1 - tanh(beta h)/(beta h), beta = 1/sqrt(K_s) [-].

    Tends to 1 for tall pillars (h >> sqrt(K_s)) and to (beta h)^2 / 3 for short pillars, which
    recovers the plane-Poiseuille film result K = h^2/3.
    """
    bh = h / math.sqrt(K_s)
    if bh < 1e-4:
        return bh**2 / 3
    return 1.0 - math.tanh(bh) / bh


def pillar_wick_permeability(d: float, l: float, h: float) -> float:
    """Effective permeability K_eff [m^2] of a pillar wick of height h (Brinkman-corrected)."""
    check_pillar_geometry(d, l, h)
    K_s = sangani_acrivos_permeability(d, l)
    return K_s * brinkman_height_factor(K_s, h)


def brinkman_velocity_profile(z: float, K_s: float, h: float, dpdx: float, mu: float) -> float:
    """Height-resolved Darcy-Brinkman velocity u(z) [m/s] (useful to compare with COMSOL)."""
    beta = 1 / math.sqrt(K_s)
    return -(K_s / mu) * dpdx * (1 - math.cosh(beta * (h - z)) / math.cosh(beta * h))


def permeability_from_flow(
    volumetric_flow: float, cross_section_area: float, length: float, pressure_drop: float, mu: float
) -> float:
    """Extract K [m^2] from a simulated or measured flow: K = mu U L / dP, U = Q / A.

    Args:
        volumetric_flow: Q through the outlet face [m^3/s]
        cross_section_area: full cross-section A normal to the flow, including pillars [m^2]
        length: flow length L over which pressure_drop is applied [m]
        pressure_drop: dP = p_in - p_out > 0 [Pa]
        mu: dynamic viscosity [Pa s]
    """
    if pressure_drop <= 0:
        raise ValueError("pressure_drop must be positive (p_in - p_out)")
    superficial_velocity = volumetric_flow / cross_section_area
    return mu * superficial_velocity * length / pressure_drop
