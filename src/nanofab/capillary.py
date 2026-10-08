"""Capillary pressure: the driving force that pulls liquid through the wick.

Equations
---------
Young-Laplace (general):
    dP = gamma * (1/R1 + 1/R2)

Pore of radius r with contact angle theta (cylindrical tube):
    dP = 2 gamma cos(theta) / r

Maximum capillary pressure of a square micropillar array (force balance on one unit cell at the
receding contact angle, i.e. when the meniscus is about to detach from the pillar sidewall):
    P_c,max = gamma cos(theta_r) * pi d / (l^2 - pi d^2 / 4)

    Contact-line force (pi d * gamma cos theta_r) balances pressure acting on the liquid area of
    the unit cell (l^2 - pi d^2/4).
    Source: vaartstra2019 Eq. (3), written there as 4 gamma cos(theta_r) / (d (4/pi (l/d)^2 - 1)),
    which is algebraically identical. The force balance is commonly credited to Xiao, Enright &
    Wang, Langmuir 26 (2010) 15070 (vaartstra2019 ref. [19]; not checked by us) and is used by
    Zhu et al., Langmuir 32 (2016) 1920.

Assumptions: static meniscus, uniform wettability, pillars tall enough that the meniscus does not
touch the substrate before reaching theta_r.
"""

from __future__ import annotations

import math

from .geometry import check_pillar_geometry


def young_laplace(surface_tension: float, R1: float, R2: float = math.inf) -> float:
    """Pressure jump across an interface with principal radii R1, R2 [m] -> [Pa].

    Use ``math.inf`` for a radius that is flat (e.g. a cylindrical interface).
    """
    curvature = (0.0 if math.isinf(R1) else 1.0 / R1) + (0.0 if math.isinf(R2) else 1.0 / R2)
    return surface_tension * curvature


def tube_capillary_pressure(surface_tension: float, pore_radius: float, contact_angle: float) -> float:
    """Capillary pressure of a cylindrical pore [Pa]. ``contact_angle`` in radians."""
    if pore_radius <= 0:
        raise ValueError("pore_radius must be positive")
    return 2 * surface_tension * math.cos(contact_angle) / pore_radius


def pillar_max_capillary_pressure(
    surface_tension: float, d: float, l: float, receding_contact_angle: float
) -> float:
    """Maximum capillary pressure of a square pillar array [Pa] (Xiao et al. 2010 force balance).

    Args:
        surface_tension: liquid-vapor surface tension gamma [N/m]
        d: pillar diameter [m]
        l: pillar pitch [m]
        receding_contact_angle: theta_r [rad]
    """
    check_pillar_geometry(d, l)
    liquid_area = l**2 - math.pi * d**2 / 4
    return surface_tension * math.cos(receding_contact_angle) * math.pi * d / liquid_area


def meniscus_radius(surface_tension: float, capillary_pressure: float) -> float:
    """Radius [m] of a spherical-cap meniscus (R1 = R2 = R) that sustains ``capillary_pressure``.

    R = 2 gamma / P_c. This is the Stage 3 COMSOL approximation; the real meniscus in a pillar array
    is not spherical (Vaartstra et al. 2019 solve the normal stress balance, their Eq. (2), as a
    PDE in COMSOL; see Sec. 2.1).
    """
    if capillary_pressure <= 0:
        return math.inf
    return 2 * surface_tension / capillary_pressure


def spherical_cap_depth(radius: float, half_span: float) -> float:
    """Sag [m] of a spherical cap of ``radius`` spanning a chord of half-width ``half_span``.

    Useful for checking whether a meniscus of a given curvature fits within the pillar height.
    """
    if math.isinf(radius):
        return 0.0
    if half_span > radius:
        raise ValueError("Meniscus radius smaller than half-span: interface cannot span the gap")
    return radius - math.sqrt(radius**2 - half_span**2)
