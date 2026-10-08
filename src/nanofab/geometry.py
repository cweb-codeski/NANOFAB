"""Geometric descriptors of a square array of cylindrical micropillars.

Notation (consistent with the Wang-group micropillar papers):
    d : pillar diameter [m]
    l : pillar pitch, center-to-center [m]
    h : pillar height [m]
"""

from __future__ import annotations

import math


def check_pillar_geometry(d: float, l: float, h: float | None = None) -> None:
    """Raise if the pillar dimensions are non-physical (pillars overlapping or non-positive)."""
    if d <= 0 or l <= 0:
        raise ValueError(f"Pillar diameter and pitch must be positive (d={d}, l={l})")
    if d >= l:
        raise ValueError(f"Pillars overlap: diameter d={d} must be smaller than pitch l={l}")
    if h is not None and h <= 0:
        raise ValueError(f"Pillar height must be positive (h={h})")


def solid_fraction(d: float, l: float) -> float:
    """Areal solid fraction phi = pi d^2 / (4 l^2) of a square pillar array [-]."""
    check_pillar_geometry(d, l)
    return math.pi * d**2 / (4 * l**2)


def porosity(d: float, l: float) -> float:
    """Areal (and, for straight pillars, volumetric) porosity 1 - phi [-]."""
    return 1.0 - solid_fraction(d, l)


def gap(d: float, l: float) -> float:
    """Minimum edge-to-edge spacing between neighboring pillars, l - d [m]."""
    check_pillar_geometry(d, l)
    return l - d


def hydraulic_diameter(d: float, l: float, h: float) -> float:
    """Hydraulic diameter of the pore space, 4 * (fluid volume) / (wetted area) per unit cell [m].

    Wetted area counts the pillar sidewall and the substrate floor (the free surface on top is not
    wetted). Used only for order-of-magnitude Reynolds numbers.
    """
    check_pillar_geometry(d, l, h)
    fluid_volume = (l**2 - math.pi * d**2 / 4) * h
    wetted_area = math.pi * d * h + (l**2 - math.pi * d**2 / 4)
    return 4 * fluid_volume / wetted_area


def pillar_ratio(d: float, l: float) -> float:
    """Diameter-to-pitch ratio d/l [-]."""
    check_pillar_geometry(d, l)
    return d / l
