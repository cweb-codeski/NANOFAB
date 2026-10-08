"""Dimensionless groups used to justify modeling assumptions.

    Re = rho U L / mu              inertia / viscous      Re << 1 -> Stokes (creeping) flow OK
    Ca = mu U / gamma              viscous / capillary    Ca << 1 -> meniscus shape set by statics
    Bo = (rho_l - rho_v) g L^2 / gamma   gravity / capillary   Bo << 1 -> gravity negligible
    Kn = lambda / L                mean free path / length     Kn << 0.01 -> continuum vapor
    Pe = U L / alpha               advection / conduction      Pe << 1 -> liquid flow can be
                                                               ignored in the heat equation

Vapor mean free path (hard-sphere kinetic theory):
    lambda = k_B T / (sqrt(2) pi d_m^2 p)
Source: Bird, Molecular Gas Dynamics (1976); Incropera et al.
"""

from __future__ import annotations

import math

from .fluids import BOLTZMANN_CONSTANT


def reynolds(rho: float, velocity: float, length: float, mu: float) -> float:
    return rho * velocity * length / mu


def capillary_number(mu: float, velocity: float, surface_tension: float) -> float:
    return mu * velocity / surface_tension


def peclet(velocity: float, length: float, thermal_diffusivity: float) -> float:
    return velocity * length / thermal_diffusivity


def bond_number(rho_l: float, rho_v: float, gravity: float, length: float, surface_tension: float) -> float:
    return (rho_l - rho_v) * gravity * length**2 / surface_tension


def mean_free_path(T: float, p: float, molecular_diameter: float) -> float:
    """Hard-sphere mean free path [m] of a gas at T [K], p [Pa], molecule diameter d_m [m]."""
    return BOLTZMANN_CONSTANT * T / (math.sqrt(2) * math.pi * molecular_diameter**2 * p)


def knudsen(mean_free_path_m: float, length: float) -> float:
    return mean_free_path_m / length


def classify_knudsen(kn: float) -> str:
    """Flow-regime label for a Knudsen number (conventional thresholds)."""
    if kn < 0.01:
        return "continuum"
    if kn < 0.1:
        return "slip"
    if kn < 10:
        return "transition"
    return "free molecular"
