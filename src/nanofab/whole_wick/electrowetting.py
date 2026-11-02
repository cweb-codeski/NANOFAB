"""Electrowetting (EWOD) as the only electrical-assistance mechanism in the whole-wick model.

Lippmann-Young (mugele_baret2005, not held, ASSUMED):

    cos(theta_V) = cos(theta_0) + eps_r eps_0 V^2 / (2 gamma t_dielectric)

The second term is the electrowetting number eta. The voltage acts on dryout **only** through the
contact angle and therefore the capillary-pressure limit:

    Pc_max(V) = 2 gamma cos(theta_V) / r_eff

With r_eff = 2 (l^2 - pi d^2/4) / (pi d), this is vaartstra2019 Eq. 3 at V = 0 (the
contact-line force over the liquid area of a square pillar cell).

Bounds (physical, not tuning):
    cos(theta_V) <= 1. The largest possible gain is Pc_max(V)/Pc_max(0) <= 1/cos(theta_0),
    which is about 3.5% for theta_0 = 15 deg.
    Optional saturation angle theta_sat: cos(theta_V) <= cos(theta_sat). If theta_0 is already
    below theta_sat, the voltage has no effect.

Not modeled here (separate mechanisms that need their own physics and experimental support):
electrohydrodynamic pumping, charge transport or leakage, direct field enhancement of evaporation,
and dielectric breakdown or charging. Electrical power is unresolved until an electrode and circuit
model exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

EPSILON_0 = 8.8541878128e-12  # F/m, CODATA 2018


@dataclass(frozen=True)
class DielectricCoating:
    name: str
    epsilon_r: float
    t_dielectric: float  # m
    saturation_angle: float = 0.0  # rad; 0 = only the cos <= 1 bound

    def __post_init__(self):
        if self.epsilon_r < 1 or self.t_dielectric <= 0:
            raise ValueError(f"coating {self.name}: need epsilon_r >= 1 and t_dielectric > 0")
        if not 0 <= self.saturation_angle < math.pi / 2:
            raise ValueError(f"coating {self.name}: saturation angle must be in [0, 90) deg")


def electrowetting_number(V, coating: DielectricCoating, gamma: float):
    """eta = eps_r eps_0 V^2 / (2 gamma t) [-]."""
    V = np.asarray(V, dtype=float)
    if np.any(V < 0):
        raise ValueError("applied voltage must be >= 0 (Lippmann term depends on V^2; give the magnitude)")
    return coating.epsilon_r * EPSILON_0 * V**2 / (2.0 * gamma * coating.t_dielectric)


def cos_theta_V(theta_0, V, coating: DielectricCoating | None, gamma: float):
    """(cos theta_V, saturated mask). Unbounded Lippmann value capped at min(1, cos theta_sat)."""
    c0 = np.cos(np.asarray(theta_0, dtype=float))
    V = np.asarray(V, dtype=float)
    if coating is None:
        if np.any(V > 0):
            raise ValueError("voltage applied to a cell type without a dielectric coating")
        return np.broadcast_to(c0, V.shape).copy(), np.zeros(V.shape, bool)
    raw = c0 + electrowetting_number(V, coating, gamma)
    cap = min(1.0, math.cos(coating.saturation_angle))
    capped = np.minimum(raw, cap)
    out = np.maximum(c0, capped)  # theta_0 below theta_sat: no change, never a decrease
    return out, (raw > cap) & (V > 0)


def r_eff_from_geometry(d: float, l: float) -> float:
    """Effective meniscus radius of a square pillar array so that 2 gamma cos/r_eff = Eq. 3."""
    return 2.0 * (l**2 - math.pi * d**2 / 4.0) / (math.pi * d)


def pc_max(gamma: float, cos_theta, r_eff) -> np.ndarray:
    """Pc_max = 2 gamma cos(theta) / r_eff [Pa]."""
    return 2.0 * gamma * np.asarray(cos_theta, dtype=float) / np.asarray(r_eff, dtype=float)


def voltage_for_full_wetting(theta_0: float, coating: DielectricCoating, gamma: float) -> float:
    """Voltage [V] at which Lippmann-Young reaches the cap (cos = min(1, cos theta_sat))."""
    gap = min(1.0, math.cos(coating.saturation_angle)) - math.cos(theta_0)
    if gap <= 0:
        return 0.0
    return math.sqrt(gap * 2.0 * gamma * coating.t_dielectric / (coating.epsilon_r * EPSILON_0))
