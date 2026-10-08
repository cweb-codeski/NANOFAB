"""Conduction resistances for the solid path from heater to evaporating interface.

Area-specific resistance R'' [m^2 K / W] is used for 1D layers so different footprints compare
directly; total resistance R [K/W] = R'' / A.

1D slab:
    R'' = t / k

Pillar wick, conduction normal to the substrate (pillars and liquid in parallel):
    k_eff = phi k_s + (1 - phi) k_l
This is an upper bound on how well the wick conducts; the real path from pillar to interface is
2D/3D and is what the COMSOL unit cell resolves.

Spreading (constriction) resistance, circular source of radius a centered on a circular plate of
radius b and thickness t, cooled on the opposite face by a uniform h. Average-temperature form of
Lee, Song, Au & Moran, "Constriction/Spreading Resistance Model for Electronics Packaging",
ASME/JSME Thermal Engineering Conf. (1995):
    eps = a/b, tau = t/b, Bi = h b / k, lambda_c = pi + 1/(sqrt(pi) eps)
    Phi_c = (tanh(lambda_c tau) + lambda_c/Bi) / (1 + (lambda_c/Bi) tanh(lambda_c tau))
    psi_avg = 0.5 (1 - eps)^1.5 Phi_c
    R_sp = psi_avg / (sqrt(pi) k a)          [K/W]
Total source-to-sink resistance = R_sp + t/(k pi b^2) + 1/(h pi b^2).
"""

from __future__ import annotations

import math


def slab_resistance_area(thickness: float, conductivity: float) -> float:
    """Area-specific 1D conduction resistance R'' = t/k [m^2 K/W]."""
    if thickness < 0 or conductivity <= 0:
        raise ValueError("thickness must be >= 0 and conductivity > 0")
    return thickness / conductivity


def series_area_resistance(*resistances: float) -> float:
    """Sum of area-specific resistances in series [m^2 K/W]."""
    return float(sum(resistances))


def parallel_conductivity(phi: float, k_solid: float, k_liquid: float) -> float:
    """Effective conductivity [W/(m K)] of pillars and liquid conducting in parallel."""
    if not 0 <= phi <= 1:
        raise ValueError("phi must be in [0, 1]")
    return phi * k_solid + (1 - phi) * k_liquid


def series_conductivity(phi: float, k_solid: float, k_liquid: float) -> float:
    """Effective conductivity [W/(m K)] of pillars and liquid in series (lower bound)."""
    return 1.0 / (phi / k_solid + (1 - phi) / k_liquid)


def spreading_resistance(
    source_radius: float, plate_radius: float, thickness: float, conductivity: float, sink_htc: float
) -> float:
    """Average spreading resistance R_sp [K/W] (Lee et al. 1995). Zero when source = plate size."""
    a, b, t, k, h = source_radius, plate_radius, thickness, conductivity, sink_htc
    if not (0 < a <= b):
        raise ValueError("Need 0 < source_radius <= plate_radius")
    eps = a / b
    tau = t / b
    bi = h * b / k
    lam = math.pi + 1 / (math.sqrt(math.pi) * eps)
    th = math.tanh(lam * tau)
    phi_c = (th + lam / bi) / (1 + (lam / bi) * th)
    psi_avg = 0.5 * (1 - eps) ** 1.5 * phi_c
    return psi_avg / (math.sqrt(math.pi) * k * a)


def equivalent_radius(area: float) -> float:
    """Radius [m] of the circle with the same area (to map square chips onto Lee et al.)."""
    return math.sqrt(area / math.pi)
