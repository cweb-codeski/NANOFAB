"""Electrokinetic relations for thin-double-layer electroosmosis in a pillar wick.

Sign convention: x points from the liquid feed toward the closed (far) end of the wick. A positive
velocity or flow rate is feed -> closed end, the direction capillary pumping drives the liquid.
E is the macroscopic applied field component along x (E > 0 points from the feed-side electrode,
the anode, to the far-end electrode, the cathode).

Helmholtz-Smoluchowski (HS) slip velocity just outside a thin electric double layer (EDL):

    u_HS = mu_EO E,     mu_EO = -eps_r eps_0 zeta / mu        [m^2/(V s)]

With a negative wall zeta potential (silica/water) mu_EO > 0, so the liquid moves along E,
toward the cathode. Source: textbook electrokinetics (e.g. Kirby, Micro- and Nanoscale Fluid
Mechanics, 2010, Ch. 6; not held, cited from memory). pujahari2022 Sec. 1 cites the HS slip
formulation (dutta_beskok2001) but solves the full Poisson-Boltzmann body force instead.

Superficial (Darcy) EO velocity of a porous medium with uniform zeta on every wetted surface,
insulating solid, thin EDL and negligible surface conduction (similitude between EO velocity and
electric field; cummings2000, not held):

    U_EO = mu_EO E / F,   F = sigma_bulk / sigma_eff (electrical formation factor, geometric)

The double-layer diagnostics below quantify when those assumptions fail.
"""

from __future__ import annotations

import math

import numpy as np

from ..fluids import BOLTZMANN_CONSTANT

EPSILON_0 = 8.8541878128e-12  # F/m, CODATA 2018
ELEMENTARY_CHARGE = 1.602176634e-19  # C, exact (SI 2019)
AVOGADRO = 6.02214076e23  # 1/mol, exact (SI 2019)


def electroosmotic_mobility(relative_permittivity: float, zeta: float, mu: float) -> float:
    """HS electroosmotic mobility mu_EO = -eps_r eps_0 zeta / mu [m^2/(V s)].

    Positive for zeta < 0: the liquid then moves in the direction of E.
    """
    if relative_permittivity < 1 or mu <= 0:
        raise ValueError("need relative_permittivity >= 1 and mu > 0")
    return -relative_permittivity * EPSILON_0 * zeta / mu


def hs_slip_velocity(mobility, E):
    """HS slip velocity outside the EDL [m/s] (unobstructed-channel EO velocity)."""
    return np.asarray(mobility, dtype=float) * np.asarray(E, dtype=float)


def superficial_eo_velocity(mobility, E, formation_factor: float):
    """Darcy (superficial, per total cross-section) EO velocity U_EO = mu_EO E / F [m/s]."""
    if formation_factor < 1:
        raise ValueError("formation factor F = sigma_bulk/sigma_eff must be >= 1 for an insulating solid")
    return hs_slip_velocity(mobility, E) / formation_factor


def thermal_voltage(T: float, valence: int = 1) -> float:
    """k_B T / (z e) [V]."""
    return BOLTZMANN_CONSTANT * T / (abs(valence) * ELEMENTARY_CHARGE)


def number_density(concentration_mol_m3: float) -> float:
    """Ion number density n0 [1/m^3] of each species of a symmetric salt at c [mol/m^3]."""
    return concentration_mol_m3 * AVOGADRO


def debye_length(relative_permittivity: float, T: float, concentration_mol_m3: float, valence: int = 1) -> float:
    """Debye length of a symmetric z:z electrolyte, sqrt(eps k T / (2 n0 z^2 e^2)) [m]."""
    n0 = number_density(concentration_mol_m3)
    eps = relative_permittivity * EPSILON_0
    return math.sqrt(eps * BOLTZMANN_CONSTANT * T / (2 * n0 * valence**2 * ELEMENTARY_CHARGE**2))


def _legendre_chi2(t: float, terms: int = 400) -> float:
    """Legendre chi function chi_2(t) = sum_k t^(2k+1) / (2k+1)^2, |t| < 1."""
    k = np.arange(terms)
    return float(np.sum(t ** (2 * k + 1) / (2 * k + 1) ** 2))


def gouy_chapman_displacement(zeta: float, T: float, valence: int = 1) -> float:
    """EO flux-displacement length of one thin Gouy-Chapman EDL, in units of the Debye length.

    For EO with no pressure gradient u(n) = u_HS (1 - psi(n)/zeta), so each charged wall removes a
    flux u_HS * delta, delta = int_0^inf psi/zeta dn. With the Gouy-Chapman potential
    tanh(z e psi / 4kT) = tanh(z e zeta / 4kT) exp(-n/lambda_D):

        delta / lambda_D = 4 chi_2(tanh(zt/4)) / zt,   zt = z e zeta / kT

    -> 1 for |zeta| << kT/e (Debye-Hueckel), and < 1 for large |zeta|. Our derivation from the
    Gouy-Chapman solution (PROJECT_DERIVED; not checked against a held source).
    """
    zt = abs(zeta) / thermal_voltage(T, valence)
    if zt < 1e-8:
        return 1.0
    return 4.0 * _legendre_chi2(math.tanh(zt / 4.0)) / zt


def thin_edl_flux_factor(debye: float, zeta: float, T: float, half_width: float, valence: int = 1) -> float:
    """First-order finite-EDL correction to the EO flux of a pore of half-width ``half_width``.

    1 - delta/half_width (one charged wall per half-width). Meaningful only when the result is
    close to 1; it is an error estimate for the HS model, not a finite-EDL model.
    """
    if half_width <= 0:
        raise ValueError("half_width must be positive")
    return 1.0 - gouy_chapman_displacement(zeta, T, valence) * debye / half_width


def dukhin_number_lower_bound(debye: float, zeta: float, T: float, length: float, valence: int = 1) -> float:
    """Lower bound of the Dukhin number Du = K_sigma / (sigma_bulk a) (diffuse-layer conduction only).

    Bikerman surface conductance of a symmetric electrolyte with equal ion diffusivities:
        Du = (2 lambda_D / a) (1 + 3m/z^2) (cosh(z e zeta / 2kT) - 1)
    with m = 2 eps (kT/ze)^2 / (3 mu D). Dropping the electroosmotic part (m = 0) gives a lower
    bound that needs no ion diffusivity. Formula cited from memory (Bikerman / Lyklema; not held).
    Du << 1 is required for the EO-velocity / electric-field similitude used for U_EO = mu_EO E / F.
    """
    zt = abs(zeta) / thermal_voltage(T, valence)
    return 2.0 * debye / length * (math.cosh(zt / 2.0) - 1.0)


def grahame_surface_charge(relative_permittivity: float, T: float, concentration_mol_m3: float, zeta: float,
                           valence: int = 1) -> float:
    """Diffuse-layer charge density [C/m^2] implied by zeta (Grahame, symmetric z:z, Gouy-Chapman).

    sigma_d = sqrt(8 eps n0 k T) sinh(z e zeta / 2kT); same sign as zeta.
    """
    eps = relative_permittivity * EPSILON_0
    n0 = number_density(concentration_mol_m3)
    zt = zeta / thermal_voltage(T, valence)
    return math.sqrt(8 * eps * n0 * BOLTZMANN_CONSTANT * T) * math.sinh(zt / 2.0)


def field_effect_charge_ratio(oxide_relative_permittivity: float, oxide_thickness: float, voltage_across_oxide: float,
                              native_surface_charge: float) -> float:
    """|induced charge| / |native diffuse charge| for a conducting substrate under a thin oxide.

    If the silicon under the oxide is conducting, the liquid-to-substrate potential difference
    induces sigma_ind = eps_ox eps_0 dV / t_ox (oxide capacitance in series with the much larger
    EDL capacitance). A ratio that is not << 1 means zeta varies along the wick with the local
    liquid potential ("field-effect" electroosmosis) and the uniform-zeta model is invalid.
    """
    if oxide_thickness <= 0 or native_surface_charge == 0:
        raise ValueError("need oxide_thickness > 0 and a nonzero native surface charge")
    induced = oxide_relative_permittivity * EPSILON_0 * abs(voltage_across_oxide) / oxide_thickness
    return induced / abs(native_surface_charge)


def joule_heat_flux(conductivity: float, E, liquid_height: float, formation_factor: float):
    """Ohmic heating per footprint area [W/m^2]: sigma_eff E^2 h with sigma_eff = sigma_bulk / F."""
    E = np.asarray(E, dtype=float)
    return conductivity / formation_factor * E**2 * liquid_height


def streaming_to_conduction_ratio(mobility: float, pressure_gradient_magnitude, conductivity: float, E):
    """Streaming-current / conduction-current density ratio (thin EDL, Onsager-reciprocal coefficients).

    i_stream = (mu_EO / F) |dp/dx|,  i_cond = (sigma / F) E. A non-small ratio means the applied
    field is not set by the electrodes alone (electro-hydraulic coupling not modeled).
    """
    E = np.asarray(E, dtype=float)
    with np.errstate(divide="ignore"):
        return mobility * np.asarray(pressure_gradient_magnitude, dtype=float) / (conductivity * E)


def electroviscous_coefficient(mobility: float, mu: float, permeability: float, conductivity: float,
                               formation_factor: float) -> float:
    """beta = mu mu_EO^2 / (F K sigma): open-circuit flow reduction U = (K/mu) G (1 - beta).

    With no external circuit the streaming potential opposes the pressure-driven flow; beta is the
    fractional reduction (thin EDL, no surface conduction). It bounds how much the passive case
    depends on whether the electrodes are open- or short-circuited.
    """
    return mu * mobility**2 / (formation_factor * permeability * conductivity)
