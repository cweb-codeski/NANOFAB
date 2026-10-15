"""Interfacial (kinetic) evaporation at a liquid-vapor interface.

Schrage relation (net evaporative mass flux):
    m'' = (2 s / (2 - s)) * sqrt(M / (2 pi R)) * ( p_sat(T_lv) / sqrt(T_lv) - p_v / sqrt(T_v) )

    m''  : evaporative mass flux [kg/(m^2 s)]
    s    : accommodation (evaporation/condensation) coefficient [-], 0 < s <= 1
    M    : molar mass [kg/mol];  R : universal gas constant [J/(mol K)]
    T_lv : interface temperature [K];  T_v, p_v : vapor temperature [K] and pressure [Pa]

Heat flux removed: q'' = m'' h_lv.

Source: Schrage, A Theoretical Study of Interphase Mass Transfer, Columbia Press (vaartstra2019
ref. [22] gives 1954; 1953 is also commonly cited; unresolved); Carey, Liquid-Vapor Phase-Change
Phenomena (edition/chapter unverified). Somasundaram et al., Int. J. Therm. Sci. 134 (2018) used
it for micropillar wicks (per vaartstra2019 Sec. 1, ref. [21]). Vaartstra et al.
(2019) instead use the linearized Boltzmann transport equation; see ``sone_bte_mass_flux``.

Linearized interfacial heat-transfer coefficient (for saturated vapor p_v = p_sat(T_v), small
superheat):
    h_i = h_lv * (2 s / (2 - s)) * sqrt(M / (2 pi R)) * (1/sqrt(T)) * (dp_sat/dT - p_sat / (2 T))
so that q'' ~= h_i (T_lv - T_v). Comparing 1/h_i to the conduction resistances shows how much of
the total resistance is kinetic.

Assumptions: no disjoining-pressure or curvature correction to p_sat (valid away from the
nanometer-scale adsorbed film), vapor near equilibrium (no Knudsen-layer resolution).
"""

from __future__ import annotations

import math

from .fluids import UNIVERSAL_GAS_CONSTANT, Fluid


def schrage_prefactor(accommodation: float, molar_mass: float) -> float:
    """(2 s/(2 - s)) sqrt(M/(2 pi R)) [s K^0.5 / m], equivalently [kg K^0.5 / (m^2 s Pa)]."""
    if not (0 < accommodation <= 1):
        raise ValueError(f"accommodation coefficient must be in (0, 1], got {accommodation}")
    return (2 * accommodation / (2 - accommodation)) * math.sqrt(
        molar_mass / (2 * math.pi * UNIVERSAL_GAS_CONSTANT)
    )


def schrage_mass_flux(
    T_interface: float,
    T_vapor: float,
    accommodation: float,
    fluid: Fluid,
    p_vapor: float | None = None,
) -> float:
    """Net evaporative mass flux m'' [kg/(m^2 s)]. ``p_vapor`` defaults to p_sat(T_vapor)."""
    if p_vapor is None:
        p_vapor = fluid.p_sat(T_vapor)
    pref = schrage_prefactor(accommodation, fluid.molar_mass)
    return pref * (fluid.p_sat(T_interface) / math.sqrt(T_interface) - p_vapor / math.sqrt(T_vapor))


def schrage_heat_flux(
    T_interface: float,
    T_vapor: float,
    accommodation: float,
    fluid: Fluid,
    p_vapor: float | None = None,
) -> float:
    """Evaporative heat flux q'' = m'' h_lv [W/m^2], h_lv evaluated at T_interface."""
    m = schrage_mass_flux(T_interface, T_vapor, accommodation, fluid, p_vapor)
    return m * fluid.h_lv(T_interface)


def interfacial_htc(T: float, accommodation: float, fluid: Fluid) -> float:
    """Linearized kinetic heat-transfer coefficient h_i [W/(m^2 K)] at temperature T [K]."""
    pref = schrage_prefactor(accommodation, fluid.molar_mass)
    p = fluid.p_sat(T)
    dpdT = fluid.dpsat_dT(T)
    return fluid.h_lv(T) * pref * (dpdT - p / (2 * T)) / math.sqrt(T)


# Magnitudes only: the thesis prints c4* = -2.13204, d4* = -0.44675; the sign is carried by the
# explicit minus in the relation below (the IJHMT PDF text drops minus signs).
SONE_C4 = 2.13204  # linearized BKW half-space evaporation coefficient (Sone 2000)
SONE_D4 = 0.44675  # temperature-jump coefficient (Sone 2000)


def sone_bte_mass_flux(
    T_interface: float,
    T_vapor: float,
    accommodation: float,
    fluid: Fluid,
    p_vapor: float | None = None,
) -> float:
    """Evaporative mass flux [kg/(m^2 s)] from a REDUCED linearized-BTE interface relation.

    Form taken from Vaartstra, Lu & Wang, IJHMT 136 (2019) 170, Eqs. (5)-(7):
        p_k / p_sat(T_l) - 1 = -C4 u_k / sqrt(2 R_s T_l)          (Sone 2000, BKW molecules)
        q''_e = rho_k u_k h_lv,   p_k = far-field vapor pressure
    Reductions relative to the paper (this is NOT a paper-verified implementation):
      * the temperature-jump relation (V2019 Eq. 6, coefficient D4) is not used;
      * rho_k is approximated by p_sat(T_l) / (R_s T_l) (linear theory), giving
            m'' = (p_sat(T_l) - p_k) * sqrt(2 / (R_s T_l)) / C4(a);
      * non-unity accommodation enters as C4 -> C4 + 2 sqrt(pi) (1 - a) / a. V2019 attributes its
        nonunity-accommodation treatment to Meland & Ytrehus, Phys. Fluids 15 (2003) 1348-1350,
        DOI 10.1063/1.1564097, which is NOT held: this correction was written from memory and is an
        open source gap (config/sources.yaml: meland_ytrehus2003).

    The flux is within ~3-5% of the Schrage relation at a = 0.052 for water at 100 C. That is a
    cross-model check between two approximate relations, not verification of either.
    """
    if not (0 < accommodation <= 1):
        raise ValueError(f"accommodation coefficient must be in (0, 1], got {accommodation}")
    if p_vapor is None:
        p_vapor = fluid.p_sat(T_vapor)
    R_s = UNIVERSAL_GAS_CONSTANT / fluid.molar_mass
    c4 = SONE_C4 + 2 * math.sqrt(math.pi) * (1 - accommodation) / accommodation
    return (fluid.p_sat(T_interface) - p_vapor) * math.sqrt(2 / (R_s * T_interface)) / c4


def sone_bte_heat_flux(
    T_interface: float,
    T_vapor: float,
    accommodation: float,
    fluid: Fluid,
    p_vapor: float | None = None,
) -> float:
    """q''_e = m'' h_lv [W/m^2] for the linearized-BTE interface model (Vaartstra et al. 2019)."""
    return sone_bte_mass_flux(T_interface, T_vapor, accommodation, fluid, p_vapor) * fluid.h_lv(T_interface)


def comsol_bte_expression(accommodation: float, molar_mass: float, T_vapor: float, h_lv: float) -> str:
    """COMSOL expression for the linearized-BTE evaporative heat flux q''_e (positive = leaving).

    Requires an Interpolation function ``psat`` (K -> Pa), as for the Schrage expression.
    """
    R_s = UNIVERSAL_GAS_CONSTANT / molar_mass
    c4 = SONE_C4 + 2 * math.sqrt(math.pi) * (1 - accommodation) / accommodation
    return (
        f"{h_lv:.6g}[J/kg]*(psat(T) - psat({T_vapor:.2f}[K]))"
        f"*sqrt(2/({R_s:.6g}[J/(kg*K)]*T))/{c4:.6g}"
    )


def clausius_clapeyron_slope(T: float, h_lv: float, molar_mass: float, p_sat: float) -> float:
    """Ideal-gas Clausius-Clapeyron estimate dp_sat/dT = h_lv p_sat M / (R T^2) [Pa/K].

    Independent check on the CoolProp slope (should agree within a few percent below ~100 C).
    """
    return h_lv * p_sat * molar_mass / (UNIVERSAL_GAS_CONSTANT * T**2)


def comsol_schrage_expression(accommodation: float, molar_mass: float, T_vapor: float, h_lv: float) -> str:
    """Text to paste into a COMSOL Boundary Heat Source on the meniscus (negative = heat removed).

    COMSOL must provide p_sat(T) itself, as an Interpolation function named ``psat`` (argument in K,
    result in Pa) built from the table written by ``analysis/baseline_estimates.py --psat-table``.
    """
    pref = schrage_prefactor(accommodation, molar_mass)
    return (
        f"-{h_lv:.6g}[J/kg]*{pref:.6g}[s*K^0.5/m]"
        f"*(psat(T)/sqrt(T) - psat({T_vapor:.2f}[K])/sqrt({T_vapor:.2f}[K]))"
    )
