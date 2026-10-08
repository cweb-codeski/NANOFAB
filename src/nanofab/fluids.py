"""Saturated working-fluid properties from CoolProp.

All properties are evaluated on the saturation curve at temperature ``T`` [K]:
liquid properties at quality Q = 0, vapor properties at Q = 1. Using CoolProp instead of typed-in
constants keeps properties consistent with temperature during sweeps.

Source: Bell et al., "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the
Open-Source Thermophysical Property Library CoolProp", Ind. Eng. Chem. Res. 53 (2014) 2498.
For water, CoolProp implements IAPWS-95 (EOS) and the IAPWS surface-tension and transport
correlations.
"""

from __future__ import annotations

from dataclasses import dataclass

from CoolProp.CoolProp import PropsSI

UNIVERSAL_GAS_CONSTANT = 8.314462618  # J/(mol K), CODATA 2018
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K, exact (SI 2019)


@dataclass(frozen=True)
class SaturatedProperties:
    """Saturated properties at one temperature. All SI."""

    T: float  # K
    p_sat: float  # Pa
    rho_l: float  # kg/m3
    rho_v: float  # kg/m3
    mu_l: float  # Pa*s
    k_l: float  # W/(m K)
    cp_l: float  # J/(kg K)
    surface_tension: float  # N/m
    h_lv: float  # J/kg
    molar_mass: float  # kg/mol
    dpsat_dT: float  # Pa/K


class Fluid:
    """Thin wrapper around CoolProp for one pure fluid (default: water)."""

    def __init__(self, name: str = "Water"):
        self.name = name
        self.molar_mass = PropsSI("M", self.name)  # kg/mol
        self.T_triple = PropsSI("Ttriple", self.name)
        self.T_critical = PropsSI("Tcrit", self.name)

    def _check_T(self, T: float) -> None:
        if not (self.T_triple <= T < self.T_critical):
            raise ValueError(
                f"T = {T:.2f} K outside saturation range of {self.name} "
                f"({self.T_triple:.2f}-{self.T_critical:.2f} K). Did you pass degC instead of K?"
            )

    def _sat(self, prop: str, T: float, quality: float) -> float:
        self._check_T(T)
        return PropsSI(prop, "T", T, "Q", quality, self.name)

    def p_sat(self, T: float) -> float:
        """Saturation pressure [Pa] at T [K]."""
        return self._sat("P", T, 0)

    def T_sat(self, p: float) -> float:
        """Saturation temperature [K] at p [Pa]."""
        return PropsSI("T", "P", p, "Q", 0, self.name)

    def rho_l(self, T: float) -> float:
        return self._sat("D", T, 0)

    def rho_v(self, T: float) -> float:
        return self._sat("D", T, 1)

    def mu_l(self, T: float) -> float:
        return self._sat("V", T, 0)

    def k_l(self, T: float) -> float:
        return self._sat("L", T, 0)

    def cp_l(self, T: float) -> float:
        return self._sat("C", T, 0)

    def surface_tension(self, T: float) -> float:
        return self._sat("I", T, 0)

    def h_lv(self, T: float) -> float:
        """Latent heat of vaporization [J/kg]."""
        return self._sat("H", T, 1) - self._sat("H", T, 0)

    def dpsat_dT(self, T: float, dT: float = 0.01) -> float:
        """Slope of the saturation curve [Pa/K], central difference."""
        return (self.p_sat(T + dT) - self.p_sat(T - dT)) / (2 * dT)

    def saturated(self, T: float) -> SaturatedProperties:
        """All saturated properties needed by the analytical models at T [K]."""
        return SaturatedProperties(
            T=T,
            p_sat=self.p_sat(T),
            rho_l=self.rho_l(T),
            rho_v=self.rho_v(T),
            mu_l=self.mu_l(T),
            k_l=self.k_l(T),
            cp_l=self.cp_l(T),
            surface_tension=self.surface_tension(T),
            h_lv=self.h_lv(T),
            molar_mass=self.molar_mass,
            dpsat_dT=self.dpsat_dT(T),
        )
