"""Baseline analytical estimates assembled from the individual physics modules.

These are the numbers a COMSOL unit-cell result should land near (permeability within tens of
percent, interfacial resistance within the same order of magnitude). Large disagreements mean
either the COMSOL setup or an analytical assumption is wrong; both are worth knowing.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import capillary, conduction, dimensionless, dryout, evaporation, geometry, permeability
from .fluids import Fluid
from .params import ParameterSet


@dataclass
class BaselineEstimates:
    # geometry
    d: float
    l: float
    h: float
    phi: float
    # fluid state
    T_v: float
    p_v: float
    # capillarity and flow
    Pc_max: float  # Pa
    K_s: float  # m^2
    K_eff: float  # m^2
    q_dry: float  # W/m^2 (footprint)
    q_target: float  # W/m^2
    dP_at_target: float  # Pa
    U_feed: float  # m/s, superficial, at target heat flux
    u_pore: float  # m/s
    Re_pore: float
    Ca: float
    Pe_pore: float  # lateral: u_pore * gap / alpha_l
    Pe_vertical: float  # V2019 A.1: u_e h / alpha_l, u_e = q''/(rho h_lv) toward the interface
    Ja: float  # V2019 A.1: cp dT / h_lv with dT = 10 K; small -> lateral convection negligible
    Re_h: float  # V2019 A.1: rho U h / mu with superficial U at the feed
    Bo_pitch: float
    gravity_over_Pc_vertical: float  # rho g L / Pc_max if the wick were vertical
    mean_free_path: float
    Kn_pitch: float
    # thermal resistances, area-specific per footprint [m^2 K/W]
    R_substrate: float
    R_wick_parallel: float
    R_wick_series: float
    R_interface: float
    h_i: float  # W/m^2K, per interface area
    accommodation: float

    @property
    def R_total_optimistic(self) -> float:
        return self.R_substrate + self.R_wick_parallel + self.R_interface

    @property
    def dT_at_target(self) -> float:
        return self.q_target * self.R_total_optimistic

    @property
    def h_evap_upper_bound(self) -> float:
        """Footprint-based h_evap if the wick conducted as well as the parallel bound [W/m^2K]."""
        return 1.0 / (self.R_wick_parallel + self.R_interface)


def interface_area_fraction(d: float, l: float) -> float:
    """Flat-meniscus interface area per footprint area (l^2 - pi d^2/4)/l^2 [-]."""
    return 1.0 - geometry.solid_fraction(d, l)


def compute_baseline(ps: ParameterSet, fluid: Fluid | None = None) -> BaselineEstimates:
    fluid = fluid or Fluid(ps["working_fluid"])
    d, l, h = ps["pillar_diameter"], ps["pillar_pitch"], ps["pillar_height"]
    L = ps["wick_feed_length"]
    T_v = ps["vapor_temperature"]
    s_acc = ps["accommodation_coefficient"]
    q = ps["target_heat_flux"]
    g = ps["gravity"]

    props = fluid.saturated(T_v)
    phi = geometry.solid_fraction(d, l)

    Pc = capillary.pillar_max_capillary_pressure(props.surface_tension, d, l, ps["receding_contact_angle"])
    K_s = permeability.sangani_acrivos_permeability(d, l)
    K_eff = permeability.pillar_wick_permeability(d, l, h)

    q_dry = dryout.dryout_heat_flux(K_eff, h, L, Pc, props.rho_l, props.mu_l, props.h_lv)
    dP = dryout.wick_pressure_drop(q, K_eff, h, L, props.rho_l, props.mu_l, props.h_lv)
    U = dryout.superficial_velocity_at_feed(q, L, h, props.rho_l, props.h_lv)
    u_pore = U / (1 - phi)

    Re = dimensionless.reynolds(props.rho_l, u_pore, geometry.gap(d, l), props.mu_l)
    Ca = dimensionless.capillary_number(props.mu_l, u_pore, props.surface_tension)
    alpha_l = props.k_l / (props.rho_l * props.cp_l)
    Pe = dimensionless.peclet(u_pore, geometry.gap(d, l), alpha_l)
    Pe_v = dimensionless.peclet(q / (props.rho_l * props.h_lv), h, alpha_l)
    Ja = props.cp_l * 10.0 / props.h_lv
    Re_h = dimensionless.reynolds(props.rho_l, U, h, props.mu_l)
    Bo = dimensionless.bond_number(props.rho_l, props.rho_v, g, l, props.surface_tension)
    grav_ratio = dryout.gravity_head(props.rho_l, g, L) / Pc
    mfp = dimensionless.mean_free_path(T_v, props.p_sat, ps["molecular_kinetic_diameter"])

    k_si = ps["silicon_thermal_conductivity"]
    R_sub = conduction.slab_resistance_area(ps["substrate_thickness"], k_si)
    R_wick_par = conduction.slab_resistance_area(h, conduction.parallel_conductivity(phi, k_si, props.k_l))
    R_wick_ser = conduction.slab_resistance_area(h, conduction.series_conductivity(phi, k_si, props.k_l))
    h_i = evaporation.interfacial_htc(T_v, s_acc, fluid)
    R_int = 1.0 / (h_i * interface_area_fraction(d, l))

    return BaselineEstimates(
        d=d, l=l, h=h, phi=phi, T_v=T_v, p_v=props.p_sat,
        Pc_max=Pc, K_s=K_s, K_eff=K_eff, q_dry=q_dry, q_target=q, dP_at_target=dP,
        U_feed=U, u_pore=u_pore, Re_pore=Re, Ca=Ca, Pe_pore=Pe, Pe_vertical=Pe_v, Ja=Ja, Re_h=Re_h,
        Bo_pitch=Bo,
        gravity_over_Pc_vertical=grav_ratio, mean_free_path=mfp,
        Kn_pitch=dimensionless.knudsen(mfp, l),
        R_substrate=R_sub, R_wick_parallel=R_wick_par, R_wick_series=R_wick_ser,
        R_interface=R_int, h_i=h_i, accommodation=s_acc,
    )


def accommodation_sensitivity(T_v: float, fluid: Fluid, d: float, l: float, values=(0.03, 0.1, 0.3, 0.5, 1.0)):
    """[(s, h_i, R''_interface per footprint)] for a sweep of accommodation coefficients."""
    frac = interface_area_fraction(d, l)
    out = []
    for s in values:
        h_i = evaporation.interfacial_htc(T_v, s, fluid)
        out.append((s, h_i, 1.0 / (h_i * frac)))
    return out


def pitch_tradeoff(ps: ParameterSet, fluid: Fluid, pitches) -> list[tuple[float, float, float, float]]:
    """[(l, Pc_max, K_eff, q_dry)] at fixed d, h, L to show the capillary/permeability tradeoff."""
    d, h, L = ps["pillar_diameter"], ps["pillar_height"], ps["wick_feed_length"]
    props = fluid.saturated(ps["vapor_temperature"])
    rows = []
    for l in pitches:
        if l <= d:
            continue
        Pc = capillary.pillar_max_capillary_pressure(props.surface_tension, d, l, ps["receding_contact_angle"])
        K = permeability.pillar_wick_permeability(d, l, h)
        qd = dryout.dryout_heat_flux(K, h, L, Pc, props.rho_l, props.mu_l, props.h_lv)
        rows.append((l, Pc, K, qd))
    return rows


def unit_cell_flow_check(K_eff: float, l: float, h: float, mu: float, dP: float) -> float:
    """Expected volumetric flow [m^3/s] through a quarter cell (width l/2, length l/2) at dP.

    Matches COMSOL Stage 1: Q = (K/mu) (dP / (l/2)) * (l/2 * h).
    """
    length = l / 2
    area = (l / 2) * h
    return K_eff / mu * dP / length * area
