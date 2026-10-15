"""Checks against Vaartstra, Lu & Wang, IJHMT 136 (2019) 170."""

import math

import pytest

from nanofab import capillary, dryout, evaporation, permeability
from nanofab.fluids import Fluid
from nanofab.params import load_parameters

T_SAT = 373.15
THETA = math.radians(15)


@pytest.fixture(scope="module")
def water():
    return Fluid("Water")


def q_dry(water, d, l, h, L=5e-3):
    s = water.saturated(T_SAT)
    Pc = capillary.pillar_max_capillary_pressure(s.surface_tension, d, l, THETA)
    K = permeability.pillar_wick_permeability(d, l, h)
    return dryout.dryout_heat_flux(K, h, L, Pc, s.rho_l, s.mu_l, s.h_lv)


def test_eq3_matches_force_balance(water):
    gamma = water.surface_tension(T_SAT)
    d, l = 10e-6, 30e-6
    eq3 = 4 * gamma * math.cos(THETA) / (d * (4 / math.pi * (l / d) ** 2 - 1))
    assert capillary.pillar_max_capillary_pressure(gamma, d, l, THETA) == pytest.approx(eq3)


def test_legacy_flat_meniscus_sanity_check(water):
    """LEGACY SANITY CHECK of the approximate analytical estimate (flat-meniscus K, uniform q).

    The 0-40% window is a plausibility range for a legacy estimate, not a source-backed bound: it
    is not established that the flat-meniscus K, constant-K transport and uniform evaporation
    together overpredict the V2019 numerical result (82 W/cm^2, Sec. 3.2) for every geometry.
    """
    q = q_dry(water, 10e-6, 30e-6, 25e-6)
    assert 82e4 < q < 1.4 * 82e4


def test_dryout_scales_with_h_squared_at_fixed_ratios(water):
    """Paper Fig. 5: q''_dry ~ h^2 at d/h = 0.4, l/h = 1.2."""
    h1, h2 = 10e-6, 40e-6
    ratio = q_dry(water, 0.4 * h2, 1.2 * h2, h2) / q_dry(water, 0.4 * h1, 1.2 * h1, h1)
    assert ratio == pytest.approx((h2 / h1) ** 2, rel=1e-6)


def test_bte_close_to_schrage_at_fitted_accommodation(water):
    a = 0.052
    bte = evaporation.sone_bte_heat_flux(T_SAT + 1, T_SAT, a, water)
    sch = evaporation.schrage_heat_flux(T_SAT + 1, T_SAT, a, water)
    assert bte == pytest.approx(sch, rel=0.05)


def test_bte_zero_at_equilibrium(water):
    assert evaporation.sone_bte_mass_flux(T_SAT, T_SAT, 0.052, water) == pytest.approx(0.0, abs=1e-12)


def test_kinetic_resistance_below_paper_total(water):
    """Consistency check only: V2019 reports h_e > 60 kW/m^2K, so R_total < 1/60000 m^2K/W, and any
    kinetic part must satisfy R_int < R_total < 1/60000. Passing is necessary, not sufficient; it
    does not bound the kinetic share (R_int / R_total = R_int h_e > 60000 R_int is a lower bound)."""
    h_i = evaporation.sone_bte_heat_flux(T_SAT + 0.1, T_SAT, 0.052, water) / 0.1
    assert 1 / h_i < 1 / 60e3


def test_parameters_hold_paper_baseline():
    ps = load_parameters()
    assert ps["pillar_height"] == pytest.approx(25e-6)
    assert ps["accommodation_coefficient"] == pytest.approx(0.052)
    assert ps["reference_dryout_heat_flux"] == pytest.approx(82e4)
