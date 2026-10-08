import math

import numpy as np
import pytest

from nanofab import conduction, dimensionless, dryout, evaporation
from nanofab.fluids import Fluid

RHO, MU, HLV = 996.0, 7.97e-4, 2.43e6
K, H, L = 2e-11, 40e-6, 5e-3


@pytest.fixture(scope="module")
def water():
    return Fluid("Water")


def test_dryout_inverts_pressure_drop():
    pc = 5000.0
    q_dry = dryout.dryout_heat_flux(K, H, L, pc, RHO, MU, HLV)
    assert dryout.wick_pressure_drop(q_dry, K, H, L, RHO, MU, HLV) == pytest.approx(pc)


def test_dryout_scales_inverse_square_with_length():
    q1 = dryout.dryout_heat_flux(K, H, L, 5000, RHO, MU, HLV)
    q2 = dryout.dryout_heat_flux(K, H, 2 * L, 5000, RHO, MU, HLV)
    assert q1 / q2 == pytest.approx(4.0)


def test_pressure_profile_mass_balance():
    """Darcy flow leaving the feed must carry all the liquid evaporated over the wick (energy balance)."""
    q = 50e4
    x = np.linspace(0, L, 20001)
    p = dryout.pressure_profile(x, q, K, H, L, RHO, MU, HLV)
    dpdx0 = (p[1] - p[0]) / (x[1] - x[0])
    U0 = -K / MU * dpdx0
    evaporated_mass_per_width = q * L / HLV
    assert RHO * U0 * H == pytest.approx(evaporated_mass_per_width, rel=1e-3)
    assert U0 == pytest.approx(dryout.superficial_velocity_at_feed(q, L, H, RHO, HLV), rel=1e-3)
    # closed end: zero gradient
    assert abs(p[-1] - p[-2]) / abs(p[1] - p[0]) < 1e-3


def test_schrage_zero_at_equilibrium(water):
    assert evaporation.schrage_mass_flux(303.15, 303.15, 0.5, water) == pytest.approx(0.0, abs=1e-12)


def test_schrage_sign(water):
    assert evaporation.schrage_mass_flux(305.0, 303.15, 0.5, water) > 0  # evaporation
    assert evaporation.schrage_mass_flux(301.0, 303.15, 0.5, water) < 0  # condensation


def test_linearized_htc_matches_small_superheat(water):
    T_v, dT, s = 303.15, 0.01, 0.5
    q = evaporation.schrage_heat_flux(T_v + dT, T_v, s, water)
    h_i = evaporation.interfacial_htc(T_v, s, water)
    assert q / dT == pytest.approx(h_i, rel=0.01)


def test_htc_order_of_magnitude(water):
    """Water near room temperature with s = 1: h_i ~ 1e6 W/m^2K (Carey, Ch. 16)."""
    h_i = evaporation.interfacial_htc(300.0, 1.0, water)
    assert 5e5 < h_i < 2e6


def test_accommodation_bounds():
    with pytest.raises(ValueError):
        evaporation.schrage_prefactor(0.0, 0.018)
    with pytest.raises(ValueError):
        evaporation.schrage_prefactor(1.5, 0.018)


def test_slab_and_parallel_conduction():
    assert conduction.slab_resistance_area(500e-6, 148) == pytest.approx(3.378e-6, rel=1e-3)
    assert conduction.parallel_conductivity(0.0, 148, 0.6) == 0.6
    assert conduction.parallel_conductivity(1.0, 148, 0.6) == 148
    assert conduction.series_conductivity(0.3, 148, 0.6) < conduction.parallel_conductivity(0.3, 148, 0.6)


def test_spreading_resistance_limits():
    assert conduction.spreading_resistance(0.01, 0.01, 1e-3, 148, 1e4) == 0.0
    small = conduction.spreading_resistance(0.002, 0.01, 1e-3, 148, 1e4)
    smaller = conduction.spreading_resistance(0.001, 0.01, 1e-3, 148, 1e4)
    assert 0 < small < smaller


def test_dimensionless_groups():
    assert dimensionless.reynolds(1000, 1e-3, 1e-5, 1e-3) == pytest.approx(0.01)
    lam = dimensionless.mean_free_path(300.0, 3500.0, 0.265e-9)
    assert 1e-6 < lam < 1e-5  # micron-scale mean free path in low-pressure water vapor
    assert dimensionless.classify_knudsen(dimensionless.knudsen(lam, 1e-3)) == "continuum"
    assert dimensionless.classify_knudsen(dimensionless.knudsen(lam, 30e-6)) == "transition"
    assert dimensionless.bond_number(1000, 0, 9.81, 30e-6, 0.072) < 1e-3
