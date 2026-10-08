import pytest

from nanofab.evaporation import clausius_clapeyron_slope
from nanofab.fluids import Fluid


@pytest.fixture(scope="module")
def water():
    return Fluid("Water")


def test_water_reference_values_at_25C(water):
    # Reference values: IAPWS / Incropera Table A.6 at 298.15 K
    s = water.saturated(298.15)
    assert s.p_sat == pytest.approx(3170, rel=0.01)
    assert s.surface_tension == pytest.approx(0.0720, rel=0.01)
    assert s.mu_l == pytest.approx(8.9e-4, rel=0.02)
    assert s.rho_l == pytest.approx(997, rel=0.005)
    assert s.h_lv == pytest.approx(2.442e6, rel=0.005)
    assert s.k_l == pytest.approx(0.607, rel=0.02)


def test_boiling_point_at_one_atm(water):
    assert water.T_sat(101325) == pytest.approx(373.12, abs=0.1)


def test_dpsat_dT_matches_clausius_clapeyron(water):
    T = 303.15
    s = water.saturated(T)
    cc = clausius_clapeyron_slope(T, s.h_lv, s.molar_mass, s.p_sat)
    assert s.dpsat_dT == pytest.approx(cc, rel=0.03)


def test_celsius_input_rejected(water):
    with pytest.raises(ValueError, match="degC"):
        water.p_sat(30.0)
