import math

import numpy as np
import pytest

from nanofab import capillary, permeability
from nanofab.geometry import solid_fraction

GAMMA = 0.0712  # N/m, water near 30 C (test constant, not a model input)
MU = 7.97e-4  # Pa s


def test_young_laplace_sphere_and_cylinder():
    R = 10e-6
    assert capillary.young_laplace(GAMMA, R, R) == pytest.approx(2 * GAMMA / R)
    assert capillary.young_laplace(GAMMA, R) == pytest.approx(GAMMA / R)
    assert capillary.young_laplace(GAMMA, math.inf, math.inf) == 0.0


def test_pillar_capillary_force_balance():
    d, l, theta = 10e-6, 30e-6, math.radians(10)
    pc = capillary.pillar_max_capillary_pressure(GAMMA, d, l, theta)
    liquid_area = l**2 - math.pi * d**2 / 4
    # pressure force on the liquid area equals contact-line force
    assert pc * liquid_area == pytest.approx(GAMMA * math.cos(theta) * math.pi * d)
    assert 1e3 < pc < 1e5  # kPa-scale for tens-of-micron pillars


def test_overlapping_pillars_rejected():
    with pytest.raises(ValueError, match="overlap"):
        capillary.pillar_max_capillary_pressure(GAMMA, 30e-6, 30e-6, 0.0)


def test_capillary_permeability_tradeoff_with_pitch():
    """Larger pitch at fixed d: lower capillary pressure, higher permeability (the central tradeoff)."""
    d, h, theta = 10e-6, 40e-6, math.radians(10)
    pitches = np.linspace(20e-6, 80e-6, 12)
    pc = [capillary.pillar_max_capillary_pressure(GAMMA, d, l, theta) for l in pitches]
    K = [permeability.pillar_wick_permeability(d, l, h) for l in pitches]
    assert np.all(np.diff(pc) < 0)
    assert np.all(np.diff(K) > 0)


def test_brinkman_tall_pillar_limit():
    d, l = 10e-6, 30e-6
    K_s = permeability.sangani_acrivos_permeability(d, l)
    assert permeability.pillar_wick_permeability(d, l, 1e-1) == pytest.approx(K_s, rel=1e-3)


def test_brinkman_short_pillar_limit_is_film_flow():
    """h << sqrt(K_s): wick behaves like a shear-free liquid film, K -> h^2/3."""
    K_s = 1e-9
    h = 1e-8
    assert K_s * permeability.brinkman_height_factor(K_s, h) == pytest.approx(h**2 / 3, rel=1e-3)


def test_brinkman_profile_average_matches_K_eff():
    d, l, h = 10e-6, 30e-6, 40e-6
    K_s = permeability.sangani_acrivos_permeability(d, l)
    dpdx = -1e6
    z = np.linspace(0, h, 4001)
    u = np.array([permeability.brinkman_velocity_profile(zi, K_s, h, dpdx, MU) for zi in z])
    u_avg = np.trapezoid(u, z) / h
    K_eff = permeability.pillar_wick_permeability(d, l, h)
    assert u_avg == pytest.approx(-K_eff / MU * dpdx, rel=1e-4)


def test_sangani_acrivos_drummond_tahir_form():
    """Spot value: phi = 0.1 -> K/a^2 = (1/(8*0.1))*(-ln 0.1 - 1.476 + 0.2 - 0.01774 + 0.004076)."""
    phi = 0.1
    a = 5e-6
    l = math.sqrt(math.pi * a**2 / phi)
    assert solid_fraction(2 * a, l) == pytest.approx(phi)
    expected = a**2 / (8 * phi) * (-math.log(phi) - 1.476 + 2 * phi - 1.774 * phi**2 + 4.076 * phi**3)
    assert permeability.sangani_acrivos_permeability(2 * a, l) == pytest.approx(expected)


def test_permeability_from_flow_recovers_darcy():
    K, A, L, dP = 2e-11, 1e-8, 15e-6, 100.0
    Q = K / MU * dP / L * A
    assert permeability.permeability_from_flow(Q, A, L, dP, MU) == pytest.approx(K)


def test_meniscus_radius_and_cap_depth():
    pc = 5000.0
    R = capillary.meniscus_radius(GAMMA, pc)
    assert capillary.young_laplace(GAMMA, R, R) == pytest.approx(pc)
    assert capillary.spherical_cap_depth(R, 0.0) == 0.0
    assert capillary.spherical_cap_depth(math.inf, 10e-6) == 0.0
