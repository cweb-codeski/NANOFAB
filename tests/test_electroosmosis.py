"""Verification of the passive-vs-electroosmotic wick model (src/nanofab/electroosmosis/)."""

import dataclasses
import math

import numpy as np
import pytest
from scipy.integrate import quad

from nanofab import dryout as dryout1d
from nanofab.electroosmosis import formation_factor as ff
from nanofab.electroosmosis import physics as ek
from nanofab.electroosmosis import wick1d as w1
from nanofab.electroosmosis.inputs import (
    MissingInputError,
    baseline_liquid,
    baseline_wick,
    load_eo_inputs,
    scenario_mobilities,
)
from nanofab.params import ParameterError, _to_si
from nanofab.permeability import pillar_wick_permeability
from nanofab.provenance import COMSOL_EXTRACTED, MODEL_CALCULATED, PRIMARY_SOURCE
from nanofab.whole_wick.darcy import FeedSpec, solve_darcy
from nanofab.whole_wick.grid import Grid

T100 = 373.15


@pytest.fixture(scope="module")
def inp():
    return load_eo_inputs()


@pytest.fixture(scope="module")
def liquid(inp):
    return baseline_liquid(inp)


@pytest.fixture(scope="module")
def wick(inp, liquid):
    return baseline_wick(inp, liquid)


@pytest.fixture(scope="module")
def mob(inp, liquid):
    return scenario_mobilities(inp, liquid)[1]["mobility"]  # zeta = -200 mV


# ------------------------------------------------------------------ units and constants

def test_new_units_convert_to_si():
    assert _to_si(-200, "mV", "z") == (pytest.approx(-0.2), "V")
    assert _to_si(0.79, "uS/cm", "s") == (pytest.approx(7.9e-5), "S/m")
    assert _to_si(1e-6, "mol/L", "c") == (pytest.approx(1e-3), "mol/m3")
    assert _to_si(2260, "kJ/kg", "h") == (pytest.approx(2.26e6), "J/kg")


def test_mobility_sign_and_dimensions():
    mu = 2.8e-4
    m = ek.electroosmotic_mobility(55.4, -0.2, mu)
    assert m > 0  # negative zeta: flow along E (toward the cathode)
    assert ek.electroosmotic_mobility(55.4, 0.2, mu) == pytest.approx(-m)
    # eps0 [C/(V m)] * zeta [V] / mu [Pa s] -> m^2/(V s); value check against hand arithmetic
    assert m == pytest.approx(55.4 * 8.8541878128e-12 * 0.2 / mu, rel=1e-12)
    with pytest.raises(ValueError):
        ek.electroosmotic_mobility(0.5, -0.1, mu)


def test_debye_length_reproduces_pujahari_values(inp):
    eps = inp.si("literature_scenarios.pujahari2022.relative_permittivity")
    for c, lam in zip(inp.si("literature_scenarios.pujahari2022.concentration"),
                      inp.si("literature_scenarios.pujahari2022.debye_length")):
        # paper quotes one significant figure; computed values must round to it
        assert ek.debye_length(eps, T100, c) == pytest.approx(lam, rel=0.1)


def test_gouy_chapman_displacement_limits_and_quadrature():
    assert ek.gouy_chapman_displacement(1e-6, T100) == pytest.approx(1.0, rel=1e-6)
    vt = ek.thermal_voltage(T100)
    for zeta in (-0.05, -0.1, -0.2):
        zt = abs(zeta) / vt
        t0 = math.tanh(zt / 4)
        psi = lambda s: 4 * math.atanh(t0 * math.exp(-s)) / zt  # psi/zeta vs n/lambda  # noqa: E731
        val, _ = quad(psi, 0, 60)
        assert ek.gouy_chapman_displacement(zeta, T100) == pytest.approx(val, rel=1e-8)
    d = [ek.gouy_chapman_displacement(z, T100) for z in (-0.01, -0.05, -0.1, -0.2)]
    assert all(a > b for a, b in zip(d, d[1:])) and d[-1] < 1


def test_dukhin_and_field_effect_formulas():
    vt = ek.thermal_voltage(T100)
    assert ek.dukhin_number_lower_bound(1e-7, 0.0, T100, 1e-5) == 0.0
    assert ek.dukhin_number_lower_bound(1e-7, -0.2, T100, 1e-5) == pytest.approx(2e-2 * (math.cosh(0.1 / vt) - 1))
    sig = ek.grahame_surface_charge(55.4, T100, 1e-3, -0.2)
    assert sig < 0
    # small-zeta Grahame -> eps zeta / lambda_D
    lam = ek.debye_length(55.4, T100, 1.0)
    assert ek.grahame_surface_charge(55.4, T100, 1.0, -1e-4) == pytest.approx(55.4 * ek.EPSILON_0 * -1e-4 / lam, rel=1e-6)
    r = ek.field_effect_charge_ratio(3.9, 1e-6, 10.0, sig)
    assert r == pytest.approx(3.9 * ek.EPSILON_0 * 10.0 / 1e-6 / abs(sig))


def test_onsager_open_circuit_reduction_matches_electroviscous_coefficient():
    K, mu, sigma, F, m, G = 3e-11, 2.8e-4, 7.9e-5, 1.19, 3.5e-7, 2e5
    # U = (K/mu) G + (m/F) E,   J = (sigma/F) E + (m/F) G = 0 at open circuit
    E = -m * G / sigma
    U = K / mu * G + m / F * E
    beta = ek.electroviscous_coefficient(m, mu, K, sigma, F)
    assert U == pytest.approx(K / mu * G * (1 - beta), rel=1e-12)


def test_joule_heating_scales_with_E_squared():
    q = ek.joule_heat_flux(1e-4, np.array([1e4, 2e4]), 25e-6, 1.2)
    assert q[1] == pytest.approx(4 * q[0])
    assert q[0] == pytest.approx(1e-4 / 1.2 * 1e8 * 25e-6)


# ------------------------------------------------------------------ formation factor

def test_maxwell_formation_factor_limits():
    assert ff.formation_factor_maxwell_2d(0.0) == 1.0
    assert ff.formation_factor_maxwell_2d(0.5) == pytest.approx(3.0)
    with pytest.raises(ValueError):
        ff.formation_factor_maxwell_2d(1.0)


def test_numerical_formation_factor_matches_maxwell_at_baseline():
    s = ff.formation_factor_numerical(10e-6, 30e-6, (128, 256, 512))
    assert abs(s.rel_diff_vs_maxwell) < 2e-4
    assert 0.7 < s.observed_order < 1.5  # sub-cell fraction staircase: first order
    assert all(a > b for a, b in zip(s.F, s.F[1:]))  # monotone convergence from above


def test_formation_factor_from_current_roundtrip():
    F = 1.19
    sigma, A, V, L = 7.9e-5, 2.5e-10, 1.0, 15e-6
    I = sigma / F * A * V / L
    assert ff.formation_factor_from_current(I, sigma, A, V, L) == pytest.approx(F)


# ------------------------------------------------------------------ cases and forced flow

def test_paired_cases_share_geometry_and_differ_only_in_forcing(wick, liquid, mob):
    a, b = w1.paired_cases(wick, liquid, w1.Forcing(mob, 1e4))
    assert a.wick is b.wick and a.liquid is b.liquid
    assert a.forcing.is_passive and not b.forcing.is_passive
    w1.assert_identical_except_forcing(a, b)
    altered = dataclasses.replace(b, wick=dataclasses.replace(wick, K=wick.K * 1.01))
    with pytest.raises(ValueError):
        w1.assert_identical_except_forcing(a, altered)


def test_passive_forced_flow_is_linear_darcy_with_zero_intercept(wick, liquid):
    a, _ = w1.paired_cases(wick, liquid, w1.Forcing())
    dP = np.linspace(0, 2000, 11)
    r = w1.forced_flow(a, dP)
    assert r["Q_total"][0] == 0.0
    assert np.allclose(r["Q_total"], wick.K * wick.area * dP / (liquid.mu * wick.L), rtol=1e-14)
    assert np.allclose(np.diff(r["Q_total"], 2), 0, atol=1e-14 * r["Q_total"].max())
    assert np.all(r["Q_eo"] == 0)


def test_eo_flow_superposes_and_has_correct_limits(wick, liquid, mob):
    dP = np.linspace(-500, 2000, 6)
    for E in (0.0, 5e3, 1e5):
        _, b = w1.paired_cases(wick, liquid, w1.Forcing(mob, E))
        r = w1.forced_flow(b, dP)
        assert np.allclose(r["Q_total"], r["Q_pressure"] + r["Q_eo"], rtol=1e-14)
        assert np.allclose(r["Q_eo"], mob * E / wick.F * wick.area)  # independent of dP
        if E == 0:
            assert np.all(r["Q_eo"] == 0)
    # Q_pressure is the same function of dP in A and B (K unchanged by the field)
    a, b = w1.paired_cases(wick, liquid, w1.Forcing(mob, 1e5))
    assert np.array_equal(w1.forced_flow(a, dP)["Q_pressure"], w1.forced_flow(b, dP)["Q_pressure"])


def test_eo_direction_and_monotonicity(wick, liquid, mob):
    E = np.linspace(-1e5, 1e5, 21)
    U = [w1.eo_superficial_velocity(w1.Case("B", wick, liquid, w1.Forcing(mob, e))) for e in E]
    assert np.all(np.diff(U) > 0)  # zeta < 0: U_EO increases with E along +x
    assert U[10] == 0.0
    pos = w1.Case("B", wick, liquid, w1.Forcing(-mob, 1e4))  # zeta > 0 reverses the flow
    assert w1.eo_superficial_velocity(pos) == pytest.approx(-U[11])


def test_eo_velocity_definitions_are_distinct(wick, liquid, mob):
    v = w1.eo_velocities(w1.Case("B", wick, liquid, w1.Forcing(mob, 1e4)))
    assert v["hs_slip"] == pytest.approx(mob * 1e4)
    assert v["superficial"] == pytest.approx(v["hs_slip"] / wick.F)
    assert v["intrinsic_pore_average"] == pytest.approx(v["superficial"] / (1 - wick.phi))
    assert v["superficial"] < v["intrinsic_pore_average"] < v["hs_slip"]


def test_apparent_permeability_is_not_a_property(wick, liquid, mob):
    a, b = w1.paired_cases(wick, liquid, w1.Forcing(mob, 1e4))
    G = np.array([1e3, 1e4, 1e5, 1e6])
    assert np.allclose(w1.apparent_permeability(a, G), wick.K)
    Kb = w1.apparent_permeability(b, G)
    assert np.all(Kb > wick.K) and np.all(np.diff(Kb) < 0)  # depends on the driving gradient
    assert Kb == pytest.approx(wick.K + liquid.mu * w1.eo_superficial_velocity(b) / G)


def test_equivalent_pressure_definition(wick, liquid, mob):
    b = w1.Case("B", wick, liquid, w1.Forcing(mob, 1e4))
    dP_eo = w1.equivalent_eo_pressure(b)
    # passive flow at dP_eo equals the EO flow at zero dP
    a = w1.Case("A", wick, liquid, w1.Forcing())
    assert w1.forced_flow(a, dP_eo)["Q_total"] == pytest.approx(w1.forced_flow(b, 0.0)["Q_total"], rel=1e-12)
    assert w1.eo_pressure_ratio(b) == pytest.approx(dP_eo / wick.Pc_max)


# ------------------------------------------------------------------ evaporating half-wick

def test_passive_regressions(wick, liquid):
    a = w1.Case("A", wick, liquid, w1.Forcing())
    q = 50e4
    assert w1.required_capillary_pressure(a, q)[0] == pytest.approx(
        dryout1d.wick_pressure_drop(q, wick.K, wick.h, wick.L, liquid.rho, liquid.mu, liquid.h_fg), rel=1e-12)
    x = np.linspace(0, wick.L, 9)
    assert np.allclose(w1.capillary_pressure_profile(a, q, x),
                       -dryout1d.pressure_profile(x, q, wick.K, wick.h, wick.L, liquid.rho, liquid.mu, liquid.h_fg))
    # analytical K reproduces the repo's legacy flat-meniscus estimate (config/model_provenance.yaml: 100.6)
    a_flat = dataclasses.replace(a, wick=dataclasses.replace(wick, K=pillar_wick_permeability(wick.d, wick.l, wick.h)))
    assert w1.passive_dryout_heat_flux(a_flat) / 1e4 == pytest.approx(100.6, abs=0.05)


def test_fv_conserves_mass_and_blocks_flow_at_closed_end(wick, liquid, mob):
    for E in (0.0, 1e4, 1e5, -1e4):
        c = w1.Case("B", wick, liquid, w1.Forcing(mob, E))
        s = w1.solve_fv(c, 80e4, n=200)
        assert abs(s.mass_residual) < 1e-10
        assert s.face_flow_per_width[-1] == 0.0
        cell_balance = s.face_flow_per_width[:-1] - s.face_flow_per_width[1:]
        assert np.allclose(cell_balance, s.evaporation_per_width / 200, rtol=1e-9)


def test_closed_form_profile_matches_fv_with_second_order_convergence(wick, liquid, mob):
    c = w1.Case("B", wick, liquid, w1.Forcing(mob, 5e4))
    errs = []
    for n in (50, 100, 200):
        s = w1.solve_fv(c, 100e4, n)
        errs.append(np.max(np.abs(s.Pc - w1.capillary_pressure_profile(c, 100e4, s.x))))
    assert errs[-1] < 1e-6 * wick.Pc_max or (errs[0] / errs[1] > 3.5 and errs[1] / errs[2] > 3.5)


def test_2d_whole_wick_solver_matches_1d_passive_profile(wick, liquid):
    """Independent solver: whole_wick.darcy (2D FV, both edges fed) vs the 1D closed form."""
    g = Grid(400, 2, 2 * wick.L, 1e-3)
    ones = np.ones(g.shape)
    q = 60e4
    res = solve_darcy(g, wick.K * ones, wick.h * ones, liquid.rho * ones, liquid.mu * ones, q / liquid.h_fg * ones,
                      FeedSpec(edges={"left": 0.0, "right": 0.0}))
    a = w1.Case("A", wick, liquid, w1.Forcing())
    xh = np.minimum(g.x, g.lx - g.x)
    exact = w1.capillary_pressure_profile(a, q, xh)
    assert np.max(np.abs(res.Pc[0] - exact)) < 1e-5 * exact.max()


def test_dryout_gain_closed_form_properties():
    assert w1.dryout_gain(0.0) == 1.0
    eps = 1e-6
    assert (w1.dryout_gain(eps) - w1.dryout_gain(-eps)) / (2 * eps) == pytest.approx(1.0, rel=1e-5)
    P = np.linspace(-0.9, 20, 400)
    assert np.all(np.diff(w1.dryout_gain(P)) > 0)
    # EO-dominated limit: U_0 -> U_EO, so G -> Pi/2 (not Pi); small Pi: G ~ 1 + Pi
    assert w1.dryout_gain(1e6) / 5e5 == pytest.approx(1.0, rel=3e-3)
    assert w1.dryout_gain(1e-3) == pytest.approx(1 + 1e-3, rel=1e-6)


@pytest.mark.parametrize("E", [0.0, 5e3, 3e4, 1e5, 3e5, -2e4])
def test_dryout_closed_form_matches_independent_fv_bisection(wick, liquid, mob, E):
    c = w1.Case("B", wick, liquid, w1.Forcing(mob, E))
    q_cf = w1.assisted_dryout_heat_flux(c)
    q_fv = w1.dryout_heat_flux_fv(c, n=400)
    # interior maximum: cell-centre max misses the peak by O(dx^2); opposing field: wall max, O(dx)
    assert q_fv == pytest.approx(q_cf, rel=2e-3 if E < 0 else 1e-5)


def test_eo_dryout_reduces_to_passive_and_increases_with_aiding_field(wick, liquid, mob):
    base = w1.passive_dryout_heat_flux(w1.Case("A", wick, liquid, w1.Forcing()))
    qs = [w1.assisted_dryout_heat_flux(w1.Case("B", wick, liquid, w1.Forcing(mob, E))) for E in (0, 1e3, 1e4, 1e5)]
    assert qs[0] == base
    assert all(b > a for a, b in zip(qs, qs[1:]))
    assert w1.assisted_dryout_heat_flux(w1.Case("B", wick, liquid, w1.Forcing(mob, -1e4))) < base


def test_flooding_threshold(wick, liquid, mob):
    c = w1.Case("B", wick, liquid, w1.Forcing(mob, 1e5))
    qf = w1.flooding_heat_flux(c)
    assert w1.closed_end_capillary_pressure(c, qf) == pytest.approx(0.0, abs=1e-9 * wick.Pc_max)
    assert w1.closed_end_capillary_pressure(c, 0.5 * qf) < 0 < w1.closed_end_capillary_pressure(c, 1.5 * qf)
    assert w1.flooding_heat_flux(w1.Case("A", wick, liquid, w1.Forcing())) == 0.0
    # at the predicted dryout point the closed end is not flooded while Pi <= 4
    assert w1.eo_pressure_ratio(c) <= 4
    assert w1.closed_end_capillary_pressure(c, w1.assisted_dryout_heat_flux(c)) >= 0


# ------------------------------------------------------------------ inputs and provenance

def test_inputs_provenance(inp, wick):
    assert inp.leaf("wick.permeability").provenance == COMSOL_EXTRACTED
    assert inp.leaf("wick.pillar_pitch").provenance == PRIMARY_SOURCE
    for key in ("zeta_potential", "relative_permittivity", "electric_field", "concentration"):
        p = inp.leaf(f"literature_scenarios.pujahari2022.{key}")
        assert p.provenance == PRIMARY_SOURCE and p.source == "pujahari2022" and p.source_location
    assert not [p for p in inp.config.registry.values() if p.provenance == MODEL_CALCULATED]
    assert wick.K == pytest.approx(2.9698e-11)


def test_requires_source_blocks_device_predictions(inp, tmp_path):
    assert "zeta_potential_device" in inp.missing
    with pytest.raises(MissingInputError, match="zeta_potential_device"):
        inp.require("zeta_potential_device", "applied_voltage")
    inp.require("not_listed")  # unknown names are not blocked
    bad = tmp_path / "bad.yaml"
    bad.write_text("requires_source:\n  zeta_x: {unit: mV, needed_for: x, value: -50}\n", encoding="utf-8")
    with pytest.raises(ParameterError, match="has a value"):
        load_eo_inputs(bad)


def test_reynolds_numbers_scale_linearly(wick, liquid):
    c = w1.Case("A", wick, liquid, w1.Forcing())
    r = w1.reynolds_numbers(c, np.array([0.01, 0.02]))
    assert r["Re_K"][1] == pytest.approx(2 * r["Re_K"][0])
    assert r["Re_gap"][0] == pytest.approx(liquid.rho * 0.01 / (1 - wick.phi) * (wick.l - wick.d) / liquid.mu)
