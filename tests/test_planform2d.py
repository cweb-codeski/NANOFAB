"""Verification of the 2D planform field-treatment solver (src/nanofab/electroosmosis/planform2d.py)."""

import numpy as np
import pytest

from nanofab.electroosmosis import planform2d as pf
from nanofab.electroosmosis.formation_factor import formation_factor_maxwell_2d
from nanofab.geometry import solid_fraction
from nanofab.permeability import sangani_acrivos_permeability

D, L = 10e-6, 30e-6  # NANOFAB baseline planform


@pytest.fixture(scope="module")
def base():
    return pf.field_treatment(D, L, 20)


def test_boundary_residuals_are_small(base):
    assert base.field.residual_equipotential < 1e-8 and base.field.residual_mirror < 1e-8
    for flow in (base.flow_solved, base.flow_uniform):
        assert max(flow.residuals.values()) < 1e-7


def test_spectral_convergence_of_ratio():
    a = pf.field_treatment(D, L, 12).uniform_over_solved
    b = pf.field_treatment(D, L, 24).uniform_over_solved
    assert abs(a - b) < 1e-5


def test_current_conservation_and_formation_factor(base):
    currents = [base.field.current(x0) for x0 in (1.0, 0.6, 0.2, 0.0)]
    assert np.ptp(currents) < 1e-10
    assert 1 / base.inv_F == pytest.approx(formation_factor_maxwell_2d(solid_fraction(D, L)), rel=1e-5)


def test_mean_gradient_and_periodicity(base):
    s = np.linspace(-1, 1, 101)
    right = base.field.phi(np.ones_like(s), s)
    left = base.field.phi(-np.ones_like(s), s)
    assert np.allclose(right, -1.0, atol=1e-9) and np.allclose(left, 1.0, atol=1e-9)  # drop = E0 l


def test_similitude_solved_field_gives_plug_flux(base):
    assert base.q_solved == pytest.approx(base.inv_F, rel=1e-9)


def test_mass_conservation_on_sections(base):
    for flow in (base.flow_solved, base.flow_uniform):
        fluxes = [flow.section_flux(x0) for x0 in (1.0, 0.5, 0.1, 0.0)]
        assert np.ptp(fluxes) < 1e-9 * abs(flow.q)
        assert fluxes[0] == pytest.approx(flow.q, rel=1e-9)


def test_reversal_and_zero_field(base):
    rev = pf.field_treatment(D, L, 20, E_sign=-1.0)
    assert rev.q_solved == pytest.approx(-base.q_solved, rel=1e-12)
    assert rev.q_uniform == pytest.approx(-base.q_uniform, rel=1e-12)
    assert pf.solve_flow(D / L, 16, slip_modes=np.zeros(16)).q == 0.0


def test_stokes_solver_matches_sangani_acrivos_permeability():
    assert pf.permeability_2d(D, L, 20) == pytest.approx(sangani_acrivos_permeability(D, L), rel=2e-3)


def test_dilute_limit_ratio_is_one_half():
    # isolated insulating cylinder: tangential field 2 E0 sin(theta) vs E0 sin(theta) -> ratio 1/2
    r = pf.field_treatment(0.01 * L, L, 12)
    assert r.inv_F == pytest.approx(1.0, abs=1e-3)
    assert r.uniform_over_solved == pytest.approx(0.5, abs=1e-3)


def test_uniform_field_gives_less_flux_than_solved_field():
    for d_over_l in (0.15, 0.33, 0.5, 0.62):
        r = pf.field_treatment(d_over_l * L, L, 20)
        assert 0.5 < r.uniform_over_solved < 0.65
