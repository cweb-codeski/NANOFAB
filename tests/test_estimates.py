import pytest

from nanofab.estimates import compute_baseline, unit_cell_flow_check
from nanofab.fluids import Fluid
from nanofab.params import load_parameters
from nanofab.permeability import permeability_from_flow


@pytest.fixture(scope="module")
def baseline():
    ps = load_parameters()
    return ps, compute_baseline(ps, Fluid(ps["working_fluid"]))


def test_baseline_runs_and_is_physical(baseline):
    _, b = baseline
    assert 0 < b.phi < 0.5
    assert b.K_eff < b.K_s
    assert b.Pc_max > 0 and b.q_dry > 0
    assert b.R_wick_parallel < b.R_wick_series


def test_dryout_consistency(baseline):
    """At q'' = q''_dry the wick pressure drop equals the capillary pressure."""
    _, b = baseline
    assert b.dP_at_target / b.Pc_max == pytest.approx(b.q_target / b.q_dry, rel=1e-9)


def test_quarter_cell_flow_round_trip(baseline):
    ps, b = baseline
    mu = Fluid(ps["working_fluid"]).mu_l(b.T_v)
    Q = unit_cell_flow_check(b.K_eff, b.l, b.h, mu, 1.0)
    K = permeability_from_flow(Q, (b.l / 2) * b.h, b.l / 2, 1.0, mu)
    assert K == pytest.approx(b.K_eff)
