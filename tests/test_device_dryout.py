"""Physical-limit search, admissibility and distinct failure statuses (VERIFICATION ONLY inputs)."""

import numpy as np
import pytest

from nanofab.device import verification as V
from nanofab.device.model import EXIT_CODES, EXIT_CONFIG_ERROR, EXIT_LOAD_ABOVE_LIMIT, SolverSettings, SolverStatus
from nanofab.device.solver import admissible_interval, search_dryout, solve_at_endpoint, solve_prescribed_load

S = SolverSettings(Nx=50, Nz=8)


def test_status_suite(tmp_path):
    res = V.suite_statuses(tmp_path, Nx=50, Nz=8)
    assert res["passed"], [g for g in res["gates"] if not g["passed"]]


def test_limit_bracket_and_monotone_continuation():
    p = V.fixture_problem()
    res = search_dryout(p, S)
    assert res.status == SolverStatus.PHYSICAL_CAPILLARY_LIMIT
    assert res.lower_operating_point < res.upper_bracket
    assert res.relative_interval <= V.GATES["limit_interval_rel"]
    assert res.monotonic_end_suction
    ends = [c["Pc_end_Pa"] for c in res.continuation]
    assert np.all(np.diff(ends) > 0) and ends[-1] < res.interval.Pc_cap
    assert res.endpoint.Pc_end == pytest.approx(res.interval.Pc_cap, rel=1e-12)


def test_endpoint_at_lower_suction_gives_lower_load():
    p = V.fixture_problem()
    a = solve_at_endpoint(p, 1000.0, S)
    b = solve_at_endpoint(p, 2000.0, S)
    assert a.converged and b.converged and a.load < b.load


def test_prescribed_overload_is_classified_not_numerical_failure():
    p = V.fixture_problem()
    sol = solve_prescribed_load(p, 2e6, S)
    assert not sol.converged
    assert sol.status == SolverStatus.PHYSICAL_CAPILLARY_LIMIT
    assert "at or above" in sol.reason


def test_geometry_limit_below_physical_limit_wins(tmp_path):
    p = V.fixture_problem(tmp_path, "geom", geometry_limit=1800.0)
    iv = admissible_interval(p, S)
    assert iv.mechanism == SolverStatus.MENISCUS_GEOMETRY_LIMIT and iv.Pc_cap == pytest.approx(1800.0)


def test_undemonstrated_geometry_limit_is_not_used(tmp_path):
    p = V.fixture_problem()
    p.closure.limits.Pc_geometry_limit = 1500.0
    p.closure.limits.geometry_limit_demonstrated = False
    iv = admissible_interval(p, S)
    assert iv.mechanism == SolverStatus.PHYSICAL_CAPILLARY_LIMIT
    assert any("not demonstrated" in n for n in iv.notes)


def test_coverage_limit_reports_lower_point_only(tmp_path):
    p = V.fixture_problem(tmp_path, "cov", pc_max=1500.0, endpoint_demonstrated=False)
    res = search_dryout(p, S)
    assert res.status == SolverStatus.TABLE_COVERAGE_LIMIT
    assert res.upper_bracket is None and res.lower_operating_point > 0


def test_exit_codes_distinct():
    nonzero = [c for c in EXIT_CODES.values() if c != 0]
    assert len(set(nonzero)) == len(nonzero)
    assert EXIT_CODES[SolverStatus.MISSING_VALIDATED_INPUTS] == 3
    assert EXIT_CONFIG_ERROR not in nonzero and EXIT_LOAD_ABOVE_LIMIT not in nonzero and EXIT_CONFIG_ERROR != 2
    assert {EXIT_CODES[s] for s in (SolverStatus.ADMISSIBLE, SolverStatus.PHYSICAL_CAPILLARY_LIMIT,
                                   SolverStatus.MENISCUS_GEOMETRY_LIMIT)} == {0}
