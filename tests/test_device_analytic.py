"""Closed-form limits, manufactured hydraulics/thermal solutions, symmetry and zero load.

All inputs are ANALYTICAL_EXAMPLE or MANUFACTURED_FIXTURE (VERIFICATION ONLY).
"""

import numpy as np
import pytest

from nanofab.device import analytic
from nanofab.device import verification as V


def _assert_suite(res):
    failed = [g for g in res["gates"] if not g["passed"]]
    assert res["passed"], failed


def test_constant_closure_matches_closed_forms():
    _assert_suite(V.suite_analytic(Nx=50, Nz=8))


def test_dryout_constant_equals_transmissivity_quadrature():
    K, h, L, rho, mu, hfg, P = 3e-11, 25e-6, 5e-3, 958.0, 2.8e-4, 2.26e6, 2176.0
    a = analytic.dryout_constant(K, h, L, rho, mu, hfg, P)
    b = analytic.dryout_from_transmissivity(lambda p: K * h, mu, rho, hfg, L, P)
    assert b == pytest.approx(a, rel=1e-12)


def test_moment_identity_uniform():
    xf = np.linspace(0, 5e-3, 101)
    q = np.full(100, 7e5)
    m = analytic.evaporation_moment(xf, q, 958.0, 2.26e6)
    assert m == pytest.approx(7e5 * 5e-3**2 / 2 / (958.0 * 2.26e6), rel=1e-13)


def test_superheats_constant():
    dTe, dTs = analytic.superheats_constant(8e5, 6e4, 100e-6, 148.0)
    assert dTe == pytest.approx(8e5 / 6e4)
    assert dTs - dTe == pytest.approx(8e5 * 100e-6 / 148.0)


def test_manufactured_table_hydraulics_second_order():
    _assert_suite(V.suite_hydraulics(levels=(40, 80, 160)))


def test_manufactured_xz_conduction_second_order():
    _assert_suite(V.suite_thermal())


def test_zero_load_exact():
    _assert_suite(V.suite_zero_load(Nx=40, Nz=8))


def test_mirrored_accounting_and_width_invariance():
    _assert_suite(V.suite_accounting(Nx=50, Nz=8))
