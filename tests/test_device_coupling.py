"""Coupled pressure / HTC / conduction iteration on MANUFACTURED_FIXTURE tables (VERIFICATION ONLY)."""

import numpy as np
import pytest

from nanofab.device import verification as V
from nanofab.device.diagnostics import energy_balance, mass_balance, regime_diagnostics
from nanofab.device.model import SolverSettings, SolverStatus
from nanofab.device.solver import search_dryout, solve_at_load
from nanofab.device.thermal import SlabMesh, reconstruct_energy_balance, solve_slab

S = SolverSettings(Nx=60, Nz=12)


@pytest.fixture(scope="module")
def problem():
    return V.fixture_problem()


def test_forward_solve_conserves_mass_and_energy(problem):
    sol = solve_at_load(problem, 5e5, S)
    assert sol.status == SolverStatus.ADMISSIBLE and sol.converged
    mb, eb = mass_balance(problem, sol), energy_balance(problem, sol)
    assert mb["global_normalized"] <= V.GATES["mass_global"]
    assert mb["local_max_normalized"] <= V.GATES["mass_local"]
    assert mb["darcy_cell_max_normalized"] <= V.GATES["mass_local"]
    assert eb["global_normalized"] <= V.GATES["energy_global"]
    assert eb["local_max_normalized"] <= V.GATES["energy_local"]
    assert sol.Q_faces[-1] == 0.0


def test_evaporation_comes_from_conduction_not_input(problem):
    sol = solve_at_load(problem, 7e5, S)
    # h_e rises with suction, so evaporation shifts toward the dry end; q_e is not q_in
    assert np.ptp(sol.q_e) / np.mean(sol.q_e) > 1e-3
    assert np.allclose(np.ptp(sol.q_in), 0.0)
    assert sol.q_e[-1] > sol.q_e[0]
    assert np.sum(sol.q_e * sol.mesh.dx) == pytest.approx(np.sum(sol.q_in * sol.mesh.dx), rel=1e-12)


def test_superheat_increases_with_load(problem):
    a, b = solve_at_load(problem, 3e5, S), solve_at_load(problem, 6e5, S)
    assert np.max(b.slab.theta_e) > np.max(a.slab.theta_e) > 0


class _Recorder:
    def __init__(self, inner):
        self.inner, self.max_query = inner, -np.inf
        self.limits, self.snap, self.name = inner.limits, inner.snap, inner.name
        self.Pc_lo, self.Pc_hi = inner.Pc_lo, inner.Pc_hi

    def _rec(self, P):
        self.max_query = max(self.max_query, float(np.max(np.atleast_1d(P))))

    def tau(self, P):
        self._rec(P)
        return self.inner.tau(P)

    def psi(self, P):
        self._rec(P)
        return self.inner.psi(P)

    def psi_inverse(self, S_):
        out = self.inner.psi_inverse(S_)
        self._rec(out)
        return out

    def h_e(self, P):
        self._rec(P)
        return self.inner.h_e(P)

    def liquid_section_area_per_width(self, P):
        return self.inner.liquid_section_area_per_width(P)


def test_closure_never_queried_above_cap():
    p = V.fixture_problem()
    rec = _Recorder(p.closure)
    p.closure = rec
    res = search_dryout(p, S)
    assert res.status == SolverStatus.PHYSICAL_CAPILLARY_LIMIT
    assert rec.max_query <= res.interval.Pc_cap + rec.snap


def test_regime_diagnostics_are_reported(problem):
    sol = solve_at_load(problem, 5e5, S)
    d = regime_diagnostics(problem, sol)
    for k in ("Re_h_superficial", "Re_h_pore_estimate", "Ca_superficial", "Bo_pitch", "We_pitch_superficial",
              "Pe_vertical_max", "Pe_axial_L_feed", "Ja_max_superheat", "vapor_Mach_max",
              "max_cell_pressure_change_over_capillary_scale", "cell_to_cell_dT_over_superheat_max"):
        assert np.isfinite(d[k]), k
    assert d["Knudsen"] is None
    assert d["within_supported_ranges"]


def test_reconstructed_energy_balance_independent_of_matrix():
    m = SlabMesh(1e-3, 2e-4, 30, 6)
    q = np.linspace(1e5, 3e5, 30)
    he = np.linspace(4e4, 8e4, 30)
    sol = solve_slab(m, 148.0, q, he)
    rec = reconstruct_energy_balance(sol)
    assert abs(rec["global_abs_W_per_m"]) / rec["input_W_per_m"] < 1e-10
    assert rec["max_local_abs_W_per_m"] / rec["local_scale_W_per_m"] < 1e-9
    assert rec["robin_max_abs_W_m2"] < 1e-6


def test_coupled_refinement_and_tolerance_gates():
    res = V.suite_coupled(levels=((25, 4), (50, 8), (100, 16)), tol_level=(50, 8))
    assert res["passed"], [g for g in res["gates"] if not g["passed"]]
