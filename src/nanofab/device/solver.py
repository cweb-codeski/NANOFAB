"""Coupled pressure / cell-HTC / substrate-conduction solver and physical-limit search.

Coupling (V2019 Sec. 2.3 and Appendix A.3, Fig. A1): starting from uniform evaporation, iterate

    q_e''(x) -> hydraulics -> Pc(x) -> h_e(Pc(x)) -> conduction -> q_e''(x)

until the scaled changes of both Pc and q_e'' are below ``rtol``. q_e'' is always the flux returned
by the conduction solution; it is never replaced by q_in''. The hydraulic and thermal problems share
one x grid, so no inter-mesh transfer is needed and integrated evaporation is conserved exactly.

Scaled field change (both fields): ||new - old||_2 / (||new||_2 + sqrt(N) * floor), with
floor = abs_floor_fraction x characteristic scale (Pc: the admissible interval width; q: the mean
input flux, or 1 W/m^2 at zero load). V2019 used a relative tolerance of 1e-6 on Euclidean norms.

Admissibility. Closure tables are never evaluated outside the admissible interval
[Pc_lo, Pc_cap]. If a proposed evaporation field would require Pc > Pc_cap, the step is
backtracked toward the last admissible field (the evaporation moment, and hence the end suction,
is linear in q_e''). Repeated backtracking is reported as NUMERICAL_FAILURE, never as a bracket.

Dryout. Following V2019 Sec. 3.3 ("minimum heat flux at which Delta p_dryout was reached on the
wick"), the limit load is solved directly with the end condition Pc(L) = Pc_cap: for fixed h_e the
conduction problem is linear in load, so the load follows from the evaporation-moment identity
(analytic.py). Every closure evaluation then lies in [Pc_feed, Pc_cap]. The result is checked
against subcritical forward continuation from low load.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .model import (AdmissibleInterval, DeviceGeometry, LiquidProperties, ModelError, SolverSettings,
                    SolverStatus)
from .thermal import SlabMesh, SlabSolution, solve_slab
from .transport import Closure, DarcyHydraulics, TableRangeError

ROUNDOFF_REL = 1e-12
GAUSS_X, GAUSS_W = np.polynomial.legendre.leggauss(8)


@dataclass
class CoupledProblem:
    geometry: DeviceGeometry
    liquid: LiquidProperties
    closure: Closure
    heat_shape: Callable[[np.ndarray], np.ndarray] = lambda xi: np.ones_like(xi)
    heat_shape_label: str = "uniform"
    Pc_feed: float = 0.0
    Pc_physical_limit: float | None = None
    physical_limit_criterion: str = ""

    def shape_cells(self, mesh: SlabMesh) -> np.ndarray:
        """Cell averages of the heat-input shape, normalized to unit mean over [0, L]."""
        f = self.x_cell_average(mesh, self.heat_shape)
        mean = float(np.sum(f * mesh.dx) / self.geometry.L)
        if not mean > 0:
            raise ModelError("heat-input shape must have a positive mean")
        return f / mean

    def x_cell_average(self, mesh: SlabMesh, fn) -> np.ndarray:
        xf = mesh.x_faces
        a, b = xf[:-1, None], xf[1:, None]
        xg = 0.5 * (a + b) + 0.5 * (b - a) * GAUSS_X[None, :]
        vals = np.asarray(fn(xg / self.geometry.L), float)
        return 0.5 * np.sum(vals * GAUSS_W[None, :], axis=1)


@dataclass
class CoupledSolution:
    status: SolverStatus
    converged: bool
    reason: str
    load: float  # mean input heat flux over the half-device footprint [W/m^2]
    iterations: int
    history: list[dict]
    mesh: SlabMesh | None = None
    q_in: np.ndarray | None = None
    q_e: np.ndarray | None = None
    h_e: np.ndarray | None = None
    Q_faces: np.ndarray | None = None
    Pc_faces: np.ndarray | None = None
    Pc_centers: np.ndarray | None = None
    slab: SlabSolution | None = None
    interval: AdmissibleInterval | None = None
    settings: dict = field(default_factory=dict)

    @property
    def Pc_end(self) -> float:
        return float(self.Pc_faces[-1])


def admissible_interval(problem: CoupledProblem, settings: SolverSettings) -> AdmissibleInterval:
    """Pressure interval in which every closure evaluation is supported, and what ends it."""
    c = problem.closure
    lo, hi = c.Pc_lo, c.Pc_hi
    snap = c.snap
    notes = []
    if not (lo - snap <= problem.Pc_feed <= hi + snap):
        raise ModelError(f"feed suction {problem.Pc_feed} Pa outside the supported interval [{lo}, {hi}]")
    phys = problem.Pc_physical_limit
    if phys is None:
        phys = c.limits.Pc_physical_limit
        if phys is not None:
            notes.append("physical limit taken from closure metadata")
    geom = c.limits.Pc_geometry_limit if c.limits.geometry_limit_demonstrated else None
    if c.limits.Pc_geometry_limit is not None and not c.limits.geometry_limit_demonstrated:
        notes.append("a meniscus geometry limit was declared but not demonstrated; it is not used")
    if geom is not None and (phys is None or geom < phys) and geom <= hi + snap:
        cap, mech = min(geom, hi), SolverStatus.MENISCUS_GEOMETRY_LIMIT
    elif phys is not None and phys <= hi + snap:
        cap, mech = min(phys, hi), SolverStatus.PHYSICAL_CAPILLARY_LIMIT
    else:
        cap, mech = hi, SolverStatus.TABLE_COVERAGE_LIMIT
        notes.append("table coverage ends below the physical limit: incomplete coverage, not dryout")
    if cap <= problem.Pc_feed:
        raise ModelError("admissible interval is empty above the feed suction")
    return AdmissibleInterval(lo, cap, mech, phys, geom, hi, notes)


def _scaled_change(new, old, floor):
    new = np.asarray(new, float)
    if old is None:
        return float("inf")
    return float(np.linalg.norm(new - old) / (np.linalg.norm(new) + np.sqrt(new.size) * floor))


def _hydraulics(problem: CoupledProblem) -> DarcyHydraulics:
    g, f = problem.geometry, problem.liquid
    return DarcyHydraulics(problem.closure, f.mu, g.W, f.rho, f.h_fg)


def _march(hyd: DarcyHydraulics, mesh, q, problem, cap):
    st = hyd.march(mesh.x_faces, q, problem.Pc_feed, cap)
    if st.exceeds:
        psi_cap = st.psi_cap
        if st.psi_required_end <= psi_cap + ROUNDOFF_REL * max(abs(psi_cap), 1e-300):
            # roundoff only: clamp to the cap through the inverse
            psi_f = np.minimum(st.psi_faces, psi_cap)
            psi_c = np.minimum(st.psi_centers, psi_cap)
            st.Pc_faces = hyd.closure.psi_inverse(psi_f)
            st.Pc_centers = hyd.closure.psi_inverse(psi_c)
            st.exceeds = False
    return st


def solve_at_load(problem: CoupledProblem, load: float, settings: SolverSettings,
                  interval: AdmissibleInterval | None = None, q_e_guess=None,
                  _inject_failure: bool = False) -> CoupledSolution:
    """Forward solve at a prescribed mean input flux ``load`` [W/m^2]."""
    g = problem.geometry
    mesh = SlabMesh(g.L, g.t_s, settings.Nx, settings.Nz)
    interval = interval or admissible_interval(problem, settings)
    hyd = _hydraulics(problem)
    if not np.isfinite(load):
        return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False, "nonfinite load", load, 0, [])
    q_in = load * problem.shape_cells(mesh) if load != 0 else np.zeros(mesh.Nx)
    if np.any(q_in < 0):
        return CoupledSolution(SolverStatus.OUTSIDE_MODEL_REGIME, False,
                               "negative heat input implies a condensing branch, which this passive "
                               "evaporating model does not support", load, 0, [], mesh=mesh, interval=interval)
    q_char = max(float(np.mean(np.abs(q_in))), 1.0)
    Pc_char = interval.Pc_cap - interval.Pc_lo
    q_floor, Pc_floor = settings.abs_floor_fraction * q_char, settings.abs_floor_fraction * Pc_char

    q_valid = np.zeros(mesh.Nx)
    q_prop = q_in.copy() if q_e_guess is None else np.asarray(q_e_guess, float).copy()
    Pc_old = None
    history: list[dict] = []
    for it in range(settings.max_iter):
        omega = 1.0
        while True:
            q_try = q_valid + omega * (q_prop - q_valid)
            st = _march(hyd, mesh, q_try, problem, interval.Pc_cap)
            if not st.exceeds:
                break
            omega *= 0.5
            if omega < settings.min_relaxation:
                return CoupledSolution(
                    SolverStatus.NUMERICAL_FAILURE, False,
                    "iterates repeatedly required suction beyond the admissible interval; this is not "
                    "treated as a physical bracket", load, it, history, mesh=mesh, interval=interval)
        try:
            h_e = problem.closure.h_e(st.Pc_centers)
            if _inject_failure:
                h_e = h_e * np.nan
            slab = solve_slab(mesh, g.k_s, q_in, h_e)
        except (TableRangeError, ValueError, np.linalg.LinAlgError) as exc:
            return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False, f"closure/conduction error: {exc}",
                                   load, it, history, mesh=mesh, interval=interval)
        q_new = slab.q_e
        if not np.all(np.isfinite(slab.theta)):
            return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False, "nonfinite temperature field",
                                   load, it, history, mesh=mesh, interval=interval)
        if np.any(q_new < -q_floor):
            return CoupledSolution(SolverStatus.OUTSIDE_MODEL_REGIME, False,
                                   "conduction returned condensation (q_e < 0); unsupported branch, not clipped",
                                   load, it, history, mesh=mesh, interval=interval)
        dq = _scaled_change(q_new, q_try, q_floor)
        dPc = _scaled_change(st.Pc_centers, Pc_old, Pc_floor)
        history.append({"iteration": it, "dq_scaled": dq, "dPc_scaled": dPc, "backtrack_omega": omega,
                        "Pc_end_Pa": float(st.Pc_faces[-1]), "linear_residual": slab.linear_residual})
        if dq <= settings.rtol and dPc <= settings.rtol:
            final = _march(hyd, mesh, q_new, problem, interval.Pc_cap)
            if final.exceeds:
                return CoupledSolution(interval.mechanism, False,
                                       "converged evaporation field requires suction beyond the admissible "
                                       "cap: load is at or above the limit", load, it + 1, history,
                                       mesh=mesh, interval=interval)
            return CoupledSolution(SolverStatus.ADMISSIBLE, True, "converged", load, it + 1, history, mesh,
                                   q_in, q_new, h_e, final.Q_faces, final.Pc_faces, final.Pc_centers, slab,
                                   interval, settings.as_dict())
        q_valid = q_try
        Pc_old = st.Pc_centers
        q_prop = q_try + settings.relaxation * (q_new - q_try)
    return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False,
                           f"no convergence in {settings.max_iter} iterations", load, settings.max_iter, history,
                           mesh=mesh, interval=interval)


def solve_at_endpoint(problem: CoupledProblem, Pc_end: float, settings: SolverSettings,
                      interval: AdmissibleInterval | None = None) -> CoupledSolution:
    """Solve for the load at which Pc(L) = Pc_end, with all closure evaluations in [Pc_feed, Pc_end]."""
    g, f, c = problem.geometry, problem.liquid, problem.closure
    mesh = SlabMesh(g.L, g.t_s, settings.Nx, settings.Nz)
    interval = interval or admissible_interval(problem, settings)
    if not (problem.Pc_feed < Pc_end <= interval.Pc_cap + c.snap):
        raise ModelError("endpoint suction must lie in (Pc_feed, Pc_cap]")
    hyd = _hydraulics(problem)
    shape = problem.shape_cells(mesh)
    psi_target = float(c.psi(Pc_end)[0] - c.psi(problem.Pc_feed)[0])
    xc, dx = mesh.xc, mesh.dx
    Pc_char = interval.Pc_cap - interval.Pc_lo
    Pc_floor = settings.abs_floor_fraction * Pc_char
    q_unit = shape.copy()  # first iteration: uniform evaporation (V2019 A.3)
    Pc_old, s_old = None, None
    history: list[dict] = []
    for it in range(settings.max_iter):
        moment = float(np.sum(xc * q_unit * dx)) / (f.rho * f.h_fg)
        if not moment > 0:
            return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False, "nonpositive evaporation moment",
                                   float("nan"), it, history, mesh=mesh, interval=interval)
        s = psi_target / (f.mu * moment)
        st = _march(hyd, mesh, s * q_unit, problem, Pc_end)
        if st.exceeds:
            return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False,
                                   "endpoint construction exceeded the target beyond roundoff", s, it, history,
                                   mesh=mesh, interval=interval)
        try:
            h_e = c.h_e(st.Pc_centers)
            slab_u = solve_slab(mesh, g.k_s, shape, h_e)
        except (TableRangeError, ValueError) as exc:
            return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False, f"closure/conduction error: {exc}",
                                   s, it, history, mesh=mesh, interval=interval)
        q_unit_new = slab_u.q_e
        dq = _scaled_change(q_unit_new, q_unit, settings.abs_floor_fraction)
        dPc = _scaled_change(st.Pc_centers, Pc_old, Pc_floor)
        ds = abs(s - s_old) / abs(s) if s_old is not None else float("inf")
        history.append({"iteration": it, "dq_scaled": dq, "dPc_scaled": dPc, "dload_rel": ds, "load_W_m2": s,
                        "linear_residual": slab_u.linear_residual})
        if dq <= settings.rtol and dPc <= settings.rtol and ds <= settings.rtol:
            # exact linear rescaling so the returned q_e has exactly the target evaporation moment
            s_star = psi_target / (f.mu * float(np.sum(xc * q_unit_new * dx)) / (f.rho * f.h_fg))
            q_in = s_star * shape
            slab = solve_slab(mesh, g.k_s, q_in, h_e)
            final = _march(hyd, mesh, slab.q_e, problem, Pc_end)
            if final.exceeds:
                return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False,
                                       "final endpoint state exceeded the target", s, it + 1, history,
                                       mesh=mesh, interval=interval)
            sol = CoupledSolution(SolverStatus.ADMISSIBLE, True, "converged (endpoint-constrained)", s_star,
                                  it + 1, history, mesh, q_in, slab.q_e, h_e, final.Q_faces, final.Pc_faces,
                                  final.Pc_centers, slab, interval, settings.as_dict())
            sol.settings["endpoint_load_consistency_rel"] = abs(s_star - s) / s_star
            return sol
        q_unit = q_unit + settings.relaxation * (q_unit_new - q_unit)
        Pc_old, s_old = st.Pc_centers, s
    return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False,
                           f"endpoint solve did not converge in {settings.max_iter} iterations",
                           float("nan"), settings.max_iter, history, mesh=mesh, interval=interval)


def solve_prescribed_load(problem: CoupledProblem, load: float, settings: SolverSettings) -> CoupledSolution:
    """Forward solve at a prescribed load, classifying loads at or above the limit by mechanism.

    The endpoint-constrained solve locates the load at which Pc(L) reaches the admissible cap. A
    prescribed load at or above it is reported with the cap's mechanism (a physical or geometry limit,
    or TABLE_COVERAGE_LIMIT), not as a numerical failure of the forward iteration.
    """
    interval = admissible_interval(problem, settings)
    if load > 0:
        end = solve_at_endpoint(problem, interval.Pc_cap, settings, interval)
        if not end.converged:
            return CoupledSolution(SolverStatus.NUMERICAL_FAILURE, False, f"endpoint solve: {end.reason}", load, 0,
                                   [], interval=interval)
        if load >= end.load:
            return CoupledSolution(interval.mechanism, False,
                                   f"prescribed load {load:.6g} W/m^2 is at or above the load {end.load:.6g} W/m^2 "
                                   f"at which Pc(L) reaches the admissible cap ({interval.mechanism.value})",
                                   load, 0, [], interval=interval)
    return solve_at_load(problem, load, settings, interval)


@dataclass
class DryoutResult:
    status: SolverStatus
    mechanism: SolverStatus | None
    reason: str
    interval: AdmissibleInterval | None
    endpoint: CoupledSolution | None
    lower_operating_point: float | None  # certified converged admissible load [W/m^2]
    upper_bracket: float | None  # load at which Pc(L) = Pc_cap [W/m^2]; None if unresolved
    relative_interval: float | None
    continuation: list[dict]
    monotonic_end_suction: bool | None
    forward_endpoint_consistency: dict
    regime: dict | None = None

    def summary(self) -> dict:
        return {
            "status": self.status.value,
            "mechanism": self.mechanism.value if self.mechanism else None,
            "reason": self.reason,
            "admissible_interval": self.interval.as_dict() if self.interval else None,
            "lower_operating_point_W_m2": self.lower_operating_point,
            "upper_bracket_W_m2": self.upper_bracket,
            "relative_interval": self.relative_interval,
            "monotonic_end_suction_along_continuation": self.monotonic_end_suction,
            "forward_endpoint_consistency": self.forward_endpoint_consistency,
            "continuation": self.continuation,
            "regime": self.regime,
        }


def search_dryout(problem: CoupledProblem, settings: SolverSettings, regime_check=None,
                  _inject_failure: bool = False) -> DryoutResult:
    """Locate the limiting load inside the admissible interval (see module docstring)."""
    try:
        interval = admissible_interval(problem, settings)
    except ModelError as exc:
        return DryoutResult(SolverStatus.OUTSIDE_MODEL_REGIME, None, str(exc), None, None, None, None, None,
                            [], None, {})
    end = solve_at_endpoint(problem, interval.Pc_cap, settings, interval)
    if not end.converged:
        return DryoutResult(SolverStatus.NUMERICAL_FAILURE, interval.mechanism, f"endpoint solve: {end.reason}",
                            interval, end, None, None, None, [], None, {})
    s_end = end.load
    loads = [fr * s_end for fr in settings.continuation_fractions] + [(1 - settings.bracket_rel_width) * s_end]
    cont, guess, prev = [], None, None
    for i, s in enumerate(loads):
        sol = solve_at_load(problem, s, settings, interval, q_e_guess=guess,
                            _inject_failure=_inject_failure and i == len(loads) - 1)
        cont.append({"load_W_m2": s, "fraction_of_endpoint": s / s_end, "status": sol.status.value,
                     "converged": sol.converged, "iterations": sol.iterations,
                     "Pc_end_Pa": sol.Pc_end if sol.converged else None})
        if not sol.converged:
            return DryoutResult(SolverStatus.NUMERICAL_FAILURE, interval.mechanism,
                                f"continuation failed at {s:.6g} W/m^2: {sol.reason}; the last converged "
                                "admissible point is a lower bound only", interval, end,
                                prev.load if prev else None, None, None, cont, None, {})
        if i + 1 < len(loads):
            guess = sol.q_e * (loads[i + 1] / s)
        prev = sol
    ends = [c["Pc_end_Pa"] for c in cont]
    monotonic = bool(np.all(np.diff(ends) > 0))
    last = prev
    consistency = {
        "endpoint_Pc_end_Pa": end.Pc_end,
        "forward_Pc_end_Pa_at_lower_bracket": last.Pc_end,
        "forward_below_cap": bool(last.Pc_end < interval.Pc_cap),
        "q_e_rel_diff_forward_vs_endpoint": float(np.linalg.norm(last.q_e - end.q_e) / np.linalg.norm(end.q_e)),
        "Pc_rel_diff_forward_vs_endpoint": float(np.linalg.norm(last.Pc_faces - end.Pc_faces)
                                                 / np.linalg.norm(end.Pc_faces)),
        "endpoint_load_consistency_rel": end.settings.get("endpoint_load_consistency_rel"),
    }
    lower = last.load
    rel = (s_end - lower) / s_end
    regime = regime_check(end) if regime_check else None
    if not monotonic or not consistency["forward_below_cap"]:
        return DryoutResult(SolverStatus.NUMERICAL_FAILURE, interval.mechanism,
                            "end suction was not monotonic along continuation or the subcritical forward solve "
                            "reached the cap; the endpoint is not certified", interval, end, None, None, None,
                            cont, monotonic, consistency, regime)
    if regime and not regime.get("within_supported_ranges", True):
        return DryoutResult(SolverStatus.OUTSIDE_MODEL_REGIME, interval.mechanism,
                            "solution leaves the closure's supported operating ranges: " + "; ".join(regime["violations"]),
                            interval, end, lower, s_end, rel, cont, monotonic, consistency, regime)
    if interval.mechanism == SolverStatus.TABLE_COVERAGE_LIMIT:
        return DryoutResult(SolverStatus.TABLE_COVERAGE_LIMIT, interval.mechanism,
                            "table coverage ends below the physical limit; the load at the coverage end is a "
                            "certified lower operating point and the physical limit is unresolved above it",
                            interval, end, s_end, None, None, cont, monotonic, consistency, regime)
    return DryoutResult(interval.mechanism, interval.mechanism, "limit located", interval, end, lower, s_end, rel,
                        cont, monotonic, consistency, regime)
