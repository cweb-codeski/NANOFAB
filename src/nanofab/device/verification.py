"""Numerical verification suites with predeclared gates (docs/passive_device_model.md, Sec. 7).

Every suite returns ``{"suite", "passed", "gates": [...], "details": {...}}``. A gate records the
metric, the observed value, the predeclared limit and the outcome; gates are never relaxed after
seeing results. All inputs here are ANALYTICAL_EXAMPLE or MANUFACTURED_FIXTURE: results are
"numerically verified for these cases" only and say nothing about physical accuracy.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

from ..provenance import REPO_ROOT
from . import analytic
from .case import load_case
from .diagnostics import energy_balance, mass_balance, reference_scales, regime_diagnostics
from .model import EXIT_CODES, DeviceGeometry, SolverSettings, SolverStatus
from .solver import CoupledProblem, search_dryout, solve_at_endpoint, solve_at_load
from .tables import load_dataset, match_case
from .thermal import SlabMesh, solve_slab
from .transport import ConstantClosure, DarcyHydraulics

FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "device"
FIXTURE_CASE = FIXTURE_DIR / "manufactured_case.yaml"
FIXTURE_TABLES = FIXTURE_DIR / "manufactured_tables"
MANUFACTURED_ENDPOINT_PA = 2200.0

GATES = {
    "closed_form_rel": 1e-9,
    "w_invariance_rel": 1e-10,
    "zero_load_abs": 1e-12,
    "mass_global": 1e-8,
    "mass_local": 1e-7,
    "energy_global": 1e-8,
    "energy_local": 1e-7,
    "hydraulic_same_input_rel": 1e-9,
    "observed_order_min": 1.8,
    "observed_order_max": 2.2,
    "thermal_finest_abs_K": 1e-3,
    "tolerance_rel": 1e-4,
    "tolerance_abs_K": 1e-3,
    "mesh_finest_pair_rel": 2e-3,
    "mesh_finest_pair_abs_K": 0.02,
    "limit_interval_rel": 1e-4,
    "moment_identity_rel": 1e-9,
}
REFINEMENT = ((50, 8), (100, 16), (200, 32))
AUDIT_SCALES_82 = {"U_feed_m_s": 0.07584, "u_pore_volume_porosity_m_s": 0.08309, "Re_h_superficial": 6.45,
                   "Re_h_pore": 7.07, "Pe_vertical": 0.0566, "Pe_axial_L_pore": 2479.0, "Ja": 0.0187}


def gate(name, value, limit, kind="max") -> dict:
    if kind == "max":
        ok = bool(np.isfinite(value) and value <= limit)
    elif kind == "min":
        ok = bool(np.isfinite(value) and value >= limit)
    elif kind == "range":
        ok = bool(np.isfinite(value) and limit[0] <= value <= limit[1])
    else:  # "equal"
        ok = value == limit
    return {"gate": name, "value": value, "limit": limit, "kind": kind, "passed": ok}


def _suite(name, gates, details) -> dict:
    return {"suite": name, "passed": all(g["passed"] for g in gates), "gates": gates, "details": details}


def _fixture_module():
    spec = importlib.util.spec_from_file_location("nanofab_fixture_tables", FIXTURE_DIR / "fixture_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def base_case():
    return load_case(FIXTURE_CASE)


def fixture_problem(workdir: Path | None = None, variant: str | None = None, W: float | None = None,
                    Pc_physical_limit: float = MANUFACTURED_ENDPOINT_PA, **build_kw) -> CoupledProblem:
    """Problem on the V2019 geometry with MANUFACTURED_FIXTURE tables (VERIFICATION ONLY)."""
    case = base_case()
    if variant is None:
        path = FIXTURE_TABLES
    else:
        path = Path(workdir) / variant
        if not (path / "manifest.yaml").exists():
            _fixture_module().write_fixture(path, **build_kw)
    ds = load_dataset(path, repo_root=REPO_ROOT)
    match_case(ds, case.case_requirements())
    if not ds.report.schema_valid or not ds.report.case_matched:
        raise RuntimeError(f"fixture {path} rejected: {ds.report.as_dict()}")
    g = case.geometry
    if W is not None:
        g = DeviceGeometry(L=g.L, W=W, h=g.h, t_s=g.t_s, k_s=g.k_s, d=g.d, l=g.l)
    return CoupledProblem(g, case.liquid, ds.to_closure(1e-9), Pc_physical_limit=Pc_physical_limit,
                          physical_limit_criterion="MANUFACTURED declared endpoint")


def constant_problem(K=3.0e-11, h_e=6.0e4) -> tuple[CoupledProblem, float]:
    """ANALYTICAL_EXAMPLE: constant K and h_e on the V2019 geometry with the Eq. (3) limit."""
    case = base_case()
    Pc_lim = case.Pc_limit_eq3
    clo = ConstantClosure(K, case.geometry.h, h_e, Pc_hi=Pc_lim)
    return CoupledProblem(case.geometry, case.liquid, clo, Pc_physical_limit=Pc_lim,
                          physical_limit_criterion="V2019 Eq. (3)"), Pc_lim


def _balance_gates(problem, sol, tag) -> list[dict]:
    mb, eb = mass_balance(problem, sol), energy_balance(problem, sol)
    return [gate(f"{tag}: mass global (normalized)", mb["global_normalized"], GATES["mass_global"]),
            gate(f"{tag}: mass local max (normalized)", mb["local_max_normalized"], GATES["mass_local"]),
            gate(f"{tag}: Darcy cell consistency (independent quadrature)", mb["darcy_cell_max_normalized"],
                 GATES["mass_local"]),
            gate(f"{tag}: energy global (normalized)", eb["global_normalized"], GATES["energy_global"]),
            gate(f"{tag}: energy local max (normalized)", eb["local_max_normalized"], GATES["energy_local"])]


# --------------------------------------------------------------------------------------------- analytic
def suite_analytic(Nx=200, Nz=32) -> dict:
    p, Pc_lim = constant_problem()
    g, f = p.geometry, p.liquid
    s = SolverSettings(Nx=Nx, Nz=Nz)
    res = search_dryout(p, s)
    q_ref = analytic.dryout_constant(p.closure.K, g.h, g.L, f.rho, f.mu, f.h_fg, Pc_lim)
    end = res.endpoint
    dTe_ref, dTs_ref = analytic.superheats_constant(end.load, p.closure.he, g.t_s, g.k_s)
    pg_ref = analytic.pressure_profile_pg(end.mesh.x_faces, end.load, g.L, p.closure.K, g.h, f.rho, f.mu, f.h_fg)
    Q_ref = analytic.flow_profile(end.mesh.x_faces, end.load, g.L, g.W, f.rho, f.h_fg)
    gates = [
        gate("status is PHYSICAL_CAPILLARY_LIMIT", res.status.value, SolverStatus.PHYSICAL_CAPILLARY_LIMIT.value, "equal"),
        gate("q_dry vs closed form (rel)", abs(end.load - q_ref) / q_ref, GATES["closed_form_rel"]),
        gate("p_g(x) vs closed form (rel to Pc_limit)", float(np.max(np.abs(-end.Pc_faces - pg_ref))) / Pc_lim,
             GATES["closed_form_rel"]),
        gate("Q(x) vs closed form (rel to Q(0))", float(np.max(np.abs(end.Q_faces - Q_ref))) / Q_ref[0],
             GATES["closed_form_rel"]),
        gate("T_e - T_v vs q/h_e (rel)", float(np.max(np.abs(end.slab.theta_e - dTe_ref))) / dTe_ref,
             GATES["closed_form_rel"]),
        gate("T_source - T_v vs q/h_e + q t_s/k_s (rel)", float(np.max(np.abs(end.slab.theta_source - dTs_ref))) / dTs_ref,
             GATES["closed_form_rel"]),
        gate("limit interval (rel)", res.relative_interval, GATES["limit_interval_rel"]),
    ]
    # pressure-dependent transmissivity, constant h_e: uniform evaporation is exact, so the quadrature
    # closed form applies to the tabulated hydraulics
    pf = fixture_problem()
    tab = pf.closure
    he_const = ConstantClosure(1.0, g.h, 6.0e4, Pc_hi=tab.Pc_hi)
    mixed = _MixedClosure(tab, he_const)
    pm = CoupledProblem(g, f, mixed, Pc_physical_limit=MANUFACTURED_ENDPOINT_PA, physical_limit_criterion="manufactured")
    endm = solve_at_endpoint(pm, MANUFACTURED_ENDPOINT_PA, s)
    q_tau = analytic.dryout_from_transmissivity(lambda P: float(tab.tau(P)[0]), f.mu, f.rho, f.h_fg, g.L,
                                                MANUFACTURED_ENDPOINT_PA, breakpoints=tab._tau.x)
    gates.append(gate("q_dry with tabulated tau(Pc) vs quadrature closed form (rel)", abs(endm.load - q_tau) / q_tau,
                      GATES["closed_form_rel"]))
    gates += _balance_gates(p, end, "constant closure endpoint")
    return _suite("analytic", gates, {"q_dry_numeric_W_m2": end.load, "q_dry_closed_form_W_m2": q_ref,
                                      "q_dry_tau_numeric_W_m2": endm.load, "q_dry_tau_quadrature_W_m2": q_tau,
                                      "Pc_limit_Pa": Pc_lim, "label": "ANALYTICAL_EXAMPLE (K = 3.0e-11 m^2, "
                                      "h_e = 6.0e4 W/m^2K chosen for verification, not fitted)"})


class _MixedClosure:
    """tau from one closure, constant h_e from another (verification helper)."""

    def __init__(self, hyd, th):
        self._h, self._t = hyd, th
        self.limits, self.snap, self.name = hyd.limits, hyd.snap, "mixed"
        self.Pc_lo, self.Pc_hi = hyd.Pc_lo, hyd.Pc_hi

    def tau(self, P):
        return self._h.tau(P)

    def psi(self, P):
        return self._h.psi(P)

    def psi_inverse(self, S):
        return self._h.psi_inverse(S)

    def h_e(self, P):
        return np.full_like(self._h.tau(P), self._t.he)

    def liquid_section_area_per_width(self, P):
        return None


# --------------------------------------------------------------------------------------------- symmetry, W
def suite_accounting(Nx=200, Nz=32) -> dict:
    s = SolverSettings(Nx=Nx, Nz=Nz)
    pa = fixture_problem()
    pb = fixture_problem(W=pa.geometry.W / 5)
    ea, eb = solve_at_endpoint(pa, MANUFACTURED_ENDPOINT_PA, s), solve_at_endpoint(pb, MANUFACTURED_ENDPOINT_PA, s)
    ratio = pa.geometry.W / pb.geometry.W
    g, f = pa.geometry, pa.liquid
    half_in = float(np.sum(ea.q_in * ea.mesh.dx)) * g.W
    half_ev = float(np.sum(ea.q_e * ea.mesh.dx)) * g.W
    full_in, full_ev = 2 * half_in, 2 * half_ev
    feed_total = 2 * ea.Q_faces[0] * f.rho * f.h_fg
    gates = [
        gate("W invariance: dryout load (rel)", abs(ea.load - eb.load) / ea.load, GATES["w_invariance_rel"]),
        gate("W invariance: Pc(x) (rel to cap)", float(np.max(np.abs(ea.Pc_faces - eb.Pc_faces))) / MANUFACTURED_ENDPOINT_PA,
             GATES["w_invariance_rel"]),
        gate("W invariance: T_e(x) (rel to max superheat)",
             float(np.max(np.abs(ea.slab.theta_e - eb.slab.theta_e)) / np.max(ea.slab.theta_e)), GATES["w_invariance_rel"]),
        gate("W scaling: Q(x) scales with W (rel)", float(np.max(np.abs(ea.Q_faces - ratio * eb.Q_faces)) / ea.Q_faces[0]),
             GATES["w_invariance_rel"]),
        gate("mirrored device: total input = 2 x half (evaporation) (rel)", abs(full_in - full_ev) / full_in,
             GATES["mass_global"]),
        gate("mirrored device: two edge feeds carry all evaporation (rel)", abs(feed_total - full_ev) / full_ev,
             GATES["mass_global"]),
        gate("mirrored device: mean flux over 2L W equals half-device mean (rel)",
             abs(full_in / g.footprint_full - half_in / g.footprint_half) / (half_in / g.footprint_half),
             GATES["closed_form_rel"]),
    ]
    return _suite("accounting", gates, {"W_m": [pa.geometry.W, pb.geometry.W], "load_W_m2": [ea.load, eb.load],
                                        "full_device_input_W": full_in, "half_device_input_W": half_in,
                                        "note": "a mirrored half-device has zero flow and zero axial heat flux "
                                                "at x = L by construction (Q(L) = 0, insulated end)"})


# --------------------------------------------------------------------------------------------- zero load
def suite_zero_load(Nx=100, Nz=16) -> dict:
    p = fixture_problem()
    sol = solve_at_load(p, 0.0, SolverSettings(Nx=Nx, Nz=Nz))
    gates = [gate("zero load converged", sol.converged, True, "equal")]
    if sol.converged:
        gates += [gate("zero load: max |Pc - Pc_feed| [Pa]", float(np.max(np.abs(sol.Pc_faces - p.Pc_feed))),
                       GATES["zero_load_abs"]),
                  gate("zero load: max |Q| [m^3/s]", float(np.max(np.abs(sol.Q_faces))), GATES["zero_load_abs"]),
                  gate("zero load: max |T - T_v| [K]", float(np.max(np.abs(sol.slab.theta))), GATES["zero_load_abs"])]
        gates += _balance_gates(p, sol, "zero load")
    return _suite("zero_load", gates, {"status": sol.status.value, "iterations": sol.iterations})


# --------------------------------------------------------------------------------------------- hydraulics
def _q_cont(x, q0, L):
    return q0 * (1 + 0.5 * np.cos(np.pi * x / L))


def _q_cell_avg(xf, q0, L):
    a, b = xf[:-1], xf[1:]
    return q0 * (1 + 0.5 * (L / np.pi) * (np.sin(np.pi * b / L) - np.sin(np.pi * a / L)) / (b - a))


def _psi_quad(closure, P):
    knots = [k for k in closure._tau.x if closure.Pc_lo < k < P]
    return integrate.quad(lambda p: float(closure.tau(p)[0]), closure.Pc_lo, P, points=knots or None,
                          limit=500, epsabs=0.0, epsrel=1e-13)[0]


def _invert_psi_quad(closure, target):
    if target <= 0:
        return closure.Pc_lo
    return optimize.brentq(lambda P: _psi_quad(closure, P) - target, closure.Pc_lo, closure.Pc_hi,
                           xtol=1e-12, rtol=1e-14)


def _ode_pressure(c, mu, W, Q_of_x, x_eval, breaks):
    """Independent reference: integrate dPc/dx = mu Q(x) / (W tau(Pc)) in the primitive variable
    (no Kirchhoff transform), segment by segment between flow-law breakpoints (DOP853)."""
    def rhs(x, P):
        return [mu * Q_of_x(x) / (W * float(c.tau(min(P[0], c.Pc_hi))[0]))]
    out = np.empty(len(x_eval))
    out[0] = c.Pc_lo
    P = c.Pc_lo
    for i in range(len(x_eval) - 1):
        a, b = x_eval[i], x_eval[i + 1]
        sol = integrate.solve_ivp(rhs, (a, b), [P], method="DOP853", rtol=1e-13, atol=1e-11)
        P = float(sol.y[0, -1])
        out[i + 1] = P
    return out


def suite_hydraulics(levels=(50, 100, 200), q0=6.0e5) -> dict:
    p = fixture_problem()
    g, f, c = p.geometry, p.liquid, p.closure
    hyd = DarcyHydraulics(c, f.mu, g.W, f.rho, f.h_fg)
    errs_same, errs_cont, rows = [], [], []
    for Nx in levels:
        xf = np.linspace(0, g.L, Nx + 1)
        qc = _q_cell_avg(xf, q0, g.L)
        st = hyd.march(xf, qc, 0.0, c.Pc_hi)
        if st.exceeds:
            raise RuntimeError("hydraulic manufactured load exceeds the fixture range")
        # (i) same piecewise-constant input (Q piecewise linear in x), independent ODE integration
        Q_pc = hyd.face_flows(xf, qc)
        P_same = _ode_pressure(c, f.mu, g.W, lambda x: float(np.interp(x, xf, Q_pc)), xf, xf)
        e_same = float(np.max(np.abs(st.Pc_faces - P_same)) / c.Pc_hi)
        # (ii) continuous q(x): exact Q(x), independent ODE integration
        Qx = lambda x: g.W * q0 * ((g.L - x) - 0.5 * (g.L / np.pi) * np.sin(np.pi * x / g.L)) / (f.rho * f.h_fg)  # noqa: E731
        P_cont = _ode_pressure(c, f.mu, g.W, Qx, xf, xf)
        e_cont = float(np.max(np.abs(st.Pc_faces - P_cont)))
        errs_same.append(e_same)
        errs_cont.append(e_cont)
        rows.append({"Nx": Nx, "err_same_input_rel": e_same, "err_continuous_q_Pa": e_cont,
                     "Pc_end_Pa": float(st.Pc_faces[-1])})
    orders = [math.log(errs_cont[i] / errs_cont[i + 1]) / math.log(levels[i + 1] / levels[i]) for i in range(len(levels) - 1)]
    gates = [gate("same piecewise-constant input vs independent ODE integration (rel, all levels)", max(errs_same),
                  GATES["hydraulic_same_input_rel"]),
             gate("continuous q(x): observed order (finest pair)", orders[-1],
                  (GATES["observed_order_min"], GATES["observed_order_max"]), "range")]
    return _suite("hydraulics", gates, {"levels": rows, "observed_orders": orders, "q0_W_m2": q0,
                                        "q_shape": "q0 (1 + 0.5 cos(pi x / L))", "tables": "MANUFACTURED_FIXTURE"})


# --------------------------------------------------------------------------------------------- thermal
MMS = {"L": 1e-3, "t_s": 0.25e-3, "k_s": 148.0, "A": 5.0, "B": 3378.0, "C": 0.04, "m": 2}


def mms_fields(x, zeta, P=MMS):
    """theta(x, zeta) = A + B (t_s - zeta) + C cos(k x) exp(k zeta), zeta = z + t_s (harmonic)."""
    k = P["m"] * math.pi / P["L"]
    return P["A"] + P["B"] * (P["t_s"] - zeta) + P["C"] * np.cos(k * x) * np.exp(k * zeta)


def mms_boundary(P=MMS):
    k = P["m"] * math.pi / P["L"]
    ks, t = P["k_s"], P["t_s"]
    q_in = lambda x: ks * P["B"] - ks * P["C"] * k * np.cos(k * x)  # noqa: E731
    q_e = lambda x: ks * P["B"] - ks * P["C"] * k * math.exp(k * t) * np.cos(k * x)  # noqa: E731
    th_top = lambda x: P["A"] + P["C"] * math.exp(k * t) * np.cos(k * x)  # noqa: E731
    th_bot = lambda x: P["A"] + P["B"] * t + P["C"] * np.cos(k * x)  # noqa: E731
    return q_in, q_e, th_top, th_bot, k


def suite_thermal(levels=REFINEMENT) -> dict:
    P = MMS
    q_in_f, q_e_f, th_top, th_bot, k = mms_boundary(P)
    rows, e_top, e_bot = [], [], []
    for Nx, Nz in levels:
        m = SlabMesh(P["L"], P["t_s"], Nx, Nz)
        xf = m.x_faces
        a, b = xf[:-1], xf[1:]
        q_avg = P["k_s"] * P["B"] - P["k_s"] * P["C"] * (np.sin(k * b) - np.sin(k * a)) / (b - a)
        h_e = q_e_f(m.xc) / th_top(m.xc)
        sol = solve_slab(m, P["k_s"], q_avg, h_e)
        et = float(np.max(np.abs(sol.theta_e - th_top(m.xc))))
        eb = float(np.max(np.abs(sol.theta_source - th_bot(m.xc))))
        X, Z = np.meshgrid(m.xc, m.zc + P["t_s"], indexing="ij")
        ef = float(np.max(np.abs(sol.theta - mms_fields(X, Z))))
        E_in, E_out = float(np.sum(q_avg * m.dx)), float(np.sum(sol.q_e * m.dx))
        e_top.append(et)
        e_bot.append(eb)
        rows.append({"Nx": Nx, "Nz": Nz, "err_T_e_K": et, "err_T_source_K": eb, "err_field_K": ef,
                     "global_energy_rel": abs(E_in - E_out) / E_in, "linear_residual": sol.linear_residual})
    o_top = [math.log(e_top[i] / e_top[i + 1]) / math.log(2) for i in range(len(levels) - 1)]
    o_bot = [math.log(e_bot[i] / e_bot[i + 1]) / math.log(2) for i in range(len(levels) - 1)]
    rng = (GATES["observed_order_min"], GATES["observed_order_max"])
    gates = [gate("T_e observed order (finest pair)", o_top[-1], rng, "range"),
             gate("T_source observed order (finest pair)", o_bot[-1], rng, "range"),
             gate("finest T_e error [K]", e_top[-1], GATES["thermal_finest_abs_K"]),
             gate("finest T_source error [K]", e_bot[-1], GATES["thermal_finest_abs_K"]),
             gate("global energy balance (rel, all levels)", max(r["global_energy_rel"] for r in rows),
                  GATES["energy_global"])]
    return _suite("thermal", gates, {"levels": rows, "orders_T_e": o_top, "orders_T_source": o_bot,
                                     "manufactured_solution": "theta = A + B (t_s - zeta) + C cos(m pi x/L) exp(m pi zeta/L)",
                                     "parameters": P})


# --------------------------------------------------------------------------------------------- coupled
def _endpoint_metrics(sol) -> dict:
    return {"load_W_m2": sol.load, "Q_feed_m3_s": float(sol.Q_faces[0]), "max_dTe_K": float(np.max(sol.slab.theta_e)),
            "max_dTsource_K": float(np.max(sol.slab.theta_source)), "iterations": sol.iterations}


FIXED_LOAD_W_M2 = 4.7e5  # forward-solve QoIs at a mesh-independent load (about half the manufactured limit)
SYSTEMATIC_FLOOR_REL = 1e-12


def _systematic(values, scale) -> bool:
    """Differences between successive refinements shrink (or are already at roundoff)."""
    d = np.abs(np.diff(values))
    if d.size < 2:
        return True
    return bool(np.all(d[1:] <= d[:-1]) or d[-1] <= SYSTEMATIC_FLOOR_REL * scale)


def suite_coupled(levels=REFINEMENT, tol_level=(100, 16)) -> dict:
    p = fixture_problem()
    rows, fwd_rows, gates = [], [], []
    for Nx, Nz in levels:
        res = search_dryout(p, SolverSettings(Nx=Nx, Nz=Nz))
        if res.endpoint is None or not res.endpoint.converged:
            gates.append(gate(f"({Nx},{Nz}) search converged", False, True, "equal"))
            continue
        e = res.endpoint
        rows.append({"Nx": Nx, "Nz": Nz, **_endpoint_metrics(e), "status": res.status.value,
                     "relative_interval": res.relative_interval, "monotonic": res.monotonic_end_suction})
        fw = solve_at_load(p, FIXED_LOAD_W_M2, SolverSettings(Nx=Nx, Nz=Nz))
        if fw.converged:
            fwd_rows.append({"Nx": Nx, "Nz": Nz, "max_suction_Pa": float(fw.Pc_end), "Q_feed_m3_s": float(fw.Q_faces[0]),
                             "max_dTe_K": float(np.max(fw.slab.theta_e)),
                             "max_dTsource_K": float(np.max(fw.slab.theta_source))})
        gates += _balance_gates(p, e, f"({Nx},{Nz}) endpoint")
        gates.append(gate(f"({Nx},{Nz}) limit interval (rel)", res.relative_interval, GATES["limit_interval_rel"]))
        gates.append(gate(f"({Nx},{Nz}) end suction monotonic along continuation", res.monotonic_end_suction, True, "equal"))
        # evaporation-moment identity with independent quadrature of tau
        f = p.liquid
        lhs = _psi_quad(p.closure, e.Pc_end) / f.mu
        rhs = analytic.evaporation_moment(e.mesh.x_faces, e.q_e, f.rho, f.h_fg)
        gates.append(gate(f"({Nx},{Nz}) moment identity (rel)", abs(lhs - rhs) / rhs, GATES["moment_identity_rel"]))
    gates.append(gate("all refinement levels solved", len(rows) == len(levels) and len(fwd_rows) == len(levels),
                      True, "equal"))
    for tag, rr, rel_keys in (("endpoint", rows, ("load_W_m2", "Q_feed_m3_s")),
                              ("fixed load", fwd_rows, ("max_suction_Pa", "Q_feed_m3_s"))):
        if len(rr) < 2:
            continue
        a, b = rr[-2], rr[-1]
        for k in rel_keys:
            gates.append(gate(f"{tag}: finest-pair change {k} (rel)", abs(a[k] - b[k]) / abs(b[k]),
                              GATES["mesh_finest_pair_rel"]))
        for k in ("max_dTe_K", "max_dTsource_K"):
            gates.append(gate(f"{tag}: finest-pair change {k} (rel to superheat)", abs(a[k] - b[k]) / b[k],
                              GATES["mesh_finest_pair_rel"]))
            gates.append(gate(f"{tag}: finest-pair change {k} [K]", abs(a[k] - b[k]), GATES["mesh_finest_pair_abs_K"]))
        for k in (*rel_keys, "max_dTe_K", "max_dTsource_K"):
            vals = [r[k] for r in rr]
            gates.append(gate(f"{tag}: systematic convergence of {k}", _systematic(vals, abs(vals[-1])), True, "equal"))
    Nx, Nz = tol_level
    lo = solve_at_endpoint(p, MANUFACTURED_ENDPOINT_PA, SolverSettings(Nx=Nx, Nz=Nz, rtol=1e-8))
    hi = solve_at_endpoint(p, MANUFACTURED_ENDPOINT_PA, SolverSettings(Nx=Nx, Nz=Nz, rtol=1e-10))
    gates.append(gate("tolerance 1e-8 vs 1e-10: dryout load (rel)", abs(lo.load - hi.load) / hi.load, GATES["tolerance_rel"]))
    gates.append(gate("tolerance 1e-8 vs 1e-10: max T_e [K]",
                      float(np.max(np.abs(lo.slab.theta_e - hi.slab.theta_e))), GATES["tolerance_abs_K"]))
    lo_f = solve_at_load(p, FIXED_LOAD_W_M2, SolverSettings(Nx=Nx, Nz=Nz, rtol=1e-8))
    hi_f = solve_at_load(p, FIXED_LOAD_W_M2, SolverSettings(Nx=Nx, Nz=Nz, rtol=1e-10))
    gates.append(gate("tolerance 1e-8 vs 1e-10 (fixed load): max suction (rel)",
                      abs(lo_f.Pc_end - hi_f.Pc_end) / hi_f.Pc_end, GATES["tolerance_rel"]))
    gates.append(gate("tolerance 1e-8 vs 1e-10 (fixed load): max T_e [K]",
                      float(np.max(np.abs(lo_f.slab.theta_e - hi_f.slab.theta_e))), GATES["tolerance_abs_K"]))
    fN = SolverSettings(Nx=levels[-1][0], Nz=levels[-1][1])
    # evaporation feedback: coupled q_e (nonuniform through h_e(Pc)) changes the pressure field
    mid = solve_at_load(p, FIXED_LOAD_W_M2, fN)
    gates.append(gate("forward solve at fixed load converged", mid.converged, True, "equal"))
    feedback = None
    if mid.converged:
        gates += _balance_gates(p, mid, "forward fixed load")
        f = p.liquid
        hyd = DarcyHydraulics(p.closure, f.mu, p.geometry.W, f.rho, f.h_fg)
        uni = hyd.march(mid.mesh.x_faces, mid.q_in, p.Pc_feed, p.closure.Pc_hi)
        feedback = {"q_e_nonuniformity_ptp_over_mean": float(np.ptp(mid.q_e) / np.mean(mid.q_e)),
                    "Pc_end_coupled_Pa": mid.Pc_end, "Pc_end_uniform_evaporation_Pa": float(uni.Pc_faces[-1])}
        gates.append(gate("nonuniform h_e redistributes evaporation (ptp/mean)",
                          feedback["q_e_nonuniformity_ptp_over_mean"], 1e-3, "min"))
        gates.append(gate("evaporation feedback changes end suction (rel)",
                          abs(mid.Pc_end - uni.Pc_faces[-1]) / mid.Pc_end, 1e-4, "min"))
    # nonuniform heating exercises lateral conduction
    pn = fixture_problem()
    pn.heat_shape = lambda xi: 1.0 + 0.5 * np.cos(np.pi * xi)
    pn.heat_shape_label = "1 + 0.5 cos(pi x/L)"
    nu = solve_at_load(pn, FIXED_LOAD_W_M2, fN)
    gates.append(gate("nonuniform heating converged", nu.converged, True, "equal"))
    if nu.converged:
        gates += _balance_gates(pn, nu, "nonuniform heating")
    # table-sampling (knot) sensitivity: same manufactured tables with every other knot dropped
    ps = fixture_problem()
    ds = load_dataset(FIXTURE_TABLES, repo_root=REPO_ROOT)
    ps.closure = ds.to_closure(1e-9, knot_stride=2)
    e2 = solve_at_endpoint(ps, MANUFACTURED_ENDPOINT_PA, fN)
    e1 = solve_at_endpoint(fixture_problem(), MANUFACTURED_ENDPOINT_PA, fN)
    sampling = {"knot_spacing_Pa": [100.0, 200.0], "load_rel_change": abs(e2.load - e1.load) / e1.load,
                "max_dTe_change_K": float(np.max(np.abs(e2.slab.theta_e - e1.slab.theta_e))),
                "note": "reported, not gated: interpolation error of the manufactured tables, distinct from device-mesh "
                        "and cell-FEM convergence"}
    return _suite("coupled", gates, {"refinement_endpoint": rows, "refinement_fixed_load": fwd_rows,
                                     "fixed_load_W_m2": FIXED_LOAD_W_M2,
                                     "tolerance_pair": {"rtol_1e-8": _endpoint_metrics(lo),
                                                        "rtol_1e-10": _endpoint_metrics(hi)},
                                     "evaporation_feedback": feedback, "table_sampling": sampling,
                                     "tables": "MANUFACTURED_FIXTURE"})


# --------------------------------------------------------------------------------------------- statuses
def suite_statuses(workdir: Path, Nx=100, Nz=16) -> dict:
    s = SolverSettings(Nx=Nx, Nz=Nz)
    out, gates = {}, []

    def run(name, problem, expect, **kw):
        res = search_dryout(problem, s, regime_check=lambda e: regime_diagnostics(problem, e), **kw)
        out[name] = {"status": res.status.value, "exit_code": EXIT_CODES[res.status], "reason": res.reason,
                     "upper_bracket_W_m2": res.upper_bracket, "lower_operating_point_W_m2": res.lower_operating_point}
        gates.append(gate(f"{name} -> {expect.value}", res.status.value, expect.value, "equal"))
        return res

    base = run("base", fixture_problem(), SolverStatus.PHYSICAL_CAPILLARY_LIMIT)
    cov = run("coverage_ends_at_1500Pa", fixture_problem(workdir, "coverage", pc_max=1500.0, endpoint_demonstrated=False),
              SolverStatus.TABLE_COVERAGE_LIMIT)
    gates.append(gate("coverage: no upper bracket reported", cov.upper_bracket is None, True, "equal"))
    gates.append(gate("coverage: lower operating point below the base limit",
                      bool(cov.lower_operating_point and cov.lower_operating_point < base.upper_bracket), True, "equal"))
    run("geometry_limit_1800Pa", fixture_problem(workdir, "geometry", geometry_limit=1800.0),
        SolverStatus.MENISCUS_GEOMETRY_LIMIT)
    run("narrow_gradient_range", fixture_problem(workdir, "narrow_gradient", gradient_hi=1e5),
        SolverStatus.OUTSIDE_MODEL_REGIME)
    run("narrow_heat_flux_range", fixture_problem(workdir, "narrow_heat_flux", heat_flux_range=(0.0, 5e5)),
        SolverStatus.OUTSIDE_MODEL_REGIME)
    run("injected_nonfinite", fixture_problem(), SolverStatus.NUMERICAL_FAILURE, _inject_failure=True)
    neg = solve_at_load(fixture_problem(), -1e4, s)
    out["negative_load"] = {"status": neg.status.value, "reason": neg.reason}
    gates.append(gate("negative (condensing) load -> OUTSIDE_MODEL_REGIME", neg.status.value,
                      SolverStatus.OUTSIDE_MODEL_REGIME.value, "equal"))
    codes = {st: EXIT_CODES[st] for st in SolverStatus}
    nonzero = [c for st, c in codes.items() if c != 0]
    gates.append(gate("nonzero exit codes are distinct", len(set(nonzero)) == len(nonzero), True, "equal"))
    return _suite("statuses", gates, {"runs": out, "exit_codes": {k.value: v for k, v in codes.items()}})


# --------------------------------------------------------------------------------------------- reference scales
def suite_reference_scales() -> dict:
    case = load_case(REPO_ROOT / "parameters" / "cases" / "vaartstra2019_fig4.yaml")
    g, f = case.geometry, case.liquid
    sc = reference_scales(82e4, g.L, g.h, g.d, g.l, f.rho, f.mu, f.h_fg, f.cp, f.k, dT=10.0)
    comp = {k: {"computed": sc[k], "audit": v, "rel_diff": (sc[k] - v) / v} for k, v in AUDIT_SCALES_82.items()}
    # informational: cross-checks the audit's scales; no physical gate is implied
    return {"suite": "reference_scales", "passed": True, "gates": [], "details": {
        "load_W_cm2": 82, "note": "82 W/cm^2 is a COMPARISON reference load here, never a fitted input; "
                                  "Ja uses an assumed 10 K superheat; estimates only",
        "scales": sc, "audit_comparison": comp}}


SUITES = {
    "analytic": lambda wd: suite_analytic(),
    "accounting": lambda wd: suite_accounting(),
    "zero_load": lambda wd: suite_zero_load(),
    "hydraulics": lambda wd: suite_hydraulics(),
    "thermal": lambda wd: suite_thermal(),
    "coupled": lambda wd: suite_coupled(),
    "statuses": lambda wd: suite_statuses(wd),
    "reference_scales": lambda wd: suite_reference_scales(),
}
