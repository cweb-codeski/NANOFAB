"""Whole-wick reduced-order model: verification tests (docs/whole_wick_assumptions.md)."""

import math
import warnings

import numpy as np
import pytest

from nanofab import dryout as dryout1d
from nanofab.params import ParameterError
from nanofab.permeability import pillar_wick_permeability
from nanofab.provenance import (
    ANALYTICAL,
    ASSUMED,
    COMSOL,
    COMSOL_EXTRACTED,
    ESTIMATE,
    LITERATURE,
    MODEL_CALCULATED,
    PRIMARY_SOURCE,
    USER,
    USER_PROVIDED,
    VALID_PROVENANCE,
    report_category,
    weakest,
)
from nanofab.whole_wick import (
    ExtrapolationError,
    ExtrapolationLog,
    ExtrapolationWarning,
    FeedSpec,
    Grid,
    LookupTable,
    SolverOptions,
    WickModel,
    find_dryout,
    heat_load_sweep,
    load_case_set,
    margins,
    solve_coupled,
    solve_darcy,
    solve_thermal,
)
from nanofab.capillary import pillar_max_capillary_pressure
from nanofab.whole_wick.config import resolve_tree
from nanofab.whole_wick.electrowetting import (
    EPSILON_0,
    DielectricCoating,
    cos_theta_V,
    electrowetting_number,
    pc_max,
    voltage_for_full_wetting,
)
from nanofab.whole_wick.maps import cell_type_map, heat_flux_map
from nanofab.whole_wick.model import ELECTRICAL_POWER_STATUS
from nanofab.whole_wick.properties import CellType, FluidState, ScalarProperty, evaluate_properties

FLUID = FluidState(T_sat=373.15, rho_l=958.4, mu_l=2.82e-4, k_l=0.679, surface_tension=0.0589, h_lv=2.257e6)
GAMMA = FLUID.surface_tension
H = 25e-6
SIO2 = DielectricCoating("sio2", epsilon_r=3.9, t_dielectric=1e-6)


def make_cell(name="c", K=3e-11, h_e=6e4, pc_max0=2000.0, k_eff=10.0, theta0_deg=15.0, coating=None,
              prov=MODEL_CALCULATED):
    theta0 = math.radians(theta0_deg)
    r_eff = 2 * GAMMA * math.cos(theta0) / pc_max0
    pc_top = 2 * GAMMA / r_eff if coating else pc_max0
    pc = np.linspace(0, pc_top, 5)
    t = {
        "K": LookupTable(f"{name}.K", "m2", {"Pc": pc}, np.full(pc.size, K), prov),
        "h_e": LookupTable(f"{name}.h_e", "W/m2/K", {"Pc": pc}, np.full(pc.size, h_e), prov),
    }
    return CellType(name, 10e-6, 30e-6, H, "silicon", "silicon", theta0, 148.0, r_eff, coating, t,
                    ScalarProperty(k_eff, "W/m/K", prov), prov)


def make_model(nx=40, ny=2, cells=None, type_map=None, V=None, q=50e4, feeds=None, edges=None, policy="warn"):
    g = Grid(nx, ny, 10e-3, 10e-3)
    cells = cells or [make_cell()]
    q_map = np.full(g.shape, q) if np.isscalar(q) else q
    return WickModel(
        name="test", grid=g, cell_types=cells,
        type_map=np.zeros(g.shape, int) if type_map is None else type_map,
        V=np.zeros(g.shape) if V is None else V, q_map=q_map, fluid=FLUID,
        substrate_k=148.0, substrate_t=100e-6,
        feeds=feeds or FeedSpec(edges={"left": 0.0, "right": 0.0}),
        thermal_edges=edges or {}, options=SolverOptions(extrapolation_policy=policy),
    )


@pytest.fixture(scope="module")
def case_set():
    return load_case_set()


# ----------------------------------------------------------------- uniform-property Darcy flow


def test_darcy_linear_pressure_without_evaporation():
    g = Grid(50, 3, 5e-3, 1e-3)
    K = np.full(g.shape, 3e-11)
    ones = np.ones(g.shape)
    res = solve_darcy(g, K, H * ones, FLUID.rho_l * ones, FLUID.mu_l * ones, np.zeros(g.shape),
                      FeedSpec(edges={"left": 0.0, "right": 1000.0}))
    exact = 1000.0 * g.x / g.lx  # Pc rises linearly from 0 to 1000 Pa
    assert np.allclose(res.Pc, exact[None, :], rtol=0, atol=1e-9)
    Q = FLUID.rho_l * 3e-11 * H / FLUID.mu_l * g.ly * 1000.0 / g.lx
    assert res.mass.inlet == pytest.approx(Q, rel=1e-10)
    assert res.mass.outlet == pytest.approx(Q, rel=1e-10)


def test_darcy_uniform_evaporation_matches_1d_profile():
    nx = 200
    g = Grid(nx, 2, 10e-3, 1e-3)
    ones = np.ones(g.shape)
    q = 50e4
    res = solve_darcy(g, 3e-11 * ones, H * ones, FLUID.rho_l * ones, FLUID.mu_l * ones, q / FLUID.h_lv * ones,
                      FeedSpec(edges={"left": 0.0, "right": 0.0}))
    L = g.lx / 2
    xh = np.minimum(g.x, g.lx - g.x)
    exact = -dryout1d.pressure_profile(xh, q, 3e-11, H, L, FLUID.rho_l, FLUID.mu_l, FLUID.h_lv)
    assert np.max(np.abs(res.Pc[0] - exact)) / exact.max() < 1e-4
    assert res.Pc.max() == pytest.approx(dryout1d.wick_pressure_drop(q, 3e-11, H, L, FLUID.rho_l, FLUID.mu_l, FLUID.h_lv), rel=1e-6)


def test_darcy_gravity_hydrostatic():
    """No evaporation and in-plane gravity: closed side, pressure grows by rho g dx (no flow)."""
    g = Grid(20, 1, 1e-2, 1e-3)
    ones = np.ones(g.shape)
    res = solve_darcy(g, 3e-11 * ones, H * ones, FLUID.rho_l * ones, FLUID.mu_l * ones, np.zeros(g.shape),
                      FeedSpec(edges={"left": 0.0}), gravity=(9.81, 0.0))
    assert np.allclose(np.diff(res.p[0]), FLUID.rho_l * 9.81 * g.dx, rtol=1e-10)
    assert abs(res.mass.inlet) < 1e-18


# ----------------------------------------------------------------- analytical permeability scaling


def test_permeability_scales_with_length_squared(case_set):
    K1 = pillar_wick_permeability(10e-6, 30e-6, 25e-6)
    for s in (0.5, 2.0, 3.0):
        assert pillar_wick_permeability(10e-6 * s, 30e-6 * s, 25e-6 * s) == pytest.approx(K1 * s**2, rel=1e-12)
    base = next(c for c in case_set.cell_types if c.name == "si_d10_l30_h25")
    assert base.tables["K"].values[0] == pytest.approx(K1, rel=1e-12)


def test_whole_wick_dryout_scales_as_h_squared():
    q = []
    for h in (10, 20, 40):
        geo = {f"cells.cell_types.si_d10_l30_h25.{k}": {"value": v, "unit": "um", "provenance": USER_PROVIDED,
                                                        "source": None, "source_location": "test"}
               for k, v in (("pillar_diameter", 0.4 * h), ("pillar_pitch", 1.2 * h), ("pillar_height", h))}
        geo["wick.grid.nx.value"] = 20
        geo["wick.grid.ny.value"] = 2
        cs = load_case_set(overrides=geo)
        q.append(find_dryout(cs.build("passive_uniform", "uniform_50")).mean_heat_flux)
    assert q[1] / q[0] == pytest.approx(4.0, rel=1e-3)
    assert q[2] / q[1] == pytest.approx(4.0, rel=1e-3)


# ----------------------------------------------------------------- conservation


def test_mass_conservation_heterogeneous():
    rng = np.random.default_rng(1)
    g = Grid(30, 25, 6e-3, 5e-3)
    K = np.exp(rng.normal(np.log(3e-11), 0.8, g.shape))
    h = rng.uniform(15e-6, 40e-6, g.shape)
    m = rng.uniform(0, 5e-4, g.shape)
    cells = np.zeros(g.shape, bool)
    cells[12, 15] = True
    res = solve_darcy(g, K, h, FLUID.rho_l * np.ones(g.shape), FLUID.mu_l * np.ones(g.shape), m,
                      FeedSpec(edges={"left": 0.0, "top": 50.0}, cells=cells, cell_pc=0.0), average="harmonic")
    assert res.mass.relative_residual < 1e-10
    assert res.mass.max_cell_residual < 1e-12 * res.mass.evaporation
    assert res.mass.evaporation == pytest.approx((m * g.cell_area).sum(), rel=1e-12)


def test_energy_conservation_with_fixed_edge():
    rng = np.random.default_rng(2)
    g = Grid(25, 20, 5e-3, 4e-3)
    G = rng.uniform(5e-3, 2e-2, g.shape)
    he = rng.uniform(0, 8e4, g.shape)
    q = rng.uniform(1e5, 1e6, g.shape)
    th = solve_thermal(g, G, he, q, 373.15, edges={"left": 370.0})
    e = th.energy
    assert e.relative_residual < 1e-10
    assert e.heat_input == pytest.approx(e.evaporative + e.edge_conduction, rel=1e-10)


def test_coupled_conservation_all_cases(case_set):
    for case in case_set.case_names:
        r = solve_coupled(case_set.build(case, "gaussian_hotspot"), 1.0)
        assert r.converged
        assert r.darcy.mass.relative_residual < 1e-10
        assert r.thermal.energy.relative_residual < 1e-10
        # adiabatic edges: all heat (incl. electrical) leaves by evaporation, and evaporation = liquid supply
        assert r.thermal.energy.evaporative == pytest.approx(r.thermal.energy.heat_input, rel=1e-10)
        assert r.darcy.mass.inlet * FLUID.h_lv == pytest.approx(r.thermal.energy.evaporative, rel=2e-3)


# ----------------------------------------------------------------- monotonicity


def test_monotonic_with_heat_load(case_set):
    rows = heat_load_sweep(case_set.build("passive_uniform", "gaussian_hotspot"), [0.2, 0.5, 1.0, 1.5, 2.0, 3.0])
    for key in ("T_max_C", "max_pressure_ratio", "Q_evap_W", "inlet_kg_s"):
        vals = [r[key] for r in rows]
        assert all(b > a for a, b in zip(vals, vals[1:])), key
    margins_ = [r["min_pressure_margin_Pa"] for r in rows]
    assert all(b < a for a, b in zip(margins_, margins_[1:]))


# ----------------------------------------------------------------- dryout detection


def test_dryout_matches_1d_analytical():
    m = make_model(nx=40, ny=2, q=50e4)
    d = find_dryout(m)
    exact = dryout1d.dryout_heat_flux(3e-11, H, 5e-3, 2000.0, FLUID.rho_l, FLUID.mu_l, FLUID.h_lv)
    assert d.mean_heat_flux == pytest.approx(exact, rel=1e-4)
    assert abs(d.location[0] - 5e-3) <= m.grid.dx
    assert abs(d.pressure_margin_at_cell) < 1e-3 * 2000.0
    assert d.criterion == "capillary_pressure"


def test_supply_criterion_equivalent_to_pressure_criterion():
    m = make_model(nx=40, ny=2, q=150e4)  # well above dryout
    r = solve_coupled(m, 1.0)
    mg = margins(r)
    evaporating = r.m_evap > 0
    assert r.dry.any()
    np.testing.assert_array_equal((mg.supply_ratio >= 1)[evaporating], (mg.pressure_ratio >= 1)[evaporating])
    np.testing.assert_array_equal(r.dry, mg.pressure_ratio > 1)


def test_feed_capacity_criterion():
    cap = 2e-5  # kg/s
    m = make_model(nx=20, ny=2, q=50e4, feeds=FeedSpec(edges={"left": 0.0, "right": 0.0}, max_mass_flow=cap))
    d = find_dryout(m)
    assert d.criterion == "feed_capacity"
    assert d.heat_load == pytest.approx(cap * FLUID.h_lv, rel=1e-4)


def test_dry_patch_mode_removes_evaporation_from_dry_cells():
    m = make_model(nx=40, ny=2, q=150e4)
    r = solve_coupled(m, 1.0, dry_model="dry_patch")
    assert r.dry.any()
    assert np.all(r.h_e_used[r.dry] == 0.0)
    assert r.thermal.energy.relative_residual < 1e-10


def test_hotspot_moves_first_dryout(case_set):
    d_uni = find_dryout(case_set.build("passive_uniform", "uniform_50"))
    d_hot = find_dryout(case_set.build("passive_uniform", "gaussian_hotspot"))
    assert d_hot.mean_heat_flux < d_uni.mean_heat_flux
    assert d_hot.peak_heat_flux > d_uni.peak_heat_flux


# ----------------------------------------------------------------- lookup-table bounds


def _table():
    return LookupTable("t", "-", {"Pc": [0.0, 100.0, 200.0], "E": [0.0, 1e6]},
                       np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]), ASSUMED)


def test_lookup_interpolates_inside_bounds_without_warning():
    t = _table()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert t(Pc=100.0, E=0.0) == pytest.approx(3.0)
        assert t(Pc=50.0, E=5e5) == pytest.approx(2.5)
        assert t(Pc=200.0, E=1e6) == pytest.approx(6.0)


def test_lookup_raise_policy():
    with pytest.raises(ExtrapolationError, match="outside table range"):
        _table()("raise", Pc=250.0, E=0.0)
    with pytest.raises(ExtrapolationError):
        _table()("raise", Pc=10.0, E=2e6)


def test_lookup_warn_policy_clips_and_logs():
    log = ExtrapolationLog()
    with pytest.warns(ExtrapolationWarning):
        v = _table()("warn", log, Pc=np.array([150.0, 400.0]), E=0.0)
    assert v[1] == pytest.approx(5.0)  # clipped to the Pc = 200 row, not extrapolated (would be 9)
    assert len(log) == 1 and log.events[0]["requested_max"] == 400.0 and log.events[0]["n_points"] == 1


def test_single_point_axis_rejects_other_values():
    t = LookupTable("z", "-", {"E": [0.0]}, np.array([7.0]), ASSUMED)
    assert t("raise", E=0.0) == pytest.approx(7.0)
    with pytest.raises(ExtrapolationError):
        t("raise", E=1.0)


def test_lookup_rejects_bad_tables():
    with pytest.raises(ValueError):
        LookupTable("bad", "-", {"Pc": [0.0, 0.0]}, np.array([1.0, 2.0]), ASSUMED)
    with pytest.raises(ValueError):
        LookupTable("bad", "-", {"Pc": [0.0, 1.0]}, np.array([1.0, np.nan]), ASSUMED)
    with pytest.raises(ValueError):
        LookupTable("bad", "-", {"Pc": [0.0, 1.0]}, np.array([1.0, 2.0]), "LITERATURE")


# ----------------------------------------------------------------- electrowetting


def test_lippmann_young_formula():
    V = 8.0
    eta = 3.9 * EPSILON_0 * V**2 / (2 * GAMMA * 1e-6)
    assert electrowetting_number(V, SIO2, GAMMA) == pytest.approx(eta, rel=1e-12)
    c, sat = cos_theta_V(math.radians(30), V, SIO2, GAMMA)
    assert float(c) == pytest.approx(math.cos(math.radians(30)) + eta, rel=1e-12)
    assert not sat


def test_lippmann_young_bounds():
    c, sat = cos_theta_V(math.radians(15), 100.0, SIO2, GAMMA)
    assert float(c) == 1.0 and sat  # cos theta_V never exceeds 1
    capped = DielectricCoating("s", 3.9, 1e-6, saturation_angle=math.radians(10))
    c, sat = cos_theta_V(math.radians(15), 100.0, capped, GAMMA)
    assert float(c) == pytest.approx(math.cos(math.radians(10))) and sat
    below = DielectricCoating("s", 3.9, 1e-6, saturation_angle=math.radians(20))
    c, _ = cos_theta_V(math.radians(15), 100.0, below, GAMMA)
    assert float(c) == pytest.approx(math.cos(math.radians(15)))  # theta_0 < theta_sat: no effect
    Vfull = voltage_for_full_wetting(math.radians(15), SIO2, GAMMA)
    assert float(cos_theta_V(math.radians(15), Vfull, SIO2, GAMMA)[0]) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        electrowetting_number(-1.0, SIO2, GAMMA)


def test_zero_voltage_pc_max_matches_eq3(case_set):
    g = case_set.gamma.value
    for ct in case_set.cell_types:
        assert ct.pc_max0(g) == pytest.approx(pillar_max_capillary_pressure(g, ct.d, ct.l, ct.theta_0), rel=1e-12)


def test_zero_voltage_coated_equals_uncoated():
    a = solve_coupled(make_model(cells=[make_cell()]), 1.0)
    b = solve_coupled(make_model(cells=[make_cell(coating=SIO2)]), 1.0)
    np.testing.assert_allclose(a.thermal.T, b.thermal.T, rtol=0, atol=1e-12)
    np.testing.assert_allclose(a.Pc, b.Pc, rtol=1e-12)
    np.testing.assert_array_equal(a.props.Pc_max, b.props.Pc_max)
    assert find_dryout(make_model(cells=[make_cell()])).scale == pytest.approx(
        find_dryout(make_model(cells=[make_cell(coating=SIO2)])).scale, rel=1e-9)


def test_global_ew_dryout_scales_with_cos_theta():
    """Constant h_e and K: 1D dryout is proportional to Pc_max, so the ratio is cos(theta_V)/cos(theta_0)."""
    cell = make_cell(coating=SIO2, theta0_deg=40.0)
    V = 12.0
    d0 = find_dryout(make_model(cells=[cell]))
    dV = find_dryout(make_model(cells=[cell], V=np.full((2, 40), V)))
    c, _ = cos_theta_V(math.radians(40), V, SIO2, GAMMA)
    assert dV.mean_heat_flux / d0.mean_heat_flux == pytest.approx(float(c) / math.cos(math.radians(40)), rel=1e-4)
    assert dV.pc_max_at_cell == pytest.approx(float(pc_max(GAMMA, c, cell.r_eff)), rel=1e-12)
    assert dV.required_dp_at_cell == pytest.approx(dV.pc_max_at_cell, rel=1e-4)


def test_voltage_changes_only_pc_max():
    cell = make_cell(coating=SIO2)
    tmap = np.zeros((2, 4), int)
    Pc = np.full((2, 4), 1500.0)
    a = evaluate_properties([cell], tmap, Pc, np.zeros((2, 4)), GAMMA)
    b = evaluate_properties([cell], tmap, Pc, np.full((2, 4), 10.0), GAMMA)
    for name in ("K", "h_e", "k_eff", "Pc_max0", "theta_0"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
    assert np.all(b.Pc_max > a.Pc_max) and np.all(b.theta_V < a.theta_V)


def test_voltage_on_uncoated_cell_raises():
    with pytest.raises(ValueError, match="dielectric coating"):
        make_model(V=np.full((2, 40), 10.0))


def test_repo_passive_cases_have_zero_voltage(case_set):
    for case in ("passive_uniform", "passive_graded"):
        r = solve_coupled(case_set.build(case, "uniform_50"), 1.0)
        assert not r.model.V.any() and r.model.actuation == "none"
        np.testing.assert_array_equal(r.props.Pc_max, r.props.Pc_max0)
        assert r.thermal.energy.electrical_heat == 0.0


# ----------------------------------------------------------------- selective actuation masking


def test_selective_voltage_only_in_vulnerable_cells(case_set):
    spec = case_set.cases.tree["voltage_maps"]["selective"]
    sel = case_set.build("selective_ew", "uniform_50")
    ref = find_dryout(case_set.build("passive_uniform", "uniform_50")).result
    mask = margins(ref).pressure_ratio >= spec["threshold"].value
    assert mask.any() and not mask.all()
    np.testing.assert_array_equal(sel.V > 0, mask)
    assert np.all(sel.V[mask] == spec["V"].value)
    assert not (sel.V[:, 0].any() or sel.V[:, -1].any())  # far from feeds
    assert sel.actuation == "selective"


def test_selective_voltage_leaves_unmasked_cells_unchanged(case_set):
    sel = case_set.build("selective_ew", "uniform_50")
    pas = case_set.build("passive_uniform", "uniform_50")
    Pc = np.full(sel.grid.shape, 500.0)
    g = case_set.gamma.value
    ps = evaluate_properties(sel.cell_types, sel.type_map, Pc, sel.V, g)
    pp = evaluate_properties(pas.cell_types, pas.type_map, Pc, pas.V, g)
    off = sel.V == 0
    for name in ("K", "h_e", "k_eff", "Pc_max", "theta_V"):
        np.testing.assert_array_equal(getattr(ps, name)[off], getattr(pp, name)[off])
    assert np.all(ps.Pc_max[~off] > pp.Pc_max[~off])
    np.testing.assert_array_equal(ps.K, pp.K)
    np.testing.assert_array_equal(ps.h_e, pp.h_e)


def test_comparison_reports_ew_outputs_and_unresolved_power(case_set):
    from nanofab.whole_wick.cases import run_case

    ref = run_case(case_set, "passive_uniform", "uniform_50")
    rows = {c: run_case(case_set, c, "uniform_50").row(ref) for c in ("global_ew", "selective_ew")}
    for c, r in rows.items():
        for key in ("V_applied_V", "actuation", "dtheta_deg", "dPc_max_Pa", "dPc_max_pct", "dryout_mean_W_cm2",
                    "first_dry_x_mm", "first_dry_y_mm", "required_dp_at_dryout_Pa", "Pc_maxV_at_dryout_Pa"):
            assert key in r, (c, key)
        assert r["electrical_power"] == ELECTRICAL_POWER_STATUS
        assert not any("P_el" in k or "per_W" in k for k in r)
        assert r["dtheta_deg"] < 0 and r["dPc_max_Pa"] > 0 and r["dQ_dryout_W"] > 0
    assert rows["global_ew"]["actuation"] == "global" and rows["selective_ew"]["actuation"] == "selective"
    assert rows["selective_ew"]["actuated_area_fraction"] < rows["global_ew"]["actuated_area_fraction"] == 1.0


# ----------------------------------------------------------------- maps, materials, provenance


def test_heat_maps():
    g = Grid(41, 41, 10e-3, 10e-3)
    leaf = lambda v, u: resolve_tree({"x": {"value": v, "unit": u, "provenance": USER_PROVIDED,  # noqa: E731
                                            "source": None, "source_location": None}}).tree["x"]
    spec = {"type": "gaussian", "background": leaf(30, "W/cm2"),
            "hotspots": [{"x": leaf(5, "mm"), "y": leaf(5, "mm"), "sigma": leaf(1, "mm"), "peak_increment": leaf(70, "W/cm2")}]}
    q = heat_flux_map(spec, g)
    assert q.max() == pytest.approx(100e4, rel=1e-9)
    assert q[0, 0] == pytest.approx(30e4, rel=1e-3)


def test_material_map_and_harmonic_faces(case_set):
    g = case_set.grid()
    names = case_set.cell_type_names
    raw = {"type": "rectangles", "background": "si_d10_l30_h25", "regions": [
        {"cell_type": "cu_d10_l30_h25", **{k: {"value": v, "unit": "mm", "provenance": USER_PROVIDED, "source": None,
                                               "source_location": None}
                                            for k, v in (("x_min", 0), ("x_max", 5), ("y_min", 0), ("y_max", 10))}}]}
    tmap = cell_type_map(resolve_tree(raw).tree, g, names, ["left", "right"])
    cu = names.index("cu_d10_l30_h25")
    assert (tmap[:, : g.nx // 2] == cu).all() and not (tmap[:, g.nx // 2:] == cu).any()
    props = evaluate_properties(case_set.cell_types, tmap, np.zeros(g.shape), np.zeros(g.shape), case_set.gamma.value)
    assert props.k_eff[0, 0] > props.k_eff[0, -1]


def test_bare_number_rejected():
    with pytest.raises(ParameterError, match="no provenance"):
        resolve_tree({"grid": {"nx": 40}})
    with pytest.raises(ParameterError, match="missing provenance"):
        resolve_tree({"q": {"value": 1, "unit": "W/cm2"}})


def test_defaults_inherited_and_ref_resolved():
    cfg = resolve_tree({"g": {"defaults": {"provenance": USER_PROVIDED, "source": None, "source_location": "x"},
                              "nx": {"value": 3, "unit": "-"}}, "d": {"ref": "pillar_diameter"}})
    assert cfg.registry["g.nx"].provenance == USER_PROVIDED
    assert cfg.tree["d"].provenance == PRIMARY_SOURCE and cfg.tree["d"].value == pytest.approx(10e-6)


def test_report_categories_and_weakest():
    assert {report_category(p) for p in VALID_PROVENANCE} == {LITERATURE, COMSOL, ANALYTICAL, ESTIMATE, USER}
    assert report_category(COMSOL_EXTRACTED) == COMSOL
    assert weakest(PRIMARY_SOURCE, ASSUMED, MODEL_CALCULATED) == ASSUMED
    assert weakest(MODEL_CALCULATED, USER_PROVIDED, PRIMARY_SOURCE) == MODEL_CALCULATED


def test_every_input_and_table_has_provenance(case_set):
    for reg in (case_set.wick.registry, case_set.cells.registry, case_set.cases.registry):
        for path, p in reg.items():
            assert p.provenance in VALID_PROVENANCE, path
    assert case_set.gamma.provenance in VALID_PROVENANCE
    for ct in case_set.cell_types:
        for t in ct.tables.values():
            assert t.provenance in VALID_PROVENANCE
        assert ct.tables["h_e"].provenance == ASSUMED  # placeholder h_e (W18)
        if ct.coating is not None:
            assert ct.pc_max_provenance == ASSUMED  # epsilon_r and saturation angle are ASSUMED
    reg = case_set.cells.registry
    assert reg["dielectric_coatings.sio2_1um.epsilon_r"].provenance == ASSUMED
