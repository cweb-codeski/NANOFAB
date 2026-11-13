"""Reference vs proposed case separation, array extent, interface-bracket permeability, Stage 1 cross-check."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest
import yaml

from nanofab.geometry_cases import (
    CASES_FILE,
    UnresolvedInputError,
    array_extent,
    load_geometry_cases,
    start_reference_interpretations,
)
from nanofab.params import ParameterError
from nanofab.permeability import (
    brinkman_height_factor,
    brinkman_no_slip_top_factor,
    pillar_wick_permeability,
    sangani_acrivos_permeability,
)

REPO = Path(__file__).resolve().parents[1]
REF, PROP = "reference_d10_l30_h25", "proposed_d20_l50_h40"


@pytest.fixture(scope="module")
def cases():
    return load_geometry_cases()


def test_two_cases_with_distinct_roles_and_geometry(cases):
    assert cases[REF].role == "REFERENCE" and cases[PROP].role == "PROPOSED_DESIGN"
    assert (cases[REF].d, cases[REF].l, cases[REF].hp) == pytest.approx((10e-6, 30e-6, 25e-6))
    assert (cases[PROP].d, cases[PROP].l, cases[PROP].hp) == pytest.approx((20e-6, 50e-6, 40e-6))


def test_proposed_geometry_is_user_provided_not_literature(cases):
    for key in ("pillar_diameter", "pillar_pitch", "pillar_height"):
        assert cases[PROP].meta(key).provenance == "USER_PROVIDED"
    assert cases[PROP].meta("receding_contact_angle").provenance == "ASSUMED"
    assert cases[REF].meta("pillar_pitch").provenance == "PRIMARY_SOURCE"


def test_no_comsol_derived_quantity_is_stored_in_either_case(cases):
    forbidden = {"permeability", "K", "h_e", "formation_factor", "accommodation_coefficient"}
    for c in cases.values():
        assert not forbidden & set(c.parameters)


@pytest.mark.parametrize("key", ["feed_length", "substrate_thickness", "heater_footprint",
                                 "reservoir_geometry", "vapor_environment", "cad_scale_and_units"])
def test_unresolved_proposed_device_inputs_refuse_a_value(cases, key):
    with pytest.raises(UnresolvedInputError):
        cases[PROP].si(key)


def test_coupon_array_and_feed_length_are_separate_keys(cases):
    p = cases[PROP]
    assert {"coupon_x_min", "coupon_x_max", "array_columns", "array_rows"} <= set(p.parameters)
    assert "feed_length" in p.unresolved


def test_array_extent_matches_stated_edge_span(cases):
    p = cases[PROP]
    ext = array_extent(int(p.si("array_columns")), int(p.si("array_rows")), p.l, p.d)
    assert ext.edge_span_x == pytest.approx(10.22e-3)
    assert ext.edge_span_y == pytest.approx(7.72e-3)
    assert ext.center_span_x == pytest.approx(204 * 50e-6)


def test_start_reference_readings_do_not_center_the_array(cases):
    rows = start_reference_interpretations(cases[PROP])
    assert {r["reading"] for r in rows} == {"first_pillar_center", "array_edge_corner"}
    for r in rows:
        assert abs(r["offset_from_coupon_center_x"]) > 4e-3
        assert r["x_max"] - r["x_min"] == pytest.approx(10.22e-3)


def test_unresolved_entry_with_a_value_is_rejected(tmp_path):
    raw = yaml.safe_load(CASES_FILE.read_text(encoding="utf-8"))
    raw["cases"][PROP]["device"]["feed_length"]["value"] = 5
    bad = tmp_path / "cases.yaml"
    bad.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ParameterError):
        load_geometry_cases(bad)


def test_no_slip_top_factor_limits():
    K_s = 1e-11
    assert brinkman_no_slip_top_factor(K_s, 1e-3) == pytest.approx(1.0, abs=1e-2)
    h = 1e-9
    assert K_s * brinkman_no_slip_top_factor(K_s, h) == pytest.approx(h**2 / 12, rel=1e-6)
    for h in (5e-6, 25e-6, 40e-6):
        assert brinkman_no_slip_top_factor(K_s, h) < brinkman_height_factor(K_s, h)


def test_no_slip_top_bracket_is_lower_for_both_cases(cases):
    for c in cases.values():
        sf = pillar_wick_permeability(c.d, c.l, c.hp)
        ns = pillar_wick_permeability(c.d, c.l, c.hp, top="no_slip")
        assert 0.5 < ns / sf < 0.8
        assert sf < sangani_acrivos_permeability(c.d, c.l)
    with pytest.raises(ValueError):
        pillar_wick_permeability(10e-6, 30e-6, 25e-6, top="partial")


def _load_script(name: str):
    sys.path.insert(0, str(REPO / "analysis"))
    spec = importlib.util.spec_from_file_location(name, REPO / "analysis" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_stage1_crosscheck_reproduces_reported_relative_changes():
    mod = _load_script("stage1_crosscheck")
    chase, stu, stu_Q = mod.load()
    changes = dict(mod.relative_changes(stu))
    assert changes[0.25] == pytest.approx(5.05e-4, abs=1e-6)
    assert changes[1.0] == pytest.approx((2.9608 - 2.9295) / 2.9295)
    for h in (2.0, 1.0, 0.5):
        assert abs(stu[h] / chase[h] - 1) < 5e-4
    assert (stu / stu_Q).max() / (stu / stu_Q).min() - 1 < 1e-4


def test_observed_order_recovers_a_known_power_law():
    mod = _load_script("stage1_crosscheck")
    hs = (1.0, 0.7, 0.5)
    ks = tuple(3.0 - 0.2 * h**2 for h in hs)
    p, kinf = mod.observed_order(hs, ks)
    assert p == pytest.approx(2.0, rel=1e-6) and kinf == pytest.approx(3.0, rel=1e-9)


def test_geometry_columns_uses_named_case_when_columns_are_absent():
    mod = _load_script("_common")
    from nanofab.comsol_io import ComsolTable
    from nanofab.params import load_parameters
    import pandas as pd

    table = ComsolTable(data=pd.DataFrame({"Q_in": [1.0]}), units={"Q_in": "W"}, metadata={})
    geo = mod.geometry_columns(table, load_parameters(), PROP)
    assert geo["l"] == pytest.approx(50e-6) and geo["h"] == pytest.approx(40e-6)
    assert math.isclose(mod.geometry_columns(table, load_parameters())["l"], 30e-6)
