import math

import pytest

from nanofab.params import ParameterError, load_parameters
from nanofab.provenance import MODEL_CALCULATED, PRIMARY_SOURCE

FIELDS = "provenance: ASSUMED, source: null, source_location: null"


def write(tmp_path, text):
    (tmp_path / "p.yaml").write_text(text, encoding="utf-8")
    return tmp_path


def test_repo_parameters_load_and_are_si():
    ps = load_parameters()
    assert ps["pillar_diameter"] == pytest.approx(10e-6)
    assert ps["vapor_temperature"] == pytest.approx(373.15)
    assert ps["target_heat_flux"] == pytest.approx(100e4)
    assert ps["receding_contact_angle"] == pytest.approx(math.radians(15))
    assert ps.meta("pillar_diameter").provenance == PRIMARY_SOURCE
    assert ps.meta("pillar_diameter").source == "vaartstra2019"


def test_every_repo_parameter_has_provenance_and_location_field():
    ps = load_parameters()
    for p in ps.parameters.values():
        assert p.provenance
        if p.source is not None:
            assert p.source.strip()


def test_missing_provenance_rejected(tmp_path):
    write(tmp_path, "x: {value: 1, unit: m, source: null, source_location: null}\n")
    with pytest.raises(ParameterError, match="provenance"):
        load_parameters(tmp_path)


def test_invalid_provenance_rejected(tmp_path):
    write(tmp_path, "x: {value: 1, unit: m, provenance: GUESS, source: null, source_location: null}\n")
    with pytest.raises(ParameterError, match="GUESS"):
        load_parameters(tmp_path)


def test_old_literature_label_rejected(tmp_path):
    write(tmp_path, "x: {value: 1, unit: m, provenance: LITERATURE, source: vaartstra2019, source_location: null}\n")
    with pytest.raises(ParameterError, match="LITERATURE"):
        load_parameters(tmp_path)


def test_primary_source_without_source_rejected(tmp_path):
    write(tmp_path, "x: {value: 1, unit: m, provenance: PRIMARY_SOURCE, source: null, source_location: null}\n")
    with pytest.raises(ParameterError, match="requires a source"):
        load_parameters(tmp_path)


def test_unregistered_source_rejected(tmp_path):
    write(tmp_path, "x: {value: 1, unit: m, provenance: PRIMARY_SOURCE, source: nobody2099, source_location: p1}\n")
    with pytest.raises(ParameterError, match="nobody2099"):
        load_parameters(tmp_path)


def test_unknown_unit_rejected(tmp_path):
    write(tmp_path, f"x: {{value: 1, unit: furlong, {FIELDS}}}\n")
    with pytest.raises(ParameterError, match="furlong"):
        load_parameters(tmp_path)


def test_sweep_converted(tmp_path):
    write(tmp_path, f"x: {{value: 10, unit: um, {FIELDS}, sweep: [5, 50]}}\n")
    ps = load_parameters(tmp_path)
    assert ps.meta("x").sweep == pytest.approx((5e-6, 50e-6))


def test_override_marks_model_calculated():
    ps = load_parameters().with_overrides(pillar_diameter=20e-6)
    assert ps["pillar_diameter"] == 20e-6
    assert ps.meta("pillar_diameter").provenance == MODEL_CALCULATED
    assert ps.meta("pillar_diameter").source is None
