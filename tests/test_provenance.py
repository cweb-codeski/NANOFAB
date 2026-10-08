"""Provenance rules from docs/SOURCE_AUTHORITY.md."""

import csv

import pytest
import yaml

from nanofab.params import DEFAULT_PARAMETER_DIR, load_parameters
from nanofab.provenance import (
    DESIGN_TARGET,
    LITERATURE_TARGET_PROVENANCE,
    MODEL_CALCULATED,
    PRIMARY_SOURCE,
    REPO_ROOT,
    VALID_PROVENANCE,
    ProvenanceError,
    check_reference,
    load_sources,
    validate_provenance,
)

TARGETS_CSV = REPO_ROOT / "data" / "literature" / "vaartstra2019_targets.csv"
MODEL_PROVENANCE = REPO_ROOT / "config" / "model_provenance.yaml"
DESIGN_FILE = "design_targets.yaml"


@pytest.fixture(scope="module")
def sources():
    return load_sources()


@pytest.fixture(scope="module")
def params():
    return load_parameters()


def target_rows():
    with open(TARGETS_CSV, encoding="utf-8") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def model_entries():
    with open(MODEL_PROVENANCE, encoding="utf-8") as fh:
        content = yaml.safe_load(fh)
    return [(f"{group}.{name}", entry) for group, items in content.items() for name, entry in items.items()]


def test_validate_provenance():
    for label in VALID_PROVENANCE:
        assert validate_provenance(label) == label
    for bad in ("LITERATURE", "CALCULATED", "TARGET", "primary_source", ""):
        with pytest.raises(ProvenanceError):
            validate_provenance(bad)


def test_all_parameters_have_valid_provenance(params):
    for p in params.parameters.values():
        assert p.provenance in VALID_PROVENANCE, p.name


def test_primary_source_parameters_have_source_and_location(params):
    for p in params.by_provenance(PRIMARY_SOURCE):
        assert p.source, p.name
        assert p.source_location, p.name


def test_all_parameter_source_ids_registered(params, sources):
    for p in params.parameters.values():
        if p.source is not None:
            assert p.source in sources, f"{p.name}: {p.source}"


def test_no_model_calculated_parameter_inputs(params):
    assert params.by_provenance(MODEL_CALCULATED) == []


def test_design_targets_separated_from_replication_parameters(params):
    for p in params.parameters.values():
        if p.provenance == DESIGN_TARGET:
            assert p.file == DESIGN_FILE, p.name
        if p.file == DESIGN_FILE:
            assert p.provenance == DESIGN_TARGET, p.name
    with open(DEFAULT_PARAMETER_DIR / DESIGN_FILE, encoding="utf-8") as fh:
        assert yaml.safe_load(fh), "design_targets.yaml should not be empty"


def test_target_heat_flux_is_design_target_without_source(params):
    meta = params.meta("target_heat_flux")
    assert meta.provenance == DESIGN_TARGET
    assert meta.source is None


def test_literature_targets_are_never_model_calculated(sources):
    rows = target_rows()
    assert rows
    for row in rows:
        where = f"{row['case']}/{row['quantity']}"
        assert row["provenance"] in LITERATURE_TARGET_PROVENANCE, where
        assert row["source"] in sources, where
        assert row["source_location"], where
        if row["provenance"] == PRIMARY_SOURCE:
            assert row["value"], f"{where}: PRIMARY_SOURCE target needs a reported value"


def test_reference_dryout_matches_targets_file(params):
    dryout = [r for r in target_rows() if r["case"] == "baseline_device" and r["quantity"] == "q_dryout"]
    assert len(dryout) == 1
    assert float(dryout[0]["value"]) == pytest.approx(params["reference_dryout_heat_flux"])


def test_model_provenance_entries_valid_and_registered(sources):
    for name, entry in model_entries():
        check_reference(entry["provenance"], entry.get("source"), sources, name)
        if entry["provenance"] == PRIMARY_SOURCE:
            assert entry.get("source_location"), name


def test_derived_results_are_not_labeled_as_literature():
    for name, entry in model_entries():
        if name.startswith("derived_results."):
            assert entry["provenance"] != PRIMARY_SOURCE, name


def test_every_source_has_valid_type_and_required_fields(sources):
    for sid, entry in sources.items():
        assert entry["type"] in VALID_PROVENANCE, sid
        for key in ("title", "authors", "year", "notes"):
            assert key in entry, f"{sid}: missing {key}"


def test_primary_model_reference_is_vaartstra2019(sources):
    roles = {sid for sid, e in sources.items() if e.get("role") == "PRIMARY_MODEL_REFERENCE"}
    assert roles == {"vaartstra2019"}
