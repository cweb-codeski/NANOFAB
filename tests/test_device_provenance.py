"""Provenance and data-gate tests: literature cases never fall back to fixtures or legacy models."""

import json
import re
import sys

import pytest
import yaml

from nanofab.device.case import CaseError, load_case, preflight
from nanofab.device.model import CaseRole, SolverStatus
from nanofab.params import load_parameters
from nanofab.provenance import REPO_ROOT, load_sources

LIT_CASE = REPO_ROOT / "parameters" / "cases" / "vaartstra2019_fig4.yaml"
FIX_CASE = REPO_ROOT / "tests" / "fixtures" / "device" / "manufactured_case.yaml"
sys.path.insert(0, str(REPO_ROOT / "analysis"))


@pytest.fixture(scope="module")
def lit():
    return load_case(LIT_CASE)


def test_literature_preflight_reports_missing_inputs(lit):
    pf = preflight(lit)
    assert pf["status"] == SolverStatus.MISSING_VALIDATED_INPUTS.value
    names = " ".join(m["input"] for m in pf["missing_inputs"])
    for t in ("hydraulic", "meniscus", "thermal"):
        assert f"cell table: {t}" in names
    assert any("meland_ytrehus2003" in g["gap"] for g in pf["blocking_source_gaps"])
    assert pf["_dataset"] is None
    assert "never substituted" in pf["fallbacks_used"]


def test_literature_search_cli_exits_3_and_solves_nothing(tmp_path):
    import run_passive_device as cli
    code = cli.main(["--case", str(LIT_CASE), "--search-dryout", "--output", str(tmp_path / "runs"),
                     "--plot-root", str(tmp_path / "plots")])
    assert code == 3
    run = next((tmp_path / "runs").iterdir())
    files = {p.name for p in run.iterdir()}
    assert {"missing_inputs.json", "run_manifest.json"} <= files
    assert not any(n.endswith("_profiles.csv") for n in files)
    man = json.loads((run / "run_manifest.json").read_text())
    assert man["status"] == "MISSING_VALIDATED_INPUTS" and man["exit_code"] == 3
    assert man["git"]["commit"] is not None and "environment" in man


def test_literature_case_rejects_fixture_dataset(tmp_path):
    raw = yaml.safe_load(LIT_CASE.read_text())
    raw["cell_tables"]["dataset"] = "tests/fixtures/device/manufactured_tables"
    p = tmp_path / "lit_with_fixture.yaml"
    p.write_text(yaml.safe_dump(raw))
    pf = preflight(load_case(p))
    assert pf["status"] == SolverStatus.MISSING_VALIDATED_INPUTS.value
    reason = " ".join(m["reason"] for m in pf["missing_inputs"])
    assert "never admitted to a literature case" in reason


def test_fixture_run_is_labelled_verification_only(tmp_path):
    import run_passive_device as cli
    code = cli.main(["--case", str(FIX_CASE), "--load-W-cm2", "40", "--Nx", "40", "--Nz", "8",
                     "--output", str(tmp_path / "runs"), "--plot-root", str(tmp_path / "plots")])
    assert code == 0
    run = next((tmp_path / "runs").iterdir())
    man = json.loads((run / "run_manifest.json").read_text())
    assert man["verification_only"] is True
    assert (run / "forward_profiles.csv").read_text().startswith("# VERIFICATION ONLY")
    assert man["parameters"] and all("provenance" in p for p in man["parameters"])


def test_design_targets_are_rejected_as_case_inputs(tmp_path):
    raw = yaml.safe_load(LIT_CASE.read_text())
    raw["parameters"]["L"] = "target_heat_flux"
    p = tmp_path / "bad.yaml"
    p.write_text(yaml.safe_dump(raw))
    with pytest.raises(CaseError, match="design target"):
        load_case(p)


def test_comparison_targets_must_be_comparison_only(tmp_path):
    raw = yaml.safe_load(LIT_CASE.read_text())
    raw["comparison_targets"][0]["use"] = "INPUT"
    p = tmp_path / "bad.yaml"
    p.write_text(yaml.safe_dump(raw))
    with pytest.raises(CaseError, match="COMPARISON_ONLY"):
        load_case(p)


def test_reference_dryout_flux_is_not_an_input(lit):
    assert "reference_dryout_heat_flux" not in lit.raw["parameters"].values()
    assert lit.role == CaseRole.LITERATURE_REPLICATION


def test_silicon_conductivity_is_provisional():
    m = load_parameters().meta("silicon_thermal_conductivity")
    assert m.provenance == "ASSUMED"
    assert "not a V2019 value" in m.notes.lower() or "not a v2019 value" in m.notes.lower()


def test_sources_registered_with_status():
    s = load_sources()
    assert "pujahari2022" in s and "10.1021/acs.langmuir.2c01048" in str(s["pujahari2022"])
    assert s["meland_ytrehus2003"].get("content_verified") in (False, None)


def test_assumption_register_complete():
    reg = yaml.safe_load((REPO_ROOT / "config" / "device_assumptions.yaml").read_text(encoding="utf-8"))
    allowed = {"SOURCE_STATED", "DIAGNOSED", "GATED", "OPEN", "PROVISIONAL_INPUT"}
    for key, e in reg.items():
        for f in ("assumption", "source", "rationale", "diagnostic", "regime", "failure_consequence",
                  "qualification_status"):
            assert e.get(f), (key, f)
        assert e["qualification_status"] in allowed, key


def test_audit_corrections_applied():
    import test_vaartstra2019 as legacy  # noqa: F401
    assert not hasattr(legacy, "test_flat_meniscus_dryout_bounds_paper_value")
    assert hasattr(legacy, "test_legacy_flat_meniscus_sanity_check")
    for rel in ("docs/assumptions.md", "src/nanofab/evaporation.py", "analysis/baseline_estimates.py",
                "docs/whole_wick_assumptions.md"):
        if not (REPO_ROOT / rel).exists():
            continue
        text = (REPO_ROOT / rel).read_text(encoding="utf-8")
        assert not re.search(r"known upper bound", text, re.I), rel
        assert not re.search(r"paper-exact", text, re.I), rel
    prov = yaml.safe_load((REPO_ROOT / "config" / "model_provenance.yaml").read_text(encoding="utf-8"))
    assert "kinetic_resistance_share_upper_bound" not in json.dumps(prov)
