"""Strict cell-table schema, interpolation, coverage and import tests (MANUFACTURED_FIXTURE only)."""

import copy
import importlib.util
import shutil

import numpy as np
import pytest
import yaml

from nanofab.device.importer import ImportError_, import_from_manifest
from nanofab.device.tables import load_dataset, match_case, sha256_file, write_dataset
from nanofab.device.transport import PiecewiseLinear, TableRangeError
from nanofab.device.verification import FIXTURE_DIR, FIXTURE_TABLES, base_case
from nanofab.provenance import REPO_ROOT


def _fixture_module():
    spec = importlib.util.spec_from_file_location("fixture_tables", FIXTURE_DIR / "fixture_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FT = _fixture_module()


@pytest.fixture(scope="module")
def built():
    return FT.build()


def _write(tmp_path, name, manifest, frames):
    write_dataset(tmp_path / name, copy.deepcopy(manifest), {k: v.copy() for k, v in frames.items()})
    return load_dataset(tmp_path / name, repo_root=REPO_ROOT)


# ------------------------------------------------------------------ interpolation
def test_piecewise_linear_exact_at_knots_and_linear_between():
    f = PiecewiseLinear([0.0, 1.0, 3.0], [2.0, 4.0, 1.0], name="t")
    assert np.allclose(f([0.0, 1.0, 3.0]), [2.0, 4.0, 1.0])
    assert f(0.5)[0] == pytest.approx(3.0)
    assert f(2.0)[0] == pytest.approx(2.5)


def test_no_extrapolation_and_narrow_snap_only():
    f = PiecewiseLinear([0.0, 10.0], [1.0, 2.0], name="t", snap=1e-8)
    assert f(10.0 + 5e-9)[0] == pytest.approx(2.0)
    with pytest.raises(TableRangeError):
        f(10.0 + 1e-6)
    with pytest.raises(TableRangeError):
        f(-1e-6)


@pytest.mark.parametrize("x,y", [([0, 2, 1], [1, 1, 1]), ([0, 1, 1], [1, 1, 1]), ([0, 1], [1, -1]),
                                 ([0, 1], [1, np.nan]), ([0, np.inf], [1, 1])])
def test_piecewise_linear_rejects_bad_tables(x, y):
    with pytest.raises(ValueError):
        PiecewiseLinear(x, y, name="bad")


def test_integral_and_inverse_round_trip():
    f = PiecewiseLinear([0.0, 1.0, 3.0], [2.0, 4.0, 1.0], name="t")
    P = np.linspace(0, 3, 31)
    S = f.integral(P)
    assert np.allclose(f.integral_inverse(S), P, rtol=0, atol=1e-12)
    assert S[-1] == pytest.approx(0.5 * (2 + 4) * 1 + 0.5 * (4 + 1) * 2)


# ------------------------------------------------------------------ fixture dataset
def test_fixture_is_schema_valid_but_verification_only():
    ds = load_dataset(FIXTURE_TABLES, repo_root=REPO_ROOT)
    r = ds.report
    assert r.schema_valid, r.errors
    assert r.verification_only
    assert not r.numerically_qualified
    match_case(ds, base_case().case_requirements())
    assert r.case_matched, r.case_mismatches


def test_validated_flag_in_manifest_is_ignored(tmp_path, built):
    m, fr = built
    m = copy.deepcopy(m)
    m["validated"] = True
    ds = _write(tmp_path, "flag", m, fr)
    assert not ds.report.numerically_qualified


def test_range_is_intersection_of_tables(tmp_path, built):
    m, fr = copy.deepcopy(built[0]), {k: v.copy() for k, v in built[1].items()}
    fr["thermal"] = fr["thermal"][fr["thermal"].Pc_Pa <= 2000.0]
    m["pressure_interval_Pa"] = [0.0, 2000.0]
    ds = _write(tmp_path, "inter", m, fr)
    assert ds.report.schema_valid, ds.report.errors
    clo = ds.to_closure()
    assert clo.Pc_hi == pytest.approx(2000.0)
    with pytest.raises(TableRangeError):
        clo.tau(2100.0)


def test_coverage_end_is_not_a_demonstrated_endpoint(tmp_path):
    m, fr = FT.build(pc_max=1500.0, endpoint_demonstrated=False)
    ds = _write(tmp_path, "cov", m, fr)
    assert ds.report.schema_valid
    assert any("physical endpoint not demonstrated" in q for q in ds.report.qualification_failures)
    assert ds.to_closure().limits.Pc_physical_limit is None


def _mutated(tmp_path, built, name, fn):
    m, fr = copy.deepcopy(built[0]), {k: v.copy() for k, v in built[1].items()}
    fn(m, fr)
    return _write(tmp_path, name, m, fr)


@pytest.mark.parametrize("name,fn,needle", [
    ("dup", lambda m, fr: fr.__setitem__("hydraulic", fr["hydraulic"].assign(
        Pc_Pa=np.r_[fr["hydraulic"].Pc_Pa.values[:-1], fr["hydraulic"].Pc_Pa.values[-2]])), "strictly increasing"),
    ("negpc", lambda m, fr: fr.__setitem__("thermal", fr["thermal"].assign(Pc_Pa=-fr["thermal"].Pc_Pa[::-1].values)),
     "negative Pc"),
    ("area", lambda m, fr: fr.__setitem__("thermal", fr["thermal"].assign(A_footprint_m2=fr["thermal"].A_footprint_m2 * 4)),
     "area basis"),
    ("knorm", lambda m, fr: fr.__setitem__("hydraulic", fr["hydraulic"].assign(K_bulk_m2=fr["hydraulic"].K_bulk_m2 * 0.75)),
     "normalization mismatch"),
    ("conv", lambda m, fr: m["conventions"].__setitem__("pressure", "p_liquid - p_v"), "conventions"),
    ("plane", lambda m, fr: m["conventions"].__setitem__("thermal_reference_plane", "substrate_bottom"), "conventions"),
    ("nan", lambda m, fr: fr.__setitem__("thermal", fr["thermal"].assign(
        h_e_W_m2K=np.r_[np.nan, fr["thermal"].h_e_W_m2K.values[1:]])), "nonfinite"),
    ("zero_load", lambda m, fr: fr.__setitem__("thermal", fr["thermal"].assign(Q_in_W=0.0, Q_evap_W=0.0)), "zero or negative"),
    ("missing_key", lambda m, fr: m.pop("gradient_study"), "gradient_study"),
    ("role", lambda m, fr: m.__setitem__("data_role", "SYNTHETIC"), "data_role"),
])
def test_malformed_tables_rejected(tmp_path, built, name, fn, needle):
    ds = _mutated(tmp_path, built, name, fn)
    assert not ds.report.schema_valid
    assert any(needle.lower() in e.lower() for e in ds.report.errors), ds.report.errors


def test_unit_and_hash_tampering_rejected(tmp_path, built):
    ds = _write(tmp_path, "ok", *built)
    assert ds.report.schema_valid
    mp = tmp_path / "ok" / "manifest.yaml"
    m = yaml.safe_load(mp.read_text())
    m["tables"]["hydraulic"]["columns"]["K_bulk_m2"] = "um2"
    mp.write_text(yaml.safe_dump(m))
    assert any("unit" in e for e in load_dataset(tmp_path / "ok").report.errors)
    ds2 = _write(tmp_path, "ok2", *built)
    with open(tmp_path / "ok2" / "hydraulic.csv", "a") as fh:
        fh.write("\n")
    assert any("sha256" in e for e in load_dataset(tmp_path / "ok2").report.errors)
    assert ds2.report.schema_valid


def test_case_mismatch_geometry_and_properties(tmp_path, built):
    m, fr = copy.deepcopy(built[0]), built[1]
    m["geometry"]["d_m"] = 12e-6
    m["fluid"]["mu_Pa_s"] *= 1.01
    for k in fr:
        fr[k] = fr[k].copy()
    fr["hydraulic"]["mu_Pa_s"] = m["fluid"]["mu_Pa_s"]
    fr["hydraulic"]["K_bulk_m2"] = (fr["hydraulic"].mu_Pa_s * fr["hydraulic"].Q_cell_m3_s * fr["hydraulic"].ell_m
                                    / (fr["hydraulic"].b_m * fr["hydraulic"].h_m * fr["hydraulic"].Delta_p_cell_Pa))
    ds = _write(tmp_path, "mismatch", m, fr)
    assert ds.report.schema_valid, ds.report.errors
    match_case(ds, base_case().case_requirements())
    assert ds.report.case_matched is False
    joined = " ".join(ds.report.case_mismatches)
    assert "d_m" in joined and "mu_Pa_s" in joined


def test_write_dataset_refuses_overwrite(tmp_path, built):
    _write(tmp_path, "once", *built)
    with pytest.raises(FileExistsError):
        write_dataset(tmp_path / "once", copy.deepcopy(built[0]), built[1])


def test_legacy_permeability_csv_is_not_a_dataset(tmp_path):
    src = REPO_ROOT / "data" / "processed" / "permeability_stage1_flow.csv"
    if not src.exists():
        pytest.skip("legacy processed table not present")
    d = tmp_path / "legacy"
    d.mkdir()
    shutil.copy(src, d / "hydraulic.csv")
    assert not load_dataset(d).report.schema_valid


# ------------------------------------------------------------------ import
def _raw_copy(tmp_path):
    raw = tmp_path / "raw"
    shutil.copytree(FIXTURE_DIR / "raw_export_fixture", raw)
    return raw


def test_import_raw_fixture(tmp_path):
    raw = _raw_copy(tmp_path)
    before = {p.name: sha256_file(p) for p in raw.iterdir()}
    res = import_from_manifest(raw / "export_manifest.yaml", tmp_path / "out", REPO_ROOT)
    r = res["report"]
    assert r["schema_valid"], r["errors"]
    assert r["verification_only"] and not r["numerically_qualified"]
    assert (tmp_path / "out" / "qualification_report.json").exists()
    assert before == {p.name: sha256_file(p) for p in raw.iterdir()}, "raw exports were modified"
    with pytest.raises(FileExistsError):
        import_from_manifest(raw / "export_manifest.yaml", tmp_path / "out", REPO_ROOT)


def test_import_refuses_unit_mismatch(tmp_path):
    raw = _raw_copy(tmp_path)
    f = raw / "flow_raw.csv"
    f.write_text(f.read_text().replace("Q, flow rate (m^3/s)", "Q, flow rate (L/min)"))
    with pytest.raises(ImportError_, match="no silent unit conversion"):
        import_from_manifest(raw / "export_manifest.yaml", tmp_path / "out", REPO_ROOT)


def test_import_refuses_fixture_under_cell_tables(tmp_path):
    raw = _raw_copy(tmp_path)
    with pytest.raises(ImportError_, match="reserved for physical"):
        import_from_manifest(raw / "export_manifest.yaml", tmp_path / "data" / "processed" / "cell_tables" / "x", REPO_ROOT)


def test_import_requires_explicit_admissibility(tmp_path):
    raw = _raw_copy(tmp_path)
    mp = raw / "export_manifest.yaml"
    m = yaml.safe_load(mp.read_text())
    m["exports"]["meniscus"]["constants"].pop("admissible")
    mp.write_text(yaml.safe_dump(m))
    with pytest.raises(ImportError_, match="never assumed"):
        import_from_manifest(mp, tmp_path / "out", REPO_ROOT)
