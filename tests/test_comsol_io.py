from pathlib import Path

import pytest

from nanofab.comsol_io import mesh_convergence_report, read_comsol_table, relative_change

EXAMPLES = Path(__file__).resolve().parents[1] / "data" / "comsol_exports" / "examples"


def test_read_csv_table_with_units_and_metadata():
    t = read_comsol_table(EXAMPLES / "stage1_flow_example.csv")
    assert list(t.data.columns) == ["hmax", "dp", "Surface Integration 1"]
    assert t.units["Surface Integration 1"] == "m^3/s"
    assert t.units["dp"] == "Pa"
    assert t.metadata["Model"] == "stage1_flow.mph"
    assert t.data["hmax"].iloc[0] == pytest.approx(4e-6)
    assert t.columns_with_unit("m^3/s") == ["Surface Integration 1"]


def test_read_comsol63_csv_with_quoted_headers(tmp_path):
    p = tmp_path / "stage1.csv"
    p.write_text(
        '% Model,permeability1.mph\n'
        '% Version,COMSOL 6.3.0.420\n'
        '% Date,"Oct 7 2026, 21:57"\n'
        '% hmax (m),"Q, flow rate (m^3/s)","K, permeability (m^2)"\n'
        "1.0E-6,2.63E-12,2.96E-11\n"
        "5.0E-7,2.64E-12,2.97E-11\n",
        encoding="utf-8",
    )
    t = read_comsol_table(p)
    assert list(t.data.columns) == ["hmax", "Q, flow rate", "K, permeability"]
    assert t.units["Q, flow rate"] == "m^3/s"
    assert t.units["K, permeability"] == "m^2"
    assert t.metadata["Model"] == "permeability1.mph"
    assert t.metadata["Date"] == "Oct 7 2026, 21:57"


def test_read_whitespace_txt_table():
    t = read_comsol_table(EXAMPLES / "stage4_thermal_example.txt")
    assert list(t.data.columns) == ["Pc", "Q_in", "Q_evap", "T_base", "T_lv"]
    assert t.units["T_base"] == "K"
    assert len(t.data) == 4
    assert t.column("q_evap") == "Q_evap"


def test_headerless_file_gets_generic_names(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("1 2\n3 4\n", encoding="utf-8")
    t = read_comsol_table(p)
    assert list(t.data.columns) == ["col0", "col1"]


def test_mesh_convergence_report():
    assert relative_change(1.0, 1.03) == pytest.approx(0.03)
    report = mesh_convergence_report({"coarse": 1.10, "medium": 1.02, "fine": 1.01}, tolerance=0.02)
    assert "converged" in report and "NOT" not in report
    report = mesh_convergence_report({"coarse": 1.5, "fine": 1.0}, tolerance=0.02)
    assert "NOT converged" in report
