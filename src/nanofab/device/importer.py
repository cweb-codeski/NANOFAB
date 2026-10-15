"""Build a strict cell-table dataset from raw COMSOL exports described by a raw-export manifest.

The raw exports are parsed with ``nanofab.comsol_io.read_comsol_table`` (generic parser, unchanged).
This module adds the strict path: every mapped column must carry the expected unit in its COMSOL
header (no silent unit conversion), derived quantities are computed with stated formulas, raw
files are hashed and never modified, and the output directory must not already hold a dataset.

Derived quantities:
    K_bulk  = mu Q_cell ell / (b h Delta_p_cell)          (bulk normalization, A_ref = b h)
    q_in    = Q_in / A_footprint
    h_e     = q_in / (T_base - T_v)                       (V2019 Eq. 4, base-plane reference)
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from ..comsol_io import read_comsol_table
from .model import DataRole
from .tables import COLUMN_UNITS, SCHEMA_VERSION, load_dataset, sha256_file, write_dataset

RAW_SCHEMA_VERSION = "nanofab.raw_export_manifest/1"
COMSOL_UNIT_MAP = {"m^3/s": "m3/s", "Pa": "Pa", "W": "W", "K": "K", "m": "m", "m^2": "m2", "deg": "deg",
                   "°": "deg", "1": "-", "": "-"}
QUALIFIED_DIR_PARTS = ("data", "processed", "cell_tables")


class ImportError_(ValueError):
    pass


def _column(table, raw_name: str) -> str:
    cols = [c for c in table.data.columns if c == raw_name]
    if len(cols) != 1:
        raise ImportError_(f"{table.path}: expected exactly one column named {raw_name!r}, have "
                           f"{list(table.data.columns)}")
    return cols[0]


def _read_mapped(base: Path, spec: dict, repo_root: Path) -> tuple[pd.DataFrame, dict]:
    fp = (base / spec["file"]).resolve()
    if not fp.is_file():
        raise ImportError_(f"raw export {fp} not found")
    tbl = read_comsol_table(fp)
    out = {}
    for target, raw in (spec.get("columns") or {}).items():
        if target not in COLUMN_UNITS:
            raise ImportError_(f"unknown target column {target!r}")
        c = _column(tbl, raw)
        unit = COMSOL_UNIT_MAP.get(tbl.units.get(c, ""), None)
        if unit != COLUMN_UNITS[target]:
            raise ImportError_(f"{fp.name}: column {raw!r} has unit {tbl.units.get(c)!r}; {target} requires "
                               f"{COLUMN_UNITS[target]!r} (no silent unit conversion)")
        out[target] = tbl.data[c].to_numpy(float)
    n = len(next(iter(out.values())))
    for target, val in (spec.get("constants") or {}).items():
        if target not in COLUMN_UNITS:
            raise ImportError_(f"unknown constant column {target!r}")
        out[target] = np.full(n, val, dtype=object if isinstance(val, str) else float)
    try:
        rel = fp.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        rel = str(fp)
    return pd.DataFrame(out), {"path": rel, "sha256": sha256_file(fp), "comsol_metadata": tbl.metadata}


def import_from_manifest(manifest_path: str | Path, out_dir: str | Path, repo_root: Path) -> dict:
    manifest_path = Path(manifest_path)
    raw = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    if raw.get("schema_version") != RAW_SCHEMA_VERSION:
        raise ImportError_(f"raw-export manifest schema_version must be {RAW_SCHEMA_VERSION!r}")
    meta = copy.deepcopy(raw.get("metadata") or {})
    role = meta.get("data_role")
    out_dir = Path(out_dir)
    resolved = out_dir.resolve()
    if role != DataRole.PHYSICAL_CELL_CALCULATION.value:
        parts = resolved.parts
        for i in range(len(parts) - 2):
            if tuple(p.lower() for p in parts[i:i + 3]) == QUALIFIED_DIR_PARTS:
                raise ImportError_(f"data_role {role} may not be written under data/processed/cell_tables "
                                   "(reserved for physical cell calculations)")
    base = manifest_path.parent
    frames, exports = {}, []
    ex = raw.get("exports") or {}
    for name in ("hydraulic", "thermal", "meniscus"):
        if name not in ex:
            raise ImportError_(f"exports.{name} missing")
        df, info = _read_mapped(base, ex[name], repo_root)
        info["table"] = name
        exports.append(info)
        frames[name] = df
    hyd = frames["hydraulic"]
    hyd["K_bulk_m2"] = hyd["mu_Pa_s"] * hyd["Q_cell_m3_s"] * hyd["ell_m"] / (
        hyd["b_m"] * hyd["h_m"] * hyd["Delta_p_cell_Pa"])
    frames["hydraulic"] = hyd[["Pc_Pa", "K_bulk_m2", "Q_cell_m3_s", "Delta_p_cell_Pa", "ell_m", "b_m", "h_m",
                               "mu_Pa_s"]]
    th = frames["thermal"]
    th["q_in_W_m2"] = th["Q_in_W"] / th["A_footprint_m2"]
    th["h_e_W_m2K"] = th["q_in_W_m2"] / (th["T_base_K"] - th["T_v_K"])
    frames["thermal"] = th[["Pc_Pa", "h_e_W_m2K", "Q_in_W", "Q_evap_W", "A_footprint_m2", "T_base_K", "T_v_K",
                            "q_in_W_m2"]]
    men = frames["meniscus"]
    spec_m = ex["meniscus"]
    men.insert(1, "state_id", [f"{spec_m.get('state_id_prefix', 'state')}{i:03d}" for i in range(len(men))])
    if "admissible" not in men.columns:
        raise ImportError_("exports.meniscus must supply 'admissible' (a column or an explicit constant); "
                           "admissibility is never assumed")
    frames["meniscus"] = men

    gen = meta.setdefault("generation", {})
    gen["export_files"] = [{"path": e["path"], "sha256": e["sha256"], "table": e["table"]} for e in exports]
    gen["raw_export_manifest"] = {"path": str(manifest_path.as_posix()), "sha256": sha256_file(manifest_path)}
    gen["comsol_export_metadata"] = {e["table"]: e["comsol_metadata"] for e in exports}
    manifest = {"schema_version": SCHEMA_VERSION, "dataset_id": raw["dataset_id"], **meta,
                "tables": {name: {} for name in frames}}
    mpath = write_dataset(out_dir, manifest, frames)
    ds = load_dataset(out_dir, repo_root=repo_root)
    report = ds.report.as_dict()
    (out_dir / "qualification_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return {"manifest": str(mpath), "report": report}
