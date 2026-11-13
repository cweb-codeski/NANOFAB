# NANOFAB — Wick Evaporator Modeling

Modeling framework for **two-phase evaporative cooling using micro/nanostructured wicks** for
high-heat-flux electronics (GPUs/TPUs). The goal is to compare a baseline periodic micropillar wick
against novel wick architectures and explain *why* one performs differently.

Causal chain the models are built around:

```text
wick geometry -> capillary pressure + permeability -> liquid transport
-> thin-film / interfacial evaporation -> heat removal -> surface temperature -> dryout limit
```

COMSOL is the primary physics tool. Python here supports it: parameter management with source
tracking, analytical sanity checks, and post-processing of COMSOL CSV exports.

## Layout

| Path | Contents |
| --- | --- |
| `docs/` | Modeling context, assumptions, equations, literature notes, COMSOL unit-cell guide |
| `parameters/` | YAML inputs; every entry has `value`, `unit`, `status`, `source` |
| `src/nanofab/` | Python library (fluids, capillary, permeability, dryout, evaporation, conduction, ...) |
| `analysis/` | Scripts: baseline estimates, permeability / heat transfer from COMSOL exports |
| `src/nanofab/device/` | Passive device-scale model (V2019 framework): table-driven Darcy transport + x-z substrate conduction, dryout search |
| `parameters/cases/` | Named device cases (manifests referencing `parameters/*.yaml` entries) |
| `src/nanofab/electroosmosis/` | 1D passive vs electroosmotically assisted wick (same K, same geometry); inputs in `parameters/electroosmosis/`, see `docs/electroosmosis_model.md` and `MODEL_STATUS.md` |
| `comsol/` | COMSOL model notes (`.mph` files are gitignored) |
| `data/` | Raw data, COMSOL exports, literature digitizations, processed results |
| `tests/` | pytest checks of analytical limits and identities |

## Setup

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

## Where to start

| Document | Purpose |
| --- | --- |
| `docs/comsol_unit_cell_guide.md` | **Canonical workflow**: geometry cases, Stages 1-7, build procedures, checks |
| `docs/modeling_audit_2026-11.md` | Prioritized audit, source/assumption register, blockers |
| `docs/validation_and_experiments.md` | Verification vs validation, literature benchmarks, experiment plan and matrix |
| `docs/stage2_undergrad_assignment.md` | Current undergraduate task (Stage 2) |
| `parameters/unit_cell/cases.yaml` | Reference (10/30/25 µm) and proposed (20/50/40 µm) cases, kept separate |

## Common commands

```powershell
# Hand-calculation equivalents of every COMSOL stage, reference vs proposed cell
python analysis/case_hand_checks.py

# Compare the two independent Stage 1 builds of the reference cell
python analysis/stage1_crosscheck.py

# Analytical sanity-check table for the baseline geometry
python analysis/baseline_estimates.py

# Analytical models vs Vaartstra, Lu & Wang (IJHMT 2019) reported values
python analysis/vaartstra_replication.py

# Permeability from a COMSOL Stage 1 export (--case when d/l/hp are not export columns)
python analysis/permeability_from_comsol.py data/comsol_exports/stage1_flow.csv --case reference_d10_l30_h25

# Effective evaporative h from a COMSOL Stage 2/4/5 export (needs an hmax column to write output)
python analysis/heat_transfer_from_comsol.py data/comsol_exports/stage4_thermal.csv --case reference_d10_l30_h25

# Passive device-scale model (docs/passive_device_model.md)
python analysis/device_benchmarks.py --suite all --output data/processed/device_runs/verification
python analysis/run_passive_device.py --case parameters/cases/vaartstra2019_fig4.yaml --preflight --output data/processed/device_runs/vaartstra_preflight
python analysis/run_passive_device.py --case parameters/cases/vaartstra2019_fig4.yaml --search-dryout --output data/processed/device_runs/vaartstra_prediction
#   (currently exits 3, MISSING_VALIDATED_INPUTS: no qualified physical cell tables exist yet)
python analysis/import_cell_tables.py --manifest PATH_TO_EXPORT_MANIFEST --output data/processed/cell_tables/DATASET_ID

# Passive vs electroosmotic wick (docs/electroosmosis_model.md) -> results/electroosmosis/
python analysis/permeability_from_comsol.py data/comsol_exports/stage1_flow.csv --dp 1 --mu 2.816e-4
python analysis/eo_assisted_wick.py
python analysis/eo_field_treatment_test.py   # uniform vs solved field (docs/eo_field_treatment_test.md)

# Tests
pytest
```

## Provenance

Every input in `parameters/*.yaml` has `value`, `unit`, `provenance`, `source` and
`source_location`. `provenance` is one of `PRIMARY_SOURCE`, `AUTHORITATIVE_DATABASE`,
`PROJECT_DERIVED`, `MODEL_CALCULATED`, `ASSUMED`, `DESIGN_TARGET`, `PROJECT_CONTEXT`,
`COMSOL_EXTRACTED` or `USER_PROVIDED`. `source`
is an id in `config/sources.yaml`. The loader rejects an invalid label, a missing source where
one is required, or an unregistered source id. Project design targets live only in
`parameters/design_targets.yaml`, kept apart from replication inputs. The definitions, the
authority order and the modeling rules are in `docs/SOURCE_AUTHORITY.md`. Anything `ASSUMED` that
affects performance should be swept rather than trusted (see `docs/assumptions.md`).

## Current status

Phase 1 (learn and reproduce). The reference case is the Vaartstra, Lu & Wang (IJHMT 136, 2019)
device: d = 10 µm, l = 30 µm, h = 25 µm, L = 5 mm, water at 100 °C (all `PRIMARY_SOURCE`).
The proposed wick (20/50/40 µm, user design inputs) is a separate case with its own runs.
Replication targets are in `data/literature/vaartstra2019_targets.csv`. Per-stage evidence status
is in `docs/comsol_unit_cell_guide.md` §0 and `MODEL_STATUS.md`.

The passive device model (`src/nanofab/device/`) is numerically verified on analytical and
manufactured inputs. A V2019 literature comparison is blocked until qualified unit-cell tables
(meniscus, hydraulic, thermal) exist and the Meland & Ytrehus accommodation source gap is closed;
see `docs/passive_device_implementation_report.md`. The legacy flat-meniscus result (100.6 W/cm²)
is an approximate analytical estimate, not a bound.
