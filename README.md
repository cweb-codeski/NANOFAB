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
| `comsol/` | COMSOL model notes (`.mph` files are gitignored) |
| `data/` | Raw data, COMSOL exports, literature digitizations, processed results |
| `tests/` | pytest checks of analytical limits and identities |

## Setup

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

## Common commands

```powershell
# Analytical sanity-check table for the baseline geometry
python analysis/baseline_estimates.py

# Analytical models vs Vaartstra, Lu & Wang (IJHMT 2019) reported values
python analysis/vaartstra_replication.py

# Permeability from a COMSOL Stage 1 export
python analysis/permeability_from_comsol.py data/comsol_exports/stage1_flow.csv

# Effective evaporative h from a COMSOL Stage 4/5 export
python analysis/heat_transfer_from_comsol.py data/comsol_exports/stage4_thermal.csv

# Tests
pytest
```

## Provenance

Every input in `parameters/*.yaml` has `value`, `unit`, `provenance`, `source` and
`source_location`. `provenance` is one of `PRIMARY_SOURCE`, `AUTHORITATIVE_DATABASE`,
`PROJECT_DERIVED`, `MODEL_CALCULATED`, `ASSUMED`, `DESIGN_TARGET` or `PROJECT_CONTEXT`. `source`
is an id in `config/sources.yaml`. The loader rejects an invalid label, a missing source where
one is required, or an unregistered source id. Project design targets live only in
`parameters/design_targets.yaml`, kept apart from replication inputs. The definitions, the
authority order and the modeling rules are in `docs/SOURCE_AUTHORITY.md`. Anything `ASSUMED` that
affects performance should be swept rather than trusted (see `docs/assumptions.md`).

## Current status

Phase 1 (learn and reproduce). The baseline is the Vaartstra, Lu & Wang (IJHMT 136, 2019)
device: d = 10 µm, l = 30 µm, h = 25 µm, L = 5 mm, water at 100 °C (all `PRIMARY_SOURCE`).
Replication targets are in `data/literature/vaartstra2019_targets.csv`. Start with
`docs/comsol_unit_cell_guide.md`.
