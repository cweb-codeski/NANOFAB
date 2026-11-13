# COMSOL exports

Put CSV/TXT files exported from COMSOL here (Results > Tables > Export, or Results > Export > Data).
Large exports are gitignored; small curated examples can go in `examples/` (those are SYNTHETIC).

Naming convention: `<case>_<stage>_<what>.csv`, with the case prefix `ref` (10/30/25 µm) or `prop`
(20/50/40 µm), e.g. `ref_stage2_flat.csv`, `prop_stage2A_eo.csv`. Older files keep their names.

Rules: never edit a raw export (add a corrected copy with a new name and say why); one row per
case × run × mesh with an `hmax` column; header lines name the model file, COMSOL version, date,
dataset, selections and expressions. Values a person reports without a raw export go to
`data/reported/` instead.

Column requirements per stage: `docs/comsol_unit_cell_guide.md` (§5 Stage 1, §6 Stage 2, §7
Stage 2A) and the docstrings of `analysis/*.py`.
