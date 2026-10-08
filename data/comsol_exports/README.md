# COMSOL exports

Put CSV/TXT files exported from COMSOL here (Results > Tables > Export, or Results > Export > Data).
Large exports are gitignored; small curated examples can go in `examples/`.

Naming convention: `<stage>_<what>_<geometry-tag>.csv`, e.g. `stage1_flow_d10_l30_h40.csv`.

What each analysis script expects is documented in `docs/comsol_unit_cell_guide.md`
(section "Exporting results") and in the docstrings of `analysis/*.py`.
