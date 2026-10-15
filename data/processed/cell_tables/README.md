# Qualified physical cell tables

Reserved for `PHYSICAL_CELL_CALCULATION` datasets imported with `analysis/import_cell_tables.py`
(`docs/cell_table_contract.md`). Fixtures and analytical examples are refused here; they live in
`tests/fixtures/device/`.

Status: **empty**. No physical meniscus, hydraulic or thermal cell tables exist yet. Each future
dataset is a new `DATASET_ID/` directory (never overwritten) containing `manifest.yaml`, the three
tables and `qualification_report.json`.
