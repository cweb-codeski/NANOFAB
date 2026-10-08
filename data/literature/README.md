# Literature data

Digitized curves and tabulated values from papers used as replication targets.
One CSV per figure/table, named `<firstauthor><year>_<fig-or-table>.csv`, with a header comment
giving the full citation, figure number, and digitization tool.

Every row carries `provenance`, `source` (an id in `config/sources.yaml`) and `source_location`.
Only values stated in a paper's text or tables are `PRIMARY_SOURCE`. Digitized points are
`PROJECT_DERIVED` with a note such as "Digitized from Fig. 3a". `MODEL_CALCULATED` values never
go here (enforced by `tests/test_provenance.py`); see `docs/SOURCE_AUTHORITY.md`.
