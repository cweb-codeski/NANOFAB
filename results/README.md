# results/

Generated output only. Regenerate rather than edit:

| Folder | Command | Contents |
| --- | --- | --- |
| `cases/` | `python whole_wick/run_cases.py` | `SUMMARY.md`; per heat map: `comparison.csv`, map PNGs (nominal and dryout onset), per-cell CSVs; `parameter_provenance.csv` |
| `sweeps/` | `python whole_wick/run_sweeps.py` | `SWEEPS.md` and one CSV per sweep |

Nothing here is validated (`validation/VALIDATION_REPORT.md`). Electrowetting numbers use
Lippmann-Young with an ASSUMED ε_r and saturation angle, and electrical power is UNRESOLVED.
Do not quote them as predictions.
