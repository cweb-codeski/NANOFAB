# Unit cell → whole wick property hand-off

The five-stage unit-cell COMSOL workflow (`docs/comsol_unit_cell_guide.md`, model notes in
`comsol/unit_cell/`) is a **property generator**. It is not the wick model. Its job is to fill
these tables for each cell type (one geometry + material):

| Table | Axes | From stage | Status |
| --- | --- | --- | --- |
| K | Pc | 1 (flat top), 3/5 (meniscus) | Stage 1 done for d10/l30/h25: K = 2.97e-11 m² (`comsol_stage1_K` recipe). K(Pc) pending. |
| h_e | Pc | 4/5 | pending; whole-wick runs use an ASSUMED placeholder (W18) |
| Pc_max(V) | closed form | 2γ cos θ_V / r_eff (Lippmann-Young); Stage 3 depinning curvature should confirm the V = 0 value | analytical |
| k_eff | - | Stage 2 (in-plane only) | analytical parallel bound; < 5% of in-plane conductance |

Electrowetting is handled in the whole-wick code, not in these tables. Voltage changes only the
receding angle, so K and h_e have **no** voltage axis, and CSV tables with a voltage or field
column are rejected (`docs/whole_wick_assumptions.md`, EW3). Electrical power is unresolved and
has no table.

## Steps

1. Run the Stage 5 parametric sweep over the capillary pressure Pc = p_v − p_l from 0 to depinning.
   For a cell type with a dielectric coating, extend the sweep to 2γ/r_eff, the Pc_max at full
   electrowetting (cos θ_V = 1). At the default 15° receding angle that is 3.5% above Pc_max(0).
2. Export the sweep table (Results > Tables > Export) to `data/comsol_exports/`.
3. Convert it:
   `python unit_cell/build_property_tables.py data/comsol_exports/stage5.csv --property h_e --pc-col dp --value-col h_e --out parameters/whole_wick/tables/si_d10_l30_h25_h_e.csv`
4. In `parameters/whole_wick/cell_types.yaml`, switch that recipe's `h_e` (or `K`) to
   `method: table` with `provenance: COMSOL_EXTRACTED` and a `source_location` naming the `.mph`
   file and export.
5. Re-run `python validation/run_validation.py`. The "unit-cell COMSOL" section then compares the
   table against the analytical values.

The table's Pc range must reach the largest Pc_max(V) used in any case. Otherwise lookups near
dryout are reported as extrapolation (clipped, never extrapolated).
