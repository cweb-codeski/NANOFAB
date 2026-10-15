# Cell-Table Contract (`nanofab.cell_tables/1`)

Implemented in `src/nanofab/device/tables.py` (validation, qualification, case match) and
`src/nanofab/device/importer.py` (raw COMSOL export -> dataset). Tests: `tests/test_cell_tables.py`.

A dataset is a directory holding `manifest.yaml` plus three CSV tables (`hydraulic.csv`,
`thermal.csv`, `meniscus.csv`). Three outcomes are reported **separately**:

1. **Schema validity**: files well formed, hashed, unit-labelled, conventions stated, values
   finite and physically signed, internal consistency (for example K_bulk recomputed from Q_cell).
2. **Numerical qualification**: evidence that the cell calculations are converged and applicable
   (meshes, gradient sweep, energy balance, demonstrated endpoint, raw exports and run log).
3. **Case match**: geometry, fluid, properties, vapor state, solid conductivity, substrate,
   kinetic model and limit criterion agree with the named case.

A user-entered `validated: true` is ignored (a warning is recorded). Relabelling a file never
promotes it.

## Data roles

| data_role | Admitted to a literature case | Notes |
| --- | --- | --- |
| PHYSICAL_CELL_CALCULATION | yes, if qualified and matched | needs registered source ids |
| LITERATURE_DIGITIZATION | no (comparison only) | needs registered source ids |
| ANALYTICAL_EXAMPLE | never (VERIFICATION ONLY) | must record `generating_expressions` |
| MANUFACTURED_FIXTURE | never (VERIFICATION ONLY) | must record `generating_expressions`; lives under `tests/fixtures/device/` |

The importer refuses to write non-physical roles under `data/processed/cell_tables/`.

## Manifest keys (all required)

`schema_version`, `dataset_id`, `data_role`, `provenance` (existing labels; COMSOL_EXTRACTED or
MODEL_CALCULATED for qualification), `source_ids`, `source_locations`, `geometry` (d_m, l_m, h_m,
geometry_id, surface_id, array, receding_contact_angle_deg), `fluid` (name, property_source,
property_version, T_ref_K, rho, mu, sigma, h_fg, k, cp), `solid` (name, k_W_mK, k_status),
`vapor` (composition = pure_saturated_vapor, T_v_K, p_v_Pa), `generation` (method, software,
software_version, model_file, export_files with sha256, run_log), `conventions`, `cell`
(flow_length_ell_m, represented_width_b_m, symmetry_factor, footprint_area_m2,
thermal_cell_substrate_m), `boundary_conditions` (hydraulic, thermal, meniscus), `kinetic_model`
(name, accommodation_coefficient, accommodation_treatment, accommodation_treatment_verified),
`pressure_interval_Pa`, `physical_endpoint` (Pc_Pa, criterion, demonstrated, evidence),
`geometry_limit` (Pc_Pa, criterion, demonstrated), `mesh_evidence` (per table: levels,
finest_pair_rel_change, systematic), `gradient_study` (Pc_values_Pa, gradient_range_Pa_m,
max_rel_K_variation, stokes_vs_navier_stokes_rel_diff), `supported_ranges`
(pressure_gradient_Pa_m, heat_flux_W_m2, T_v_K), `energy_balance_tolerance`, `limitations`,
`tables` (per table: file, sha256, columns -> unit; optional pore_to_bulk_area_ratio).

Required conventions (exact strings; no automatic reconciliation):

| key | value |
| --- | --- |
| pressure | `Pc = p_v - p_liquid [Pa], positive suction` |
| hydraulic_normalization | `bulk: Q = (W h K_bulk / mu) dPc/dx, A_ref = W h` |
| thermal_reference_plane | `pillar_base_plane` |
| area_basis | `footprint` |
| curvature_normal | free text stating the V2019 Eq. (2) normal convention |

## Columns and units

| Table | Column | Unit | Required |
| --- | --- | --- | --- |
| hydraulic | Pc_Pa | Pa | yes |
| hydraulic | K_bulk_m2 or tau_m3 (or K_pore_m2 with a stated pore_to_bulk_area_ratio) | m2 / m3 | one of |
| hydraulic | Q_cell_m3_s, Delta_p_cell_Pa, ell_m, b_m, h_m, mu_Pa_s | m3/s, Pa, m, m, m, Pa*s | all, when COMSOL-generated |
| thermal | Pc_Pa, h_e_W_m2K, Q_in_W, Q_evap_W, A_footprint_m2, T_base_K, T_v_K, q_in_W_m2 | Pa, W/m2/K, W, W, m2, K, K, W/m2 | yes |
| meniscus | Pc_Pa, state_id, branch_id, min_clearance_m, interface_area_m2, contact_angle_measure_deg, admissible | Pa, -, -, m, m2, deg, - | yes |
| meniscus | A_liquid_section_m2 + section_location | m2, - | optional (needed for intrinsic velocity) |

Derived definitions checked on load (relative tolerance 1e-6):
`K_bulk = mu Q_cell ell / (b h Delta_p_cell)`, `tau = h K_bulk`, `q_in = Q_in / A_footprint`,
`h_e = Q_in / (A_footprint (T_base - T_v))`. A bulk-normalized K is never multiplied by a further
liquid-area or saturation factor (that would count meniscus recession twice). Multiple gradient or
mesh studies are kept out of the final one-dimensional table; they are summarized in
`mesh_evidence` and `gradient_study`, never averaged silently.

## Rejections (schema errors)

Nonfinite or nonnumeric values; nonpositive K, tau, h_e, Q_in, interface area or admissible
clearance; negative Pc (sign convention); unsorted or duplicate pressures; duplicate state ids;
unknown columns or units; undeclared columns; hash mismatch; missing keys; convention mismatch;
area-basis mismatch; COMSOL-column inconsistency with K; T_base <= T_v; energy imbalance above the
declared tolerance; a zero-load thermal row (h_e would be 0/0: use a documented small-load limit);
rows not spanning `pressure_interval_Pa`; unknown data_role or provenance; unregistered source ids.

## Qualification targets (predeclared; `QUALIFICATION_TARGETS`)

- At least 3 mesh levels per table; finest-pair relative change <= 1e-3; systematic convergence.
  (V2019 Table A1 reports finest-pair changes of 1.1e-4 interface minimum, 1.4e-4 cell HTC,
  5.5e-4 cell K, 2.3e-6 device pressure, 9.9e-7 device temperature range: reported observations,
  not guaranteed error bounds.)
- Gradient sweep at >= 3 representative Pc values spanning the interval, covering device
  gradients, with max K variation <= 1%, plus a Stokes vs Navier-Stokes comparison. If the
  variation exceeds 1%, a one-dimensional K(Pc) table is not qualified and Q(Pc, gradient) data
  are required.
- Thermal energy-balance tolerance <= 1e-3.
- Raw exports present and hash-matched; `.mph` model file named; actual solver run log.
- Physical endpoint demonstrated with evidence (the last CSV row is not evidence).

## Case match (`match_case`)

Geometry within 1e-9 relative (mixed geometries are never interpolated); same fluid; T_v and
property reference temperature within 0.01 K; rho, mu, sigma, h_fg within 1e-3 (property-source
mismatch); solid k within 1e-3 of the case substrate; same thermal cell substrate thickness;
for V2019 replication, `kinetic_model.name = linearized_BTE_V2019_Eqs5-7` and
`accommodation_treatment_verified: true`; a demonstrated endpoint within 1e-3 of the case's Eq. (3)
limit.

## Interpolation and range

Piecewise-linear in tau and h_e (positivity preserved), exact at knots, no extrapolation, no
clamping. A narrow snap (1e-9 x pressure span) absorbs floating-point error at the ends only. The
usable range is the intersection of the hydraulic, thermal and meniscus ranges, cut at the last
admissible meniscus state.

## Raw-export manifest (`nanofab.raw_export_manifest/1`)

```yaml
schema_version: nanofab.raw_export_manifest/1
dataset_id: DATASET_ID
metadata: { ...all dataset keys except schema_version, dataset_id, tables... }
exports:
  hydraulic: {file: flow.csv, columns: {Pc_Pa: "Pc", Q_cell_m3_s: "Q", Delta_p_cell_Pa: "dp"},
              constants: {ell_m: ..., b_m: ..., h_m: ..., mu_Pa_s: ...}}
  thermal:   {file: thermal.csv, columns: {Pc_Pa: ..., Q_in_W: ..., Q_evap_W: ..., T_base_K: ...},
              constants: {A_footprint_m2: ..., T_v_K: ...}}
  meniscus:  {file: meniscus.csv, state_id_prefix: s,
              columns: {Pc_Pa: ..., min_clearance_m: ..., interface_area_m2: ..., contact_angle_measure_deg: ...},
              constants: {branch_id: ..., admissible: 1}}
```

Column headers are parsed by `nanofab.comsol_io.read_comsol_table` (unchanged generic parser).
The unit in each COMSOL header must equal the contract unit (`m^3/s` -> m3/s, `Pa`, `W`, `K`, `m`,
`m^2` -> m2, `deg`); there is no silent conversion. Admissibility must be supplied explicitly.
Raw files are hashed and never modified; an existing dataset directory is never overwritten.

```powershell
python analysis/import_cell_tables.py --manifest PATH_TO_EXPORT_MANIFEST --output data/processed/cell_tables/DATASET_ID
# runnable demonstration with the MANUFACTURED raw-export fixture (VERIFICATION ONLY; refused under cell_tables/):
python analysis/import_cell_tables.py --manifest tests/fixtures/device/raw_export_fixture/export_manifest.yaml --output data/processed/device_runs/import_demo/fixture_raw_import
```

## Existing data and why none of it qualifies

See `config/data_inventory.yaml` (reported by every literature preflight). In short:
`data/processed/permeability_stage1_flow.csv` is a stored derived COMSOL result with incomplete
provenance (flat interface Pc = 0 only, one low gradient ~6.7e4 Pa/m vs device gradients of order
1e5-1e6 Pa/m, no Stokes/Navier-Stokes comparison, `.mph` not tracked); the `examples/` exports
are SYNTHETIC; no thermal or meniscus tables exist.
