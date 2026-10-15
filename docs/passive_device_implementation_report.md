# Passive device-scale wick model: implementation report

Branch `cursor/source-provenance-system`, base commit `b12675f4` (changes uncommitted, not pushed).
Run date 2026-10-14 (local) / 2026-10-15 UTC. The specification is `docs/passive_device_model.md`;
the cell-table contract is `docs/cell_table_contract.md`; assumptions are in
`config/device_assumptions.yaml` (D1-D19).

Claim vocabulary used here, kept distinct:

- **numerically verified for these cases** -- the discretized equations reproduce exact or
  independently computed solutions of the *same* equations within predeclared gates (fixtures only).
- **approximate analytical estimate** -- closed-form or legacy estimates (`analysis/baseline_estimates.py`,
  `analysis/vaartstra_replication.py`).
- **literature comparison** -- comparison against V2019 reported values; not performed for the
  device model because validated inputs are missing (below).
- **experimentally validated** -- not claimed for anything in this repository.

## 1. Outcome

| item | outcome |
|---|---|
| device solver (Kirchhoff hydraulics + x-z conduction + Picard coupling + endpoint dryout solve) | implemented; numerically verified for these cases (109/109 predeclared gates pass) |
| V2019 Fig. 4 device prediction | **BLOCKED: MISSING_VALIDATED_INPUTS** (exit 3); nothing solved, no fallback used |
| COMSOL | not run (no software access or execution evidence); checklist in `docs/comsol_unit_cell_guide.md` Sec. 10 |
| audit corrections (a)-(f) | applied (Sec. 6) |
| 82 W/cm^2 | used only as a comparison reference scale; never an input or fit target |

## 2. Changed files

New:

- `src/nanofab/device/` -- `model.py` (statuses, exit codes, dataclasses), `tables.py` (strict cell-table
  loading/validation, piecewise-linear closures, `knot_stride`), `importer.py` (raw-export import; no unit
  conversion, admissibility never assumed), `hydraulics.py`, `conduction.py`, `solver.py` (Picard
  coupling, endpoint solve, continuation, `solve_prescribed_load`), `case.py` (case manifests, preflight,
  CoolProp liquid properties, V2019 Eq. 3), `output.py` (run directories, manifests, CSV/plots),
  `verification.py` (predeclared gates and suites), `__init__.py`.
- `analysis/run_passive_device.py`, `analysis/device_benchmarks.py`, `analysis/import_cell_tables.py`.
- `parameters/cases/vaartstra2019_fig4.yaml`; `tests/fixtures/device/` (manufactured case, tables, raw export).
- `config/device_assumptions.yaml`, `config/data_inventory.yaml`.
- `docs/passive_device_model.md`, `docs/cell_table_contract.md`, this report,
  `data/processed/cell_tables/README.md`, `data/processed/device_runs/README.md`.
- Tests: `tests/test_device_analytic.py`, `test_device_coupling.py`, `test_device_dryout.py`,
  `test_cell_tables.py`, `test_device_provenance.py`.

Modified (audit corrections and documentation): `README.md`, `analysis/baseline_estimates.py`,
`analysis/vaartstra_replication.py`, `analysis/heat_transfer_from_comsol.py`, `config/model_provenance.yaml`,
`config/sources.yaml`, `docs/SOURCE_AUTHORITY.md`, `docs/assumptions.md`, `docs/comsol_unit_cell_guide.md`,
`docs/equations.md`, `docs/literature_notes.md`, `parameters/materials.yaml`, `src/nanofab/evaporation.py`,
`tests/test_vaartstra2019.py`; one phrase in the untracked user file `docs/whole_wick_assumptions.md` (W3).
`README.md`, `config/sources.yaml` and `docs/SOURCE_AUTHORITY.md` already carried uncommitted user edits
at session start (git snapshot 2026-10-15T02:42Z); this work added to them without reverting anything.

Not touched (pre-existing user work): `.gitignore`, `src/nanofab/params.py`, `src/nanofab/provenance.py`,
`src/nanofab/whole_wick/`, `whole_wick/`, `unit_cell/`, `validation/`, `results/`, `parameters/whole_wick/`,
`tests/test_whole_wick.py`. `src/nanofab/comsol_io.py` is unchanged; the strict import path is
`src/nanofab/device/importer.py`.

## 3. Equations and closures actually implemented

Half device 0 <= x <= L (inlet x = 0 insulated in conduction, symmetry x = L), Pc = p_v - p_liquid,
tau(Pc) = h K_bulk(Pc).

- Hydraulics (Darcy with superficial velocity, V2019 framework): Q = -(W tau/mu) dp_l/dx = (W tau/mu) dPc/dx,
  dQ/dx = W q_e(x)/(rho_l h_fg). Kirchhoff variable Psi(Pc) = int tau dPc, exact for piecewise-linear tau
  (quadratic inverse per segment). Evaporation-moment identity Psi(Pc(L)) - Psi(Pc_feed) =
  (mu/(rho h_fg)) int_0^L int_x^L q_e dx' dx gives the endpoint solve with an exact rescale.
- Interface: q_e = h_e(Pc) (T_e - T_v), h_e from the cell thermal table (dataset-defined reference).
- Substrate: steady conduction k_s grad^2 T = 0 on x-z, bottom flux q''_in(x), top Robin with h_e(Pc(x)),
  cell-centred finite volumes, half-cell resistance at the Robin face, face-reconstructed T_e and T_source.
- Coupling: Picard on h_e(Pc(x)) <-> q_e(x), with predeclared tolerances.
- Dryout: the endpoint-constrained load search for Pc(L) = Pc_cap, where Pc_cap is the minimum of the
  table-demonstrated meniscus geometry limit, the physical limit (V2019 Eq. 3 for the literature case,
  2176.58 Pa from CoolProp water at 373.15 K), and table coverage. Every evaluation stays within
  [Pc_feed, Pc_cap]. Forward continuation fractions certify the bracket through monotonicity.
- Status taxonomy and exit codes: ADMISSIBLE / PHYSICAL_CAPILLARY_LIMIT / MENISCUS_GEOMETRY_LIMIT 0,
  config error 1, argparse 2, MISSING_VALIDATED_INPUTS 3, OUTSIDE_MODEL_REGIME 4,
  TABLE_COVERAGE_LIMIT 5, NUMERICAL_FAILURE 6, prescribed load above limit 7.
- The transport interface takes closures K(Pc), h_e(Pc) and source terms only; no electrical,
  boiling, graded-geometry or moving-interface terms exist.

## 4. Commands actually run and outcomes

All commands used `.\.venv\Scripts\python.exe` (Python 3.12, Windows).

| command | exit | outcome |
|---|---|---|
| `python -m pip install -e ".[dev]"` | 0 | installed |
| `python -m pytest -q` | 0 | 163 passed (53 s) |
| `python analysis/baseline_estimates.py` | 0 | approximate analytical estimates printed |
| `python analysis/vaartstra_replication.py` | 0 | legacy approximate analytical comparison; output CSV unchanged byte-for-byte |
| `python analysis/device_benchmarks.py --suite all --output data/processed/device_runs/verification` | 0 | ALL GATES PASSED, 109 gates; run `20261015T031857870284Z_b12675f4_bench_all` |
| `python analysis/run_passive_device.py --case parameters/cases/vaartstra2019_fig4.yaml --preflight --output data/processed/device_runs/vaartstra_preflight` | 3 | MISSING_VALIDATED_INPUTS; `missing_inputs.json` written |
| `python analysis/run_passive_device.py --case parameters/cases/vaartstra2019_fig4.yaml --search-dryout --output data/processed/device_runs/vaartstra_prediction` | 3 | MISSING_VALIDATED_INPUTS; nothing solved |
| `python analysis/run_passive_device.py --case tests/fixtures/device/manufactured_case.yaml --search-dryout --output data/processed/device_runs/verification_fixture_run` | 0 | VERIFICATION ONLY; PHYSICAL_CAPILLARY_LIMIT in [94.2631, 94.2678] W/cm^2 |
| `python analysis/import_cell_tables.py --manifest tests/fixtures/device/raw_export_fixture/export_manifest.yaml --output data/processed/device_runs/import_demo/fixture_raw_import` | 0 | schema-valid, verification only, not qualified (15 listed reasons) |
| `python analysis/import_cell_tables.py --manifest PATH_TO_EXPORT_MANIFEST --output data/processed/cell_tables/DATASET_ID` | 1 | refused: manifest not found (placeholder template) |
| same fixture manifest with `--output data/processed/cell_tables/fixture_should_refuse` | 1 | refused: fixtures may not be written under `data/processed/cell_tables`; nothing created |

A second identical benchmark run (`20261015T031803007999Z_b12675f4_bench_all`, started 03:18:03 UTC,
also all gates passed) exists in the same directory; it was not started by this session's command
sequence and was left in place (historical runs are not deleted).

Template command for real COMSOL data (not runnable until exports exist):

```
python analysis/import_cell_tables.py --manifest data/raw/comsol/<EXPORT_ID>/export_manifest.yaml --output data/processed/cell_tables/<DATASET_ID>
```

## 5. Numerical evidence (numerically verified for these cases)

All on manufactured or analytical fixtures (labelled VERIFICATION ONLY / ANALYTICAL_EXAMPLE).

- Analytic (constant K, h_e = 6e4): dryout load 1002907.9221088 vs closed form 1002907.9221079 W/m^2
  (rel 8.7e-13); tabulated tau with constant h_e vs quadrature, rel 8.7e-13; profile errors ~1e-12;
  balances 1e-13 to 1e-16.
- Accounting: results invariant to the width W (W/5) to ~1e-15; mirrored-domain accounting consistent.
- Zero load: exact zeros.
- Hydraulics vs an independent DOP853 ODE in Pc with q = q0(1 + 0.5 cos(pi x/L)): same-input relative
  error <= 3.0e-12; continuous-q errors 0.0910 / 0.0228 / 0.00569 Pa, observed order 2.0001 / 2.00002.
- Thermal MMS (theta = A + B(t_s - zeta) + C cos(m pi x/L) exp(m pi zeta/L)): T_e errors 1.22e-3 /
  3.06e-4 / 7.66e-5 K (orders 1.9955, 1.9989); T_source orders 1.9980, 1.9995; energy residual <= 1.5e-12.
- Coupled fixture, grids (50,8)/(100,16)/(200,32): endpoint loads 942682.08 / 942678.57 / 942677.69 W/m^2
  (finest-pair change 9.3e-7, monotone); fixed load 4.7e5 W/m^2 max suction 1037.4466 / 1037.4487 /
  1037.4492 Pa; tolerance pair rtol 1e-8 vs 1e-10 differs by 1.4e-14 in load (the endpoint iteration already
  contracts well below 1e-8); evaporation feedback q_e ptp/mean 0.0254, Pc_end 1037.449 Pa coupled vs
  1035.385 Pa with uniform evaporation; knot_stride 2 sampling changes the load by 2.2e-4 (relative) and
  T_e by 3.3e-3 K (reported, not gated).
- Status suite: base PHYSICAL_CAPILLARY_LIMIT; coverage ends at 1500 Pa -> TABLE_COVERAGE_LIMIT
  (certified lower point 66.79 W/cm^2); geometry limit at 1800 Pa -> MENISCUS_GEOMETRY_LIMIT (78.99 W/cm^2);
  narrow gradient range and narrow heat-flux range -> OUTSIDE_MODEL_REGIME; injected NaN -> NUMERICAL_FAILURE;
  negative load -> OUTSIDE_MODEL_REGIME.
- Reference scales at 82 W/cm^2 (comparison only): U_feed 0.07584 m/s, pore velocity 0.08309 m/s,
  Re_h 6.45 / 7.07, Pe_vertical 0.0566, Pe_L 2479 (pore-velocity based; superficial 2262), Ja 0.0187
  (assumed 10 K superheat). All match the audit within 0.1%.

These results verify that the code solves its stated equations. They say nothing about whether
those equations, or any closure, represent the V2019 device.

## 6. Audit corrections applied

(a) The "flat-meniscus upper bound" claim was removed; the legacy test is renamed
`test_legacy_flat_meniscus_sanity_check`. The moment identity shows the uniform-flux model is not a
bound in general. (b) The reversed inequality was corrected in `docs/SOURCE_AUTHORITY.md`. (c) Kinetic
resistance share is reported as a conditional lower bound (">= ~16%"), not a fraction.
(d) "Paper-exact" wording was replaced by "reduced BTE (V2019 Eqs. 5, 7; Eq. 6 unused)"; Schrage is a
cross-model check only. (e) Silicon k = 148 W/m K is relabelled ASSUMED (not a V2019 value).
(f) COMSOL post-processing reports *apparent* resistances with stated averaging definitions.

## 7. Discrepancies from V2019, approximations, unresolved sources

- No device-level comparison with V2019 is possible yet; Fig. 4 profiles are not digitized
  (non-blocking for the solver).
- Approximations (register D1-D19): steady state; uniform geometry; Darcy with superficial velocity and
  K(Pc) from cells (inertia only diagnosed, Re_h ~ 6-7); a 2D x-z substrate with a Robin top condition
  using a cell h_e; no lateral conduction through the liquid layer beyond the cell closure; constant
  liquid properties at T_v; inlet insulated; table interpolation piecewise-linear.
- Unresolved sources: Meland & Ytrehus 2003 (DOI 10.1063/1.1564097; NOT_HELD) for the accommodation
  correction used in the kinetic model; it is recorded as a blocking source gap for the literature case.
  Pujahari et al. 2022 (Langmuir 38, 8442) is held but not read; its SI is not held (future audit only).
- Silicon k at the operating temperature is ASSUMED (148 W/m K) and marked PROVISIONAL_INPUT.

## 8. Missing physical inputs (block the literature prediction)

1. Hydraulic cell table K_bulk(Pc) (or Q(Pc, gradient)) for the V2019 pillar geometry with mesh evidence.
2. Meniscus table: shape, admissibility and the geometry limit vs Pc.
3. Thermal table h_e(Pc) with the kinetic boundary condition (requires item 4).
4. Verified accommodation treatment (Meland & Ytrehus 2003) or an explicit decision to use V2019's
   stated form only, recorded as a source decision.

## 9. Next COMSOL action

Run the Sec. 10 checklist in `docs/comsol_unit_cell_guide.md` for the V2019 Fig. 4 geometry: compute the
meniscus family over Pc in [0, Pc_limit], solve the Stokes unit-cell flow at three mesh levels per Pc
(plus a Navier-Stokes spot check and a gradient study), and the thermal cell with the reduced-BTE interface.
Export the raw CSVs with `export_manifest.yaml` under `data/raw/comsol/<EXPORT_ID>/`, run the import
command above, then point `cell_tables.dataset` in `parameters/cases/vaartstra2019_fig4.yaml` at the
qualified dataset and rerun `--preflight` and `--search-dryout`.
