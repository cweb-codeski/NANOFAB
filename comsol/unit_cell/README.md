# Unit-cell COMSOL models

`.mph` files are gitignored (large binaries). For every model, keep a row here, a text copy
(File > Compact History, then Save As > Model File for Java, `.java`, committed in this folder)
and the raw export in `data/comsol_exports/`. Build instructions: `docs/comsol_unit_cell_guide.md`.
Status labels are defined in guide §0.

## Models known to exist (reported; files not in the repository)

| File | Stage | Case | Status | What is known | Missing evidence |
| --- | --- | --- | --- | --- | --- |
| `permeability1.mph` (Chase) | 1 | reference | STORED_DERIVED | COMSOL 6.3.0.420; quarter cell, dp = 1 Pa, μ = 2.816e-4 Pa s, flat symmetry top; hmax 2/1/0.7/0.5 µm → K 2.9698e-11 m² (`data/processed/permeability_stage1_flow.csv`). Earlier runs in the same export gave K ~ 7.9e-12 ("setup error"; cause not recorded) | `.mph`/`.java`, selections, Stokes vs NS setting |
| student Stage 1 model | 1 | reference | REPORTED_UNVERIFIED | hmax 2/1/0.5/0.25 µm → K 2.9711e-11 m² (`data/reported/`) | file, export, μ, dp, BC options |
| `stage2_flat_interface.mph` (Chase) | 2 | reference | REPORTED_UNVERIFIED | COMSOL 6.3.0.420; two Q_in (1.125e-4 W) and two T_base (381.19639, 381.21757 K) columns without identifiers | hi value; pillar trimming (volume 490.9 µm³?); T_base selection and area; mesh/study per column; Q_evap |

## Planned models (one file per stage and case)

| File | Stage | Physics (base-package route) | Status |
| --- | --- | --- | --- |
| `ref_stage1_flow.mph`, `prop_stage1_flow.mph` | 1 | Laminar Flow, Stokes option | ref: see above; prop: PROPOSED |
| `ref_stage1q_periodic.mph`, `prop_stage1q_periodic.mph` | 1Q | Laminar Flow, periodic flow condition, NS vs Stokes, interface bracket | PROPOSED |
| `ref_stage2_flat.mph`, `prop_stage2_flat.mph` | 2 | Heat Transfer in Solids and Fluids | PROPOSED (assignment) |
| `ref_stage2A_eo.mph`, `prop_stage2A_eo.mph` | 2A | Electric Currents (or Laplace Equation) + Laminar Flow, EO slip | PROPOSED |
| `*_stage3_meniscus.mph` | 3, 3H, 3E | General Form Boundary PDE + deformed mesh + Stage 1/2A physics | PROPOSED |
| `*_stage4_evap_<closure>.mph` | 4 | Stage 2 physics on menisci, named closure | PROPOSED; V2019 replication blocked by the accommodation-treatment source gap |
| `*_stage5_sweep.mph` | 5 | sweeps for tables | PROPOSED |
| `*_stage6_device.mph` | 6 | Heat Transfer in Solids + boundary PDEs | PROPOSED; proposed case blocked by unresolved device inputs |

Record for each new row: COMSOL version, licensed modules used, case id, `hmax` list, DOF, the
export file name, and the checks of the stage (pass/fail with numbers).
