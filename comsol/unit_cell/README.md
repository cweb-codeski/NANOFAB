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
| `stage2_flat_interface.mph` (Chase) | 2 | reference | STORED_DERIVED | COMSOL 6.3.0.420, sha256 `fe4cb245…3de43`. Model tree read 2026-11-13: hi = 4.2e5; trimmed pillar; heater bnd 3, meniscus bnd 11, T_base bnds 6 + 10; hmax sweep 2/1/0.5/0.25 µm solved; exported column 1 = automatic mesh, column 2 = probably 0.25 µm; `dp_dry` typo (`1/d`) (`data/reported/stage2_chase_identified_from_model_file.md`) | evaluation of the four sweep solutions on the parametric dataset; A_base; interior-plane flux; `.java` copy |
| `model 1.mph` (Chase, D:\NANOFAB\COMSOL) | sketch | — | not a result | unsolved geometry sketch with parameters r = 10 µm, h = 50 µm, pitch 50 µm. Note h = 50, not the stated 40 µm design height | none needed; confirm the intended proposed height |

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
