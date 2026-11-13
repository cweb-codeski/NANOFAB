# Reported results (transcribed, not raw exports)

Values that a team member **reported** (in a message, slide or notebook) but whose raw COMSOL
export and `.mph` file are not in the repository. They are kept so the comparison is reproducible
and traceable. They are not evidence of a particular model setup.

Rules:

- Copy the numbers exactly as reported. Never edit a value; add a new file if a correction arrives.
- The header says who reported it, when, what is known about the setup, and what is **not** known.
- Status is always `REPORTED_UNVERIFIED` until the raw export, the `.mph` file, and the selection
  and BC record are available. Then the run moves to `data/comsol_exports/` and
  `comsol/unit_cell/README.md` like any other run.
- Never average reported columns, assign them to meshes, or derive an energy balance unless the
  report itself identifies them.

| File | Reported by | Stage / case | Status |
| --- | --- | --- | --- |
| `stage1_student_reference_d10_l30_h25.csv` | undergraduate (independent Stage 1 build) | Stage 1, reference | REPORTED_UNVERIFIED |
| `stage2_chase_reported_columns.csv` | Chase (text exports of `stage2_flat_interface.mph`, COMSOL 6.3.0.420) | Stage 2, reference | REPORTED_UNVERIFIED; columns lack mesh and selection identifiers in the export itself |
| `stage2_chase_identified_from_model_file.md` | model-file inspection (2026-11-13) | Stage 2, reference | settings and column-to-solution mapping from the `.mph` tree and action history; column 1 = automatic mesh, column 2 = probably hmax 0.25 µm |

Comparison script: `python analysis/stage1_crosscheck.py`.
