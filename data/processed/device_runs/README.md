# Device-model run records

Written by `analysis/run_passive_device.py` and `analysis/device_benchmarks.py`. Each invocation
creates a new `<output>/<run_id>/` (UTC timestamp + commit prefix + mode; never overwritten) with
`run_manifest.json` (command line, git commit and dirty state, environment and package versions,
SI parameter snapshot with provenance, table ids and hashes, qualification reports, settings,
convergence, balances, regime diagnostics, boundary conditions, status and exit code) plus CSV
profiles and fields. Plots go to `plots/device/<run_id>/` (git-ignored).

- `verification/`: benchmark suites (VERIFICATION ONLY inputs).
- `vaartstra_preflight/`, `vaartstra_prediction/`: literature case; currently blocked with
  MISSING_VALIDATED_INPUTS (`missing_inputs.json` names what is missing).
- `verification_fixture_run/`: full pipeline on manufactured tables (VERIFICATION ONLY).
- `import_demo/`: import CLI demonstration on the manufactured raw-export fixture.
