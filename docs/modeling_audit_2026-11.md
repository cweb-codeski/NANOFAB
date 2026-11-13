# Modeling Workflow Audit, 2026-11-13

Scope: the staged COMSOL workflow, its documentation, the stored and reported results, and the
Python support code, checked against held primary papers (PDF text read on 2026-11-13), official
COMSOL guidance (EO blog, KB 1261), dimensional consistency and the repository's own data.
Changes made in response are in `docs/comsol_unit_cell_guide.md` (canonical plan),
`docs/validation_and_experiments.md`, `docs/stage2_undergrad_assignment.md` and the files listed
in §4.

Environment limits during the audit: COMSOL was not installed, so no model was built, opened or
solved; licensed modules are unknown. The local COMSOL folder (`D:\NANOFAB`) was not accessible.
The CAD files `Wick_V2.cat` / `Coupon_V2.cat` were not available, and their stated values are
taken as given.

## 1. Prioritized findings

P1 = could materially change a reported number or conclusion. P2 = could mislead a reader or a
later stage. P3 = hygiene.

| # | P | Finding | Evidence | Consequence | Correction |
| --- | --- | --- | --- | --- | --- |
| 1 | P1 | Stage 2 geometry step used **Union of the block with the full cylinder**, leaving 3/4 of the pillar outside the cell. The "3 domains" check does not catch it. | old guide §2 steps 2, 4, 6 | insulated pillar extension and extra z = 0 faces bias T_base; Chase's Stage 2 model may contain it (unknown) | trim by Intersection; volume and A_base checks (guide §4, §6); inspect `stage2_flat_interface.mph` |
| 2 | P1 | Stage 2 evidence over-read: a gitignored assembled export labels rows "coarser/finer" by **inference** and picks the Q_evap column that balances to 4e-12; `data/processed/h_evap_stage2_flat_d10_l30_h25.csv` derives h = 62 kW/m²K from it. Q_in is the `q0*(l/2)^2` echo. | export header text; reported columns lack identifiers | false convergence and energy-balance evidence; 62 kW/m²K sits next to V2019's "h_e > 60" and invites a false validation claim (flat interface, linearized closure: not comparable) | transcribed as REPORTED_UNVERIFIED (`data/reported/`); classified unqualified in `config/data_inventory.yaml`; post-processor now refuses unlabeled rows and checks A_base; Q_in must be a heater-flux integral |
| 3 | P1 | **Interface shear condition is uncertain.** V2019 uses a shear-free meniscus (Sec. 2.2). Fleming et al. 2020 measured curvature in pure vapor that agreed much better with a **no-slip** interface (Sec. 3.1, Figs. 4-5). | held PDFs | flat-top K falls ~31% (ref) / ~28% (proposed) under no slip; q''_dry falls in proportion; EO flux also changes. The reference flat q''_dry brackets V2019's 82 W/cm² between 101 (shear-free) and 69 (no-slip) W/cm², so the "−19% gap is curvature" story is not unique | bracket runs in 1Q, 2A (R4), 3H; Exp-1/Exp-2 discriminate; never choose one to match a target |
| 4 | P1 | **Linear Darcy unverified at device gradients, especially for the proposed cell.** Hand estimates at the capillary limit: reference Re_K ≈ 1.7; proposed Re_K ≈ 2-8, Re_h ≈ 10-40 (illustrative L = 10-2.5 mm). The quarter-cell pressure planes are exact only for Stokes flow. | `case_hand_checks.py` Sec. C; symmetry argument (guide §4) | K and q''_dry of the proposed design could be overestimated near the feed | Stage 1Q: periodic half cell, Navier-Stokes vs Stokes gradient sweep, predeclared 1% rule (contract) |
| 5 | P1 | **Evaporation closure and accommodation.** The Meland-Ytrehus treatment V2019 cites (Sec. 2.2, ref. 28) is not held; our correction form is from memory. a = 0.052 was fitted (V2019 Sec. 3.4) to Adera A1, measured at ~24 °C / 3 kPa after subtracting a 62% heat loss (adera2016 Sec. 4). | held PDFs | no Stage 4 output can be called a V2019 replication; `a` from 0.052 to 1 moves the Stage 2 reference superheat by ~1.3 K of ~8 K. The fitted `a` absorbs the loss correction and model-form errors | closure gate and named closures (guide §9); sweep `a`; calibrate only on own pure-vapor data, validate on a different condition |
| 6 | P1 | **Vapor environment undefined.** Unit-cell h_e assumes pure saturated vapor. In open air, evaporation is diffusion-limited. | V2019 Sec. 3.4 (air data used only for dryout); docs/assumptions A6 | heated tests in air cannot validate h_e or superheat predictions | `vapor_environment` UNRESOLVED in the case file; Exp-4 states what is comparable in each environment |
| 7 | P1 | **Dryout logic.** Reduced-order models with prescribed uniform evaporation, or an interface that stays wet regardless of supply, cannot establish dryout. Other limits can come first: Adera et al. report boiling before dryout. EO adds flooding and electrochemical limits. | adera2016 Sec. 4; `docs/electroosmosis_model.md` | a "dryout gain" may be a limit that is never reached | Stage 6 couples supply and evaporation with conservation checks and reports which limit is reached first |
| 8 | P1 | **Reference → proposed contamination.** Post-processors fell back silently to the reference geometry when an export lacked d/l/hp columns. All stored K, F, EO and poster numbers are for the reference cell. | `analysis/_common.py` (before fix) | a proposed-cell export would be normalized with the l = 30 µm footprint | case file + loader with UNRESOLVED refusal; `--case` flags; warning on fallback; tests |
| 9 | P1 | **Downstream (cathode) reservoir hydraulic BC is undefined.** | CAD brief (separate cathode compartment, gas port) | sealed: EO raises end pressure (flooding risk); vented: EO drives net throughflow (subcooled feed adds sensible heat); gas-filled: pressure drifts | must be chosen per run (guide §11); `reservoir_geometry` UNRESOLVED |
| 10 | P1 | **Array geometry unresolved.** 205 × 155 at 50 µm pitch spans 10.22 × 7.72 mm (outer edges). Under either reading of the (−10, −7.5) mm start reference, the array sits ~4.9 mm and ~3.6 mm off the coupon center. | `start_reference_interpretations` | feed length, heater placement and any device model depend on the true layout; a 20 × 15 mm wick must not be substituted | recorded as UNRESOLVED; resolve in the native CAD application (scale, units, `Coupon_V2.cat` reference) |
| 11 | P2 | **Field treatment and electrical accounting.** The uniform imposed field gives ~0.5-0.6× the solved-field EO flux in the 2D limit (done in Python); 3D pending. The field in the wick is not supply voltage / length; DC current needs Faradaic reactions. | `docs/eo_field_treatment_test.md`; `MODEL_STATUS.md` | EO magnitudes and power can be off by ~2× or more; bulk Joule power is not the electrical input | Stage 2A detailed procedure with signed current integrals, P_joule = I·ΔV check, and E0 inferred from measured current |
| 12 | P2 | **EO stall pressure is small compared with reservoir head effects.** Illustratively, tens of Pa over 10 mm at 1e4 V/m and ζ = −50 mV, i.e. a few mm of water. | `case_hand_checks.py` Sec. E | unheated EO tests are dominated by head drift unless controlled | Exp-3 controls; matrix |
| 13 | P2 | **Stage 3 endpoint.** V2019 Eq. 3 is the mean force balance (exact identity of the Young-Laplace equation). The local edge angle reaches θ_rec earlier where the slope is steepest. Floor contact is not expected (sag ~10-12% of hp for both cells). | derivation in guide §8; hand estimate | dryout Pc may be overestimated if local depinning governs | report mean and local criteria; endpoint is a recorded choice |
| 14 | P2 | **Measured vs model temperature planes.** Experiments measure heater/backside temperatures through the substrate with losses; the model reports the pillar-base plane. | adera2016 (650 µm substrate, 1D correction, 62% loss) | direct comparison would mix substrate, spreading and loss errors into h_e | Exp-4 inference chain; Stage 7 trigger |
| 15 | P2 | **Thin-EDL validity** holds in the gaps for all screened chemistries (λ_D/half-gap ≤ 0.1, finite-EDL factor ≥ 0.91), but the Dukhin lower bound reaches ~0.1 for pure water, and the contact-line film is locally thinner than λ_D. | `case_hand_checks.py` Sec. E | ~10% EO error from surface conduction in very dilute water; local film region unmodeled | prefer a defined dilute electrolyte; resolve EDL only if the chosen chemistry fails the screen |
| 16 | P2 | **EO and the conduction approximation.** Joule ≤ 1e-5 of q''; vertical Pe ~1e-2; net flow is fixed by evaporation unless the cathode compartment allows throughflow. | `MODEL_STATUS.md`; V2019 App. A.1 | EO thermal coupling at the unit cell is likely negligible | one 4E check instead of an EO thermal stage |
| 17 | P2 | **Silicon conductivity** 148 W/m/K is a 300 K value from memory; the 373 K value is lower (also from memory, ~100-120). The sweep lower bound of 120 may be too high. Pillar conduction is a large share of the reference R''_wick. | `parameters/materials.yaml` | Stage 2/4 superheat biased low | keep ASSUMED; obtain a held source; widen the sweep when sourced (not changed silently) |
| 18 | P2 | **Three Python device-level models** (`device/`, `whole_wick/` with electrowetting and an ASSUMED h_e placeholder, `electroosmosis/wick1d`) with different assumptions. | repo | risk of quoting a reduced-order number as a prediction | roles fixed as cross-checks of COMSOL Stage 6; no extension planned |
| 19 | P2 | **Convergence threshold inconsistency:** modeling context said 1-5%; the contract requires 1e-3. Chase's stored Stage 1 finest pair is 1.02e-3; the student's is 5.05e-4. | docs | "converged" meant different things | preliminary (≤ 1e-2, labeled) vs qualified (≤ 1e-3) in the guide; a Richardson-based alternative is proposed, not adopted |
| 20 | P3 | Old Stage 2 mesh table listed a "Corner Refinement" node (a 2D feature). | old guide §4 | build confusion | replaced by edge/boundary Size nodes |
| 21 | P3 | `analysis/permeability_from_comsol.py` rewrites `data/processed/permeability_<stem>.csv` on every run. | code | a rerun with different flags could overwrite a stored result | noted; run with a different export stem or commit first |

## 2. Stage 1 and Stage 2 evidence, as it stands

- **Stage 1 reference.** Chase (`permeability1.mph`, stored): K = 2.96979e-11 m² at 0.5 µm;
  changes 1.02%, 0.20%, 0.10% (coarser denominator). Student (reported): K = 2.9711e-11 m² at
  0.25 µm; changes 1.068%, 0.297%, 0.0505%. Matching meshes differ by −0.047%, −0.0007% and
  −0.006%; finest-available values differ by +0.044%. Both are 1.3% below the analytical flat-top
  estimate. This is reproducibility evidence for one flat-top, low-gradient case. It is not
  qualification (no gradient sweep, no interface bracket, `.mph` not tracked).
- **Stage 2 reference.** Two Q_in values (1.125e-4 W, the input echo) and two T_base values
  (381.19639 K, 381.21757 K) without identifiers. Nothing beyond "a Stage 2 model exists and
  produced superheats of ~8.05-8.07 K under unknown mesh/selection" can be concluded.

## 3. Source and assumption register (consequential items only)

| Item | Current treatment | Provenance | If wrong | Resolution |
| --- | --- | --- | --- | --- |
| Interface shear condition | shear-free baseline, no-slip bracket | V2019 Sec. 2.2 vs fleming2020 Sec. 3.1 | K and q''_dry ±30%; EO flux changes | Exp-1 + Exp-2 on our surface |
| Accommodation treatment (Meland-Ytrehus) | reduced form from memory, labeled ASSUMED | meland_ytrehus2003 not held | blocks V2019 replication claims | obtain Phys. Fluids 15 (2003) 1348 |
| Accommodation coefficient | 0.052, swept 0.03-1 | fitted (V2019 Sec. 3.4) | superheat ±~15% (reference, flat) | pure-vapor calibration on own wick, then validation |
| θ_rec, proposed surface | 15° ASSUMED, swept 5-40° | V2019 value transferred | Pc_max ∝ cos θ_rec | co-processed coupon measurement |
| ζ (pillar, floor), electrolyte, σ | parameters only; no device value | REQUIRES_SOURCE (`eo_inputs.yaml`) | all EO magnitudes linear in ζ | streaming potential/current |
| Silicon k at operating T | 148 W/m/K ASSUMED | memory | superheat bias | held property table |
| Proposed geometry | 20/50/40 µm USER_PROVIDED | design input | every proposed-cell result | SEM/profilometry of the fabricated wick |
| Array layout, CAD scale/units | UNRESOLVED | CAD not inspected | feed length, heater placement | native CAD check |
| Feed length, substrate, heater, reservoirs, vapor environment (proposed) | UNRESOLVED | — | Stage 6 cannot run | design decisions |
| Cathode compartment hydraulic BC | UNRESOLVED | — | EO sign of effect (throughflow vs pressure rise) | design decision |
| Linear Darcy at device gradients | unverified | Re estimates | K overestimated near feed | Stage 1Q |
| Water permittivity at 25 °C | 78.4 ASSUMED (screening only) | memory | Debye-length screening only | IAPWS release |

## 4. Change summary (this revision)

Added:
- `parameters/unit_cell/cases.yaml`, `src/nanofab/geometry_cases.py`: separate reference and
  proposed cases; unresolved inputs refuse to return values; array extent and start-reference
  readings.
- `analysis/case_hand_checks.py`: hand-calculation equivalents for Stages 1-2A, both cases.
- `analysis/stage1_crosscheck.py`, `data/reported/` (README, student Stage 1, Chase Stage 2):
  reported results preserved verbatim with status.
- `nanofab.permeability.brinkman_no_slip_top_factor` and `pillar_wick_permeability(top=...)`.
- `tests/test_geometry_cases.py` (18 tests).
- `docs/validation_and_experiments.md`, `docs/stage2_undergrad_assignment.md`, this audit.
- Sources: `fleming2020`, `comsol_blog_eof_edl`, `comsol_kb_1261`; `adera2016` now held with
  checked details.

Changed:
- `docs/comsol_unit_cell_guide.md` rewritten as the canonical staged workflow (stage cards,
  Stage 1Q, corrected Stage 2, detailed Stage 2A, Stages 3-7 specifications).
- `analysis/_common.py`, `permeability_from_comsol.py`, `heat_transfer_from_comsol.py`: `--case`,
  fallback warning, mesh-identifier guard, A_base check.
- `docs/modeling_context.md`, `MODEL_STATUS.md`, `comsol/unit_cell/README.md`,
  `config/data_inventory.yaml`, `docs/SOURCE_AUTHORITY.md`, `docs/cell_table_contract.md`,
  `docs/electroosmosis_model.md`, `docs/eo_field_treatment_test.md`, `README.md`: pointers,
  status and classifications.

Not changed: raw exports, stored results, existing equations and parameter values, the Python
reduced-order models.

## 5. Remaining blockers

1. COMSOL access and licensed-module list (every COMSOL step here is unexecuted).
2. `stage2_flat_interface.mph`: hi value, pillar trimming, T_base selection, mesh identity of each
   reported column.
3. Meland & Ytrehus (2003) for any V2019 replication claim.
4. Proposed-device definitions: CAD scale/units and array placement, feed arrangement and length,
   substrate thickness, heater footprint, reservoir levels, cathode compartment BC, vapor
   environment, electrolyte.
5. Surface data on our oxide: θ_adv, θ_rec, ζ (streaming), aging behavior.

## 6. Next three concrete actions

1. **Chase:** open `stage2_flat_interface.mph`, run the geometry Measure check (pillar volume
   490.9 µm³) and record hi, the T_base selection with its area, and which mesh/study produced
   each reported column. Save a `.java` copy and register the model in `comsol/unit_cell/README.md`.
2. **Undergraduate:** do `docs/stage2_undergrad_assignment.md` on the reference cell (hand in
   model, selections, Q_in, Q_evap_bc, Q_evap_flux, T_base, A_base, T_lv, mesh study), then
   duplicate for the proposed cell.
3. **Either, in parallel:** Stage 2A runs R0-R2 on the reference cell (guide §7), with Ũ·F = 1
   as the pass/fail check, followed by Stage 1 + 1Q on the proposed cell (inertia sweep and
   interface bracket).
