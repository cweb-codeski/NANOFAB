# Modeling Workflow Audit, 2026-11-13

Scope: the staged COMSOL workflow, its documentation, the stored and reported results, and the
Python support code, checked against held primary papers (PDF text read on 2026-11-13), official
COMSOL guidance (EO blog, KB 1261), dimensional consistency and the repository's own data.
Changes made in response are in `docs/comsol_unit_cell_guide.md` (canonical plan),
`docs/validation_and_experiments.md`, `docs/stage2_undergrad_assignment.md` and the files listed
in §4.

Environment limits during the audit: COMSOL was not installed, so no model was built, opened or
solved; licensed modules are unknown. In revision A the local COMSOL folder (`D:\NANOFAB`) was not
accessible. In revision B (§7) it was read: the raw Stage 2 tables, plus the model tree and
action history inside `stage2_flat_interface.mph` (a ZIP container), read-only. The CAD files
`Wick_V2.cat` / `Coupon_V2.cat` were not available, and their stated values are taken as given.

## 0. Revision B (2026-11-13, later the same day): source-level re-audit

Every consequential claim was re-checked against the held PDFs (page numbers are journal pages; the
thesis uses its printed page numbers) and against Chase's Stage 2 model file. Labels:
**M** = measured by the authors, **I** = inferred by the authors from measurements through a model,
**A** = assumed/chosen by the authors, **C** = calculated by us. "Affects" says whether the finding
changes the **reference replication** (REF), the **proposed-device model** (PROP), or only an
**uncertainty/sensitivity study** (SENS). Agreement between papers, or a successful fit, is never
counted as experimental validation.

### 0.1 No-slip vs shear-free liquid-vapor interface

| # | Claim (earlier wording → verdict) | Evidence | Type | Affects |
| --- | --- | --- | --- | --- |
| B1 | V2019 uses a shear-free interface for K → **stands** | "importing the corresponding interface shape and setting it as a symmetry boundary, which enforces no flux and no shear stress at the interface" (V2019 Sec. 2.2, p. 172); thesis Fig. 4c caption (p. 20) | A | REF: keep shear-free for replication |
| B2 | "Fleming measured that the interface is no-slip" → **corrected** | They measured curvature by thin-film interference (Sec. 2.2, p. 43268). No-slip is the better-fitting of two model bounds: "When using these two bounds for permeability, we found much better agreement with the no-slip condition when comparing experiments vs a semianalytical Darcian model with curvature-dependent permeability" (p. 43270). Model: Eqs. 1-3 (p. 43268), COMSOL creeping flow with "a periodic pressure boundary condition" and "The impact of evaporation on the flow dynamics is neglected" (p. 43269) | curvature M; interface condition I | PROP and SENS: bracket, not a replacement |
| B3 | Fleming's conditions | room temperature, degassed DI water, "~20 mbar" (p. 43268); Piranha-treated just before testing (p. 43267); vertical wick with gravity (Eqs. 1-2); d/pitch/height 60/90/50 and 80/120/50 µm, solid fraction 0.35 (Table 1, p. 43267) | M | SENS: geometry is far from ours (φ 0.087 / 0.126); not a matched benchmark |
| B4 | Fleming's proposed mechanism | surface-tension gradients from "surfactant contamination or mobile surface charge" and heating (p. 43269), and possible depinning near the top and near dryout (pp. 43269-43270) | A (hypotheses) | experiment design: surface cleanliness and pinning must be recorded |
| B5 | "Flat K falls 31% / 28% under no-slip" → **stands as our number only** | Brinkman closed-channel factor, `nanofab.permeability.brinkman_no_slip_top_factor` | C | SENS |
| B6 | "V2019's 82 W/cm² lies between our flat 101 (shear-free) and 69 (no-slip), so the −19% gap is not uniquely curvature" → **withdrawn** | V2019 used a shear-free interface (B1), so its 82 W/cm² must be reproduced within its own model: curvature-dependent K, non-uniform evaporation ("first calculated for uniform evaporation ... iterated once more—now accounting for non-uniform evaporative flux", Sec. 3.3, p. 173) and any error in our analytical K. The no-slip bracket is about real wicks, not about V2019 | C | REF: interpretation corrected; no number changes |

### 0.2 Origin and applicability of a = 0.052 and the 62% heat loss

| # | Claim → verdict | Evidence | Type | Affects |
| --- | --- | --- | --- | --- |
| B7 | a = 0.052 was fitted to Adera "A1" → **stands, sharpened** | "We first determined the accommodation coefficient a = 0.052 by fitting to one set of data from [12] (A1) and then verified that our model matched well to their other data (A2)" (V2019 Sec. 3.4, p. 174). A1 = d/l/h 5/12/82 µm, A2 = 12/20/90 µm (Fig. 6 caption, p. 175) = Adera devices 1 and 5 (Adera Table 2, p. 284). Thesis: fitted "to the last point on (A1)" because the low-flux points "are offset (the intercept is not zero)" (p. 32) | fit = C by V2019 | REF: use as V2019's input, labeled fitted |
| B8 | a is a property → **no; the authors say so** | "This accommodation coefficient is not intended to represent the theoretical evaporation/condensation coefficients because it absorbs the non-idealities of the experiments from which it was calculated" (thesis p. 22) | A | PROP: a is not transferable; sweep and calibrate on own data |
| B9 | The A2 agreement validates the model → **no** | A2 is another device from the same apparatus, chamber and loss correction. The thesis itself calls the mean-superheat comparison "not an entirely satisfactory validation" (p. 33) | — | SENS |
| B10 | Same data, different closure, different a (**new**) | Adera's own model uses "the Schrage equation [36] based on the experimentally obtained average accommodation coefficient (0.08)" (Eq. 16, p. 289) | I | SENS: a depends on the closure and the rest of the model |
| B11 | Measured receding angle on A1 (**new**) | interferometry vs Surface Evolver gave 65° and 75°, "The average receding contact angle 70° was used", attributed to C4F8 residue from DRIE (Sec. 5.1.3, Fig. 9, p. 288). V2019 uses 15° "for water on silicon dioxide [29]" (p. 172) and does not state the angle used for its A1/A2 runs | I (Adera) / A (V2019) | SENS: the fit of a is conditional on unstated V2019 inputs (θ_rec, substrate, property temperature, 2D feed). Also a warning for our own DRIE wicks: C4F8 residue can raise θ_rec |
| B12 | "62% heat loss" → **stands, scoped** | "The combined heat loss via sensible, convection, and ohmic was 62% of the input heating power. This heat loss was deducted from the input power in all the heat fluxes reported in this study" (Adera Sec. 4, p. 286); Fig. 6 shows it for device 1 (= A1). Components: sensible from measured flow and temperatures (M); convection "estimated conservatively for a 10 W/m² K heat transfer coefficient" (A); ohmic from measured resistance (C) | mixed | SENS: a was fitted to loss-corrected fluxes; V2019 never mentions the correction. Whether the 62% fraction was applied unchanged to every device is not stated |
| B13 | Superheat is measured → **inferred** | ΔT from backside RTDs by 1D conduction through 1 µm SiO₂ + 650 µm Si (Adera Eq. 2, p. 285); chamber evacuated "to remove all noncondensable gases" (p. 282), held at "24 C and corresponding saturation pressure 3 kPa" (p. 284) | I | SENS |
| B14 | "Adera reports boiling before dryout" → **wrong, corrected** | those entries are in Adera Table 1 (p. 282), a summary of other groups' studies. Adera's own experiments are "in the absence of nucleate boiling" (Sec. 4, p. 285) | — | none; text corrected in guide §11, sources and §1 item 7 |
| B15 | a fitted at 24 °C applied at 100 °C → **stands as an assumption** | the V2019 parametric and device studies are at 100 °C (Sec. 3.1-3.3, pp. 172-173) with the 3 kPa-fitted a | A | REF: replicate as V2019 did; PROP and SENS: not transferable |

### 0.3 Meland-Ytrehus correction and our evaporation law

| # | Claim → verdict | Evidence | Type | Affects |
| --- | --- | --- | --- | --- |
| B16 | V2019 does not print the correction → **stands** | "was accounted for using the method outlined by Meland and Ytrehus [28]" (V2019 p. 172); thesis: Eqs. 6-7 "correspond to the 'total accommodation' case and must be corrected for fluids such as water" (p. 21) and cites the method (p. 22) without a formula. Ref. [28] = Phys. Fluids 15 (2003) 1348 (p. 177); still not held | — | REF: source gap remains for the general form |
| B17 | "Our correction form is unverified" → **functionally checked at a = 0.052** | thesis Fig. 7 (p. 23) plots q''_e vs ΔT_int "for a = 0.052". Digitized (`analysis/digitize_vaartstra_thesis_fig7.py` → `data/literature/vaartstra2019_thesis_fig7_digitized.csv`). Our `sone_bte_heat_flux` (C4 → C4 + 2√π(1−a)/a, ρ_k from the liquid-side saturation state, T_sat = 100 °C) matches within 0.5% for ΔT = 3-9.5 K and within 5 W/cm² (1% of axis) everywhere; a least-squares fit of a returns 0.052. C4/a overshoots by ~60%; using the far-field vapor density undershoots by 28% at 9.5 K (`tests/test_evaporation_thesis_fig7.py`) | digitization = PROJECT_DERIVED; comparison = C | REF: our law reproduces V2019's own curve at a = 0.052. PROP/SENS: the form for other a is still unverified |
| B18 | What controls the closure | at a = 0.052 the accommodation term 2√π(1−a)/a = 64.6 dominates C4 = 2.13, so the law is close to Hertz-Knudsen-Schrage scaling; Schrage at the same a is 2.7% lower (432.9 vs 445.0 W/cm² at 9.5 K). Dropping V2019 Eq. 6 (temperature jump) changes the mass flux only at second order through ρ_k | C | SENS: the value of a matters far more than BTE vs Schrage |
| B19 | Fig. 7 conditions | the figure does not print T_sat. At 3 kPa the same law gives ~20× smaller fluxes, so 100 °C is the only consistent reading | inference (ours) | REF |

Contract consequence: `docs/cell_table_contract.md` still requires `accommodation_treatment_verified:
true` for a V2019 replication. The functional match is recorded as evidence. Proposed revision,
**not adopted**: accept `accommodation_treatment_verified: functional_match_thesis_fig7`
for a = 0.052 at T_sat = 100 °C only. Effect: allows a V2019 reference replication now. It does not
allow any other a.

### 0.4 Geometry, averaging and normalization needed to reproduce V2019

| # | Item | Evidence | Type | Affects |
| --- | --- | --- | --- | --- |
| B20 | Thermal/meniscus domain is the 1/8 triangle | "Only the region enclosed by the dashed lines is treated as the domain" (V2019 Fig. 2a, p. 172); triangle drawn in thesis Fig. 4a (p. 20). Our quarter cell is exactly equivalent for these symmetric problems | A | REF: no change |
| B21 | Flow-study domain and BCs are not specified | "A pressure difference was imposed across the cell and the mass flux (Γ) was reported" with laminar Navier-Stokes, "Re_h ~ 10¹" (p. 172); Fig. 2c shows p₁, p₂ across one cell | A (unstated details) | REF: our Stokes quarter cell equals their setup only if inertia is negligible at their Δp (unknown); Stage 1Q covers this |
| B22 | T_e definition | "average temperature at the surface (red line)" (V2019 Fig. 2b, p. 172). In thesis Fig. 4b (p. 20) the red line spans the full cell width at the pillar-base plane. Area weighting is not stated (our PROJECT_DERIVED reading) | A | REF |
| B23 | h_e is per footprint, heat input 1 µm below the pillar base | Eq. 4, "a uniform heat flux is supplied to the bottom of the cell, 1 µm below the base of the pillars" (p. 171) | A | REF |
| B24 | K is superficial (full wick height) | Eq. 8 does not state the area basis. Device Eq. 9 uses Q_m = q''_e/(h_lv h) over the full height h (p. 172), which implies a superficial basis | inference (ours) | REF: our bulk convention is consistent |
| B25 | Device BCs | "Reference pressure is set at the inlet while x = L is a symmetry boundary" and thermal symmetry at x = 0 and L, "reasonable so long as liquid is supplied at T_sat" (p. 172); dryout iteration as in B6 | A | REF |
| B26 | Not stated by V2019 | silicon conductivity, water property source, flow-study BCs, θ_rec and substrate used for Fig. 6a | — | REF: sweep k_Si; record our choices |

### 0.5 Electric-field and EO assumptions in Pujahari et al. 2022

| # | Claim → verdict | Evidence | Type | Affects |
| --- | --- | --- | --- | --- |
| B27 | "Uniform imposed field" → **stands** | "As the double derivative of external potential is zero ... it is assumed to be constant along the wick direction. Therefore, a constant electric field is considered to be acting along each unit cell to calculate the volumetric force" (Sec. 2.1.2, Eq. 8, p. 8446); BC "E = dφ/dx at x = 0 and l plane" (p. 8446) | A | PROP: our solved-field baseline remains the physical choice for insulating pillars |
| B28 | Their EO mechanism is HS slip → **corrected** | "Unlike many studies, here, the effect of electro-osmosis is not just added ... by using the Helmholtz−Smoluchowski slip velocity boundary condition, but the EDL distribution of ions is calculated and then the volumetric Coulombic force is added to the Navier−Stokes equation" (pp. 8444-8445). Poisson-Boltzmann (Eqs. 5-6) with "Ψ = ζ at the bottom wall (z = 0) and at the pillar walls", zero normal derivative at the meniscus; "u, v, w = 0 at the bottom wall and at the pillar walls" (p. 8446) | A | SENS: our uniform-field (Case A) comparison is the thin-EDL equivalent. For 1e-6 M (λ_D 0.3 µm) and d = 5 µm at small pitch, λ_D/half-gap reaches ~0.1, so that equivalence carries a finite-EDL error of order 10% |
| B29 | Consistent with COMSOL guidance | resolved EDL + body force + no-slip is the blog's alternative to slip when λ_D is not small ("use the traditional no-slip velocity wall boundary condition and resolve the flow in the EDL"); no double counting | — | none |
| B30 | ζ values | "10⁻⁶ and 10⁻⁴ M ... in agreement with the experimental work of Mala et al., resulting in zeta potentials of −200 and −100 mV" (Sec. 2.1.1, p. 8446); ε(T) = 305.7 exp(−T/219) (p. 8446) | A (from another paper) | PROP: no ζ for our surface |
| B31 | Electrolysis statement | "for the considered range of electric field and wick length, the voltage values will be very less compared to the voltage that causes electrolysis" (pp. 8446-8447); "the applied voltage or electric field should not be high enough to result in electrolysis or water splitting" (p. 8453). E × 5 mm = 25-500 V is our calculation | A / C | PROP: electrode chemistry is unaddressed in the paper |
| B32 | Power claim | "in the order of milliwatts (Supplementary Information)" (p. 8452); Joule heat "a maximum of 4.5 × 10⁻⁸" W (p. 8453) | C by authors; SI not held | SENS |
| B33 | Their thermal closure (**new**) | Schrage (Eq. 9) is "valid because the values of the accommodation coefficient (0.052)..." (p. 8447): V2019's BTE-fitted a inside a different closure (~3% at that a, B18) | A | SENS |
| B34 | Their temperature "validation" is independent → **no (new)** | Fig. 3b uses d/l/h = 5/12/82 µm (p. 8449), Adera device 1, the data set V2019 used to fit a. Dryout check (Fig. 3a) uses Zhu's air-ambient B1 device | — | SENS |
| B35 | Velocity definition | U_i is the integral of u over the liquid midplane area A (Sec. 2.1.2, p. 8447), with mass M_i = ρ A_i U_i (Eq. 13). It is an intrinsic (liquid-area) velocity, not superficial | A | SENS: our Fig. 6 reproduction used dryout fluxes, so no number changes |

### 0.6 Chase's Stage 2 model file (inspected read-only)

Details and hashes: `data/reported/stage2_chase_identified_from_model_file.md`.

| # | Finding | Evidence | Affects |
| --- | --- | --- | --- |
| B36 | hi = 4.2e5 W/(m²K), not 4.17e5 | parameter `hi`; on column 2, ∫ht.ntflux / ∫(T − T_v) dA = 4.200000e5 | REF: assignment updated |
| B37 | Pillar **is** trimmed | geometry: Cylinder ∩ block, liquid block minus the trimmed pillar, Union with substrate (history actions 139-165) | REF: §1 item 1 resolved for this model; the guide check stays for new builds |
| B38 | Columns are not one mesh sweep | column 1 = automatic physics-controlled mesh (actions 204-330); column 2 = "Study 2/Solution 1" after the hmax sweep, probably 0.25 µm (actions 351-399). The 2, 1 and 0.5 µm solutions are stored but were never evaluated | REF: convergence evidence must be re-evaluated from the stored sweep |
| B39 | Column-2 "energy balance" is an identity | ht.ntflux on the Robin boundary equals hi∫(T − T_v) dA to 1e-6 | REF |
| B40 | Column-1 inconsistency | ∫ht.ntflux = +9.3e-4 vs Q_in, but hi∫(T − T_v) dA from the same table = −7.4e-4; unexplained without COMSOL | REF |
| B41 | `dp_dry` typo | `(1/d)` instead of `(l/d)`; evaluates to 1.8e-6 Pa | Stage 3 reuse would fail |
| B42 | T_base selection is two faces (boundaries 6 and 10) | `av2` selection | REF: area check still needed (expect 2.25e-10 m²) |
| B43 | Heater selection is the full bottom face | on column 2, ∫ht.ntflux over boundary 11 equals q0(l/2)² to 4e-12 with all other faces insulated, so boundary 3 has area (l/2)² | REF |
| B44 | "Corner Refinement is a 2D-only node" → **wrong** | the 3D mesh sequence contains `Corner Refinement 1` | guide text corrected |

## 1. Prioritized findings

P1 = could materially change a reported number or conclusion. P2 = could mislead a reader or a
later stage. P3 = hygiene.

| # | P | Finding | Evidence | Consequence | Correction |
| --- | --- | --- | --- | --- | --- |
| 1 | P1 | Stage 2 geometry step used **Union of the block with the full cylinder**, leaving 3/4 of the pillar outside the cell. The "3 domains" check does not catch it. | old guide §2 steps 2, 4, 6 | insulated pillar extension and extra z = 0 faces bias T_base | trim by Intersection; volume and A_base checks (guide §4, §6). **Rev. B:** Chase's model is **not** affected; its geometry sequence trims the pillar (Cylinder ∩ block, then Difference, then Union; `data/reported/stage2_chase_identified_from_model_file.md`). The old-guide risk remains for new builds |
| 2 | P1 | Stage 2 evidence over-read: a gitignored assembled export labels rows "coarser/finer" by **inference** and picks the Q_evap column that balances to 4e-12; `data/processed/h_evap_stage2_flat_d10_l30_h25.csv` derives h = 62 kW/m²K from it. Q_in is the `q0*(l/2)^2` echo. | export header text | false convergence and energy-balance evidence; 62 kW/m²K invites a false comparison with V2019's "h_e > 60" | as before, plus **rev. B** from the model file: column 1 is the **automatic-mesh** solution, column 2 is "Study 2/Solution 1" after the hmax sweep (probably 0.25 µm). They are not a convergence pair. The 4e-12 balance is `ht.ntflux` returning the Robin BC flux (Q_evap/∫(T−T_v)dA = 4.200000e5 = hi). hi = 4.2e5 exactly. The 2/1/0.5 µm solutions exist but were never evaluated |
| 3 | P1 | **Interface shear condition is uncertain.** V2019 uses a shear-free meniscus (Sec. 2.2, p. 172). Fleming et al. 2020 measured local curvature (room temperature, ~20 mbar vapor, φ = 0.35 pillars) and **inferred**, through a Darcian model with COMSOL K(κ), that a no-slip interface fits much better (abstract; Sec. 3.1, Figs. 4-5, pp. 43269-43270). | held PDFs | flat-top K falls ~31% (ref) / ~28% (proposed) under no slip (our Brinkman estimate for our cells, not a Fleming number); q''_dry falls in proportion; EO flux also changes | bracket in 1Q, 2A (R4), 3H; Exp-1/Exp-2 discriminate. **Rev. B correction:** the earlier remark that the 101/69 W/cm² bracket makes V2019's 82 W/cm² "not uniquely a curvature effect" was a non sequitur. V2019 used shear-free, so a replication must explain its 82 W/cm² with shear-free ingredients (curvature, nonuniform evaporation, K accuracy). The bracket affects only the proposed-device prediction and the sensitivity study |
| 4 | P1 | **Linear Darcy unverified at device gradients, especially for the proposed cell.** Hand estimates at the capillary limit: reference Re_K ≈ 1.7; proposed Re_K ≈ 2-8, Re_h ≈ 10-40 (illustrative L = 10-2.5 mm). The quarter-cell pressure planes are exact only for Stokes flow. | `case_hand_checks.py` Sec. C; symmetry argument (guide §4) | K and q''_dry of the proposed design could be overestimated near the feed | Stage 1Q: periodic half cell, Navier-Stokes vs Stokes gradient sweep, predeclared 1% rule (contract) |
| 5 | P1 | **Evaporation closure and accommodation.** The Meland-Ytrehus treatment V2019 cites (Sec. 2.2, ref. 28, p. 172) is not held; our correction form is from memory. a = 0.052 was fitted (V2019 Sec. 3.4, p. 174) to the **last point** of Adera device 1 ("A1"; thesis p. 32), whose reported fluxes already have Adera's 62% loss deducted (adera2016 p. 286), at 24 °C / 3 kPa (p. 284) | held PDFs | `a` from 0.052 to 1 moves the Stage 2 reference superheat by ~1.3 K of ~8 K. The fitted `a` absorbs the loss correction, the 1D substrate inference and model-form errors ("absorbs the non-idealities", thesis p. 22) | **Rev. B:** our implementation reproduces the authors' own q''_e(ΔT) curve for a = 0.052 (thesis Fig. 7, p. 23) within ≤ 0.5% for ΔT > 3 K and ≤ 5 W/cm² everywhere. A least-squares fit of `a` to that curve with our form returns 0.052. Alternative forms (C4/a; far-field Knudsen-layer density) miss by > 25% (`tests/test_evaporation_thesis_fig7.py`). For the **reference replication** the closure gap is reduced to "functional match at a = 0.052 and 100 °C". The form at other `a` (sensitivity sweeps, proposed device) remains ASSUMED. The contract flag `accommodation_treatment_verified` stays false until the paper is held (no silent relaxation) |
| 6 | P1 | **Vapor environment undefined.** Unit-cell h_e assumes pure saturated vapor. In open air, evaporation is diffusion-limited. | V2019 Sec. 3.4 (air data used only for dryout); docs/assumptions A6 | heated tests in air cannot validate h_e or superheat predictions | `vapor_environment` UNRESOLVED in the case file; Exp-4 states what is comparable in each environment |
| 7 | P1 | **Dryout logic.** Reduced-order models with prescribed uniform evaporation, or an interface that stays wet regardless of supply, cannot establish dryout. Other limits can come first. **Rev. B correction:** the "boiling before dryout" entries are in adera2016 Table 1 (p. 282), a summary of other groups' wicks at 1 atm; Adera's own runs were "in the absence of nucleate boiling" (p. 285). EO adds flooding and electrochemical limits. | adera2016 pp. 282, 285; `docs/electroosmosis_model.md` | a "dryout gain" may be a limit that is never reached | Stage 6 couples supply and evaporation with conservation checks and reports which limit is reached first |
| 8 | P1 | **Reference → proposed contamination.** Post-processors fell back silently to the reference geometry when an export lacked d/l/hp columns. All stored K, F, EO and poster numbers are for the reference cell. | `analysis/_common.py` (before fix) | a proposed-cell export would be normalized with the l = 30 µm footprint | case file + loader with UNRESOLVED refusal; `--case` flags; warning on fallback; tests |
| 9 | P1 | **Downstream (cathode) reservoir hydraulic BC is undefined.** | CAD brief (separate cathode compartment, gas port) | sealed: EO raises end pressure (flooding risk); vented: EO drives net throughflow (subcooled feed adds sensible heat); gas-filled: pressure drifts | must be chosen per run (guide §11); `reservoir_geometry` UNRESOLVED |
| 10 | P1 | **Array geometry unresolved.** 205 × 155 at 50 µm pitch spans 10.22 × 7.72 mm (outer edges). Under either reading of the (−10, −7.5) mm start reference, the array sits ~4.9 mm and ~3.6 mm off the coupon center. | `start_reference_interpretations` | feed length, heater placement and any device model depend on the true layout; a 20 × 15 mm wick must not be substituted | recorded as UNRESOLVED; resolve in the native CAD application (scale, units, `Coupon_V2.cat` reference) |
| 11 | P2 | **Field treatment and electrical accounting.** Pujahari resolves the EDL (Poisson-Boltzmann, ψ = ζ on floor and pillars, Eqs. 5-6) and applies a ρ_e E body force with no-slip walls and a constant x-field per cell (Eqs. 7-8, p. 8446). It uses no slip velocity and solves no field around the pillars. Our uniform-vs-solved comparison (thin-EDL slip, Python 2D: ~0.5-0.6×) is the thin-EDL equivalent, approximate where λ_D/half-gap ~ 0.1 (their 1e-6 M cases at small pitch). The field in the wick is not supply voltage / length; DC current needs Faradaic reactions. | pujahari2022 pp. 8446-8447; `docs/eo_field_treatment_test.md` | EO magnitudes and power can be off by ~2× or more; bulk Joule power is not the electrical input | Stage 2A procedure (solved field, thin-EDL slip; resolved-EDL only if the screen fails); compare with Pujahari only as context |
| 12 | P2 | **EO stall pressure is small compared with reservoir head effects.** Illustratively, tens of Pa over 10 mm at 1e4 V/m and ζ = −50 mV, i.e. a few mm of water. | `case_hand_checks.py` Sec. E | unheated EO tests are dominated by head drift unless controlled | Exp-3 controls; matrix |
| 13 | P2 | **Stage 3 endpoint.** V2019 Eq. 3 is the mean force balance (exact identity of the Young-Laplace equation). The local edge angle reaches θ_rec earlier where the slope is steepest. Floor contact is not expected (sag ~10-12% of hp for both cells). | derivation in guide §8; hand estimate | dryout Pc may be overestimated if local depinning governs | report mean and local criteria; endpoint is a recorded choice |
| 14 | P2 | **Measured vs model temperature planes.** Experiments measure heater/backside temperatures through the substrate with losses; the model reports the pillar-base plane. | adera2016 (650 µm substrate, 1D correction, 62% loss) | direct comparison would mix substrate, spreading and loss errors into h_e | Exp-4 inference chain; Stage 7 trigger |
| 15 | P2 | **Thin-EDL validity** holds in the gaps for all screened chemistries (λ_D/half-gap ≤ 0.1, finite-EDL factor ≥ 0.91), but the Dukhin lower bound reaches ~0.1 for pure water, and the contact-line film is locally thinner than λ_D. | `case_hand_checks.py` Sec. E | ~10% EO error from surface conduction in very dilute water; local film region unmodeled | prefer a defined dilute electrolyte; resolve EDL only if the chosen chemistry fails the screen |
| 16 | P2 | **EO and the conduction approximation.** Joule ≤ 1e-5 of q''; vertical Pe ~1e-2; net flow is fixed by evaporation unless the cathode compartment allows throughflow. | `MODEL_STATUS.md`; V2019 App. A.1 | EO thermal coupling at the unit cell is likely negligible | one 4E check instead of an EO thermal stage |
| 17 | P2 | **Silicon conductivity** 148 W/m/K is a 300 K value from memory; the 373 K value is lower (also from memory, ~100-120). The sweep lower bound of 120 may be too high. Pillar conduction is a large share of the reference R''_wick. | `parameters/materials.yaml` | Stage 2/4 superheat biased low | keep ASSUMED; obtain a held source; widen the sweep when sourced (not changed silently) |
| 18 | P2 | **Three Python device-level models** (`device/`, `whole_wick/` with electrowetting and an ASSUMED h_e placeholder, `electroosmosis/wick1d`) with different assumptions. | repo | risk of quoting a reduced-order number as a prediction | roles fixed as cross-checks of COMSOL Stage 6; no extension planned |
| 19 | P2 | **Convergence threshold inconsistency:** modeling context said 1-5%; the contract requires 1e-3. Chase's stored Stage 1 finest pair is 1.02e-3; the student's is 5.05e-4. | docs | "converged" meant different things | preliminary (≤ 1e-2, labeled) vs qualified (≤ 1e-3) in the guide; a Richardson-based alternative is proposed, not adopted |
| 20 | P3 | ~~Old Stage 2 mesh table listed a "Corner Refinement" node (a 2D feature).~~ **Rev. B correction: wrong.** Chase's 3D Stage 2 mesh contains a working Corner Refinement node (`cr1`, scaling 0.35), so it exists in 3D in COMSOL 6.3. | `stage2_flat_interface.mph` tree | none | the guide's edge Size node is an alternative, not a fix; Corner Refinement may be kept |
| 21 | P3 | `analysis/permeability_from_comsol.py` rewrites `data/processed/permeability_<stem>.csv` on every run. | code | a rerun with different flags could overwrite a stored result | noted; run with a different export stem or commit first |
| 22 | P2 | **`dp_dry` typo in Chase's Stage 2 model (rev. B):** `4*gam*cos(theta_r)/(d*(4/pi*(1/d)^2-1))` uses `1/d` for `l/d` and evaluates to 1.8e-6 Pa. | model parameters | none at Pc = 0; wrong for any Stage 3 copy of this model | fix before reuse; compare with 2177 Pa |
| 23 | P2 | **The a = 0.052 fit is conditional on V2019 inputs it does not state (rev. B).** V2019 modeled Adera A1/A2 with "2D flow" (Fig. 6 caption, p. 175) but gives no θ_rec, property set or substrate model for that comparison. Adera measured θ_rec ≈ 70° on those devices (C4F8 residue; p. 288), not the 15° V2019 uses for SiO₂. Adera's own model fitted a = 0.08 with Schrage (Eq. 16, p. 289). | held PDFs | the same experiment yields a = 0.052 (BTE, V2019) or 0.08 (Schrage, Adera): `a` is closure- and model-dependent | treat `a` as a closure-specific fitting parameter; never transfer it between closures (Pujahari uses V2019's 0.052 inside Schrage, p. 8447) |

## 2. Stage 1 and Stage 2 evidence, as it stands

- **Stage 1 reference.** Chase (`permeability1.mph`, stored): K = 2.96979e-11 m² at 0.5 µm;
  changes 1.02%, 0.20%, 0.10% (coarser denominator). Student (reported): K = 2.9711e-11 m² at
  0.25 µm; changes 1.068%, 0.297%, 0.0505%. Matching meshes differ by −0.047%, −0.0007% and
  −0.006%; finest-available values differ by +0.044%. Both are 1.3% below the analytical flat-top
  estimate. This is reproducibility evidence for one flat-top, low-gradient case. It is not
  qualification (no gradient sweep, no interface bracket, `.mph` not tracked).
- **Stage 2 reference.** Two Q_in values (1.125e-4 W, the input echo) and two T_base values
  (381.19639 K, 381.21757 K). **Rev. B (model file inspected):** hi = 4.2e5 W/(m²K); trimmed
  pillar; heater on boundary 3, T_base over boundaries 6 and 10, meniscus boundary 11. Column 1 is
  the automatic-mesh solution; column 2 is probably hmax = 0.25 µm. The heater area is confirmed
  indirectly: with every other boundary insulated, ∫ht.ntflux over boundary 11 equals q0·(l/2)² to
  4e-12 on column 2. The T_base selection area is still unchecked. Superheat: 8.046 K (automatic
  mesh) and 8.068 K (probably 0.25 µm), 0.26% apart. That difference is between two different
  mesh families, not a finest-pair change.

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
2. `stage2_flat_interface.mph` (rev. B: hi, trimming, selections and column mapping resolved): still
   needed are the area of boundaries 6 + 10 (= 2.25e-10 m²?), confirmation that "Study 2/Solution 1"
   holds hmax = 0.25 µm, and an evaluation of all four sweep solutions.
3. Meland & Ytrehus (2003). Reduced by rev. B to: confirming the general form (other `a` values).
   At a = 0.052 our law matches the authors' own Fig. 7.
4. Proposed-device definitions: CAD scale/units and array placement, feed arrangement and length,
   substrate thickness, heater footprint, reservoir levels, cathode compartment BC, vapor
   environment, electrolyte.
5. Surface data on our oxide: θ_adv, θ_rec, ζ (streaming), aging behavior.

## 6. Next three concrete actions

1. **Chase** (revised after rev. B): in `stage2_flat_interface.mph`, set every Derived Value to
   Dataset "Study 2/Parametric Solutions 1", parameter selection All. Add a Surface Integration of 1
   over boundaries 6 + 10 (A_base) and of the inward heater flux (`-ht.ntflux`, sign-checked) over
   boundary 3 (Q_in). Export one table with an hmax column. Fix `dp_dry` (`l/d`). Save a `.java`
   copy. Screenshots are only needed if A_base is not (l/2)².
2. **Undergraduate:** do `docs/stage2_undergrad_assignment.md` on the reference cell (hand in
   model, selections, Q_in, Q_evap_bc, Q_evap_flux, T_base, A_base, T_lv, mesh study), then
   duplicate for the proposed cell.
3. **Either, in parallel:** Stage 2A runs R0-R2 on the reference cell (guide §7), with Ũ·F = 1
   as the pass/fail check, followed by Stage 1 + 1Q on the proposed cell (inertia sweep and
   interface bracket).

## 7. Revision B (2026-11-13): source-level re-audit of the consequential claims

Page numbers are journal pages. PDF-to-journal offsets: vaartstra2019 +169, adera2016 +279,
pujahari2022 +8441, fleming2020 +43265; the thesis uses its printed page numbers. Labels:
**M** = measured by the authors; **A** = assumed or chosen by the authors; **I** = inferred by the
authors from a model fit or comparison; **C** = calculated by us. "Affects" says what the finding
changes: R = reference replication of V2019, P = proposed-device model, S = only an uncertainty or
sensitivity study.

### 7.1 No-slip vs shear-free liquid-vapor interface

| Claim | Evidence | Label | Affects | Status vs rev. A |
| --- | --- | --- | --- | --- |
| V2019 uses a shear-free interface for K(Pc) | "setting it as a symmetry boundary, which enforces no flux and no shear stress at the interface" (Sec. 2.2, p. 172; also Fig. 2c caption) | A | R (must keep shear-free) | unchanged |
| Fleming found no-slip fits better | "we find much better agreement using a no-slip condition" (abstract, p. 43266); "This finding contradicts typical analytical or numerical methods for estimating permeability" (p. 43270) | **I**, from M curvature | P, S | **sharpened**: an inference through their Darcian model (Eqs. 1-3, p. 43268), with evaporation neglected in the flow (p. 43269); not a direct measurement of interface velocity |
| Fleming conditions | ~20 mbar vapor at room temperature (p. 43268); Piranha-cleaned (p. 43267); d/pitch/height 60/90/50 and 80/120/50 µm, φ = 0.35 (Table 1, p. 43267) | M | — | **new**: their geometry and state differ from both our cells and from V2019's 100 °C |
| Proposed mechanism | surface-tension gradients, "surfactant contamination or mobile surface charge", thermal Marangoni (p. 43269) | A (hypothesis) | S | unchanged |
| −31% / −28% flat-K bracket for our cells | Brinkman closed-channel factor (`brinkman_no_slip_top_factor`) | **C** | P, S | relabeled: our number, not Fleming's |
| "82 W/cm² is not uniquely curvature" | — | — | — | **withdrawn** (see item 3) |

### 7.2 Origin and applicability of a = 0.052 and the 62% heat loss

| Claim | Evidence | Label | Affects | Status |
| --- | --- | --- | --- | --- |
| a = 0.052 fitted to Adera A1, checked on A2 | "determined the accommodation coefficient a = 0.052 by fitting to one set of data from [12] (A1) and then verified ... (A2)" (Sec. 3.4, p. 174) | I | R | unchanged |
| Fit used one point | "We chose to fit to the last point on (A1) ... the lowest heat fluxes in this experiment are offset (the intercept is not zero)" (thesis p. 32) | I | R, S | **new** |
| `a` is not a physical coefficient | "not intended to represent the theoretical evaporation/condensation coefficients because it absorbs the non-idealities of the experiments" (thesis p. 22) | A (authors' statement) | P (do not transfer), S | **new**; strengthens item 5 |
| A1/A2 = Adera devices 1 and 5 | V2019 Fig. 6 caption (p. 175): A1 d/l/h = 5/12/82, A2 = 12/20/90; adera2016 Table 2 (p. 284) | M (geometry) | R | **new** |
| 62% heat loss | "The combined heat loss via sensible, convection, and ohmic was 62% of the input heating power. This heat loss was deducted from the input power in all the heat fluxes reported" (Sec. 4, p. 286); characterization shown "for device 1" (Fig. 6 caption, p. 286) | M (sensible, ohmic) + **A** (convective, "conservatively for a 10 W/m² K heat transfer coefficient", p. 286) | R (inherited by the a fit), S | **sharpened**: mixed measured/assumed; whether one 62% figure was applied to all devices is not stated explicitly beyond "all the heat fluxes" |
| Wall superheat from 1D conduction | Eq. 2, 1 µm SiO₂ + 650 µm Si (p. 285) | I | R (the fit target), S | unchanged |
| Chamber state | "maintained at 24 C and corresponding saturation pressure 3 kPa"; noncondensables pumped out (pp. 282, 284) | M | R (fit at 24 °C, used at 100 °C), S | unchanged |
| Adera's own θ_rec ≈ 70° | interferometry vs Surface Evolver, 65°/75°, "average receding contact angle 70°"; C4F8 residue (Sec. 5.1.3, p. 288) | I | R (V2019's fit inputs unstated), P (our oxide state matters) | **new** (item 23) |
| Adera's own a = 0.08 (Schrage) | "the experimentally obtained average accommodation coefficient (0.08)" (Eq. 16, p. 289) | I | S | **new**: `a` depends on the closure |
| "Boiling before dryout" in Adera | Table 1 (p. 282) summarizes other groups; Adera: "in the absence of nucleate boiling" (p. 285) | — | — | **corrected** (item 7) |

What this does **not** establish: that V2019's h_e model is experimentally validated. The A2 check
(thesis Fig. 13, p. 32) is one independent device. It is compared on mean superheat, which the
authors call "not an entirely satisfactory validation" (thesis p. 33). It inherits Adera's loss and
1D-conduction corrections.

### 7.3 Meland-Ytrehus correction and our evaporation law

| Claim | Evidence | Label | Affects | Status |
| --- | --- | --- | --- | --- |
| V2019 uses the linearized BTE (Sone 2000) with p_k = p_∞ | Eqs. 5-7, "constraining the pressure at the Knudsen layer boundary to equal the far field pressure" (p. 172) | A | R | unchanged |
| The a = 1 solution must be corrected | "The macroscopic flow parameters determined by Equations 6,7 correspond to the 'total accommodation' case and must be corrected for fluids such as water" (thesis p. 21); method of Meland & Ytrehus (Phys. Fluids 15 (2003) 1348) cited but not printed (p. 172; thesis p. 22) | A | R | unchanged: form not printed anywhere held |
| Our law reproduces the authors' curve at a = 0.052 | thesis Fig. 7 (p. 23) digitized (`data/literature/vaartstra2019_thesis_fig7_digitized.csv`); our `sone_bte_heat_flux` within ≤ 0.5% for ΔT > 3 K, ≤ 5 W/cm² (1% of span) below; fitted a = 0.052; C4/a and far-field ρ_k variants are > 25% off (`tests/test_evaporation_thesis_fig7.py`) | **C** vs authors' plotted **calculation** | R | **new**: the closure gap for the reference replication is reduced to a functional match. Not experimental evidence. T_sat = 100 °C is inferred (figure unlabeled) |
| The correction form at other `a` | not testable from held sources | A (ours) | P, S | unchanged |
| Dropping Eq. 6 (temperature jump) | at u_k/√(2RT) ~ 1e-4 the T_k correction to ρ_k is negligible; the Fig. 7 match supports this | C | R | **new** |

### 7.4 Geometry, averaging and normalization needed to reproduce V2019

| Item | Evidence | Label | Our treatment |
| --- | --- | --- | --- |
| Domain | "Only the region enclosed by the dashed lines is treated as the domain" (Fig. 2a, p. 172); the thesis drawing shows a 1/8 triangle (Fig. 4a, p. 20). No separate flow domain is stated; Fig. 2c shows a side view between two pillars with p1, p2 | A | quarter cell (exact for Stokes, §4); triangle also valid for thermal/meniscus; flow-domain choice of V2019 unknown |
| Heat input | "1 µm below the base of the pillars to account for constriction resistance" (Sec. 2.2, p. 171) | A | `tsub = 1 µm` |
| T_e | "average temperature at the surface (red line)" (Fig. 2b, p. 172); red line spans the full width at the pillar-base plane (thesis Fig. 4b, p. 20) | A | area-weighted average over the whole z = 0 plane (PROJECT_DERIVED reading; "area-weighted" is COMSOL's Surface Average definition, not stated by V2019) |
| h_e | Eq. 4, h_e = q''_in/(T_e − T_sat), q''_in at the cell bottom (p. 171) | A | footprint basis |
| K normalization | Darcy Eq. 8 (p. 172); "the mass flux (Γ) was reported" without stating the area; device Eq. 9 uses Q_m = q''_e/(h_lv h) over the full wick height (p. 173) | A + **C** (our reading) | superficial over (l/2)·hp. Eq. 9's use of h supports a bulk (superficial) K, but this is our inference |
| Re | "Reynolds number Re_h = ρuh/μ ~ 10¹" (p. 172; App. A.1, p. 175) | C (authors) | Stage 1Q gradient sweep |
| Device | L = 5 mm, t_s = 100 µm, pure vapor at 100 °C, θ_rec = 15° (Sec. 3.2, p. 173); symmetry at x = 0 and x = L "so long as liquid is supplied at T_sat" (p. 173) | A | reference case file |
| Mesh | Table A1 (p. 176), evaluated "at the most extreme interfacial curvature" | C (authors) | contract thresholds |
| Not stated | silicon conductivity, water property source, flow-domain BCs, θ_rec and property set for the Fig. 6a Adera comparison | — | swept or blocked |

### 7.5 Electric-field and EO assumptions in Pujahari (simulation only)

| Claim | Evidence | Label | Affects | Status |
| --- | --- | --- | --- | --- |
| EDL resolved, not slip | "Unlike many studies, here, the effect of electro-osmosis is not just added ... by using the Helmholtz-Smoluchowski slip velocity boundary condition, but the EDL distribution of ions is calculated and then the volumetric Coulombic force is added" (pp. 8444-8445); PB Eqs. 5-6, "Ψ = ζ at the bottom wall (z = 0) and at the pillar walls", zero normal derivative at the meniscus (p. 8446) | A | S (comparison context) | **clarified**: their method is the COMSOL blog's "resolve the EDL" route; our Stage 2A uses the thin-EDL slip route; the two are equivalent only for λ_D ≪ gap |
| Uniform field | "As the double derivative of external potential is zero, the external electric field (E) obeys a linear profile ... a constant electric field is considered to be acting along each unit cell" (Eq. 8, p. 8446); BC "E = dφ/dx at x = 0 and l plane" | A | S | confirmed: no field solved around insulating pillars |
| Walls | "u, v, w = 0 at the bottom wall and at the pillar walls" (p. 8446) | A | S | no slip with body force, so no double counting |
| ζ and chemistry | 1e-6 and 1e-4 M KCl "in agreement with the experimental work of Mala et al., resulting in zeta potentials of −200 and −100 mV" (p. 8446) | A (borrowed from another surface) | P (cannot transfer), S | unchanged |
| ε(T) | ε(T) = 305.7 exp(−T/219) (p. 8446) | A | S | new detail |
| Electrolysis | "the voltage values will be very less compared to the voltage that causes electrolysis" (pp. 8446-8447); "should not be high enough to result in electrolysis" (p. 8453) | A (unquantified) | P | unchanged: E × 5 mm = 25-500 V is our calculation (C) |
| Power | "order of milliwatts (Supplementary Information)" (p. 8452); Joule ≤ 4.5e-8 W per domain (p. 8453) | C (authors; SI not held) | P | unchecked |
| Thermal closure | Schrage (Eq. 9) with "accommodation coefficient (0.052)" (p. 8447) | A | S | **new**: V2019's BTE-fitted `a` reused in Schrage, a closure mismatch (~3% flux, C) |
| "Validation" | EO vs an analytical 2D channel (Dutta & Beskok) (Sec. 3.1.1, p. 8448); thermal vs Zhu (dryout) and Adera via V2019 (Fig. 3, p. 8449), where Fig. 3b uses d/l/h = 5/12/82 µm = Adera device 1 = the A1 set V2019 used to **fit** `a` | C | — | **new**: the temperature check is not independent of the fit; no EO experiment exists in the paper |
| Velocity definition | U_i = ∬ u dy dz over the midplane liquid area A (Sec. 2.1.2, p. 8447); M_i = ρ A_i U_i (Eq. 13) | A | S | **new**: their U is a liquid-area (intrinsic) average, not superficial |

### 7.6 COMSOL guidance (re-read)

- comsol_blog_eof_edl: HS slip "u = ε_w ζ₀ ∇_T V / η"; "first calculate the Debye length ... If this
  length is much smaller than your geometrical length scale, use the electroosmotic velocity
  boundary condition. If not, use the traditional no-slip velocity wall boundary condition and
  resolve the flow in the EDL." Pujahari follows the second route; Stage 2A the first.
- comsol_kb_1261: mesh-size parameter plus Parametric Sweep; "integrals over the entire model will
  converge the most quickly, and local values the most slowly"; derivatives at sharp corners are
  "locally non-convergent". Chase's Stage 2 sweep follows this. Only the evaluation step
  (dataset selection) was missed.
