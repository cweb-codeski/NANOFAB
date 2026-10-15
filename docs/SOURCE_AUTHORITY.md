# Source Authority and Provenance

This document defines how NANOFAB decides which source to trust for which question, and how every
value records where it came from. The rules are enforced in code by `src/nanofab/provenance.py`,
`src/nanofab/params.py` and `tests/test_provenance.py`.

Files involved:

| File | Contents |
| --- | --- |
| `config/sources.yaml` | Registry of every source id (papers, databases, documentation, project notes) |
| `parameters/*.yaml` | Model inputs: `value`, `unit`, `provenance`, `source`, `source_location`, `notes` |
| `parameters/design_targets.yaml` | Project design targets only (kept apart from replication inputs) |
| `config/model_provenance.yaml` | Provenance of equations, boundary conditions, interpretations, derived results |
| `data/literature/vaartstra2019_targets.csv` | Replication targets, each with provenance and location |

## Provenance classes

| Class | Meaning | Needs a `source` id? |
| --- | --- | --- |
| `PRIMARY_SOURCE` | Stated explicitly in the text, a table or an equation of a paper we hold (or one whose citation we copied from a held paper's reference list). `source_location` gives the section, equation, table or figure. | yes |
| `AUTHORITATIVE_DATABASE` | Taken from a reference database or official documentation: IAPWS/NIST properties, CODATA constants, COMSOL documentation. Textbook property tables are filed here as secondary compilations. | yes |
| `PROJECT_DERIVED` | Inferred by us from a source's stated values: digitized from a figure, read from an ambiguous passage, or simple arithmetic on reported numbers. The source it was derived from goes in `source`. | yes |
| `MODEL_CALCULATED` | Produced by this repository's models or COMSOL runs. **Never presented as literature**, and never allowed in the literature targets file. | no (`null`) |
| `ASSUMED` | Not supported by a held source: a modeling choice, an engineering guess, or a value from a citation we have not checked. If an unverified citation exists, it goes in `source` with a note. | optional |
| `DESIGN_TARGET` | A project goal (for example 100 W/cm²). Lives only in `parameters/design_targets.yaml`. | no |
| `PROJECT_CONTEXT` | Chase's notes, the project brief, Drive documents. These explain intent and priorities; they are not scientific evidence. | yes |

A value's provenance describes how **this value** was obtained, not how reputable the source is.
For example, the 19% permeability reduction uses the PRIMARY_SOURCE value of 82 W/cm², but the
reduction itself is MODEL_CALCULATED because it comes from our 1D model.

## Authority order for the baseline replication

When sources conflict on a question within their scope, prefer them in this order:

1. **Vaartstra, Lu & Wang, IJHMT 136 (2019) 170** (`vaartstra2019`, role `PRIMARY_MODEL_REFERENCE`).
2. Other primary thin-film / micropillar papers (`zhou_adera2026`, `vaartstra2019_thesis`,
   `vaartstra2020`, `ravi2014`, and the papers vaartstra2019 cites).
3. COMSOL documentation (`comsol_docs`) for implementation questions.
4. NIST / IAPWS (`nist_webbook`, `iapws` via `coolprop`) for fluid properties where the paper does
   not specify them.
5. Project-derived values.
6. Assumptions.

Scope of each authority:

- **vaartstra2019 governs** the baseline geometry, governing equations, boundary conditions,
  modeling assumptions (App. A.1) and replication targets. It does **not** govern COMSOL
  workflow details beyond what it states (it used COMSOL 5.3a). It also says nothing on material
  properties it does not state (silicon conductivity, the water property source).
- **COMSOL documentation governs** how to implement a physics choice in COMSOL: interfaces,
  boundary-condition features, solver and mesh settings. It does not override the paper's
  physics choices.
- **NIST / IAPWS govern** fluid properties in new predictive models. Use a paper's exact property
  set only when reproducing its numbers exactly; vaartstra2019 does not state one.
- The **thesis** ranks below the IJHMT paper (not peer reviewed) but may clarify ambiguities. When
  it disagrees with the paper, record both (rule 9).
- **Review papers** (`vaartstra2020`) are primary only for their own new results. For values they
  compile, cite the original paper.

## Parameter schema

```yaml
target_heat_flux:
  value: 100
  unit: W/cm2
  provenance: DESIGN_TARGET
  source: null
  source_location: null
  notes: "Project performance target; not a value from vaartstra2019."
```

`value`, `unit`, `provenance`, `source` and `source_location` are required (use `null` when not
applicable). `notes` and `sweep` are optional. Overriding a value in a script
(`ParameterSet.with_overrides`) relabels it `MODEL_CALCULATED`.

## Script output

Scripts label results by provenance, for example (`analysis/vaartstra_replication.py`):

```text
   [MODEL_CALCULATED] q''_dry = 100.6 W/cm^2 (1D Darcy, flat-meniscus K ...)
   [PRIMARY_SOURCE] vaartstra2019 target = 82 W/cm^2 (Sec. 3.2)
   relative error vs vaartstra2019: +22.7%
```

## Audit of values attributed to Vaartstra et al. (2019)

The audit was done against the held PDF, reading the full text.

### Confirmed as PRIMARY_SOURCE (explicitly stated)

| Value | Location |
| --- | --- |
| d = 10, l = 30, h = 25 µm, L = 5 mm, t_s = 100 µm | Sec. 3.2, Fig. 4 caption |
| 1 cm × 1 cm wick (→ wick_width 10 mm) | Sec. 1 |
| Heat input 1 µm below the pillar base | Sec. 2.2 |
| θ_rec = 15° (water on SiO₂, citing raj2013) | Sec. 3.1 |
| T_sat = 100 °C, pure vapor ambient | Sec. 3.1, 3.2 |
| q''_in = 50 W/cm² for the unit-cell sweeps | Sec. 3.1 |
| q''_dryout = 82 W/cm² | Sec. 3.2 |
| h_e > 60 kW/m²K; h_e varies 30%; ΔT_e 3 °C higher at the inlet | Sec. 3.2, Fig. 4 |
| a = 0.052, fitted to Adera A1 | Sec. 2.2, 3.4 |
| C4* = 2.13204, d4* = 0.44675 | Sec. 2.2 |
| 0.082% h_e change for a 10% change in q''_in | Sec. 3.1 |
| "scaling with h²"; d/h = 0.4, l/h = 1.2 | Sec. 3.3 |
| Validation geometries; Adera at 3 kPa, 650 µm substrate | Fig. 6 caption, Sec. 3.4 |
| Table A1 mesh statistics (minus signs lost in PDF text, restored) | App. A.2 |
| Symmetry check 0.19% (uniform), 1% (hotspot) | App. A.1, A.4 |

### Attributed to the paper but not stated by it

| Value | Was presented as | Actual provenance | Finding |
| --- | --- | --- | --- |
| Dryout heat flux 100.6 W/cm² | replication result | MODEL_CALCULATED | Our flat-meniscus 1D model. +22.7% vs the paper's 82. |
| "K about 19% below flat" | "the paper implies" | MODEL_CALCULATED; the cause is a PROJECT_DERIVED hypothesis | The paper states no K reduction, only that curvature lowers K. Our docs said "this is the curvature effect", but analytical K error and non-uniform evaporation may contribute. Reworded everywhere. |
| "Kinetic resistance ≤ 16% of total" | "V2019 reports … so the kinetic part is at most ~16%" | MODEL_CALCULATED; **inequality reversed (corrected)** | R_int/R_total = R_int h_e and h_e > 60 kW/m²K, so the share is > R_int × 60000 ≈ 16%: a conditional *lower* bound, not an upper bound. Conditional on our R''_int (flat interface area, unverified accommodation correction). The paper gives no resistance breakdown. Corrected in code, docs and `config/model_provenance.yaml`. |
| "Flat K is a known upper bound" on dryout | assumption status (A2) | Unsupported claim (removed) | V2019 says curvature lowers K; it does not follow that the legacy constant-K, uniform-q estimate bounds the coupled result. Now an approximate analytical estimate. |
| Reduced BTE called "paper-exact" | docstring / printed label | Unsupported claim (removed) | The implementation drops the temperature-jump relation (Eq. 6), approximates rho_k, and uses an accommodation correction written from memory (Meland & Ytrehus 2003 not held). BTE-vs-Schrage agreement (~3-5%) is a cross-model check only. |
| Silicon k = 148 W/(m·K) | AUTHORITATIVE_DATABASE | ASSUMED (provisional) | A textbook value recalled for 300 K, never re-checked against the table, not verified at 373.15 K, and not a V2019 value. |
| T_e = mean temperature of the pillar-base plane | the paper's definition | PROJECT_DERIVED (well supported) | The text says only "average temperature at the surface (red line)". The thesis drawing (Fig. 4b, p. 20) puts the red line along the substrate top at the pillar base. |
| Re_h ≈ 10 (exponent +1) | "From the stated values Re_h ≈ 16, so the exponent is +1" | Now PRIMARY_SOURCE (thesis p. 23, 51) | The IJHMT text was ambiguous; the thesis page images print "~ 10¹". Our earlier reasoning was not proof, but the conclusion holds. |
| C4* = 2.13204, d4* = 0.44675 | as printed | PRIMARY_SOURCE, sign convention noted | The thesis prints them as negative. `evaporation.py` uses magnitudes with an explicit minus sign, so the physics is equivalent. |
| Δp_dryout = 2177 Pa | note on a target row | MODEL_CALCULATED | Eq. (3) is primary; the number is ours. Moved out of the targets file. |
| Adera case T = 297.2 K | target condition | MODEL_CALCULATED | T_sat(3 kPa) via IAPWS. The 3 kPa is primary. |
| Zhu case T = 298.15 K | target condition | ASSUMED | Not stated (air ambient). |
| Meniscus "computed with Surface Evolver in Vaartstra et al." | docstring in `capillary.py` | **Misattribution** | The paper solves Eq. (2) in COMSOL (Sec. 2.1). Docstring corrected. No physics affected. |
| Pitch optimum ≈ 40 µm | finding | MODEL_CALCULATED | The paper calls l = 30 µm "near-optimal", citing zhu2016. |

### Other parameters whose classification changed

| Parameter | Old label | New label | Reason |
| --- | --- | --- | --- |
| silicon_thermal_conductivity (148) | LITERATURE | ASSUMED, provisional (intended source `incropera_fhmt`) | Not in vaartstra2019. Textbook value from memory for 300 K, table not re-checked, not verified at 373.15 K; relabelled from AUTHORITATIVE_DATABASE. vaartstra2020 uses 153 at 20 °C; both are kept, and the sweep was widened to 120-153. |
| silicon_density, silicon_specific_heat | LITERATURE | AUTHORITATIVE_DATABASE (`incropera_fhmt`) | Same as above. |
| gravity (9.81) | LITERATURE | AUTHORITATIVE_DATABASE (`nist_si_constants`) | Rounded from g_n = 9.80665 m/s². Value unchanged (negligible effect, Bo ~1e-5). |
| molecular_kinetic_diameter (0.265 nm) | LITERATURE | ASSUMED (`breck1974`, unverified) | Citation from memory, not held. |
| advancing_contact_angle (30°) | ASSUMED | ASSUMED | No source. Not used by the baseline. |
| heat_flux (100 W/cm²) | TARGET | DESIGN_TARGET, renamed `target_heat_flux` and moved to `design_targets.yaml` | Separated from replication inputs. |
| accommodation_coefficient | LITERATURE | PRIMARY_SOURCE | Confirmed in the paper. Range now cited to vaartstra2020 Sec. III C instead of an unverified Marek & Straub citation. |
| Geometry, θ_rec, T_v, unit-cell and dryout heat flux | LITERATURE | PRIMARY_SOURCE | Confirmed with section numbers. |

### Citation problems found (documented before changing)

1. **"Byon et al., IJHMT 79 (2014)"** for the Brinkman height correction (`permeability.py`,
   `equations.md`, `literature_notes.md`) was written from memory and cannot be verified. The held
   papers cite Byon & Kim, J. Micromech. Microeng. 21 (2011) for meniscus effects on K, and
   Brinkman (1949) / Adera et al. (2016) for the Brinkman approach. The citation was replaced; the
   formula is unchanged and now labeled as our own derivation.
2. **Schrage year.** Our docs said 1953; vaartstra2019 ref. [22] prints 1954. Both are recorded
   (`schrage1953` notes); unresolved.
3. **Marek & Straub (2001)** for the accommodation range: not held. Replaced with vaartstra2020
   Sec. III C (held), which gives 0.001-1 for water.
4. **Surface Evolver** misattribution (see above).
5. **Meland & Ytrehus correction form** `C4 + 2√π(1−a)/a`: vaartstra2019 cites the method without
   printing the form. The form is ASSUMED until checked.

No equation, constant or numerical result changed as a result of this audit.

## Unresolved source questions

- Silicon conductivity and water property source used by vaartstra2019. Neither the paper nor
  the thesis states them, so they are probably not recoverable from held sources. Keep sweeping k_Si.
- Exact Meland-Ytrehus correction form. Neither the paper nor the thesis prints it, so it needs
  meland_ytrehus2003 itself.

Resolved by the thesis page images: Re_h exponent (+1); T_e surface (substrate top at pillar
base, from the figure); heat input depth 1 µm below the pillar base.
- Sangani-Acrivos / Drummond-Tahir coefficients, Lee et al. (1995) spreading formula,
  Carey/Faghri editions, Breck kinetic diameter: all cited from memory.
- Which IAPWS releases CoolProp 8.0 implements for each water property.
- The COMSOL guide steps have not been checked against COMSOL documentation.
- Year of SAND92-2347C; DOI of narayanan2024_interpack.

## Rules for Future Modeling Work

1. Never present an inferred or calculated value as though it were directly reported by a paper.
2. Prefer primary papers over project notes for scientific claims.
3. Prefer COMSOL documentation for COMSOL-specific implementation questions.
4. Prefer NIST/IAPWS for fluid properties in new predictions unless reproducing a paper's exact
   property set.
5. Use Chase's Drive notes and project docs as context, not as primary scientific evidence.
6. Any unsupported value must be marked ASSUMED.
7. Any assumption that materially affects performance should be included in a sensitivity
   analysis.
8. Preserve exact replication values separately from project design targets.
9. When multiple sources disagree, preserve both and document the discrepancy rather than
   silently choosing one.
10. Every important model output should be distinguishable as published, derived, calculated, or
    assumed.
