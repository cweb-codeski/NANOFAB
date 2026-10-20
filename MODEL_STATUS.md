# MODEL_STATUS — Passive vs Electroosmotically Assisted Micropillar Wick

Revision 2 (2026-10-18): corrected audit after re-reading pujahari2022 and digitizing its Fig. 6.
Equations and assumptions: `docs/electroosmosis_model.md`. Regenerate with
`python analysis/eo_assisted_wick.py` (≈3 min) → `results/electroosmosis/` (gitignored).
Re-digitize with `python analysis/digitize_pujahari2022_fig6.py` (needs dev-only PyMuPDF).
Tests: `pytest` → **211 passed** (163 pre-existing + 38 electroosmosis + 10 planform field-treatment).

**Classification of all EO results:** *idealized/conditional electroosmotic hydraulic predictions
under an imposed electric field.* They are not device predictions.

## Corrected audit

### Conclusions that changed

| Previous statement | Correction |
| --- | --- |
| Figures and text described the liquid only as "water at 100 °C" | The EO working fluid is **dilute aqueous KCl** (1e-6 M with ζ = −200 mV; 1e-4 M with ζ = −100 mV). Water-like thermophysical properties (IAPWS) are used on pujahari2022's dilute-solution basis (end of Sec. 2.1.4). Labels corrected. |
| "Electrolysis voltage exceeded 20–400×"; 1.23 V drawn as a limit | Withdrawn as a validity cutoff. 1.23 V is an electrode-interface thermodynamic quantity, not a limit on the bulk drop. The defensible statement is current-based: a steady DC current of 1.7 µA (1e-6 M, Table 1 σ) to ~100 µA (1e-4 M, KCl-corrected σ) per half-wick at 1e5 V/m must be carried by electrode reactions unless the electrodes are capacitive, AC or isolated. That is a major open feasibility question, not an invalidation of the hydraulics. |
| "pujahari2022's electrolysis claim is inconsistent" | Restated. The paper says the voltages "will be very less compared to the voltage that causes electrolysis" (Sec. 2.1.2) and that the field "should not be high enough to result in electrolysis" (Sec. 3.5). It **does not model electrodes**: E is imposed uniformly per unit cell (Eq. 8, 1D). V_bulk = E × 5 mm = 25–500 V is **our engineering extrapolation**. If that is the right length, it is hard to reconcile with Sec. 2.1.2, but the paper does not specify the configuration. |
| "Salt accumulates without bound; no steady state; enrichment ~10^1000" | Wrong in substance; withdrawn. With Pe ≫ 1, salt-flux conservation gives a **quasi-steady advective enrichment** c/c_feed = 1/(1 − x/L) over most of the wick. Unbounded accumulation is confined to a diffusive end layer ~6–22 µm thick. At the Case B dryout-critical location (0.04–1.0 mm from the closed end), c/c_feed ≈ 5–130. The issue is a non-uniform electrolyte (λ_D, σ, ζ at x*), not the absence of any steady state. |
| Cross-check: "agreement for l ≥ 20 µm, unresolved over-prediction at small pitch" | The text ranges ("1.4–1.45", "1.74–1.8") are the **maxima** of U-shaped curves, and I had overlaid them as bands across every pitch. With the digitized Fig. 6 compared per pitch, our model over-predicts their gain at **every** pitch (by up to 32% at h = 25 µm). The cause is decomposed below. |
| Conductivity used as given | Table 1's single 0.79 µS/cm is plausible for 1e-6 M (KCl ≈ 0.4–0.5 µS/cm plus intrinsic H⁺/OH⁻) but ~50× low for 1e-4 M (≈40–48 µS/cm). Hydraulics are unaffected; current, Faradaic rate and Joule heat are ×50 for that case (Joule still ≤ 8e-4 of q''). |
| Debye length taken as stated | The values reproduce from c (0.286 / 0.0286 µm). For 1e-6 M, intrinsic water ions at 100 °C (pK_w ≈ 12.27, ASSUMED) shorten λ_D to ≈ 0.22 µm. |

### Conclusions that survived

- Q_total = Q_pressure + Q_EO is exact under thin-EDL, uniform-ζ, linear-Stokes assumptions. K is
  geometric and identical in both cases. F = 1.1912 is verified numerically (0.008%).
- The dryout gain depends only on Pi = U_EO/U_cap, with G = (1 + Pi + √(1 + 2Pi))/2. This structure
  is now **confirmed against pujahari2022's own full PB-NS results** (below).
- The required bulk potential drop V_cap = F K Pc,max/(ε|ζ|) is independent of L (≈ 785 V for
  Pi = 1 at −200 mV for our wick), now framed as a bulk drop.
- A new over-pumping (flooding) limit exists at low heat flux in Case B.
- Diagnostics unchanged in substance:
  - Re_K 1.7–2.7 at the feed in both cases.
  - Dukhin lower bound 0.61 for the 1e-6 M / −200 mV case.
  - Field-effect charge through a thin oxide comparable to or larger than the native charge at
    high E, if the substrate conducts.
  - Joule heating negligible.

## A. Reproduction of pujahari2022 (d = 5 µm, L = 5 mm, h = 5 µm at l = 8–18 µm, h = 25 µm at l = 10–45 µm, ζ = −100/−200 mV, E = 5e3–1e5 V/m)

Source: digitized Fig. 6 (`data/literature/pujahari2022_fig6_digitized.csv`, PROJECT_DERIVED,
±0.5% of axis span). Their duplicated passive curves agree within 0.2%.

| Factor | Result |
| --- | --- |
| Field dependence | **Reproduced.** Pi from their 1e5 V/m curve, scaled ∝ E through our G(Pi), predicts their 5e3, 1e4, 5e4 V/m curves: median error 0.35%, max 5.4%. |
| Trends with pitch | **Reproduced in shape.** h = 25 µm gains are U-shaped in both models. For h = 5 µm their gain rises with pitch; it follows once their lower passive K is used. |
| Passive K | Their passive dryout implies K/K_flat = 0.75–0.97 (h = 25) and 0.44–0.67 (h = 5), falling with pitch. Consistent with meniscus curvature, which our flat-top K omits. |
| EO-flux magnitude | **Not reproduced.** Their EO flux is R = 0.52–0.66 (h = 25) / 0.60–0.69 (h = 5) × mu_EO E/F. |
| Gain magnitude at 1e5 V/m | Our flat model is off by up to 24% (h = 5) / 32% (h = 25). The zero-parameter uniform-imposed-field hypothesis with their implied K gets within 21% / 9%. A single fitted R (0.67 / 0.58) gets within 5%. |

**Small-pitch discrepancy (priority item).**

- **Part 1, misreading:** I compared against range maxima.
- **Part 2, two separable physical differences:**
  1. Their passive K is lower, which raises their Pi, mostly at large pitch and short pillars.
  2. Their EO flux is ~0.5–0.7× ours, which lowers their Pi everywhere and dominates at small
     pitch and tall pillars.
- **Hypothesis for difference 2:** a uniform imposed x-field drives sidewall slip E sinθ instead
  of the ≈2E sinθ an insulating wall requires. That is exactly half in the dilute 2D limit; the
  floor slip is screened by √K_s. It gives R = 0.54–0.80. It matches the inferred R at small pitch
  but over-predicts it at large pitch, so it is a strong partial explanation, **not** confirmed.
- **Which baseline is physical is unresolved:** our solved-field value or their ~0.6× value.
- **Independent 2D test (revision 3, `docs/eo_field_treatment_test.md`):** a resolved Laplace +
  Stokes planform calculation with identical HS slip, ζ, μ and mean field gives
  Q_uniform/Q_solved = 0.50–0.63 (R_field = Q_solved/Q_uniform = 1.59–1.98) in the tall-pillar
  limit, with no fitting.
  - It agrees with R_inferred within 0.03 for h = 25 µm at l = 12.5–25 µm.
  - 23 of 26 inferred points lie between the exact 2D limit and a heuristic floor bridge.
  - Verified: current and mass conservation (~1e-13), similitude (Case B = 1/F), F vs the FV
    solver (7e-5), K vs Sangani–Acrivos (0.05%), reversal, E → 0, and spectral convergence.
  - Conclusion: field treatment explains the **magnitude**. The pitch and height trend at large
    pitch and short pillars needs 3D (floor and meniscus).
  - Physically, insulating pillars under a DC field require the solved field, so mu_EO E/F remains
    our baseline.

## B. Our current wick (d/l/h = 10/30/25 µm, L = 5 mm, COMSOL K), at E = 1e5 V/m

The band spans the two unresolved factors:
- **EO-flux model:** solved field (R = 1) vs pujahari-equivalent R = 0.52–0.66.
- **Passive K:** COMSOL flat-top K vs the V2019-implied K (82/99.3 = 0.83 × K, PROJECT_DERIVED).

| ζ | Pi | G = q''_B/q''_A | q''_A [W/cm²] | q''_B [W/cm²] |
| --- | --- | --- | --- | --- |
| −100 mV (1e-4 M) | 0.17–0.39 | 1.16–1.36 | 82–99 | 98–129 |
| −200 mV (1e-6 M) | 0.33–0.77 | 1.31–1.68 | 82–99 | 112–156 |

At 5e3–1e4 V/m the gains are 1.6–6% (solved field) and smaller if R ≈ 0.6.

## C. Extrapolation to a proposed fabricated geometry

**Not done.** No fabricated geometry is defined in the repository. The geometry maps
(`fig5_regime_map`, `fig6_required_field_vs_geometry`, `fig11_voltage_scale`) give the regime
for any (d, l, h, L), with the B-band caveats (R and passive K) applied on top.

## Electrochemical and E-to-voltage status

- **What pujahari2022 actually models:**
  - a uniform imposed field per unit cell;
  - an external potential written as 1D Laplace, with BC dφ/dx = E on the cell faces;
  - an electrode sign convention (cathode opposite the reservoir).

  Electrode spacing, material, area, overpotential and current density are not in the main text;
  the Supporting Information (power estimate) is not held.
- **Hence:** translating E into a device voltage is an engineering extrapolation. V_bulk = E × 5 mm
  is the bulk drop *if* the field spans the wick; the electrode voltage adds interfacial drops.
- **Retained as a major uncertainty, not a cutoff:**
  - DC Faradaic products: up to 1.4 µL/h of gas per half-wick (1e-6 M) or ~90 µL/h (1e-4 M) at
    1e5 V/m, if all current splits water (upper bound). Plus pH fronts that change ζ.
  - Field effect through the oxide.
  - Advective enrichment at the dryout-critical location.

## Parameter provenance (full table: `results/electroosmosis/parameter_provenance.csv`)

| Quantity | Provenance |
| --- | --- |
| d, l, h, L, W, θ_r (our wick) | PRIMARY_SOURCE (vaartstra2019) |
| K = 2.9698e-11 m² | COMSOL_EXTRACTED |
| F | MODEL_CALCULATED (Maxwell 2D, numerically verified) |
| Liquid properties | AUTHORITATIVE_DATABASE (IAPWS water; dilute-KCl approximation per pujahari2022) |
| KCl c, ζ, λ_D, ε_r, σ, E range, their property set | PRIMARY_SOURCE (pujahari2022, "stated in paper") |
| Fig. 6 curves | PROJECT_DERIVED (digitized) |
| Field application length 5 mm | ASSUMED (engineering extrapolation) |
| Λ°_KCl(100 °C), pK_w(100 °C), oxide ε_r, K band | ASSUMED (from memory; consistency checks only) |
| REQUIRES_SOURCE (19) | device ζ, ζ(T), ζ(c); electrolyte; σ; ε_r(T); K_w(T); Λ_KCl(T); KCl solubility; applied voltage; electrode spacing, material, reactions, area and overpotential; oxide thickness and substrate potential; breakdown field; ion diffusivities; interface ζ |

## Recommended poster figures

1. `fig2_Q_vs_pressure_gradient.png`: passive vs assisted Q through the same wick.
2. `fig9_pujahari_reproduction.png`: the published model reproduced factor by factor, with the
   disagreement shown.
3. `fig5_regime_map.png`: universal Pi map.
4. `fig8_preliminary_dryout.png`, labelled "idealized/conditional", with the B-band numbers.

## Single most important next calculation

Status: the 2D Python version is done (above). The remaining step is the 3D COMSOL version
specified in `docs/eo_field_treatment_test.md`, which adds the floor and finite height.

**A unit-cell EO-flux calculation with the two field treatments, on identical geometry.** Use the
same quarter cell, flat top, zero ΔP and thin-EDL slip (or PB body force) on pillars and floor.
Compare (a) the field solved with insulating walls against (b) a uniform imposed E_x. Run it for
the pujahari2022 cells (d = 5; h = 25 at l = 10 and 45; h = 5 at l = 8 and 18 µm) and for our
wick.

This yields R = U_b/U_a directly. It decides whether the ~0.6× difference is the field treatment,
and therefore which EO magnitude (Pi differs by 1.5–1.9×) should feed every poster number. It is
cheap in COMSOL (Electric Currents + Creeping Flow with the electroosmotic wall velocity) and needs
no new physical inputs.

Second, after that: a quasi-steady 1D calculation along the wick that couples the advective KCl
enrichment, current conservation with σ(c), and a parameterized ζ(c). Not the full transient
Nernst-Planck model.
