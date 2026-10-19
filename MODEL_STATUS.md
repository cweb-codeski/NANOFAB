# MODEL_STATUS — Passive vs Electroosmotically Assisted Micropillar Wick

Date: 2026-10-18. Scope: first pass of a reduced-order passive-vs-EO liquid-transport framework.
Equations and assumptions: `docs/electroosmosis_model.md`. Regenerate everything with
`python analysis/eo_assisted_wick.py` (≈2 min) → `results/electroosmosis/` (gitignored).
Tests: `pytest` → **196 passed** (163 pre-existing + 33 new in `tests/test_electroosmosis.py`).

## Bottom line

1. **The transport comparison is defensible as a model.** Under thin-EDL, uniform-zeta conditions
   the total flow is exactly Q_total = K A (-dp/dx)/mu + (mu_EO E / F) A. K is the same COMSOL
   value in both cases; F (formation factor) is geometric and verified to 0.008%.
2. **Hydrodynamically, EO can matter.** With pujahari2022's inputs, the baseline wick (d/l/h =
   10/30/25 µm, L = 5 mm) gains +1.6% to +57% in capillary-limited dryout heat flux over
   E = 5e3–1e5 V/m. It is never dominant there (Pi ≤ 0.64 < 1).
3. **The gain is set by voltage, not field.** Pi = V ε|ζ| / (F K Pc,max), independent of wick
   length. A 5% gain needs ≈ 40 V (−200 mV) or ≈ 80 V (−100 mV) across the wick. No
   geometrically similar pillar wick with a valid thin EDL reaches 5% below ~1.2 V.
4. **The binding constraints are electrochemical, not hydraulic.** At the fields that give useful
   gains, the model's own diagnostics show the following:
   - electrolysis voltage is exceeded by 20–400×;
   - field-effect charge through a thin oxide exceeds the native charge;
   - Dukhin number ≥ 0.6 for the −200 mV case;
   - salt accumulates without bound at the evaporating end of a sealed wick.

   A dryout improvement is therefore **not** claimable for a real device yet.

## What was verified

| Check | Result |
| --- | --- |
| COMSOL K reproduced from the raw export | 2.9698e-11 m² (needs `--dp 1 --mu 2.816e-4`; the export has no Pa column) |
| K mesh convergence | finest-pair change 0.10%; three-mesh extrapolation +0.13% (observed order 1.77) |
| K vs analytical (Sangani-Acrivos + Brinkman) | −1.33% |
| COMSOL viscosity vs IAPWS at 100 °C | 2.816e-4 vs 2.8158e-4 Pa s |
| Formation factor F, FV Laplace solve vs Maxwell 2D | 1.19131 vs 1.19122 (0.008%, first-order convergence 128/256/512) |
| Passive Darcy forced flow | linear, Q(0) = 0, matches K A G/mu to 1e-14 |
| Passive dryout regression | analytical K → 100.6 W/cm² (repo legacy value); COMSOL K → 99.28 W/cm² |
| 1D closed form vs existing 2D whole-wick Darcy solver | max difference < 1e-5 of the peak Pc |
| 1D closed form vs independent FV solve (with EO) | 2nd-order convergence; dryout bisection within 1e-5 (aiding) / 2e-3 (opposing, wall maximum) |
| Mass conservation (FV) | feed inflow = evaporation to 1e-10; zero total flux at the closed end |
| Q_EO → 0 as E → 0; sign | ζ < 0 drives flow along E; reversing E or ζ reverses it |
| Monotonicity | U_EO and q''_dry,B strictly increase with an aiding E |
| Identical geometry | Cases share one `Wick` object; any difference raises |
| Onsager / open-circuit reduction | formula for β checked against the 2×2 Onsager system |
| Gouy-Chapman flux-displacement length | matches quadrature to 1e-8; → λ_D at small ζ |
| Debye length | reproduces pujahari2022's 0.03 / 0.3 µm |
| Units | new units (mV, µS/cm, mol/L, kJ/kg) convert correctly; mobility checked dimensionally |

## What was implemented

- `src/nanofab/electroosmosis/physics.py`: HS mobility; thermal voltage; Debye length;
  Gouy-Chapman displacement; thin-EDL flux factor; Dukhin lower bound; Grahame charge;
  field-effect ratio; Joule heating; streaming/conduction ratio; electroviscous coefficient.
- `.../formation_factor.py`: Maxwell 2D, FV Laplace solver with Richardson extrapolation, and
  extraction from a COMSOL current.
- `.../wick1d.py`: `Wick`/`Liquid`/`Forcing`/`Case`; forced flow; apparent K; equivalent EO
  pressure; Pi; evaporating half-wick Pc(x); closed-form dryout gain G(Pi); over-pumping
  (flooding) limit; independent FV solver and dryout bisection; Reynolds numbers.
- `.../inputs.py` + `parameters/electroosmosis/eo_inputs.yaml`: provenance-checked loader. Device
  inputs live in a `requires_source` section (`value: null`), and `require()` blocks on them.
- `analysis/eo_assisted_wick.py`: benchmarks, sweeps, sensitivity, falsification maps, the
  pujahari2022 cross-check, validity diagnostics, 11 figures, 10 CSVs, `summary.json`.
- Registry and docs: pujahari2022 upgraded to read/verified, with audit findings recorded. Added
  mala1997, dutta_beskok2001 (from reference lists), and cummings2000, kirby2010, rayleigh1892,
  lyklema_fics (from memory, marked unverified). Added `docs/electroosmosis_model.md` and new
  equation/derived-result entries in `config/model_provenance.yaml`. Updated README,
  `docs/equations.md` and `docs/whole_wick_assumptions.md` (D1 note). `params.py` gained units only.
- Not modified: `src/nanofab/device/`, the 2D whole-wick solver, electrowetting, COMSOL import.

## Quantitative results (baseline wick, IAPWS water at 100 °C; MODEL_CALCULATED)

Pc,max = 2177 Pa. U_cap = K Pc,max/(mu L) = 0.0459 m/s. Passive dryout q''_dry,A = 99.3 W/cm².
vaartstra2019 reports 82 W/cm² (PRIMARY_SOURCE); the gap is the known flat-meniscus limitation.

| ζ (pujahari2022) | E [V/m] | V over 5 mm | Pi | G = q''_B/q''_A | q''_dry,B [W/cm²] | flooding below [W/cm²] |
| --- | --- | --- | --- | --- | --- | --- |
| −100 mV | 5e3 | 25 V | 0.016 | 1.016 | 100.9 | 1.6 |
| −100 mV | 1e4 | 50 V | 0.032 | 1.032 | 102.4 | 3.2 |
| −100 mV | 1e5 | 500 V | 0.319 | 1.299 | 129.0 | 31.6 |
| −200 mV | 5e3 | 25 V | 0.032 | 1.032 | 102.4 | 3.2 |
| −200 mV | 1e4 | 50 V | 0.064 | 1.063 | 105.5 | 6.3 |
| −200 mV | 1e5 | 500 V | 0.637 | 1.573 | 156.1 | 63.2 |

- **Velocities at 1e5 V/m, −200 mV:** HS slip 34.8 mm/s; superficial U_EO 29.2 mm/s; pore-average
  32.0 mm/s. For comparison, the passive feed velocity at dryout is 91.8 mm/s.
- **Sensitivity at −200 mV, 1e5 V/m.** These are one-at-a-time bounding intervals, not
  probabilities:

  | Parameter varied | Range of G |
  | --- | --- |
  | ζ pair | 1.30–1.57 |
  | K × 0.63–1.0 | 1.57–1.87 |
  | θ_r 0–40° | 1.55–1.71 |
  | L 1–10 mm at fixed E | 1.12–2.08 |
  | pitch 15–60 µm | 1.55–1.86 (minimum near the passive-optimal pitch) |
  | h 5–50 µm | 1.49–3.21 |
  | F | negligible |

  The corner bounds over ζ, K and θ_r give G = 1.29–2.08.
- **Regimes.** Negligible means Pi < 0.05; dominant means Pi > 1.
  - Negligible: ζ = −100 mV at E ≤ 1.6e4 V/m.
  - Useful: the rest of the literature window at the baseline geometry.
  - Dominant (Pi > 1 at 1e5 V/m): only for short pillars (h ≲ 10 µm at −200 mV, ≲ 6 µm at
    −100 mV) or long wicks (L ≳ 8 mm at −200 mV).
  - The field needed for a 5% gain ranges from 1.8e3 V/m (h = 5 µm, −200 mV) to 7.9e4 V/m
    (L = 1 mm, −100 mV).
- **Cross-check vs pujahari2022.** This compares two models; it is not validation. It is at
  E = 1e5 V/m, d = 5 µm, using their property set:
  - h = 25 µm, l = 20–45 µm: our G is 1.40–1.47 / 1.76–1.89 for −100 / −200 mV, against their
    stated 1.40–1.45 / 1.74–1.80. Agreement with no fitted parameters.
  - h = 25 µm, l = 10–15 µm: we predict higher gains (up to 2.31), while they report that gains
    *decrease* at small pitch. **Unresolved**; the causes have not been isolated.
  - h = 5 µm: we fall inside their stated ranges (1.88–2.16 vs 1.9–2.6; 2.65–3.16 vs 1.7–4.18),
    trending low at large pitch. That direction is consistent with flat-meniscus K overestimating
    K when the meniscus sag is 19–27% of h.
- **Diagnostics at the predicted Case B dryout point** (`validity_diagnostics.csv`):

  | Diagnostic | Value | Verdict |
  | --- | --- | --- |
  | Re_K at the feed | 1.7–2.7 (Re_gap 7–11), Case A included | near Darcy-Forchheimer onset |
  | Finite-EDL flux error | 0.3% (−100 mV) / 2.0% (−200 mV) | OK |
  | Dukhin lower bound | 0.009 (−100 mV) / 0.61 (−200 mV) | the −200 mV case violates the similitude |
  | Streaming/conduction current | 0.39–0.78 at 5e3 V/m; 0.02–0.05 at 1e5 V/m | strong coupling at low E |
  | Joule heating | ≤ 1.3e-5 of q'' | negligible |
  | Open- vs short-circuit β | ≤ 1.2% | negligible |
  | Field-effect charge / native charge | 0.09–3.5 (2 µm oxide); 1.7–70 (0.1 µm oxide) | uniform-ζ invalid unless the substrate is isolated |
  | Salt Péclet number | ≥ 4.7e3 even at an extreme D = 1e-7 m²/s | sealed-system enrichment at the closed end ~10^1000 |

## What remains provisional

- **Everything electrokinetic.** ζ, ε_r, σ, concentration and E come from a simulation paper,
  not from our surface. ζ was measured by mala1997 at unchecked conditions.
- **The assisted dryout (Fig. 8)** is conditional on EO1–EO12. It is a hydraulic upper estimate,
  not a device prediction.
- **Passive dryout** is the flat-meniscus estimate: +21% vs vaartstra2019. That gap carries into
  absolute values; the gain ratio G is more robust than absolute q''.
- **Off-baseline geometries** use analytical K and Maxwell F.
- **Sources cited from memory and not verified:** the HS relation, the similitude, the Dukhin
  formula, and the 1.23 V decomposition voltage (used only as a reference line).

## Parameter provenance (full table: `results/electroosmosis/parameter_provenance.csv`)

| Quantity | Value | Provenance |
| --- | --- | --- |
| d, l, h, L, W, θ_r | 10, 30, 25 µm; 5, 10 mm; 15° | PRIMARY_SOURCE (vaartstra2019) |
| K | 2.9698e-11 m² | COMSOL_EXTRACTED (Stage 1, flat top) |
| F | 1.1912 | MODEL_CALCULATED (Maxwell 2D, numerically verified) |
| Water properties at 100 °C | CoolProp | AUTHORITATIVE_DATABASE (IAPWS) |
| ζ, c, λ_D, ε_r, σ, E range | −100/−200 mV, 1e-4/1e-6 M, 0.03/0.3 µm, 55.4, 0.79 µS/cm, 5e3–1e5 V/m | PRIMARY_SOURCE (pujahari2022; "stated in paper") |
| Electrode spacing = 5 mm | for V = E·L | PROJECT_DERIVED (pujahari2022) |
| K band 0.63–1.0; oxide ε_r 3.9 | | ASSUMED |
| Sweep axes (mobility, L, pitch, h, t_ox, q'') | | USER_PROVIDED (not property values) |
| **REQUIRES_SOURCE** (12 inputs) | device ζ (at T, on our oxide); electrolyte and ionic strength; conductivity; ε_r at T (IAPWS); applied voltage; electrode spacing; electrode material and gas management; oxide thickness and substrate potential; breakdown field; ion diffusivities; liquid–vapor interface ζ; dζ/dT | — |

## Failed approaches and corrections made during this pass

- The first FV right-hand side had the wrong sign. Caught by the closed-form comparison, fixed
  before any results were produced.
- Formation-factor extrapolation on 64/128/256 grids was not asymptotic (observed order 0.54,
  −0.14% error). Switched to 128/256/512 (order 1.06, 0.008%).
- My stated large-Pi asymptote G → Pi was wrong; the exact limit is G → Pi/2. The test caught it.
- The endpoint-based tornado hid the non-monotonic pitch effect; replaced with dense-sweep min/max.
- The first regime map placed geometry markers at arbitrary y-values, visually implying
  "negligible". Replaced with each wick's actual U_EO at the strongest and weakest literature
  conditions.
- `analysis/permeability_from_comsol.py` fails on the curated export without `--dp/--mu`. Left
  unchanged; the command is documented in the README.

## Recommended poster figures (`results/electroosmosis/`)

1. `fig2_Q_vs_pressure_gradient.png`: passive vs assisted Q through the same wick. The EO part is
   the zero-gradient intercept, and the slope (K) is unchanged.
2. `fig5_regime_map.png`: the single falsification figure. G depends only on
   Pi = U_EO/U_cap; negligible, useful and dominant regions; literature band; geometry sweeps.
3. `fig11_voltage_scale.png`: the voltage needed for a meaningful gain vs geometry, and why
   down-scaling doesn't escape electrolysis while the thin-EDL condition holds.
4. Either `fig8_preliminary_dryout.png` (marked PRELIMINARY: operating window between over-pumping
   and dryout) **or** `fig10_validity_diagnostics.png`.
   Supporting: `fig1` (passive benchmark), `fig3`/`fig4` (velocities; K_app is not a property),
   `fig9` (cross-check).

## Is a reduced-order dryout prediction defensible?

- **Passive:** yes, as the flat-meniscus capillary-limit estimate (99.3 W/cm²), clearly labelled
  as +21% vs vaartstra2019.
- **Assisted:** only as a conditional hydraulic estimate. The equations do show EO lowering the
  pressure gradient needed to deliver the feed flow, so dryout rises by G(Pi). But in a sealed
  evaporating wick the uniform-electrolyte assumption (EO5) cannot hold, and the required voltages
  exceed electrolysis (EO6). An improved-dryout claim for a device is not supported.

## Exact next modeling step

Add 1D ion conservation to the evaporating half-wick: Nernst-Planck for K⁺/Cl⁻ with advection
by U(x), electromigration in E(x), and zero ion flux into the vapor. Couple it to charge
conservation, ∂x(σ(c) E + streaming current) = 0 at fixed electrode voltage. Run it transiently
from a uniform fill.

Output: the time for closed-end enrichment to change λ_D, σ and ζ by a given factor, as a
function of ion inventory. This decides whether *any* quasi-steady EO operating window exists.
It needs one new REQUIRES_SOURCE relation, ζ(c); parameterize it and sweep.

The cheapest secondary step is the COMSOL hand-off in `docs/electroosmosis_model.md`: Electric
Currents → F, EO-slip verification, and a Laminar-Flow inertia check on the same Stage 1 cell.
