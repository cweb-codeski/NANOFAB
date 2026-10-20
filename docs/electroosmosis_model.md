# Passive vs Electroosmotically Assisted Wick: Reduced-Order Model

Code: `src/nanofab/electroosmosis/` (physics, formation_factor, wick1d, inputs).
Inputs: `parameters/electroosmosis/eo_inputs.yaml`. Tests: `tests/test_electroosmosis.py`.
Runner: `python analysis/eo_assisted_wick.py` → `results/electroosmosis/` (figures, CSVs, `summary.json`).
Status and results: `MODEL_STATUS.md`.

## Comparison being made

| | Case A (passive) | Case B (assisted) |
| --- | --- | --- |
| Geometry d, l, h, L, W | shared `Wick` object | the **same** object |
| Intrinsic permeability K | COMSOL Stage 1, 2.9698e-11 m² | the **same** K |
| Formation factor F | geometric, 1.1912 | the **same** F |
| Pc,max (V2019 Eq. 3) | 2177 Pa | the **same** |
| Liquid: dilute aqueous KCl, water-like properties (IAPWS, 100 °C) | shared `Liquid` object | the **same** |
| Electrical forcing | mu_EO E = 0 | mu_EO, imposed E |

Results are **idealized/conditional electroosmotic hydraulic predictions under an imposed electric
field**. They are not device predictions; the electrochemistry that would supply the field is not
modeled.

`wick1d.paired_cases` builds both from one wick and one liquid object, and
`assert_identical_except_forcing` raises if anything else differs. K is never a function of E.

## Governing equations

Sign convention: x runs from the liquid feed (x = 0) toward the closed end (x = L). Positive
velocity = toward the closed end. E is the macroscopic field along +x.

**Pore scale (thin EDL).** In the electroneutral bulk outside a thin double layer:

```text
mu ∇²u = ∇p,   ∇·u = 0                           Stokes (Re << 1 assumed; see diagnostics)
∇²φ = 0,  E_local = -∇φ                            Ohmic conduction, uniform bulk conductivity
u = mu_EO E_local,t  on charged walls (pillars AND floor)     Helmholtz-Smoluchowski slip
mu_EO = -eps_r eps_0 zeta / mu                      [m²/(V s)], > 0 for zeta < 0
n·∇φ = 0 on insulating walls and on the meniscus;  meniscus shear-free, uncharged
```

**Superposition (exact under the assumptions).** The problem is linear, so
u = u_pressure (no-slip, imposed ∇p, E = 0) + u_EO (slip, ∇p = 0). With uniform zeta on every
wetted surface, u_EO = mu_EO E_local satisfies Stokes with constant pressure, the slip condition
and, for a flat meniscus with z-independent field, zero shear at the free surface (similitude,
`cummings2000`, not held). Volume-averaging (current conservation gives
⟨E_local⟩_superficial = E / F):

```text
U = (K/mu)(-dp/dx) + U_EO,        U_EO = mu_EO E / F
Q = U W h                          [m³/s] volumetric flow through the full cross-section W h
```

- U is the **superficial** (Darcy) velocity, per total cross-section including pillars. This
  matches how K was extracted from COMSOL (A = (l/2) h in the quarter cell).
- Pore-average (intrinsic) EO velocity: U_EO / (1 - φ). HS slip velocity: mu_EO E.
- F = σ_bulk / σ_eff is the electrical formation factor of the liquid-filled pore space, a
  geometric property like K. Maxwell 2D: F = (1 + φ)/(1 - φ) = 1.1912 at φ = 0.0873; the repo's
  FV Laplace solver gives 1.19131 (Richardson, 128/256/512 grids), a 0.008% difference.

A Darcy-Brinkman equation with a ρ_e E body force (`pujahari2022`'s approach) is needed only if
the EDL is not thin, zeta is non-uniform, or surface conduction matters. Under thin EDL and
uniform zeta, the HS-slip + formation-factor result is the exact outer limit of that formulation,
so the simpler model is used and the diagnostics test when it fails.

**Apparent (effective) hydraulic response.** K_app = mu U_total / (-dp/dx) = K + mu U_EO/(-dp/dx).
It is NOT a property: it diverges as dp/dx → 0 (Fig. 4b). pujahari2022 Sec. 3.3 calls the same
quantity "effective permeability". This repository reports K (geometric) and U_EO separately, and
uses the **equivalent EO pressure** ΔP_EO = mu U_EO L / K (the pressure drop that would drive U_EO
through the same K) as the comparable quantity.

## Forced flow (no evaporation)

```text
Q_pressure = K W h ΔP / (mu L),   Q_EO = (mu_EO E / F) W h,   Q_total = Q_pressure + Q_EO
Q_total / Q_pressure = 1 + ΔP_EO / ΔP
```

## Evaporating half-wick and dryout (PRELIMINARY)

Feed at x = 0 with Pc = p_v - p_l = 0; closed end at x = L blocks **all** flow (pressure-driven
and electroosmotic); uniform q'' fully converted to evaporation; flat-meniscus K.

```text
U(x) = U_0 (1 - x/L),  U_0 = q'' L / (rho h_fg h)            mass conservation
dPc/dx = (mu/K)(U(x) - U_EO)                                   momentum
Pc(x) = (mu/K)[U_0 (x - x²/2L) - U_EO x]
max Pc at x* = L (1 - U_EO/U_0):  Pc* = (mu U_0 L / 2K)(1 - U_EO/U_0)²     (0 <= U_EO <= U_0)
```

Dryout when Pc* = Pc,max. With **Pi = mu U_EO L / (K Pc,max) = U_EO / U_cap** and
U_cap = K Pc,max / (mu L):

```text
G(Pi) = q''_dry,B / q''_dry,A = (1 + Pi + sqrt(1 + 2 Pi)) / 2     Pi >= 0
G(Pi) = 1 + Pi                                                    Pi < 0 (opposing field)
G ≈ 1 + Pi (small Pi);  G → Pi/2 (EO-dominated: U_0 → U_EO)
```

Derivation is ours (PROJECT_DERIVED); it reduces exactly to `dryout.dryout_heat_flux` at Pi = 0.
Verified against an independent finite-volume solve with a dryout bisection (relative difference
< 1e-5 for aiding fields).

**New constraint introduced by the field: over-pumping.** Pc(L) = (mu U_0 L/2K)(1 - 2 U_EO/U_0)
is negative when U_0 < 2 U_EO, i.e. below q''_flood = 2 U_EO rho h_fg h / L. Then the liquid
pressure at the closed end exceeds the vapor pressure (convex meniscus, pillar-top flooding),
outside the validity of the flat-top K. The assisted wick therefore has an operating window
q''_flood < q'' < q''_dry,B, and the field must be throttled at low load. At the dryout point
the closed end is not flooded while Pi <= 4.

**The bulk potential drop, not the field, controls the gain.** If the field acts over the whole
feed length, the bulk (ohmic) liquid potential drop is V_bulk = E L and

```text
Pi = V_bulk eps |zeta| / (F K Pc,max)  →  V_cap = F K Pc,max / (eps |zeta|),   Pi = V_bulk / V_cap
```

independent of L. For geometrically similar wicks K Pc,max ∝ scale, so V_cap ∝ scale while the
passive dryout flux ∝ scale² (at fixed L). V_bulk is not the electrode-to-electrode voltage (see
"Electrochemical system").

## Electrochemical system (dilute aqueous KCl)

The pujahari2022 working fluid is a dilute KCl-water electrolyte (Table 1), with water-like
thermophysical properties because the solution is dilute (end of Sec. 2.1.4). This repository
uses IAPWS water properties on exactly that basis.

| Case | c_KCl | ζ (paper) | λ_D stated | λ_D computed, KCl only | λ_D incl. intrinsic H⁺/OH⁻ | σ in Table 1 | σ_KCl estimate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1e-6 M | −200 mV | 0.3 µm | 0.286 µm | ≈0.22 µm | 0.79 µS/cm | ≈0.4–0.5 µS/cm |
| 2 | 1e-4 M | −100 mV | 0.03 µm | 0.0286 µm | 0.0285 µm | 0.79 µS/cm | ≈40–48 µS/cm |

ε_r = 55.4 and T = 373 K. The intrinsic-ion and σ_KCl columns use ASSUMED inputs from memory
(pK_w ≈ 12.27 at 100 °C; Λ°_KCl ≈ 400–475 S cm²/mol), so they are consistency checks only.

- **Debye length:** the stated values reproduce from c (to rounding).
- **Intrinsic ions:** at 1e-6 M, the intrinsic water ions are comparable to the KCl and shorten
  λ_D by ~25%.
- **Conductivity:** Table 1's single 0.79 µS/cm is plausible for case 1 but ~50× low for case 2.
  It enters the paper only through Eq. 8, which doesn't affect the flow under a uniform field. So
  it changes neither their nor our hydraulics; it changes current, Joule heating and Faradaic
  rates.
- **Zeta:** ζ values come from mala1997 (not held). Their temperature dependence and the ζ(c)
  relation needed along an enriching wick are REQUIRES_SOURCE.

**Field vs voltage vs electrolysis.** These are distinct quantities:
- **Imposed field E:** in pujahari2022, a uniform x-field per unit cell. Eq. 8 is written in 1D,
  the cell faces have dφ/dx = E, and the text says "a constant electric field is considered to be
  acting along each unit cell to calculate the volumetric force".
- **Bulk liquid potential drop:** E × (length the field acts over). That length is not specified.
  The only electrode statement is that the negative electrode is on the side opposite the liquid
  reservoir (Sec. 2). The electrolysis sentence refers to "the considered range of electric field
  and wick length", and the Supporting Information (not held) contains the power estimate.
  **V_bulk = E × 5 mm is therefore an engineering extrapolation, not a reproduction of the paper.**
- **Electrode-to-electrode voltage:** the bulk drop plus the interfacial drop at each electrode.
- **Thermodynamic decomposition potential and electrode overpotentials:** these apply at the
  electrode interfaces. They depend on electrode material and current density, and they are not
  limits on the bulk drop.

The defensible statement is current-based. A steady DC current
I = (σ/F) E W h (1.7 µA per half-wick at 1e5 V/m with Table 1's σ; ~100 µA for case 2 with the
KCl-corrected σ) must cross each electrode interface. With inert electrodes in dilute KCl, that
requires a Faradaic reaction (water splitting or chloride oxidation) at a rate I/(nF), unless the
design uses capacitive, AC or membrane-isolated electrodes. pujahari2022 states that the applied
voltage or field "should not be high enough to result in electrolysis or water splitting"
(Sec. 3.5) and that the voltages "will be very less compared to the voltage that causes
electrolysis" (Sec. 2.1.2). How that is achieved depends on electrode geometry and electrochemistry
the paper does not specify. Electrochemical feasibility is a major open uncertainty, not a
validity cutoff on the hydraulic model.

## Salt enrichment (dilute KCl, solvent evaporates)

Ions don't leave with the vapor and Pe = U L / D ≫ 1. Salt-flux conservation then gives a
quasi-steady **advective enrichment** profile, established within about one liquid transit time:

```text
c(x) / c_feed = Q_0 / Q(x) = 1 / (1 - x/L)        (same in Case A and Case B at the same q'')
```

It holds where L − x ≫ δ = sqrt(D ε L / U_0) ≈ 6–22 µm (D = 1e-9–1e-8 m²/s). Inside that end layer,
salt accumulates with time; with a constant-concentration reservoir there is no steady state
there until precipitation (pujahari2022 Sec. 3.5 mentions KCl precipitation only after a complete
loss of liquid). The Case B dryout-critical location x* = L (1 − U_EO/U_0) lies 0.04–1.0 mm from
the closed end, where c/c_feed ≈ 5–130. The electrolyte at the critical location therefore is
**not** the feed composition: λ_D falls by sqrt of that factor, σ rises, and ζ changes by an
unknown amount.

The global balance and Pe ≫ 1 are enough for this profile and for "no steady state in the end
layer". They are not enough for ζ(x), E(x) (current conservation with σ(c)), electrode-generated
pH fronts, or precipitation timing. Those need Nernst-Planck with electroneutrality, electrode
reactions and a sourced ζ(c).

## Reproduction of pujahari2022 Fig. 6

Their Fig. 6 was digitized from the held PDF (`analysis/digitize_pujahari2022_fig6.py` →
`data/literature/pujahari2022_fig6_digitized.csv`, PROJECT_DERIVED). The passive curve is drawn
twice in the figure, and the two copies agree within 0.2% at interior pitches. Comparing factor by
factor at the same d, h, l, L, ζ and E (`results/electroosmosis/fig9_pujahari_reproduction.png`):

1. **Field dependence (structure):** Pi inferred from their 1e5 V/m curve predicts their
   5e3–5e4 V/m curves through G(Pi) with Pi ∝ E. Median error 0.35%, maximum 5.4%. Their model
   behaves exactly as Darcy plus a field-linear EO flux in the 1D capillary limit.
2. **Passive K:** their passive dryout implies K/K_flat = 0.75–0.97 for h = 25 µm and 0.44–0.67
   for h = 5 µm, decreasing with pitch. This is consistent with meniscus curvature, which our
   flat-top K omits.
3. **EO flux magnitude:** their EO flux is R = 0.52–0.66 (h = 25 µm) and 0.60–0.69 (h = 5 µm)
   times mu_EO E / F. This is the real discrepancy.
4. **Hypothesis for item 3 (PROJECT_DERIVED, not confirmed):** a uniform imposed x-field drives
   slip mu_EO E sinθ on an insulating pillar sidewall, half the ≈2 mu_EO E sinθ of the field that
   satisfies n·E = 0 (exact factor ½ in the dilute 2D limit by linearity). The floor sees the full
   field but is screened by the pillars (Brinkman length sqrt(K_s)). This gives R ≈ (1 + φ) f +
   (1 − f)/2 = 0.54–0.80. It matches the inferred R at small pitch but overpredicts it at large
   pitch.
   - With no fitted parameters, it reproduces their 1e5 V/m gains within 9% (h = 25 µm) and 21%
     (h = 5 µm).
   - A single fitted R per height (0.58 / 0.67) reproduces them within 5%.

The earlier "smaller pitch gives smaller gain" disagreement was partly a misreading: their stated
gain ranges are the extremes of U-shaped curves (h = 25 µm), not ranges over every pitch.

**Independent test of item 4** (`docs/eo_field_treatment_test.md`): a resolved 2D planform
Laplace + Stokes calculation finds Q_uniform/Q_solved = 0.50–0.63 in the tall-pillar limit, with no
fitting. That is the observed magnitude; for h = 25 µm at l = 12.5–25 µm it agrees with R_inferred
within 0.03. Height dependence (floor and meniscus) needs the 3D COMSOL check specified there.

## Assumption register

| # | Assumption | Consequence if wrong | Diagnostic in the code |
| --- | --- | --- | --- |
| EO1 | Thin EDL: λ_D << pore half-gap and h | EO flux overestimated by ~δ/a | `thin_edl_flux_factor` (Gouy-Chapman displacement, exact for thin EDL) |
| EO2 | Uniform zeta on pillars **and** floor; uncharged, shear-free meniscus | similitude fails; U_EO ≠ mu_EO E/F | none (needs data); `requires_source.liquid_vapor_interface_zeta` |
| EO3 | No surface conduction (Du << 1) | field distorted near walls; similitude fails | `dukhin_number_lower_bound` (Bikerman, m = 0) |
| EO4 | Field set by electrodes only; streaming current negligible | E(x) varies with the flow | `streaming_to_conduction_ratio` |
| EO5 | Feed-composition electrolyte at the dryout-critical location | KCl enriches as 1/(1 − x/L): c/c_feed ≈ 5–130 at x*; λ_D, σ, ζ change there | `reproduction.advective_enrichment`, enrichment-layer width |
| EO6 | Imposed uniform field; electrode reactions not modeled | a steady DC current I must be carried by Faradaic reactions (gas, pH fronts) unless the electrodes are capacitive, AC or membrane-isolated | current and Faradaic upper-bound gas rate vs E |
| EO7 | Oxide is only an insulator (no field effect) | a conducting substrate under thin oxide induces charge that changes zeta along x | `field_effect_charge_ratio` (t_ox swept, REQUIRES_SOURCE) |
| EO8 | Flat meniscus at the pillar top (K, F flat-top) | receding meniscus lowers K and changes F | K × 0.63–1.0 band (ASSUMED) |
| EO9 | Creeping flow (Darcy) | Re_K at the feed is 1.7–2.7 → Darcy-Forchheimer onset | `reynolds_numbers` |
| EO10 | Joule and viscous heating negligible | — | `joule_heat_flux` (≤ 1.3e-5 of q'') |
| EO11 | Open- vs short-circuit electrodes do not change Case A | electroviscous reduction β | `electroviscous_coefficient` (≤ 1.2%) |
| EO12 | No electrowetting cross-coupling | V/2 across the oxide is far above EWOD actuation voltages (whole_wick model: ~11 V for 1 µm SiO₂) | flagged only |

## Where the electrokinetic numbers come from

`pujahari2022` (held, read, Fig. 6 digitized) supplies:
- dilute aqueous KCl at 1e-4 / 1e-6 M, with ζ = −100 / −200 mV (taken from `mala1997`, not held);
- λ_D = 0.03 / 0.3 µm;
- ε_r = 55.4 and σ = 0.79 µS/cm at 100 °C;
- an imposed field E = 5e3–1e5 V/m along the wick direction.

These are inputs of a simulation paper, not measurements on our surface. The model runs them as a
**literature-bounded scenario** plus exploration sweeps. Device-specific inputs are
`requires_source` entries with `value: null`, and `EOInputs.require()` blocks on them.

Audit notes on `pujahari2022` (also in `config/sources.yaml`):
1. The field is imposed as a uniform x-field in the body force, and electrodes are not modeled.
   V_bulk = E × 5 mm (25–500 V) is our extrapolation. If the field does span the wick, it is
   difficult to reconcile with the paper's statement that the voltages are far below the
   electrolysis voltage. Their electrode configuration, which isn't specified, would decide it.
2. One conductivity (0.79 µS/cm) is given for both concentrations. It is plausible for 1e-6 M and
   ~50× low for 1e-4 M, and it does not affect the hydraulics.
3. Ion conservation is not modeled: KCl enriches toward the evaporating end. The paper discusses
   precipitation only after a complete loss of liquid.
4. The oxide is treated only as an insulator, with no field-effect change of ζ.
5. A uniform imposed field gives about half the physically consistent EO slip on insulating
   pillar sidewalls (hypothesis; see the reproduction section).

## COMSOL hand-off (same quarter cell as Stage 1)

- **Priority: field treatment (resolves the R = 0.52–0.69 discrepancy).** Use the pujahari2022
  cells (d = 5 µm; h = 25 µm with l = 10 and 45 µm; h = 5 µm with l = 8 and 18 µm) and our
  baseline. Use a flat top, zero ΔP, and the same thin-EDL slip (or the PB body force) on pillars
  and floor. Run two cases: (a) field from Electric Currents with insulating pillars/floor/top, and
  (b) a uniform imposed E_x. The superficial velocities give R = U_b / U_a directly. Expect
  U_a = mu_EO E / F and R ≈ 0.5–0.8.

- **F:** Electric Currents, σ = 1 S/m in the liquid, insulating pillars/floor/top, ΔV = 1 V over
  l/2. `F = σ A ΔV / (I L)` with A = (l/2) h, L = l/2 (`formation_factor_from_current`). Expect 1.191.
- **EO verification:** Creeping Flow with the Electroosmotic Velocity wall condition on pillars and
  floor, slip top, zero ΔP, field from the Electric Currents study. Superficial velocity should be
  mu_EO E / F to within the mesh error.
- **Inertia:** rerun Stage 1 with Laminar Flow at ΔP that gives U ≈ 0.1 m/s (≈ 14 Pa over l/2) and
  compare K with the creeping-flow value.
