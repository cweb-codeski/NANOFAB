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
| Liquid (IAPWS water, 100 °C) | shared `Liquid` object | the **same** |
| Electrical forcing | mu_EO E = 0 | mu_EO, E |

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

**Voltage, not field, controls the gain.** With electrodes at the feed and at the closed end,
V = E L and

```text
Pi = V eps |zeta| / (F K Pc,max)  →  V_cap = F K Pc,max / (eps |zeta|),   Pi = V / V_cap
```

independent of L. For geometrically similar wicks K Pc,max ∝ scale, so V_cap ∝ scale while the
passive dryout flux ∝ scale² (at fixed L).

## Assumption register

| # | Assumption | Consequence if wrong | Diagnostic in the code |
| --- | --- | --- | --- |
| EO1 | Thin EDL: λ_D << pore half-gap and h | EO flux overestimated by ~δ/a | `thin_edl_flux_factor` (Gouy-Chapman displacement, exact for thin EDL) |
| EO2 | Uniform zeta on pillars **and** floor; uncharged, shear-free meniscus | similitude fails; U_EO ≠ mu_EO E/F | none (needs data); `requires_source.liquid_vapor_interface_zeta` |
| EO3 | No surface conduction (Du << 1) | field distorted near walls; similitude fails | `dukhin_number_lower_bound` (Bikerman, m = 0) |
| EO4 | Field set by electrodes only; streaming current negligible | E(x) varies with the flow | `streaming_to_conduction_ratio` |
| EO5 | Uniform electrolyte (c, σ, zeta) along the wick | **fails in a sealed evaporating wick**: ions cannot leave with the vapor | salt Péclet number; zero-flux enrichment exp(Pe/2) |
| EO6 | No electrode reactions | DC voltage >> ~1.2 V electrolyzes water → H₂/O₂ non-condensables in a vapor chamber | applied voltage V = E × spacing |
| EO7 | Oxide is only an insulator (no field effect) | a conducting substrate under thin oxide induces charge that changes zeta along x | `field_effect_charge_ratio` (t_ox swept, REQUIRES_SOURCE) |
| EO8 | Flat meniscus at the pillar top (K, F flat-top) | receding meniscus lowers K and changes F | K × 0.63–1.0 band (ASSUMED) |
| EO9 | Creeping flow (Darcy) | Re_K at the feed is 1.7–2.7 → Darcy-Forchheimer onset | `reynolds_numbers` |
| EO10 | Joule and viscous heating negligible | — | `joule_heat_flux` (≤ 1.3e-5 of q'') |
| EO11 | Open- vs short-circuit electrodes do not change Case A | electroviscous reduction β | `electroviscous_coefficient` (≤ 1.2%) |
| EO12 | No electrowetting cross-coupling | V/2 across the oxide is far above EWOD actuation voltages (whole_wick model: ~11 V for 1 µm SiO₂) | flagged only |

## Where the electrokinetic numbers come from

`pujahari2022` (held, read): KCl-water at 1e-4 / 1e-6 M with zeta = -100 / -200 mV (taken from
`mala1997`, not held), λ_D = 0.03 / 0.3 µm, ε_r = 55.4 and σ = 0.79 µS/cm at 100 °C, E = 5e3–1e5
V/m along a 5 mm wick. These are inputs of a simulation paper, not measurements on our surface.
The model is run as a **literature-bounded scenario** plus exploration sweeps; device-specific
inputs are `requires_source` entries with `value: null`, and `EOInputs.require()` blocks on them.

Audit findings on `pujahari2022` (recorded in `config/sources.yaml`):
1. Sec. 2.1.2 says the voltages are far below electrolysis, but E = 5e3–1e5 V/m over the 5 mm
   wick is 25–500 V.
2. One conductivity is used for both concentrations.
3. No ion conservation: the evaporating end accumulates salt.
4. No field-effect change of zeta through the thin oxide.

## COMSOL hand-off (same quarter cell as Stage 1)

- **F:** Electric Currents, σ = 1 S/m in the liquid, insulating pillars/floor/top, ΔV = 1 V over
  l/2. `F = σ A ΔV / (I L)` with A = (l/2) h, L = l/2 (`formation_factor_from_current`). Expect 1.191.
- **EO verification:** Creeping Flow with the Electroosmotic Velocity wall condition on pillars and
  floor, slip top, zero ΔP, field from the Electric Currents study. Superficial velocity should be
  mu_EO E / F to within the mesh error.
- **Inertia:** rerun Stage 1 with Laminar Flow at ΔP that gives U ≈ 0.1 m/s (≈ 14 Pa over l/2) and
  compare K with the creeping-flow value.
