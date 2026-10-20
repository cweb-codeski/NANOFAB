# Field-Treatment Test of the pujahari2022 EO-Flux Discrepancy (2D Python)

Code: `src/nanofab/electroosmosis/planform2d.py`, runner `analysis/eo_field_treatment_test.py`
→ `results/electroosmosis/field_treatment/` (`R_field_by_geometry.csv`, `convergence.csv`,
`verification_checks.csv`, `tables.md`, `fig_R_vs_pitch.png`, `fig_baseline_cell_fields.png`).
Tests: `tests/test_planform2d.py`. Dependencies: numpy/scipy only.

## Question and the direction of the ratio

The reduced-order model gives U_EO = mu_EO E / F. pujahari2022's Fig. 6 implies an EO flux
R_inferred = 0.52–0.69 times that. The hypothesis under test: pujahari2022 imposes a uniform field,
while a field solved around insulating pillars gives a different net flux.

Our reduced-order model **is** the solved-field result: by similitude, a field solved with insulating
walls gives exactly mu_EO E/F. A uniform field drives only the wall-tangential projection
E0 sin θ on a pillar sidewall, whereas the solved field wraps around the pillar (≈2 E0 sin θ). So the
uniform field gives the **smaller** flux. The quantity comparable to R_inferred is
Q_uniform/Q_solved = **1/R_field** in the definition R_field = Q_solved/Q_uniform. Both are reported.

## Model

- **Cases:** A, uniform E = E0 x̂ in the liquid; B, ∇²φ = 0 with n·∇φ = 0 on the pillar, equipotential
  faces midway between pillar columns (exact for a square array), mirror faces between rows, and mean
  gradient −E0.
- **Common to both:** the same HS slip u_t = mu_EO E_t on the pillar wall, the same ζ, ε, μ, geometry
  and mean field. Results are in units of mu_EO E0, so R is purely geometric.
- **Flow:** Stokes, no-penetration plus prescribed slip on the pillar (exact, mode by mode), zero
  mean pressure gradient (zero disturbance pressure on the fore-aft symmetry face), v = 0 and zero
  vorticity on the mirror face.
- **Method:** boundary collocation with odd-harmonic multipole series (Laplace) and the four
  biharmonic families per mode (Stokes). Spectral convergence.
- **Limitations (by construction):**
  - z-invariant: the pillars are infinitely tall and there is no floor, so **pillar height does
    not enter**.
  - No curved meniscus, no finite EDL, no PB body force.
  - The 2D result is the exact h → ∞ limit of the 3D flat-top problem, not a model of
    pujahari2022's curved-meniscus cells.

## Verification (baseline cell d/l = 1/3, plus d/l = 0.11, 0.4, 0.625)

| Check | Result |
| --- | --- |
| Current conservation, 4 sections | spread 2e-13 relative |
| Mean imposed gradient | potential drop across the cell = E0 l (1 ± 1e-13) |
| Periodicity / symmetry faces | φ + E0 x periodic to 5e-12; cell-face BC residuals ≤ 1e-10 (24 modes) |
| F vs independent FV Laplace solver | 7e-5 relative; vs Maxwell 2D 2e-6 (φ = 0.087) |
| Stokes solver vs Sangani–Acrivos K (pressure-driven, same code) | 0.05% (φ = 0.087), 1.3% (φ = 0.196); 14% at φ = 0.31, where the truncated series degrades |
| Similitude: Case B flux = 1/F | 3e-13 |
| Mass conservation, 4 sections, both cases | 9e-13 relative |
| Uniform-field slip imposed exactly | wall-slip residual 1e-12 |
| Reversing E | flux exactly antisymmetric (both cases) |
| E → 0 | Q_EO = 0 exactly |
| Resolution study | modes 8/12/16/24/32 and 48/96/192 collocation points: ratio converged to < 1e-6 by 16 modes |

**Vanishing obstruction (d/l → 0):** 1/F → 1, but Q_uniform/Q_solved → **0.5, not 1**. This is
correct physics, not a numerical defect. In the planform model the pillar wall is the only charged
surface. A uniform field supplies E0 sin θ there against 2 E0 sin θ for the solved field around any
isolated insulating cylinder, so by linearity the flux ratio stays ½ however small the pillar is.

The limit in which the ratio does go to 1 needs the floor and is exact. When h → 0, floor slip
dominates. The depth-averaged flow is then U = mu_EO E0 x̂ − (h²/3μ)∇p with no penetration at the
pillars, which is the same potential problem as the solved field, so Q_uniform = Q_solved = mu_EO E0/F.

## Results

| geometry (d/l/h, µm) | Q_uniform/Q_solved, 2D (R_field) | height bridge (heuristic) | R inferred (pujahari2022) |
| --- | --- | --- | --- |
| 5/10/25 | 0.592 (1.69) | 0.62–0.63 | 0.521 |
| 5/12.5/25 | 0.561 (1.78) | 0.60–0.61 | 0.558 |
| 5/15/25 | 0.543 (1.84) | 0.60–0.61 | 0.549 |
| 5/20/25 | 0.524 (1.91) | 0.62–0.63 | 0.534 |
| 5/25/25 | 0.516 (1.94) | 0.65–0.66 | 0.545 |
| 5/30/25 | 0.511 (1.96) | 0.69–0.70 | 0.579 |
| 5/35/25 | 0.508 (1.97) | 0.73 | 0.616 |
| 5/45/25 | 0.505 (1.98) | 0.80 | 0.657 |
| 5/8/5 | 0.631 (1.59) | 0.69–0.74 | 0.610 |
| 5/10/5 | 0.592 (1.69) | 0.71–0.76 | 0.666 |
| 5/14/5 | 0.549 (1.82) | 0.79–0.84 | 0.674 |
| 5/18/5 | 0.530 (1.89) | 0.86–0.91 | 0.663 |
| 10/30/25 (NANOFAB baseline) | 0.543 (1.84) | 0.66–0.68 | – |
| 20/50/40 (mask r = 10 µm, user-provided) | 0.561 (1.78) | 0.66–0.69 | – |

The full table is in `tables.md`. The **height bridge is a heuristic, not a 3D solution.** It
exploits two facts: Case B is exactly z-invariant, and the floor slip of Case A reaches a depth
fraction f = tanh(h/√K_s)/(h/√K_s). It is evaluated with two closures, exact-endpoint and
Brinkman-shear. It equals the 2D value as h → ∞ and 1 as h → 0.

## Interpretation

1. **Does field distortion plausibly explain the 0.52–0.69 factor?** Yes in magnitude.
   - With no fitting, the uniform-field treatment alone gives 0.50–0.63 (R_field = 1.59–1.98) in the
     tall-pillar limit. That is the same size as the inferred 0.52–0.69, and R_field matches the
     "1.5–1.9× larger" overprediction.
   - For h = 25 µm at pitch 12.5–25 µm, where pillar walls dominate, 2D and inferred agree within
     0.003–0.03.
   - 23 of 26 inferred points lie between the exact tall-pillar limit and the heuristic floor bridge.
     The three exceptions are the densest cells (d/l ≥ 0.4 and d/l = 0.625), 0.003–0.07 below the 2D
     limit.
2. **Does R depend strongly on pitch?**
   - In 2D, weakly: 0.505 → 0.63 as d/l rises from 0.11 to 0.63, with a slight decrease with pitch.
   - Once the floor is included, the dependence on pitch and height becomes strong through f, which
     is set by h/√K_s.
3. **Does the trend agree with the digitized data?** Partly.
   - h = 25 µm: their ratio rises from ≈ 0.53 to 0.66 with pitch. The 2D limit falls and the
     heuristic bridge rises, but the bridge overshoots (0.80 at l = 45 µm). The data sit between the
     two, consistent with a floor contribution smaller than the Brinkman bridge estimates.
   - h = 5 µm: their ratio is nearly flat at 0.61–0.69, while the bridge rises to 0.86–0.91. Their
     lower values are consistent with effects absent here, chiefly curved menisci that remove
     liquid near the floor (sag ≈ 20–27% of h) and reduce the floor-driven flux.
4. **What cannot be concluded:**
   - The quantitative 3D floor contribution, which needs a 3D solve.
   - Effects of the curved meniscus on the EO flux.
   - Finite-EDL or Poisson–Boltzmann body-force differences.
   - Whether pujahari2022 actually used a uniform field. Their text and the 1D form of Eq. 8 point
     that way, and this test shows the treatment produces the observed magnitude, but it is not proof.
   - Which treatment is physical. For DC fields and insulating (oxide-covered) pillars, n·J = 0
     forces the solved field, so our reduced-order mu_EO E/F is the physically consistent baseline.

**Outcome category:**
- **A in magnitude:** consistent with the observed ratio.
- **C for geometry:** the floor and meniscus introduce a pitch and height dependence that 2D cannot
  capture.
- **D for the quantitative 3D answer:** the minimum COMSOL calculation is below.

## COMSOL verification specification

**Geometry (one per case below):** full cell x, y ∈ [−l/2, l/2], z ∈ [0, h]. The cylindrical pillar
(diameter d) is centred and spans 0 ≤ z ≤ h, with a flat top at z = h. Remove the pillar from the
liquid domain.

**Geometries to run:**
- pujahari2022, d = 5 µm: h = 25 µm with l = 10, 15, 25, 45 µm; h = 5 µm with l = 8, 12, 18 µm.
- Baseline: d/l/h = 10/30/25 µm.
- Mask: d/l/h = 20/50/40 µm.

**Electric Currents (Case B only):**
- σ = 1 S/m (R is independent of σ).
- Electric insulation on the pillar wall, floor, top and y = ±l/2.
- Electric potential ΔV = E0·l on x = −l/2 and ground on x = +l/2, with E0 = 1e4 V/m (any value).

**Creeping Flow (both cases):**
- μ = 2.822e-4 Pa·s.
- Periodic flow condition between x = ±l/2 with pressure difference 0, plus a pressure point
  constraint.
- Symmetry on y = ±l/2. Slip (symmetry) on the top z = h.
- Pillar wall and floor: **Electroosmotic Velocity** wall condition with ε_r = 55.4, ζ = −0.1 V.
  - Case A: user-defined field (E0, 0, 0).
  - Case B: the Electric Currents field (ec.Ex, ec.Ey, ec.Ez).
- COMSOL projects the field onto the wall tangent in both cases; confirm this in the condition
  settings.

**Extra runs for cross-checking the Python model:**
- 2D-limit check: for the baseline and d/l = 0.5, replace the floor EO condition with slip and set
  h ≥ 10 l. Q_A/Q_B should reproduce the Python 2D values (0.543 and 0.592).
- Thin-film check: one run with h = 0.05 l and floor EO. Q_A/Q_B should be close to 1.

**Mesh:** three meshes (for example hmax = l/10, l/20, l/40, refined on the pillar wall). Report DOF.

**Exports, one row per geometry × case × mesh** (`case, floor_bc, d_um, l_um, h_um, hmax_um, dof,
Q_out, Q_mid, U_sup, I_out, F, Utilde`):

1. **Q_out:** ∫ u dA over x = +l/2 [m³/s].
2. **Q_mid:** ∫ u dA over the liquid part of x = 0. This is a mass check; it must equal Q_out.
3. **U_sup:** Q_out/(l h) [m/s]. Also report Ũ = U_sup / (mu_EO E0), with mu_EO = −ε_r ε0 ζ / μ.
4. **I_out (Case B):** ∫ ec.normJ over x = +l/2, from which F = σ (l h) ΔV / (I_out l). Expect
   1.191 for the baseline.
5. **Similitude check (Case B):** Ũ·F, which should equal 1.
6. **Wall field (Case B):** tangential field magnitude on the pillar wall at z = h/2 for θ = 0–180°.
   Compare with `fig_baseline_cell_fields.png` panel (b).
7. **Velocity profiles:** u(z) at (0, l/2) and (l/2, l/2) for both cases. The floor-layer thickness
   tests the height bridge.

**Primary result:** Q_A/Q_B per geometry, compared with the Python 2D value, the heuristic bridge and
R_inferred in `R_field_by_geometry.csv`.
