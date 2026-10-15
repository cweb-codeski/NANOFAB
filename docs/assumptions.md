# Assumptions Register

Each modeling assumption, why it was made, how it could fail, and how to test it. "Baseline
check" values come from `python analysis/baseline_estimates.py`, using the V2019 baseline
(d = 10, l = 30, h = 25 µm, L = 5 mm, water at 100 °C) at the project target of 100 W/cm²
(DESIGN_TARGET, `parameters/design_targets.yaml`). All "Baseline check" numbers are
MODEL_CALCULATED; "V2019" entries are PRIMARY_SOURCE statements. Provenance classes are defined in
`docs/SOURCE_AUTHORITY.md`. Rule 7 there: an assumption that materially affects performance
must be swept.

The device-scale model has its own register with regime diagnostics, failure consequences and
qualification status: `config/device_assumptions.yaml` (summarized in `docs/passive_device_model.md`).

| # | Assumption | Where used | Baseline check | Status / how to test |
| --- | --- | --- | --- | --- |
| A1 | Laminar flow; inertia small in the unit cell | Stage 1, 3 | Re_h ≈ 8 at the feed; Re ≈ 0.5 in the Stage 1 test cell | Use Laminar Flow as V2019 does. Compare against Creeping Flow once. |
| A2 | Flat meniscus for K (analytical, Stage 1) | analytical K, dryout | flat q''_dry is 1.23x the V2019 value | **Not established as a bound.** V2019 Sec. 2.2 says curvature lowers K, so the flat K is expected to exceed the curved-meniscus K (PROJECT_DERIVED expectation). Even so, the legacy constant-K, uniform-q dryout is not shown to bound the coupled V2019 result: that would also need nonuniform evaporation to lower the limit and the analytical K coefficients to be right, neither of which is verified. Treat it as an approximate analytical estimate. Stage 3 provides K(Pc). |
| A3 | Liquid heat transfer by conduction only | Stage 2, 4 | Pe_h (vertical) ≈ 0.07; Ja ≈ 0.02 | Valid (V2019 A.1). A large lateral Pe is harmless while Ja << 1. |
| A4 | Gravity negligible | all | Bo(pitch) ≈ 1e-4; ρgL/Pc ≈ 0.02 | OK for horizontal wicks. |
| A5 | Vapor not resolved; kinetic BC at the interface | Stage 4 | Ma ~ 1e-3 (V2019) | Valid, near-equilibrium. |
| A6 | Pure saturated vapor at T_v | Stage 4 | – | Matches V2019 and a vacuum-chamber experiment. **Not valid in air**, where evaporation becomes diffusion-limited. |
| A7 | Accommodation coefficient a = 0.052 | Stage 2, 4 | R_int ≈ 2.6 mm²K/W. Combined with V2019's h_e > 60 kW/m²K this gives R_int/R_total = R_int h_e > ~16%, a conditional **lower** bound on the kinetic share (valid only if our R_int is right; it uses an ASSUMED flat interface area and an unverified accommodation correction) | Fitted by V2019 to one data set. The nonunity-accommodation treatment (Meland & Ytrehus 2003) is a source gap. Sweep 0.03-1 for our own surfaces. |
| A8 | Receding contact angle 15° | dryout Δp | – | V2019 value for water on SiO₂. Measure on our surfaces; sweep 0-40°. |
| A9 | Meniscus pinned at pillar top edge | Stage 3 | – | Confirmed as the V2019 assumption (Sec. 2.1). |
| A10 | Insulated pillar tops | Stage 2-4 | – | Confirmed in V2019: natural convection is 2-3 orders of magnitude weaker. |
| A11 | Constant curvature within a cell | Stage 3 | Ca ≈ 5e-4 | Valid (V2019 A.1). |
| A12 | Symmetric unit cell (no lateral heat flow between cells) | all unit-cell models | – | V2019: 0.19% (uniform load), 1% (hotspot). Re-check for hotspots or novel geometries. |
| A13 | Uniform heat flux, uniform wick | analytical dryout | – | Device model (Phase 7) handles non-uniform loads. |
| A14 | Steady state | all | – | Transient only for power-step studies later. |
| A15 | Properties evaluated at T_v | all | – | Superheat is a few K, so OK. |
| A16 | Silicon k = 148 W/(m·K) | Stage 2-4, device substrate | substrate ΔT < 1 K at 1 µm; ≈ 0.55 K across t_s = 100 µm at 82 W/cm² | **ASSUMED, provisional.** A textbook value recalled for 300 K, not re-checked against the table, not verified at 373.15 K, and not a V2019 value (V2019 does not state its silicon conductivity; vaartstra2020 uses 153 at 20 °C). Silicon k falls with temperature, so the 373 K value is likely lower. Non-blocking for the device model; sweep 120-153. |
| A17 | Advancing contact angle 30° | not used in baseline | – | ASSUMED (no source). |
| A18 | Water kinetic diameter 0.265 nm | Knudsen estimate only | – | ASSUMED (citation unverified). Affects only the Kn sanity check. |
| A19 | T_e = mean temperature of the pillar-base plane | Stage 2, 4 | – | PROJECT_DERIVED reading of V2019 Fig. 2(b) ("surface (red line)"), supported by the thesis drawing (Fig. 4b, p. 20), where the red line runs along the substrate top. Text never names the surface. |

## Consequences already visible

- **Dryout.** With the V2019 geometry, the MODEL_CALCULATED flat-meniscus capillary limit
  (100.6 W/cm²) is +22.7% above the paper's 82 W/cm² (PRIMARY_SOURCE). Curvature-dependent
  permeability (Stage 3, 5) is the likely main cause (PROJECT_DERIVED hypothesis), but error in
  the analytical K model is not yet ruled out.
- **Where the resistance sits.** Not settled. Our approximate R_int and V2019's h_e > 60 kW/m²K
  give a conditional lower bound of ~16% on the kinetic share (A7); the actual share needs the
  unit-cell model (Stage 4). Liquid conduction between pillar and interface is expected to be the
  larger part (PROJECT_DERIVED expectation, not demonstrated), in which case geometry (h, l,
  meniscus shape) controls h_e and thin-film area near the pillar is the lever for novel designs.
- **Pitch optimum.** The MODEL_CALCULATED pitch sweep at d = 10, h = 25 µm peaks near
  l ≈ 40 µm for flat-meniscus dryout (V2019 calls l = 30 µm "near-optimal", citing zhu2016). That's the capillary-pressure vs permeability tradeoff in one table.
- **Feed length.** q''_dry ∝ 1/L² and, at fixed shape ratios, ∝ h². Distributed feeds (smaller
  L) and taller pillars both raise dryout, but taller pillars lower h_e (V2019 Fig. 5).
