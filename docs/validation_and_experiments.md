# Verification, Validation and Experiments

Companion to `docs/comsol_unit_cell_guide.md` (what to build) and
`docs/modeling_audit_2026-11.md` (why). This file says how we will know whether a model is
right, and what measurement would support or challenge each prediction.

## 1. Five kinds of evidence (kept separate)

| Kind | Question it answers | Current evidence | Next |
| --- | --- | --- | --- |
| Numerical verification | Does the code solve the stated equations accurately? | Stage 1 mesh studies (two builds); Python 2D EO solver verified (spectral convergence, conservation ~1e-13, `docs/eo_field_treatment_test.md`); Python device-model verification suites | Stage 2A similitude check Ũ·F = 1 (exact); Stage 2 energy balance; Stage 3 force identity |
| Reproducibility | Does an independent rebuild with the same inputs give the same output? | Stage 1 reference: Chase vs student within 0.05% at matching meshes (student result REPORTED_UNVERIFIED) | Stage 2 rebuild (`docs/stage2_undergrad_assignment.md`) |
| Analytical checks | Do tractable limits and orders of magnitude agree? | `analysis/case_hand_checks.py`: COMSOL K is 1.3% below the Sangani-Acrivos + Brinkman flat-top estimate; Stage 2 superheat inside conduction bounds | every stage has a hand bound in the guide |
| Literature benchmarking | Does our model reproduce a paper **under matched geometry, fluid, temperature, BCs and output definitions**? | none quantitative yet | V2019 Fig. 3 h_e(Pc) for the reference cell (after the closure gate); a matched-geometry K vs ravi2014 measurements (see §2) |
| Experimental validation | Does the model predict **our** measurements, with stated uncertainty, on conditions not used for calibration? | none | §3 progression |

A model that passes verification can still be wrong physics. A literature match can still be a
shared assumption. Only the last row tests the physics of our device.

## 2. What the papers can and cannot benchmark

| Source | Experiment or simulation | Use for us |
| --- | --- | --- |
| Vaartstra et al. 2019 (held) | Simulation. Its h_e "validation" uses `a` **fitted** to Adera sample A1 and checked on A2 (Sec. 3.4, Fig. 6a); its dryout check is a trend against Zhu et al. air-ambient experiments (Fig. 6b) | quantitative benchmark for the **reference cell only**, with matched closure, T_v = 100 °C, θ_rec = 15°, pure vapor; the fitted `a` is not a property |
| Adera et al. 2016 (held) | Experiment + model. Environmental chamber evacuated, then held near saturation at 24 °C / 3 kPa; RTD heater temperature; 62% of input power lost (sensible + convective + ohmic) and subtracted; several conditions boil before dryout (Sec. 4, Fig. 6) | experimental benchmark only if our model runs their geometry, substrate (650 µm), pressure and loss-corrected flux definition; shows how large heat-loss corrections are |
| Pujahari et al. 2022 (held) | Simulation only. EO model checked against an analytical 2D channel solution (Dutta & Beskok); thermal model checked against Zhu (dryout) and Adera (via V2019) data (Sec. 3.1). Uniform imposed field per unit cell; dilute KCl with ζ from Mala et al. | structural reproduction (done, `MODEL_STATUS.md`); its EO magnitudes are contextual, not a benchmark for a solved field or our chemistry |
| Fleming et al. 2020 (held) | Experiment: local meniscus curvature in pure vapor by thin-film interference; no-slip interface fits much better than shear-free (Sec. 3.1, Figs. 4-5) | evidence that the interface shear condition is uncertain; bracket, not a matched benchmark |
| Ravi et al. 2014 (held) | Experiment: K and capillary pressure determined **separately** for six Si micropillar wicks from dryout-threshold mass flow at several wicking lengths against gravity | candidate experimental benchmark for Stage 1/3 K if one of their geometries and fluids is run exactly (check their Table 3 and conditions before use) |
| Elkholy et al. 2023 (held) | Experiment/method: compares Lucas-Washburn and gravity-based rate-of-rise analyses | method reference: which analysis gives K/r_eff only and which separates them |
| Zhou & Adera 2026, Vaartstra 2020 (held, partial) | Simulation / review | context |

Rule: a paper number is a quantitative benchmark only when geometry, fluid, temperature,
vapor state, boundary conditions and output definition all match. Otherwise label it contextual
or order of magnitude. Never tune `a`, θ_rec, K or ζ to reproduce an expected number. If a
parameter must be fitted, fit it on one data set and report the comparison on a different one.

## 3. Experimental progression (small, ordered)

Each step isolates one physical question. Equipment is not prescribed. The bench must only
control the boundary conditions listed.

### Exp-1: unheated hydraulic transport (permeability)

- **Measure:** flow rate Q vs imposed pressure drop Δp through the wick with a **cover sealed on
  the pillar tops** (closed channel). Liquid temperature (μ changes ~3× between 25 and 100 °C).
- **Model to compare:** Stage 1Q item 3 (no-slip top) on metrology-corrected geometry. This
  comparison has unambiguous boundary conditions, so it tests geometry and fabrication (diameter,
  taper, scalloping, height), not meniscus physics.
- **Controls:** known length between pressure taps, entrance/exit losses measured with a blank or
  short sample, bubble-free degassed liquid.
- **Confounders:** gap or leak between cover and pillar tops (a thin gap carries a lot of flow);
  cover deflection; bubbles.

### Exp-2: capillary behavior, separately from hydraulic resistance

- **Directly measured:** advancing and receding contact angles on a flat oxide coupon processed
  with the wick (same oxidation, cleaning and storage); wicking-front position x(t) on the open
  wick (horizontal, reservoir level at the wick plane); the suction or reservoir offset at which a
  wetted wick can no longer hold liquid (receding limit).
- **Why rate-of-rise is not enough:** for horizontal or early-time wicking, x² ≈ 2 K P_front t /
  (μ ε), so the front speed gives only the **product** of permeability and the front capillary
  pressure (ravi2014 abstract; fleming2020 Sec. 1; elkholy2023 Sec. 4.1, Lucas-Washburn method). The advancing front's
  capillary pressure is set by the advancing angle, not by the receding Pc_max of V2019 Eq. 3 that
  limits dryout. Vertical rise separates K and P only when the height approaches the equilibrium
  height P/(ρg) ≈ 20 cm (proposed cell, 25 °C), far beyond a 20-25 mm coupon. Separation therefore
  needs K from Exp-1 (with the interface condition bracketed), or the dryout-threshold method of
  ravi2014, or a direct receding-limit (suction) test.
- **Confounders:** surface contamination drift (angles change within hours to days on oxide),
  evaporation during wicking, pinning at the feed edge, entrance meniscus.

### Exp-3: unheated electroosmosis

- **Directly measured:** supply voltage and current (I-V curve); flow rate (meniscus displacement
  in a calibrated capillary or mass change); stall pressure (the head difference at zero net
  flow); liquid conductivity σ (meter); liquid temperature; gas bubbles and pH at both electrodes.
- **Inferred, not measured:** mean field in the wick E0 = I F / (σ A) (F from Stage 2A/3E;
  avoids electrode overpotentials); effective ζ.
- **Calibration vs validation:** calibrate ζ_eff from **streaming potential or streaming current**
  under pressure-driven flow (zero net current, so no electrode reactions or overpotentials enter),
  then predict EO flow and stall pressure at other currents without refitting. Alternatively,
  calibrate on one current and predict the others; never report the calibration point as
  validation.
- **Scale (ILLUSTRATIVE: ζ = −50 mV, E0 = 1e4 V/m, 25 °C, proposed cell, `case_hand_checks.py`
  Sec. E):** U_EO ~ 3e-4 m/s, i.e. µL/min per mm of wick width; stall gradient ~4e3 Pa/m, i.e.
  tens of Pa (a few mm of water head) over a 10 mm wick; bulk drop ~100 V over 10 mm; current of
  order µA for a 1e-4 M electrolyte. **Reservoir head must be controlled to well below 1 mm**, or
  the hydrostatic signal swamps the EO signal.
- **Hydraulic boundary of the cathode compartment** (sealed / vented / gas headspace) must be set
  and recorded; it decides whether EO produces throughflow or a pressure rise.
- **Confounders:** electrolytic gas blocking pores, pH fronts changing ζ during the run, Joule
  heating, evaporation from the reservoirs, CO₂ uptake changing σ and pH of "DI" water, surface
  conduction at low ionic strength (Dukhin up to ~0.1 for pure water).

### Exp-4: heated passive, then heated EO, at matched conditions

- **Directly measured:** heater electrical power; heater or backside temperatures (RTD/
  thermocouples) at stated locations; reservoir mass loss rate; feed liquid temperature; vapor or
  ambient temperature and pressure; onset of dryout (visual plus temperature runaway).
- **Inferred:** heat reaching evaporation = power − losses (losses calibrated in a separate test;
  Adera et al. lost 62%); latent mass balance check: ṁ_measured ≈ (power − losses)/h_fg; T at the
  pillar-base plane by subtracting substrate conduction (1D only if spreading is negligible).
- **Vapor environment decides what can be compared.** Open laboratory air: dryout load and the
  capillary limit are comparable (V2019 Sec. 3.4 used air data only for dryout); h_e is **not**,
  because evaporation is diffusion-limited. Degassed pure-vapor chamber: h_e is comparable, and
  `a` may be calibrated on one condition and tested on another.
- **EO comparison:** identical sample, liquid, reservoir levels, heater power steps and vapor
  state; only the current changes. Record electrical input power separately (wick bulk, reservoir,
  electrodes, supply). A gross dryout gain is not a net benefit until electrical power and
  electrochemical degradation are counted.
- **Confounders:** heat losses, subcooled feed (sensible heating near the feed edge),
  noncondensable gas, boiling before dryout, nonuniform heater, wick aging between runs.

## 4. Model-to-experiment matrix

| Model output (stage) | Feasible measurement | Boundary condition that must be controlled | Major confounder | What disagreement would teach |
| --- | --- | --- | --- | --- |
| K, flat no-slip top (1Q) | Exp-1 closed-channel Q-Δp | sealed cover, tap spacing, liquid T | cover gap or leak | actual geometry differs from design, or Stokes model/normalization error |
| K, shear-free vs no-slip open interface (1Q, 3H) | Exp-2 wicking with K from Exp-1 | horizontal, reservoir level | advancing-angle uncertainty, evaporation | whether the interface is mobile (V2019) or immobilized (Fleming) on our surfaces |
| Pc_max, endpoint (3) | Exp-2 receding-limit test; receding angle on co-processed coupon | known hydrostatic offset | hysteresis, edge rounding, contamination | whether Eq. 3 (mean) or local depinning governs; actual θ_rec |
| F (2A, 3E) | AC conductance of the liquid-filled wick vs bulk σ | four-electrode or AC measurement, full wetting | electrode polarization; partial wetting | geometry/wetting differs; surface conduction |
| Ũ, EO flow and stall pressure (2A, 3E) | Exp-3 flow, current, stall head | reservoir levels, cathode compartment BC, temperature | gas evolution, ζ drift, head drift | thin-EDL/uniform-ζ model fails, or floor and pillar charge differ |
| ζ_eff (input, not output) | streaming potential/current | zero net current, known Δp | surface conduction, asymmetry potentials | chemistry differs from assumption; sets all EO magnitudes |
| T_base − T_v, h_e (2, 4) | Exp-4 in pure vapor only | vapor pressure, noncondensables, losses | heat loss, substrate spreading | closure or `a` wrong for our surface; conduction path wrong |
| q''_dry and location (6) | Exp-4 (air or vapor) | feed head, feed edge, heater footprint, uniformity | boiling before dryout; loss calibration | K(Pc), endpoint or supply model wrong; another limit dominates |
| EO dryout gain, flooding (6) | Exp-4 EO vs passive at matched conditions | everything in Exp-3 and Exp-4 | electrochemistry over time | gain not realized (field, ζ or gas), or flooding at low load |
| Mass and energy closure (6) | reservoir mass loss vs power | sealed reservoirs, evaporation from reservoirs | reservoir evaporation, leaks | losses or latent fraction misjudged |

## 5. Real conditions to record for every test

Liquid temperature (μ, σ, ζ, γ), surface state (oxide process, cleaning, time since cleaning,
measured angles), fluid chemistry (electrolyte, concentration, conductivity, pH, degassing),
reservoir heads relative to the wick (mm), cathode compartment state, vapor environment (air or
vapor, pressure, humidity), and the metrology of the tested sample (d, pitch, height, sidewall).
A model comparison uses these measured values, not the design values in
`parameters/unit_cell/cases.yaml`.
