# Whole-wick reduced-order model

A 2D grid of homogenized cells representing the full capillary-fed micropillar wick, with graded
geometry, mixed materials, non-uniform heat flux and electrowetting actuation by a prescribed
voltage map. The unit-cell COMSOL model is **not** the wick model. It generates property tables
(`unit_cell/README.md`).

**Status: runs and conserves mass and energy to machine precision, but is NOT validated.**
h_e(Pc) is a placeholder, K(Pc) is the flat-meniscus analytical value, and the electrowetting
parameters ε_r and θ_sat are ASSUMED (`validation/VALIDATION_REPORT.md`).

Code: `src/nanofab/whole_wick/`. Assumptions and open decisions:
`docs/whole_wick_assumptions.md`. Read it before trusting any number.

## Equations (2D planform, averaged through the wick height h)

| | Equation | Notes |
| --- | --- | --- |
| Darcy | u_l = −(K/μ_l)(∇p_l − ρ_l g) | superficial velocity; K = K(Pc) per cell, independent of V |
| Liquid mass | ∇·(ρ_l h u_l) = −m''_evap | = ∇·(ρ_l u_l) = −ṁ_evap with ṁ_evap = m''_evap/h |
| Evaporation | m''_evap = q''_evap / h_lv | |
| Evaporative flux | q''_evap = h_e(Pc, geometry, material)·(T − T_sat) | h_e per vaartstra2019 Eq. 4; **no voltage dependence** |
| Heat | ∇·(G∇T) + q''_in − q''_evap = 0 | G = k_sub t_sub + k_eff h [W/K]; no electrical heat (D14) |
| Electrowetting | cos θ_V = cos θ_0 + ε_r ε_0 V² / (2 γ t_dielectric) | capped at min(1, cos θ_sat) |
| Capillary limit | Pc_max(V) = 2 γ cos θ_V / r_eff | r_eff from geometry reproduces vaartstra2019 Eq. 3 at V = 0 |

T is the pillar-base plane temperature (T_e of vaartstra2019). h_e already includes conduction
through the pillar/liquid layer, so k_eff enters **only** the in-plane conductance (W6).

**Pressure convention.** p_v is uniform and is the reference. Pc = p_v − p_l ≥ 0 is the pressure
difference the wick must sustain. The feed edges are held at Pc = Pc_feed (0). Liquid flows toward
higher Pc.

**Dryout.** A cell dries when the required pressure difference p_v − p_l reaches Pc_max(V) for
that cell. The liquid-supply margin (Pc_max − Pc)·ΣT_faces is the same criterion expressed as
flow headroom. The independent supply limit is an optional finite feed capacity. The dryout load is
the first heat-load scale at which any cell reaches either limit (Brent root-finding).

**Bound on electrowetting.** cos θ_V ≤ 1, so Pc_max can rise by at most 1/cos θ_0 − 1. That is
3.5% for the 15° receding angle of the baseline wick, reached at about 11 V with 1 µm SiO₂.

**Not modeled:** EHD pumping, charge transport, and direct field enhancement of evaporation.
These are separate mechanisms that need their own physics and experimental justification
(`docs/whole_wick_assumptions.md`, electrowetting section).

**Numerics.** Cell-centered finite volume with harmonic face averaging (D5), a sparse direct
solve, and Picard iteration on Pc for the K(Pc) and h_e(Pc) coupling.

## Inputs (`parameters/whole_wick/`)

| File | Contents |
| --- | --- |
| `wick.yaml` | domain, grid, substrate, T_v, γ (from the fluid), feeds, thermal edges, gravity, solver |
| `cell_types.yaml` | materials; dielectric coatings (ε_r, t_dielectric, θ_sat); cell types (d, l, h, material, θ_0, r_eff, coating); property recipes (analytical, COMSOL table, placeholder) |
| `cases.yaml` | heat maps (uniform, Gaussian hotspots, patches, CSV file); geometry/material maps (uniform, feed-distance bands, rectangles, file); voltage maps (zero, uniform, rectangles, file, vulnerability mask); the four cases |

Every number carries `provenance`. Bare numbers are rejected. Reports map the provenance to
*literature value / extracted from COMSOL / analytical value / estimate / user-provided value*.
Lookup tables are never extrapolated: an out-of-range query raises an error or is clipped, warned
about, and logged. Applying V > 0 to an uncoated cell type is an error.

## Cases

| Case | Geometry | Voltage |
| --- | --- | --- |
| `passive_uniform` | d10/l30/h25 everywhere | V = 0 |
| `global_ew` | uniform | 10 V everywhere |
| `selective_ew` | uniform | 10 V where Pc/Pc_max ≥ 0.7 in `passive_uniform` at its dryout |
| `passive_graded` | l = 40 / 30 / 22 µm in bands by distance from the feed | V = 0 |

## Outputs

```powershell
python whole_wick/run_cases.py        # results/cases/: Stage C-F reports, maps, per-cell CSVs
python whole_wick/run_sweeps.py       # results/sweeps/: load, V, θ_0, θ_sat, t_dielectric, threshold, K(Pc), h_e, pitch
python validation/run_validation.py   # validation/VALIDATION_REPORT.md
pytest tests/test_whole_wick.py
```

Per case: inlet/outlet/evaporation flow and mass residual; heat input, evaporative removal, edge
loss and energy residual; T map and T_max; first dry location, dryout heat load (mean and peak
local flux), and the required p_v − p_l vs Pc_max(V) at that cell. The comparison table reports
the actuation (none, global or selective), applied voltage, actuated area, θ_0 → θ_V, Δθ,
Pc_max(0) → Pc_max(V), dryout heat flux and location, and its change relative to
`passive_uniform`. **Electrical power is reported as UNRESOLVED.** No power or efficiency claim is
made until an electrode and circuit model exists.

## Validation procedure

1. Verification: `pytest` (uniform Darcy, permeability scaling, conservation, monotonicity,
   dryout vs exact 1D, lookup bounds, Lippmann-Young and its caps, V = 0 equals passive,
   dryout ∝ cos θ_V, selective masking).
2. Grid refinement, solver sensitivity and electrowetting verification:
   `validation/run_validation.py` sections 1-3b.
3. Replace the placeholders with COMSOL Stage 3-5 tables (`unit_cell/`), then rerun.
4. Explicit small-array COMSOL vs homogenized model across a graded transition (report section 6).
5. Device data (adera2016, zhu2016 via vaartstra2019 Fig. 6). Only this step validates the model.
   Electrowetting additionally needs measured receding angle vs V on the coated pillars.
