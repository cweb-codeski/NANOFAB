# NANOFAB Modeling Context

Condensed project brief. It's the reference for anyone (human or AI assistant) working in this
repo.

## Goal

Two-phase evaporative cooling with micro/nanostructured wicks for high-heat-flux electronics
(GPUs/TPUs). Quantify:

> How much thermal performance improvement does a proposed wick architecture provide, what
> physical mechanism produces it, and what tradeoffs does it introduce?

The chain to connect: wick geometry -> fluid/thermal physics -> chip temperature and dryout ->
experimental validation -> system-level usefulness.

## Causal model

```text
wick geometry -> pore / interface geometry -> capillary pressure + permeability
-> liquid transport -> liquid availability at evaporation sites
-> thin-film / interfacial evaporation -> heat removal -> surface temperature -> dryout limit

material + geometry -> solid conduction -> temperature distribution -> local evaporation
```

These two branches are coupled. That coupling is the core modeling problem.

## Modeling hierarchy

COMSOL is the primary environment for the physical simulations. Python keeps hand-calculation
equivalents, limiting-case checks, parameter management, export processing, plotting and
targeted cross-checks. The staged plan (Stages 1, 1Q, 2, 2A, 3-7) with prerequisites, equations,
boundary conditions, checks and handoffs is **`docs/comsol_unit_cell_guide.md`**. That guide is
the only workflow plan; this section only summarizes it.

1. **Unit cell (Stages 1-5):** one pillar cell in COMSOL -> K, F, EO flux, interface shape,
   h_e, and their mesh evidence, per geometry case.
2. **Device (Stage 6):** COMSOL effective-medium model (boundary PDE for the liquid supply on the
   substrate, Heat Transfer in Solids, optional current conservation for EO) using the Stage 5
   tables. Not explicitly meshed pillars. The Python reduced-order models (`src/nanofab/device/`,
   `src/nanofab/whole_wick/`, `src/nanofab/electroosmosis/`) are cross-checks of it, not
   predictions.
3. **Cold plate (Stage 7):** only when spreading, contacts, reservoirs or vapor transport change
   Stage 6 predictions.

The focus is the evaporator/wick. The full loop (vapor, condenser, return, pump) comes later.

## Geometry cases

Two cases are kept separate (`parameters/unit_cell/cases.yaml`):

- **Reference** d/l/hp = 10/30/25 µm (Vaartstra et al. 2019): learning, analytical checks and
  matched literature comparison.
- **Proposed** d/l/hp = 20/50/40 µm: user-specified design inputs, not measured fabrication
  dimensions. Array 205 × 155 (outer pillar-edge span 10.22 × 7.72 mm); coupon placement, CAD
  scale, feed length, substrate, heater, reservoirs and vapor environment are UNRESOLVED.

No COMSOL result (K, F, h_e, menisci) or fitted parameter is transferred between them.

## Fidelity ladder for the interface

| Model | Interface treatment | Where it is appropriate |
| --- | --- | --- |
| A | prescribed evaporative heat flux | device-scale energy-balance checks |
| B | prescribed effective h_evap | device scale only (h_evap comes *from* the unit cell) |
| C | kinetic flux on a fixed meniscus under a **named** closure (linearized BTE as in V2019 once the Meland-Ytrehus treatment is verified, or Schrage; never mixed) | unit cell (Stage 4); the model used to extract h_evap; pure saturated vapor only |
| D | thin-film / disjoining pressure / Knudsen layer | only if C fails to match the reference |

Using a prescribed h_evap inside the unit cell would be circular, because extracting h_evap is the
unit cell's purpose. Two-phase interface tracking (Level Set or Phase Field) is not used; the
meniscus is a fixed shape computed separately, as in the reference work.

## Central design tradeoff

Small pores give high capillary pressure but low permeability. Large pores give the reverse. For a
1D wick, q''_dry ∝ K·Pc / L². So:

- geometry sets the product K·Pc (with an optimum in pitch for fixed diameter);
- architecture sets L (distributed feeds reduce it quadratically);
- thin-film area, contact-line length and solid conduction set h_evap.

Each proposed geometry must come with a hypothesis tied to one of these. For example: *radial
feed channels reduce L, delaying dryout at the evaporator center.*

## Metrics

| Metric | Symbol | Notes |
| --- | --- | --- |
| dryout heat flux | q''_dry | primary comparison metric |
| max temperature / superheat | T_max, ΔT | above T_v |
| thermal resistance | R''_th | reported as substrate + wick + interface |
| permeability | K | Stage 1 |
| capillary pressure | P_c | at depinning |
| evaporative HTC | h_evap | footprint-based |
| wick pressure drop | ΔP | vs P_c |

## Fair comparison rules

Same heat input, footprint, feed length, working fluid, saturation condition, base material,
substrate thickness, comparable wick height and fabrication constraints. Change one variable at a
time.

## Validation hierarchy

Five kinds of evidence, kept separate (`docs/validation_and_experiments.md`): numerical
verification, reproducibility (independent rebuilds), analytical checks
(`analysis/case_hand_checks.py`, `tests/`), literature benchmarking (only under matched
geometry, fluid, temperature, BCs and output definitions), and experimental validation (our
measurements, with calibration kept separate from validation conditions).

Mesh convergence is required for every reported number. A finest-pair change ≤ 1e-2 is a
labeled PRELIMINARY learning check; qualified tables and reported numbers need ≤ 1e-3 with at
least three systematic levels (`docs/cell_table_contract.md`). Track the physical outputs (K, EO
flux, superheat, h_e, Q_evap, meniscus clearance, dryout load), not residuals.

## Source tracking

Every input carries a provenance class and a registered source id (enforced by
`src/nanofab/params.py`; see `docs/SOURCE_AUTHORITY.md`). This replaces the original brief's
`MEASURED`/`LITERATURE`/`CALCULATED`/`ASSUMED`/`DESIGN`/`TARGET` tags. Uncertain inputs are swept, not trusted: contact angle,
accommodation coefficient, interface geometry, contact resistance, inlet pressure.

## Phases

The phases map onto the stages of `docs/comsol_unit_cell_guide.md`; the guide's §0 table is the
current status.

1. Learn and reproduce the reference unit cell (Stages 1, 2, 2A). **Current.**
2. Qualify: Stage 1Q, mesh convergence, conservation and analytical checks; reproduce V2019 only
   under matched conditions and a verified closure.
3. Menisci, evaporation and sweeps per case (Stages 3-5).
4. The proposed geometry, built as its own case, never by transferring reference results.
5. Device-scale COMSOL model (Stage 6), cross-checked against the Python reduced-order models.
6. Experiments in the order of `docs/validation_and_experiments.md` §3.

## Coding conventions

- SI units internally. Units documented in docstrings and parameter files.
- No hard-coded physical constants in analysis code: parameters come from YAML, fluid properties
  from CoolProp.
- Readable equations with cited sources. Every analytical model has a test of a known limit.
- Modular scripts over notebooks. Python supports COMSOL; it doesn't replace it.

## How an assistant should help

For any calculation: name the physical quantity, state the equation, define variables and units,
list assumptions and required inputs, implement cleanly with sanity checks, and explain how the
output connects to the cooling design. Prefer the simplest model that captures the physics needed
for the current question.

The goal is to answer *why does this wick perform differently?*, not just *which one is cooler?*
