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

1. **Unit cell (Level 1):** one pillar cell in COMSOL -> K, h_evap, interface shape, local flux.
2. **Reduced order (Level 2):** lookup tables K(geometry), h_evap(Pc, geometry); wick treated as
   an effective medium.
3. **Device / cold plate (Level 3):** effective properties -> chip temperature, dryout location,
   pressure requirements.

The focus is the evaporator/wick. The full loop (vapor, condenser, return, pump) comes later.

## Fidelity ladder for the interface

| Model | Interface treatment | Where it is appropriate |
| --- | --- | --- |
| A | prescribed evaporative heat flux | device-scale energy-balance checks |
| B | prescribed effective h_evap | device scale only (h_evap comes *from* the unit cell) |
| C | kinetic flux (linearized BTE as in V2019, or Schrage) on a fixed meniscus | unit cell (Stage 4); the model used to extract h_evap |
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

1. Analytical sanity checks (`analysis/baseline_estimates.py`, `tests/`).
2. Literature replication (Vaartstra et al.; targets in `docs/literature_notes.md`).
3. Experiment (thermal test vehicle -> wick -> vapor -> condenser; heater power and temperature,
   wick temperature, pressure, flow, dryout onset).

Mesh convergence is required for every reported number (change below 1-5% between the two finest
meshes).

## Source tracking

Every input carries a provenance class and a registered source id (enforced by
`src/nanofab/params.py`; see `docs/SOURCE_AUTHORITY.md`). This replaces the original brief's
`MEASURED`/`LITERATURE`/`CALCULATED`/`ASSUMED`/`DESIGN`/`TARGET` tags. Uncertain inputs are swept, not trusted: contact angle,
accommodation coefficient, interface geometry, contact resistance, inlet pressure.

## Phases

1. Learn and reproduce the reference unit cell. **Current.**
2. Validate: mesh convergence, analytical checks, literature comparison.
3. Parameter sweeps (d, l, h, q'', Pc, θ).
4. Validated baseline micropillar wick (control case).
5. Novel geometries, one change at a time.
6. Reduced-order relations K(geometry), h_evap(Pc, geometry).
7. Device-scale analysis at realistic chip sizes and heat loads.

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
