# Whole-Wick Model: Assumptions Register and Unresolved Decisions

Scope: the 2D reduced-order wick model in `src/nanofab/whole_wick/`. The unit-cell COMSOL model
(`docs/comsol_unit_cell_guide.md`) is treated as a **property generator** only. Its outputs
K(Pc), h_e(Pc), Pc_max and k_eff feed lookup tables, and the whole-wick model never re-solves
unit-cell physics. Unit-cell assumptions A1-A19 are in `docs/assumptions.md` and still apply to
every table extracted from COMSOL.

Each entry gives the default the code uses now, so runs are reproducible while the decision is
open. "Switch" names the config key that changes it.

## Pressure reference and sign convention (decision D7, adopted)

```text
p_v      vapor pressure, uniform, = p_sat(T_v). This is the pressure REFERENCE (gauge zero).
p_l      liquid pressure in the wick, reported as gauge: p_l - p_v  (<= 0 in a working wick).
Pc       capillary pressure = p_v - p_l  (>= 0 when the meniscus is concave, liquid in suction).
Pc_max   largest Pc the local pillar array can sustain before the meniscus depins at the
         receding angle (vaartstra2019 Eq. 3). A cell is pressure-dry when Pc > Pc_max.
Feed     p_l = p_v - Pc_feed at feed edges; Pc_feed = 0 for a flat reservoir meniscus.
```

Liquid flows toward lower p_l, that is toward higher Pc. Pc therefore increases with distance
from the feed and peaks at the point farthest from it.

## Assumptions register

| # | Assumption | Default | Consequence / when it fails | How to test |
| --- | --- | --- | --- | --- |
| W1 | 2D planform model, averaged through the wick thickness (D6) | 2D | Misses through-thickness K variation and 3D hotspot spreading. OK when lateral scales >> t_s + h (here 125 µm). | Compare with explicit small-array COMSOL (validation item 6). |
| W2 | Homogenization: cell size >= several pitches | warn if dx < 3 l | Below about 3 pitches a "cell" is not a representative volume, so grid refinement past that converges the *equations*, not the *physics*. | Grid study in `validation/`. |
| W3 | Darcy flow with tabulated K(Pc) | flat-meniscus analytical K | Flat K is expected to exceed the curved-meniscus K (V2019 Sec. 2.2 says curvature lowers K), but it is not established as a bound on dryout: the legacy estimate is +22.7% vs vaartstra2019 (`docs/assumptions.md` A2). | Replace with COMSOL Stage 3/5 K(Pc). |
| W4 | Uniform vapor pressure p_v = p_sat(T_v) (D2) | uniform | Fails in air ambient (diffusion-limited) or a vapor space with significant pressure drop. | Vapor-space model, if needed later. |
| W5 | Evaporation enters as a sink per unit footprint area (D4) | areal flux = depth-averaged volumetric sink | Identical in 2D. Differs only in a 3D model. | n/a in 2D |
| W6 | T is the pillar-base plane temperature (T_e of vaartstra2019) | | h_e already contains the through-wick conduction resistance, so k_eff is used **only** for in-plane conduction. Adding a through-plane k_eff resistance in series would double count. | Check that the COMSOL h_e definition uses the same T_e surface (A19). |
| W7 | In-plane conductance G = k_sub t_sub + k_eff h | Si 100 µm dominates (> 95% of G) | The k_eff model barely matters in-plane. | Sweep k_eff 0-2x. |
| W8 | Liquid sensible heat and heat advection neglected | | Ja ≈ 0.02 (A3). | Energy audit in the Stage D report. |
| W9 | Fluid properties uniform at T_v; per-cell override allowed | CoolProp at T_v | Superheats are a few K (A15). | |
| W10 | Steady state | | No transient dryout or rewetting. | Phase-later transient. |
| W11 | Dryout = onset: the first wet cell with Pc >= Pc_max (D3) | `dry_model: onset` | Conservative. Real wicks can run with partial dry patches. | `dry_model: dry_patch` gives a post-onset estimate (heuristic). |
| W12 | Pc_max is the static receding-angle limit of the local unit cell | vaartstra2019 Eq. 3 | No contact-angle hysteresis, no meniscus-touches-floor limit for short pillars. | COMSOL Stage 3 should report the depinning curvature directly. |
| W13 | Feed is a fixed-pressure edge, Pc_feed = 0 | left and right edges | Makes the liquid-supply criterion algebraically equal to the pressure criterion (see D3). | Optional finite feed capacity `max_feed_mass_flow`. |
| W14 | Thermal edges adiabatic | vaartstra2019 Sec. 2.3 symmetry | Heat leaves only by evaporation (or by an edge sink if one is configured). | Energy balance report. |
| W15 | Electrical assistance is electrowetting only, driven by a prescribed voltage map V(x, y) with no field or circuit solve (D1) | | V acts only through θ_V and Pc_max(V). See the electrowetting section below. | |
| W16 | No electrical power and no Joule heat enter the model (D14) | electrical power UNRESOLVED | An ideal dielectric draws no DC current. Leakage, AC displacement current and driver losses are unknown. | Electrode and circuit model. |
| W17 | In-plane gravity component | 0 (horizontal wick) | Bo ~ 1e-4 (A4). | Set `gravity_inplane`. |
| W18 | h_e(Pc) for the baseline is a **placeholder** (60 kW/m²K, rising 30% from Pc = 0 to Pc_max) | ASSUMED | The 30% variation and center-high trend from vaartstra2019 become **inputs, not predictions**. Do not cite them as validation. | Replace with COMSOL Stage 4/5 h_e(Pc). |
| W19 | Property transitions happen at cell faces with harmonic face averaging (D5) | harmonic | Unit-cell tables assume periodicity, which fails within about one pitch of a geometry change. | Explicit-array COMSOL across one transition. |
| W20 | Dry cells evaporate nothing (`h_dry = 0`) in `dry_patch` mode | 0 | Pessimistic for post-onset spreading. | |
| W21 | Graded wicks are built from discrete cell types, each with its own Pc tables (D9) | | Avoids a 4-D (d, l, h, Pc) table the COMSOL sweep cannot fill. | |

## Electrowetting model (D1, adopted)

```text
eta        = epsilon_r epsilon_0 V^2 / (2 gamma t_dielectric)          electrowetting number
cos θ_V    = min( cos θ_0 + eta , cos θ_sat , 1 ),  never below cos θ_0  (Lippmann-Young + saturation)
Pc_max(V)  = 2 gamma cos θ_V / r_eff
r_eff      = 2 (l^2 - pi d^2/4) / (pi d)    so that Pc_max(0) equals vaartstra2019 Eq. 3
dryout     first wet cell where required (p_v - p_l) >= Pc_max(V)
```

Voltage changes **only** the receding contact angle and therefore Pc_max. It does not change
h_e, K, k_eff or the heat input. Because cos θ_V ≤ 1, the largest possible gain in Pc_max (and,
with constant h_e, in dryout load) is **1/cos θ_0 − 1**: 3.5% at θ_0 = 15°, 15% at 30°, 41% at
45°. The 1 µm SiO₂ coating reaches that limit at about 11 V.

| # | Assumption | Default | Consequence / when it fails | How to test |
| --- | --- | --- | --- | --- |
| EW1 | Lippmann-Young applies to the **receding** angle at the pinning edge | θ_0 = receding_contact_angle (15°) | Lippmann describes the static apparent angle. Hysteresis under voltage is not modeled. | Measured receding angle vs V on the coated pillar surface. |
| EW2 | Saturation is a hard cap at θ_sat | θ_sat = 0° (ASSUMED, i.e. only the cos ≤ 1 cap) | Real saturation is typically 30-80° on flat surfaces (mugele_baret2005, cited from memory). If θ_sat ≥ θ_0, electrowetting has **no** effect. | Sweep θ_sat (`run_sweeps.py`); measure. |
| EW3 | The meniscus stays pinned at the pillar top; voltage moves only the angle limit | | K(Pc) and h_e(Pc) are V-independent at a given Pc. Higher Pc_max lets the wick reach higher Pc, where K and h_e are taken from the same tables, extended up to Pc_max at the maximum V. | COMSOL Stage 3 at θ_V. |
| EW4 | Dielectric is ideal: uniform ε_r and t_dielectric, no leakage, no charge trapping, no breakdown | ε_r 3.9 (ASSUMED), t 1 µm (USER_PROVIDED) | The run reports V/t (V/µm). It does **not** check breakdown strength or charging. | Dielectric characterization. |
| EW5 | The dielectric layer adds no thermal resistance and no diameter | | 1 µm SiO₂ (k ≈ 1.4 W/m/K) adds about 0.7 µK·m²/W to every heat path, and it thickens the pillars, which lowers K and raises Pc_max. Neither is included. | Add to the h_e and geometry inputs once the coating is fixed. |
| EW6 | Every pillar in an actuated cell sees the full V | | Requires an electrode layout that addresses each cell. The electrode geometry is undefined. | Electrode design. |
| EW7 | r_eff comes from geometry (`r_eff: from_geometry`) | | Can be overridden per cell type with a value. | |
| EW8 | Selective actuation uses the passive-uniform pressure ratio at that case's dryout load (D11) | threshold 0.7 | The mask is fixed before actuation, so it does not adapt as the dry spot moves. | Threshold sweep. |

**Mechanisms deliberately NOT modeled.** These are separate physical mechanisms. Each would need
its own governing equations, material data and experimental evidence before it enters this model:

- **Electrohydrodynamic (EHD) pumping.** A body force on the liquid. It would need a new term in
  Darcy's law plus charge-relaxation physics, not a change in K or Pc_max.
- **Charge transport.** Ionic conduction, leakage through the dielectric, double-layer charging.
  This needs a circuit and electrochemistry model, and it sets the electrical power.
- **Direct field enhancement of evaporation.** Any change of h_e or interfacial kinetics with
  field. The code no longer lets h_e depend on V, and CSV tables with a voltage or field column
  are rejected. With adiabatic edges an h_e-only boost changes temperature, not dryout load.

## Unresolved decisions

Each decision lists the default the code uses now and what would resolve it.

**D1 - Electrical-assistance mechanism. (Adopted: electrowetting.)** See the electrowetting
section above. Electroosmosis is modeled separately in the 1D passive-vs-assisted model
(`src/nanofab/electroosmosis/`, `docs/electroosmosis_model.md`); it is NOT part of this 2D
whole-wick model. EHD pumping, charge transport and direct evaporation enhancement are excluded
until each has its own physics and experimental justification. Still open inside this decision:
the coating material and thickness, θ_sat on the real pillar surface, and whether the receding
angle follows Lippmann at all (EW1, EW2).

**D2 - Vapor-pressure treatment.** Default: uniform p_v = p_sat(T_v), pure vapor. Needs a
decision when the target environment is known (vacuum chamber, vapor chamber, or air).

**D3 - Dryout definition.** Default: onset, meaning the first wet cell where Pc >= Pc_max.
The requested liquid-supply criterion needs care. With fixed-pressure feeds and Darcy flow,
a cell's maximum net inflow (its pressure lowered to p_v - Pc_max, neighbors fixed) is
`m_evap·A + (Pc_max - Pc)·ΣT_faces`. So "demand exceeds local supply" is *algebraically the same*
as "Pc exceeds Pc_max". The code reports both margins (the test suite checks they agree) and adds
one criterion that is independent: total evaporation exceeding a finite feed capacity
`max_feed_mass_flow`. Other candidates still to decide: dry-area fraction > X%, or temperature
runaway.

**D4 - Boundary flux or volumetric source.** Default: areal sink (W5). Equivalent in 2D. It
matters only if we move to 3D or let K vary through the wick thickness.

**D5 - Averaging across material transitions.** Default: harmonic mean of the cell mobility
ρKh/μ and of the conductance G at faces. This keeps flux continuous for flow normal to the
transition. Arithmetic averaging is available (`face_average: arithmetic`) for sensitivity runs.

**D6 - 2D or 3D.** Default: 2D (W1). Move to 3D only if the explicit-array comparison shows
errors above about 10% near transitions or hotspots.

**D7 - Pressure reference and sign.** Adopted as above.

**D8 - Provenance labels.** The five requested labels are mapped onto the repo's enforced classes
(`docs/SOURCE_AUTHORITY.md`). Two classes were added: `COMSOL_EXTRACTED` and `USER_PROVIDED`.

**D9 - Graded-geometry parameterization.** Default: discrete cell types (fabricated graded wicks
are zoned anyway). If continuous grading is needed, restrict it to one geometry axis (pitch at
fixed d and h) so the COMSOL sweep stays tractable.

**D10 - Base wick for the electrowetting cases.** Default: the uniform wick, so the
electrowetting effect is isolated. A graded wick with selective actuation is the obvious fifth case.

**D11 - Selective-actuation vulnerability criterion.** Default: cells with Pc/Pc_max >= 0.7 in the
passive uniform case at that case's dryout load. Alternatives: a fixed operating load, the top-N%
of cells, or a temperature criterion.

**D12 - Figure of merit.** Default reports: dryout-load gain, Δθ, ΔPc_max, applied voltage,
actuated area fraction, and whether actuation is global or selective. No per-watt or efficiency
metric is reported until D14 is resolved. Actuated area is **not** a proxy for power.

**D13 - Heater footprint vs wick footprint.** Default: the heat map covers the whole wick, and
"dryout heat flux" is the total heat input at dryout divided by the wick footprint. The peak local
heat flux is also reported.

**D14 - Electrical power. (Open.)** Reported as UNRESOLVED. It needs an electrode layout,
dielectric leakage and capacitance, the drive waveform (DC or AC frequency) and driver losses.
Until then the model adds no Joule heat and makes no power-saving claim, including for selective
vs global actuation.
