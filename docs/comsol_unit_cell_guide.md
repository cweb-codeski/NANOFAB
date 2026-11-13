# COMSOL Modeling Workflow: Unit Cell to Device (canonical guide)

Revision 2026-11-13. This guide is the **single plan** for NANOFAB physical simulations. COMSOL
6.3 runs the physics; Python keeps the hand calculations, limiting-case checks, parameters,
export processing and plots. The audit behind this revision, with evidence and consequences, is
`docs/modeling_audit_2026-11.md`. Validation and experiments: `docs/validation_and_experiments.md`.
Undergraduate Stage 2 task: `docs/stage2_undergrad_assignment.md`.

**Nothing in this guide was executed in COMSOL by its author.** COMSOL was not installed in the
environment that wrote it, and the licensed modules are unknown. Menu paths follow COMSOL 6.3
naming as far as known. Where a node may differ by version or license, the text says so and gives
a fallback that uses only the base package. Check every step against the Help page of the node
(F1) and record any correction here.

---

## 0. Evidence vocabulary and current status

Every result is one of:

| Label | Meaning | Can it feed a later stage? |
| --- | --- | --- |
| `EXECUTED_EVIDENCE` | raw export, `.mph` (or `.java` model file / report), selections and run log are in the repo or archived with checksums | yes, after qualification |
| `STORED_DERIVED` | numbers from one of our runs, but setup evidence incomplete | learning and cross-checks only |
| `REPORTED_UNVERIFIED` | a team member reports having run it; numbers transcribed in `data/reported/` | comparison only |
| `PROPOSED` | instructions in this guide, not yet run | no |
| `ANALYTICAL_ESTIMATE` | Python hand calculation (`analysis/case_hand_checks.py`) | sanity bounds only |
| `SYNTHETIC` | made-up parser or fixture data | never |

A missing export does not prove a run was never done. It only means the repository cannot use it
yet.

| Stage | Reference cell (10/30/25 µm) | Proposed cell (20/50/40 µm) |
| --- | --- | --- |
| 1 flat-top flow | Chase `STORED_DERIVED` (K = 2.9698e-11 m² at hmax 0.5 µm); student `REPORTED_UNVERIFIED` (K = 2.9711e-11 m² at 0.25 µm). Builds agree within 0.05% at matching meshes (`analysis/stage1_crosscheck.py`) | `PROPOSED` |
| 1Q qualification (gradient/inertia, interface bracket, reduced-cell check) | `PROPOSED` | `PROPOSED` (high priority: Re_K ≈ 2-8 at device gradients) |
| 2 flat-top conduction | Chase `REPORTED_UNVERIFIED` (two unlabeled T_base values, see §6) | `PROPOSED` |
| 2A isothermal EO, flat top | 2D Python limit done (`docs/eo_field_treatment_test.md`); COMSOL `PROPOSED` | 2D Python limit done; COMSOL `PROPOSED` |
| 3 menisci + transport | `PROPOSED` | `PROPOSED` |
| 4 evaporation on menisci | `PROPOSED`, blocked for V2019 replication by the accommodation-treatment source gap | `PROPOSED` |
| 5 property tables | none qualified | none |
| 6 device (COMSOL effective medium) | `PROPOSED`; Python reduced-order models exist as cross-checks | `PROPOSED`, blocked by unresolved feed, heater, reservoir, substrate |
| 7 cold-plate integration | not triggered | not triggered |

---

## 1. Cases: reference and proposed design are separate

Definitions: `parameters/unit_cell/cases.yaml`, loaded by `nanofab.geometry_cases`.

| Quantity | Reference `reference_d10_l30_h25` | Proposed `proposed_d20_l50_h40` |
| --- | --- | --- |
| Role | learning, analytical checks, matched literature comparison (V2019) | user-specified design input, **not** measured fabrication dimensions |
| d / l / hp | 10 / 30 / 25 µm (PRIMARY_SOURCE, V2019 Sec. 3.2) | 20 / 50 / 40 µm (USER_PROVIDED) |
| Solid fraction φ, porosity | 0.0873, 0.913 | 0.1257, 0.874 |
| θ_rec | 15° (V2019 Sec. 3.1) | 15° ASSUMED, transferred without evidence; sweep 5-40° |
| Pc_max, V2019 Eq. 3 (100 °C) | 2177 Pa | 1636 Pa |
| Feed length | 5 mm (V2019) | **UNRESOLVED** |
| Substrate, heater, reservoirs, vapor environment | V2019 idealization | **UNRESOLVED** |
| Array | n/a | 205 columns × 155 rows, outer pillar-edge span 10.22 × 7.72 mm |
| Coupon | n/a | corners (−12.5, −10) and (12.5, 10) mm, as reported |

Rules:

1. **No transfer.** K, F, h_e, menisci, fitted accommodation and contact angle belong to one case.
   The proposed cell gets its own Stage 1-5 results. Scripts take `--case` so an export without
   d/l/hp columns is never normalized with the wrong footprint
   (`analysis/permeability_from_comsol.py`, `analysis/heat_transfer_from_comsol.py`).
2. **Coupon ≠ array ≠ feed distance.** Coupon size, active-array size, heater footprint,
   substrate thickness, reservoir geometry and feed distance are separate keys. None is inferred
   from another. Do not substitute a 20 × 15 mm active wick.
3. **Array placement is unresolved.** With the start reference (−10, −7.5) mm read as either the
   first pillar center or the array's edge corner, the array spans about x −10.0..0.2 mm and
   y −7.5..0.2 mm. It is offset ~4.9 mm and ~3.6 mm from the coupon center. That could be the
   intended layout, a different meaning of the start reference, or a CAD scale or unit problem
   (`Wick_V2.cat` references `Coupon_V2.cat`; unit metadata unchecked). Resolve it in the native
   CAD application before any device-scale model of the proposed wick.
4. **One `.mph` per stage and case.** Name files `ref_stage1_flow.mph`,
   `prop_stage1_flow.mph`, and so on, and register them in `comsol/unit_cell/README.md`.

Hand-calculation equivalents for every stage and both cases:

```powershell
python analysis/case_hand_checks.py
```

---

## 2. Stage roadmap

Every stage answers one physical question and passes its verification checks before its outputs
feed another stage.

```text
Stage 1  flat-top passive flow ─┬─> 1Q qualification (gradient/inertia, interface bracket, cell size)
                                ├─> Stage 2A isothermal EO, flat top (needs Stage 1 only)
Stage 2  flat-top conduction ───┘     │
         (needs Stage 1 geometry)     │
Stage 3  menisci (Young-Laplace) ─> 3H passive K(Pc) ─> 3E EO on menisci
Stage 4  evaporation on menisci: 4P passive h_e(Pc) (closure-labeled) ─> 4E one EO check
Stage 5  sweeps -> qualified property tables, per case
Stage 6  device: COMSOL effective medium (liquid supply + conduction + optional EO current)
Stage 7  cold-plate integration, only when its triggers apply
```

| Stage | Physical question | Prerequisites | Main outputs | Changed in this revision |
| --- | --- | --- | --- | --- |
| 1 | How much liquid flows through the array per unit pressure gradient with a flat, shear-free top? | geometry checks | K (superficial), mesh study | quarter cell is exact only for Stokes flow; see 1Q |
| 1Q | Is K linear over device gradients, independent of the reduced cell, and how much does an immobile interface lower it? | Stage 1 | K(G) for Stokes and Navier-Stokes; K(no-slip top)/K(shear-free) | **new**; required before tables, not before Stage 2 |
| 2 | How much superheat does conduction plus a linearized interfacial resistance produce with a flat interface? | Stage 1 geometry, trimmed pillar | T_base − T_v, Q balance, apparent h | Q_in from heater-flux integral; A_base check; closure named |
| 2A | For a given mean field in the liquid, how much superficial EO flow, current and Joule power does the array produce, and how do field treatment, floor/pillar charge and interface shear change it? | Stage 1 converged on the same case | Ũ = U_EO/(μ_EO E0), F, R = Q_uniform/Q_solved, shares | depends on Stage 1 only; detailed procedure in §7 |
| 3 | What meniscus shape and physical endpoint does each Pc produce, and how do K and EO flux change with it? | Stages 1, 2A | f(x,y;Pc), clearance, mean and local edge angle, K(Pc), F(Pc), Ũ(Pc) | local depinning criterion added; interface bracket |
| 4 | How much heat leaves a curved meniscus at a given superheat, under a named closure? | Stage 3; closure gate | h_e(Pc), Q balance | passive first; EO enters only as a single check (4E) |
| 5 | Tables for the device model | Stages 1Q, 3, 4 qualified per case | `docs/cell_table_contract.md` datasets | per case; closure and interface condition recorded |
| 6 | Where does the liquid supply fail first, at what heat load, and what are the temperatures? | Stage 5 tables; resolved device inputs | q''_dry, location, T field, mass and energy balance | COMSOL becomes primary; Python models are cross-checks |
| 7 | Do contacts, spreading, reservoirs or vapor transport change Stage 6 predictions? | Stage 6 | corrected temperatures / limits | only when triggered (§12) |

Why this order, and what was not added:

- **Stage 2A after Stage 1, parallel to Stage 2.** EO slip in the thin-EDL limit is a
  boundary condition on the same Stokes problem. It reuses Stage 1 geometry, mesh and selections
  and needs nothing from Stage 2.
- **EO does not get its own thermal stage.** At the unit cell, EO changes the heat equation only
  through Joule heating (≤ 1e-5 of q'' in `MODEL_STATUS.md`) and advection. Advection stays small
  while the vertical Péclet number is ~1e-2 and the liquid arrives near T_sat (V2019 App. A.1).
  One check at the largest planned field (4E) decides whether h_e needs a field axis. It is
  expected not to.
- **No moving interfaces, transient physics or resolved EDL by default.** Each needs a stated
  trigger: Ca and Bo stay small (V2019 App. A.1), steady operation is the target, and λ_D/half-gap
  ≤ 0.1 for every chemistry screened (`case_hand_checks.py`, Sec. E). Resolve the EDL only if
  the chosen electrolyte breaks that, or in a separate thin-film study.
- **Stage 6 is an effective medium**, not explicitly meshed pillars (≈ 32 000 pillars on the
  proposed array).

---

## 3. Conventions shared by all stages

**Parameters.** Global Definitions > Parameters, one table per `.mph`. Never use COMSOL reserved
names (`h` is the mesh element size; also `t`, `x`, `y`, `z`, `pi`, `freq`, `lambda`, `phase`).
Values (reference case, water at 100 °C from CoolProp/IAPWS):

| Name | Reference | Proposed | Notes |
| --- | --- | --- | --- |
| `d` | `10[um]` | `20[um]` | pillar diameter |
| `l` | `30[um]` | `50[um]` | center-to-center pitch |
| `hp` | `25[um]` | `40[um]` | pillar height (not `h`) |
| `tsub` | `1[um]` | `1[um]` | heated plane 1 µm below the pillar base (V2019 Sec. 2.2 convention) |
| `Tv` | `373.15[K]` | `373.15[K]` | saturated vapor (proposed: matched for comparison, ASSUMED) |
| `q0` | `50[W/cm^2]` | `50[W/cm^2]` | unit-cell heat flux |
| `gam` | `0.05892[N/m]` | same | surface tension at 100 °C |
| `mu0` | `2.816e-4[Pa*s]` | same | viscosity at 100 °C |
| `rho0` | `958.3[kg/m^3]` | same | density at 100 °C |
| `cpl` | `4216[J/(kg*K)]` | same | only needed so COMSOL accepts the fluid |
| `kl` | `0.677[W/(m*K)]` | same | liquid conductivity |
| `ks` | `148[W/(m*K)]` | same | silicon, ASSUMED (300 K textbook value from memory; 373 K value is lower; sweep) |
| `theta_r` | `15[deg]` | `15[deg]` | proposed value ASSUMED |
| `dp_dry` | `4*gam*cos(theta_r)/(d*(4/pi*(l/d)^2-1))` | same expression | V2019 Eq. 3 (2177 / 1636 Pa) |
| `hmax` | `1[um]` | `1.6[um]` | mesh-size parameter; proposed scaled by l |

Properties at another temperature (for example a room-temperature experiment) are a separate
parameter set, not an edit of this one. K, F and Ũ are geometric and do not depend on μ if μ is
consistent, but every dimensional flow, current or heat value does.

**Signs.**

- Capillary pressure: `Pc = p_v − p_liquid`, positive when the liquid is in suction.
- Flow and current: positive in +x, the direction from the feed (anode reservoir) toward the
  closed or cathode end. Integrate the x-components (`u`, `ec.Jx`) over x-planes, never
  magnitudes such as `spf.U` or `ec.normJ`.
- Heat: Q_in and Q_evap are both reported positive. Check the sign of every flux integral once
  and write the expression used into the export header.

**Normalization.** Darcy (superficial) velocity U = Q/A_ref with A_ref the full cross-section
including pillars: quarter cell `(l/2)·hp`, flow length `l/2`. Pore velocity = U/porosity is a
different quantity and never enters K. Heat fluxes are per footprint `(l/2)^2`.

**Thermal reference plane.** T_base is the **area-weighted** average temperature over the whole
pillar-base plane z = 0 (pillar footprint + liquid floor), V2019's T_e ("average temperature at
the surface (red line)", Fig. 2b; the thesis drawing puts it at the substrate top). In COMSOL,
Surface Average over a set of boundaries is ∫T dA / ∫dA, i.e. area-weighted. Always export the
selection's area too: `A_base = ∫1 dA` must equal `(l/2)^2`.

**Mesh studies.** A parameter `hmax` controls every Size node. A **Parametric Sweep** (not an
auxiliary sweep, which keeps one mesh) rebuilds the mesh for each value (COMSOL KB 1261).
At least three levels, coarse to fine. Track the physical output, not residuals: K, Ũ, F,
T_base − T_v, Q_evap, z_min, mean edge angle, q''_dry. Integrals converge faster than local
values. Point values at the pillar-top edge or the meniscus-pillar corner may never converge
(KB 1261), so never use a point maximum as a convergence metric.

Thresholds (predeclared; do not relax after seeing results):

| Use | Finest-pair relative change | Other requirements |
| --- | --- | --- |
| Preliminary learning check | ≤ 1e-2 | label PRELIMINARY |
| Qualified table / reported number | ≤ 1e-3 (`docs/cell_table_contract.md`) | ≥ 3 systematic levels; conservation checks below |

Proposed revision, **not adopted**: where a 1e-3 finest-pair change is impractical (curved
menisci near the endpoint), accept a Richardson-extrapolated discretization estimate ≤ 1e-3,
provided three levels give an observed order between 0.5 and 4. Effect on claims: the same
tolerance, but resting on an asymptotic-range assumption. Adopt it only by editing the contract
and saying so.

**Exports.** One table per study. One row per (case, run, mesh). Required identifier columns:
`hmax` (and DOF from the solver log). Raw files go to `data/comsol_exports/` unedited. The header
states model file, COMSOL version, date, dataset, selections and expressions. Name columns as the
scripts expect, or pass flags.

**Model record.** For each saved `.mph`: (1) File > Compact History, then save a text copy with
File > Save As > Model File for Java (`.java`), small enough to commit next to the README in
`comsol/unit_cell/`; (2) a Report (Intermediate) listing selections, boundary conditions, mesh
statistics and solver settings; (3) the row in `comsol/unit_cell/README.md`. If your version's
save dialog differs, record what it offers.

**Licenses.** Before starting a stage, open the Licensed and Used Products dialog (in 6.x on
Windows: File > Help) and record which modules are available. The Add Physics window only
lists licensed interfaces. Stages 1-6 below have a base-package route.

---

## 4. Geometry (all unit-cell stages)

A square array is mirror-symmetric about planes through pillar centers and planes midway between
pillars (V2019 Fig. 2a).

```text
 y
 ^  (0,l/2) ________ (l/2,l/2)
 |         |        |
 |  pillar |  liquid|            quarter cell: x, y in [0, l/2]
 |  quarter|        |            flow and mean field along +x
 |  _/     |        |
 +-------------------> x
 (0,0)            (l/2,0)
```

Which reduced cell is exact:

| Problem | Quarter cell | 1/8 triangle | Needs periodic cell in x |
| --- | --- | --- | --- |
| Stokes flow (Stage 1, 3H) | exact: x = 0 and x = l/2 are planes of constant pressure and zero tangential velocity by reversibility plus mirror symmetry | no (directed flow) | — |
| Navier-Stokes at finite Re (1Q) | **not exact**: inertia breaks fore-aft symmetry | no | yes (y-mirror planes still exact) |
| Thin-EDL EO + Laplace field (2A, 3E) | exact: the x-planes are equipotentials and the EO flow is fore-aft symmetric (linear problem) | no | — |
| Meniscus shape, uniform-load conduction (2, 3, 4) | exact | exact | — |

Build (Geometry node, 3D, unit µm or m consistently):

1. **Block** `blk_cell`: `l/2` × `l/2` × `hp`, corner (0, 0, 0).
2. **Cylinder** `cyl_full`: radius `d/2`, height `hp`, at (0, 0, 0), axis z. Three quarters of it
   lie outside the cell.
3. **Trim the pillar to the cell.** Intersection of `cyl_full` with a copy of `blk_cell` (check
   *Keep input objects* for the block), giving `pillar_q`. An equivalent route: Work Plane at z = 0,
   Circle with sector angle 90°, Extrude by `hp`.
4. Flow and EO models: **Difference** `blk_cell` minus `pillar_q` → one liquid domain.
5. Thermal models: **Block** `blk_sub` `l/2` × `l/2` × `tsub` at (0, 0, −`tsub`), then **Union** of
   `blk_cell`, `pillar_q`, `blk_sub` with *Keep interior boundaries* → three domains.
6. Build All and check with **Measure** (Geometry toolbar):
   - pillar domain volume = π d² hp / 16 (reference 490.9 µm³, proposed 3141.6 µm³);
   - liquid domain volume = (l/2)² hp − π d² hp / 16;
   - every z = 0 boundary together has area (l/2)².
   A larger pillar volume means step 3 was skipped. **The domain count alone does not catch
   this:** an untrimmed cylinder still gives three domains, but adds an insulated pillar
   extension outside the cell and extra z = 0 faces that bias T_base.

Named selections (Definitions > Selections > Explicit; build them once, reuse everywhere):
`inlet` (x = 0, liquid part), `outlet` (x = l/2), `sym_y` (y = 0 and y = l/2), `floor` (z = 0,
liquid part), `pillar_wall`, `top_liquid` (z = hp, liquid part; the flat meniscus),
`pillar_top` (z = hp, pillar part), `base_plane` (all z = 0 faces: floor + pillar footprint),
`heater` (z = −tsub). Assigning a condition to the wrong face is the most common error. The
Report lists every selection by boundary number; keep it.

---

## 5. Stage 1: passive flow, flat interface

**Question.** Permeability K of the array with a flat, shear-free liquid top at Pc = 0.

**Equations.** Steady incompressible Stokes (or laminar Navier-Stokes):
μ∇²u = ∇p, ∇·u = 0. Darcy (superficial): U = (K/μ)(Δp/L), so in the quarter cell
`K = μ Q (l/2) / ((l/2) hp Δp) = μ Q / (hp Δp)`.

**COMSOL.** Laminar Flow (`spf`) with *Neglect inertial term (Stokes flow)* checked (or Creeping
Flow if licensed). V2019 used laminar Navier-Stokes (Sec. 2.2, Re_h ~ 10¹). At Δp = 1 Pa over
`l/2` both give the same K within the mesh error; the inertia question belongs to 1Q.

| Selection | Condition |
| --- | --- |
| `inlet` | Inlet > Pressure, p0 = `dp`, *Normal flow* checked |
| `outlet` | Outlet > Pressure, p0 = 0, *Normal flow* checked, *Suppress backflow* unchecked |
| `sym_y` | Symmetry |
| `floor`, `pillar_wall` | Wall, no slip |
| `top_liquid` | Symmetry (= slip; V2019 Sec. 2.2: "no flux and no shear stress") |

Parameter `dp = 1[Pa]`. Initial values u = 0, p = 0. Mesh: Free Tetrahedral, Size max `hmax`,
min `hmax/10`; optionally 2-3 boundary layers on `floor` and `pillar_wall`. Study: Stationary
inside a Parametric Sweep `hmax = 2 1 0.5 0.25 um` (reference; proposed `3.2 1.6 0.8 0.4 um`).

**Outputs.** Surface Integration of `u` on `outlet` (Q) and on `inlet` (mass check);
Global Evaluation of `dp`; DOF.

```powershell
python analysis/permeability_from_comsol.py data/comsol_exports/ref_stage1_flow.csv --mu 2.816e-4 --case reference_d10_l30_h25
```

**Checks.** |Q_in − Q_out|/Q_out < 1e-6 (Stokes, exact discrete conservation is typical);
K within a few percent of the analytical flat-top value (reference 3.01e-11, proposed
6.38e-11 m²; `case_hand_checks.py` Sec. B), finest-pair change ≤ 1e-3.

**Status (reference).** Chase: K = 2.96979e-11 m² at hmax 0.5 µm; finest-pair change 1.02e-3
(0.7 → 0.5 µm), just above the 1e-3 target. Student (reported): K = 2.9711e-11 m² at 0.25 µm,
finest-pair change 5.05e-4; K/Q is constant at μ/(hp Δp), consistent with the same normalization.
At matching meshes the builds differ by −0.047% (2 µm), −0.0007% (1 µm) and −0.006% (0.5 µm).
Richardson estimates from the three finest levels are 2.9735e-11 (Chase) and 2.9714e-11 (student).
Free-mesh refinement is not geometrically similar, so these are indicative. The comparison is
independent-rebuild evidence for Stage 1, **conditional on the student's unverified setup**.
Neither result is qualified: one low gradient (6.7e4 Pa/m), flat top only, `.mph` not tracked.

**Proposed case.** Save As `prop_stage1_flow.mph`; change only `d`, `l`, `hp`, `hmax` list.
Expected Q ≈ 9.06e-12 m³/s at Δp = 1 Pa (analytical; COMSOL a few percent lower, as for the
reference).

**Wrong-result symptoms.** K ≈ 7.9e-12 (Chase's original export contains such "setup-error"
runs; their cause was not recorded, so record the cause of any rejected run); K a factor 2 off (full-cell area with quarter-cell flow); Q_in ≠ Q_out (a selection
that misses part of a face); K far above K_s = 4.0e-11 (ref) (top or floor not no-slip).

### 5.1 Stage 1Q: qualification before any table

Run on both cases. These are the gaps that block `permeability_stage1_flow.csv` from
qualification (`config/data_inventory.yaml`).

1. **Reduced-cell check.** Build a periodic half cell: x ∈ [−l/2, l/2], y ∈ [0, l/2], half pillar
   at the origin (trimmed). Laminar Flow with *Periodic Flow Condition* between x = −l/2 and
   x = +l/2 (pressure difference `dp`), a *Pressure Point Constraint* at one point, `sym_y`
   symmetry, slip top. Here K = μ Q l / ((l/2) hp Δp). At Δp = 1 Pa with Stokes it must equal the
   quarter-cell K within 1e-3.
2. **Gradient and inertia sweep.** Same half cell, Navier-Stokes (inertia on), auxiliary sweep of
   `dp` from 0.1 Pa up to the device gradient times `l`. Device gradient at the capillary limit:
   2 Pc_max / L (reference 8.7e5 Pa/m at L = 5 mm, so Δp ≈ 26 Pa; proposed 3.3e5-1.3e6 Pa/m for
   illustrative L = 10-2.5 mm, so Δp ≈ 16-65 Pa). Report K_app(G) = μ U/G. Predeclared: if
   max|K_app/K_Stokes − 1| > 1% the 1D table must be Q(Pc, G) or carry a Forchheimer-type
   correction (`docs/cell_table_contract.md`). Hand estimate at the limit: reference Re_K ≈ 1.7;
   proposed Re_K ≈ 2-8, Re_h ≈ 10-40. **Inertia is likely not negligible for the proposed cell.**
3. **Interface bracket.** Repeat Δp = 1 Pa with `top_liquid` = no-slip wall. Analytical estimate
   K_no-slip / K_shear-free ≈ 0.69 (ref), 0.72 (proposed). Fleming et al. (2020, Sec. 3.1,
   Figs. 4-5) measured meniscus curvature that matched a no-slip interface much better than a
   shear-free one. V2019 used shear-free. Carry both until an experiment discriminates; never
   pick one to match a target.

---

## 6. Stage 2: conduction with a flat interface

**Question.** Superheat T_base − T_v produced by conduction through substrate, pillar and
stagnant liquid plus a linearized interfacial resistance at a flat meniscus. This is a learning
and verification stage. Its h is a **flat-interface, linearized-closure, conduction-only metric**,
not V2019's h_e(Pc) and not a design prediction.

**Equations.** ∇·(k∇T) = 0 in substrate and pillar (k = `ks`) and liquid (k = `kl`, u = 0).
Meniscus: −n·q = hi (T − T_v). Heater: inward flux q0 at z = −tsub. Everything else insulated
(symmetry planes; pillar top per V2019 App. A.1, natural convection ~3 orders weaker).

**Why stagnant liquid is acceptable.** Vertical Péclet number ~1e-2 and Jakob number ~1e-2
(V2019 App. A.1; `docs/assumptions.md` A3). The liquid arrives near T_sat, so lateral advection
carries little heat. In an experiment with a **subcooled** feed this fails near the feed edge.

**Interfacial closure.** `hi` is a named input, not a property. Record which value you used:

| Closure | hi at 100 °C [W/(m²K)] | Footprint R''_int, ref [mm²K/W] | Status |
| --- | --- | --- | --- |
| reduced linearized BTE, accommodation correction C4 + 2√π(1−a)/a written from memory, a = 0.052 | 4.17e5 (the old guide's "about 4.2e5") | 2.63 | ASSUMED form (Meland & Ytrehus 2003 not held) |
| Schrage, a = 0.052 | 4.03e5 | 2.72 | different closure; never a silent substitute |
| Schrage, a = 1 | 1.51e7 | 0.07 | bound: kinetic resistance nearly removed |

The closure choice matters more through `a` than through the BTE-vs-Schrage form: going from
a = 0.052 to 1 lowers the reference superheat by ~1.3 K out of ~8 K. a = 0.052 was fitted by
V2019 to Adera et al. sample A1, measured at ~24 °C / 3 kPa with a 62% heat-loss correction
(adera2016 Sec. 4). Carrying it to 100 °C or to our surfaces is an assumption.

**COMSOL.** Heat Transfer in Solids and Fluids (`ht`). Liquid: Fluid with k = `kl`, ρ = `rho0`,
cp = `cpl`, γ = 1, velocity 0 (ρ and cp drop out of a stationary problem with zero velocity).

| Selection | Condition |
| --- | --- |
| `heater` | Heat Flux, general inward, q0 = `q0` |
| side faces | Thermal Insulation (default) |
| `pillar_top` | Thermal Insulation |
| `top_liquid` | Heat Flux, convective, h = `hi`, T_ext = `Tv` |

Initial T = `Tv + 2[K]`. Mesh: global Size max `hmax`, min `hmax/10`, curvature factor 0.3; a
boundary Size node max `hmax/2` on `top_liquid` + `pillar_top`; a domain Size node on the
substrate with max `min(hmax, tsub/3)`; an edge Size node on the pillar-top circular edge
(heat flux crowds where the liquid surface meets the insulated pillar top). No boundary layers
(no flow). Study: Parametric Sweep over `hmax` as in Stage 1. Check that the element count grows
roughly as (hmax ratio)³ between levels; if not, a Size node ignores `hmax`.

**Exports (one row per `hmax`).**

| Column | Expression | Purpose |
| --- | --- | --- |
| `hmax`, `dof` | sweep parameter; solver log | identifies the row |
| `Q_in` | Surface Integration over `heater` of the inward normal flux (sign-checked) | proves the heater selection carries q0·A |
| `Q_in_param` | Global Evaluation `q0*(l/2)^2` | the input echo; must equal `Q_in` |
| `Q_evap_bc` | Surface Integration over `top_liquid` of `hi*(T-Tv)` | implementation check; equals Q_in to solver tolerance |
| `Q_evap_flux` | Surface Integration over `top_liquid` of the outward normal total flux (`ht.ntflux`, sign-checked) | discretization-sensitive energy balance |
| `T_base` | Surface Average of T over `base_plane` | V2019 T_e definition |
| `A_base` | Surface Integration of 1 over `base_plane` | must equal (l/2)² |
| `T_lv`, `A_lv` | Surface Average of T / Integration of 1 over `top_liquid` | diagnostic only |

```powershell
python analysis/heat_transfer_from_comsol.py data/comsol_exports/ref_stage2_flat.csv --qevap-col Q_evap_flux --case reference_d10_l30_h25
```

The script refuses to write a processed table if rows carry no `hmax`, and checks `A_base`.

**Checks.** A_base = (l/2)² within 1e-6; pillar volume per §4; Q_in = Q_in_param within 1e-6;
|Q_evap_flux − Q_in|/Q_in ≤ 1e-3 at the finest mesh; finest-pair change of T_base − T_v ≤ 1e-3;
T_base − T_v inside the hand bounds (reference 2.2-18.2 K, proposed 2.4-27.2 K at the BTE hi).

**Status (reference).** Chase reports two Q_in columns (both 1.125e-4 W, the q0·A echo) and two
T_base columns, 381.19639 K and 381.21757 K (superheats 8.046 and 8.068 K, 0.26% apart). The
columns carry no mesh or selection identifiers and no Q_evap is part of the report
(`data/reported/stage2_chase_reported_columns.csv`). They are **not** treated as two meshes,
averaged, or used for an energy balance. A separately assembled file labels the rows "coarser /
finer" by inference and picks the Q_evap column that balances exactly. That labelling is not
evidence, and `data/processed/h_evap_stage2_flat_d10_l30_h25.csv` derived from it is unqualified
(`config/data_inventory.yaml`). Before reuse, open `stage2_flat_interface.mph` and record: the
`hi` value and formula, whether the pillar was trimmed (Measure), which faces form the T_base
selection (and their area), and which study/mesh produced each column.

**Next.** `docs/stage2_undergrad_assignment.md` (rebuild on the reference cell, compare with
Chase, then duplicate for the proposed cell).

---

## 7. Stage 2A: isothermal electroosmotic transport, flat interface (detailed procedure)

**Question.** For a given mean electric field E0 in the liquid (along +x), what superficial EO
velocity, current and Joule power does the pillar array produce? How do the field treatment
(solved vs uniform imposed), floor vs pillar surface charge and the interface shear condition
change that answer? The answers are dimensionless and geometric, so they hold for any ζ, μ, σ
within the thin-EDL limit.

**Prerequisites.**

1. Stage 1 converged for the same case (same geometry, selections, mesh sequence).
2. Thin-EDL screening for the intended electrolyte: `python analysis/case_hand_checks.py`,
   Sec. E. λ_D/half-gap ≤ 0.1 and the finite-EDL factor ≥ 0.9 for every chemistry listed (pure
   water at 25 °C is the worst case: λ_D ≈ 0.96 µm, 0.096 of the reference half-gap). The Dukhin
   lower bound reaches ~0.1 for pure water. Surface conduction is then not negligible; prefer a
   defined dilute electrolyte (≥ 1e-5 M) for experiments.
3. ζ is a **parameter**, not a property. No value for our oxide, electrolyte and temperature is
   held (`parameters/electroosmosis/eo_inputs.yaml`, `requires_source`). Results are reported
   per μ_EO E0.

**Equations (thin EDL, outer solution).** Following the COMSOL guidance (comsol_blog_eof_edl):
compute the Debye length first, and if it is much smaller than the geometry, replace the EDL by
the Helmholtz-Smoluchowski slip.

```text
Electric:   ∇·(σ ∇V) = 0 in the liquid;  n·J = 0 on pillar wall, floor, top, y-planes
            V = dV at x = 0 (liquid face), V = 0 at x = l/2;  E0 = dV/(l/2)
Flow:       μ ∇²u = ∇p,  ∇·u = 0  (Stokes)
Slip:       u = μ_EO (E − (E·n) n) on charged walls,  μ_EO = −ε_r ε0 ζ / μ
            (equivalently u = ε ζ ∇_T V / μ;  ζ < 0 gives flow along +E)
Top:        shear-free, impermeable (baseline) or no-slip (bracket)
x-planes:   p = 0 on both (no imposed pressure drop), normal flow
```

Never combine the slip condition with a ρ_e E body force: they represent the same physics, so
using both double-counts it. A body force needs no-slip walls and a resolved EDL
(Poisson + Nernst-Planck), which is not justified here.

**Exact result to verify against (similitude).** With uniform ζ on pillar and floor, a flat
insulating top and insulating walls, the field is z-invariant and u = μ_EO E solves the problem
exactly, including the shear-free top. Then U_sup = μ_EO E0 / F, where F is the 2D formation
factor: **Ũ·F = 1**, with F = 1.1912 (reference) and 1.2875 (proposed) from the Python 2D
solver. Any deviation beyond the mesh error is a setup error.

**COMSOL setup.**

1. Open the converged Stage 1 model, Save As `ref_stage2A_eo.mph` (proposed: `prop_stage2A_eo.mph`).
   Keep geometry, selections and mesh. Split the wall selection into `pillar_wall` and `floor`
   if it is not already.
2. Parameters to add:

   | Name | Expression | Note |
   | --- | --- | --- |
   | `E0` | `1e4[V/m]` | any value; results are normalized |
   | `dV` | `E0*l/2` | potential drop across the quarter cell |
   | `sig` | `1[S/m]` | any value; F and Ũ do not depend on it |
   | `epsr` | `55.4` | value only scales μ_EO; pujahari2022 value at 100 °C |
   | `zeta_p`, `zeta_f` | `-50[mV]` | pillar and floor ζ, ILLUSTRATIVE |
   | `mueo_p`, `mueo_f` | `-epsr*epsilon0_const*zeta_p/mu0` (and `_f`) | HS mobilities |
   | `Euni` | `1` | 1 = uniform imposed field (Case A), 0 = solved field (Case B) |
   | `dp` | `0[Pa]` | imposed pressure drop (nonzero only for the superposition run) |

3. **Electric field.** Electric Currents (`ec`) on the liquid domain (in the base package;
   fallback: Mathematics > Classical PDEs > Laplace Equation with dependent variable `V`).
   Conductivity `sig`, relative permittivity 1 (irrelevant for a stationary current problem).
   Default Electric Insulation everywhere, then Electric Potential `dV` on `inlet`, Ground on
   `outlet`.
4. **Wall field variables.** Definitions > Variables, geometric entity Boundary, selection
   `pillar_wall` + `floor`:

   ```text
   Exw = Euni*E0 + (1-Euni)*ec.Ex
   Eyw = (1-Euni)*ec.Ey
   Ezw = (1-Euni)*ec.Ez
   Enw = Exw*nx + Eyw*ny + Ezw*nz
   ```

   With the Laplace fallback use `-Vx`, `-Vy`, `-Vz` in place of `ec.Ex`, `ec.Ey`, `ec.Ez`.
5. **Flow.** Laminar Flow from Stage 1 (Stokes option on). Inlet and outlet: Pressure, p0 = `dp`
   and 0, *Normal flow*, *Suppress backflow* unchecked. `sym_y`: Symmetry. `top_liquid`:
   Symmetry (baseline).
   Wall on `pillar_wall` with the EO slip. Either:
   - (Microfluidics Module) the electroosmotic-velocity wall option, with ζ = `zeta_p`,
     ε_r = `epsr`, and electric field (`Exw`, `Eyw`, `Ezw`); or
   - (any license) Wall with a prescribed wall velocity (the *Translational velocity*/*Sliding
     wall* option your version offers; see the Wall node Help) set to the explicitly tangential
     vector:

     ```text
     mueo_p*(Exw-Enw*nx),  mueo_p*(Eyw-Enw*ny),  mueo_p*(Ezw-Enw*nz)
     ```

   Repeat for `floor` with `mueo_f`. The vector is already tangential, so the result is the same
   whether or not COMSOL projects it.
6. **Study.** Step 1 Stationary: `ec` only. Step 2 Stationary: `spf` only, with *Values of
   variables not solved for* taken from Step 1. Both problems are linear and one-way coupled.
   Wrap in a Parametric Sweep over `hmax` (Stage 1 list), and run the run list below with
   additional sweep parameters or separate studies.
7. **Mesh.** Reuse Stage 1. The EDL is not resolved (slip model), so no λ_D-scale mesh is needed.
   Keep the Stage 1 boundary layers if they exist; they matter for Case A and for the no-slip top.

**Run list (each on ≥ 3 meshes unless noted).**

| Run | Settings | Purpose / expected result |
| --- | --- | --- |
| R0 | `E0 = 0`, `dp = 1[Pa]` | K must reproduce Stage 1 (catches accidental edits) |
| R1 | Case B (`Euni = 0`), ζ_p = ζ_f, slip top | **Ũ·F = 1**; F = 1.1912 / 1.2875; ∫\|ec.Ez\| ≈ 0 |
| R2 | Case A (`Euni = 1`), ζ_p = ζ_f, slip top | R = Q_A/Q_B; between the 2D tall-pillar limit (0.543 / 0.561) and 1; heuristic bridge 0.66-0.69 |
| R3a, R3b | Case B with (ζ_p, 0) and (0, ζ_f) | pillar and floor shares; Q(R3a) + Q(R3b) = Q(R1) (linearity) |
| R4 | R1 with `top_liquid` no-slip wall | EO reduction by an immobile interface (Fleming 2020 bracket) |
| R5 | R1 with `dp = 1[Pa]` | Q = Q(R0) + Q(R1) within solver tolerance (superposition) |
| R6 | R1 with `E0 → −E0` (one mesh) | Q exactly reversed |
| R7 | R1 with `zeta_p = zeta_f = 0` (one mesh) | Q = 0 |
| R8 (optional) | R5 with inertia on and `dp` at the device gradient (1Q list) | departure from superposition when inertia matters |

**Exports (one row per run × mesh).** Header: model file, version, date, case.

| Column | Expression (quarter cell) |
| --- | --- |
| `case_id`, `run_id`, `Euni`, `zeta_p`, `zeta_f`, `top_bc`, `hmax`, `dof`, `E0`, `dV`, `dp` | parameters |
| `Q_out` | ∫ u dA over `outlet` [m³/s] |
| `Q_in` | ∫ u dA over `inlet` |
| `I_out` | ∫ ec.Jx dA over `outlet` [A], signed |
| `I_in` | ∫ ec.Jx dA over `inlet` |
| `P_joule` | ∫ ec.Qrh dV over the liquid (or ∫ sig*(Vx²+Vy²+Vz²) dV) [W] |
| `A_out` | ∫ 1 dA over `outlet`, must be (l/2)·hp |
| `wall_normal_u` | max of \|u·n\| over `pillar_wall` + `floor`, must be ~0 |

Derived in Python (by hand or a short script):

```text
U_sup  = Q_out / ((l/2) hp)
Ũ      = U_sup / (mueo E0)            (use mueo_p for R1, R2, R4-R8)
F      = sig * hp * dV / I_out        (quarter cell: F = σ A ΔV / (I L), A = (l/2) hp, L = l/2)
R      = Q_out(R2) / Q_out(R1)
```

**Verification checks (predeclared).**

| Check | Tolerance |
| --- | --- |
| Ũ·F − 1 (R1) | ≤ finest-pair change of F, and ≤ 1e-3 at the finest mesh |
| F vs Python 2D (1.1912 / 1.2875) | ≤ 1e-3 |
| \|I_in − I_out\| / I_out, \|Q_in − Q_out\| / Q_out | ≤ 1e-3 |
| P_joule − I_out·dV | ≤ 1e-3 relative (signed: both positive for dV > 0) |
| Linearity (R3a + R3b vs R1; R5 vs R0 + R1) | ≤ 1e-6 relative (same mesh) |
| Reversal and zero-ζ (R6, R7) | exact to solver tolerance |
| Mesh: finest-pair change of Ũ (R1, R2) and of R | ≤ 1e-3 |

**What a wrong result looks like.**

| Symptom | Likely cause |
| --- | --- |
| Ũ·F ≈ 2 or 0.5 | full-cell area with quarter-cell flow, or `dV` across `l` instead of `l/2` |
| Ũ < 0 with ζ < 0 | sign of μ_EO or of the potential BC swapped |
| Ũ·F < 1 in R1 by more than the mesh error | slip missing on the floor, a no-slip top, or slip applied to the top |
| R > 1 | uniform field applied without tangential projection (normal velocity into the wall) |
| F ≈ 1 | pillar not excluded from the liquid domain, or the pillar wall not insulating |
| `wall_normal_u` not ~0 | the wall-velocity option enforces the full vector including a normal part |

**Limitations.** Thin EDL with uniform σ (no surface conduction, no concentration polarization);
flat interface; isothermal; uniform ζ per surface; no electrodes. E0 is the **mean field in the
wick liquid**. It is not the supply voltage divided by a length: the supply voltage adds the
reservoir and entrance resistances and both electrode overpotentials, and DC current requires
Faradaic reactions at the electrodes (gas, pH change). Power accounting must keep the wick bulk
term (P_joule, from this stage) separate from electrode and reservoir losses (not modeled).

**Dimensional translation.** For any ζ, μ and E0 at the same geometry: U_EO = Ũ μ_EO E0, and the
mean field follows from a **measured** current: E0 = I F / (σ A). In an experiment this avoids
electrode overpotentials (§ Exp-3 in `docs/validation_and_experiments.md`). Equivalent
EO pressure gradient: μ U_EO / K (reference ~1e4 Pa/m at ζ = −50 mV, E0 = 1e4 V/m, 25 °C; device
capillary gradients are ~1e5-1e6 Pa/m).

**Handoff.** Ũ and F per case to the 1D EO model (`src/nanofab/electroosmosis/`, which assumes
Ũ = 1/F) and to Stage 3E (curved menisci, where similitude no longer holds). The 2D Python results
(`docs/eo_field_treatment_test.md`) are the tall-pillar limit that R2 must approach as hp grows.

---

## 8. Stage 3: menisci, then transport on them

**Question.** For each Pc from 0 to the physical endpoint: the meniscus shape, the liquid
clearance above the floor, the edge angle, and how K, F and the EO flux change.

**Equation.** Static normal-stress balance for z = hp + f(x, y), constant curvature over the cell
(V2019 Sec. 2.1, Eq. 2; App. A.1: Bo ~ 1e-5, Ca ~ 1e-3):

```text
∇·( ∇f / sqrt(1 + |∇f|²) ) = Pc / gam,   f = 0 on the pinned pillar-top edge,
zero normal flux on symmetry edges
```

Sign test that does not depend on any paper's normal convention: **Pc > 0 must give f < 0**
(the meniscus sags). V2019 Eq. 2 is printed with H positive for a convex surface; the PDF text
loses a minus sign, so check by the sign test, not by copying the printed form.

**COMSOL (recommended: one 3D model).**

1. In the 3D flow model, General Form Boundary PDE on `top_liquid`, dependent variable `f` [m],
   conservative flux −∇_T f / sqrt(1 + |∇_T f|²), source −`Pc/gam`. Dirichlet f = 0 on the pillar
   edge, zero flux on the symmetry edges. (A 2D General Form PDE on the planform, as in V2019, is
   equivalent; then pass f to 3D with a General Extrusion.)
2. Deform the liquid domain with the deformed-geometry / moving-mesh features (in recent 6.x
   versions these live under Component > Definitions as a Deforming Domain plus boundary
   Prescribed Mesh Displacement; older versions use a Deformed Geometry interface; check Help):
   dz = f on `top_liquid`; dz = 0 on `floor`; dx = 0 on x-planes, dy = 0 on y-planes;
   dx = dy = 0 on `pillar_wall`. Physics must evaluate on the deformed (spatial) frame; confirm
   that the liquid volume integral decreases as Pc rises.
   Sag is ~10-12% of hp at the endpoint for both cases (`case_hand_checks.py` Sec. A), so the
   mesh deformation stays mild. V2019 imported the shape as a mesh; that route is fine too.
3. Parametric continuation in Pc from 0 upward (the PDE is nonlinear), with a finer step near
   the endpoint.

**Outputs per Pc.** z_min of the surface (clearance = hp + min f), liquid volume, interface area,
the mean edge slope and the **maximum local** edge slope.

**Endpoint.** Integrating the PDE over the planform gives the exact identity
`∮_edge sin β ds = Pc · A_liquid / gam` (β = meniscus angle below horizontal at the pinned edge).
V2019 Eq. 3 is this identity with the **mean** sin β = cos θ_rec. So:

- check the identity numerically at every Pc (it verifies units, sign and selections);
- the mean criterion reproduces Eq. 3 by construction; the **local** angle reaches θ_rec earlier
  where the edge slope is steepest, so local depinning can start below Eq. 3;
- report both, and floor contact (clearance → 0; not expected for either case).

The endpoint used in Stage 6 is a recorded choice with this evidence, not "the last row".
Real pinning also depends on edge sharpness after oxidation and on contact-angle hysteresis;
both are experimental questions.

**3H: passive transport on menisci.** Stage 1 (and 1Q item 3) on the deformed geometry:
K(Pc) with a shear-free and a no-slip interface. K must fall as Pc rises (V2019 Sec. 2.2).

**3E: EO on menisci.** Stage 2A R1-R4 on the deformed geometry: F(Pc), Ũ(Pc), floor/pillar
shares. On a curved shear-free top the similitude result no longer holds exactly, so Ũ·F ≠ 1
is expected here, which is why this stage exists. Pujahari et al. (2022) treated curved menisci
with an imposed uniform field. Their numbers are context, not a benchmark for the solved field.

**Checks.** Pc = 0 gives f ≡ 0 and recovers Stage 1/2A; force identity ≤ 1e-4 relative;
convergence of z_min (V2019 Table A1: 1.1e-4), K and Ũ at the endpoint, where meshes struggle
most.

---

## 9. Stage 4: evaporation on curved interfaces

**Question.** h_e(Pc) = q0 / (T_base − T_v) for a meniscus evaporating into pure saturated vapor
under a **named** kinetic closure.

**Closure gate.**

- "V2019 replication" requires V2019's closure: linearized BTE for BKW molecules (Eqs. 5-7,
  C4* = 2.13204, d4* = 0.44675), Knudsen-layer pressure equal to the far-field pressure, and
  the Meland & Ytrehus (2003) accommodation treatment. That paper is not held and our correction
  form is from memory. **Until it is obtained and checked, no Stage 4 output is labeled a V2019
  replication** (`docs/cell_table_contract.md`: `accommodation_treatment_verified`).
- Meanwhile, run named closures as separate datasets: `reduced_BTE_assumed_MY` (the
  `comsol_bte_expression` in `src/nanofab/evaporation.py`) and `Schrage`. Never mix closures in
  one table or switch closure to match a number.
- `a` is a fitted parameter (V2019 Sec. 3.4, to Adera A1 at 3 kPa). Sweep it (0.03-1). For our own
  surfaces it can only be calibrated on our own pure-vapor data and then tested on a different
  condition.

**COMSOL.** Stage 2 thermal model on the Stage 3 deformed geometry. Interpolation function `psat`
(K → Pa) from `data/processed/psat_water.txt` (`python analysis/baseline_estimates.py
--psat-table`). Boundary variable `q_evap` = the closure expression; on `top_liquid` Heat Flux,
general inward, −`q_evap`. Initial T = Tv + 1 K; if Newton fails, ramp q0 with an auxiliary
sweep. Refine at the meniscus-pillar edge; track integrals (h_e, Q_evap), not the edge flux peak.

**4P passive first.** h_e(Pc) per case. Checks: energy balance ≤ 1e-3; h_e at Pc = 0 near
the Stage 2 value with the matching linearized hi; h_e rises with Pc (V2019 Fig. 3); a 10% change
in q0 changes h_e by ~0.1% (V2019 Sec. 3.1); linearization error (nonlinear psat vs hi) reported.

**4E one EO check.** Take the 3E flow field and current at the largest planned field, add
convection u·∇T and the Joule source σ|E|² to the liquid, and recompute h_e at two Pc values.
Predeclared: if |Δh_e/h_e| < 1e-3, h_e tables carry no field axis.

**Limitations to state with every h_e.** Pure saturated vapor only. **In open air, evaporation is
limited by vapor diffusion, and this unit-cell h_e does not apply** (V2019 Sec. 3.4 uses an
air-ambient experiment only for dryout). No disjoining pressure, no thermal Marangoni (Fleming
2020 suggests surface-tension gradients may immobilize the interface), no nucleation. Report
the maximum liquid superheat so a boiling-incipience comparison is possible.

---

## 10. Stage 5: sweeps and property tables

Per case: Pc from 0 to the recorded endpoint (denser near it); K(Pc) for both interface
conditions; h_e(Pc) per closure; F(Pc), Ũ(Pc) for EO cases; gradient study at ≥ 3 Pc values
(1Q method). Import with `analysis/import_cell_tables.py` under `docs/cell_table_contract.md`. The
manifest records the case id from `parameters/unit_cell/cases.yaml`, the interface condition
(in `boundary_conditions.hydraulic`), the closure and `a`. A table that misses a qualification
target stays a preliminary learning table; relabeling never promotes it.

V2019 replication targets (reference case only, matched geometry and conditions): Fig. 3(a)
pitch sweep and Fig. 3(b) height sweep of h_e(Pc) at q0 = 50 W/cm²; Sec. 3.1 q0 sensitivity;
Table A1 convergence levels. Digitize figures into `data/literature/` before comparing.

---

## 11. Stage 6: device model in COMSOL (effective medium)

**Question.** Where and at what heat load does liquid supply fail first, what are the substrate
temperatures, and how does EO assistance change both, with mass and energy conserved?

**Prerequisites.** Stage 5 tables for the case; resolved feed arrangement, heater footprint,
substrate thickness, reservoir levels and hydraulic BC of every reservoir, vapor environment
(for the proposed case all are UNRESOLVED today).

**Formulation (V2019 Sec. 2.3 generalized to 2D planform).** Unknowns Pc(x, y) on the wick
boundary, T in the substrate, optionally V(x, y) in the wick layer.

```text
liquid:   ∇_T·[ (ρ hp K(Pc)/μ) ∇_T Pc + ρ hp U_EO ] = −m''_evap ,  m''_evap = h_e(Pc)(T − T_v)/h_fg
          (superficial U = −(K/μ)∇p_l = (K/μ)∇Pc, so liquid flows toward higher Pc;
           U_EO = Ũ(Pc) μ_EO E_T from Stages 2A/3E; the bracket is the liquid mass flux per width)
current:  ∇_T·( (σ hp / F(Pc)) ∇_T V ) = 0        (only for EO cases)
heat:     ∇·(k_s ∇T) = 0 in the substrate;  heater flux map on the bottom/backside;
          top (pillar-base plane, wetted area): −n·q = h_e(Pc)(T − T_v) + (Joule, if EO)
```

T at the top surface is the pillar-base-plane temperature, matching the h_e definition. Do not
add a through-wick resistance on top of h_e (it already contains it; `docs/whole_wick_assumptions.md` W6).

**COMSOL (base-package route).** Heat Transfer in Solids on the substrate domain (3D, or 2D x-z
for a two-edge 1D feed as in V2019); General Form Boundary PDE for Pc on the active-array top
boundary; a second boundary PDE (or Electric Currents, Shell if licensed) for V; Interpolation
functions K(Pc), h_e(Pc), F(Pc), Ũ(Pc) from Stage 5 with **no extrapolation** (error outside the
table range). Darcy's Law (Porous Media Flow / Subsurface Flow modules) is an alternative only if
licensed; the boundary-PDE route is identical in physics.

**Boundary conditions that must be stated.**

- Feed edge: Pc = Pc_feed = ρ g (z_wick − z_reservoir) + entrance loss. A reservoir 5 mm below
  the wick costs ~47 Pa (~3% of the proposed Pc_max). The same 5 mm is comparable to or larger
  than the EO stall pressure of the wick at modest fields (§ Exp-3).
- Closed edges: zero liquid flux.
- Cathode-side edge (EO): sealed compartment (zero net flux; EO raises the end pressure, which
  can flood), vented free surface (fixed pressure; EO drives a net throughflow, which adds
  sensible-heat advection if the feed is subcooled), or gas-filled (pressure rises with gas
  evolution). Each gives a different prediction; pick one per run and record it.
- Electrodes: potential on the reservoir-side wick edges is the **bulk** potential there, not the
  supply voltage.

**Solver.** Fully coupled stationary Newton. Continuation in total heater power from a low load
with the passive solution as the initial guess. Dryout search: a Stop Condition (or Global
Equation) on `maxop(Pc) ≥ Pc_end`, with Pc_end from Stage 3, then refine the load by bisection.
Beyond onset the formulation is invalid; post-onset dry-patch modeling is a separate decision.

**Outputs.** q''_dry (total heater power at onset / heater footprint, plus peak local flux),
onset location, T_max, Pc and h_e maps, feed mass flow, and for EO the current, P_joule and the
bulk voltage drop.

**Verification.**

- Mass: feed inflow = ∫ m''_evap dA (+ cathode outflow) within 1e-4. Energy: heater power =
  ∫ q_evap dA + modeled losses within 1e-4. Latent heat is the only sink; a model where q_evap is
  prescribed, or the interface stays wet regardless of supply, cannot establish dryout.
- Analytic limits: constant K and h_e, uniform flux, one feed edge → `dryout.dryout_heat_flux`
  within 1e-3; with constant U_EO → G(Π) of `src/nanofab/electroosmosis/wick1d.py` within 1e-3.
- Cross-model: the V2019 case with identical tables against `src/nanofab/device/` (Python), and
  the planform case against `src/nanofab/whole_wick/`. Disagreements are investigated, not
  averaged. The Python models are cross-checks and have no experimental standing.
- Mesh: q''_dry and T_max converged ≤ 1e-3; cell size ≥ ~3 pitches (homogenization; W2).

**Limits other than capillary dryout** (report which one is reached first): boiling incipience
(Adera et al. 2016 report "boiling before dryout" for several conditions); a T_max limit of the
chip; flooding (Pc < 0 at a closed end under EO; `docs/electroosmosis_model.md`); and
electrical/chemical limits (gas evolution blocking pores, pH-driven ζ change, salt enrichment
toward the closed end, oxide field). Stage 6 detects the hydraulic ones. The others are flagged
with their diagnostics, not modeled.

---

## 12. Stage 7: simplified cold-plate integration (trigger-based)

Add only what a trigger requires:

| Trigger | Addition |
| --- | --- |
| heater or die smaller than the wick, or a non-1D spreading path | die, TIMs (Thin Layer feature) and copper base as solids |
| feed liquid subcooled by more than a few K, or EO throughflow | sensible-heat advection term in the wick equation and reservoir temperatures |
| vapor pressure not uniform (long plenum, noncondensables) or open-air operation | vapor-space model (uniform p_v replaced by a diffusion/flow model) |
| measured heat losses comparable to latent heat | calibrated loss boundary conditions |

The CAD assembly in `NANOFAB_Cursor_CAD_Handoff.zip` is a conceptual illustration with exaggerated
display geometry; never mesh it as simulation geometry.

---

## 13. Recognizing a wrong result (all stages)

| Observation | Check first |
| --- | --- |
| Integral does not change between meshes at all | the Size nodes do not use `hmax`, or an auxiliary sweep was used |
| K, F or Ũ off by exactly 2 or 4 | quarter vs half vs full cell area or length |
| T_base selection area ≠ (l/2)² | wrong faces, or untrimmed pillar (§4) |
| Q_evap equals Q_in to 1e-12 | you integrated the boundary-condition expression; also export the flux-based value |
| Superheat outside the hand bounds | k units, hi value, heater on the wrong face |
| h_e decreases with Pc | sign of Pc in the meniscus PDE (f > 0) |
| Result matches a literature number "too well" after a change | check that the change was not made to match it |

---

## 14. File and export checklist (every run)

1. `.mph` saved under the case-specific name; `.java` text copy (after Compact History) and a
   Report.
2. Raw export in `data/comsol_exports/`, unedited, with header lines for model, version, date,
   dataset, selections, expressions.
3. Row in `comsol/unit_cell/README.md` with status label (§0).
4. Post-processing command and its output kept (`results/` is regenerable; numbers that matter go
   into the README row).
5. Mesh study with ≥ 3 levels and DOF, and the conservation checks of the stage.
