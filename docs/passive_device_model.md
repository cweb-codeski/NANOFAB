# Passive Device-Scale Micropillar-Wick Model

Code: `src/nanofab/device/`. Commands: `analysis/run_passive_device.py`,
`analysis/device_benchmarks.py`, `analysis/import_cell_tables.py`. Table contract:
`docs/cell_table_contract.md`. Assumption register: `config/device_assumptions.yaml`.
Implementation status and evidence: `docs/passive_device_implementation_report.md`.

The model follows the device framework of Vaartstra, Lu & Wang, IJHMT 136 (2019) 170 (V2019),
Sec. 2.3, Eq. (9) and App. A.3. Cell properties K(Pc) and h_e(Pc) come from tables computed on
the pinned-meniscus unit cell (V2019 Sec. 2.1-2.2). No electro-osmosis, electrode physics, graded
geometry, boiling or moving interface is implemented. The transport interface (`Closure` in
`transport.py`) takes Pc and returns tau and h_e; a later provider can add terms without changing
the solver.

Equation labels: **V2019** (stated in the paper), **REDUCTION** (algebraically equivalent
rearrangement), **PROJECT_DERIVED** (our derivation, used for verification).

## 1. Geometry, conventions, units

- Half device 0 <= x <= L, uniform across width W (invariant in y). x = 0 is a reservoir edge,
  x = L the central symmetry plane. The full two-edge-fed device is 2L x W (footprint 2LW); only
  LW is simulated. Heat fluxes are never doubled; only integrated heat and feed rates are doubled
  for full-device totals.
- Substrate -t_s <= z <= 0; pillar bases at z = 0 (thermal reference plane); pillar tops at z = h.
- p_g = p_liquid - p_v; **Pc = p_v - p_liquid = -p_g > 0 is suction**. V2019 writes
  Delta p = p_liquid - p_0 = 2 sigma H (negative curvature for a concave meniscus); here
  Pc = -2 sigma H. V2019 Eq. (3) gives the positive limit
  `Pc_limit = sigma pi d cos(theta_rec) / (l^2 - pi d^2/4)`.
- Q(x) [m^3/s] is total liquid flow across W; A_ref = W h; U = Q/A_ref is the superficial
  velocity. An intrinsic section velocity Q/(W A_liquid/b) is reported only when the meniscus
  table supplies A_liquid_section; otherwise it is reported as UNAVAILABLE. Volume-porosity
  "pore" velocities are labelled as estimates.
- SI throughout; q'' values are per footprint area.

## 2. Equations

Liquid transport (V2019 Eq. 8, depth-integrated; REDUCTION):

    Q/W = (tau(Pc)/mu) dPc/dx,     tau = h K_bulk(Pc)  [m^3]

Mass (V2019 Eq. 9 with j_e = q_e''/h_fg; REDUCTION): `dQ/dx = -W q_e''(x) / (rho h_fg)`.
Boundary conditions (V2019 Sec. 2.3): `Pc(0) = Pc_feed` (0 for the literature case; a documented
offset is allowed in verification cases), `Q(L) = 0`.

Substrate conduction (V2019 Sec. 2.3, App. A.3):

    div(k_s grad T) = 0,   q_c = -k_s grad T
    bottom  z = -t_s:  q_c . n = -q_in''(x)
    top     z = 0:     q_c . n =  q_e''(x) = h_e(Pc(x)) (T_e(x) - T_v)      (V2019 Eq. 4)
    x = 0, L:          q_c . n = 0  (insulated; the inlet does not clamp T)

h_e is the base-plane-to-vapor coefficient of the cell; it already contains pillar/liquid
conduction and the interfacial kinetics, so no extra series resistance is added. T_e is the
top-face (pillar-base plane) temperature and T_source the bottom-face temperature.

Kirchhoff transform (PROJECT_DERIVED): `Psi(Pc) = integral tau dPc`, so `Q/W = (1/mu) dPsi/dx`.
For piecewise-constant q_e'' on cells, Q is piecewise linear and Psi at faces and centres follows by
exact integration; Pc = Psi^-1 (exact for piecewise-linear tau: a quadratic per segment, solved in a
cancellation-free form).

Evaporation-moment identity (PROJECT_DERIVED, any q_e'' >= 0 with Q(L) = 0):

    integral_{Pc(0)}^{Pc(L)} tau/mu dPc = integral_0^L x q_e''(x) / (rho h_fg) dx

Closed forms (PROJECT_DERIVED; `analytic.py`). Constant K, constant h_e, uniform
q_in'' = q_e'' = q, insulated ends, saturated inlet:

    Q(x) = W q (L - x) / (rho h_fg)
    p_g(x) = p_g(0) - mu q (L x - x^2/2) / (rho K h h_fg)
    q_dry = 2 rho K h h_fg (Pc_limit - Pc_feed) / (mu L^2)
    T_e = T_v + q / h_e,   T_source = T_e + q t_s / k_s

Pressure-dependent tau, prescribed uniform evaporation:
`q_dry = (2 rho h_fg / L^2) integral_{Pc_feed}^{Pc_limit} tau/mu dPc`.
The uniform-q forms are local identities only because the constant-property thermal problem is
one-dimensional. The moment identity also shows why any "upper bound" argument between two models
needs matching evaporation distributions.

## 3. Numerics

- One axial grid (Nx cells) shared by hydraulics and conduction, so evaporation is transferred
  exactly. Cell-centred finite volumes in x-z (Nz rows); unknowns theta = T - T_v; sparse direct
  solve (SciPy). The top Robin condition uses the face temperature through the half-cell
  resistance `dz/(2 k_s)`: `q_e = theta_c / (1/h_e + dz/(2 k_s))`; T_e and T_source are face
  reconstructions, never cell-centre values.
- Coupling (V2019 App. A.3): start from uniform evaporation; iterate
  q_e -> Q -> Pc -> h_e(Pc) -> conduction -> q_e. Stop when both scaled field changes
  `||new - old|| / (||new|| + sqrt(N) floor)` are below rtol (default 1e-10; V2019 used 1e-6),
  with floor = 1e-6 x characteristic scale (Pc: admissible interval width; q: mean input flux,
  1 W/m^2 at zero load). q_e is always the conduction result, never replaced by q_in.
- Admissibility: tables are never evaluated outside [Pc_lo, Pc_cap]. If a trial evaporation field
  would need Pc > Pc_cap, the step is backtracked toward the last admissible field (the end
  suction is linear in q_e through the moment identity). Repeated backtracking is
  NUMERICAL_FAILURE, never a bracket. Condensation (q_e < 0) is OUTSIDE_MODEL_REGIME, never clipped.
- Interpolation: positive-preserving piecewise-linear; exact at knots; no extrapolation or
  clamping; a narrow endpoint snap of 1e-9 x (table pressure span) absorbs roundoff only. The
  supported interval is the intersection of the hydraulic, thermal and meniscus ranges, further
  limited by the last admissible meniscus state.

## 4. Physical-limit search and statuses

Admissible cap: a **demonstrated** meniscus geometry limit below the physical limit
(MENISCUS_GEOMETRY_LIMIT); else the physical capillary limit if the tables reach it
(PHYSICAL_CAPILLARY_LIMIT; V2019 Eq. 3 for the literature case); else the end of table coverage
(TABLE_COVERAGE_LIMIT: incomplete coverage, **not** dryout; a final CSV row is never taken as a
physical endpoint).

Dryout (V2019 Sec. 3.3, "minimum heat flux at which Delta p_dryout was reached") is solved directly
with the end condition Pc(L) = Pc_cap. For fixed h_e the conduction problem is linear in load, so
the load follows from the moment identity each iteration; on convergence the load is rescaled
exactly so the returned field has the target moment. Every table evaluation lies in
[Pc_feed, Pc_cap]. The result is certified by forward continuation from 10% of the endpoint load
through 99.9% and (1 - 5e-5): end suction must increase monotonically and stay below the cap. The
reported bracket is [last converged admissible forward load, endpoint load].

| Status | Meaning | CLI exit |
| --- | --- | --- |
| ADMISSIBLE | converged inside the admissible interval | 0 |
| PHYSICAL_CAPILLARY_LIMIT | limit located at the physical capillary limit | 0 |
| MENISCUS_GEOMETRY_LIMIT | limit located at an earlier demonstrated meniscus limit | 0 |
| TABLE_COVERAGE_LIMIT | tables end below the physical limit; lower operating point only | 5 |
| MISSING_VALIDATED_INPUTS | required qualified inputs absent; nothing solved | 3 |
| OUTSIDE_MODEL_REGIME | solution leaves supported ranges (gradient, load) or condenses | 4 |
| NUMERICAL_FAILURE | nonconvergence, nonfinite fields, unverifiable endpoint | 6 |

Exit 1 = invalid case file, 2 = usage error, 7 = a prescribed forward load at or above the located
limit. Solver status is separate from qualification status: a fixture run can be numerically
ADMISSIBLE and still VERIFICATION ONLY.

## 5. Regime diagnostics (every run)

From the solution (`diagnostics.regime_diagnostics`): Re (superficial; volume-porosity pore
estimate; gap; flat hydraulic diameter; intrinsic section if available), Ca, Bo, We, maximum
dPc/dx and the cell pressure change over the nonzero capillary scale (Pc_cap - Pc_feed; the local
Pc is never a divisor), u_e, Pe_vertical, axial Pe (superficial and pore, length L), Ja, omitted
sensible/latent and local axial-enthalpy estimates (liquid temperature approximated by T_e),
vapor Mach, cell-to-cell temperature change over local superheat. Knudsen is reported as null:
the molecular diameter is unverified. These are calculated estimates, not rigorous thresholds.
Only the closure's documented supported ranges (pressure gradient, local heat flux) stop a run.

Audit reference scales at 82 W/cm^2 (legacy flat geometry; reproduced by
`device_benchmarks.py --suite reference_scales` to within 0.1%): U_feed 0.07584 m/s, pore
velocity 0.08309 m/s, Re_h 6.45 / 7.07, Pe_vertical 0.0566, Pe_L 2479 (pore-velocity based;
the superficial value is 2262), Ja 0.0187 for 10 K. Laminar is not creeping: V2019 states Re_h of
order 10 and uses laminar Navier-Stokes, so the cell tables must include a gradient sweep
(<= 1% K variation; otherwise Q(Pc, gradient) data are needed).

## 6. Verification gates (predeclared; `nanofab.device.verification.GATES`)

| Check | Gate |
| --- | --- |
| Closed forms (constant closure): q_dry, p_g(x), Q(x), T_e, T_source; tau-quadrature dryout | <= 1e-9 rel (stricter than the 1e-4 / 1e-8 targets) |
| Global mass / energy imbalance (normalized by actual balance terms or a documented nonzero scale) | <= 1e-8 |
| Local mass / energy residual; Darcy cell consistency (independent quadrature) | <= 1e-7 |
| W invariance of intensive outputs; mirrored full-device accounting | <= 1e-10 rel |
| Zero load: Pc, Q, T - T_v | <= 1e-12 abs |
| Manufactured-table hydraulics vs independent ODE integration (same input) | <= 1e-9 rel |
| Hydraulics with continuous q(x) and manufactured x-z conduction: observed order | 1.8-2.2; finest thermal error <= 1e-3 K |
| Coupled refinement (Nx, Nz) = (50, 8), (100, 16), (200, 32): finest-pair change of dryout, end suction, feed flow, superheats | <= 0.2 % and <= 0.02 K; changes must shrink systematically |
| Solver tolerance 1e-8 vs 1e-10 | <= 1e-4 rel, <= 0.001 K |
| Limit bracket | <= 1e-4 rel |
| Moment identity (independent quadrature of tau) | <= 1e-9 rel |
| Status injections: coverage, geometry limit, gradient range, load range, NaN, condensation | distinct statuses and exit codes |

Table-knot sensitivity (manufactured tables at 100 vs 200 Pa spacing) is reported, not gated; it
is distinct from device-mesh convergence and from cell-FEM convergence.

## 7. Literature comparison (separate from verification)

Only after the checks in `parameters/cases/vaartstra2019_fig4.yaml` (geometry, feed topology,
properties, vapor conditions, h_e definition and plane, kinetic closure, load dependence, pressure
normalization, dryout criterion) and only for a located physical limit with qualified physical
tables. 82 W/cm^2 is a reported numerical result: COMPARISON ONLY, never an input; 5% is a
predeclared investigation threshold, not a validation tolerance; no experimental uncertainty is
assigned to it. Fig. 4 profiles need digitizations with documented procedure and uncertainty (not
held). The ~30% HTC variation and ~3 K inlet-centre difference in the V2019 prose lack a fully
specified normalization and are not treated as exact scalar targets. a, theta, K, h_e and k_s are
never tuned toward agreement.

## 8. Future audit items: Pujahari, DasGupta & Bhattacharya, Langmuir 38 (2022) 8442-8455 (P2022)

Listed for a later audit only (from the project brief; the article PDF is held but was not read
for this task, and the Supporting Information is not held). Nothing below is implemented.

1. P2022 Sec. 2.1.4, Eqs. (13)-(16) fix cell mass loss from the local imposed heat input and then
   compute substrate temperature; its SI flowchart does not show V2019's evaporation-feedback
   iteration. That coupling is not substituted for V2019's.
2. Its mid-plane flow area and velocity conventions differ from ours; the velocity integral near
   Eq. (13) needs area normalization to have units of velocity.
3. P2022 Sec. 2 and Eq. (1) give conflicting pressure-sign descriptions; we use Pc = p_v - p_liquid.
4. Its Schrage closure, electrolyte properties, field assumptions and power estimates are not
   passive V2019 inputs.
5. SI Sec. 3 gives 500 V for 1e5 V/m across 5 mm; a blanket "electrolysis negligible" claim is
   not inherited. SI Sec. 2's heating arithmetic needs an independent check.

## 9. What is and is not established

- **Numerically verified for these cases**: the solver reproduces the closed forms and manufactured
  solutions above and meets every gate on ANALYTICAL_EXAMPLE and MANUFACTURED_FIXTURE inputs.
- **Approximate analytical estimate**: the legacy flat-meniscus model (`dryout.py`,
  `analysis/vaartstra_replication.py`, 100.6 W/cm^2). It is not a bound on the coupled result.
- **Literature comparison**: not yet possible (no qualified physical cell tables).
- **Experimentally validated**: nothing in this model.
