# COMSOL Micropillar Unit-Cell Guide (Beginner, Staged)

Goal: reproduce the unit-cell parametric studies of **Vaartstra, Lu & Wang, IJHMT 136 (2019)
170** (called V2019 below) step by step. Check each stage against an analytical number before
adding the next piece of physics.

```text
Stage 1  Laminar flow, flat top                 -> permeability K (flat interface)
Stage 2  Conduction, flat top, linear kinetic BC -> wick thermal resistance
Stage 3  Meniscus shape from Young-Laplace       -> interface surfaces for each capillary pressure
Stage 4  Kinetic evaporation BC on the meniscus  -> h_e at one capillary pressure
Stage 5  Sweep capillary pressure and geometry   -> lookup tables K(dp, geometry), h_e(dp, geometry)
```

V2019 does Stages 3-5 directly. Stages 1-2 are extra stepping stones with flat interfaces, where
analytical checks exist.

Before you start, run:

```powershell
python analysis/baseline_estimates.py --psat-table
python analysis/vaartstra_replication.py
```

Keep the output open. Every stage below says which number to compare against.

COMSOL version assumed: 6.x (V2019 used 5.3a). Menu names may differ slightly.

---

## 0. Vocabulary (COMSOL terms used below)

| Term | Meaning |
| --- | --- |
| Component | One model geometry plus its physics and mesh. |
| Global Definitions > Parameters | Named constants (`d`, `l`, `hp`, ...) usable everywhere and sweepable. Avoid COMSOL's reserved names (`h` = mesh element size, `t`, `x`, `y`, `z`, `pi`, `freq`, `lambda`, `phase`). |
| Domain | A 3D region (liquid, pillar, substrate). Materials and volume physics live on domains. |
| Boundary | A face of a domain (an edge in 2D). Boundary conditions live on boundaries. |
| Physics interface | A set of equations, e.g. Laminar Flow (`spf`), Heat Transfer in Solids and Fluids (`ht`), General Form PDE (`g`). |
| Study | The solver setup (Stationary = steady state). |
| Derived Values | Numbers computed from a solution (integrals, averages). They fill a Table you can export. |
| Selection | A named set of domains/boundaries. Create them once, reuse them everywhere. |

Tip: create named Selections (Definitions > Selections > Explicit) for `inlet`, `outlet`,
`symmetry`, `floor`, `pillar_wall`, `pillar_top`, `meniscus`, `heater`. Assigning a condition to
the wrong face is the most common beginner error.

---

## 1. Parameters (all stages)

Global Definitions > Parameters. Values match `parameters/*.yaml` (V2019 baseline, water at
100 °C from CoolProp), so COMSOL and the Python checks use identical inputs.

| Name | Expression | Description | Source |
| --- | --- | --- | --- |
| `d` | `10[um]` | pillar diameter | V2019 Fig. 4 |
| `l` | `30[um]` | pillar pitch | V2019 Fig. 4 |
| `hp` | `25[um]` | pillar height (not `h`: COMSOL reserves `h` for the mesh element size) | V2019 Fig. 4 |
| `tsub` | `1[um]` | silicon below the pillar base in the cell | V2019 Sec. 2.2 |
| `Tv` | `373.15[K]` | saturated vapor temperature | V2019 Sec. 3 |
| `q0` | `50[W/cm^2]` | heat flux at the cell bottom | V2019 Sec. 3.1 |
| `acc` | `0.052` | accommodation coefficient | V2019 (fitted) |
| `theta_r` | `15[deg]` | receding contact angle | V2019 Sec. 3.1 |
| `gam` | `0.05892[N/m]` | surface tension | CoolProp, 100 °C |
| `mu0` | `2.816e-4[Pa*s]` | liquid viscosity | CoolProp, 100 °C |
| `rho0` | `958.3[kg/m^3]` | liquid density | CoolProp, 100 °C |
| `kl` | `0.677[W/(m*K)]` | liquid conductivity | CoolProp, 100 °C |
| `ks` | `148[W/(m*K)]` | silicon conductivity | Incropera, 300 K (V2019 doesn't state; sweep 120-148) |
| `dp_dry` | `4*gam*cos(theta_r)/(d*(4/pi*(l/d)^2-1))` | dryout capillary pressure, V2019 Eq. 3 (≈ 2177 Pa) | |
| `Pc` | `0[Pa]` | local capillary pressure (swept 0 to `dp_dry`) | |
| `dp` | `1[Pa]` | pressure drop across the cell (flow studies) | |
| `hmax` | `1[um]` | maximum mesh element size (convergence sweep) | |

Why user-defined properties instead of COMSOL's built-in Water material: the built-in material
uses its own correlations. Using the CoolProp values makes COMSOL and the Python checks agree
exactly.

---

## 2. Geometry: the unit-cell domain

A square pillar array repeats with period `l`. Planes through pillar centers and planes midway
between pillars are mirror planes (V2019 Fig. 2a), so only part of one cell is modeled:

```text
 y
 ^  (0,l/2) ________ (l/2,l/2)
 |         |        |
 |         |  liquid|
 |  pillar |        |
 |  quarter|        |
 |  _/     |        |
 +-------------------> x      flow direction for flow studies: +x
 (0,0)            (l/2,0)
```

- **Flow studies** (pressure-driven in x) need this **quarter cell**: the diagonal is not a
  symmetry plane once flow has a direction.
- **Meniscus and heat-transfer studies** are also symmetric about the diagonal, so a 1/8 cell
  (triangle) works. Use the quarter cell anyway while learning: one geometry serves every stage,
  at only 2x the cost of the triangle.

Steps (Geometry node, 3D component):

1. **Block** `blk_liquid`: Width `l/2`, Depth `l/2`, Height `hp`, Corner (0, 0, 0).
2. **Cylinder** `cyl_pillar`: Radius `d/2`, Height `hp`, Position (0, 0, 0), Axis z.
3. Flow-only models: **Booleans and Partitions > Difference**, add `blk_liquid`, subtract
   `cyl_pillar`. That leaves one liquid domain.
4. Heat-transfer models: use **Union** with **Keep interior boundaries** checked (or leave both
   objects and let *Form Union* finalize). That gives two domains, liquid and pillar, which share
   the pillar sidewall where heat passes from solid to liquid.
5. Heat-transfer models: **Block** `blk_sub`: Width `l/2`, Depth `l/2`, Height `tsub`, Corner
   (0, 0, `-tsub`). V2019 applies the heat flux 1 µm below the pillar base so that constriction
   (heat crowding into the pillar) is captured.
6. **Build All**. Check the domain count: 1 for flow models, 3 for thermal models (substrate,
   pillar, liquid).

Why cutting at symmetry planes is exact: mirror symmetry of the geometry and loading means zero
normal flux across those planes. V2019 A.1 checks this for the device: the temperature gradient
between neighboring cells was 0.19% of the evaporator-to-vapor gradient.

---

## 3. Stage 1: Permeability with a flat interface (Laminar Flow)

**Question answered:** how easily does liquid flow through the pillars? Output: K [m²].

**Physics interface:** Fluid Flow > Single-Phase Flow > **Laminar Flow** (`spf`), as in V2019.
It solves the steady incompressible Navier-Stokes equations. V2019 A.1 gives Re_h ~ 10 (h and
the superficial velocity at the feed), so the flow is laminar but inertia is not strictly
negligible at the device feed. At the small cell pressure drops used here (Re ≈ 0.5 at
dp = 1 Pa), Laminar Flow and Creeping Flow should agree. Run both once; if K differs by more than
~1%, inertia matters and you keep Laminar Flow.

**Fluid properties:** density `rho0`, dynamic viscosity `mu0` (User defined).

**Boundary conditions:**

| Boundary | Location | Condition | Meaning |
| --- | --- | --- | --- |
| `inlet` | x = 0 face | Inlet > Pressure, p0 = `dp`, check **Normal flow** | p = dp, tangential velocity 0 |
| `outlet` | x = l/2 face | Outlet > Pressure, p0 = 0, check **Normal flow** | p = 0 |
| `symmetry` | y = 0 and y = l/2 faces | Symmetry | no normal flow, zero shear |
| `floor` | z = 0 | Wall, No slip (default) | u = 0 |
| `pillar_wall` | curved pillar face | Wall, No slip (default) | u = 0 |
| `meniscus` | z = h (flat for now) | **Symmetry** | no flux, no shear (V2019 Sec. 2.2) |

The inlet and outlet planes pass through a pillar center and a mid-gap, both mirror planes, so
pressure plus normal flow there is exact for slow flow.

**Initial values:** u = 0, p = 0.

**Mesh:** Free Tetrahedral, Size > Custom, Maximum element size = `hmax`. Add Boundary Layers
(2-3 layers) on `floor` and `pillar_wall` once the basic run works.

**Study:** Stationary, plus **Parametric Sweep** over `hmax` = `2[um] 1[um] 0.5[um] 0.25[um]`
(coarse to fine) for mesh convergence.

**Results:**

1. Derived Values > **Surface Integration** on `outlet`, expression `u`, giving Q in m³/s.
2. Derived Values > **Global Evaluation**, expression `dp`.
3. Evaluate both into the same Table, then Export as CSV to `data/comsol_exports/stage1_flow.csv`.

**Check:**

```powershell
python analysis/permeability_from_comsol.py data/comsol_exports/stage1_flow.csv --mu 2.816e-4
```

Expected: Q ≈ 2.7e-12 m³/s at dp = 1 Pa, i.e. K ≈ 3.0e-11 m² (Sangani-Acrivos + Brinkman,
flat top). Agreement within 10-30% is plausible; the fine-mesh COMSOL value becomes the reference.
V2019 Table A1 reached a relative mesh change of 5.5e-4 for permeability. Aim for < 1e-3
between your two finest meshes.

**Common mistakes:** forgetting Normal flow on the inlet and outlet; using No slip on the top
(that makes a closed channel); using the pore area instead of the full cross-section (l/2)·h when
computing U (the script uses the full area, consistent with Darcy's law).

---

## 4. Stage 2: Conduction with a flat interface

**Question answered:** how much temperature rise does conduction through substrate, pillar and
liquid cause?

**Physics interface:** Heat Transfer > **Heat Transfer in Solids and Fluids** (`ht`), solving
`div(k grad T) = 0`. The liquid is a stagnant conductor (no velocity). V2019 A.1 justifies this:

- **Vertical Péclet number.** Pe_h = u_e·h/α ~ 1e-2, where u_e is the liquid velocity toward the
  interface. The baseline script prints it.
- **Jakob number.** Ja = c_p·ΔT/h_lv ~ 1e-2. Liquid enters at T_sat and barely heats up, so
  lateral advection carries almost no heat, even though the lateral Péclet number is large.

**Domains:** substrate and pillar are Solid with k = `ks`; liquid is Fluid with k = `kl` and the
velocity field set to zero.

**Boundary conditions:**

| Boundary | Condition | Meaning |
| --- | --- | --- |
| `heater` (z = -tsub) | Heat Flux, General inward heat flux, q0 = `q0` | heat entering the cell |
| side faces (x = 0, l/2; y = 0, l/2) | Thermal Insulation (default) | symmetry plane |
| `pillar_top` (z = h, pillar part) | Thermal Insulation | V2019: natural convection 2-3 orders of magnitude weaker than evaporation |
| `meniscus` (z = h, liquid part) | Heat Flux, Convective, h = `hi`, T_ext = `Tv` | linearized kinetic resistance |

`hi` is the linearized interfacial coefficient: about 4.2e5 W/(m²K) at a = 0.052 and 100 °C (see
`vaartstra_replication.py`, section 3). It is a kinetic-theory value, not a prescribed h_e. Stage
4 replaces it with the full relation.

**Initial values:** T = `Tv + 2[K]`.

**Results:** export one table with columns named:

| Column | How to compute it |
| --- | --- |
| `Q_in` | Global Evaluation of `q0*(l/2)^2` |
| `Q_evap` | Surface Integration of `-ht.ntflux` on `meniscus` |
| `T_base` | Surface Average of T on the z = 0 plane (liquid floor and pillar base) |
| `T_lv` | Surface Average of T on `meniscus` |

`T_base` is V2019's evaporator temperature T_e. It is the temperature of the flat "evaporator
surface" to which h_e is applied in their device model (Fig. 2b, red line). Check this reading
against the figure; if the red line is elsewhere, export that surface instead and pass
`--base-col`.

**Check:**

```powershell
python analysis/heat_transfer_from_comsol.py data/comsol_exports/stage2_flat.csv
```

- Energy balance error below 1%.
- R_wick between the parallel and series bounds printed by `baseline_estimates.py`.
- Kinetic fraction around 0.15 or less. V2019 reports h_e > 60 kW/(m²K) (PRIMARY_SOURCE). Our
  analytical R''_int combined with that bound puts the kinetic part at most ~16% of the total
  (MODEL_CALCULATED, not a paper value); the rest is expected to be conduction through the liquid.

---

## 5. Stage 3: Meniscus shape (V2019 Sec. 2.1)

**Question answered:** what does the liquid surface look like at each capillary pressure?

**Physics:** the static normal-stress balance (Young-Laplace for a surface z = f(x, y)):

```text
div( grad f / sqrt(1 + |grad f|^2) ) = Pc / gam
```

The left side is twice the mean curvature. With `Pc` > 0 (liquid pressure below vapor) and f = 0
at the pillar edge, f is negative inside the cell, so the meniscus sags below the pillar tops.
V2019 neglects gravity (Bo ~ 1e-5) and treats the curvature as constant over one cell (Ca ~ 1e-3).

**Model (separate 2D component):**

1. Geometry (2D): Rectangle `l/2` × `l/2` at the origin, minus a Circle of radius `d/2` at the
   origin (Difference). This is the top view of the liquid area.
2. Physics: Mathematics > PDE Interfaces > **General Form PDE**, dependent variable `f`
   (unit m).
   - Conservative flux Γ: `-fx/sqrt(1+fx^2+fy^2)`, `-fy/sqrt(1+fx^2+fy^2)`
   - Source term F: `-Pc/gam`
   - Damping/mass coefficients: 0 (stationary)
   These settings give `div(grad f/sqrt(...)) = Pc/gam`.
3. Boundary conditions: **Dirichlet** f = 0 on the quarter-circle arc (pinned at the pillar top
   edge, V2019). **Zero Flux** (the default) on the four straight edges, which are the symmetry
   lines.
4. Study: Stationary with Parametric Sweep `Pc` = `range(0, dp_dry/10, dp_dry)`. The equation is
   nonlinear, so sweeping from 0 lets each solution start from the previous one.
5. Checks:
   - At `Pc = 0`, f = 0 everywhere (flat).
   - The minimum of f (at the cell center) should stay above `-hp`. V2019 tracked z_min for mesh
     convergence (Table A1: 1.1e-4).
   - At `dp_dry`, the average interface angle at the pillar should match the receding angle,
     about 15° from the pillar wall. Equivalently, the slope is about tan(75°) ≈ 3.7. A big
     mismatch means a sign or unit error.

**Using the shape in 3D:** export f on a grid (Results > Export > Data, regular grid) for each
`Pc`. In the 3D component, load it as an Interpolation function `fint(x,y)`. Build a Parametric
Surface (x = s1, y = s2, z = `hp + fint(s1,s2)`) and use it with **Partition Objects** to split
the liquid block; delete the part above the surface. V2019 imported the meniscus as a mesh. Either
way, the curved face becomes the `meniscus` selection.

The quick alternative (rough, for learning only) is to subtract a sphere of radius
`2*gam/Pc` centered above the cell center. It is not symmetric about the x = 0 and y = 0 planes,
so use it only for small `Pc`.

**Flow with the curved surface:** repeat Stage 1 with the curved `meniscus` as a Symmetry
boundary (V2019 Sec. 2.2). K drops as `Pc` rises (stated qualitatively by V2019). For our 1D
analytical model to reproduce the paper's 82 W/cm², K would have to be ~19% below the analytical
flat value (MODEL_CALCULATED). Part of that gap may be error in the analytical K rather than
curvature, so compare your flat Stage 1 K with the analytical K before attributing the gap.

---

## 6. Stage 4: Kinetic evaporation boundary condition

**Question answered:** how much heat does the meniscus remove at a given superheat, and where?

**V2019 model:** the linearized Boltzmann transport equation (Sone 2000), with the Meland &
Ytrehus correction for accommodation coefficient a. The vapor pressure at the Knudsen-layer edge
is set to the far-field pressure. See `sone_bte_mass_flux` in `src/nanofab/evaporation.py`.

```text
q_evap = h_lv * (psat(T) - psat(Tv)) * sqrt(2/(R_s T)) / C4(a),   C4(a) = 2.13204 + 2 sqrt(pi)(1-a)/a
```

At a = 0.052 this is within ~3% of the Schrage relation; either works.

**Steps:**

1. Definitions > Functions > **Interpolation**, named `psat`, Data source: File,
   `data/processed/psat_water.txt` (written by `baseline_estimates.py --psat-table`). Argument
   unit `K`, function unit `Pa`. Extrapolation: none.
2. Component Definitions > **Variables**, geometric entity: Boundary, selection `meniscus`.
   Define `q_evap` = the "linearized BTE" expression printed by `baseline_estimates.py`.
3. On `meniscus`: **Heat Flux > General inward heat flux**, q0 = `-q_evap`.

**Initial values:** T = `Tv + 1[K]`. The condition is nonlinear (psat is exponential). If the
solver fails, ramp `q0` with an Auxiliary Sweep from 10% to 100%.

**Results:** same four columns as Stage 2, with Q_evap now the Surface Integration of `q_evap`
on `meniscus`. Plot `q_evap` on the meniscus surface: it should peak near the pillar, in the
thin-film region.

**Checks:**

- Energy balance below 1%.
- h_e = q0/(T_base − Tv) at `Pc = 0` should be close to the Stage 2 result.
- h_e must increase with `Pc` (V2019 Fig. 3: a thinner effective film and more interface area).
- Mesh convergence of h_e: V2019 Table A1 gives 1.4e-4. Refine at the meniscus-pillar edge.

---

## 7. Stage 5: Sweeps and lookup tables

Parametric Sweep over `Pc` from 0 to `dp_dry` for each geometry, using the matching meniscus from
Stage 3. Export:

- thermal: `Pc, d, l, hp, Q_in, Q_evap, T_base, T_lv` to `data/comsol_exports/stage5_thermal.csv`
- flow: `Pc, d, l, hp, dp, Q` to `data/comsol_exports/stage5_flow.csv`

```powershell
python analysis/heat_transfer_from_comsol.py data/comsol_exports/stage5_thermal.csv
python analysis/permeability_from_comsol.py  data/comsol_exports/stage5_flow.csv --mu 2.816e-4
```

**Replication targets (V2019):**

| Figure | Sweep | Fixed | What to compare |
| --- | --- | --- | --- |
| Fig. 3(a) | l, and Pc from 0 to dp_dry | d = 10, h = 25 µm, q0 = 50 W/cm² | h_e(Pc); higher h_e at smaller pitch |
| Fig. 3(b) | h, and Pc from 0 to dp_dry | d = 10, l = 30 µm | h_e drops with h; flatter vs Pc for tall pillars |
| Sec. 3.1 | q0 ± 10% | baseline | h_e changes by only ~0.08% |
| Fig. 4 | device model (Phase 6-7) | baseline, L = 5 mm, t_s = 100 µm | q''_dry = 82 W/cm²; h_e > 60 kW/(m²K) and varies ~30%; ΔT_e is 3 °C higher at the inlet |

Digitize Fig. 3 curves into `data/literature/` to compare quantitatively.

**Fair comparison rule:** change one variable at a time and record the fixed ones in the table.

---

## 8. Mesh convergence table (fill in for every stage)

V2019 Eq. (A1): U = |Y_{i+1} − Y_i| / Y_i between the mesh used and one level finer.

| Calculation | Metric | V2019 U (Table A1) | Your DOF | Your U |
| --- | --- | --- | --- | --- |
| interface shape | z_min | 1.1e-4 | | |
| cell HTC | h_e | 1.4e-4 | | |
| cell permeability | K | 5.5e-4 | | |

`nanofab.comsol_io.mesh_convergence_report` computes this column. V2019 evaluated U at the
most extreme curvature (dryout), where meshes struggle most. Do the same.

## 9. Exporting results

- Tables: Table window toolbar > Export, as `.csv` or `.txt`. The `%` header lines carry column
  names and units, which `read_comsol_table` parses.
- Give each Derived Value the column name the scripts expect (`Q_in`, `Q_evap`, `T_base`,
  `T_lv`), or pass `--qin-col` etc.
- Save each `.mph` in `comsol/unit_cell/` and log it in `comsol/unit_cell/README.md`.
