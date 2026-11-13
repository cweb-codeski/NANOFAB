# Assignment: Stage 2, reference cell (then the proposed cell)

Build guide: `docs/comsol_unit_cell_guide.md` §4 (geometry) and §6 (Stage 2). Hand bounds:
`python analysis/case_hand_checks.py` (Sec. D).

## Goal

Rebuild Stage 2 (conduction with a flat interface and a linearized interfacial resistance) on the
**unchanged reference cell** (d/l/hp = 10/30/25 µm), independently of Chase's
`stage2_flat_interface.mph`, so that the two builds can be compared. Then duplicate it for the
proposed cell (20/50/40 µm).

This is a verification exercise. The resulting h = q0/(T_base − T_v) is a flat-interface,
linearized-closure learning metric. Do not compare it with V2019's h_e or present it as a design
prediction.

## Before you start

1. Get the exact `hi` value Chase used and write it in your parameter table. If it is not
   available, use `hi = 4.17e5[W/(m^2*K)]` and say so. A 1% change in `hi` moves T_base by about
   0.013 K, comparable to the 0.021 K between the two values Chase reported, so the comparison
   only means something with identical inputs.
2. Use the parameter values of guide §3 (reference column) exactly: `ks = 148`, `kl = 0.677`,
   `q0 = 50[W/cm^2]`, `tsub = 1[um]`, `Tv = 373.15[K]`.
3. Start from your converged Stage 1 geometry.

## Steps

1. Add the substrate block and **trim the pillar to the quarter cell** (guide §4 step 3). Check
   with Geometry > Measure: pillar volume 490.9 µm³. Union with *Keep interior boundaries* →
   3 domains.
2. Make the named selections, especially `heater` (z = −tsub), `base_plane` (**both** z = 0 faces:
   pillar footprint and liquid floor), `top_liquid`, `pillar_top`.
3. Physics and mesh as in guide §6. Parametric Sweep `hmax = 2 1 0.5 0.25 um` (if 0.25 µm runs out
   of memory, use `1 0.5 0.35`).
4. Export one table, one row per `hmax`, with these columns:
   `hmax, dof, Q_in, Q_in_param, Q_evap_bc, Q_evap_flux, T_base, A_base, T_lv, A_lv`
   (definitions in guide §6). Save it unedited as `data/comsol_exports/ref_stage2_flat_<initials>.csv`.

## Hand in

- The saved model `ref_stage2_flat_<initials>.mph`, plus its `.java` text copy (after File >
  Compact History) or an Intermediate Report.
- A selection record: which boundary numbers are in `heater`, `base_plane`, `top_liquid`,
  `pillar_top` (screenshot or the Report).
- The export file and the output of
  `python analysis/heat_transfer_from_comsol.py <export> --qevap-col Q_evap_flux --case reference_d10_l30_h25`.
- A convergence table: T_base − T_v and Q_evap_flux per mesh, relative change vs the next
  coarser mesh (coarser value in the denominator), DOF.

## Pass criteria (decided before you run)

| Check | Requirement |
| --- | --- |
| `A_base` | = 2.25e-10 m² within 1e-6 relative. This confirms T_base is the area-weighted average over the whole pillar-base plane |
| Pillar volume | 490.9 µm³ |
| `Q_in` vs `Q_in_param` | equal within 1e-6 (1.125e-4 W) |
| `Q_evap_bc` vs `Q_in` | equal to solver tolerance (implementation check) |
| \|Q_evap_flux − Q_in\| / Q_in | ≤ 1e-3 at the finest mesh |
| T_base − T_v, finest-pair change | ≤ 1e-3 |
| T_base − T_v | between the hand bounds 2.2 and 18.2 K |

If a check fails, report it. Do not adjust anything to make it pass.

## Comparison with Chase

Only after your table passes: compare T_base − T_v with Chase's model **at the same `hmax`**, once
Chase has identified which mesh and selection produced each of his two reported T_base values
(381.19639 K and 381.21757 K). Until then his numbers are unlabeled and must not be matched to
your rows. Report the difference at each matching mesh and at the finest mesh of each build.

## Then: the proposed cell

Save As `prop_stage2_flat_<initials>.mph`. Change **only** `d = 20[um]`, `l = 50[um]`, `hp = 40[um]`
and the `hmax` list (`3.2 1.6 0.8 0.4 um`). Keep `hi`, `ks`, `kl`, `q0`, `tsub`, `Tv` unchanged so
the comparison isolates geometry. Expected checks: pillar volume 3141.6 µm³, `A_base` = 6.25e-10 m²,
`Q_in` = 3.125e-4 W, T_base − T_v between 2.4 and 27.2 K. Run the script with
`--case proposed_d20_l50_h40`.

## Questions to answer (a few sentences each)

1. Why must `Q_in` be an integral over the heater boundary rather than `q0*(l/2)^2`?
2. Why does `Q_evap_bc` match `Q_in` almost exactly while `Q_evap_flux` does not?
3. Your apparent interfacial share (T_lv − T_v)/(T_base − T_v) depends on averaging T_lv over
   area. Why is that a convention and not a unique resistance split?
4. What would change in T_base − T_v if the accommodation coefficient were 1 instead of 0.052
   (use the hand-check table), and why does that matter for experiments?
