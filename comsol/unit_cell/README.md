# Unit-cell COMSOL models

`.mph` files are gitignored (large binaries). Keep one model file per stage and record here what
each contains, so a result can always be traced to a model state.

| File | Stage | Physics | Geometry params | Notes |
| --- | --- | --- | --- | --- |
| `stage1_flow.mph` | 1 | Laminar Flow | d, l, h from `parameters/geometry.yaml` | Permeability, flat symmetry top |
| `stage2_conduction.mph` | 2 | Heat Transfer in Solids and Fluids | + 1 um substrate | Flat top, linearized kinetic BC |
| `stage3_meniscus.mph` | 3 | General Form PDE (2D) + Laminar Flow | meniscus f(x,y) per Pc | K(Pc) |
| `stage4_evaporation.mph` | 4 | + linearized-BTE heat-flux BC | | h_e at one Pc |
| `stage5_sweep.mph` | 5 | Parametric sweep | | h_e(Pc, geometry), K(Pc, geometry) |

Step-by-step build instructions: `docs/comsol_unit_cell_guide.md`.
