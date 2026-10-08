# Literature Notes

For each paper, record: what was modeled, the geometry and fluid, the boundary conditions, the key
results that can serve as replication targets, and the values to copy into `parameters/` (with
`provenance: PRIMARY_SOURCE`, the source id from `config/sources.yaml`, and the section, figure or
table in `source_location`). Provenance classes are defined in `docs/SOURCE_AUTHORITY.md`. In
these notes, anything not stated by the paper is labeled PROJECT_DERIVED, MODEL_CALCULATED or
ASSUMED.

---

## Vaartstra, Lu & Wang (2019): primary reference (V2019)

G. Vaartstra, Z. Lu, E. N. Wang, "Simultaneous prediction of dryout heat flux and local
temperature for thin film evaporation in micropillar wicks", Int. J. Heat Mass Transfer 136
(2019) 170-177. doi:10.1016/j.ijheatmasstransfer.2019.02.074.
Local copy: `Downloads/NANOFAB PROJECT RESEARCH/research documents/research papers/`.

### Framework

1. **Meniscus shape** (Sec. 2.1). The contact line is pinned at the pillar tops. The shape comes
   from the static normal-stress balance σ ∇·n̂ = 2H (Eq. 2), solved in COMSOL with z = 0 at the
   pillar edges and symmetry at the cell boundaries. Curvature is constant per cell. Δp = 2σH
   (Eq. 1).
2. **Dryout capillary pressure** (Eq. 3). Δp_dryout = 4σ cos θ_rec / [d (4/π (l/d)² − 1)],
   algebraically identical to `capillary.pillar_max_capillary_pressure`. The paper gives the
   equation but no numeric Δp_dryout; our ~2.18 kPa at 100 °C is MODEL_CALCULATED.
3. **Heat-transfer parametric study** (Sec. 2.2, Fig. 2b).
   - Uniform q''_in is applied 1 µm below the pillar base.
   - The lateral faces and pillar tops are symmetry (insulated) boundaries.
   - Heat leaves only through the meniscus.
   - h_e = q''_in / (T_e − T_sat) (Eq. 4). PRIMARY_SOURCE (Fig. 2b caption): T_e is "the average
     temperature at the surface (red line)". In the thesis drawing (Fig. 4b, p. 20) the red line
     runs along the top of the substrate at the pillar base. The text never names the surface,
     so the pillar-base reading is PROJECT_DERIVED from the figure, but well supported.
   - Water property source and silicon conductivity are stated in neither the paper nor the thesis.
4. **Interface model** (Eqs. 5-7). Linearized BTE for BKW molecules (Sone 2000), C4* = 2.13204,
   d4* = 0.44675. p_k is set to the far-field pressure, and q''_e = ρ_k u_k h_lv. The
   accommodation coefficient a = 0.052 is fitted to Adera et al. data set A1, using the Meland &
   Ytrehus (2003) method. Implemented as `evaporation.sone_bte_mass_flux`.
5. **Flow parametric study** (Sec. 2.2, Fig. 2c). Laminar Navier-Stokes (Re_h order 10; see note below) with an
   imposed Δp across the cell. The interface is a symmetry boundary (no flux, no shear). K is
   back-calculated from Darcy's law, giving κ(Δp, geometry).
6. **Device model** (Sec. 2.3, App. A.3).
   - Darcy flow in the wick with pressure-dependent κ, with evaporation as a volumetric sink
     Q_m = q''_e/(h_lv h).
   - Conduction in the substrate (thickness t_s) with h_e(Δp(x)) on the evaporator surface.
   - The two are iterated until both converge to 1e-6.
   - Symmetry at x = 0 and x = L, with liquid supplied at T_sat.

### Key assumptions (App. A.1, with characteristic values l = h = 10 µm, L = 1 cm, 100 W/cm², ΔT = 10 K, 100 °C)

| Group | Value | Conclusion |
| --- | --- | --- |
| Bo | ~1e-5 | gravity negligible for the meniscus |
| We | ~1e-2 | inertia negligible for the meniscus |
| Ca | ~1e-3 | curvature constant within a cell |
| Re_h | ~10 | laminar (not creeping), so they use full Navier-Stokes |
| Pe_h (vertical) | ~1e-2 | conduction-only heat transfer in the liquid |
| Ja | ~1e-2 | negligible liquid superheat, hence negligible lateral convection |
| Ma (vapor) | ~1e-3 | near-equilibrium evaporation, so the linearized BTE is valid |
| R_nat / R_int | 1e2-1e3 | natural convection from pillar tops negligible |

All exponents in this table are PRIMARY_SOURCE. The IJHMT PDF text drops minus signs (Bo prints
as "105", Table A1 as "1.1e4") and renders Re_h ambiguously as "10 1". The thesis page images
settle it: "Re_h = ρuh/μ ~ 10¹" (vaartstra2019_thesis Sec. 2.2.1 p. 23, Appendix p. 51). This
agrees with our MODEL_CALCULATED ≈ 16.

### Results to replicate (also in `data/literature/vaartstra2019_targets.csv`)

- **Fig. 3.** h_e vs Δp at q''_in = 50 W/cm²:
  - (a) l sweep at d = 10, h = 25 µm;
  - (b) h sweep at d = 10, l = 30 µm.
  h_e rises with Δp and falls with h. A ±10% change in q''_in shifts h_e by only 0.082%.
- **Fig. 4.** Baseline device: d = 10, l = 30, h = 25 µm, L = 5 mm, t_s = 100 µm, 100 °C.
  q''_dryout = 82 W/cm², h_e > 60 kW/(m²K) and varies ~30%, ΔT_e is 3 °C higher at the inlet.
- **Fig. 5.** h = 5-50 µm with d/h = 0.4 and l/h = 1.2. q''_dryout ∝ h², while h_e decreases
  with h.
- **Fig. 6.** Validation cases:
  - Adera A1 (d = 5, l = 12, h = 82 µm) and A2 (d = 12, l = 20, h = 90 µm): 3 kPa, 650 µm
    substrate, fed on all sides.
  - Zhu B1 (d = 7, l = 20, h = 20 µm) and B2 (d = 6, l = 50, h = 19 µm): air ambient,
    q''_dry vs L.
- **Table A1.** Mesh convergence: interface z_min 1.1e-4, cell h_e 1.4e-4, cell K 5.5e-4,
  device p_min 2.3e-6, device ΔT 9.9e-7 (minus signs restored; see above).
- **Cross-check (zhou_adera2026, PRIMARY_SOURCE).** An independent model gives 84 W/cm² for the
  same geometry at 1 atm.

### What our analytical models give (`analysis/vaartstra_replication.py`)

Everything in this subsection is MODEL_CALCULATED or PROJECT_DERIVED. None of it is reported by V2019.

- **Dryout.** MODEL_CALCULATED flat-meniscus q''_dry = 100.6 W/cm², against the PRIMARY_SOURCE
  value of 82 W/cm² (+22.7%). For our 1D model to give 82 W/cm², K would have to be ~19% lower
  (MODEL_CALCULATED). V2019 does not state a K reduction. It states qualitatively that curvature
  lowers K (Sec. 2.2). Attributing the 19% to curvature alone is a PROJECT_DERIVED hypothesis.
  Other contributors are the accuracy of our analytical K (Sangani-Acrivos + Brinkman,
  coefficients unverified), the paper's non-uniform evaporation, and property differences.
  COMSOL Stage 3 (flat-meniscus K) separates the first two.
- **Height scaling.** Our model gives exactly h² by geometric similarity (MODEL_CALCULATED);
  V2019 states "scaling with h²" (PRIMARY_SOURCE, Sec. 3.3).
- **Interface resistance.** MODEL_CALCULATED R''_int ≈ 2.6 mm²K/W at a = 0.052 (linearized BTE,
  ASSUMED flat interface area). Combined with the PRIMARY_SOURCE bound h_e > 60 kW/m²K, the
  kinetic share is at most ~16% of the total unit-cell resistance (MODEL_CALCULATED; the paper
  gives no resistance breakdown).
- **Schrage vs BTE.** They agree within ~3% at a = 0.052 and are ~14% apart at a = 1 (MODEL_CALCULATED).

### Open items

- [ ] Digitize Fig. 3(a), 3(b), 5, 6(a), 6(b), A2 into `data/literature/` (PROJECT_DERIVED,
      note the figure).
- [x] T_e surface: pillar-base plane per thesis Fig. 4(b) drawing (PROJECT_DERIVED, figure only).
- [x] Re_h exponent: +1 (thesis p. 23, 51).
- [x] Silicon conductivity: stated in neither paper nor thesis. We use 148 W/m/K
      (AUTHORITATIVE_DATABASE); vaartstra2020 uses 153 W/m/K.
- [x] Water property source: stated in neither (vaartstra2020's own Fig. 8 model uses REFPROP).
- [ ] Check Meland & Ytrehus (2003) for the exact accommodation correction form.
- [ ] Contact-line treatment: V2019 does not mention a film cut-off or disjoining pressure. Check
      whether h_e is mesh-converged without one.

## Analytical models used in `src/nanofab/`

Source ids refer to `config/sources.yaml`; `config/model_provenance.yaml` records the provenance
of each model element. "Held" means we have the PDF and checked the relevant part.

| Topic | Source id(s) | Checked? | Used in |
| --- | --- | --- | --- |
| Pillar P_c,max | vaartstra2019 Eq. 3 (held); xiao2010 (origin, not held) | Eq. 3 yes | `capillary.py` |
| Array permeability | sangani_acrivos1982, drummond_tahir1984 | no (citation from memory) | `permeability.py` |
| Brinkman height correction | brinkman1949, adera2016 (via vaartstra2019 Sec. 1) | no | `permeability.py` |
| Capillary-limited dryout | faghri_heatpipe; adera2016, zhu2016 (via vaartstra2019 Sec. 1) | no | `dryout.py` |
| BTE interface model | vaartstra2019 Eqs. 5-7 (held); sone2000; meland_ytrehus2003 | Eqs. yes; correction form no | `evaporation.py` |
| Schrage relation | schrage1953 (year discrepancy, 1953 vs 1954); carey_lvpcp | no | `evaporation.py` |
| Accommodation coefficient range | vaartstra2020 Sec. III C (held) | yes | `operating_conditions.yaml` |
| Spreading resistance | lee1995_spreading | no | `conduction.py` |
| Water properties | iapws via coolprop | spot checks only | `fluids.py` |

Removed: "Byon et al., IJHMT (2014)" for the Brinkman correction. It was written from memory and
could not be verified; the held papers cite byon_kim2011 (J. Micromech. Microeng.) instead.

## Architecture inspiration (for Phase 5, not yet modeled)

| Architecture | Bottleneck it addresses | Notes |
| --- | --- | --- |
| Periodic micropillars | baseline | well characterized, easy to simulate |
| Fiber / membrane wicks | interconnected pores, high K·Pc | |
| Inverse opals | high Pc and area | hydraulic resistance, fabrication |
| Embedded / distributed liquid feeds | long feed length L (q''_dry ∝ 1/L²) | strongest lever in the baseline numbers |
| Hierarchical wicks | separate transport (large channels) from pumping (small pores) | |
