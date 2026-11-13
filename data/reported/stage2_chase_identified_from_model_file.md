# Stage 2 (Chase): what the model file shows about the exported columns

Inspected 2026-11-13, read-only, from `D:\NANOFAB\COMSOL\stage2_flat_interface.mph` (COMSOL 6.3.0.420,
sha256 `fe4cb24649c2e6a0477e15206138de24f9780c3de30f5eeeb3214c745753de43`). Without COMSOL, the
inspection read the model tree (`smodel.json`) and the stored action history (`dmodel.xml`) inside
the `.mph` ZIP container. No solution data were read. Raw tables were copied unmodified to
`data/comsol_exports/stage2_raw_chase/` (gitignored):

| File | sha256 |
| --- | --- |
| `Qin table.txt` | `8958413cb509ac066d029b7b489a116bb13837813bdb924d05313cf518cb0ead` |
| `Qevap table.txt` | `c499e599fbe79335b2ef101f8f5eed89b62bf36699eba7c51725f32a8a46165a` |
| `Tlv table.txt` | `6136a739eecaf2357afeb5fb96fc3719b2452f86b31825a1335779510afb250a` |
| `Tbase table.txt` | `f7e718f11b9c8f44ae28f5d54152e3dff51b12a08fc18b2d3f7825f151c8fdf8` |

## Established from the model file (stated settings)

| Item | Finding | Where |
| --- | --- | --- |
| Interfacial coefficient | `hi = 4.2e5 [W/(m^2*K)]`, convective Heat Flux on boundary 11, `Text = tv` (parameter `tv = 373.15[K]`) | Parameters; ht/hf2; history actions 195-202 |
| Heater | General inward heat flux `q0` on boundary 3 | ht/hf1; action 192 |
| Pillar trimming | Cylinder ∩ block (`int1`), then a liquid block minus the trimmed pillar (`dif2`, keep subtracted), then Union with the substrate block. The pillar **is** trimmed to the cell | geometry sequence; actions 139-165 |
| T_base selection | Surface Average over boundaries **6 and 10** (two faces) | `av2`; actions 323-325 |
| T_lv, Q_evap selection | boundary 11 | `av1`, `int1`; actions 280, 294 |
| Q_evap expression | first `-ht.ntflux` (column pair 1, negative value), then `ht.ntflux` (pairs 2-4). `+ht.ntflux` is positive (outward) in this model | actions 281, 287 |
| Q_in | Global Evaluation of `q0*(l/2)^2` (input echo) | `gev1`; action 270 |
| Mesh sweep | Parametric Sweep `hmax = 2[um] 1[um] 0.5[um] 0.25[um]`; four solutions stored in "Parametric Solutions 1" | std2/param; actions 358-375 |
| Laminar Flow | present but not solved in Study 2 (`setSolveFor spf false`) | action 206 |
| Typo | `dp_dry = 4*gam*cos(theta_r)/(d*(4/pi*(1/d)^2-1))` uses `1/d`, not `l/d`; it evaluates to 1.8e-6 Pa. Unused at Pc = 0, wrong for any Stage 3 reuse | Parameters |

## Established from the action order (which solution each column came from)

1. Actions 204-251: mesh set to physics-controlled automatic size (`autoMeshSize 3`), Solution 1 solved.
2. Actions 269-330: Q_in column 1, Q_evap column pairs 1-3, T_lv column 1, T_base column 1 evaluated
   on dataset `dset1` (Study 2/Solution 1) = **the automatic-mesh solution**.
3. Actions 331-375: user-controlled `hmax` mesh, Parametric Sweep run.
4. Actions 392-399: one more value appended to each table, still on `dset1` (the data setting was never
   changed to the parametric dataset). This is column 2 of Q_in, T_lv, T_base and pair 4 of Q_evap.

Interpretation:

- Column 1 is the **automatic mesh**, not `hmax = 2 µm`.
- Column 2 is "Study 2/Solution 1" after the sweep. In COMSOL this normally holds the last swept value
  (`hmax = 0.25 µm`); that is expected behavior, not verified here. Supporting evidence: the
  boundary-11 area implied by ∫T dA / T_lv is 7.4e-5 above the exact trimmed-cell value
  (l/2)² − πd²/16 for column 1, and 1.0e-5 for column 2.
- The `hmax = 2, 1, 0.5 µm` solutions exist in the file but were **never evaluated or exported**.
- On column 2, ∫ht.ntflux dA = 1.12500000000050e-4 W, and that value divided by ∫(T − T_v) dA gives
  exactly 4.200000e5. On this boundary `ht.ntflux` returned the Robin boundary-condition flux, so the
  4e-12 "balance" is an identity, not a discretization check.
- On column 1, ∫ht.ntflux dA = 1.126051e-4 W (+9.3e-4 vs Q_in), while hi·∫(T − T_v) dA computed from
  the same table's ∫T dA and T_lv is 1.12417e-4 W (−7.4e-4). That inconsistency is unexplained without
  opening the model in COMSOL.

## Consequence

The two T_base values (8.046 and 8.068 K superheat) are an automatic-mesh result and (probably) the
0.25 µm result. They are not a convergence pair from one sweep. A proper study needs the existing
four sweep solutions evaluated with Dataset = "Study 2/Parametric Solutions 1", parameter selection
= All, plus Q_in from the heater-flux integral and `A_base` = ∫1 dA over boundaries 6 and 10.
