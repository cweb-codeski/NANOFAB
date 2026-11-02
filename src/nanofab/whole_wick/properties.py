"""Cell types and their properties (the unit-cell -> whole-wick interface).

A *cell type* is one unit-cell design: pillar geometry, pillar/substrate material, zero-voltage
receding contact angle theta_0, effective meniscus radius r_eff, an optional dielectric coating
(electrowetting), and a *property recipe* saying where K and h_e come from:

    K(Pc)        permeability [m^2]                         lookup table
    h_e(Pc)      evaporative coefficient [W/m^2/K]          lookup table (vaartstra2019 Eq. 4 definition)
    k_eff        in-plane conductivity of the pillar/liquid layer [W/m/K]   scalar
    Pc_max(V)    2 gamma cos(theta_V) / r_eff [Pa]          closed form (electrowetting.py)

K and h_e do NOT depend on voltage. The meniscus is pinned at the pillar top edge (A9), so its
shape, and with it K and h_e, is set by Pc alone. The contact angle only decides when the
meniscus depins, which is Pc_max. Voltage therefore changes dryout only through Pc_max.

Recipe methods (``property_recipes`` in cell_types.yaml):

    K:     analytical_flat | analytical_curvature | constant | table
    h_e:   placeholder_linear | constant | table
    k_eff: parallel | constant

``table`` reads a long-format CSV with columns ``Pc_Pa`` and the value column (SI), for example
from ``unit_cell/build_property_tables.py``. Tables are tabulated up to the largest Pc_max any
voltage can reach (cos theta_V = 1, or the saturation cap) so electrowetting cells stay inside the
tabulated range. COMSOL tables must cover that range too; otherwise lookups are reported.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

import numpy as np

from ..geometry import check_pillar_geometry, solid_fraction
from ..params import Parameter
from ..permeability import pillar_wick_permeability
from ..provenance import MODEL_CALCULATED, REPO_ROOT, check_reference, load_sources, weakest
from .config import si
from .electrowetting import DielectricCoating, cos_theta_V, pc_max, r_eff_from_geometry
from .lookup import ExtrapolationLog, LookupTable

PC_TABLES = ("K", "h_e")


@dataclass
class FluidState:
    """Saturated liquid properties used by the whole-wick model (SI)."""

    T_sat: float
    rho_l: float
    mu_l: float
    k_l: float
    surface_tension: float
    h_lv: float
    provenance: str = "AUTHORITATIVE_DATABASE"
    source: str = "iapws"

    @classmethod
    def from_coolprop(cls, T_sat: float, fluid: str = "Water") -> "FluidState":
        from ..fluids import Fluid

        f = Fluid(fluid)
        return cls(T_sat, f.rho_l(T_sat), f.mu_l(T_sat), f.k_l(T_sat), f.surface_tension(T_sat), f.h_lv(T_sat))


@dataclass
class ScalarProperty:
    value: float
    unit: str
    provenance: str
    source_location: str | None = None


@dataclass
class CellType:
    name: str
    d: float
    l: float
    h: float
    pillar_material: str
    substrate_material: str
    theta_0: float  # zero-voltage receding contact angle [rad]
    k_pillar: float
    r_eff: float  # m
    coating: DielectricCoating | None
    tables: dict[str, LookupTable]
    k_eff: ScalarProperty
    pc_max_provenance: str
    records: dict[str, Parameter] = field(default_factory=dict)

    @property
    def phi(self) -> float:
        return solid_fraction(self.d, self.l)

    @property
    def can_actuate(self) -> bool:
        return self.coating is not None

    def pc_max0(self, gamma: float) -> float:
        return float(pc_max(gamma, math.cos(self.theta_0), self.r_eff))

    def pc_upper(self, gamma: float) -> float:
        """Largest Pc_max this cell can reach at any voltage."""
        if self.coating is None:
            return self.pc_max0(gamma)
        cap = max(math.cos(self.theta_0), min(1.0, math.cos(self.coating.saturation_angle)))
        return float(pc_max(gamma, cap, self.r_eff))


def _prov(*nodes) -> str:
    return weakest(*[n.provenance if isinstance(n, Parameter) else n for n in nodes])


def _read_table_csv(spec: dict) -> tuple[np.ndarray, np.ndarray]:
    path = REPO_ROOT / spec["file"]
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    if any(k in rows[0] for k in ("E_V_per_m", "V_V")):
        raise ValueError(
            f"{spec['file']}: voltage/field-dependent K or h_e tables are not supported. Under the "
            "electrowetting model voltage acts only through Pc_max (pinned meniscus, A9)."
        )
    pcs = np.array([float(r["Pc_Pa"]) for r in rows])
    vals = np.array([float(r[spec["column"]]) for r in rows])
    order = np.argsort(pcs)
    return pcs[order], vals[order]


def _table_meta(spec: dict, name: str) -> tuple[str, str]:
    prov, src, loc = spec.get("provenance"), spec.get("source"), spec.get("source_location")
    check_reference(prov, src, load_sources(), f"{name} table")
    if not loc:
        raise ValueError(f"{name}: a table file needs source_location (COMSOL model/export and stage)")
    return prov, loc


def build_coatings(spec: dict | None) -> dict[str, DielectricCoating]:
    out = {}
    for name, c in (spec or {}).items():
        sat = c.get("saturation_angle")
        out[name] = DielectricCoating(name, si(c["epsilon_r"]), si(c["t_dielectric"]), si(sat) if sat is not None else 0.0)
    return out


def build_cell_type(
    name: str,
    spec: dict,
    materials: dict,
    recipes: dict,
    coatings: dict[str, DielectricCoating],
    coating_specs: dict,
    fluid: FluidState,
    gamma: Parameter,
    n_pc: int = 21,
) -> CellType:
    """Build tables and the Pc_max(V) description for one cell type."""
    d, l, h = (spec[k] for k in ("pillar_diameter", "pillar_pitch", "pillar_height"))
    theta0 = spec["theta_0"]
    check_pillar_geometry(si(d), si(l), si(h))
    pmat = spec["pillar_material"]
    k_pillar = materials[pmat]["thermal_conductivity"]
    recipe = recipes[spec["property_recipe"]]
    geo = (d, l, h)

    r = spec.get("r_eff", "from_geometry")
    if r == "from_geometry":
        r_eff, r_prov = r_eff_from_geometry(si(d), si(l)), _prov(MODEL_CALCULATED, d, l)
    elif isinstance(r, Parameter):
        r_eff, r_prov = si(r), r.provenance
    else:
        raise ValueError(f"{name}: r_eff must be 'from_geometry' or a value leaf")

    coating = None
    pc_prov = _prov(MODEL_CALCULATED, theta0, gamma, r_prov)
    if spec.get("dielectric_coating"):
        cname = spec["dielectric_coating"]
        coating = coatings[cname]
        cs = coating_specs[cname]
        pc_prov = _prov(pc_prov, cs["epsilon_r"], cs["t_dielectric"], *( [cs["saturation_angle"]] if cs.get("saturation_angle") else []))

    gam = si(gamma)
    pc0 = float(pc_max(gam, math.cos(si(theta0)), r_eff))
    cap = 1.0 if coating is None else max(math.cos(si(theta0)), min(1.0, math.cos(coating.saturation_angle)))
    pc_top = float(pc_max(gam, cap if coating is not None else math.cos(si(theta0)), r_eff))
    pc_axis = np.linspace(0.0, pc_top, n_pc)
    s = pc_axis / pc0  # normalized by the zero-voltage limit; > 1 only for electrowetting cells
    tables: dict[str, LookupTable] = {}

    rk = recipe["K"]
    if rk["method"] in ("analytical_flat", "analytical_curvature"):
        K0 = pillar_wick_permeability(si(d), si(l), si(h))
        a = si(rk["reduction_at_pc_max"]) if rk["method"] == "analytical_curvature" else 0.0
        vals = K0 * (1.0 - a * s)
        extra = [rk["reduction_at_pc_max"]] if rk["method"] == "analytical_curvature" else []
        prov, loc = _prov(MODEL_CALCULATED, *geo, *extra), "permeability.pillar_wick_permeability (Sangani-Acrivos + Brinkman)"
        if np.any(vals <= 0):
            raise ValueError(f"{name}: K(Pc) placeholder becomes non-positive within the tabulated Pc range")
        axis = pc_axis
    elif rk["method"] == "constant":
        vals, prov, loc, axis = np.full(pc_axis.size, si(rk["value"])), rk["value"].provenance, rk["value"].source_location, pc_axis
    elif rk["method"] == "table":
        axis, vals = _read_table_csv(rk)
        prov, loc = _table_meta(rk, f"{name}.K")
    else:
        raise ValueError(f"{name}: unknown K method {rk['method']}")
    tables["K"] = LookupTable(f"{name}.K", "m2", {"Pc": axis}, vals, prov, None, loc)

    rh = recipe["h_e"]
    if rh["method"] == "placeholder_linear":
        vals = si(rh["h_e_at_zero_pc"]) * (1.0 + si(rh["rise_at_pc_max"]) * s)
        prov = _prov(rh["h_e_at_zero_pc"], rh["rise_at_pc_max"])
        loc, axis = "placeholder h_e0 (1 + rise Pc/Pc_max(V=0)) - assumption W18", pc_axis
    elif rh["method"] == "constant":
        vals, prov, loc, axis = np.full(pc_axis.size, si(rh["value"])), rh["value"].provenance, rh["value"].source_location, pc_axis
    elif rh["method"] == "table":
        axis, vals = _read_table_csv(rh)
        prov, loc = _table_meta(rh, f"{name}.h_e")
    else:
        raise ValueError(f"{name}: unknown h_e method {rh['method']}")
    tables["h_e"] = LookupTable(f"{name}.h_e", "W/m2/K", {"Pc": axis}, vals, prov, None, loc)

    rke = recipe["k_eff"]
    phi = solid_fraction(si(d), si(l))
    if rke["method"] == "parallel":
        k_eff = ScalarProperty(phi * si(k_pillar) + (1 - phi) * fluid.k_l, "W/m/K",
                               _prov(MODEL_CALCULATED, d, l, k_pillar, fluid.provenance),
                               "conduction.parallel_conductivity (upper bound; in-plane only, W6/W7)")
    elif rke["method"] == "constant":
        k_eff = ScalarProperty(si(rke["value"]), "W/m/K", rke["value"].provenance, rke["value"].source_location)
    else:
        raise ValueError(f"{name}: unknown k_eff method {rke['method']}")

    return CellType(
        name=name, d=si(d), l=si(l), h=si(h), pillar_material=pmat,
        substrate_material=spec.get("substrate_material", pmat), theta_0=si(theta0), k_pillar=si(k_pillar),
        r_eff=r_eff, coating=coating, tables=tables, k_eff=k_eff, pc_max_provenance=pc_prov,
        records={"d": d, "l": l, "h": h, "theta_0": theta0, "k_pillar": k_pillar},
    )


@dataclass
class PropertyMaps:
    K: np.ndarray
    h_e: np.ndarray
    k_eff: np.ndarray
    Pc_max: np.ndarray  # Pc_max(V) [Pa]
    Pc_max0: np.ndarray  # Pc_max(V = 0) [Pa]
    theta_0: np.ndarray  # rad
    theta_V: np.ndarray  # rad
    ew_saturated: np.ndarray  # Lippmann value capped by cos <= 1 or the saturation angle
    beyond_pc_max: np.ndarray  # cells whose Pc exceeded Pc_max(V) (looked up at Pc_max; they are dry)


def evaluate_properties(
    cell_types: list[CellType],
    type_map: np.ndarray,
    Pc: np.ndarray,
    V: np.ndarray,
    gamma: float,
    policy: str = "warn",
    log: ExtrapolationLog | None = None,
) -> PropertyMaps:
    """Per-cell properties at capillary pressure ``Pc`` and applied voltage ``V``.

    A cell whose Pc exceeds its Pc_max(V) is past local dryout. Its K and h_e are looked up at
    Pc_max(V) and the cell is flagged in ``beyond_pc_max``, so it is reported as dry, never as a
    valid wet state.
    """
    shape = type_map.shape
    out = {k: np.zeros(shape) for k in ("K", "h_e", "k_eff", "Pc_max", "Pc_max0", "theta_0", "theta_V")}
    sat = np.zeros(shape, bool)
    beyond = np.zeros(shape, bool)
    for idx, ct in enumerate(cell_types):
        mask = type_map == idx
        if not mask.any():
            continue
        v = V[mask]
        if ct.coating is None and np.any(v > 0):
            raise ValueError(f"voltage applied to cell type '{ct.name}', which has no dielectric coating")
        c, s = cos_theta_V(ct.theta_0, v, ct.coating, gamma)
        out["theta_0"][mask] = ct.theta_0
        out["theta_V"][mask] = np.arccos(np.clip(c, -1.0, 1.0))
        out["Pc_max"][mask] = pc_max(gamma, c, ct.r_eff)
        out["Pc_max0"][mask] = ct.pc_max0(gamma)
        out["k_eff"][mask] = ct.k_eff.value
        sat[mask] = s
        pcm = out["Pc_max"][mask]
        pc = Pc[mask]
        beyond[mask] = pc > pcm
        pc_q = np.minimum(pc, pcm)
        for p in PC_TABLES:
            out[p][mask] = ct.tables[p](policy, log, ct.name, Pc=pc_q)
    return PropertyMaps(**out, ew_saturated=sat, beyond_pc_max=beyond)


def table_summary(cell_types: list[CellType], gamma: float) -> list[dict]:
    """Provenance and range of every property, for reports."""
    from ..provenance import report_category

    rows = []
    for ct in cell_types:
        for p, t in ct.tables.items():
            rows.append({"cell_type": ct.name, "property": p, "unit": t.unit, "min": float(t.values.min()),
                         "max": float(t.values.max()),
                         "axes": ", ".join(f"{k}[{a[0]:.3g}..{a[-1]:.3g}]x{a.size}" for k, a in t.axes.items()),
                         "provenance": t.provenance, "category": report_category(t.provenance),
                         "source_location": t.source_location})
        rows.append({"cell_type": ct.name, "property": "k_eff", "unit": "W/m/K", "min": ct.k_eff.value,
                     "max": ct.k_eff.value, "axes": "-", "provenance": ct.k_eff.provenance,
                     "category": report_category(ct.k_eff.provenance), "source_location": ct.k_eff.source_location})
        rows.append({"cell_type": ct.name, "property": "Pc_max(V)", "unit": "Pa", "min": ct.pc_max0(gamma),
                     "max": ct.pc_upper(gamma), "axes": "V (closed form)" if ct.coating else "none (uncoated)",
                     "provenance": ct.pc_max_provenance, "category": report_category(ct.pc_max_provenance),
                     "source_location": "electrowetting.pc_max: 2 gamma cos(theta_V) / r_eff"})
    return rows
