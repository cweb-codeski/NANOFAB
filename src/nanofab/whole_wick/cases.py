"""Stage A/F: build whole-wick models from parameter files, compare cases, run sweeps.

Parameter files (``parameters/whole_wick/``):

    wick.yaml        domain, grid, substrate, operating point (T_v, gamma), feeds, thermal edges, solver
    cell_types.yaml  materials, dielectric coatings, unit-cell designs, property recipes
    cases.yaml       heat-flux maps, geometry/material maps, voltage maps, cases, comparison set

``load_case_set(overrides={"wick.grid.nx.value": 80})`` edits the raw YAML before validation,
which is how sweeps and grid studies vary inputs without bypassing provenance checks. An
overridden leaf keeps its provenance label but gets a note saying it was overridden.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..params import Parameter, load_parameters
from ..provenance import AUTHORITATIVE_DATABASE, REPO_ROOT, check_reference, load_sources
from .config import Config, load_raw, resolve_tree, si
from .darcy import FeedSpec
from .dryout import DryoutResult, find_dryout, margins
from .grid import Grid
from .maps import ACTUATION, cell_type_map, heat_flux_map, voltage_map, vulnerability_mask
from .model import ELECTRICAL_POWER_STATUS, CoupledResult, SolverOptions, WickModel, solve_coupled
from .properties import CellType, FluidState, build_cell_type, build_coatings

DEFAULT_DIR = REPO_ROOT / "parameters" / "whole_wick"
FILES = {"wick": "wick.yaml", "cells": "cell_types.yaml", "cases": "cases.yaml"}


def _set_dotted(raw: dict, dotted: str, value: Any) -> None:
    keys = dotted.split(".")
    node = raw
    for k in keys[:-1]:
        node = node[int(k)] if isinstance(node, list) else node[k]
    last = keys[-1]
    node[last] = value
    if last == "value" and isinstance(node, dict):
        node["notes"] = (str(node.get("notes", "")) + f" [overridden to {value!r} by a sweep]").strip()


@dataclass
class CaseSet:
    wick: Config
    cells: Config
    cases: Config
    fluid: FluidState
    gamma: Parameter
    cell_types: list[CellType]
    directory: Path
    _mask_cache: dict = field(default_factory=dict)

    @property
    def case_names(self) -> list[str]:
        return list(self.cases.tree["cases"])

    @property
    def heat_map_names(self) -> list[str]:
        return list(self.cases.tree["heat_maps"])

    @property
    def cell_type_names(self) -> list[str]:
        return [c.name for c in self.cell_types]

    @property
    def registry(self) -> dict[str, Parameter]:
        return {**self.wick.registry, **self.cells.registry, **self.cases.registry}

    def grid(self) -> Grid:
        g, d = self.wick.tree["grid"], self.wick.tree["domain"]
        return Grid(int(si(g["nx"])), int(si(g["ny"])), si(d["length_x"]), si(d["length_y"]))

    def feeds(self, grid: Grid) -> FeedSpec:
        f = self.wick.tree["feeds"]
        edges = {e: si(v["Pc_feed"]) for e, v in (f.get("edges") or {}).items()}
        cap = f.get("max_mass_flow")
        return FeedSpec(edges=edges, max_mass_flow=si(cap) if cap is not None else None)

    def options(self) -> SolverOptions:
        s = self.wick.tree["solver"]
        return SolverOptions(
            face_average=s["face_average"],
            extrapolation_policy=s["extrapolation_policy"],
            dry_model=s["dry_model"],
            h_dry=si(s["h_dry"]),
            tol=si(s["tol"]),
            max_iter=int(si(s["max_iter"])),
            relax=si(s["relax"]),
        )

    def heat_map(self, name: str, grid: Grid) -> np.ndarray:
        spec = self.cases.tree["heat_maps"][name]
        if spec["type"] == "file":
            _check_file_provenance(spec, f"heat_maps.{name}")
        return heat_flux_map(spec, grid)

    def build(self, case: str, heat_map: str, V_override: np.ndarray | None = None) -> WickModel:
        grid = self.grid()
        feeds = self.feeds(grid)
        cdef = self.cases.tree["cases"][case]
        gspec = self.cases.tree["geometry_maps"][cdef["geometry_map"]]
        tmap = cell_type_map(gspec, grid, self.cell_type_names, list(feeds.edges))
        for idx in np.unique(tmap):
            grid.check_homogenization(self.cell_types[idx].l)
        vspec = self.cases.tree["voltage_maps"][cdef["voltage_map"]]
        if V_override is not None:
            V, actuation = V_override, "custom"
        elif vspec["type"] == "vulnerability_mask":
            V, actuation = voltage_map(vspec, grid, self.selective_mask(vspec, heat_map)), "selective"
        else:
            if vspec["type"] == "file":
                _check_file_provenance(vspec, f"voltage_maps.{cdef['voltage_map']}")
            V, actuation = voltage_map(vspec, grid), ACTUATION[vspec["type"]]
        if not np.any(V > 0):
            actuation = "none"
        sub = self.wick.tree["substrate"]
        k_sub = si(self.cells.tree["materials"][sub["material"]]["thermal_conductivity"])
        te = {e: si(v["T"]) for e, v in (self.wick.tree.get("thermal_edges") or {}).items()}
        gv = self.wick.tree["gravity_inplane"]
        return WickModel(
            name=case, grid=grid, cell_types=self.cell_types, type_map=tmap, V=V,
            q_map=self.heat_map(heat_map, grid), fluid=self.fluid, substrate_k=k_sub,
            substrate_t=si(sub["thickness"]), feeds=feeds, gamma=si(self.gamma), actuation=actuation,
            thermal_edges=te, gravity=(si(gv["gx"]), si(gv["gy"])), options=self.options(),
            provenance=self.registry,
        )

    def selective_mask(self, vspec: dict, heat_map: str) -> np.ndarray:
        """Vulnerable cells predicted by the passive reference case (decision D11)."""
        key = (vspec["reference_case"], heat_map, vspec["evaluate_at"], si(vspec["threshold"]))
        if key not in self._mask_cache:
            ref = self.build(vspec["reference_case"], heat_map)
            if not ref.is_passive:
                raise ValueError("the selective-actuation reference case must be passive (V = 0)")
            if vspec["evaluate_at"] == "dryout":
                res = find_dryout(ref).result
            elif vspec["evaluate_at"] == "nominal":
                res = solve_coupled(ref, 1.0)
            else:
                raise ValueError("evaluate_at must be 'dryout' or 'nominal'")
            ratio = margins(res).pressure_ratio
            self._mask_cache[key] = vulnerability_mask(ratio, si(vspec["threshold"]), res.darcy.feed_cells)
        return self._mask_cache[key]


def _check_file_provenance(spec: dict, where: str) -> None:
    check_reference(spec.get("provenance"), spec.get("source"), load_sources(), where)
    if not spec.get("source_location"):
        raise ValueError(f"{where}: file maps need a source_location")


def _surface_tension(op: dict, fluid: FluidState) -> Parameter:
    st = op.get("surface_tension", "from_fluid")
    if isinstance(st, Parameter):
        return st
    if st != "from_fluid":
        raise ValueError("operating.surface_tension must be 'from_fluid' or a value leaf")
    return Parameter(
        name="operating.surface_tension", value=fluid.surface_tension, si_unit="N/m",
        provenance=AUTHORITATIVE_DATABASE, source="iapws",
        source_location=f"IAPWS surface tension via CoolProp at T_v = {fluid.T_sat:.2f} K",
        raw_value=round(fluid.surface_tension * 1e3, 4), raw_unit="mN/m", file="wick.yaml",
    )


def load_case_set(directory: str | Path | None = None, overrides: dict[str, Any] | None = None) -> CaseSet:
    directory = Path(directory) if directory else DEFAULT_DIR
    raw = {k: load_raw(directory / f) for k, f in FILES.items()}
    for dotted, value in (overrides or {}).items():
        group, rest = dotted.split(".", 1)
        _set_dotted(raw[group], rest, value)
    base = load_parameters()
    cfg = {k: resolve_tree(raw[k], FILES[k], base) for k in FILES}
    op = cfg["wick"].tree["operating"]
    fluid = FluidState.from_coolprop(si(op["vapor_temperature"]), op.get("fluid", "Water"))
    gamma = _surface_tension(op, fluid)
    cfg["wick"].registry["operating.surface_tension"] = gamma
    ct = cfg["cells"].tree
    coating_specs = ct.get("dielectric_coatings") or {}
    coatings = build_coatings(coating_specs)
    cell_types = [
        build_cell_type(name, spec, ct["materials"], ct["property_recipes"], coatings, coating_specs, fluid, gamma,
                        n_pc=int(si(ct["table_resolution"]["n_pc"])))
        for name, spec in ct["cell_types"].items()
    ]
    return CaseSet(cfg["wick"], cfg["cells"], cfg["cases"], fluid, gamma, cell_types, directory)


# --------------------------------------------------------------------------- comparison


def electrowetting_summary(model: WickModel, result: CoupledResult) -> dict:
    """Contact-angle change and capillary-pressure increase over the actuated cells."""
    p = result.props
    act = model.actuated
    sel = act if act.any() else np.ones(model.grid.shape, bool)
    t_by_cell = np.array([ct.coating.t_dielectric if ct.coating else np.nan for ct in model.cell_types])[model.type_map]
    field_v_um = (model.V / t_by_cell)[act] * 1e-6 if act.any() else np.array([0.0])
    return {
        "actuation": model.actuation,
        "V_applied_V": float(model.V.max()),
        "actuated_area_fraction": float(act.mean()),
        "theta_0_deg": float(np.degrees(p.theta_0[sel]).mean()),
        "theta_V_deg": float(np.degrees(p.theta_V[sel]).mean()),
        "dtheta_deg": float(np.degrees(p.theta_V[sel] - p.theta_0[sel]).mean()) + 0.0,
        "Pc_max0_Pa": float(p.Pc_max0[sel].mean()),
        "Pc_maxV_Pa": float(p.Pc_max[sel].mean()),
        "dPc_max_Pa": float((p.Pc_max - p.Pc_max0)[sel].mean()),
        "dPc_max_pct": float(100 * ((p.Pc_max / p.Pc_max0) - 1)[sel].mean()),
        "dielectric_field_V_um": float(np.nanmax(field_v_um)),
        "ew_saturated": bool(p.ew_saturated.any()),
    }


@dataclass
class CaseSummary:
    case: str
    heat_map: str
    dryout: DryoutResult
    nominal: CoupledResult

    def row(self, ref: "CaseSummary | None" = None) -> dict:
        d, n = self.dryout, self.nominal
        m = n.model
        g = m.grid
        row = {"case": self.case, "heat_map": self.heat_map, **electrowetting_summary(m, n)}
        row.update({
            "dryout_mean_W_cm2": d.mean_heat_flux / 1e4,
            "dryout_peak_W_cm2": d.peak_heat_flux / 1e4,
            "dryout_heat_load_W": d.heat_load,
            "first_dry_x_mm": d.location[0] * 1e3,
            "first_dry_y_mm": d.location[1] * 1e3,
            "first_dry_cell_type": d.cell_type,
            "first_dry_V": d.voltage_at_cell,
            "first_dry_actuated": d.voltage_at_cell > 0,
            "required_dp_at_dryout_Pa": d.required_dp_at_cell,
            "Pc_maxV_at_dryout_Pa": d.pc_max_at_cell,
            "criterion": d.criterion,
            "nominal_mean_W_cm2": float((n.q_external * g.cell_area).sum()) / g.footprint / 1e4,
            "nominal_wet": not n.dry.any(),
            "T_max_C": n.thermal.T_max - 273.15,
            "Q_evap_W": n.thermal.energy.evaporative,
            "max_pressure_ratio_nominal": float(margins(n).pressure_ratio.max()),
            "electrical_power": ELECTRICAL_POWER_STATUS,
            "mass_rel_residual": n.darcy.mass.relative_residual,
            "energy_rel_residual": n.thermal.energy.relative_residual,
            "converged": n.converged and d.result.converged,
            "extrapolation_events": len(n.log) + len(d.result.log),
        })
        if ref is not None:
            row["dQ_dryout_W"] = d.heat_load - ref.dryout.heat_load
            row["dQ_dryout_pct"] = 100 * (d.heat_load / ref.dryout.heat_load - 1)
            row["dT_max_K"] = n.thermal.T_max - ref.nominal.thermal.T_max
        return row


def run_case(cs: CaseSet, case: str, heat_map: str) -> CaseSummary:
    model = cs.build(case, heat_map)
    return CaseSummary(case, heat_map, find_dryout(model), solve_coupled(model, 1.0))


def compare_cases(cs: CaseSet, heat_map: str, cases: list[str] | None = None, reference: str | None = None) -> list[dict]:
    """Run every case under the same heat map and operating point (Stage F)."""
    cases = cases or cs.case_names
    reference = reference or cs.cases.tree["comparison"]["reference_case"]
    summaries = {c: run_case(cs, c, heat_map) for c in [reference] + [c for c in cases if c != reference]}
    ref = summaries[reference]
    return [summaries[c].row(ref) for c in cases]


def heat_load_sweep(model: WickModel, scales: list[float]) -> list[dict]:
    rows = []
    for s in scales:
        r = solve_coupled(model, s, dry_model="onset")
        mg = margins(r)
        g = model.grid
        rows.append({
            "scale": s,
            "mean_heat_flux_W_cm2": float((r.q_external * g.cell_area).sum()) / g.footprint / 1e4,
            "T_max_C": r.thermal.T_max - 273.15,
            "max_pressure_ratio": float(mg.pressure_ratio.max()),
            "min_pressure_margin_Pa": float(np.where(r.darcy.feed_cells, np.inf, mg.pressure_margin).min()),
            "Q_evap_W": r.thermal.energy.evaporative,
            "inlet_kg_s": r.darcy.mass.inlet,
            "any_dry": bool(r.dry.any()),
        })
    return rows


def parameter_sweep(
    path: str, values: list[Any], case: str, heat_map: str, directory: str | Path | None = None, extra: dict | None = None
) -> list[dict]:
    """Vary one raw YAML entry (dotted path, e.g. 'cases.voltage_maps.global.V.value')."""
    rows = []
    for v in values:
        cs = load_case_set(directory, {**(extra or {}), path: v})
        s = run_case(cs, case, heat_map)
        rows.append({"parameter": path, "value": v, **s.row()})
    return rows


def cell_table(result: CoupledResult) -> dict[str, np.ndarray]:
    """Everything stored per cell, flattened (for CSV export)."""
    m = result.model
    X, Y = m.grid.XY
    mg = margins(result)
    p = result.props
    cols = {
        "x_m": X, "y_m": Y,
        "cell_type": np.array([c.name for c in m.cell_types], dtype=object)[m.type_map],
        "pillar_material": m.per_cell("pillar_material"),
        "substrate_material": m.per_cell("substrate_material"),
        "d_m": m.per_cell("d"), "l_m": m.per_cell("l"), "h_m": m.per_cell("h"), "r_eff_m": m.per_cell("r_eff"),
        "rho_l_kg_m3": m.rho, "mu_l_Pa_s": m.mu,
        "V_applied_V": m.V,
        "theta_0_deg": np.degrees(p.theta_0), "theta_V_deg": np.degrees(p.theta_V),
        "ew_saturated": p.ew_saturated,
        "Pc_max0_Pa": p.Pc_max0, "Pc_maxV_Pa": p.Pc_max,
        "required_dp_Pa": result.Pc,
        "pressure_margin_Pa": mg.pressure_margin, "pressure_ratio": mg.pressure_ratio,
        "supply_margin_kg_s": mg.supply_margin, "supply_ratio": mg.supply_ratio,
        "q_in_W_m2": result.q_in,
        "K_m2": p.K, "k_eff_W_mK": p.k_eff, "h_e_W_m2K": result.h_e_used,
        "T_K": result.thermal.T, "q_evap_W_m2": result.thermal.q_evap, "m_evap_kg_m2s": result.m_evap,
        "ux_m_s": result.darcy.ux, "uy_m_s": result.darcy.uy,
        "dry": result.dry,
    }
    return {k: np.asarray(v).ravel() for k, v in cols.items()}
