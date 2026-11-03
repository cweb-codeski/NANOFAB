"""Stage E: local dryout criteria and the dryout heat load.

Pressure criterion (sign convention D7). The *required* liquid-vapor pressure difference from the
Darcy solve, Pc = p_v - p_l, is compared with the *available* electrowetting-dependent limit
Pc_max(V) = 2 gamma cos(theta_V) / r_eff of the same cell:

    pressure margin  = Pc_max(V) - Pc        [Pa]   dry when <= 0
    pressure ratio   = Pc / Pc_max(V)        [-]    dry when >= 1

Liquid-supply criterion. The largest net inflow a cell can draw, with its own pressure lowered to
its limit p_v - Pc_max and its neighbors held fixed, is

    capacity = m_evap A + (Pc_max(V) - Pc) * sum(face transmissibilities)

so  supply margin = capacity - m_evap A = (Pc_max(V) - Pc) * sum(T_f)  [kg/s],
    supply ratio  = m_evap A / capacity                                 dry when >= 1.

With fixed-pressure feeds this is algebraically the pressure criterion scaled by the local
conductance (decision D3). It is reported because it shows *how much* flow headroom a cell has.
The independent supply limit is the optional finite feed capacity ``FeedSpec.max_mass_flow``:

    feed supply ratio = total inlet flow / max_mass_flow                dry when >= 1.

Dryout onset = the first heat-load scale at which any non-feed cell reaches pressure ratio 1, or
the feed supply ratio reaches 1. Found by root bracketing on the heat-load scale.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .model import CoupledResult, WickModel, solve_coupled


@dataclass
class Margins:
    pressure_margin: np.ndarray  # Pa
    pressure_ratio: np.ndarray  # -
    supply_margin: np.ndarray  # kg/s
    supply_ratio: np.ndarray  # -
    feed_supply_ratio: float  # - (nan if no finite feed capacity)


def margins(result: CoupledResult) -> Margins:
    da, props = result.darcy, result.props
    A = result.model.grid.cell_area
    pm = props.Pc_max - da.Pc
    ratio = np.divide(da.Pc, props.Pc_max, out=np.full_like(da.Pc, np.inf), where=props.Pc_max > 0)
    demand = result.m_evap * A
    smargin = pm * da.conductance_total
    capacity = demand + smargin
    sratio = np.divide(demand, capacity, out=np.full_like(demand, np.inf), where=capacity > 0)
    sratio = np.where(demand <= 0, 0.0, sratio)
    cap = result.model.feeds.max_mass_flow
    feed_ratio = da.mass.inlet / cap if cap else math.nan
    fc = da.feed_cells
    ratio = np.where(fc, 0.0, ratio)
    sratio = np.where(fc, 0.0, sratio)
    return Margins(pm, ratio, smargin, sratio, feed_ratio)


def onset_metric(result: CoupledResult) -> tuple[float, str]:
    """max(pressure ratio over cells, feed supply ratio) and which one governs."""
    mg = margins(result)
    p = float(mg.pressure_ratio.max())
    f = mg.feed_supply_ratio
    if not math.isnan(f) and f > p:
        return f, "feed_capacity"
    return p, "capillary_pressure"


@dataclass
class DryoutResult:
    scale: float
    heat_load: float  # W, external heat at dryout
    mean_heat_flux: float  # W/m^2 over the wick footprint
    peak_heat_flux: float  # W/m^2, max local external heat flux at dryout
    location: tuple[float, float]  # (x, y) of the first dry cell [m]
    cell_index: tuple[int, int]  # (j, i)
    cell_type: str
    criterion: str
    required_dp_at_cell: float  # Pa, p_v - p_l required by the flow
    pc_max_at_cell: float  # Pa, Pc_max(V) available there
    voltage_at_cell: float  # V
    pressure_margin_at_cell: float  # Pa, ~0 at onset
    supply_margin_at_cell: float  # kg/s
    result: CoupledResult

    def lines(self) -> list[str]:
        x, y = self.location
        return [
            f"dryout heat load      {self.heat_load:.4g} W (mean {self.mean_heat_flux / 1e4:.4g} W/cm^2, "
            f"peak local {self.peak_heat_flux / 1e4:.4g} W/cm^2)",
            f"first dry cell        x = {x * 1e3:.3f} mm, y = {y * 1e3:.3f} mm, cell type {self.cell_type}, "
            f"V = {self.voltage_at_cell:.3g} V",
            f"governing criterion   {self.criterion}",
            f"required p_v - p_l    {self.required_dp_at_cell:.4g} Pa vs Pc_max(V) {self.pc_max_at_cell:.4g} Pa",
            f"pressure margin there {self.pressure_margin_at_cell:.3g} Pa",
            f"supply margin there   {self.supply_margin_at_cell:.3g} kg/s",
        ]


def find_dryout(model: WickModel, rtol: float = 1e-5, max_doublings: int = 40) -> DryoutResult:
    """Heat-load scale at dryout onset (always evaluated with dry_model='onset')."""

    def f(s: float) -> float:
        return onset_metric(solve_coupled(model, s, dry_model="onset"))[0] - 1.0

    lo, hi = 0.0, 1.0
    f_hi = f(hi)
    if f_hi < 0:
        for _ in range(max_doublings):
            lo, hi = hi, 2 * hi
            f_hi = f(hi)
            if f_hi >= 0:
                break
        else:
            raise RuntimeError("no dryout found up to 2^40 x the heat map")
    else:
        if f(0.0) >= 0:
            s_star = 0.0
            return _package(model, s_star)
    s_star = brentq(f, lo, hi, rtol=rtol, xtol=1e-12)
    return _package(model, s_star)


def _package(model: WickModel, s: float) -> DryoutResult:
    res = solve_coupled(model, s, dry_model="onset")
    mg = margins(res)
    _, crit = onset_metric(res)
    j, i = np.unravel_index(int(np.argmax(mg.pressure_ratio)), model.grid.shape)
    A = model.grid.cell_area
    Q = float((res.q_external * A).sum())
    return DryoutResult(
        scale=s,
        heat_load=Q,
        mean_heat_flux=Q / model.grid.footprint,
        peak_heat_flux=float(res.q_external.max()),
        location=(float(model.grid.x[i]), float(model.grid.y[j])),
        cell_index=(int(j), int(i)),
        cell_type=model.cell_types[model.type_map[j, i]].name,
        criterion=crit,
        required_dp_at_cell=float(res.Pc[j, i]),
        pc_max_at_cell=float(res.props.Pc_max[j, i]),
        voltage_at_cell=float(model.V[j, i]),
        pressure_margin_at_cell=float(mg.pressure_margin[j, i]),
        supply_margin_at_cell=float(mg.supply_margin[j, i]),
        result=res,
    )
