"""Whole-wick model container and the coupled Darcy-thermal solver.

Coupling loop (Picard iteration on the capillary pressure field):

    Pc -> K(Pc), h_e(Pc), k_eff; Pc_max(V) from electrowetting (independent of Pc)
       -> thermal solve -> q''_evap = h_e (T - T_sat) -> m''_evap = q''_evap / h_lv
       -> Darcy solve -> required Pc = p_v - p_l -> relax -> repeat until max|dPc| / max(Pc_max) < tol

The applied voltage V acts only through Pc_max(V) (electrowetting.py). It does not change K, h_e or
the heat input. Electrical power, and any Joule or leakage heating, is UNRESOLVED until an
electrode and circuit model exists (decision D14), so no electrical heat is added.

Dry cells (required Pc > Pc_max(V)):
    dry_model = "onset"     cells keep evaporating and are only flagged. Used for onset detection;
                            the state past onset is outside model validity.
    dry_model = "dry_patch" flagged cells get h_e = h_dry (default 0) and stay dry (no rewetting).
                            Post-onset estimate only (assumption W20).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..params import Parameter
from .darcy import DarcyResult, FeedSpec, solve_darcy
from .grid import Grid
from .lookup import ExtrapolationLog
from .properties import CellType, FluidState, PropertyMaps, evaluate_properties
from .thermal import ThermalResult, solve_thermal

DRY_MODELS = ("onset", "dry_patch")
ELECTRICAL_POWER_STATUS = "UNRESOLVED (no electrode/circuit model)"


@dataclass
class SolverOptions:
    face_average: str = "harmonic"
    extrapolation_policy: str = "warn"
    dry_model: str = "onset"
    h_dry: float = 0.0
    tol: float = 1e-8
    max_iter: int = 200
    relax: float = 1.0


@dataclass
class WickModel:
    name: str
    grid: Grid
    cell_types: list[CellType]
    type_map: np.ndarray
    V: np.ndarray  # applied electrowetting voltage magnitude [V]
    q_map: np.ndarray  # external heat flux at load scale 1 [W/m^2]
    fluid: FluidState
    substrate_k: float
    substrate_t: float
    feeds: FeedSpec
    gamma: float | None = None  # liquid-vapor surface tension [N/m]; default fluid value at T_sat
    actuation: str = "none"  # none | global | selective | custom
    thermal_edges: dict[str, float] = field(default_factory=dict)
    gravity: tuple[float, float] = (0.0, 0.0)
    options: SolverOptions = field(default_factory=SolverOptions)
    rho: np.ndarray | None = None
    mu: np.ndarray | None = None
    provenance: dict[str, Parameter] = field(default_factory=dict)

    def __post_init__(self):
        if self.gamma is None:
            self.gamma = self.fluid.surface_tension
        if self.rho is None:
            self.rho = np.full(self.grid.shape, self.fluid.rho_l)
        if self.mu is None:
            self.mu = np.full(self.grid.shape, self.fluid.mu_l)
        for arr, nm in ((self.type_map, "type_map"), (self.V, "V"), (self.q_map, "q_map")):
            if arr.shape != self.grid.shape:
                raise ValueError(f"{nm} shape {arr.shape} != grid {self.grid.shape}")
        if np.any(self.V < 0):
            raise ValueError("V map must be >= 0 (voltage magnitude)")
        for idx in np.unique(self.type_map[self.V > 0]):
            if not self.cell_types[idx].can_actuate:
                raise ValueError(f"voltage applied to cell type '{self.cell_types[idx].name}' without a dielectric coating")

    def per_cell(self, attr: str) -> np.ndarray:
        """Per-cell array of a cell-type attribute (d, l, h, theta_0, r_eff, pillar_material, ...)."""
        values = [getattr(ct, attr) for ct in self.cell_types]
        dtype = object if isinstance(values[0], str) else float
        return np.array(values, dtype=dtype)[self.type_map]

    @property
    def is_passive(self) -> bool:
        return not np.any(self.V > 0)

    @property
    def actuated(self) -> np.ndarray:
        return self.V > 0


@dataclass
class CoupledResult:
    model: WickModel
    scale: float
    q_external: np.ndarray  # W/m^2
    q_in: np.ndarray  # W/m^2 (= q_external; no electrical heating while power is unresolved)
    props: PropertyMaps
    h_e_used: np.ndarray
    thermal: ThermalResult
    darcy: DarcyResult
    m_evap: np.ndarray  # kg/m^2/s
    dry: np.ndarray  # cells with required Pc > Pc_max(V)
    iterations: int
    converged: bool
    history: list[float]
    log: ExtrapolationLog
    electrical_power: str = ELECTRICAL_POWER_STATUS

    @property
    def Pc(self) -> np.ndarray:
        """Required liquid-vapor pressure difference p_v - p_l [Pa]."""
        return self.darcy.Pc


def solve_coupled(model: WickModel, scale: float = 1.0, log: ExtrapolationLog | None = None, dry_model: str | None = None) -> CoupledResult:
    """Solve the coupled flow/thermal problem at heat-load ``scale`` x ``model.q_map``."""
    opt = model.options
    dry_model = dry_model or opt.dry_model
    if dry_model not in DRY_MODELS:
        raise ValueError(f"dry_model must be one of {DRY_MODELS}")
    log = log if log is not None else ExtrapolationLog()
    g = model.grid
    h_map = model.per_cell("h")
    Pc = np.zeros(g.shape)
    sticky_dry = np.zeros(g.shape, bool)
    q_in = scale * model.q_map
    history: list[float] = []
    converged = False
    for it in range(1, opt.max_iter + 1):
        props = evaluate_properties(model.cell_types, model.type_map, Pc, model.V, model.gamma, opt.extrapolation_policy, log)
        if dry_model == "dry_patch":
            sticky_dry |= props.beyond_pc_max
            h_e = np.where(sticky_dry, opt.h_dry, props.h_e)
        else:
            h_e = props.h_e
        G = model.substrate_k * model.substrate_t + props.k_eff * h_map
        th = solve_thermal(g, G, h_e, q_in, model.fluid.T_sat, model.thermal_edges, None, opt.face_average)
        m = th.q_evap / model.fluid.h_lv
        da = solve_darcy(g, props.K, h_map, model.rho, model.mu, m, model.feeds, model.gravity, opt.face_average)
        change = float(np.abs(da.Pc - Pc).max()) / max(float(props.Pc_max.max()), 1.0)
        history.append(change)
        Pc = Pc + opt.relax * (da.Pc - Pc)
        if change < opt.tol:
            converged = True
            break
    dry = (da.Pc > props.Pc_max) & ~da.feed_cells
    if dry_model == "dry_patch":
        dry |= sticky_dry
    return CoupledResult(
        model=model, scale=scale, q_external=q_in, q_in=q_in, props=props, h_e_used=h_e, thermal=th,
        darcy=da, m_evap=m, dry=dry, iterations=it, converged=converged, history=history, log=log,
    )
