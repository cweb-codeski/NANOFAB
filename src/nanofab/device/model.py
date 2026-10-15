"""Definitions shared by the passive device-scale micropillar-wick model.

Coordinates and conventions (docs/passive_device_model.md):

* Half device 0 <= x <= L, uniform across width W. x = 0 is a reservoir edge, x = L the central
  symmetry plane. The full two-edge-fed device is 2L x W; only L x W is simulated.
* Substrate -t_s <= z <= 0; pillar bases at z = 0 (the thermal reference plane), tops at z = h.
* Pc = p_v - p_liquid [Pa], positive suction. p_g = p_liquid - p_v = -Pc.
* Q(x) [m^3/s] is the total liquid flow across W; U = Q / (W h) is the superficial velocity.
* q'' values are per footprint area [W/m^2].

Solver status (what happened numerically/physically) is kept separate from qualification status
(whether the inputs are admissible for a given claim). A run can be numerically ADMISSIBLE and still
VERIFICATION ONLY because its tables are manufactured fixtures.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum


class SolverStatus(str, Enum):
    ADMISSIBLE = "ADMISSIBLE"
    PHYSICAL_CAPILLARY_LIMIT = "PHYSICAL_CAPILLARY_LIMIT"
    MENISCUS_GEOMETRY_LIMIT = "MENISCUS_GEOMETRY_LIMIT"
    TABLE_COVERAGE_LIMIT = "TABLE_COVERAGE_LIMIT"
    MISSING_VALIDATED_INPUTS = "MISSING_VALIDATED_INPUTS"
    OUTSIDE_MODEL_REGIME = "OUTSIDE_MODEL_REGIME"
    NUMERICAL_FAILURE = "NUMERICAL_FAILURE"


# Process exit codes for the CLIs. A located physical limit is the successful outcome of a dryout
# search; incomplete coverage, missing inputs, regime violations and numerical failures each have
# their own nonzero code so scripts can tell them apart. 1 = invalid case/configuration, 2 = usage
# error (argparse).
EXIT_CONFIG_ERROR = 1
EXIT_LOAD_ABOVE_LIMIT = 7  # forward solve only: the prescribed load is not admissible
EXIT_CODES = {
    SolverStatus.ADMISSIBLE: 0,
    SolverStatus.PHYSICAL_CAPILLARY_LIMIT: 0,
    SolverStatus.MENISCUS_GEOMETRY_LIMIT: 0,
    SolverStatus.NUMERICAL_FAILURE: 6,
    SolverStatus.MISSING_VALIDATED_INPUTS: 3,
    SolverStatus.OUTSIDE_MODEL_REGIME: 4,
    SolverStatus.TABLE_COVERAGE_LIMIT: 5,
}


class DataRole(str, Enum):
    PHYSICAL_CELL_CALCULATION = "PHYSICAL_CELL_CALCULATION"
    LITERATURE_DIGITIZATION = "LITERATURE_DIGITIZATION"
    ANALYTICAL_EXAMPLE = "ANALYTICAL_EXAMPLE"
    MANUFACTURED_FIXTURE = "MANUFACTURED_FIXTURE"


VERIFICATION_ONLY_ROLES = frozenset({DataRole.ANALYTICAL_EXAMPLE, DataRole.MANUFACTURED_FIXTURE})


class CaseRole(str, Enum):
    LITERATURE_REPLICATION = "LITERATURE_REPLICATION"
    VERIFICATION_FIXTURE = "VERIFICATION_FIXTURE"


class ModelError(RuntimeError):
    """Raised when the solver is asked for something outside its admissible domain."""


@dataclass(frozen=True)
class LiquidProperties:
    """Saturated-liquid (and vapor) properties at T_v, all SI. Constant over the device."""

    T_v: float  # K, vapor (saturation) temperature
    p_v: float  # Pa
    rho: float  # kg/m3, liquid
    mu: float  # Pa s
    h_fg: float  # J/kg
    sigma: float  # N/m
    cp: float  # J/(kg K)
    k: float  # W/(m K), liquid
    rho_v: float  # kg/m3, saturated vapor
    gamma_v: float  # cp/cv of saturated vapor
    R_specific: float  # J/(kg K)
    source: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class DeviceGeometry:
    """Half-device geometry. ``k_s`` is the substrate conductivity [W/(m K)]."""

    L: float  # m, edge-to-center wicking length
    W: float  # m, width (intensive outputs do not depend on it)
    h: float  # m, pillar height = reference wick thickness for A_ref = W h
    t_s: float  # m, substrate thickness
    k_s: float  # W/(m K)
    d: float | None = None  # m, pillar diameter (diagnostics only)
    l: float | None = None  # m, pillar pitch (diagnostics only)

    def __post_init__(self) -> None:
        for name in ("L", "W", "h", "t_s", "k_s"):
            v = getattr(self, name)
            if not (v > 0 and v < float("inf")):
                raise ValueError(f"DeviceGeometry.{name} must be positive and finite, got {v}")
        if self.d is not None and self.l is not None and not (0 < self.d < self.l):
            raise ValueError(f"invalid pillar geometry d={self.d}, l={self.l}")

    @property
    def footprint_half(self) -> float:
        return self.L * self.W

    @property
    def footprint_full(self) -> float:
        return 2 * self.L * self.W

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class SolverSettings:
    """Numerical settings. ``rtol`` is a scaled field-change tolerance (see solver docstring)."""

    Nx: int = 200
    Nz: int = 32
    rtol: float = 1e-10
    max_iter: int = 300
    relaxation: float = 1.0
    min_relaxation: float = 1e-6
    abs_floor_fraction: float = 1e-6  # absolute floor = fraction x characteristic scale
    continuation_fractions: tuple[float, ...] = (0.1, 0.25, 0.5, 0.75, 0.9, 0.97, 0.99, 0.999)
    bracket_rel_width: float = 5e-5
    endpoint_snap_rel: float = 1e-9

    def as_dict(self) -> dict:
        d = asdict(self)
        d["continuation_fractions"] = list(self.continuation_fractions)
        return d


@dataclass
class AdmissibleInterval:
    """Pressure interval in which every closure evaluation is supported, and what ends it."""

    Pc_lo: float
    Pc_cap: float
    mechanism: SolverStatus  # PHYSICAL_CAPILLARY_LIMIT, MENISCUS_GEOMETRY_LIMIT or TABLE_COVERAGE_LIMIT
    Pc_physical_limit: float | None
    Pc_geometry_limit: float | None
    Pc_coverage_hi: float
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        d = asdict(self)
        d["mechanism"] = self.mechanism.value
        return d
