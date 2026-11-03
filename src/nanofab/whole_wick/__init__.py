"""Whole-wick reduced-order model: a 2D grid of homogenized cells fed by unit-cell property tables.

See ``whole_wick/README.md`` (equations, inputs, outputs, validation) and
``docs/whole_wick_assumptions.md`` (assumptions register and unresolved decisions).
"""

from .cases import CaseSet, compare_cases, heat_load_sweep, load_case_set, parameter_sweep
from .darcy import FeedSpec, solve_darcy
from .electrowetting import DielectricCoating, cos_theta_V, electrowetting_number, r_eff_from_geometry
from .dryout import find_dryout, margins
from .grid import Grid
from .lookup import ExtrapolationError, ExtrapolationLog, ExtrapolationWarning, LookupTable
from .model import SolverOptions, WickModel, solve_coupled
from .thermal import solve_thermal

__all__ = [
    "CaseSet",
    "DielectricCoating",
    "ExtrapolationError",
    "cos_theta_V",
    "electrowetting_number",
    "r_eff_from_geometry",
    "ExtrapolationLog",
    "ExtrapolationWarning",
    "FeedSpec",
    "Grid",
    "LookupTable",
    "SolverOptions",
    "WickModel",
    "compare_cases",
    "find_dryout",
    "heat_load_sweep",
    "load_case_set",
    "margins",
    "parameter_sweep",
    "solve_coupled",
    "solve_darcy",
    "solve_thermal",
]
