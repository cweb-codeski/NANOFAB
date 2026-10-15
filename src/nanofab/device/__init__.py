"""Passive device-scale micropillar-wick model (V2019 framework, no electrical terms).

Modules:
    model        statuses, geometry, properties, settings
    analytic     closed-form verification limits (PROJECT_DERIVED)
    transport    cell closures (constant, tabulated) and Darcy/Kirchhoff liquid transport
    thermal      x-z substrate conduction (finite volumes)
    solver       coupled iteration, endpoint-constrained limit solve, dryout search
    diagnostics  conservation checks and regime scales
    tables       strict cell-table schema, validation and qualification
    importer     raw COMSOL export manifest -> strict dataset
    case         named cases, preflight, problem construction
    output       reproducible run directories

See docs/passive_device_model.md and docs/cell_table_contract.md.
"""

from .model import CaseRole, DataRole, DeviceGeometry, LiquidProperties, SolverSettings, SolverStatus

__all__ = ["CaseRole", "DataRole", "DeviceGeometry", "LiquidProperties", "SolverSettings", "SolverStatus"]
