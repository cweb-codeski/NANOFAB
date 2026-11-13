"""Reference and proposed-design geometry cases (``parameters/unit_cell/cases.yaml``).

The two cases are kept apart on purpose: COMSOL results (K, h_e, F, meniscus shapes) belong to one
geometry and are never transferred to the other. Device-level quantities that are not yet defined
(feed length, substrate thickness, heater footprint, reservoirs, vapor environment, CAD scale) are
stored as UNRESOLVED and :meth:`GeometryCase.si` refuses to return a number for them.

Array extent of a square array with ``n`` pillars per row at pitch ``l`` and diameter ``d``:
    center-to-center span   (n - 1) l
    outer pillar-edge span  (n - 1) l + d
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .params import DEFAULT_PARAMETER_DIR, Parameter, ParameterError, load_parameters, parse_entry
from .provenance import load_sources

CASES_FILE = DEFAULT_PARAMETER_DIR / "unit_cell" / "cases.yaml"
ROLES = ("REFERENCE", "PROPOSED_DESIGN")
SECTIONS = ("unit_cell", "device")


class UnresolvedInputError(KeyError):
    """Raised when a calculation asks for an input that the case file marks UNRESOLVED."""


@dataclass(frozen=True)
class GeometryCase:
    name: str
    role: str
    description: str
    parameters: dict[str, Parameter] = field(default_factory=dict)
    unresolved: dict[str, str] = field(default_factory=dict)

    def si(self, key: str) -> float:
        """SI value of ``key``. Raises UnresolvedInputError for UNRESOLVED inputs."""
        if key in self.unresolved:
            raise UnresolvedInputError(f"{self.name}: '{key}' is UNRESOLVED: {self.unresolved[key]}")
        return self.parameters[key].value

    def meta(self, key: str) -> Parameter:
        return self.parameters[key]

    @property
    def d(self) -> float:
        return self.si("pillar_diameter")

    @property
    def l(self) -> float:  # noqa: E743 - pitch symbol used throughout the repo
        return self.si("pillar_pitch")

    @property
    def hp(self) -> float:
        return self.si("pillar_height")


def _resolve_leaf(name: str, entry, base, file: str, sources: dict) -> tuple[Parameter | None, str | None]:
    if isinstance(entry, dict) and set(entry) == {"ref"}:
        if entry["ref"] not in base:
            raise ParameterError(f"{file}:{name}: ref '{entry['ref']}' not found in parameters/*.yaml")
        return base.meta(entry["ref"]), None
    if isinstance(entry, dict) and entry.get("status") == "UNRESOLVED":
        if entry.get("value") is not None:
            raise ParameterError(f"{file}:{name}: an UNRESOLVED entry must have value: null")
        return None, str(entry.get("notes", "")).strip()
    return parse_entry(name, entry, file, sources), None


def load_geometry_cases(path: str | Path | None = None) -> dict[str, GeometryCase]:
    """All cases in the case file, keyed by case name."""
    path = Path(path) if path is not None else CASES_FILE
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    base = load_parameters()
    sources = load_sources()
    out: dict[str, GeometryCase] = {}
    for cname, spec in raw["cases"].items():
        role = spec.get("role")
        if role not in ROLES:
            raise ParameterError(f"{path.name}:{cname}: role must be one of {ROLES}, got {role!r}")
        params: dict[str, Parameter] = {}
        unresolved: dict[str, str] = {}
        for section in SECTIONS:
            for key, entry in (spec.get(section) or {}).items():
                if key in params or key in unresolved:
                    raise ParameterError(f"{path.name}:{cname}: '{key}' defined twice")
                param, reason = _resolve_leaf(f"{cname}.{key}", entry, base, path.name, sources)
                if param is None:
                    unresolved[key] = reason or "UNRESOLVED"
                else:
                    params[key] = param
        out[cname] = GeometryCase(cname, role, str(spec.get("description", "")).strip(), params, unresolved)
    return out


def load_screening_inputs(path: str | Path | None = None) -> dict[str, Parameter]:
    """Shared screening inputs (e.g. room-temperature permittivity) from the same file."""
    path = Path(path) if path is not None else CASES_FILE
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    sources = load_sources()
    return {k: parse_entry(k, v, path.name, sources) for k, v in (raw.get("screening_checks") or {}).items()}


@dataclass(frozen=True)
class ArrayExtent:
    """Spans [m] of an n_cols x n_rows square array."""

    center_span_x: float
    center_span_y: float
    edge_span_x: float
    edge_span_y: float


def array_extent(n_cols: int, n_rows: int, pitch: float, diameter: float) -> ArrayExtent:
    if n_cols < 1 or n_rows < 1:
        raise ValueError("array needs at least one column and one row")
    cx, cy = (n_cols - 1) * pitch, (n_rows - 1) * pitch
    return ArrayExtent(cx, cy, cx + diameter, cy + diameter)


def start_reference_interpretations(case: GeometryCase) -> list[dict]:
    """Outer pillar-edge bounding box [m] of the array under each reading of the start reference.

    Readings (none is adopted; the CAD files decide):
      first_pillar_center  the start point is the center of the first pillar
      array_edge_corner    the start point is the lower-left corner of the outer pillar edges
    Each row also reports whether the box lies inside the coupon and how far its center is from
    the coupon center.
    """
    ext = array_extent(int(case.si("array_columns")), int(case.si("array_rows")), case.l, case.d)
    x0, y0 = case.si("array_start_reference_x"), case.si("array_start_reference_y")
    cxmin, cxmax = case.si("coupon_x_min"), case.si("coupon_x_max")
    cymin, cymax = case.si("coupon_y_min"), case.si("coupon_y_max")
    boxes = {
        "first_pillar_center": (x0 - case.d / 2, y0 - case.d / 2),
        "array_edge_corner": (x0, y0),
    }
    rows = []
    for reading, (xmin, ymin) in boxes.items():
        xmax, ymax = xmin + ext.edge_span_x, ymin + ext.edge_span_y
        rows.append(
            {
                "reading": reading,
                "x_min": xmin, "x_max": xmax, "y_min": ymin, "y_max": ymax,
                "inside_coupon": cxmin <= xmin and xmax <= cxmax and cymin <= ymin and ymax <= cymax,
                "offset_from_coupon_center_x": (xmin + xmax) / 2 - (cxmin + cxmax) / 2,
                "offset_from_coupon_center_y": (ymin + ymax) / 2 - (cymin + cymax) / 2,
            }
        )
    return rows
