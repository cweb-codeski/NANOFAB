"""Load provenance-tracked parameters from ``parameters/*.yaml`` and convert them to SI.

Each YAML entry must look like::

    pillar_diameter:
      value: 10
      unit: um
      provenance: PRIMARY_SOURCE      # see nanofab.provenance / docs/SOURCE_AUTHORITY.md
      source: vaartstra2019           # id in config/sources.yaml, or null
      source_location: "Sec. 3.2"     # where in the source, or null
      sweep: [5, 50]                  # optional, same unit as value
      notes: "..."                    # optional

Entries with an invalid ``provenance``, a missing source where one is required, or a source id
that is not registered are rejected.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .provenance import MODEL_CALCULATED, check_reference, load_sources

DEFAULT_PARAMETER_DIR = Path(__file__).resolve().parents[2] / "parameters"
REQUIRED_FIELDS = ("value", "unit", "provenance", "source", "source_location")

# unit string -> (multiplier to SI, SI unit string). Temperatures with offsets are handled
# separately in _to_si.
_UNIT_TABLE: dict[str, tuple[float, str]] = {
    "-": (1.0, "-"),
    # length
    "m": (1.0, "m"),
    "mm": (1e-3, "m"),
    "um": (1e-6, "m"),
    "nm": (1e-9, "m"),
    "cm": (1e-2, "m"),
    # angle
    "rad": (1.0, "rad"),
    "deg": (3.141592653589793 / 180.0, "rad"),
    # pressure
    "Pa": (1.0, "Pa"),
    "kPa": (1e3, "Pa"),
    "bar": (1e5, "Pa"),
    # heat flux
    "W/m2": (1.0, "W/m2"),
    "W/cm2": (1e4, "W/m2"),
    # power
    "W": (1.0, "W"),
    # temperature (absolute)
    "K": (1.0, "K"),
    # material properties
    "W/m/K": (1.0, "W/m/K"),
    "kg/m3": (1.0, "kg/m3"),
    "J/kg/K": (1.0, "J/kg/K"),
    "J/kg": (1.0, "J/kg"),
    "N/m": (1.0, "N/m"),
    "Pa*s": (1.0, "Pa*s"),
    "m2": (1.0, "m2"),
    "m/s2": (1.0, "m/s2"),
    "kg/s": (1.0, "kg/s"),
}


class ParameterError(ValueError):
    """Raised when a parameter file entry is malformed."""


@dataclass(frozen=True)
class Parameter:
    """One source-tracked input. ``value`` is in SI (``si_unit``); ``raw_*`` keep the file entry."""

    name: str
    value: Any
    si_unit: str
    provenance: str
    source: str | None
    source_location: str | None
    raw_value: Any
    raw_unit: str
    sweep: tuple[float, float] | None = None
    notes: str = ""
    file: str = ""


def _to_si(value: Any, unit: str, name: str) -> tuple[Any, str]:
    if isinstance(value, str):
        if unit != "-":
            raise ParameterError(f"{name}: string values must have unit '-'")
        return value, "-"
    if unit == "degC":
        return float(value) + 273.15, "K"
    if unit not in _UNIT_TABLE:
        raise ParameterError(
            f"{name}: unknown unit '{unit}'. Add it to _UNIT_TABLE in params.py with its SI factor."
        )
    factor, si_unit = _UNIT_TABLE[unit]
    return float(value) * factor, si_unit


def _parse_entry(name: str, entry: Any, file: str, sources: dict) -> Parameter:
    if not isinstance(entry, dict):
        raise ParameterError(f"{file}:{name}: entry must be a mapping with {REQUIRED_FIELDS}")
    missing = [k for k in REQUIRED_FIELDS if k not in entry]
    if missing:
        raise ParameterError(f"{file}:{name}: missing required field(s) {missing}")

    provenance = str(entry["provenance"])
    source = entry["source"]
    try:
        check_reference(provenance, source, sources, f"{file}:{name}")
    except ValueError as exc:
        raise ParameterError(str(exc)) from exc

    unit = str(entry["unit"])
    value, si_unit = _to_si(entry["value"], unit, name)

    sweep = None
    if "sweep" in entry:
        lo, hi = entry["sweep"]
        sweep = (_to_si(lo, unit, name)[0], _to_si(hi, unit, name)[0])

    return Parameter(
        name=name,
        value=value,
        si_unit=si_unit,
        provenance=provenance,
        source=source,
        source_location=entry["source_location"],
        raw_value=entry["value"],
        raw_unit=unit,
        sweep=sweep,
        notes=str(entry.get("notes", "")).strip(),
        file=file,
    )


@dataclass
class ParameterSet:
    """Collection of parameters. ``ps["name"]`` returns the SI value; ``ps.meta("name")`` the record."""

    parameters: dict[str, Parameter] = field(default_factory=dict)

    def __getitem__(self, name: str) -> Any:
        return self.parameters[name].value

    def __contains__(self, name: str) -> bool:
        return name in self.parameters

    def meta(self, name: str) -> Parameter:
        return self.parameters[name]

    def with_overrides(self, **si_values: float) -> "ParameterSet":
        """Copy with some SI values replaced (provenance becomes MODEL_CALCULATED, no source)."""
        new = dict(self.parameters)
        for name, value in si_values.items():
            old = self.parameters[name]
            new[name] = Parameter(
                name=name,
                value=value,
                si_unit=old.si_unit,
                provenance=MODEL_CALCULATED,
                source=None,
                source_location="override in script",
                raw_value=value,
                raw_unit=old.si_unit,
                sweep=old.sweep,
                notes=old.notes,
                file=old.file,
            )
        return ParameterSet(new)

    def by_provenance(self, provenance: str) -> list[Parameter]:
        return [p for p in self.parameters.values() if p.provenance == provenance]

    def summary_table(self) -> str:
        """Plain-text traceability table (name, value, unit, provenance, source, location)."""
        rows = [("Parameter", "Value", "Unit", "Provenance", "Source", "Location")]
        for p in self.parameters.values():
            rows.append(
                (p.name, str(p.raw_value), p.raw_unit, p.provenance, str(p.source), str(p.source_location))
            )
        widths = [max(len(r[i]) for r in rows) for i in range(5)]
        lines = []
        for r in rows:
            lines.append("  ".join(r[i].ljust(widths[i]) for i in range(5)) + "  " + r[5])
        return "\n".join(lines)


def load_parameters(directory: str | Path | None = None, sources_file: str | Path | None = None) -> ParameterSet:
    """Load every ``*.yaml`` file in ``directory`` (default: repo ``parameters/``)."""
    directory = Path(directory) if directory is not None else DEFAULT_PARAMETER_DIR
    files = sorted(directory.glob("*.yaml"))
    if not files:
        raise ParameterError(f"No .yaml parameter files found in {directory}")
    sources = load_sources(sources_file)

    params: dict[str, Parameter] = {}
    for path in files:
        with open(path, encoding="utf-8") as fh:
            content = yaml.safe_load(fh) or {}
        for name, entry in content.items():
            if name in params:
                raise ParameterError(f"{name} defined in both {params[name].file} and {path.name}")
            params[name] = _parse_entry(name, entry, path.name, sources)
    return ParameterSet(params)
