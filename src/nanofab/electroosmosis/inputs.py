"""Load parameters/electroosmosis/eo_inputs.yaml with provenance checks.

Every numeric leaf goes through the whole-wick loader (``whole_wick.config.resolve_tree``), so it
carries a provenance label and a registered source. The ``requires_source`` section is parsed
separately into ``MissingInput`` records: those entries must have ``value: null`` and are never
converted to numbers. ``require()`` raises ``MissingInputError`` listing them, which is how
device-specific EO predictions are blocked until the inputs exist.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from pathlib import Path

from ..capillary import pillar_max_capillary_pressure
from ..fluids import Fluid
from ..geometry import solid_fraction
from ..params import Parameter, ParameterError
from ..provenance import REPO_ROOT
from ..whole_wick.config import Config, load_raw, resolve_tree
from .formation_factor import formation_factor_maxwell_2d
from .physics import electroosmotic_mobility
from .wick1d import Liquid, Wick

DEFAULT_FILE = REPO_ROOT / "parameters" / "electroosmosis" / "eo_inputs.yaml"


class MissingInputError(RuntimeError):
    """A device-specific input marked REQUIRES_SOURCE was needed."""


@dataclass(frozen=True)
class MissingInput:
    name: str
    unit: str
    needed_for: str
    notes: str = ""
    status: str = "REQUIRES_SOURCE"


@dataclass
class EOInputs:
    config: Config
    missing: dict[str, MissingInput]
    file: str

    def leaf(self, dotted: str) -> Parameter:
        node = self.config.get(dotted)
        if not isinstance(node, Parameter):
            raise KeyError(f"{dotted} is not a numeric leaf")
        return node

    def si(self, dotted: str):
        return self.leaf(dotted).value

    def require(self, *names: str) -> None:
        """Raise unless none of ``names`` is still REQUIRES_SOURCE."""
        blocked = [self.missing[n] for n in names if n in self.missing]
        if blocked:
            lines = [f"  {m.name} [{m.unit}]: needed for {m.needed_for}" for m in blocked]
            raise MissingInputError("REQUIRES_SOURCE inputs are not available:\n" + "\n".join(lines))

    def provenance_rows(self) -> list[dict]:
        rows = []
        for path, p in sorted(self.config.registry.items()):
            rows.append({"parameter": path, "value": p.raw_value, "unit": p.raw_unit, "provenance": p.provenance,
                         "source": p.source, "source_location": p.source_location})
        for m in self.missing.values():
            rows.append({"parameter": f"requires_source.{m.name}", "value": None, "unit": m.unit,
                         "provenance": m.status, "source": None, "source_location": m.needed_for})
        return rows


def _parse_missing(section: dict | None, file: str) -> dict[str, MissingInput]:
    out: dict[str, MissingInput] = {}
    for name, entry in (section or {}).items():
        if not isinstance(entry, dict) or "unit" not in entry or "needed_for" not in entry:
            raise ParameterError(f"{file}:requires_source.{name}: needs 'unit' and 'needed_for'")
        if entry.get("value") is not None:
            raise ParameterError(
                f"{file}:requires_source.{name}: has a value; move it to a scenario with full provenance instead"
            )
        out[name] = MissingInput(name=name, unit=str(entry["unit"]), needed_for=str(entry["needed_for"]),
                                 notes=str(entry.get("notes", "")).strip())
    return out


def load_eo_inputs(path: str | Path | None = None) -> EOInputs:
    path = Path(path) if path is not None else DEFAULT_FILE
    raw = copy.deepcopy(load_raw(path))
    missing = _parse_missing(raw.pop("requires_source", None), path.name)
    return EOInputs(config=resolve_tree(raw, path.name), missing=missing, file=path.name)


def baseline_liquid(inp: EOInputs, scenario: str = "pujahari2022") -> Liquid:
    """IAPWS (CoolProp) liquid at T_v with the scenario's relative permittivity."""
    fluid = Fluid(inp.config.get("operating.fluid"))
    T = inp.si("operating.vapor_temperature")
    return Liquid.from_fluid(fluid, T, inp.si(f"literature_scenarios.{scenario}.relative_permittivity"))


def baseline_wick(inp: EOInputs, liquid: Liquid) -> Wick:
    """The single wick shared by Case A and Case B (COMSOL K, geometric F, V2019 Eq. 3 Pc_max)."""
    d, l, h = (inp.si(f"wick.{k}") for k in ("pillar_diameter", "pillar_pitch", "pillar_height"))
    F_spec = inp.config.get("wick.formation_factor")
    if F_spec == "maxwell_2d":
        F = formation_factor_maxwell_2d(solid_fraction(d, l))
    elif isinstance(F_spec, Parameter):
        F = F_spec.value
    else:
        raise ParameterError(f"wick.formation_factor: expected 'maxwell_2d' or a leaf, got {F_spec!r}")
    return Wick(
        d=d, l=l, h=h,
        L=inp.si("wick.feed_length"),
        W=inp.si("wick.width"),
        K=inp.si("wick.permeability"),
        F=F,
        Pc_max=pillar_max_capillary_pressure(liquid.surface_tension, d, l, inp.si("wick.receding_contact_angle")),
    )


def scenario_mobilities(inp: EOInputs, liquid: Liquid, scenario: str = "pujahari2022") -> list[dict]:
    """mu_EO for each (concentration, zeta, Debye length) triple of a literature scenario."""
    base = f"literature_scenarios.{scenario}"
    zetas = list(inp.si(f"{base}.zeta_potential"))
    concs = list(inp.si(f"{base}.concentration"))
    debyes = list(inp.si(f"{base}.debye_length"))
    if not len(zetas) == len(concs) == len(debyes):
        raise ParameterError(f"{base}: zeta, concentration and debye_length lists must pair up")
    eps_r = inp.si(f"{base}.relative_permittivity")
    return [
        {"zeta": z, "concentration": c, "debye_length": lam,
         "mobility": electroosmotic_mobility(eps_r, z, liquid.mu),
         "label": f"zeta = {z*1e3:.0f} mV, {c/1e3:.0e} M KCl"}
        for z, c, lam in zip(zetas, concs, debyes)
    ]


def theta(inp: EOInputs, dotted: str) -> float:
    v = inp.si(dotted)
    if not 0 <= v < math.pi / 2:
        raise ParameterError(f"{dotted}: contact angle out of range")
    return v
