"""Provenance-checked loading of nested whole-wick YAML files.

A *leaf* is a mapping with a ``value`` key (scalar or list) or a ``ref`` key::

    nx: {value: 40, unit: "-", provenance: USER_PROVIDED, source: null, source_location: null}
    pitch: {ref: pillar_pitch}          # copies an entry from parameters/*.yaml, provenance included
    V: {value: [0, 5, 10], unit: V}   # list values are converted element-wise

``provenance``, ``source`` and ``source_location`` may be inherited from a ``defaults`` mapping
on any enclosing group. A bare number anywhere in the file is rejected, so every numeric input
carries a provenance label. Strings and booleans are structural (map types, names, switches).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from ..params import Parameter, ParameterError, ParameterSet, _to_si, load_parameters
from ..provenance import check_reference, load_sources

PROVENANCE_KEYS = ("provenance", "source", "source_location")


@dataclass
class Config:
    """Resolved config tree. Leaves are ``Parameter`` objects; ``registry`` maps dotted path -> leaf."""

    tree: dict[str, Any]
    registry: dict[str, Parameter] = field(default_factory=dict)
    file: str = ""

    def get(self, dotted: str) -> Any:
        node: Any = self.tree
        for key in dotted.split("."):
            node = node[key]
        return node


def si(node: Any) -> Any:
    """SI value of a leaf (float or ndarray); passes non-leaves through."""
    if isinstance(node, Parameter):
        return node.value
    return node


def _is_leaf(node: Any) -> bool:
    if not isinstance(node, dict):
        return False
    return "ref" in node or ("value" in node and not isinstance(node["value"], dict))


def _parse_leaf(
    path: str, node: dict, inherited: dict, file: str, sources: dict, base: ParameterSet | None
) -> Parameter:
    if "ref" in node:
        if base is None:
            raise ParameterError(f"{file}:{path}: 'ref' used but no base parameter set loaded")
        name = node["ref"]
        if name not in base:
            raise ParameterError(f"{file}:{path}: ref '{name}' not found in parameters/*.yaml")
        return base.meta(name)

    merged = {k: inherited.get(k) for k in PROVENANCE_KEYS}
    merged.update({k: node[k] for k in PROVENANCE_KEYS if k in node})
    if merged["provenance"] is None:
        raise ParameterError(f"{file}:{path}: missing provenance (set it on the leaf or in 'defaults')")
    if "unit" not in node:
        raise ParameterError(f"{file}:{path}: missing unit")
    try:
        check_reference(merged["provenance"], merged["source"], sources, f"{file}:{path}")
    except ValueError as exc:
        raise ParameterError(str(exc)) from exc

    unit = str(node["unit"])
    raw = node["value"]
    if isinstance(raw, list):
        converted = [_to_si(v, unit, path) for v in raw]
        value: Any = np.array([c[0] for c in converted], dtype=float)
        si_unit = converted[0][1] if converted else unit
    else:
        value, si_unit = _to_si(raw, unit, path)
    return Parameter(
        name=path,
        value=value,
        si_unit=si_unit,
        provenance=str(merged["provenance"]),
        source=merged["source"],
        source_location=merged["source_location"],
        raw_value=raw,
        raw_unit=unit,
        notes=str(node.get("notes", "")).strip(),
        file=file,
    )


def _resolve(node: Any, path: str, inherited: dict, file: str, sources: dict, base, registry: dict) -> Any:
    if _is_leaf(node):
        leaf = _parse_leaf(path, node, inherited, file, sources, base)
        registry[path] = leaf
        return leaf
    if isinstance(node, dict):
        local = dict(inherited)
        if isinstance(node.get("defaults"), dict):
            local.update({k: v for k, v in node["defaults"].items() if k in PROVENANCE_KEYS})
        return {
            k: _resolve(v, f"{path}.{k}" if path else str(k), local, file, sources, base, registry)
            for k, v in node.items()
            if k != "defaults"
        }
    if isinstance(node, list):
        return [_resolve(v, f"{path}[{i}]", inherited, file, sources, base, registry) for i, v in enumerate(node)]
    if isinstance(node, bool) or node is None or isinstance(node, str):
        return node
    if isinstance(node, (int, float)):
        raise ParameterError(
            f"{file}:{path}: bare number {node!r} has no provenance. "
            "Write {value: ..., unit: ..., provenance: ...} or inherit provenance from 'defaults'."
        )
    raise ParameterError(f"{file}:{path}: unsupported YAML node {type(node).__name__}")


def resolve_tree(raw: dict, file: str = "<memory>", base: ParameterSet | None = None) -> Config:
    """Resolve an already-parsed YAML mapping (used by sweeps that edit the raw dict)."""
    sources = load_sources()
    if base is None:
        base = load_parameters()
    registry: dict[str, Parameter] = {}
    tree = _resolve(raw, "", {}, file, sources, base, registry)
    return Config(tree=tree, registry=registry, file=file)


def load_raw(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def load_config(path: str | Path, base: ParameterSet | None = None) -> Config:
    path = Path(path)
    return resolve_tree(load_raw(path), path.name, base)
