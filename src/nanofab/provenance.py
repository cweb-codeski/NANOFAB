"""Provenance categories and the source registry (``config/sources.yaml``).

See ``docs/SOURCE_AUTHORITY.md`` for what each category means and the authority order.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

PRIMARY_SOURCE = "PRIMARY_SOURCE"
AUTHORITATIVE_DATABASE = "AUTHORITATIVE_DATABASE"
PROJECT_DERIVED = "PROJECT_DERIVED"
MODEL_CALCULATED = "MODEL_CALCULATED"
ASSUMED = "ASSUMED"
DESIGN_TARGET = "DESIGN_TARGET"
PROJECT_CONTEXT = "PROJECT_CONTEXT"
COMSOL_EXTRACTED = "COMSOL_EXTRACTED"
USER_PROVIDED = "USER_PROVIDED"

VALID_PROVENANCE = frozenset(
    {
        PRIMARY_SOURCE,
        AUTHORITATIVE_DATABASE,
        PROJECT_DERIVED,
        MODEL_CALCULATED,
        ASSUMED,
        DESIGN_TARGET,
        PROJECT_CONTEXT,
        COMSOL_EXTRACTED,
        USER_PROVIDED,
    }
)

# Categories that must name a registered source.
REQUIRES_SOURCE = frozenset({PRIMARY_SOURCE, AUTHORITATIVE_DATABASE, PROJECT_DERIVED, PROJECT_CONTEXT})

# Five reporting categories used by the whole-wick model (docs/whole_wick_assumptions.md, D8).
LITERATURE = "literature value"
COMSOL = "extracted from COMSOL"
ANALYTICAL = "analytical value"
ESTIMATE = "estimate"
USER = "user-provided value"

REPORT_CATEGORY = {
    PRIMARY_SOURCE: LITERATURE,
    AUTHORITATIVE_DATABASE: LITERATURE,
    PROJECT_DERIVED: LITERATURE,
    COMSOL_EXTRACTED: COMSOL,
    MODEL_CALCULATED: ANALYTICAL,
    ASSUMED: ESTIMATE,
    DESIGN_TARGET: USER,
    PROJECT_CONTEXT: USER,
    USER_PROVIDED: USER,
}

# Weakest-first order used when a derived quantity inherits provenance from several inputs.
# User design choices (e.g. a chosen pitch) do not downgrade a calculated result; ASSUMED inputs do.
_WEAKNESS_ORDER = (
    ASSUMED,
    MODEL_CALCULATED,
    COMSOL_EXTRACTED,
    USER_PROVIDED,
    DESIGN_TARGET,
    PROJECT_CONTEXT,
    PROJECT_DERIVED,
    AUTHORITATIVE_DATABASE,
    PRIMARY_SOURCE,
)


def report_category(label: str) -> str:
    """Map a provenance class to one of the five whole-wick reporting categories."""
    return REPORT_CATEGORY[validate_provenance(label)]


def weakest(*labels: str) -> str:
    """Return the least authoritative of ``labels`` (an ASSUMED input makes the result ASSUMED)."""
    for label in labels:
        validate_provenance(label)
    return min(labels, key=_WEAKNESS_ORDER.index)

# Categories that may appear as literature replication targets.
LITERATURE_TARGET_PROVENANCE = frozenset({PRIMARY_SOURCE, PROJECT_DERIVED})

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCES_FILE = REPO_ROOT / "config" / "sources.yaml"


class ProvenanceError(ValueError):
    """Raised for an invalid provenance label or a missing/unknown source."""


def validate_provenance(label: str, context: str = "") -> str:
    """Return ``label`` if it is a valid category, else raise ``ProvenanceError``."""
    if label not in VALID_PROVENANCE:
        where = f"{context}: " if context else ""
        raise ProvenanceError(f"{where}provenance '{label}' not in {sorted(VALID_PROVENANCE)}")
    return label


def load_sources(path: str | Path | None = None) -> dict[str, dict[str, Any]]:
    """Load the source registry, validating each entry's ``type``."""
    path = Path(path) if path is not None else SOURCES_FILE
    with open(path, encoding="utf-8") as fh:
        sources = yaml.safe_load(fh) or {}
    for sid, entry in sources.items():
        validate_provenance(entry.get("type", ""), f"sources.yaml:{sid}")
    return sources


def check_reference(
    provenance: str,
    source: str | None,
    sources: dict[str, dict[str, Any]],
    context: str = "",
) -> None:
    """Validate a (provenance, source) pair against the registry."""
    validate_provenance(provenance, context)
    where = f"{context}: " if context else ""
    if provenance in REQUIRES_SOURCE and not source:
        raise ProvenanceError(f"{where}{provenance} requires a source id from config/sources.yaml")
    if source and source not in sources:
        raise ProvenanceError(f"{where}source '{source}' is not registered in config/sources.yaml")
