"""Lookup tables for unit-cell properties with explicit out-of-range handling.

Tables are multilinear on a regular grid of named axes (``Pc`` [Pa], ``E`` [V/m]). Linear
interpolation never overshoots the tabulated values. Queries outside an axis range are **never
extrapolated**. Depending on ``policy``:

* ``"raise"``: raise ``ExtrapolationError``.
* ``"warn"``: clip to the nearest bound, emit a ``warnings.warn`` and record the event in an
  ``ExtrapolationLog`` so it appears in reports.

An axis with a single point (for example ``E = [0]`` for a table with no field data) accepts only
that exact value. Asking for E > 0 from a zero-field table is therefore reported, not silently
answered with zero-field properties.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from ..provenance import validate_provenance

POLICIES = ("raise", "warn")


class ExtrapolationError(ValueError):
    """Raised when a lookup would need values outside the table and the policy is 'raise'."""


class ExtrapolationWarning(UserWarning):
    """Emitted when a lookup is clipped to the table bounds."""


@dataclass
class ExtrapolationLog:
    events: list[dict] = field(default_factory=list)

    def record(self, **event) -> None:
        self.events.append(event)

    def __len__(self) -> int:
        return len(self.events)

    def summary(self) -> list[dict]:
        """One row per (table, axis), merging repeated events from nonlinear iterations."""
        merged: dict[tuple, dict] = {}
        for e in self.events:
            key = (e["table"], e["axis"])
            row = merged.setdefault(key, dict(e))
            row["requested_min"] = min(row["requested_min"], e["requested_min"])
            row["requested_max"] = max(row["requested_max"], e["requested_max"])
            row["n_points"] = max(row["n_points"], e["n_points"])
        return list(merged.values())


@dataclass
class LookupTable:
    """``values`` has one dimension per entry of ``axes`` (in order). All SI."""

    name: str
    unit: str
    axes: dict[str, np.ndarray]
    values: np.ndarray
    provenance: str
    source: str | None = None
    source_location: str | None = None
    notes: str = ""

    def __post_init__(self):
        validate_provenance(self.provenance, f"table {self.name}")
        self.axes = {k: np.asarray(v, dtype=float) for k, v in self.axes.items()}
        self.values = np.asarray(self.values, dtype=float)
        expected = tuple(len(a) for a in self.axes.values())
        if self.values.shape != expected:
            raise ValueError(f"{self.name}: values shape {self.values.shape} != axes {expected}")
        for k, a in self.axes.items():
            if len(a) > 1 and np.any(np.diff(a) <= 0):
                raise ValueError(f"{self.name}: axis {k} must be strictly increasing")
        if not np.all(np.isfinite(self.values)):
            raise ValueError(f"{self.name}: non-finite table values")
        live = [k for k, a in self.axes.items() if len(a) > 1]
        live_values = self.values.reshape([len(self.axes[k]) for k in live]) if live else self.values.ravel()
        self._live = live
        self._interp = (
            RegularGridInterpolator([self.axes[k] for k in live], live_values, method="linear") if live else None
        )

    def bounds(self, axis: str) -> tuple[float, float]:
        a = self.axes[axis]
        return float(a[0]), float(a[-1])

    def __call__(self, policy: str = "warn", log: ExtrapolationLog | None = None, context: str = "", **coords) -> np.ndarray:
        if policy not in POLICIES:
            raise ValueError(f"policy must be one of {POLICIES}")
        missing = set(self.axes) - set(coords)
        if missing:
            raise ValueError(f"{self.name}: missing coordinates {sorted(missing)}")
        arrays = np.broadcast_arrays(*[np.asarray(coords[k], dtype=float) for k in self.axes])
        shape = arrays[0].shape
        clipped = {}
        for k, arr in zip(self.axes, arrays):
            lo, hi = self.bounds(k)
            tol = 1e-9 * max(abs(lo), abs(hi), 1.0)
            out = (arr < lo - tol) | (arr > hi + tol)
            if out.any():
                msg = (
                    f"{self.name}{' [' + context + ']' if context else ''}: axis {k} requested "
                    f"[{arr[out].min():.4g}, {arr[out].max():.4g}] outside table range [{lo:.4g}, {hi:.4g}] "
                    f"at {int(out.sum())} point(s)"
                )
                if policy == "raise":
                    raise ExtrapolationError(msg)
                warnings.warn(msg + "; clipped to bounds (no extrapolation).", ExtrapolationWarning, stacklevel=2)
                if log is not None:
                    log.record(
                        table=self.name,
                        context=context,
                        axis=k,
                        requested_min=float(arr[out].min()),
                        requested_max=float(arr[out].max()),
                        table_min=lo,
                        table_max=hi,
                        n_points=int(out.sum()),
                    )
            clipped[k] = np.clip(arr, lo, hi)
        if self._interp is None:
            return np.full(shape, float(self.values.ravel()[0]))
        pts = np.column_stack([clipped[k].ravel() for k in self._live])
        return self._interp(pts).reshape(shape)
