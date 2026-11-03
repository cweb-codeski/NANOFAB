"""Uniform cell-centered 2D grid. Arrays are indexed ``[j, i]`` = ``[y, x]`` with shape (ny, nx)."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from functools import cached_property

import numpy as np

EDGES = ("left", "right", "bottom", "top")


@dataclass(frozen=True)
class Grid:
    nx: int
    ny: int
    lx: float  # m
    ly: float  # m

    def __post_init__(self):
        if self.nx < 1 or self.ny < 1 or self.lx <= 0 or self.ly <= 0:
            raise ValueError(f"Invalid grid {self}")

    @property
    def dx(self) -> float:
        return self.lx / self.nx

    @property
    def dy(self) -> float:
        return self.ly / self.ny

    @property
    def cell_area(self) -> float:
        return self.dx * self.dy

    @property
    def footprint(self) -> float:
        return self.lx * self.ly

    @property
    def shape(self) -> tuple[int, int]:
        return (self.ny, self.nx)

    @cached_property
    def x(self) -> np.ndarray:
        return (np.arange(self.nx) + 0.5) * self.dx

    @cached_property
    def y(self) -> np.ndarray:
        return (np.arange(self.ny) + 0.5) * self.dy

    @cached_property
    def XY(self) -> tuple[np.ndarray, np.ndarray]:
        return np.meshgrid(self.x, self.y)

    def distance_to_edges(self, edges: list[str]) -> np.ndarray:
        """Distance [m] from each cell center to the nearest of the named edges."""
        X, Y = self.XY
        dist = {"left": X, "right": self.lx - X, "bottom": Y, "top": self.ly - Y}
        if not edges:
            raise ValueError("need at least one edge")
        return np.min([dist[e] for e in edges], axis=0)

    def check_homogenization(self, pitch: float, min_pitches: float = 3.0) -> None:
        """Warn when a cell spans fewer than ``min_pitches`` pillar pitches (assumption W2)."""
        cells = min(self.dx, self.dy) / pitch
        if cells < min_pitches:
            warnings.warn(
                f"Grid cell spans only {cells:.1f} pillar pitches (< {min_pitches}). Homogenized "
                "properties are not representative at this resolution (assumption W2).",
                stacklevel=2,
            )
