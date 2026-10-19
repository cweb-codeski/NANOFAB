"""Electrical formation factor F = sigma_bulk / sigma_eff of a square pillar array (in-plane field).

F is a geometric property of the liquid-filled pore space, like the intrinsic permeability K: it
does not depend on the field strength or on zeta. With insulating (oxide-covered) pillars that span
the full liquid height and a flat meniscus at the pillar tops, the potential is independent of z
and F is the 2D effective conductivity of a square array of insulating cylinders.

Three routes, in decreasing order of authority once the COMSOL run exists:

1. COMSOL (Electric Currents in the Stage 1 unit cell, same mesh):
       F = sigma_bulk A dV / (I L)          (formation_factor_from_current)
2. Numerical (this module): cell-centered finite volumes on the full unit cell, sub-cell fluid
   fractions, harmonic face conductances, Richardson extrapolation (formation_factor_numerical).
3. Maxwell (Maxwell-Garnett) 2D estimate for insulating cylinders (formation_factor_maxwell_2d):
       F = (1 + phi) / (1 - phi)
   Rayleigh's square-array corrections enter at O(phi^4) (cited from memory, not held).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..geometry import check_pillar_geometry, solid_fraction


def formation_factor_maxwell_2d(phi: float) -> float:
    """Maxwell 2D formation factor of insulating cylinders at areal solid fraction phi [-]."""
    if not 0 <= phi < 1:
        raise ValueError("solid fraction must be in [0, 1)")
    return (1 + phi) / (1 - phi)


def formation_factor_from_current(current: float, conductivity: float, area: float, voltage: float,
                                  length: float) -> float:
    """F from a simulated or measured current I [A] through a cell of cross-section A, length L.

    Same bookkeeping as permeability.permeability_from_flow: A is the full cross-section including
    pillars, and voltage is the potential drop applied over length.
    """
    if current <= 0 or voltage <= 0:
        raise ValueError("current and voltage must be positive")
    return conductivity * area * voltage / (current * length)


def _fluid_fraction(d: float, l: float, n: int, subsample: int) -> np.ndarray:
    s = (np.arange(subsample) + 0.5) / subsample
    centers = ((np.arange(n)[:, None] + s[None, :]).ravel()) * (l / n)
    X, Y = np.meshgrid(centers, centers, indexing="xy")
    solid = (X - l / 2) ** 2 + (Y - l / 2) ** 2 < (d / 2) ** 2
    return 1.0 - solid.reshape(n, subsample, n, subsample).mean(axis=(1, 3))


def _harmonic(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.zeros_like(a)
    m = (a > 0) & (b > 0)
    out[m] = 2 * a[m] * b[m] / (a[m] + b[m])
    return out


def formation_factor_grid(d: float, l: float, n: int, subsample: int = 8) -> float:
    """F on one n x n grid of the unit cell [0, l]^2 with the pillar at the center.

    The planes x = 0 and x = l (midway between pillar columns) are equipotentials and y = 0, l are
    mirror planes, so phi = 1 / 0 on the x faces and no flux on the y faces is exact for the
    periodic array. Conductivity per cell = fluid fraction (bulk = 1).
    """
    check_pillar_geometry(d, l)
    sig = _fluid_fraction(d, l, n, subsample)
    sig = np.where(sig < 1e-12, 0.0, sig)
    Cx = _harmonic(sig[:, :-1], sig[:, 1:])
    Cy = _harmonic(sig[:-1, :], sig[1:, :])
    Cl = 2 * sig[:, 0]
    Cr = 2 * sig[:, -1]

    active = sig > 0
    idx = -np.ones(sig.shape, dtype=int)
    idx[active] = np.arange(active.sum())
    rows, cols, vals = [], [], []
    diag = np.zeros(active.sum())
    rhs = np.zeros(active.sum())

    def couple(C, ia, ib):
        m = (C > 0) & (ia >= 0) & (ib >= 0)
        a, b, c = ia[m], ib[m], C[m]
        rows.extend([a, b])
        cols.extend([b, a])
        vals.extend([-c, -c])
        np.add.at(diag, a, c)
        np.add.at(diag, b, c)

    couple(Cx, idx[:, :-1], idx[:, 1:])
    couple(Cy, idx[:-1, :], idx[1:, :])
    left, right = idx[:, 0], idx[:, -1]
    ml, mr = left >= 0, right >= 0
    np.add.at(diag, left[ml], Cl[ml])
    np.add.at(rhs, left[ml], Cl[ml] * 1.0)
    np.add.at(diag, right[mr], Cr[mr])
    n_act = active.sum()
    A = sp.csr_matrix(
        (np.concatenate(vals + [diag]), (np.concatenate(rows + [np.arange(n_act)]), np.concatenate(cols + [np.arange(n_act)]))),
        shape=(n_act, n_act),
    )
    phi = spla.spsolve(A, rhs)
    current = float(np.sum(Cl[ml] * (1.0 - phi[left[ml]])))
    current_out = float(np.sum(Cr[mr] * phi[right[mr]]))
    if abs(current - current_out) > 1e-8 * current:
        raise RuntimeError(f"current not conserved: in {current}, out {current_out}")
    # unit depth, square cell: sigma_eff / sigma_bulk = I / dV
    return 1.0 / current


@dataclass(frozen=True)
class FormationFactorStudy:
    n: tuple[int, ...]
    F: tuple[float, ...]
    F_extrapolated: float
    observed_order: float
    maxwell_2d: float

    @property
    def rel_diff_vs_maxwell(self) -> float:
        return (self.F_extrapolated - self.maxwell_2d) / self.maxwell_2d


def formation_factor_numerical(d: float, l: float, n_levels=(64, 128, 256), subsample: int = 8) -> FormationFactorStudy:
    """Grid study of F with Richardson extrapolation from the three finest levels (ratio 2)."""
    if len(n_levels) < 3:
        raise ValueError("need at least three grid levels")
    n_levels = tuple(sorted(n_levels))
    F = tuple(formation_factor_grid(d, l, n, subsample) for n in n_levels)
    f1, f2, f3 = F[-3:]  # coarse -> fine
    r = n_levels[-1] / n_levels[-2]
    num, den = f1 - f2, f2 - f3
    if den == 0 or num / den <= 0:
        p, Fx = float("nan"), f3
    else:
        p = float(np.log(num / den) / np.log(r))
        Fx = f3 - den / (r**p - 1)
    return FormationFactorStudy(n=n_levels, F=F, F_extrapolated=float(Fx), observed_order=p,
                                maxwell_2d=formation_factor_maxwell_2d(solid_fraction(d, l)))
