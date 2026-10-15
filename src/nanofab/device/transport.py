"""Cell-property closures and the one-dimensional liquid-transport law.

A *closure* supplies, as functions of the local suction Pc [Pa]:

* ``tau(Pc)``  bulk transmissivity h K_bulk(Pc) [m^3], with Q/W = tau/mu dPc/dx;
* ``h_e(Pc)``  effective base-plane-to-vapor HTC [W/(m^2 K)] (V2019 Eq. 4 definition);
* the supported pressure interval and the limits that end it.

Every evaluation outside the supported interval raises ``TableRangeError``. There is no
extrapolation and no silent clamping. The only tolerance is a narrow endpoint snap: a query within
``snap`` [Pa] of an end knot is evaluated at that knot (it absorbs floating-point roundoff from
cumulative sums, nothing more).

Transport law. For bulk-normalized Darcy flow (a depth-integrated form of V2019 Eq. 8),
Q/W = (tau(Pc)/mu) dPc/dx. With the Kirchhoff transform Psi(Pc) = integral tau dPc this is
Q/W = (1/mu) dPsi/dx, so for a piecewise-constant evaporation field the face flows (from mass
conservation) and Psi (from exact integration of the linear-in-cell Q) are obtained exactly. A
future gradient-dependent law Q(Pc, dPc/dx) would replace ``DarcyHydraulics.march`` with an ODE
integration using the same face flows; the coupled solver only calls ``march``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


class TableRangeError(ValueError):
    """A closure was queried outside its supported pressure interval."""


class PiecewiseLinear:
    """Positive piecewise-linear interpolant with strict range checks and exact integral/inverse."""

    def __init__(self, x, y, name: str = "", snap: float = 0.0):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.ndim != 1 or x.shape != y.shape or x.size < 2:
            raise ValueError(f"{name}: need matching 1-D arrays with at least two knots")
        if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
            raise ValueError(f"{name}: nonfinite knot")
        if np.any(np.diff(x) <= 0):
            raise ValueError(f"{name}: knots must be strictly increasing (duplicate or unsorted state)")
        if np.any(y <= 0):
            raise ValueError(f"{name}: values must be positive")
        if snap < 0:
            raise ValueError("snap must be nonnegative")
        self.x, self.y, self.name, self.snap = x, y, name, float(snap)
        seg = 0.5 * (y[1:] + y[:-1]) * np.diff(x)
        self._cum = np.concatenate([[0.0], np.cumsum(seg)])

    @property
    def lo(self) -> float:
        return float(self.x[0])

    @property
    def hi(self) -> float:
        return float(self.x[-1])

    def _checked(self, p) -> np.ndarray:
        p = np.atleast_1d(np.asarray(p, dtype=float)).copy()
        if not np.all(np.isfinite(p)):
            raise TableRangeError(f"{self.name}: nonfinite query")
        lo_snap = (p < self.lo) & (p >= self.lo - self.snap)
        hi_snap = (p > self.hi) & (p <= self.hi + self.snap)
        p[lo_snap] = self.lo
        p[hi_snap] = self.hi
        bad = (p < self.lo) | (p > self.hi)
        if np.any(bad):
            raise TableRangeError(
                f"{self.name}: query {p[bad].min() if np.any(p < self.lo) else p[bad].max():.9g} outside "
                f"supported interval [{self.lo:.9g}, {self.hi:.9g}] (no extrapolation)"
            )
        return p

    def __call__(self, p):
        p = self._checked(p)
        return np.interp(p, self.x, self.y)

    def integral(self, p):
        """integral_{lo}^{p} y dP (exact for the piecewise-linear interpolant)."""
        p = self._checked(p)
        i = np.clip(np.searchsorted(self.x, p, side="right") - 1, 0, self.x.size - 2)
        s = p - self.x[i]
        m = (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
        return self._cum[i] + self.y[i] * s + 0.5 * m * s**2

    def integral_inverse(self, F):
        """P such that integral(P) = F (exact; positive y makes the integral strictly increasing)."""
        F = np.atleast_1d(np.asarray(F, dtype=float)).copy()
        total = self._cum[-1]
        tolF = self.snap * float(np.max(self.y))
        F[(F < 0) & (F >= -tolF)] = 0.0
        F[(F > total) & (F <= total + tolF)] = total
        if np.any((F < 0) | (F > total)) or not np.all(np.isfinite(F)):
            raise TableRangeError(f"{self.name}: integral value outside the supported interval")
        i = np.clip(np.searchsorted(self._cum, F, side="right") - 1, 0, self.x.size - 2)
        dF = F - self._cum[i]
        a = self.y[i]
        m = (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
        # 0.5 m s^2 + a s - dF = 0, stable root for either sign of m
        s = 2 * dF / (a + np.sqrt(np.maximum(a * a + 2 * m * dF, 0.0)))
        return np.minimum(self.x[i] + s, self.x[i + 1])


@dataclass
class ClosureLimits:
    """Pressure limits attached to a closure. Values are None where not established."""

    Pc_physical_limit: float | None = None
    physical_limit_criterion: str = ""
    physical_limit_demonstrated: bool = False
    Pc_geometry_limit: float | None = None
    geometry_limit_criterion: str = ""
    geometry_limit_demonstrated: bool = False
    supported_gradient_Pa_m: tuple[float, float] = (0.0, float("inf"))
    supported_heat_flux_W_m2: tuple[float, float] = (0.0, float("inf"))
    notes: list[str] = field(default_factory=list)


class Closure:
    """Base class: tau(Pc), its Kirchhoff transform and inverse, h_e(Pc), and the supported range."""

    name = "closure"
    limits: ClosureLimits
    has_intrinsic_area: bool = False
    snap: float = 0.0

    @property
    def Pc_lo(self) -> float:
        raise NotImplementedError

    @property
    def Pc_hi(self) -> float:
        raise NotImplementedError

    def tau(self, Pc):
        raise NotImplementedError

    def psi(self, Pc):
        """Kirchhoff transform integral_{Pc_lo}^{Pc} tau dP [m^3 Pa]."""
        raise NotImplementedError

    def psi_inverse(self, Psi):
        raise NotImplementedError

    def h_e(self, Pc):
        raise NotImplementedError

    def liquid_section_area_per_width(self, Pc):
        """A_liquid_section / W [m] at the identified section, if the meniscus table provides it."""
        return None

    def describe(self) -> dict:
        return {"name": self.name, "Pc_lo": self.Pc_lo, "Pc_hi": self.Pc_hi}


class ConstantClosure(Closure):
    """Constant K_bulk and h_e over [Pc_lo, Pc_hi]. Used for analytical verification limits."""

    def __init__(self, K_bulk: float, h: float, h_e: float, Pc_hi: float, Pc_lo: float = 0.0,
                 limits: ClosureLimits | None = None, snap: float = 0.0):
        if not (K_bulk > 0 and h > 0 and h_e > 0 and Pc_hi > Pc_lo):
            raise ValueError("ConstantClosure needs positive K, h, h_e and Pc_hi > Pc_lo")
        self.K, self.hh, self.he = float(K_bulk), float(h), float(h_e)
        self._lo, self._hi, self.snap = float(Pc_lo), float(Pc_hi), float(snap)
        self.limits = limits or ClosureLimits(Pc_physical_limit=Pc_hi, physical_limit_criterion="declared",
                                              physical_limit_demonstrated=True)
        self.name = "constant"

    @property
    def Pc_lo(self):
        return self._lo

    @property
    def Pc_hi(self):
        return self._hi

    def _checked(self, Pc):
        p = np.atleast_1d(np.asarray(Pc, float)).copy()
        p[(p < self._lo) & (p >= self._lo - self.snap)] = self._lo
        p[(p > self._hi) & (p <= self._hi + self.snap)] = self._hi
        if not np.all(np.isfinite(p)) or np.any((p < self._lo) | (p > self._hi)):
            raise TableRangeError(f"constant closure: query outside [{self._lo}, {self._hi}]")
        return p

    def tau(self, Pc):
        return np.full_like(self._checked(Pc), self.K * self.hh)

    def psi(self, Pc):
        return self.K * self.hh * (self._checked(Pc) - self._lo)

    def psi_inverse(self, Psi):
        Psi = np.atleast_1d(np.asarray(Psi, float))
        return self._checked(self._lo + Psi / (self.K * self.hh))

    def h_e(self, Pc):
        return np.full_like(self._checked(Pc), self.he)

    def describe(self):
        return {**super().describe(), "K_bulk_m2": self.K, "h_m": self.hh, "h_e_W_m2K": self.he}


class TabulatedClosure(Closure):
    """Piecewise-linear tau(Pc) and h_e(Pc) from imported cell tables.

    The supported interval is the intersection of the hydraulic, thermal and (if present) meniscus
    table intervals. Interpolation is linear in tau and h_e, which preserves positivity.
    """

    def __init__(self, Pc_tau, tau, Pc_he, he, limits: ClosureLimits, snap: float = 0.0,
                 Pc_area=None, area_per_width=None, name: str = "tabulated"):
        self._tau = PiecewiseLinear(Pc_tau, tau, name=f"{name}:tau", snap=snap)
        self._he = PiecewiseLinear(Pc_he, he, name=f"{name}:h_e", snap=snap)
        self._area = None
        if Pc_area is not None and area_per_width is not None:
            self._area = PiecewiseLinear(Pc_area, area_per_width, name=f"{name}:A_liquid", snap=snap)
            self.has_intrinsic_area = True
        self.limits = limits
        self.name = name
        self.snap = snap
        self._lo = max(self._tau.lo, self._he.lo)
        self._hi = min(self._tau.hi, self._he.hi)
        if self._area is not None:
            self._lo = max(self._lo, self._area.lo)
            self._hi = min(self._hi, self._area.hi)
        if not self._hi > self._lo:
            raise ValueError(f"{name}: hydraulic and thermal tables have no overlapping pressure range")
        self._psi_lo = float(self._tau.integral(self._lo)[0])

    @property
    def Pc_lo(self):
        return self._lo

    @property
    def Pc_hi(self):
        return self._hi

    def _in_range(self, Pc):
        p = np.atleast_1d(np.asarray(Pc, float)).copy()
        p[(p < self._lo) & (p >= self._lo - self.snap)] = self._lo
        p[(p > self._hi) & (p <= self._hi + self.snap)] = self._hi
        if not np.all(np.isfinite(p)) or np.any((p < self._lo) | (p > self._hi)):
            raise TableRangeError(
                f"{self.name}: query outside the intersected supported interval [{self._lo:.9g}, {self._hi:.9g}]")
        return p

    def tau(self, Pc):
        return self._tau(self._in_range(Pc))

    def psi(self, Pc):
        return self._tau.integral(self._in_range(Pc)) - self._psi_lo

    def psi_inverse(self, Psi):
        Psi = np.atleast_1d(np.asarray(Psi, float))
        return self._in_range(self._tau.integral_inverse(Psi + self._psi_lo))

    def h_e(self, Pc):
        return self._he(self._in_range(Pc))

    def liquid_section_area_per_width(self, Pc):
        if self._area is None:
            return None
        return self._area(self._in_range(Pc))

    def describe(self):
        return {**super().describe(), "tau_knots": int(self._tau.x.size), "h_e_knots": int(self._he.x.size),
                "intrinsic_area_available": self._area is not None}


@dataclass
class HydraulicState:
    """Result of one hydraulic march for a given piecewise-constant evaporation field."""

    Q_faces: np.ndarray  # m^3/s at x_faces
    Pc_faces: np.ndarray | None
    Pc_centers: np.ndarray | None
    psi_faces: np.ndarray  # m^3 Pa, relative to Pc_lo
    psi_centers: np.ndarray
    exceeds: bool
    psi_required_end: float
    psi_cap: float


class DarcyHydraulics:
    """Bulk-normalized Darcy transport with exact Kirchhoff integration (see module docstring)."""

    def __init__(self, closure: Closure, mu: float, W: float, rho: float, h_fg: float):
        self.closure, self.mu, self.W, self.rho, self.h_fg = closure, float(mu), float(W), float(rho), float(h_fg)

    def face_flows(self, x_faces, q_cells) -> np.ndarray:
        """Q at faces from mass conservation with Q(L) = 0 [m^3/s]."""
        dx = np.diff(x_faces)
        sink = self.W * dx * np.asarray(q_cells, float) / (self.rho * self.h_fg)
        Q = np.zeros(len(x_faces))
        Q[:-1] = np.cumsum(sink[::-1])[::-1]
        return Q

    def march(self, x_faces, q_cells, Pc_feed: float, Pc_cap: float) -> HydraulicState:
        """Pc at faces and centers for evaporation q_cells; no evaluation above Pc_cap."""
        x_faces = np.asarray(x_faces, float)
        q_cells = np.asarray(q_cells, float)
        dx = np.diff(x_faces)
        Q = self.face_flows(x_faces, q_cells)
        sink = Q[:-1] - Q[1:]
        # integral of Q/W over a cell and over its left half (Q linear inside the cell)
        full = dx * (Q[:-1] + Q[1:]) / (2 * self.W)
        half = dx / (2 * self.W) * (Q[:-1] - sink / 4)
        psi0 = float(self.closure.psi(Pc_feed)[0])
        psi_faces = psi0 + self.mu * np.concatenate([[0.0], np.cumsum(full)])
        psi_centers = psi_faces[:-1] + self.mu * half
        psi_cap = float(self.closure.psi(Pc_cap)[0])
        tol = self.closure.snap * float(np.max(self.closure.tau(Pc_cap)))
        exceeds = bool(np.max(psi_faces) > psi_cap + tol)
        if exceeds:
            return HydraulicState(Q, None, None, psi_faces, psi_centers, True, float(psi_faces[-1]), psi_cap)
        psi_faces_c = np.minimum(psi_faces, psi_cap)
        psi_centers_c = np.minimum(psi_centers, psi_cap)
        return HydraulicState(Q, self.closure.psi_inverse(psi_faces_c), self.closure.psi_inverse(psi_centers_c),
                              psi_faces, psi_centers, False, float(psi_faces[-1]), psi_cap)

    def gradient(self, Q_faces, Pc_faces) -> np.ndarray:
        """dPc/dx at faces [Pa/m] from the local transport law."""
        return self.mu * np.asarray(Q_faces) / (self.W * self.closure.tau(Pc_faces))
