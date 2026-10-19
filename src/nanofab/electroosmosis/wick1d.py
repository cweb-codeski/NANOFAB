"""Case A (passive) vs Case B (identical wick + electroosmosis): 1D superficial-velocity model.

Momentum (Darcy + thin-EDL electroosmosis; superposition of two linear Stokes problems):

    U = (K / mu) (-dp_l/dx) + U_EO,      U_EO = mu_EO E / F

    U      total superficial velocity [m/s] (volumetric flow per unit TOTAL cross-section W h)
    K      intrinsic permeability [m^2]   geometric; identical in both cases, never field-dependent
    F      electrical formation factor [-] geometric; identical in both cases
    mu_EO  HS mobility -eps zeta / mu [m^2/(V s)];  E macroscopic field along +x [V/m]

Volumetric flow rate Q = U W h [m^3/s]. Case A is mu_EO E = 0.

Forced-flow configuration (no evaporation, imposed dP over length L):
    Q_pressure = K W h dP / (mu L),   Q_EO = U_EO W h,   Q_total = Q_pressure + Q_EO

Evaporating half-wick (feed at x = 0 with Pc = 0, closed end at x = L, uniform q''):
    Pc = p_v - p_l,  dPc/dx = (mu / K) (U(x) - U_EO),  U(x) = U_0 (1 - x/L),  U_0 = q'' L / (rho h_fg h)
    Pc(x) = (mu / K) [U_0 (x - x^2 / 2L) - U_EO x]
The capillary (dryout) limit is max_x Pc(x) = Pc_max. With Pi = mu U_EO L / (K Pc_max) (EO pressure
over capillary pressure) the dryout heat flux ratio is

    G(Pi) = q''_dry,B / q''_dry,A = (1 + Pi + sqrt(1 + 2 Pi)) / 2    (Pi >= 0)
                                  = 1 + Pi                            (Pi < 0, opposing field)

and the closed end floods (Pc(L) < 0, liquid pressure above vapor pressure) when U_0 < 2 U_EO.
Derivation: docs/electroosmosis_model.md (PROJECT_DERIVED).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..capillary import pillar_max_capillary_pressure
from ..dryout import dryout_heat_flux
from ..geometry import check_pillar_geometry, porosity, solid_fraction
from ..permeability import pillar_wick_permeability
from .formation_factor import formation_factor_maxwell_2d


@dataclass(frozen=True)
class Liquid:
    mu: float  # Pa s
    rho: float  # kg/m^3
    h_fg: float  # J/kg
    surface_tension: float  # N/m
    relative_permittivity: float  # -

    @classmethod
    def from_fluid(cls, fluid, T: float, relative_permittivity: float) -> "Liquid":
        p = fluid.saturated(T)
        return cls(mu=p.mu_l, rho=p.rho_l, h_fg=p.h_lv, surface_tension=p.surface_tension,
                   relative_permittivity=relative_permittivity)


@dataclass(frozen=True)
class Wick:
    """Geometry and geometric transport properties. One instance is shared by Case A and Case B."""

    d: float  # pillar diameter [m]
    l: float  # pitch [m]
    h: float  # pillar (liquid) height [m]
    L: float  # feed length (feed to closed end) [m]
    W: float  # wick width normal to the flow [m]
    K: float  # intrinsic permeability [m^2]
    F: float  # electrical formation factor [-]
    Pc_max: float  # maximum capillary pressure [Pa]

    def __post_init__(self):
        check_pillar_geometry(self.d, self.l, self.h)
        if min(self.L, self.W, self.K, self.Pc_max) <= 0 or self.F < 1:
            raise ValueError("need L, W, K, Pc_max > 0 and F >= 1")

    @property
    def phi(self) -> float:
        return solid_fraction(self.d, self.l)

    @property
    def porosity(self) -> float:
        return porosity(self.d, self.l)

    @property
    def area(self) -> float:
        """Full cross-section normal to the flow, including pillars [m^2] (the K/Darcy area)."""
        return self.W * self.h

    @classmethod
    def analytical(cls, d, l, h, L, W, theta_r, surface_tension, K=None, F=None) -> "Wick":
        """Wick with analytical K (Sangani-Acrivos + Brinkman) and Maxwell F unless given."""
        return cls(
            d=d, l=l, h=h, L=L, W=W,
            K=pillar_wick_permeability(d, l, h) if K is None else K,
            F=formation_factor_maxwell_2d(solid_fraction(d, l)) if F is None else F,
            Pc_max=pillar_max_capillary_pressure(surface_tension, d, l, theta_r),
        )


@dataclass(frozen=True)
class Forcing:
    """Electrical forcing. Case A is Forcing.passive(); Case B adds mu_EO and E only."""

    mobility: float = 0.0  # HS mobility mu_EO [m^2/(V s)]
    E: float = 0.0  # macroscopic field along +x [V/m]

    @classmethod
    def passive(cls) -> "Forcing":
        return cls(0.0, 0.0)

    @property
    def is_passive(self) -> bool:
        return self.mobility * self.E == 0.0


@dataclass(frozen=True)
class Case:
    name: str
    wick: Wick
    liquid: Liquid
    forcing: Forcing


def paired_cases(wick: Wick, liquid: Liquid, forcing: Forcing) -> tuple[Case, Case]:
    """Case A (passive) and Case B (assisted) sharing the SAME wick and liquid objects."""
    return Case("A_passive", wick, liquid, Forcing.passive()), Case("B_electroosmotic", wick, liquid, forcing)


def assert_identical_except_forcing(a: Case, b: Case) -> None:
    """Raise if the two cases differ in anything other than the electrical forcing."""
    if a.wick != b.wick:
        raise ValueError(f"geometry/hydraulics differ between {a.name} and {b.name}: {a.wick} vs {b.wick}")
    if a.liquid != b.liquid:
        raise ValueError(f"liquid properties differ between {a.name} and {b.name}")


# ---------------------------------------------------------------------------------------------
# velocities and flow rates
# ---------------------------------------------------------------------------------------------

def eo_superficial_velocity(case: Case) -> float:
    """U_EO = mu_EO E / F [m/s] (superficial)."""
    return case.forcing.mobility * case.forcing.E / case.wick.F


def eo_velocities(case: Case) -> dict[str, float]:
    """The three EO velocities that must not be confused [m/s]."""
    u_hs = case.forcing.mobility * case.forcing.E
    U = u_hs / case.wick.F
    return {"hs_slip": u_hs, "superficial": U, "intrinsic_pore_average": U / case.wick.porosity}


def pressure_superficial_velocity(case: Case, pressure_gradient):
    """Darcy superficial velocity (K/mu)(-dp/dx) [m/s]; ``pressure_gradient`` = -dp/dx [Pa/m]."""
    return case.wick.K / case.liquid.mu * np.asarray(pressure_gradient, dtype=float)


def forced_flow(case: Case, pressure_drop, length: float | None = None) -> dict[str, np.ndarray]:
    """Forced-flow test: imposed dP = p(0) - p(length) [Pa], no evaporation.

    Returns volumetric flow rates [m^3/s] and superficial velocities [m/s].
    """
    L = case.wick.L if length is None else length
    G = np.asarray(pressure_drop, dtype=float) / L
    A = case.wick.area
    U_p = pressure_superficial_velocity(case, G)
    U_eo = np.full_like(U_p, eo_superficial_velocity(case))
    return {
        "pressure_gradient": G,
        "U_pressure": U_p,
        "U_eo": U_eo,
        "U_total": U_p + U_eo,
        "Q_pressure": U_p * A,
        "Q_eo": U_eo * A,
        "Q_total": (U_p + U_eo) * A,
    }


def apparent_permeability(case: Case, pressure_gradient):
    """K_app = mu U_total / (-dp/dx) [m^2]. NOT a material property: it diverges as dp/dx -> 0."""
    G = np.asarray(pressure_gradient, dtype=float)
    U = pressure_superficial_velocity(case, G) + eo_superficial_velocity(case)
    with np.errstate(divide="ignore", invalid="ignore"):
        return case.liquid.mu * U / G


def equivalent_eo_pressure(case: Case, length: float | None = None) -> float:
    """Pressure drop that would drive U_EO through the same K over ``length`` [Pa]: mu U_EO L / K."""
    L = case.wick.L if length is None else length
    return case.liquid.mu * eo_superficial_velocity(case) * L / case.wick.K


def eo_pressure_ratio(case: Case) -> float:
    """Pi = mu U_EO L / (K Pc_max): EO pumping pressure over the capillary pumping capacity [-]."""
    return equivalent_eo_pressure(case) / case.wick.Pc_max


# ---------------------------------------------------------------------------------------------
# evaporating half-wick: closed form
# ---------------------------------------------------------------------------------------------

def feed_velocity(case: Case, heat_flux) -> np.ndarray:
    """Superficial velocity entering the wick at x = 0 [m/s]: q'' L / (rho h_fg h)."""
    w, f = case.wick, case.liquid
    return np.asarray(heat_flux, dtype=float) * w.L / (f.rho * f.h_fg * w.h)


def capillary_pressure_profile(case: Case, heat_flux: float, x) -> np.ndarray:
    """Pc(x) = p_v - p_l [Pa] along the evaporating half-wick (feed Pc = 0)."""
    w, f = case.wick, case.liquid
    x = np.asarray(x, dtype=float)
    U0 = float(feed_velocity(case, heat_flux))
    return f.mu / w.K * (U0 * (x - x**2 / (2 * w.L)) - eo_superficial_velocity(case) * x)


def required_capillary_pressure(case: Case, heat_flux: float) -> tuple[float, float]:
    """(max Pc over the wick [Pa], location [m]). Liquid pressure is lowest there."""
    w, f = case.wick, case.liquid
    U0 = float(feed_velocity(case, heat_flux))
    Ueo = eo_superficial_velocity(case)
    if U0 <= 0:
        return 0.0, 0.0
    x_star = min(max(w.L * (1 - Ueo / U0), 0.0), w.L)
    return float(capillary_pressure_profile(case, heat_flux, x_star)), x_star


def closed_end_capillary_pressure(case: Case, heat_flux: float) -> float:
    """Pc(L) [Pa]. Negative = liquid above vapor pressure at the closed end (over-pumping)."""
    return float(capillary_pressure_profile(case, heat_flux, case.wick.L))


def passive_dryout_heat_flux(case: Case) -> float:
    """Case A capillary limit [W/m^2] (dryout.dryout_heat_flux with this wick's K and Pc_max)."""
    w, f = case.wick, case.liquid
    return dryout_heat_flux(w.K, w.h, w.L, w.Pc_max, f.rho, f.mu, f.h_fg)


def dryout_gain(Pi):
    """G(Pi) = q''_dry with EO / q''_dry passive for the uniform-flux half-wick [-]."""
    Pi = np.asarray(Pi, dtype=float)
    pos = (1 + Pi + np.sqrt(np.clip(1 + 2 * Pi, 0, None))) / 2
    return np.where(Pi >= 0, pos, np.clip(1 + Pi, 0, None))


def assisted_dryout_heat_flux(case: Case) -> float:
    """Capillary-limited dryout heat flux [W/m^2] including the EO flux (closed form)."""
    return passive_dryout_heat_flux(case) * float(dryout_gain(eo_pressure_ratio(case)))


def flooding_heat_flux(case: Case, pc_flood: float = 0.0) -> float:
    """Lowest q'' [W/m^2] at which Pc(L) >= pc_flood. Below it the field over-pumps the closed end.

    pc_flood = 0 is a flat meniscus at the closed end (liquid pressure = vapor pressure). A pinned
    meniscus can bulge somewhat further; that limit is not modeled (it needs the pinning angle).
    """
    w, f = case.wick, case.liquid
    U0_min = 2 * eo_superficial_velocity(case) + 2 * w.K * pc_flood / (f.mu * w.L)
    return max(U0_min, 0.0) * f.rho * f.h_fg * w.h / w.L


# ---------------------------------------------------------------------------------------------
# evaporating half-wick: independent finite-volume solution (verification of the closed form)
# ---------------------------------------------------------------------------------------------

@dataclass
class FVSolution:
    x: np.ndarray  # cell centers [m]
    Pc: np.ndarray  # [Pa]
    face_x: np.ndarray
    face_flow_per_width: np.ndarray  # volumetric flow per unit width at faces [m^2/s], + toward closed end
    evaporation_per_width: float  # [m^2/s]
    feed_inflow_per_width: float  # [m^2/s]

    @property
    def mass_residual(self) -> float:
        return (self.feed_inflow_per_width - self.evaporation_per_width) / self.evaporation_per_width


def solve_fv(case: Case, heat_flux: float, n: int = 400) -> FVSolution:
    """Cell-centered FV solve of d/dx[h (K/mu) dPc/dx + h U_EO] ... in conservative flux form.

    Face flow per width: f = h [(K/mu) (Pc_R - Pc_L)/dx + U_EO]; feed face Pc = 0 (half cell);
    closed end f = 0 (blocks pressure-driven AND electroosmotic flow). Sink per cell q'' dx/(rho h_fg).
    """
    w, f = case.wick, case.liquid
    dx = w.L / n
    x = (np.arange(n) + 0.5) * dx
    c = w.h * w.K / f.mu / dx  # face conductance per width
    ueo = w.h * eo_superficial_velocity(case)
    sink = heat_flux * dx / (f.rho * f.h_fg)

    main = np.full(n, 2 * c)
    main[0] = c + 2 * c  # interior face to cell 1 + half-cell feed face
    main[-1] = c  # closed end
    off = np.full(n - 1, -c)
    A = sp.diags([off, main, off], [-1, 0, 1], format="csc")
    # f_in - f_out = sink with f = c (Pc_R - Pc_L) + ueo: the EO terms cancel except at the closed end
    rhs = np.full(n, sink)
    rhs[-1] -= ueo
    Pc = spla.spsolve(A, rhs)
    faces = np.concatenate([[0.0], (np.arange(1, n)) * dx, [w.L]])
    flow = np.empty(n + 1)
    flow[0] = 2 * c * (Pc[0] - 0.0) + ueo
    flow[1:-1] = c * (Pc[1:] - Pc[:-1]) + ueo
    flow[-1] = 0.0
    return FVSolution(x=x, Pc=Pc, face_x=faces, face_flow_per_width=flow,
                      evaporation_per_width=sink * n, feed_inflow_per_width=float(flow[0]))


def dryout_heat_flux_fv(case: Case, n: int = 400, rtol: float = 1e-10) -> float:
    """Dryout q'' [W/m^2] from the FV model by bisection on max Pc = Pc_max (independent of G(Pi))."""
    lo, hi = 0.0, 1e4
    while solve_fv(case, hi, n).Pc.max() < case.wick.Pc_max:
        hi *= 2
        if hi > 1e12:
            raise RuntimeError("no dryout found below 1e12 W/m^2")
    while hi - lo > rtol * hi:
        mid = 0.5 * (lo + hi)
        if solve_fv(case, mid, n).Pc.max() < case.wick.Pc_max:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def reynolds_numbers(case: Case, superficial_velocity) -> dict[str, np.ndarray]:
    """Re_K = rho U sqrt(K)/mu (Darcy-Forchheimer onset ~ 1-10) and Re_gap on pore velocity and l - d."""
    w, f = case.wick, case.liquid
    U = np.abs(np.asarray(superficial_velocity, dtype=float))
    return {
        "Re_K": f.rho * U * math.sqrt(w.K) / f.mu,
        "Re_gap": f.rho * (U / w.porosity) * (w.l - w.d) / f.mu,
    }
