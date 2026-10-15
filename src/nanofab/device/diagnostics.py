"""Conservation checks and regime diagnostics computed from a coupled solution.

Balances are recomputed from the returned fields, independently of the solver's convergence flag.
Normalizations use magnitudes of the actual balance terms; at zero load (where every term
vanishes) they fall back to documented nonzero characteristic scales:

* mass:   Q_char = max(Q(0), W L q_char / (rho h_fg)) with q_char = 1e4 W/m^2 (1 W/cm^2);
* energy: E_char = max(sum |q_in| dx, h_e_mean x 1 K x L), per unit width.

Regime numbers are calculated estimates, not rigorous thresholds (docs/passive_device_model.md).
Only the closure's documented supported ranges (pressure gradient, local heat flux) are gated.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate

from .solver import CoupledProblem, CoupledSolution
from .thermal import reconstruct_energy_balance

G_STANDARD = 9.80665  # m/s^2 (CGPM 1901, exact)
Q_CHAR_FLUX = 1e4  # W/m^2


def mass_balance(problem: CoupledProblem, sol: CoupledSolution, independent_quadrature: bool = True) -> dict:
    g, f = problem.geometry, problem.liquid
    m = sol.mesh
    dx = m.dx
    Q = sol.Q_faces
    sink = g.W * dx * sol.q_e / (f.rho * f.h_fg)
    local = Q[:-1] - Q[1:] - sink
    Q_char = max(abs(Q[0]), g.W * g.L * Q_CHAR_FLUX / (f.rho * f.h_fg))
    glob = Q[0] - float(np.sum(sink))
    out = {
        "Q_feed_m3_s": float(Q[0]),
        "Q_end_m3_s": float(Q[-1]),
        "integrated_evaporation_m3_s": float(np.sum(sink)),
        "global_abs_m3_s": glob,
        "global_normalized": abs(glob) / Q_char,
        "local_max_abs_m3_s": float(np.max(np.abs(local))),
        "local_max_normalized": float(np.max(np.abs(local)) / Q_char),
        "end_flow_normalized": abs(float(Q[-1])) / Q_char,
        "normalization_Q_char_m3_s": Q_char,
    }
    if independent_quadrature:
        # Darcy/Kirchhoff consistency: integral_{Pc_l}^{Pc_r} tau dP (adaptive quadrature, not the
        # closure's own antiderivative) vs mu * integral Q/W dx over each cell.
        c = problem.closure
        Pf = sol.Pc_faces
        lhs = np.array([integrate.quad(lambda p: float(c.tau(p)[0]), Pf[i], Pf[i + 1], limit=200,
                                       epsabs=0.0, epsrel=1e-12)[0] for i in range(m.Nx)])
        rhs = f.mu * dx * (Q[:-1] + Q[1:]) / (2 * g.W)
        scale = max(float(np.sum(np.abs(rhs))), float(c.tau(Pf[0])[0]) * 1e-6 * max(c.Pc_hi - c.Pc_lo, 1.0))
        out["darcy_cell_max_abs_m3Pa"] = float(np.max(np.abs(lhs - rhs)))
        out["darcy_cell_max_normalized"] = float(np.max(np.abs(lhs - rhs)) / scale)
    return out


def energy_balance(problem: CoupledProblem, sol: CoupledSolution) -> dict:
    g = problem.geometry
    rec = reconstruct_energy_balance(sol.slab)
    dx = sol.mesh.dx
    E_char = max(float(np.sum(np.abs(sol.q_in) * dx)), float(np.mean(sol.h_e)) * 1.0 * g.L)
    local_scale = max(rec["local_scale_W_per_m"], float(np.mean(sol.h_e)) * 1.0 * float(np.mean(dx)))
    return {
        "input_W_half_device": rec["input_W_per_m"] * g.W,
        "evaporation_W_half_device": rec["evaporation_W_per_m"] * g.W,
        "other_losses_W": 0.0,
        "global_abs_W": rec["global_abs_W_per_m"] * g.W,
        "global_normalized": abs(rec["global_abs_W_per_m"]) / E_char,
        "local_max_abs_W": rec["max_local_abs_W_per_m"] * g.W,
        "local_max_normalized": rec["max_local_abs_W_per_m"] / local_scale,
        "robin_face_max_abs_W_m2": rec["robin_max_abs_W_m2"],
        "normalization_E_char_W": E_char * g.W,
        "linear_solve_residual": sol.slab.linear_residual,
    }


def regime_diagnostics(problem: CoupledProblem, sol: CoupledSolution) -> dict:
    """Solution-based dimensionless scales and supported-range checks."""
    g, f, c = problem.geometry, problem.liquid, problem.closure
    m = sol.mesh
    Q = sol.Q_faces
    A_ref = g.W * g.h
    U = Q / A_ref
    U_feed = float(U[0])
    theta_e = sol.slab.theta_e
    dT_max = float(np.max(theta_e))
    out: dict = {"definitions": {}}
    d, l = g.d, g.l
    porosity = None
    if d is not None and l is not None:
        porosity = 1 - math.pi * d**2 / (4 * l**2)
        gap = l - d
        D_h = 4 * (l**2 - math.pi * d**2 / 4) * g.h / (math.pi * d * g.h + (l**2 - math.pi * d**2 / 4))
    out["U_superficial_feed_m_s"] = U_feed
    out["Re_h_superficial"] = f.rho * U_feed * g.h / f.mu
    out["definitions"]["Re_h_superficial"] = "rho U_feed h / mu, U = Q/(W h) (V2019 A.1 form)"
    if porosity is not None:
        u_pore = U_feed / porosity
        out["u_pore_volume_porosity_estimate_m_s"] = u_pore
        out["Re_h_pore_estimate"] = f.rho * u_pore * g.h / f.mu
        out["Re_gap_pore_estimate"] = f.rho * u_pore * gap / f.mu
        out["Re_Dh_pore_estimate"] = f.rho * u_pore * D_h / f.mu
        out["Pe_axial_L_pore_estimate"] = f.rho * f.cp * u_pore * g.L / f.k
        out["hydraulic_diameter_flat_m"] = D_h
        out["definitions"]["u_pore_volume_porosity_estimate"] = (
            "U / (1 - pi d^2/(4 l^2)); flat-interface volume porosity, NOT an intrinsic section velocity")
    area = c.liquid_section_area_per_width(sol.Pc_faces[0]) if hasattr(c, "liquid_section_area_per_width") else None
    if area is not None:
        u_int = Q[0] / (g.W * float(area[0]))
        out["u_intrinsic_section_feed_m_s"] = u_int
        out["Re_h_intrinsic_section"] = f.rho * u_int * g.h / f.mu
    else:
        out["u_intrinsic_section_feed_m_s"] = None
        out["definitions"]["u_intrinsic_section"] = "UNAVAILABLE: no A_liquid_section in the closure tables"
    out["Ca_superficial"] = f.mu * U_feed / f.sigma
    if l is not None:
        out["Bo_pitch"] = (f.rho - f.rho_v) * G_STANDARD * l**2 / f.sigma
        out["We_pitch_superficial"] = f.rho * U_feed**2 * l / f.sigma
        grad = f.mu * Q / (g.W * c.tau(sol.Pc_faces))
        out["max_dPc_dx_Pa_m"] = float(np.max(grad))
        out["max_cell_pressure_change_Pa"] = float(np.max(grad) * l)
        cap = sol.interval.Pc_cap - problem.Pc_feed if sol.interval else None
        if cap:
            out["max_cell_pressure_change_over_capillary_scale"] = float(np.max(grad) * l / cap)
            out["definitions"]["capillary_scale"] = "Pc_cap - Pc_feed (nonzero; local Pc is not used as a divisor)"
    q_e = sol.q_e
    u_e = q_e / (f.rho * f.h_fg)
    out["u_e_max_m_s"] = float(np.max(u_e))
    out["Pe_vertical_max"] = float(f.rho * f.cp * np.max(u_e) * g.h / f.k)
    out["Pe_axial_L_feed"] = float(f.rho * f.cp * U_feed * g.L / f.k)
    out["definitions"]["Pe_axial_L_feed"] = "rho cp U_feed L / k_l (superficial velocity, half length L)"
    out["Ja_max_superheat"] = f.cp * dT_max / f.h_fg
    # Omitted sensible heat: liquid warming from T_v to roughly the local T_e (explicit approximation).
    Q_in_total = float(np.sum(sol.q_in * m.dx)) * g.W
    if Q_in_total > 0:
        mdot_feed = f.rho * Q[0]
        out["omitted_sensible_over_latent_global_estimate"] = mdot_feed * f.cp * dT_max / Q_in_total
        dTe_dx = np.gradient(theta_e, m.xc)
        Qc = 0.5 * (Q[:-1] + Q[1:])
        adv = f.rho * f.cp * Qc * np.abs(dTe_dx)
        evap = g.W * np.maximum(q_e, 1e-300)
        out["axial_enthalpy_over_evaporation_local_max_estimate"] = float(np.max(adv / evap))
        out["definitions"]["axial_enthalpy_estimate"] = (
            "rho cp Q |dT_e/dx| / (W q_e) with liquid temperature approximated by T_e; approximate")
    u_v = float(np.max(q_e)) / (f.rho_v * f.h_fg)
    out["vapor_Mach_max"] = u_v / math.sqrt(f.gamma_v * f.R_specific * f.T_v)
    out["Knudsen"] = None
    out["definitions"]["Knudsen"] = "not computed: the molecular diameter is unverified (breck1974, ASSUMED)"
    if l is not None:
        dTe_dx = np.gradient(theta_e, m.xc)
        out["cell_to_cell_dT_over_superheat_max"] = float(np.max(np.abs(dTe_dx) * l / np.maximum(theta_e, 1e-300)))
    out["max_superheat_T_e_K"] = dT_max
    out["max_superheat_T_source_K"] = float(np.max(sol.slab.theta_source))
    # supported operating ranges of the closure
    viol = []
    lim = c.limits
    if "max_dPc_dx_Pa_m" not in out:
        grad = f.mu * Q / (g.W * c.tau(sol.Pc_faces))
        out["max_dPc_dx_Pa_m"] = float(np.max(grad))
    glo, ghi = lim.supported_gradient_Pa_m
    if out["max_dPc_dx_Pa_m"] > ghi:
        viol.append(f"max dPc/dx {out['max_dPc_dx_Pa_m']:.4g} Pa/m above supported {ghi:.4g} Pa/m")
    qlo, qhi = lim.supported_heat_flux_W_m2
    if float(np.max(q_e)) > qhi or float(np.min(q_e)) < qlo:
        viol.append(f"local q_e in [{np.min(q_e):.4g}, {np.max(q_e):.4g}] W/m^2 outside supported "
                    f"[{qlo:.4g}, {qhi:.4g}] W/m^2 (closure load applicability)")
    out["supported_ranges"] = {"pressure_gradient_Pa_m": list(lim.supported_gradient_Pa_m),
                               "heat_flux_W_m2": list(lim.supported_heat_flux_W_m2)}
    out["violations"] = viol
    out["within_supported_ranges"] = not viol
    return out


def reference_scales(q, L, h, d, l, rho, mu, h_fg, cp, k, dT=10.0) -> dict:
    """Audit reference scales for a uniform load q with the legacy flat geometry (estimates only)."""
    U = q * L / (rho * h_fg * h)
    porosity = 1 - math.pi * d**2 / (4 * l**2)
    u_pore = U / porosity
    u_e = q / (rho * h_fg)
    return {
        "U_feed_m_s": U,
        "u_pore_volume_porosity_m_s": u_pore,
        "Re_h_superficial": rho * U * h / mu,
        "Re_h_pore": rho * u_pore * h / mu,
        "Pe_vertical": rho * cp * u_e * h / k,
        "Pe_axial_L_superficial": rho * cp * U * L / k,
        "Pe_axial_L_pore": rho * cp * u_pore * L / k,
        "Ja": cp * dT / h_fg,
    }
