"""Reproducible run directories: environment, provenance snapshot, profiles, fields and plots.

Every run is written to ``<output>/<run_id>/`` (plots to ``plots/device/<run_id>/``) where run_id is
a UTC timestamp plus a short commit hash, made unique if needed; historical runs are never
overwritten. Diagnostics that are unavailable are written as null with a reason, never filled.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from ..provenance import REPO_ROOT

VERIFICATION_BANNER = "VERIFICATION ONLY - manufactured or analytical inputs; not a physical prediction"
PACKAGES = ("numpy", "scipy", "pandas", "matplotlib", "PyYAML", "CoolProp", "nanofab", "pytest")


def rel_path(p) -> str:
    """Repository-relative POSIX path when inside the repository, else the absolute path."""
    p = Path(p).resolve()
    return p.relative_to(REPO_ROOT).as_posix() if p.is_relative_to(REPO_ROOT) else str(p)


def git_state() -> dict:
    def run(*args):
        try:
            return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:  # noqa: BLE001
            return None
    status = run("status", "--porcelain")
    return {"commit": run("rev-parse", "HEAD"), "branch": run("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(status), "status_porcelain": status.splitlines() if status else []}


def environment() -> dict:
    versions = {}
    for p in PACKAGES:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            versions[p] = None
    return {"python": sys.version, "executable": sys.executable, "platform": platform.platform(),
            "packages": versions}


def new_run_dir(output: Path, label: str = "", plot_root: Path | None = None) -> tuple[str, Path, Path]:
    g = git_state()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    base = f"{stamp}_{(g['commit'] or 'nogit')[:8]}{'_' + label if label else ''}"
    run_id, k = base, 1
    output = Path(output)
    while (output / run_id).exists():
        run_id = f"{base}_{k}"
        k += 1
    run_dir = output / run_id
    run_dir.mkdir(parents=True)
    plot_dir = Path(plot_root or REPO_ROOT / "plots" / "device") / run_id
    plot_dir.mkdir(parents=True, exist_ok=True)
    return run_id, run_dir, plot_dir


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and not np.isfinite(obj):
        return str(obj)
    if hasattr(obj, "value") and obj.__class__.__name__.endswith("Status"):
        return obj.value
    if isinstance(obj, Path):
        return obj.as_posix()
    return obj


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(to_jsonable(obj), indent=2), encoding="utf-8")


def run_header(argv, model_mode: str, verification_only: bool) -> dict:
    return {
        "command_line": " ".join([Path(sys.executable).name] + list(argv)),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git": git_state(),
        "environment": environment(),
        "model_mode": model_mode,
        "verification_only": verification_only,
        "banner": VERIFICATION_BANNER if verification_only else None,
    }


def profiles_frame(problem, sol) -> pd.DataFrame:
    m = sol.mesh
    g = problem.geometry
    xf = m.x_faces
    Pc_c = sol.Pc_centers
    Qc = 0.5 * (sol.Q_faces[:-1] + sol.Q_faces[1:])
    df = pd.DataFrame({
        "x_m": m.xc, "x_over_L": m.xc / g.L,
        "Pc_Pa": Pc_c, "p_g_Pa": -Pc_c,
        "Q_cell_avg_m3_s": Qc, "U_superficial_m_s": Qc / (g.W * g.h),
        "q_in_W_m2": sol.q_in, "q_e_W_m2": sol.q_e, "h_e_W_m2K": sol.h_e,
        "T_e_K": problem.liquid.T_v + sol.slab.theta_e, "T_source_K": problem.liquid.T_v + sol.slab.theta_source,
        "dT_e_K": sol.slab.theta_e, "dT_source_K": sol.slab.theta_source,
        "admissibility_margin_Pa": sol.interval.Pc_cap - Pc_c,
    })
    area = problem.closure.liquid_section_area_per_width(Pc_c)
    df["u_intrinsic_section_m_s"] = Qc / (g.W * area) if area is not None else np.nan
    faces = pd.DataFrame({"x_face_m": xf, "Q_face_m3_s": sol.Q_faces, "Pc_face_Pa": sol.Pc_faces})
    return df, faces


def write_solution(run_dir: Path, plot_dir: Path, problem, sol, tag: str, verification_only: bool,
                   write_field: bool = True) -> dict:
    prof, faces = profiles_frame(problem, sol)
    banner = f"# {VERIFICATION_BANNER}\n" if verification_only else ""
    files = {}
    for name, df in ((f"{tag}_profiles.csv", prof), (f"{tag}_faces.csv", faces)):
        p = run_dir / name
        p.write_text(banner + df.to_csv(index=False, float_format="%.12g"), encoding="utf-8")
        files[name] = p.name
    if prof["u_intrinsic_section_m_s"].isna().all():
        files["u_intrinsic_note"] = "u_intrinsic_section_m_s is UNAVAILABLE (no A_liquid_section in the tables)"
    hist = pd.DataFrame(sol.history)
    hist.to_csv(run_dir / f"{tag}_convergence_history.csv", index=False)
    files["history"] = f"{tag}_convergence_history.csv"
    if write_field:
        m = sol.mesh
        X, Z = np.meshgrid(m.xc, m.zc, indexing="ij")
        fld = pd.DataFrame({"x_m": X.ravel(), "z_m": Z.ravel(),
                            "T_K": (problem.liquid.T_v + sol.slab.theta).ravel()})
        (run_dir / f"{tag}_substrate_T_field.csv").write_text(banner + fld.to_csv(index=False, float_format="%.10g"),
                                                               encoding="utf-8")
        files["field"] = f"{tag}_substrate_T_field.csv"
    files["plots"] = plot_solution(plot_dir, problem, sol, tag, verification_only)
    return files


def plot_solution(plot_dir: Path, problem, sol, tag: str, verification_only: bool) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    prof, _ = profiles_frame(problem, sol)
    xi = prof["x_over_L"]
    title_sfx = "\nVERIFICATION ONLY" if verification_only else ""
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    ax[0, 0].plot(xi, prof["Pc_Pa"])
    ax[0, 0].axhline(sol.interval.Pc_cap, ls="--", c="r", lw=0.8, label="admissible cap")
    ax[0, 0].set(xlabel="x/L [-]", ylabel="Pc = p_v - p_l [Pa]")
    ax[0, 0].legend()
    ax[0, 1].plot(xi, prof["h_e_W_m2K"] / 1e3)
    ax[0, 1].set(xlabel="x/L [-]", ylabel="h_e [kW/(m^2 K)]")
    ax[1, 0].plot(xi, prof["dT_e_K"], label="T_e - T_v (base plane)")
    ax[1, 0].plot(xi, prof["dT_source_K"], label="T_source - T_v")
    ax[1, 0].set(xlabel="x/L [-]", ylabel="superheat [K]")
    ax[1, 0].legend()
    ax[1, 1].plot(xi, prof["q_in_W_m2"] / 1e4, label="q_in''")
    ax[1, 1].plot(xi, prof["q_e_W_m2"] / 1e4, label="q_e''")
    ax[1, 1].set(xlabel="x/L [-]", ylabel="heat flux [W/cm^2]")
    ax[1, 1].legend()
    fig.suptitle(f"{tag}: load {sol.load/1e4:.6g} W/cm^2 (half-device mean){title_sfx}")
    fig.tight_layout()
    p1 = plot_dir / f"{tag}_profiles.png"
    fig.savefig(p1, dpi=120)
    plt.close(fig)
    m = sol.mesh
    fig, axf = plt.subplots(figsize=(10, 3))
    im = axf.pcolormesh(m.x_faces / problem.geometry.L, m.z_faces * 1e6, sol.slab.theta.T, shading="flat")
    fig.colorbar(im, ax=axf, label="T - T_v [K]")
    axf.set(xlabel="x/L [-]", ylabel="z [um]", title=f"{tag}: substrate temperature{title_sfx}")
    fig.tight_layout()
    p2 = plot_dir / f"{tag}_substrate_field.png"
    fig.savefig(p2, dpi=120)
    plt.close(fig)
    return [rel_path(p1), rel_path(p2)]
