"""Digitize Fig. 7 of the Vaartstra thesis: evaporative flux q''_e vs interface superheat, a = 0.052.

Source: vaartstra2019_thesis, p. 23, Fig. 7 ("Evaporative heat flux as a function of superheat at the
liquid-vapor interface as calculated by solving the linearized BTE for a = 0.052"). The figure is a
raster image. The plot frame is located from the long dark frame lines, then the curve is read as the
mean dark-pixel row per column. Axes: 0-10 C on x, 0-500 W/cm^2 on y (as printed).

The figure does not state T_sat. The thesis unit-cell studies are at 100 C (pure vapor), and at the
3 kPa / 24 C state of the Adera experiments the same law would give fluxes ~20x lower, so 100 C is the
only reading consistent with the plotted magnitudes (PROJECT_DERIVED inference).

Output: data/literature/vaartstra2019_thesis_fig7_digitized.csv (PROJECT_DERIVED). Needs the
dev-only PyMuPDF package and the held PDF (config/sources.yaml: vaartstra2019_thesis).

Usage:
    python analysis/digitize_vaartstra_thesis_fig7.py
"""

from __future__ import annotations

import numpy as np
import yaml

from _common import REPO

PAGE_INDEX = 22  # 0-based index of thesis page 23
DPI = 300
DARK = 160
OUT = REPO / "data" / "literature" / "vaartstra2019_thesis_fig7_digitized.csv"


def frame(dark: np.ndarray, x_lo: int, x_hi: int, y_lo: int, y_hi: int, min_count: int = 300):
    """Mean row/column index of the frame lines (thick dark lines inside the search window)."""
    rows = [i for i in range(y_lo, y_hi) if dark[i, x_lo:x_hi].sum() > min_count]
    cols = [j for j in range(x_lo, x_hi) if dark[y_lo:y_hi, j].sum() > min_count]
    top = np.mean([r for r in rows if r < (y_lo + y_hi) / 2])
    bottom = np.mean([r for r in rows if r > (y_lo + y_hi) / 2])
    left = np.mean([c for c in cols if c < (x_lo + x_hi) / 2])
    right = np.mean([c for c in cols if c > (x_lo + x_hi) / 2])
    return top, bottom, left, right


def main() -> None:
    import fitz  # PyMuPDF, dev-only

    sources = yaml.safe_load((REPO / "config" / "sources.yaml").read_text(encoding="utf-8"))
    pdf = sources["vaartstra2019_thesis"]["local_path"]
    pix = fitz.open(pdf)[PAGE_INDEX].get_pixmap(dpi=DPI, colorspace=fitz.csGRAY)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    dark = img < DARK
    top, bottom, left, right = frame(dark, 700, 1800, 150, 1150)
    margin = 20
    rows = []
    for x in range(int(left) + margin, int(right) - margin, 4):
        ys = np.where(dark[int(top) + margin:int(bottom) - margin, x])[0]
        if len(ys) and ys.max() - ys.min() < 25:  # a single thin curve in this column
            y = ys.mean() + int(top) + margin
            rows.append((10 * (x - left) / (right - left), 500 * (bottom - y) / (bottom - top)))
    header = (
        "# vaartstra2019_thesis p. 23, Fig. 7, digitized by analysis/digitize_vaartstra_thesis_fig7.py. "
        "provenance: PROJECT_DERIVED (raster figure; reading uncertainty ~1% of the axis span). "
        "Linearized BTE with a = 0.052 and the Meland-Ytrehus treatment as implemented by the authors; "
        "T_sat not printed (100 C inferred).\n"
        "dT_int_K,q_e_W_cm2\n"
    )
    OUT.write_text(header + "".join(f"{a:.4f},{b:.2f}\n" for a, b in rows), encoding="utf-8")
    print(f"frame px: top {top:.1f} bottom {bottom:.1f} left {left:.1f} right {right:.1f}; {len(rows)} points -> {OUT}")


if __name__ == "__main__":
    main()
