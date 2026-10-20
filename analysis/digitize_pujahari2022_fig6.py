"""Digitize pujahari2022 Fig. 6 (dryout heat flux vs pitch) from the held PDF.

Output: data/literature/pujahari2022_fig6_digitized.csv (provenance PROJECT_DERIVED: digitized from
a figure; see docs/SOURCE_AUTHORITY.md). Requires PyMuPDF (`pip install pymupdf`, dev-only).

Method: extract the embedded Fig. 6 bitmap (page 10), locate each panel's axes frame from long dark
pixel runs, calibrate with the frame = printed axis limits, classify pixels by curve colour, and
track each curve column-by-column from a seed point in a clean region (nearest pixel cluster to the
previous point), so overlapping text and legend samples are not picked up. Dashed gaps are linearly
interpolated. Reading uncertainty is about +-1 px, i.e. +-0.5% of the axis span.

Usage:
    python analysis/digitize_pujahari2022_fig6.py [--plot]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import REPO

from nanofab.provenance import load_sources

OUT = REPO / "data" / "literature" / "pujahari2022_fig6_digitized.csv"

# panel: frame search box (y0, y1, x0, x1) in the 1566 x 1712 bitmap, axis limits, seeds (x, y) in data units
PANELS = {
    "a": {"box": (0, 760, 0, 783), "xlim": (7, 19), "ylim": (2, 8), "h_um": 5, "zeta_mV": -100, "pitch_range": (8, 18),
          "seed_x": 11.0, "seed": {"none": 3.25, "5e3": 3.45, "1e4": 3.62, "5e4": 4.95, "1e5": 6.5}},
    "b": {"box": (0, 760, 783, 1566), "xlim": (7, 19), "ylim": (2, 12), "h_um": 5, "zeta_mV": -200, "pitch_range": (8, 18),
          "seed_x": 11.0, "seed": {"none": 3.25, "5e3": 3.62, "1e4": 4.0, "5e4": 6.5, "1e5": 9.3}},
    "c": {"box": (820, 1712, 0, 783), "xlim": (5, 50), "ylim": (30, 90), "h_um": 25, "zeta_mV": -100, "pitch_range": (10, 45),
          "seed_x": 27.0, "seed": {"none": 69.0, "5e3": 70.0, "1e4": 71.5, "5e4": 79.0, "1e5": 88.0},
          "seed2_x": 11.0, "seed2": {"none": 35.3, "5e3": 36.2, "1e4": 37.5, "5e4": 43.2, "1e5": 49.0}},
    "d": {"box": (820, 1712, 783, 1566), "xlim": (5, 50), "ylim": (20, 112), "h_um": 25, "zeta_mV": -200, "pitch_range": (10, 45),
          "seed_x": 27.0, "seed": {"none": 69.0, "5e3": 71.0, "1e4": 73.0, "5e4": 88.0, "1e5": 106.0},
          "seed2_x": 11.0, "seed2": {"none": 35.6, "5e3": 36.5, "1e4": 39.9, "5e4": 49.5, "1e5": 61.4}},
}
FIELDS = {"none": 0.0, "5e3": 5e3, "1e4": 1e4, "5e4": 5e4, "1e5": 1e5}


def colour_masks(im):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    return {
        "none": (r + g + b) < 200,
        "5e3": (r > 170) & (g < 60) & (b < 60),
        "1e4": (b > 170) & (r < 90) & (g < 90),
        "5e4": (r > 200) & (g > 40) & (g < 140) & (b > 130),
        "1e5": (r > 200) & (g > 50) & (g < 140) & (b < 50),
    }


def frame(dark, y0, y1, x0, x1):
    sub = dark[y0:y1, x0:x1]
    rows = [y0 + i for i, v in enumerate(sub.sum(axis=1)) if v > 0.6 * (x1 - x0)]
    cols = [x0 + j for j, v in enumerate(sub.sum(axis=0)) if v > 0.6 * (y1 - y0) - 100]
    return min(rows) + 1, max(rows) - 1, min(cols) + 1, max(cols) - 1  # inner edges


def track(mask, top, bot, left, right, xs, ys, seed_px, max_jump=6, margin=6):
    """Follow one curve from the seed column in both directions; returns {col: row}."""
    m = mask[top + margin:bot - margin, left + margin:right - margin]
    out = {}
    c0, y_seed = seed_px[0] - left - margin, seed_px[1] - top - margin
    best = None
    for dc in sorted(range(-25, 26), key=abs):
        rows = np.flatnonzero(m[:, c0 + dc])
        if rows.size:
            d = np.abs(rows - y_seed)
            if d.min() <= 30 and (best is None or d.min() < best[2]):
                best = (c0 + dc, float(rows[np.argmin(d)]), d.min())
    if best is None:
        raise RuntimeError(f"seed not found near {seed_px}")
    for direction in (1, -1):
        hist = [(best[0], best[1])]
        c = best[0] + direction
        miss = 0
        while 0 <= c < m.shape[1]:
            recent = hist[-12:]
            if len(recent) >= 4:
                xs_, ys_ = np.array(recent).T
                slope = np.polyfit(xs_, ys_, 1)[0]
            else:
                slope = 0.0
            pred = hist[-1][1] + slope * (c - hist[-1][0])
            rows = np.flatnonzero(m[:, c])
            if rows.size:
                runs = np.split(rows, np.flatnonzero(np.diff(rows) > 1) + 1)
                centers = np.array([r.mean() for r in runs])
                k = int(np.argmin(np.abs(centers - pred)))
                if abs(centers[k] - pred) <= max_jump + 0.25 * miss:
                    hist.append((c, centers[k]))
                    out[c + left + margin] = centers[k] + top + margin
                    miss = 0
                else:
                    miss += 1
            else:
                miss += 1
            if miss > 150:
                break
            c += direction
    return out


def clipped_fit(x, y, deg=6, clip=3.0, iters=8):
    """Polynomial fit with iterative sigma clipping (rejects legend samples / stray text)."""
    keep = np.ones(x.size, bool)
    for _ in range(iters):
        p = np.polyfit(x[keep], y[keep], deg)
        r = y - np.polyval(p, x)
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        new = np.abs(r) <= max(clip * s, 1e-9)
        if new.sum() == keep.sum() and np.all(new == keep):
            break
        keep = new
    return p, keep


def digitize(pdf_path):
    import pymupdf
    from PIL import Image
    import io

    doc = pymupdf.open(pdf_path)
    xref = doc[9].get_images()[0][0]
    pix = pymupdf.Pixmap(doc, xref)
    if pix.n - pix.alpha >= 4:
        pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
    im = np.asarray(Image.open(io.BytesIO(pix.tobytes("png"))).convert("RGB")).astype(int)
    dark = im.sum(axis=2) < 250
    masks = colour_masks(im)
    rows = []
    curves = {}
    for name, P in PANELS.items():
        top, bot, left, right = frame(dark, *P["box"])
        (x0, x1), (y0, y1) = P["xlim"], P["ylim"]
        to_px = lambda x, y: (int(round(left + (x - x0) / (x1 - x0) * (right - left))),  # noqa: E731
                              top + (y1 - y) / (y1 - y0) * (bot - top))
        for key, ys in P["seed"].items():
            pts = track(masks[key], top, bot, left, right, None, None, to_px(P["seed_x"], ys))
            if "seed2" in P:
                extra = track(masks[key], top, bot, left, right, None, None, to_px(P["seed2_x"], P["seed2"][key]))
                pts.update({c: v for c, v in extra.items() if c not in pts})
            cols = np.array(sorted(pts))
            x = x0 + (cols - left) / (right - left) * (x1 - x0)
            y = y1 - (np.array([pts[c] for c in cols]) - top) / (bot - top) * (y1 - y0)
            lo, hi = P["pitch_range"]
            inside = (x >= lo - 0.3) & (x <= hi + 0.3)
            x, y = x[inside], y[inside]
            p, keep = clipped_fit(x, y)

            def direct(l, p=p, mask=masks[key]):
                """Pixel read at column l (+-3 px): cluster nearest the extrapolated fit, within 3% of span."""
                c = int(round(left + (l - x0) / (x1 - x0) * (right - left)))
                guess = np.polyval(p, l)
                vals = []
                for cc in range(c - 3, c + 4):
                    rows_ = np.flatnonzero(mask[top + 8:bot - 8, cc]) + top + 8
                    if rows_.size:
                        for r_ in np.split(rows_, np.flatnonzero(np.diff(rows_) > 1) + 1):
                            vals.append(y1 - (r_.mean() - top) / (bot - top) * (y1 - y0))
                vals = np.array(vals)
                ok = vals[np.abs(vals - guess) < 0.03 * (y1 - y0)] if vals.size else vals
                return float(np.median(ok)) if ok.size else float("nan")

            curves[(name, key)] = (x[keep], y[keep], p, direct)
            rows.append({"panel": name, "h_um": P["h_um"], "zeta_mV": P["zeta_mV"], "E_V_per_m": FIELDS[key],
                         "x_min_um": x[keep].min(), "x_max_um": x[keep].max(), "n_points": int(keep.sum()),
                         "n_rejected": int((~keep).sum()),
                         "fit_rms_Wcm2": float(np.sqrt(np.mean((y[keep] - np.polyval(p, x[keep])) ** 2)))})
    return curves, pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    pdf = load_sources()["pujahari2022"]["local_path"]
    curves, coverage = digitize(pdf)
    print(coverage.to_string(index=False))
    sample = {"a": np.arange(8, 18.01, 1.0), "b": np.arange(8, 18.01, 1.0),
              "c": np.arange(10, 45.01, 2.5), "d": np.arange(10, 45.01, 2.5)}
    out = []
    for (name, key), (x, y, p, direct) in curves.items():
        P = PANELS[name]
        for l in sample[name]:
            if x.min() - 0.15 <= l <= x.max() + 0.15:
                q, method = float(np.polyval(p, l)), "fit"
            elif l >= x.min() - 2.0:
                q, method = direct(l), "direct"
            else:
                continue
            if np.isfinite(q):
                out.append({"panel": name, "d_um": 5, "h_um": P["h_um"], "L_mm": 5, "zeta_mV": P["zeta_mV"],
                            "E_V_per_m": FIELDS[key], "pitch_um": l, "q_dry_Wcm2": q, "method": method,
                            "reading_uncertainty_Wcm2": 0.005 * (P["ylim"][1] - P["ylim"][0])})
    df = pd.DataFrame(out).sort_values(["panel", "E_V_per_m", "pitch_um"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write("# pujahari2022 Fig. 6, digitized by analysis/digitize_pujahari2022_fig6.py. provenance: PROJECT_DERIVED "
                 "(source pujahari2022, Fig. 6 a-d). Values are read from a raster figure, not reported numbers.\n")
        df.to_csv(fh, index=False)
    print(f"wrote {OUT} ({len(df)} rows)")
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2, figsize=(9, 8))
        for a, name in zip(ax.ravel(), "abcd"):
            for key in FIELDS:
                x, y, p, _ = curves[(name, key)]
                a.plot(x, y, ".", ms=1.5, label=key)
                xx = np.linspace(x.min(), x.max(), 200)
                a.plot(xx, np.polyval(p, xx), "k-", lw=0.5)
            a.set_xlim(*PANELS[name]["xlim"])
            a.set_ylim(*PANELS[name]["ylim"])
            a.set_title(name)
            a.legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(REPO / "results" / "electroosmosis" / "pujahari2022_fig6_redigitized.png", dpi=150)


if __name__ == "__main__":
    main()
