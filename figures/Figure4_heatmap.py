#!/usr/bin/env python3
"""Render Figure 4b from saved aggregate estimates; perform no statistical tests.

Required CSV columns:
  --heatmap: model_id, domain_id, outcome, point_per_1000, q_BH
  --joint: model_id, point_per_1000, CI95_low_per_1000,
           CI95_high_per_1000, mediated_proportion_percent

The three model IDs and six domain IDs below define the manuscript layout.
The heatmap has three effect columns (IHD, T2D, Delta); DEP is not displayed.
Saved q values determine cell stars only. Missing joint q values are not
interpreted as evidence of non-significance; no joint stars are drawn.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.patches import Rectangle
import numpy as np

MODELS = [("S44_BASE11", "11 Risk Factors", "#555555"),
          ("S44_NMR_PLUS40", "+29 Metabolites", "#A46A2E"),
          ("S44_PROT_PLUS15", "+4 Proteins", "#325B87")]
DOMAINS = [("SMOKING", "Smoking", "#C49A9E"),
           ("BODY_COMPOSITION", "Body composition", "#A4BEBC"),
           ("CLINICAL_LIPIDS", "Clinical lipids", "#B3BFCF"),
           ("GLYCAEMIA", "Baseline glycaemia", "#D5B88E"),
           ("NMR_ALL", "Metabolites", "#8CAB9F"),
           ("CANDIDATE_PROTEINS", "Proteins", "#A899B7")]
MODEL_DOMAINS = {
    "S44_BASE11": tuple(d[0] for d in DOMAINS[:4]),
    "S44_NMR_PLUS40": tuple(d[0] for d in DOMAINS[:5]),
    "S44_PROT_PLUS15": tuple(d[0] for d in DOMAINS[:4]) + ("CANDIDATE_PROTEINS",),
}
W, H, GRID_W, CELL_GAP, RH = 180., 98., 39., .45, 10.8
CELL_W = (GRID_W - 2 * CELL_GAP) / 3
INK, BLUE, RED, REMAIN = "#202124", "#3E7198", "#A64B43", "#DEE2E5"


def read_rows(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        return list(reader)


def number(value):
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError("Plotting requires finite saved estimates")
    return parsed


def fmt(value, digits=2):
    return f"{number(value):.{digits}f}".replace("-", "−")


def saved_q(value):
    if value is None or str(value).strip().lower() in ("", "na", "nan", "none"):
        return None
    parsed = number(value)
    if not 0 <= parsed <= 1:
        raise ValueError("Saved q must be between zero and one")
    return parsed


def validate(data, joint):
    keys = [(r["model_id"], r["domain_id"], r["outcome"]) for r in data]
    expected = {(m, d, e) for m, ds in MODEL_DOMAINS.items()
                for d in ds for e in ("IHD", "T2D", "Delta")}
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError("Heatmap rows must uniquely cover the manuscript model/domain/effect layout")
    models = [r["model_id"] for r in joint]
    if len(models) != len(set(models)) or set(models) != set(MODEL_DOMAINS):
        raise ValueError("Joint CSV must contain one row per manuscript model")
    for row in data:
        if abs(number(row["point_per_1000"])) > 1.5:
            raise ValueError("Saved heatmap effect exceeds the fixed colour-scale limits")
        q = saved_q(row["q_BH"])
        if row.get("significance_marker", "") not in ("", "NA"):
            marker = row["significance_marker"].strip().lower()
            if marker not in ("true", "false") or (q is not None and q < .05) != (marker == "true"):
                raise ValueError("Saved significance marker disagrees with saved q")
    for row in joint:
        lo, point, hi = (number(row[k]) for k in
                         ("CI95_low_per_1000", "point_per_1000", "CI95_high_per_1000"))
        if not 0 <= lo <= point <= hi <= 5:
            raise ValueError("Joint interval must fit the fixed common 0–5 forest scale")
        if not 0 <= number(row["mediated_proportion_percent"]) <= 100:
            raise ValueError("Mediated proportion cannot be drawn as a 0–100% bar")


def render(data, joint, out, font_family="Arial"):
    """Keep the manuscript's categorical-band geometry and colour mapping."""
    validate(data, joint)
    plt.rcParams.update({"font.family": font_family, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "mathtext.fontset": "custom",
                         "mathtext.rm": font_family, "mathtext.it": font_family + ":italic",
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(W / 25.4, H / 25.4), dpi=180)
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, W), ylim=(H, 0))
    ax.set_axis_off()
    texts, stars, geometry, star_cells = [], [], [], []

    def text(x, y, label, size=8, bold=False, ha="left"):
        artist = ax.text(x, y, label, fontsize=size, fontweight="bold" if bold else "normal",
                         ha=ha, va="center", color=INK, zorder=5)
        texts.append(artist)
        return artist

    def rect(x, y, width, height, colour):
        ax.add_patch(Rectangle((x, y), width, height, facecolor=colour, edgecolor="none", zorder=1))

    cmap = colors.LinearSegmentedColormap.from_list("Blue_White_Brick", [BLUE, "#FFFFFF", RED], N=513)
    norm = colors.TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5)
    cells = {(r["model_id"], r["domain_id"], r["outcome"]): r for r in data}
    joint_by_model = {r["model_id"]: r for r in joint}
    for k, (mid, title, cohort_colour) in enumerate(MODELS):
        gx, gy = 7 + 45 * k, 15.8
        centre = gx + GRID_W / 2
        rect(gx, 13.5, GRID_W, 2.3, cohort_colour)
        for j, label in enumerate(("IHD", "T2D", "ΔIE")):
            text(gx + j * (CELL_W + CELL_GAP) + CELL_W / 2, 9.6, label, 8.2, True, ha="center")
        gs = [d for d in DOMAINS if d[0] in MODEL_DOMAINS[mid]]
        for i, (gid, _, category_colour) in enumerate(gs):
            y = gy + i * RH
            rect(gx - 2.8, y, 2.5, RH - .4, category_colour)
            for j, endpoint in enumerate(("IHD", "T2D", "Delta")):
                row = cells[mid, gid, endpoint]
                x = gx + j * (CELL_W + CELL_GAP)
                val, q = number(row["point_per_1000"]), saved_q(row["q_BH"])
                rect(x, y, CELL_W, RH - .4, cmap(norm(val)))
                text(x + CELL_W / 2, y + (RH - .4) / 2, fmt(val), 8.2, ha="center")
                if q is not None and q < .05:
                    mark = text(x + CELL_W - 1.15, y + (RH - .4) / 2 - 1.1, "*", 12.8, True, ha="center")
                    star_cells.append((mark, x, y, CELL_W, RH - .4))
                    stars.append({"model_id": mid, "domain_id": gid, "outcome": endpoint})
        row = joint_by_model[mid]
        fraction = number(row["mediated_proportion_percent"]) / 100
        bar_y, bar_h = gy + len(gs) * RH + .4, 4.8
        rect(gx, bar_y, GRID_W, bar_h, REMAIN)
        rect(gx, bar_y, GRID_W * fraction, bar_h, cohort_colour)
        text(centre, bar_y + bar_h / 2, f"{100 * fraction:.1f}%", 8.2, True, ha="center")
        text(gx, bar_y + 9.0, "Joint ΔIE", 7.7, True)
        forest_y = bar_y + 15.0
        estimate = number(row["point_per_1000"])
        ci_low, ci_high = number(row["CI95_low_per_1000"]), number(row["CI95_high_per_1000"])

        def fx(value):
            return gx + value / 5 * GRID_W

        ax.plot([fx(0), fx(0)], [forest_y - 2.4, forest_y + 2.6], color="#A7ADB2", lw=.55, ls=(0, (2, 2)), zorder=1)
        ax.plot([fx(ci_low), fx(ci_high)], [forest_y, forest_y], color=cohort_colour, lw=.95, zorder=2)
        for value in (ci_low, ci_high):
            ax.plot([fx(value), fx(value)], [forest_y - .8, forest_y + .8], color=cohort_colour, lw=.8)
        ax.plot([fx(estimate)], [forest_y], marker="o", markersize=4.5,
                markerfacecolor=cohort_colour, markeredgecolor=cohort_colour, zorder=3)
        ax.plot([gx, gx + GRID_W], [forest_y + 2.6, forest_y + 2.6], color="#C7CCD0", lw=.4)
        for tick in (0, 5):
            ax.plot([fx(tick), fx(tick)], [forest_y + 2.6, forest_y + 3.1], color="#A7ADB2", lw=.4)
            text(fx(tick), forest_y + 5.0, str(tick), 6.2, ha="center")
        text(centre, bar_y + 23.2, fmt(estimate, 3) + " (95% CI " + fmt(ci_low, 3) + "–" + fmt(ci_high, 3) + ")", 7.1, ha="center")
        geometry.append({"model_id": mid, "heatmap_x_mm": gx, "heatmap_width_mm": GRID_W,
                         "proportion_bar_x_mm": gx, "proportion_bar_width_mm": GRID_W,
                         "heatmap_y_mm": gy, "top_band_bottom_mm": 15.8,
                         "row_height_mm": RH, "proportion_bar_y_mm": bar_y})

    lx = 143.5
    text(lx, 6.0, "Population", 8.1, True)
    for i, (_, label, colour) in enumerate(MODELS):
        cy = 11.5 + i * 4.7
        rect(lx, cy - 1.25, 2.7, 2.5, colour)
        text(lx + 4.2, cy, label, 7.1)
    text(lx, 27.0, "Domains", 8.1, True)
    for i, (_, label, colour) in enumerate(DOMAINS):
        cy = 32.2 + i * 4.2
        rect(lx, cy - 1.25, 2.7, 2.5, colour)
        text(lx + 4.2, cy, label, 7.0)
    text(lx, 60.0, "IE (per 1,000)", 8.1, True)
    cbx, cby, cbw, cbh = lx, 66.6, 3.4, 15.
    for i, val in enumerate(np.linspace(1.5, -1.5, 256)):
        rect(cbx, cby + i * cbh / 256, cbw, cbh / 256 + .002, cmap(norm(val)))
    for val in (1.5, 0, -1.5):
        yy = cby + (1.5 - val) / 3 * cbh
        ax.plot([cbx + cbw, cbx + cbw + .8], [yy, yy], color=INK, lw=.4)
        text(cbx + cbw + 1.8, yy, f"{val:g}".replace("-", "−"), 7.0)
    rx = 156.0
    text(rx - 1.0, 65.3, "Proportion mediated", 7.0, True)
    for cy, label, colour in ((70.6, "Mediated", None), (75.1, "Remaining", REMAIN)):
        if colour is None:
            for j, (_, _, cohort_colour) in enumerate(MODELS):
                rect(rx + j * 2.4 / 3, cy - 1.1, 2.4 / 3, 2.2, cohort_colour)
        else:
            rect(rx, cy - 1.1, 2.4, 2.2, colour)
        text(rx + 3.8, cy, label, 6.9)
    text(rx + 1.0, 83.8, "*", 11.5, True, ha="center")
    text(rx + 4.0, 83.8, r"BH $\mathit{q}$ < 0.05", 6.7)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for mark, x, y, width, height in star_cells:
        box = mark.get_window_extent(renderer).transformed(ax.transData.inverted())
        bx0, bx1 = sorted((box.x0, box.x1))
        by0, by1 = sorted((box.y0, box.y1))
        if min(bx0 - x, x + width - bx1, by0 - y, y + height - by1) <= .25:
            raise ValueError("A star is too close to its cell boundary")
    for i, artist in enumerate(texts):
        box = artist.get_window_extent(renderer)
        if box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height:
            raise ValueError(f"Clipped text: {artist.get_text()}")
        for other in texts[i + 1:]:
            if box.overlaps(other.get_window_extent(renderer)):
                raise ValueError(f"Text overlap: {artist.get_text()} / {other.get_text()}")
    out.mkdir(parents=True, exist_ok=False)
    for ext in ("pdf", "svg"):
        fig.savefig(out / ("Figure4_Heatmap." + ext))
    fig.savefig(out / "Figure4_Heatmap.png", dpi=600)
    fig.savefig(out / "Figure4_Heatmap_preview.png", dpi=180)
    plt.close(fig)
    return {"displayed_cells": len(data), "stars": stars, "geometry": geometry,
            "width_mm": W, "height_mm": H, "font": font_family,
            "text_overlap": [], "text_clipped": [], "new_P_q_calculations": 0,
            "new_model_fits": 0, "new_bootstrap_draws": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heatmap", type=Path, required=True)
    parser.add_argument("--joint", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="New output directory; never replaced")
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--font", type=Path, action="append", default=[], help="Optional font file; repeat for bold")
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Output directory must be new")
    for font in args.font:
        font_manager.fontManager.addfont(str(font))
    # Require an explicit available font to avoid silently changing geometry.
    font_manager.findfont(args.font_family, fallback_to_default=False)
    inputs = (args.heatmap, args.joint)
    pins = [hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs]
    data = read_rows(args.heatmap, ("model_id", "domain_id", "outcome", "point_per_1000", "q_BH"))
    joint = read_rows(args.joint, ("model_id", "point_per_1000", "CI95_low_per_1000",
                                  "CI95_high_per_1000", "mediated_proportion_percent"))
    audit = render(data, joint, args.out, args.font_family)
    if pins != [hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs]:
        raise RuntimeError("Input changed during rendering")
    audit["input_sha256"] = dict(zip(("heatmap", "joint"), pins))
    with (args.out / "Figure4_Heatmap_QC.json").open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(json.dumps({"output": str(args.out), "displayed_cells": len(data), "new_P_q_calculations": 0}))


if __name__ == "__main__":
    main()
