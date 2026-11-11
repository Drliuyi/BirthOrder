#!/usr/bin/env python3
"""Render Figure 3d from saved single-mediator six-contrast aggregates.

--values columns: indicator_id, indicator_name_en, domain (Clinical/NMR/Olink),
  mediator_order, contrast, Delta_IE_per_1000, Delta_TE_per_1000, Pcontrast_percent.
--a-directions columns: indicator_id, sign (-1/+1), layer.

Exactly 18 mediators and the six contrast IDs below define this panel. Widths
encode absolute saved Delta_IE on one linear scale; inline numbers are saved
100*Delta_IE/Delta_TE percentages, not P values. Left-band colour encodes the
supplied a-reference sign. Right-ribbon colour encodes sign(a)*sign(Delta_IE),
a derived graphical code, NOT an independently fitted b coefficient. The
supplied layer must describe the provenance of each a-reference sign. No
models, resampling, P values or q values are calculated by this renderer.
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
from matplotlib.colors import to_rgba
from matplotlib import font_manager
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle

POS, NEG, CORONARY = "#E67B8D", "#20A9AE", "#B65F64"
DOMAIN_COLOURS = {"Clinical": "#E39D3D", "NMR": "#8CAB9F", "Olink": "#A899B7"}
PRIMARY = "411_minus_250.2"
CONTRASTS = [("411_minus_250", "IHD−DM"), ("411_minus_250.2", "IHD−T2D"),
             ("411.3_minus_250", "Angina−DM"), ("411.3_minus_250.2", "Angina−T2D"),
             ("411.4_minus_250", "CAS−DM"), ("411.4_minus_250.2", "CAS−T2D")]
RIBBON_TOTAL_IN = 9.0
MX, MW, TX, TW, SX, SW, W0 = 6.59, .18, 11.59, .16, 1.59, .18, 15.50


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


def load(values, directions):
    cells, mediators = {}, {}
    for row in values:
        iid, contrast = row["indicator_id"], row["contrast"]
        key = (iid, contrast)
        if key in cells:
            raise ValueError("Duplicate mediator/contrast")
        cells[key] = row
        mediator = (int(row["mediator_order"]), iid, row["indicator_name_en"], row["domain"])
        if iid in mediators and mediator != mediators[iid]:
            raise ValueError("Inconsistent mediator label, domain or order")
        if mediator[3] not in DOMAIN_COLOURS:
            raise ValueError("Unknown mediator domain")
        mediators[iid] = mediator
        if row.get("horizon_years") and number(row["horizon_years"]) != 11:
            raise ValueError("This panel is defined for the saved 11-year estimand")
        die, dte, share = (number(row[c]) for c in
                            ("Delta_IE_per_1000", "Delta_TE_per_1000", "Pcontrast_percent"))
        if die == 0 or dte == 0:
            raise ValueError("Zero-width ribbons or undefined shares cannot use this layout")
        if abs(100 * die / dte - share) >= 1e-4:
            raise ValueError("Saved contrast share does not match its saved numerator/denominator")
    meds = sorted(mediators.values())
    if len(meds) != 18 or len({m[0] for m in meds}) != 18:
        raise ValueError("Expected 18 uniquely ordered mediators")
    expected = {(m[1], c) for m in meds for c, _ in CONTRASTS}
    if set(cells) != expected:
        raise ValueError("Expected one saved estimate per mediator and manuscript contrast")
    a_dir, layers = {}, {}
    for row in directions:
        iid, sign = row["indicator_id"], int(row["sign"])
        if iid in a_dir or sign not in (-1, 1) or not row["layer"].strip():
            raise ValueError("Directions must be unique, signed and provenance-labelled")
        a_dir[iid], layers[iid] = sign, row["layer"]
    if set(a_dir) != set(mediators):
        raise ValueError("A-reference directions must cover exactly the plotted mediators")
    return meds, cells, a_dir, layers


def band(ax, points, widths, colour, alpha, gid):
    """Piecewise Bezier ribbon; exact vertical width at every attachment."""
    if len(points) != len(widths) or any(w < 0 for w in widths):
        raise ValueError("Invalid ribbon geometry")
    top = [(x, y + w / 2) for (x, y), w in zip(points, widths)]
    bottom = [(x, y - w / 2) for (x, y), w in zip(points, widths)]
    vertices, codes = [top[0]], [MPath.MOVETO]
    for start, end in zip(top, top[1:]):
        distance = (end[0] - start[0]) * .43
        vertices.extend([(start[0] + distance, start[1]), (end[0] - distance, end[1]), end])
        codes.extend([MPath.CURVE4] * 3)
    target_top = len(vertices) - 1
    vertices.append(bottom[-1]); codes.append(MPath.LINETO)
    target_bottom = len(vertices) - 1
    for start, end in zip(bottom[::-1], bottom[-2::-1]):
        distance = (start[0] - end[0]) * .43
        vertices.extend([(start[0] - distance, start[1]), (end[0] + distance, end[1]), end])
        codes.extend([MPath.CURVE4] * 3)
    source_bottom = len(vertices) - 1
    vertices.append(top[0]); codes.append(MPath.CLOSEPOLY)
    patch = PathPatch(MPath(vertices, codes), fc=to_rgba(colour, alpha), ec="none", lw=0, zorder=2)
    patch.set_gid(gid)
    patch.attachment_indices = (0, source_bottom, target_top, target_bottom)
    ax.add_patch(patch)
    return patch


def render(meds, cells, a_dir, layers, out, font_family="Arial"):
    ids = [m[1] for m in meds]
    labels, groups = {m[1]: m[2] for m in meds}, {m[1]: m[3] for m in meds}
    pc = {(iid, c): number(cells[iid, c]["Pcontrast_percent"]) for iid in ids for c, _ in CONTRASTS}
    die = {(iid, c): number(cells[iid, c]["Delta_IE_per_1000"]) for iid in ids for c, _ in CONTRASTS}
    scale = RIBBON_TOTAL_IN / sum(abs(v) for v in die.values())
    widths = {key: abs(value) * scale for key, value in die.items()}
    nh = {iid: sum(widths[iid, c] for c, _ in CONTRASTS) for iid in ids}
    rh = {iid: max(.34, nh[iid] + .12) for iid in ids}
    totals = {c: sum(widths[iid, c] for iid in ids) for c, _ in CONTRASTS}
    span = sum(rh.values()) + .05 * 17
    work = max(span, sum(totals.values()) + .16 * 5)
    height, ys = work + 2.35, {}
    top = height - .80
    cursor = top - (work - span) / 2
    for iid in ids:
        ys[iid] = cursor - rh[iid] / 2
        cursor -= rh[iid] + .05
    target = {}
    cursor = top - (work - sum(totals.values()) - .16 * 5) / 2
    for contrast, _ in CONTRASTS:
        target[contrast] = (cursor, cursor - totals[contrast])
        cursor -= totals[contrast] + .16
    plt.rcParams.update({"font.family": font_family, "text.color": "black",
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(W0, height), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, W0), ylim=(0, height)); ax.axis("off")
    artists = []

    def label(x, y, text, **kwargs):
        kwargs["fontsize"] = kwargs.get("fontsize", 10) * 1.5
        artist = ax.text(x, y, text, color=kwargs.pop("color", "black"), **kwargs)
        artists.append(artist)
        return artist

    label(MX - .45, height - .68, "Mediator", fontsize=11.5, weight="bold", ha="center")
    label(MX + MW + .07, height - .68, "Pcontrast\n(IHD−T2D), pp", fontsize=9.6, weight="bold", ha="left")
    label(TX + TW + .13, height - .68, "Contrast", fontsize=11.5, weight="bold", ha="left")
    centre = (max(ys.values()) + min(ys.values())) / 2
    source_height = sum(nh.values())
    ax.add_patch(Rectangle((SX, centre - source_height / 2), SW, source_height, fc="#9DA3A9", ec="none", zorder=5))
    label(SX - .13, centre, "Relative\nbirth order", ha="right", va="center", fontsize=12.5, weight="bold")
    label(SX - .13, centre - .40, "B = 0 → 1", ha="right", va="center", fontsize=9.5)
    ends = {c: target[c][0] for c, _ in CONTRASTS}
    max_width, start_x = max(widths.values()), MX + MW + .07
    lc = centre + source_height / 2
    for iid in ids:
        centre_y, left_width = ys[iid], nh[iid]
        a_colour = POS if a_dir[iid] > 0 else NEG
        band(ax, [(SX + SW, lc - left_width / 2), (MX, centre_y)],
             [left_width, left_width], a_colour, .36, "a-path-" + iid)
        lc -= left_width
        ax.add_patch(Rectangle((MX, centre_y - nh[iid] / 2), MW, nh[iid],
                               fc=DOMAIN_COLOURS[groups[iid]], ec="none", zorder=5))
        label(MX - .11, centre_y, labels[iid], ha="right", va="center", fontsize=11.5, weight="bold", zorder=10)
        root = centre_y + nh[iid] / 2
        for contrast, _ in CONTRASTS:
            key = (iid, contrast)
            width, ribbon_top = widths[key], root
            root -= width
            b_sign = a_dir[iid] * (1 if die[key] > 0 else -1)
            colour = POS if b_sign > 0 else NEG
            patch = band(ax, [(MX + MW, ribbon_top - width / 2),
                              (MX + MW + .16, ribbon_top - width / 2), (TX, ends[contrast] - width / 2)],
                         [width, width, width], colour, .46, "b-contrast-" + iid + "-" + contrast)
            patch.set_zorder(2 + .2 * (1 - width / max_width))
            ends[contrast] -= width
            vertices = ax.transData.transform(patch.get_path().vertices)
            st, sb, tt, tb = patch.attachment_indices
            if not all(math.isclose(abs(vertices[a, 1] - vertices[b, 1]) / fig.dpi, width, abs_tol=1e-10)
                       for a, b in ((st, sb), (tt, tb))):
                raise ValueError("Ribbon attachment width changed")
            if contrast == PRIMARY:
                artist = label(start_x, centre_y, f"{pc[key]:+.2f}", ha="left", va="center",
                               fontsize=12.4, weight="bold", zorder=12)
                artist.set_gid("PCONTRAST-" + iid + "-" + contrast)
        if not math.isclose(root, centre_y - nh[iid] / 2, abs_tol=1e-12):
            raise ValueError("Ribbon widths do not fill mediator node")
    for contrast, name in CONTRASTS:
        target_top, bottom = target[contrast]
        if not math.isclose(ends[contrast], bottom, abs_tol=1e-10):
            raise ValueError("Ribbon widths do not fill contrast node")
        ax.add_patch(Rectangle((TX, bottom), TW, target_top - bottom, fc=CORONARY, ec="none", lw=0, zorder=5))
        label(TX + TW + .13, (target_top + bottom) / 2, name, ha="left", va="center", fontsize=12.5, weight="bold")

    lx, tx, y = W0 - 2.05, W0 - 2.05 + .33, height - .95

    def head(text):
        nonlocal y
        label(lx, y, text, fontsize=7.8, weight="bold", ha="left", va="center")
        y -= .28

    def key(colour, text, square=False, alpha=None):
        nonlocal y
        h, width = (.10, .10) if square else (.045, .22)
        ax.add_patch(Rectangle((lx, y - h / 2), width, h,
                               fc=to_rgba(colour, 1 if square else .65 if alpha is None else alpha), ec="none"))
        label(tx, y, text, fontsize=7.5, ha="left", va="center")
        y -= .27

    head("Colour")
    label(lx, y, "Left band = a path\nRibbon = b_contrast", fontsize=5.6, ha="left", va="top", linespacing=1.2)
    y -= .40
    key(POS, "Positive direction", alpha=.46); key(NEG, "Negative direction", alpha=.46)
    y -= .08
    head("Mediator domain")
    for name, colour in DOMAIN_COLOURS.items():
        key(colour, name, square=True)
    y -= .12
    head("Width  |ΔIE|")
    label(lx, y - .02, "per 1,000", fontsize=6.7, ha="left", va="top", linespacing=1.2)
    y -= .60
    for value in (.2, .4, .6):
        h = value * scale
        ax.add_patch(Rectangle((lx, y - h / 2), .22, h, fc="#A5A8AC", ec="none"))
        label(tx, y, f"{value:g}", fontsize=7.5, ha="left", va="center")
        y -= max(.27, h + .09)
    y -= .10
    head("Inline numbers")
    label(lx, y, "Pcontrast (IHD−T2D), pp\nbold", fontsize=6.7, ha="left", va="top", linespacing=1.3)
    y -= .80
    head("Note")
    label(lx, y, "sign(b_contrast) = sign(a) × sign(ΔIE).\nsign(ΔIE) = sign(a) × sign(b_contrast).\n"
          "Number = Pcontrast = ΔIE/ΔTE.\nSingle-mediator; not additive.",
          fontsize=5.6, ha="left", va="top", linespacing=1.25)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    area = ax.get_window_extent(renderer)
    boxes = [(artist, artist.get_window_extent(renderer)) for artist in artists]
    for index, (artist, box) in enumerate(boxes):
        if not (area.contains(box.x0, box.y0) and area.contains(box.x1, box.y1)):
            raise ValueError(f"Clipped text: {artist.get_text()}")
        for other, other_box in boxes[index + 1:]:
            if box.overlaps(other_box):
                raise ValueError(f"Text overlap: {artist.get_text()} / {other.get_text()}")
    out.mkdir(parents=True, exist_ok=False)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / ("Figure3_SixContrasts." + ext), dpi=350, facecolor="white")
    fig.savefig(out / "Figure3_SixContrasts_preview.png", dpi=120, facecolor="white")
    plt.close(fig)
    return {"mediators": len(ids), "contrasts": len(CONTRASTS), "cells": len(cells),
            "size_inches": [W0, height], "width_total_inches": RIBBON_TOTAL_IN,
            "font": font_family, "a_direction_layers": sorted(set(layers.values())),
            "width_encodes": "abs(Delta_IE_per_1000)",
            "left_colour_encodes": "a-reference sign",
            "right_colour_encodes": "sign(a-reference) * sign(Delta_IE)",
            "inline_numbers": "saved Pcontrast_percent for IHD−T2D; not additive",
            "text_overlap": [], "text_clipped": [],
            "new_model_fits": 0, "new_bootstrap_draws": 0, "new_P_q_calculations": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--values", type=Path, required=True)
    parser.add_argument("--a-directions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="New output directory; never replaced")
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--font", type=Path, action="append", default=[], help="Optional font file; repeat for bold")
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Output directory must be new")
    for font in args.font:
        font_manager.fontManager.addfont(str(font))
    font_manager.findfont(args.font_family, fallback_to_default=False)
    paths = (args.values, args.a_directions)
    pins = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    values = read_rows(args.values, ("indicator_id", "indicator_name_en", "domain", "mediator_order",
                                    "contrast", "Delta_IE_per_1000", "Delta_TE_per_1000", "Pcontrast_percent"))
    directions = read_rows(args.a_directions, ("indicator_id", "sign", "layer"))
    audit = render(*load(values, directions), args.out, args.font_family)
    if pins != [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]:
        raise RuntimeError("Input changed during rendering")
    audit["input_sha256"] = dict(zip(("values", "a_directions"), pins))
    with (args.out / "Figure3_SixContrasts_QC.json").open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(json.dumps({"output": str(args.out), "cells": len(values), "new_P_q_calculations": 0}))


if __name__ == "__main__":
    main()
