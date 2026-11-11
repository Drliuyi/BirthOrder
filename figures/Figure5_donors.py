"""Render donor-level expression using supplied values and significance labels.

The module performs no normalization, model fitting or significance calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

GENES = ("FCRL6", "GZMA", "GZMH", "KLRD1")
CTS = ("NK", "CD8_T")
DATASETS = ("GSE289422", "GSE255566")
DISEASE = dict(zip(DATASETS, ("CAD", "T2D")))
GROUP_PALETTE = {"CTRL": "#9B7BC6", "CAD": "#D96C80", "T2D": "#5B9FC8"}
GROUP_COLORS = {"Control": GROUP_PALETTE["CTRL"], "CAD": GROUP_PALETTE["CAD"], "T2D": GROUP_PALETTE["T2D"]}
COHORT_COLORS = dict(zip(DATASETS, (GROUP_PALETTE["CAD"], GROUP_PALETTE["T2D"])))
MM = 1 / 25.4


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def jitter(sample, ct, gene):
    value = int(hashlib.sha256(f"{sample}|{ct}|{gene}".encode()).hexdigest()[:12], 16)
    return (value / (16**12 - 1) - .5) * .20


def configure_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    family = "DejaVu Sans"
    for candidate in ("Arial", "Helvetica"):
        try:
            font_manager.findfont(font_manager.FontProperties(family=candidate), fallback_to_default=False)
            family = candidate
            break
        except ValueError:
            pass
    plt.rcParams.update({
        "font.family": family, "font.size": 7.3, "axes.labelsize": 7.6,
        "axes.labelweight": "bold", "axes.titlesize": 8.3,
        "axes.linewidth": .65, "axes.edgecolor": "#202020",
        "xtick.labelsize": 7.0, "ytick.labelsize": 7.1,
        "xtick.major.size": 2.7, "ytick.major.size": 2.7,
        "xtick.major.width": .65, "ytick.major.width": .65,
        "text.color": "#202020", "axes.labelcolor": "#202020",
        "xtick.color": "#202020", "ytick.color": "#202020",
        "figure.facecolor": "white", "axes.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.bbox": None,
    })
    for key in ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize"):
        plt.rcParams[key] *= 1.30
    for key in ("text.color", "axes.labelcolor", "xtick.color", "ytick.color"):
        plt.rcParams[key] = "#000000"
    return family


def save_figure(fig, axes, out, stem, expected_mm):
    import matplotlib.pyplot as plt
    import numpy as np
    fig.canvas.draw()
    actual_mm = fig.get_size_inches() / MM
    require(np.allclose(actual_mm, expected_mm, atol=1e-9), "Unexpected canvas size")
    bounds = [list(map(float, ax.get_position().bounds)) for ax in axes]
    for dpi in (300, 600):
        path = out / f"{stem}_{dpi}dpi.png"
        require(not path.exists(), f"Refusing overwrite: {path}")
        fig.savefig(path, dpi=dpi, facecolor="white")
    for extension in ("pdf", "svg"):
        path = out / f"{stem}.{extension}"
        require(not path.exists(), f"Refusing overwrite: {path}")
        fig.savefig(path, facecolor="white")
    plt.close(fig)
    return {"figure": stem, "canvas_mm": actual_mm.tolist(), "axes_bounds": bounds}


def donor_plot(data, annotations, out, layout="vertical"):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle
    import numpy as np
    import pandas as pd
    if layout == "vertical":
        canvas_mm = (180, 151)
        fig, axes = plt.subplots(2, 1, figsize=(canvas_mm[0] * MM, canvas_mm[1] * MM))
        fig.subplots_adjust(left=.095, right=.985, bottom=.115, top=.895, hspace=.63)
        legend_center = .55
    elif layout == "horizontal":
        # Preserve each panel's physical width and height when arranging columns.
        canvas_mm = (360, 75.5)
        fig = plt.figure(figsize=(canvas_mm[0] * MM, canvas_mm[1] * MM))
        panel_height = (.895 - .115) / (2 + .63)
        axes = [fig.add_axes([(.095 + index) / 2, .115 * 2, (.985 - .095) / 2,
                              panel_height * 2]) for index in range(2)]
        legend_center = .55 / 2
    else:
        raise ValueError("Layout must be vertical or horizontal")
    xpositions = {ds: np.arange(4) + 4.5 * i for i, ds in enumerate(DATASETS)}
    for rownum, (ax, ct) in enumerate(zip(axes, CTS)):
        for ds in DATASETS:
            for gene, x in zip(GENES, xpositions[ds]):
                part = data.loc[data.dataset.eq(ds) & data.cell_type.eq(ct) & data.gene.eq(gene)]
                for group, offset in (("Control", -.18), (DISEASE[ds], .18)):
                    group_data = part.loc[part.group.eq(group)].sort_values("sample")
                    values = group_data.log2_one_plus_TMM_CPM.to_numpy()
                    color = GROUP_COLORS[group]
                    is_control = group == 'Control'
                    ax.boxplot([values], positions=[x + offset], widths=.29, manage_ticks=False,
                               patch_artist=True, showfliers=False, whis=1.5, zorder=2,
                               boxprops={"facecolor": matplotlib.colors.to_rgba(color, .22), "edgecolor": color, "linewidth": .7},
                               whiskerprops={"color": color, "linewidth": .65},
                               capprops={"color": color, "linewidth": .65},
                               medianprops={"color": color, "linewidth": 1.05})
                    xs = x + offset + np.array([jitter(s, ct, gene) for s in group_data["sample"]])
                    ax.scatter(xs, values, s=9, facecolors="none" if is_control else color,
                               edgecolors=color, alpha=.80, linewidth=.60 if is_control else 0, zorder=3)
                annotation = annotations.loc[annotations.dataset.eq(ds) & annotations.cell_type.eq(ct) & annotations.gene.eq(gene)].iloc[0]
                y = 11.60
                sig = str(annotation.raw_P_stars) if pd.notna(annotation.raw_P_stars) else ""
                require(sig in ("", "*", "**", "***"), "Invalid supplied significance label")
                if sig:
                    ax.plot([x - .34, x - .34, x + .34, x + .34],
                            [y - .14, y, y, y - .14], color="#3C3C3C", lw=.55, clip_on=False)
                    ax.text(x, y + .15, sig, ha="center", va="bottom", fontsize=12.75, fontweight="bold")
        ax.set_xlim(-.65, 8.15)
        ax.set_ylim(3, 12.35)
        ax.set_yticks([4, 6, 8, 10, 12])
        ax.set_xticks(np.concatenate(list(xpositions.values())), list(GENES) * 2, fontstyle="italic")
        ax.tick_params(axis="x", length=0, pad=4)
        ax.set_ylabel(r"TMM expression [log$_2$(1 + CPM)]", labelpad=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(0, 1.075, "NK cells" if ct == "NK" else "CD8 T cells", transform=ax.transAxes,
                fontsize=11.05, fontweight="bold", ha="left", va="bottom")
        for i, ds in enumerate(DATASETS):
            ax.add_patch(Rectangle((.5 * i, -.27), .5, .12, transform=ax.transAxes,
                                   fc=COHORT_COLORS[ds], ec="none", clip_on=False))
            label = ("CAD · GSE289422  |  15 CTRL / 29 CAD" if i == 0
                     else "T2D · GSE255566  |  3 CTRL / 3 T2D")
            ax.text(.25 + .5 * i, -.21, label, transform=ax.transAxes, color="white",
                    fontsize=8.645, fontweight="bold", va="center", ha="center", clip_on=False)
    fig.legend(handles=[Patch(facecolor=matplotlib.colors.to_rgba(GROUP_PALETTE[group], .22),
                             edgecolor=GROUP_PALETTE[group], label=group)
                        for group in ("CTRL", "CAD", "T2D")],
               loc="upper center", bbox_to_anchor=(legend_center, .995), ncol=3, frameon=False,
               fontsize=9.62, handlelength=1.6, columnspacing=2.0)
    audit = save_figure(fig, axes, out, "Figure5_donors", canvas_mm)
    b = np.array(audit["axes_bounds"])
    shared_edge = 0 if layout == "vertical" else 1
    require(np.ptp(b[:, shared_edge]) < 1e-10 and np.ptp(b[:, 2]) < 1e-10 and np.ptp(b[:, 3]) < 1e-10,
            "Donor panels do not share baselines or outer edges and dimensions")
    audit.update(layout=layout, single_coordinate_per_celltype=True, panel_widths_heights_aligned=True,
                 shared_y_limits=[3, 12.35], gene_order=list(GENES),
                 group_colors=GROUP_COLORS, cohort_colors=COHORT_COLORS,
                 zero_reference_line=False, model_refits=0, new_tests=0)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--layout", choices=["vertical", "horizontal"], default="vertical",
                        help="Arrange NK and CD8 T panels in rows or columns")
    args = parser.parse_args()
    inp, out = [Path(path).expanduser().resolve() for path in (args.input_root, args.output_root)]
    files = [inp / "donor_expression_plot_data.csv", inp / "frozen_P_annotations.csv"]
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    if out.exists():
        raise FileExistsError(f"Output root must be new: {out}")
    if out.is_relative_to(inp):
        raise ValueError("Output root must be outside the input directory")
    hashes = {path.name: sha(path) for path in files}
    import pandas as pd
    data = pd.read_csv(files[0])
    annotations = pd.read_csv(files[1])
    require(len(data) == 400 and len(annotations) == 16, "Unexpected plot input grid")
    require("raw_P_stars" in annotations, "Supplied significance labels are required")
    require(not annotations.duplicated(["dataset", "cell_type", "gene"]).any(), "Repeated annotation")
    out.mkdir(parents=True)
    os.environ.setdefault("MPLCONFIGDIR", str(out / "matplotlib_cache"))
    configure_style()
    audit = donor_plot(data, annotations, out, layout=args.layout)
    require(hashes == {path.name: sha(path) for path in files}, "Input changed during rendering")
    with (out / "RENDER_MANIFEST.json").open("x") as handle:
        json.dump({"input_sha256": hashes, "palette": GROUP_PALETTE, "geometry": audit,
                   "significance_labels": "supplied unchanged", "statistics_recomputed": 0}, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
