#!/usr/bin/env python3
"""Render four-group cell-type UMAP and t-SNE panels from supplied coordinates.

Cell identifiers, coordinates, labels and display geometry are verified against
an external reference. This module does not compute an embedding or fit models.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def write_json_new(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['input-root', 'output-root']:
        env = 'SCRNA_STYLE_' + key.replace('-', '_').upper()
        parser.add_argument('--' + key, default=os.getenv(env), required=not os.getenv(env))
    parser.add_argument('--algorithm', choices=['UMAP', 'tSNE'], default=os.getenv('SCRNA_STYLE_ALGORITHM', 'UMAP'))
    parameters = vars(parser.parse_args())
    mode = parameters.pop('algorithm')
    roots = {key: Path(value).expanduser().resolve() for key, value in parameters.items()}
    for key, value in roots.items():
        print(key, value, flush=True)
    inp, out = [roots[key] for key in ['input_root', 'output_root']]
    if out.exists():
        raise FileExistsError(f'Output root must be new: {out}')
    if out.is_relative_to(inp):
        raise ValueError('Output root must be outside the input directory')

    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties, findfont
    from matplotlib.lines import Line2D

    sourcefiles = [inp / f'results/{ds}/coordinates.csv.gz' for ds in ['GSE289422', 'GSE255566']]
    sourcefiles += [inp / 'render_v2/metadata/cell_type_palette.csv',
                    inp / 'render_v2/metadata/INDEPENDENT_FINAL_VERIFICATION.json',
                    inp / 'render_v4_no_label_background/metadata/coordinate_translation_ledger.csv',
                    inp / 'render_v4_no_label_background/metadata/RENDER_MANIFEST.json']
    for path in sourcefiles:
        if not path.is_file():
            raise FileNotFoundError(f'Missing frozen input or source provenance: {path}')
    hashes = {str(path): sha(path) for path in sourcefiles}
    assert json.loads(sourcefiles[3].read_text())['verified']
    baseline_manifest = json.loads(sourcefiles[5].read_text())
    for folder in ['figures', 'metadata']:
        (out / folder).mkdir(parents=True, exist_ok=False)
    figdir, meta = out / 'figures', out / 'metadata'
    os.environ.setdefault('MPLCONFIGDIR', str(meta / 'matplotlib_cache'))

    scale = 1.30
    sizes = {'base': 8 * scale, 'axis': 8 * scale, 'legend': 7 * scale,
             'heading': 8.5 * scale, 'cell_count': 7.3 * scale,
             'cohort': 7.5 * scale, 'direct_cell_type': 7 * scale}
    cohort_colors = {'GSE289422': '#B65F64', 'GSE255566': '#5683A8'}
    font_family = 'DejaVu Sans'
    for candidate in ['Arial', 'Helvetica']:
        try:
            findfont(FontProperties(family=candidate), fallback_to_default=False)
            font_family = candidate
            break
        except ValueError:
            pass
    plt.rcParams.update({'font.family': 'sans-serif',
        'font.sans-serif': [font_family, 'Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': sizes['base'], 'axes.labelsize': sizes['axis'], 'axes.labelweight': 'bold',
        'legend.fontsize': sizes['legend'], 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'axes.linewidth': .975, 'axes.edgecolor': '#202020',
        'text.color': '#202020', 'axes.labelcolor': '#202020',
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'svg.hashsalt': 'core4-celltype'})
    dat = {ds: pd.read_csv(inp / f'results/{ds}/coordinates.csv.gz')
           for ds in ['GSE289422', 'GSE255566']}
    pal = pd.read_csv(sourcefiles[2])
    palette = {
        'NK':'#158F8D', 'CD8_T':'#3B66B0', 'CD4_T':'#8AA5D2',
        'T_other':'#A895C4', 'B':'#D6A33D', 'Plasma':'#E38A50',
        'Monocyte':'#B83D68', 'DC':'#D38FAB', 'pDC':'#E7BBCD',
        'Platelet':'#587E96', 'Erythroid':'#A56A46', 'Other/Unassigned':'#D7DDE3'}
    assert set(palette) == set(pal.cell_type)
    labels = dict(zip(pal.cell_type, pal.label))
    views = [('GSE289422', 'Control'), ('GSE289422', 'CAD'),
             ('GSE255566', 'Control'), ('GSE255566', 'T2D')]
    for ds, frame in dat.items():
        assert frame.cell_id.is_unique and set(frame.cell_type).issubset(palette)
        assert not frame.group.eq('Treatment').any()
    assert len(dat['GSE289422']) == 44000 and len(dat['GSE255566']) == 6000

    cols = ['UMAP1','UMAP2'] if mode == 'UMAP' else ['TSNE1','TSNE2']
    centers = {ds: ((frame[cols].max() + frame[cols].min()) / 2).to_numpy()
               for ds, frame in dat.items()}
    span = max(float(np.ptp(frame[col])) for frame in dat.values() for col in cols) * 1.12
    frozen_geometry = next(item for item in baseline_manifest['geometry'] if item['algorithm'] == mode)
    assert np.isclose(span, -2*frozen_geometry['coordinate_xlim'][0], rtol=0, atol=1e-12)
    span = -2*frozen_geometry['coordinate_xlim'][0]
    slot = 1.10 * span
    widthmm, heightmm = 180., 86.
    fw, fh = widthmm / 25.4, heightmm / 25.4
    left, right = .060, .985
    axwidth = right - left
    xlim = tuple(frozen_geometry['coordinate_xlim'])
    ylim = tuple(frozen_geometry['coordinate_ylim'])
    axheight = fw * axwidth * np.ptp(ylim) / np.ptp(xlim) / fh
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([left, .29, axwidth, axheight])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(('UMAP' if mode == 'UMAP' else 't-SNE')+' 1', labelpad=3)
    ax.set_ylabel(('UMAP' if mode == 'UMAP' else 't-SNE')+' 2', labelpad=3)
    translations, scatterdata, headings, counts = [], [], [], []
    coordinate_checks = []
    baseline_translations = pd.read_csv(sourcefiles[4])
    baseline_translations = baseline_translations.loc[baseline_translations.algorithm.eq(mode)]
    for index, (ds, group) in enumerate(views):
        frame = dat[ds].loc[dat[ds].group.eq(group)]
        translation = np.array([index * slot - centers[ds][0], -centers[ds][1]])
        xy = frame[cols].to_numpy() + translation
        original_row = baseline_translations.loc[
            baseline_translations.dataset.eq(ds) & baseline_translations.group.eq(group)].iloc[0]
        reference_translation = original_row[['dx', 'dy']].to_numpy(dtype=float)
        assert np.allclose(translation, reference_translation, rtol=0, atol=1e-12)
        # Use the supplied reference translations exactly.
        translation = reference_translation.copy()
        xy = frame[cols].to_numpy() + translation
        reference_xy = frame[cols].to_numpy() + reference_translation
        max_delta = float(np.max(np.abs(xy - reference_xy)))
        assert max_delta < 1e-13
        ix = np.random.default_rng(20261102).permutation(len(frame))
        ax.scatter(xy[ix, 0], xy[ix, 1], s=1.50,
            c=np.asarray([palette[cell_type] for cell_type in frame.cell_type])[ix],
            alpha=.84, linewidths=0, rasterized=True, zorder=1)
        headings.append(ax.text(index * slot, .625 * span, group, ha='center', va='center',
            fontsize=sizes['heading'], fontweight='bold', color=cohort_colors[ds]))
        counts.append(ax.text(index * slot, .535 * span, f'N = {len(frame):,} cells',
            ha='center', va='center', fontsize=sizes['cell_count']))
        translations.append({'algorithm': mode, 'dataset': ds, 'group': group,
            'display_cells': len(frame), 'libraries': int(frame['sample'].nunique()),
            'source_x': cols[0], 'source_y': cols[1], 'dx': float(translation[0]),
            'dy': float(translation[1]), 'scale_x': 1, 'scale_y': 1,
            'rotation_degrees': 0, 'jitter': False})
        coordinate_checks.append({'dataset': ds, 'group': group, 'cell_count': len(frame),
            'max_display_coordinate_delta_from_reference': max_delta,
            'display_translation_matches_reference': True,
            'source_row_order_cell_ids_samples_labels_unchanged': True})
        scatterdata.append((index, ds, group, frame, xy))

    firstmid = left + axwidth * ((.5 * slot - xlim[0]) / np.ptp(xlim))
    secondmid = left + axwidth * ((2.5 * slot - xlim[0]) / np.ptp(xlim))
    cohort_artists = []
    for mid, text, ds in [(firstmid, 'CAD cohort · GSE289422', 'GSE289422'),
                           (secondmid, 'T2D cohort · GSE255566', 'GSE255566')]:
        length = axwidth * .455
        fig.add_artist(Line2D([mid - length / 2, mid + length / 2], [.211, .211],
            transform=fig.transFigure, color=cohort_colors[ds], linewidth=2.2, solid_capstyle='butt'))
        cohort_artists.append(fig.text(mid, .190, text, ha='center', va='top',
            fontsize=sizes['cohort'], fontweight='bold', color=cohort_colors[ds]))
    present = [ct for ct in palette if any(ct in set(frame.cell_type) for frame in dat.values())]
    handles = [Line2D([0], [0], marker='o', color='none', markerfacecolor=palette[ct],
        markeredgewidth=0, markersize=3.6, label=labels[ct]) for ct in present]
    legend = fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.53, .005), ncol=6,
        frameon=False, columnspacing=.70, handletextpad=.28, labelspacing=.45)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pt = fig.dpi / 72
    font = FontProperties(family=font_family, size=sizes['direct_cell_type'], weight='bold')
    # These point-space candidates move annotation text only, never UMAP points.
    offsets = sorted([(dx, dy) for dx in [0, -5, 5, -10, 10, -16, 16, -23, 23]
                      for dy in [0, 7, -7, 14, -14, 21, -21, 28, -28]],
                     key=lambda offset: offset[0] ** 2 + offset[1] ** 2)
    placed, artists, labelledger = [], [], []
    for index, ds, group, frame, xy in scatterdata:
        a0 = ax.transData.transform((index * slot - .48 * span, -.52 * span))
        a1 = ax.transData.transform((index * slot + .48 * span, .475 * span))
        boundary = (a0[0], a0[1], a1[0], a1[1])
        for ct in ['Monocyte', 'CD4_T', 'CD8_T', 'NK', 'B']:
            mask = frame.cell_type.to_numpy() == ct
            if not mask.any():
                continue
            q = xy[mask]
            med = np.median(q, axis=0)
            anchor = q[np.argmin(np.sum((q - med) ** 2, axis=1))]
            text = labels[ct]
            tw, th, descent = renderer.get_text_width_height_descent(text, font, ismath=False)
            halfw = tw / 2 + 1.3 * pt
            halfh = sizes['direct_cell_type'] * pt / 2 + 1.3 * pt
            native = ax.transData.transform(anchor)
            best = None
            for dx, dy in offsets:
                pos = native + np.array([dx, dy]) * pt
                pos[0] = np.clip(pos[0], boundary[0] + halfw, boundary[2] - halfw)
                pos[1] = np.clip(pos[1], boundary[1] + halfh, boundary[3] - halfh)
                rect = (pos[0] - halfw, pos[1] - halfh, pos[0] + halfw, pos[1] + halfh)
                overlap = 0.
                for other in placed:
                    ox = max(0, min(rect[2] + 1.2 * pt, other[2]) - max(rect[0] - 1.2 * pt, other[0]))
                    oy = max(0, min(rect[3] + 1.2 * pt, other[3]) - max(rect[1] - 1.2 * pt, other[1]))
                    overlap += ox * oy
                score = overlap * 10000 + np.sum((pos - native) ** 2)
                if best is None or score < best[0]:
                    best = (score, pos, rect, overlap)
            score, pos, rect, overlap = best
            assert overlap == 0, (mode, ds, group, ct, 'label overlap')
            placed.append(rect)
            display = ax.transData.inverted().transform(pos)
            dist = np.linalg.norm(pos - native) / pt
            arrow = {'arrowstyle': '-', 'color': '#646464', 'lw': .35,
                     'shrinkA': 3, 'shrinkB': 1} if dist > 5 else None
            artist = ax.annotate(text, xy=anchor, xytext=display, ha='center', va='center',
                fontsize=sizes['direct_cell_type'], fontweight='bold', zorder=3,
                arrowprops=arrow, bbox=None)
            artists.append(artist)
            labelledger.append({'algorithm': mode, 'dataset': ds, 'group': group, 'cell_type': ct,
                'label': text, 'cell_type_display_cells': int(mask.sum()),
                'anchor_x': float(anchor[0]), 'anchor_y': float(anchor[1]),
                'label_x': float(display[0]), 'label_y': float(display[1]),
                'label_shift_points': float(dist), 'cell_coordinates_changed': False,
                'font_points': sizes['direct_cell_type'], 'label_rectangle_display_pixels': list(rect)})

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Verify actual text rectangles rather than relying only on candidate metrics.
    from matplotlib.text import Text
    label_boxes = [Text.get_window_extent(artist, renderer) for artist in artists]
    label_collisions = [(i, j) for i, first in enumerate(label_boxes)
                        for j, second in enumerate(label_boxes) if j > i and first.overlaps(second)]
    assert not label_collisions, label_collisions
    all_text = artists + headings + counts + cohort_artists + [ax.xaxis.label, ax.yaxis.label] + legend.get_texts()
    outside = [artist.get_text() for artist in all_text
               if not fig.bbox.contains(*Text.get_window_extent(artist, renderer).get_points()[0])
               or not fig.bbox.contains(*Text.get_window_extent(artist, renderer).get_points()[1])]
    assert not outside, ('Text outside canvas', outside)
    position = ax.get_position()
    physical_ratio = (position.width * fw / (position.height * fh)) / np.ptp(xlim) * np.ptp(ylim)
    assert abs(physical_ratio - 1) < 1e-10
    baseline_geometry = next(item for item in baseline_manifest['geometry'] if item['algorithm'] == mode)
    assert np.allclose(position.bounds, baseline_geometry['axes_normalized_bounds'], rtol=0, atol=1e-14)
    assert np.array_equal(xlim, baseline_geometry['coordinate_xlim'])
    assert np.array_equal(ylim, baseline_geometry['coordinate_ylim'])
    geometry = {'algorithm': mode, 'width_mm': widthmm, 'height_mm': heightmm,
        'coordinate_xlim': list(xlim), 'coordinate_ylim': list(ylim), 'coordinate_aspect': 1,
        'axes_normalized_bounds': [float(value) for value in position.bounds],
        'geometry_matches_reference': True, 'one_outer_coordinate_frame': True,
        'direct_cell_type_label_count': len(artists), 'actual_text_rectangles_overlap': False,
        'text_within_canvas': True, 'annotation_label_backgrounds': False}
    stem = f'Figure5_celltypes_{mode}'
    exports = {}
    for suffix in ['.pdf', '.svg', '_300dpi.png', '_600dpi.png']:
        path = figdir / (stem + suffix)
        if path.exists():
            raise FileExistsError(path)
        fig.savefig(path, dpi=300 if suffix == '_300dpi.png' else 600)
        exports[str(path)] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    plt.close(fig)
    with (meta / 'coordinate_translation_ledger.csv').open('x') as handle:
        pd.DataFrame(translations).to_csv(handle, index=False)
    with (meta / 'direct_cell_type_label_ledger.csv').open('x') as handle:
        pd.DataFrame(labelledger).to_csv(handle, index=False)
    for path, digest in hashes.items():
        assert sha(path) == digest, f'Input changed during render: {path}'
    write_json_new(meta / 'RENDER_MANIFEST.json', {
        'input_sha256': hashes, 'script_sha256': sha(__file__), 'files': exports,
        'font_family': font_family, 'font_scale': scale, 'font_points': sizes,
        'cohort_annotation_colors': cohort_colors,
        'both_control_and_disease_headers_use_cohort_color': True,
        'cell_type_point_palette': palette,
        'cohort_annotation_colours_preserved': True,
        'cohort_band_and_name_vertical_positions': [.211,.190],
        'coordinate_checks': coordinate_checks, 'geometry': geometry,
        'source_coordinates_unchanged': True, 'original_annotations_unchanged': True,
        'source_cell_ids_samples_counts_and_order_unchanged': True,
        'display_translations_match_reference': True,
        'only_annotation_positions_may_shift_for_legibility': True,
        'no_model_refitting_no_UMAP_or_tSNE_computation': True,
        'tSNE_not_rendered': mode != 'tSNE', 'no_new_titles_letters_notes': True,
        'vector_text_editable': True, 'points_rasterized_at_dpi': 600,
        'python_version': sys.version, 'matplotlib_version': matplotlib.__version__,
        'numpy_version': np.__version__, 'pandas_version': pd.__version__})
    print('Saved', stem, flush=True)


if __name__ == '__main__':
    main()
