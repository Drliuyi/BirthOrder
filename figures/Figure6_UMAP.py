"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

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
    for key in ['input-root', 'script-root', 'output-root']:
        env = 'SCRNA_STYLE_' + key.replace('-', '_').upper()
        parser.add_argument('--' + key, default=os.getenv(env), required=not os.getenv(env))
    roots = {key: Path(value).expanduser().resolve() for key, value in vars(parser.parse_args()).items()}
    for key, value in roots.items():
        print(key, value, flush=True)
    inp, code, out = [roots[key] for key in ['input_root', 'script_root', 'output_root']]
    if not (code / Path(__file__).name).is_file():
        raise FileNotFoundError(f'Missing script-root entrypoint: {code / Path(__file__).name}')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Output root must be new or empty: {out}')
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties, fontManager
    from matplotlib.lines import Line2D
    sourcefiles = [inp / f'results/{ds}/coordinates.csv.gz' for ds in ['GSE289422', 'GSE255566']]
    sourcefiles += [inp / 'render_v2/metadata/cell_type_palette.csv', inp / 'render_v2/metadata/INDEPENDENT_FINAL_VERIFICATION.json', inp / 'render_v4_no_label_background/metadata/coordinate_translation_ledger.csv', inp / 'render_v4_no_label_background/metadata/RENDER_MANIFEST.json', inp / 'scripts/draw_primary_single_row_v4_no_background.py']
    for path in sourcefiles:
        if not path.is_file():
            raise FileNotFoundError(f'Missing frozen input or source provenance: {path}')
    hashes = {str(path): sha(path) for path in sourcefiles}
    assert json.loads(sourcefiles[3].read_text())['verified']
    baseline_manifest = json.loads(sourcefiles[5].read_text())
    for path, digest in baseline_manifest['input_sha256'].items():
        assert sha(path) == digest, f'Frozen source changed since original render: {path}'
    for folder in ['figures', 'metadata']:
        (out / folder).mkdir(parents=True, exist_ok=False)
    figdir, meta = (out / 'figures', out / 'metadata')
    os.environ.setdefault('MPLCONFIGDIR', str(meta / 'matplotlib_cache'))
    scale = 1.3
    sizes = {'base': 8 * scale, 'axis': 8 * scale, 'legend': 7 * scale, 'heading': 8.5 * scale, 'cell_count': 7.3 * scale, 'cohort': 7.5 * scale, 'direct_cell_type': 7 * scale}
    cohort_colors = {'GSE289422': '#B65F64', 'GSE255566': '#5683A8'}
    fontfile = Path('/System/Library/Fonts/Helvetica.ttc')
    if fontfile.is_file():
        fontManager.addfont(str(fontfile))
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'], 'font.size': sizes['base'], 'axes.labelsize': sizes['axis'], 'axes.labelweight': 'bold', 'legend.fontsize': sizes['legend'], 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.975, 'axes.edgecolor': '#202020', 'text.color': '#202020', 'axes.labelcolor': '#202020', 'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'svg.hashsalt': 'umap-expression-style-20261105-v1'})
    dat = {ds: pd.read_csv(inp / f'results/{ds}/coordinates.csv.gz') for ds in ['GSE289422', 'GSE255566']}
    pal = pd.read_csv(sourcefiles[2])
    palette = dict(zip(pal.cell_type, pal.colour))
    labels = dict(zip(pal.cell_type, pal.label))
    views = [('GSE289422', 'Control'), ('GSE289422', 'CAD'), ('GSE255566', 'Control'), ('GSE255566', 'T2D')]
    for ds, frame in dat.items():
        assert frame.cell_id.is_unique and set(frame.cell_type).issubset(palette)
        assert not frame.group.eq('Treatment').any()
    assert len(dat['GSE289422']) == 44000 and len(dat['GSE255566']) == 6000
    mode, cols = ('UMAP', ['UMAP1', 'UMAP2'])
    centers = {ds: ((frame[cols].max() + frame[cols].min()) / 2).to_numpy() for ds, frame in dat.items()}
    span = max((float(np.ptp(frame[col])) for frame in dat.values() for col in cols)) * 1.12
    slot = 1.1 * span
    widthmm, heightmm = (180.0, 86.0)
    fw, fh = (widthmm / 25.4, heightmm / 25.4)
    left, right = (0.06, 0.985)
    axwidth = right - left
    xlim = (-0.5 * span, 0.5 * span + 3 * slot)
    ylim = (-0.58 * span, 0.72 * span)
    axheight = fw * axwidth * np.ptp(ylim) / np.ptp(xlim) / fh
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([left, 0.29, axwidth, axheight])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP 1', labelpad=3)
    ax.set_ylabel('UMAP 2', labelpad=3)
    translations, scatterdata, headings, counts = ([], [], [], [])
    coordinate_checks = []
    baseline_translations = pd.read_csv(sourcefiles[4])
    baseline_translations = baseline_translations.loc[baseline_translations.algorithm.eq('UMAP')]
    for index, (ds, group) in enumerate(views):
        frame = dat[ds].loc[dat[ds].group.eq(group)]
        translation = np.array([index * slot - centers[ds][0], -centers[ds][1]])
        xy = frame[cols].to_numpy() + translation
        original_row = baseline_translations.loc[baseline_translations.dataset.eq(ds) & baseline_translations.group.eq(group)].iloc[0]
        reference_translation = original_row[['dx', 'dy']].to_numpy(dtype=float)
        assert np.allclose(translation, reference_translation, rtol=0, atol=1e-13)
        reference_xy = frame[cols].to_numpy() + reference_translation
        max_delta = float(np.max(np.abs(xy - reference_xy)))
        assert max_delta < 1e-13
        ix = np.random.default_rng(20261102).permutation(len(frame))
        ax.scatter(xy[ix, 0], xy[ix, 1], s=1.5, c=np.asarray([palette[cell_type] for cell_type in frame.cell_type])[ix], alpha=0.84, linewidths=0, rasterized=True, zorder=1)
        headings.append(ax.text(index * slot, 0.625 * span, group, ha='center', va='center', fontsize=sizes['heading'], fontweight='bold', color=cohort_colors[ds]))
        counts.append(ax.text(index * slot, 0.535 * span, f'N = {len(frame):,} cells', ha='center', va='center', fontsize=sizes['cell_count']))
        translations.append({'algorithm': mode, 'dataset': ds, 'group': group, 'display_cells': len(frame), 'libraries': int(frame['sample'].nunique()), 'source_x': cols[0], 'source_y': cols[1], 'dx': float(translation[0]), 'dy': float(translation[1]), 'scale_x': 1, 'scale_y': 1, 'rotation_degrees': 0, 'jitter': False})
        coordinate_checks.append({'dataset': ds, 'group': group, 'cell_count': len(frame), 'max_display_coordinate_delta_from_v4': max_delta, 'display_translation_matches_v4': True, 'source_row_order_cell_ids_samples_labels_unchanged': True})
        scatterdata.append((index, ds, group, frame, xy))
    firstmid = left + axwidth * ((0.5 * slot - xlim[0]) / np.ptp(xlim))
    secondmid = left + axwidth * ((2.5 * slot - xlim[0]) / np.ptp(xlim))
    cohort_artists = []
    for mid, text, ds in [(firstmid, 'CAD cohort · GSE289422', 'GSE289422'), (secondmid, 'T2D cohort · GSE255566', 'GSE255566')]:
        length = axwidth * 0.455
        fig.add_artist(Line2D([mid - length / 2, mid + length / 2], [0.181, 0.181], transform=fig.transFigure, color=cohort_colors[ds], linewidth=2.2, solid_capstyle='butt'))
        cohort_artists.append(fig.text(mid, 0.16, text, ha='center', va='top', fontsize=sizes['cohort'], fontweight='bold', color=cohort_colors[ds]))
    present = [ct for ct in palette if any((ct in set(frame.cell_type) for frame in dat.values()))]
    handles = [Line2D([0], [0], marker='o', color='none', markerfacecolor=palette[ct], markeredgewidth=0, markersize=3.6, label=labels[ct]) for ct in present]
    legend = fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.53, 0.005), ncol=6, frameon=False, columnspacing=0.7, handletextpad=0.28, labelspacing=0.45)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pt = fig.dpi / 72
    font = FontProperties(family='Helvetica', size=sizes['direct_cell_type'], weight='bold')
    offsets = sorted([(dx, dy) for dx in [0, -5, 5, -10, 10, -16, 16, -23, 23] for dy in [0, 7, -7, 14, -14, 21, -21, 28, -28]], key=lambda offset: offset[0] ** 2 + offset[1] ** 2)
    placed, artists, labelledger = ([], [], [])
    for index, ds, group, frame, xy in scatterdata:
        a0 = ax.transData.transform((index * slot - 0.48 * span, -0.52 * span))
        a1 = ax.transData.transform((index * slot + 0.48 * span, 0.475 * span))
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
                overlap = 0.0
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
            arrow = {'arrowstyle': '-', 'color': '#646464', 'lw': 0.35, 'shrinkA': 3, 'shrinkB': 1} if dist > 5 else None
            artist = ax.annotate(text, xy=anchor, xytext=display, ha='center', va='center', fontsize=sizes['direct_cell_type'], fontweight='bold', zorder=3, arrowprops=arrow, bbox=None)
            artists.append(artist)
            labelledger.append({'algorithm': mode, 'dataset': ds, 'group': group, 'cell_type': ct, 'label': text, 'cell_type_display_cells': int(mask.sum()), 'anchor_x': float(anchor[0]), 'anchor_y': float(anchor[1]), 'label_x': float(display[0]), 'label_y': float(display[1]), 'label_shift_points': float(dist), 'cell_coordinates_changed': False, 'font_points': sizes['direct_cell_type'], 'label_rectangle_display_pixels': list(rect)})
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    from matplotlib.text import Text
    label_boxes = [Text.get_window_extent(artist, renderer) for artist in artists]
    label_collisions = [(i, j) for i, first in enumerate(label_boxes) for j, second in enumerate(label_boxes) if j > i and first.overlaps(second)]
    assert not label_collisions, label_collisions
    all_text = artists + headings + counts + cohort_artists + [ax.xaxis.label, ax.yaxis.label] + legend.get_texts()
    outside = [artist.get_text() for artist in all_text if not fig.bbox.contains(*Text.get_window_extent(artist, renderer).get_points()[0]) or not fig.bbox.contains(*Text.get_window_extent(artist, renderer).get_points()[1])]
    assert not outside, ('Text outside canvas', outside)
    position = ax.get_position()
    physical_ratio = position.width * fw / (position.height * fh) / np.ptp(xlim) * np.ptp(ylim)
    assert abs(physical_ratio - 1) < 1e-10
    baseline_geometry = next((item for item in baseline_manifest['geometry'] if item['algorithm'] == mode))
    assert np.allclose(position.bounds, baseline_geometry['axes_normalized_bounds'], rtol=0, atol=1e-14)
    assert np.array_equal(xlim, baseline_geometry['coordinate_xlim'])
    assert np.array_equal(ylim, baseline_geometry['coordinate_ylim'])
    geometry = {'algorithm': mode, 'width_mm': widthmm, 'height_mm': heightmm, 'coordinate_xlim': list(xlim), 'coordinate_ylim': list(ylim), 'coordinate_aspect': 1, 'axes_normalized_bounds': [float(value) for value in position.bounds], 'geometry_identical_to_v4': True, 'one_outer_coordinate_frame': True, 'direct_cell_type_label_count': len(artists), 'actual_text_rectangles_overlap': False, 'text_within_canvas': True, 'annotation_label_backgrounds': False}
    stem = 'Figure_primary_UMAP_single_row_expression_style_v1'
    exports = {}
    for suffix in ['.pdf', '.svg', '_600dpi.png']:
        path = figdir / (stem + suffix)
        if path.exists():
            raise FileExistsError(path)
        fig.savefig(path, dpi=600)
        exports[str(path)] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    plt.close(fig)
    with (meta / 'coordinate_translation_ledger.csv').open('x') as handle:
        pd.DataFrame(translations).to_csv(handle, index=False)
    with (meta / 'direct_cell_type_label_ledger.csv').open('x') as handle:
        pd.DataFrame(labelledger).to_csv(handle, index=False)
    for path, digest in hashes.items():
        assert sha(path) == digest, f'Input changed during render: {path}'
    write_json_new(meta / 'RENDER_MANIFEST.json', {'input_sha256': hashes, 'script_sha256': sha(__file__), 'files': exports, 'font_family': 'Helvetica', 'font_scale_from_v4': scale, 'font_points': sizes, 'cohort_annotation_colors': cohort_colors, 'both_control_and_disease_headers_use_cohort_color': True, 'cell_type_point_palette_unchanged': palette, 'coordinate_checks': coordinate_checks, 'geometry': geometry, 'source_coordinates_unchanged': True, 'original_annotations_unchanged': True, 'source_cell_ids_samples_counts_and_order_unchanged': True, 'display_translations_identical_to_v4': True, 'only_annotation_positions_may_shift_for_legibility': True, 'no_model_refitting_no_UMAP_or_tSNE_computation': True, 'tSNE_not_rendered': True, 'no_new_titles_letters_notes': True, 'vector_text_editable': True, 'points_rasterized_at_dpi': 600, 'python_version': sys.version, 'matplotlib_version': matplotlib.__version__, 'numpy_version': np.__version__, 'pandas_version': pd.__version__})
    print('Saved', stem, flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Restyle the frozen four-group UMAP; no embedding or data recomputation.\n\nThe four source maps, display translations, point palette, cell IDs, sample\nlabels, counts, axis geometry, and dimensions are retained. Only type size,\ncohort annotation colors, and collision-aware direct label positions change.\n"""\nimport argparse\nimport hashlib\nimport json\nimport os\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\n\ndef sha(path):\n    h = hashlib.sha256()\n    with Path(path).open(\'rb\') as handle:\n        for block in iter(lambda: handle.read(1048576), b\'\'):\n            h.update(block)\n    return h.hexdigest()\n\ndef write_json_new(path, value):\n    with Path(path).open(\'x\') as handle:\n        json.dump(value, handle, indent=2)\n        handle.write(\'\\n\')\n\ndef main():\n    parser = argparse.ArgumentParser(description=__doc__)\n    for key in [\'input-root\', \'script-root\', \'output-root\']:\n        env = \'SCRNA_STYLE_\' + key.replace(\'-\', \'_\').upper()\n        parser.add_argument(\'--\' + key, default=os.getenv(env), required=not os.getenv(env))\n    roots = {key: Path(value).expanduser().resolve() for key, value in vars(parser.parse_args()).items()}\n    for key, value in roots.items():\n        print(key, value, flush=True)\n    inp, code, out = [roots[key] for key in [\'input_root\', \'script_root\', \'output_root\']]\n    if not (code / Path(__file__).name).is_file():\n        raise FileNotFoundError(f\'Missing script-root entrypoint: {code / Path(__file__).name}\')\n    if out.exists() and any(out.iterdir()):\n        raise FileExistsError(f\'Output root must be new or empty: {out}\')\n    import numpy as np\n    import pandas as pd\n    import matplotlib\n    matplotlib.use(\'Agg\')\n    import matplotlib.pyplot as plt\n    from matplotlib.font_manager import FontProperties, fontManager\n    from matplotlib.lines import Line2D\n    sourcefiles = [inp / f\'results/{ds}/coordinates.csv.gz\' for ds in [\'GSE289422\', \'GSE255566\']]\n    sourcefiles += [inp / \'render_v2/metadata/cell_type_palette.csv\', inp / \'render_v2/metadata/INDEPENDENT_FINAL_VERIFICATION.json\', inp / \'render_v4_no_label_background/metadata/coordinate_translation_ledger.csv\', inp / \'render_v4_no_label_background/metadata/RENDER_MANIFEST.json\', inp / \'scripts/draw_primary_single_row_v4_no_background.py\']\n    for path in sourcefiles:\n        if not path.is_file():\n            raise FileNotFoundError(f\'Missing frozen input or source provenance: {path}\')\n    hashes = {str(path): sha(path) for path in sourcefiles}\n    assert json.loads(sourcefiles[3].read_text())[\'verified\']\n    baseline_manifest = json.loads(sourcefiles[5].read_text())\n    for path, digest in baseline_manifest[\'input_sha256\'].items():\n        assert sha(path) == digest, f\'Frozen source changed since original render: {path}\'\n    for folder in [\'figures\', \'metadata\']:\n        (out / folder).mkdir(parents=True, exist_ok=False)\n    figdir, meta = (out / \'figures\', out / \'metadata\')\n    os.environ.setdefault(\'MPLCONFIGDIR\', str(meta / \'matplotlib_cache\'))\n    scale = 1.3\n    sizes = {\'base\': 8 * scale, \'axis\': 8 * scale, \'legend\': 7 * scale, \'heading\': 8.5 * scale, \'cell_count\': 7.3 * scale, \'cohort\': 7.5 * scale, \'direct_cell_type\': 7 * scale}\n    cohort_colors = {\'GSE289422\': \'#B65F64\', \'GSE255566\': \'#5683A8\'}\n    fontfile = Path(\'/System/Library/Fonts/Helvetica.ttc\')\n    if fontfile.is_file():\n        fontManager.addfont(str(fontfile))\n    plt.rcParams.update({\'font.family\': \'sans-serif\', \'font.sans-serif\': [\'Helvetica\', \'Arial\', \'DejaVu Sans\'], \'font.size\': sizes[\'base\'], \'axes.labelsize\': sizes[\'axis\'], \'axes.labelweight\': \'bold\', \'legend.fontsize\': sizes[\'legend\'], \'pdf.fonttype\': 42, \'ps.fonttype\': 42, \'svg.fonttype\': \'none\', \'axes.linewidth\': 0.975, \'axes.edgecolor\': \'#202020\', \'text.color\': \'#202020\', \'axes.labelcolor\': \'#202020\', \'figure.facecolor\': \'white\', \'savefig.facecolor\': \'white\', \'svg.hashsalt\': \'umap-expression-style-20261105-v1\'})\n    dat = {ds: pd.read_csv(inp / f\'results/{ds}/coordinates.csv.gz\') for ds in [\'GSE289422\', \'GSE255566\']}\n    pal = pd.read_csv(sourcefiles[2])\n    palette = dict(zip(pal.cell_type, pal.colour))\n    labels = dict(zip(pal.cell_type, pal.label))\n    views = [(\'GSE289422\', \'Control\'), (\'GSE289422\', \'CAD\'), (\'GSE255566\', \'Control\'), (\'GSE255566\', \'T2D\')]\n    for ds, frame in dat.items():\n        assert frame.cell_id.is_unique and set(frame.cell_type).issubset(palette)\n        assert not frame.group.eq(\'Treatment\').any()\n    assert len(dat[\'GSE289422\']) == 44000 and len(dat[\'GSE255566\']) == 6000\n    mode, cols = (\'UMAP\', [\'UMAP1\', \'UMAP2\'])\n    centers = {ds: ((frame[cols].max() + frame[cols].min()) / 2).to_numpy() for ds, frame in dat.items()}\n    span = max((float(np.ptp(frame[col])) for frame in dat.values() for col in cols)) * 1.12\n    slot = 1.1 * span\n    widthmm, heightmm = (180.0, 86.0)\n    fw, fh = (widthmm / 25.4, heightmm / 25.4)\n    left, right = (0.06, 0.985)\n    axwidth = right - left\n    xlim = (-0.5 * span, 0.5 * span + 3 * slot)\n    ylim = (-0.58 * span, 0.72 * span)\n    axheight = fw * axwidth * np.ptp(ylim) / np.ptp(xlim) / fh\n    fig = plt.figure(figsize=(fw, fh))\n    ax = fig.add_axes([left, 0.29, axwidth, axheight])\n    ax.set_xlim(xlim)\n    ax.set_ylim(ylim)\n    ax.set_aspect(\'equal\', adjustable=\'box\')\n    ax.set_xticks([])\n    ax.set_yticks([])\n    ax.set_xlabel(\'UMAP 1\', labelpad=3)\n    ax.set_ylabel(\'UMAP 2\', labelpad=3)\n    translations, scatterdata, headings, counts = ([], [], [], [])\n    coordinate_checks = []\n    baseline_translations = pd.read_csv(sourcefiles[4])\n    baseline_translations = baseline_translations.loc[baseline_translations.algorithm.eq(\'UMAP\')]\n    for index, (ds, group) in enumerate(views):\n        frame = dat[ds].loc[dat[ds].group.eq(group)]\n        translation = np.array([index * slot - centers[ds][0], -centers[ds][1]])\n        xy = frame[cols].to_numpy() + translation\n        original_row = baseline_translations.loc[baseline_translations.dataset.eq(ds) & baseline_translations.group.eq(group)].iloc[0]\n        reference_translation = original_row[[\'dx\', \'dy\']].to_numpy(dtype=float)\n        assert np.allclose(translation, reference_translation, rtol=0, atol=1e-13)\n        reference_xy = frame[cols].to_numpy() + reference_translation\n        max_delta = float(np.max(np.abs(xy - reference_xy)))\n        assert max_delta < 1e-13\n        ix = np.random.default_rng(20261102).permutation(len(frame))\n        ax.scatter(xy[ix, 0], xy[ix, 1], s=1.5, c=np.asarray([palette[cell_type] for cell_type in frame.cell_type])[ix], alpha=0.84, linewidths=0, rasterized=True, zorder=1)\n        headings.append(ax.text(index * slot, 0.625 * span, group, ha=\'center\', va=\'center\', fontsize=sizes[\'heading\'], fontweight=\'bold\', color=cohort_colors[ds]))\n        counts.append(ax.text(index * slot, 0.535 * span, f\'N = {len(frame):,} cells\', ha=\'center\', va=\'center\', fontsize=sizes[\'cell_count\']))\n        translations.append({\'algorithm\': mode, \'dataset\': ds, \'group\': group, \'display_cells\': len(frame), \'libraries\': int(frame[\'sample\'].nunique()), \'source_x\': cols[0], \'source_y\': cols[1], \'dx\': float(translation[0]), \'dy\': float(translation[1]), \'scale_x\': 1, \'scale_y\': 1, \'rotation_degrees\': 0, \'jitter\': False})\n        coordinate_checks.append({\'dataset\': ds, \'group\': group, \'cell_count\': len(frame), \'max_display_coordinate_delta_from_v4\': max_delta, \'display_translation_matches_v4\': True, \'source_row_order_cell_ids_samples_labels_unchanged\': True})\n        scatterdata.append((index, ds, group, frame, xy))\n    firstmid = left + axwidth * ((0.5 * slot - xlim[0]) / np.ptp(xlim))\n    secondmid = left + axwidth * ((2.5 * slot - xlim[0]) / np.ptp(xlim))\n    cohort_artists = []\n    for mid, text, ds in [(firstmid, \'CAD cohort · GSE289422\', \'GSE289422\'), (secondmid, \'T2D cohort · GSE255566\', \'GSE255566\')]:\n        length = axwidth * 0.455\n        fig.add_artist(Line2D([mid - length / 2, mid + length / 2], [0.181, 0.181], transform=fig.transFigure, color=cohort_colors[ds], linewidth=2.2, solid_capstyle=\'butt\'))\n        cohort_artists.append(fig.text(mid, 0.16, text, ha=\'center\', va=\'top\', fontsize=sizes[\'cohort\'], fontweight=\'bold\', color=cohort_colors[ds]))\n    present = [ct for ct in palette if any((ct in set(frame.cell_type) for frame in dat.values()))]\n    handles = [Line2D([0], [0], marker=\'o\', color=\'none\', markerfacecolor=palette[ct], markeredgewidth=0, markersize=3.6, label=labels[ct]) for ct in present]\n    legend = fig.legend(handles=handles, loc=\'lower center\', bbox_to_anchor=(0.53, 0.005), ncol=6, frameon=False, columnspacing=0.7, handletextpad=0.28, labelspacing=0.45)\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    pt = fig.dpi / 72\n    font = FontProperties(family=\'Helvetica\', size=sizes[\'direct_cell_type\'], weight=\'bold\')\n    offsets = sorted([(dx, dy) for dx in [0, -5, 5, -10, 10, -16, 16, -23, 23] for dy in [0, 7, -7, 14, -14, 21, -21, 28, -28]], key=lambda offset: offset[0] ** 2 + offset[1] ** 2)\n    placed, artists, labelledger = ([], [], [])\n    for index, ds, group, frame, xy in scatterdata:\n        a0 = ax.transData.transform((index * slot - 0.48 * span, -0.52 * span))\n        a1 = ax.transData.transform((index * slot + 0.48 * span, 0.475 * span))\n        boundary = (a0[0], a0[1], a1[0], a1[1])\n        for ct in [\'Monocyte\', \'CD4_T\', \'CD8_T\', \'NK\', \'B\']:\n            mask = frame.cell_type.to_numpy() == ct\n            if not mask.any():\n                continue\n            q = xy[mask]\n            med = np.median(q, axis=0)\n            anchor = q[np.argmin(np.sum((q - med) ** 2, axis=1))]\n            text = labels[ct]\n            tw, th, descent = renderer.get_text_width_height_descent(text, font, ismath=False)\n            halfw = tw / 2 + 1.3 * pt\n            halfh = sizes[\'direct_cell_type\'] * pt / 2 + 1.3 * pt\n            native = ax.transData.transform(anchor)\n            best = None\n            for dx, dy in offsets:\n                pos = native + np.array([dx, dy]) * pt\n                pos[0] = np.clip(pos[0], boundary[0] + halfw, boundary[2] - halfw)\n                pos[1] = np.clip(pos[1], boundary[1] + halfh, boundary[3] - halfh)\n                rect = (pos[0] - halfw, pos[1] - halfh, pos[0] + halfw, pos[1] + halfh)\n                overlap = 0.0\n                for other in placed:\n                    ox = max(0, min(rect[2] + 1.2 * pt, other[2]) - max(rect[0] - 1.2 * pt, other[0]))\n                    oy = max(0, min(rect[3] + 1.2 * pt, other[3]) - max(rect[1] - 1.2 * pt, other[1]))\n                    overlap += ox * oy\n                score = overlap * 10000 + np.sum((pos - native) ** 2)\n                if best is None or score < best[0]:\n                    best = (score, pos, rect, overlap)\n            score, pos, rect, overlap = best\n            assert overlap == 0, (mode, ds, group, ct, \'label overlap\')\n            placed.append(rect)\n            display = ax.transData.inverted().transform(pos)\n            dist = np.linalg.norm(pos - native) / pt\n            arrow = {\'arrowstyle\': \'-\', \'color\': \'#646464\', \'lw\': 0.35, \'shrinkA\': 3, \'shrinkB\': 1} if dist > 5 else None\n            artist = ax.annotate(text, xy=anchor, xytext=display, ha=\'center\', va=\'center\', fontsize=sizes[\'direct_cell_type\'], fontweight=\'bold\', zorder=3, arrowprops=arrow, bbox=None)\n            artists.append(artist)\n            labelledger.append({\'algorithm\': mode, \'dataset\': ds, \'group\': group, \'cell_type\': ct, \'label\': text, \'cell_type_display_cells\': int(mask.sum()), \'anchor_x\': float(anchor[0]), \'anchor_y\': float(anchor[1]), \'label_x\': float(display[0]), \'label_y\': float(display[1]), \'label_shift_points\': float(dist), \'cell_coordinates_changed\': False, \'font_points\': sizes[\'direct_cell_type\'], \'label_rectangle_display_pixels\': list(rect)})\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    from matplotlib.text import Text\n    label_boxes = [Text.get_window_extent(artist, renderer) for artist in artists]\n    label_collisions = [(i, j) for i, first in enumerate(label_boxes) for j, second in enumerate(label_boxes) if j > i and first.overlaps(second)]\n    assert not label_collisions, label_collisions\n    all_text = artists + headings + counts + cohort_artists + [ax.xaxis.label, ax.yaxis.label] + legend.get_texts()\n    outside = [artist.get_text() for artist in all_text if not fig.bbox.contains(*Text.get_window_extent(artist, renderer).get_points()[0]) or not fig.bbox.contains(*Text.get_window_extent(artist, renderer).get_points()[1])]\n    assert not outside, (\'Text outside canvas\', outside)\n    position = ax.get_position()\n    physical_ratio = position.width * fw / (position.height * fh) / np.ptp(xlim) * np.ptp(ylim)\n    assert abs(physical_ratio - 1) < 1e-10\n    baseline_geometry = next((item for item in baseline_manifest[\'geometry\'] if item[\'algorithm\'] == mode))\n    assert np.allclose(position.bounds, baseline_geometry[\'axes_normalized_bounds\'], rtol=0, atol=1e-14)\n    assert np.array_equal(xlim, baseline_geometry[\'coordinate_xlim\'])\n    assert np.array_equal(ylim, baseline_geometry[\'coordinate_ylim\'])\n    geometry = {\'algorithm\': mode, \'width_mm\': widthmm, \'height_mm\': heightmm, \'coordinate_xlim\': list(xlim), \'coordinate_ylim\': list(ylim), \'coordinate_aspect\': 1, \'axes_normalized_bounds\': [float(value) for value in position.bounds], \'geometry_identical_to_v4\': True, \'one_outer_coordinate_frame\': True, \'direct_cell_type_label_count\': len(artists), \'actual_text_rectangles_overlap\': False, \'text_within_canvas\': True, \'annotation_label_backgrounds\': False}\n    stem = \'Figure_primary_UMAP_single_row_expression_style_v1\'\n    exports = {}\n    for suffix in [\'.pdf\', \'.svg\', \'_600dpi.png\']:\n        path = figdir / (stem + suffix)\n        if path.exists():\n            raise FileExistsError(path)\n        fig.savefig(path, dpi=600)\n        exports[str(path)] = {\'sha256\': sha(path), \'bytes\': path.stat().st_size}\n    plt.close(fig)\n    with (meta / \'coordinate_translation_ledger.csv\').open(\'x\') as handle:\n        pd.DataFrame(translations).to_csv(handle, index=False)\n    with (meta / \'direct_cell_type_label_ledger.csv\').open(\'x\') as handle:\n        pd.DataFrame(labelledger).to_csv(handle, index=False)\n    for path, digest in hashes.items():\n        assert sha(path) == digest, f\'Input changed during render: {path}\'\n    write_json_new(meta / \'RENDER_MANIFEST.json\', {\'input_sha256\': hashes, \'script_sha256\': sha(__file__), \'files\': exports, \'font_family\': \'Helvetica\', \'font_scale_from_v4\': scale, \'font_points\': sizes, \'cohort_annotation_colors\': cohort_colors, \'both_control_and_disease_headers_use_cohort_color\': True, \'cell_type_point_palette_unchanged\': palette, \'coordinate_checks\': coordinate_checks, \'geometry\': geometry, \'source_coordinates_unchanged\': True, \'original_annotations_unchanged\': True, \'source_cell_ids_samples_counts_and_order_unchanged\': True, \'display_translations_identical_to_v4\': True, \'only_annotation_positions_may_shift_for_legibility\': True, \'no_model_refitting_no_UMAP_or_tSNE_computation\': True, \'tSNE_not_rendered\': True, \'no_new_titles_letters_notes\': True, \'vector_text_editable\': True, \'points_rasterized_at_dpi\': 600, \'python_version\': sys.version, \'matplotlib_version\': matplotlib.__version__, \'numpy_version\': np.__version__, \'pandas_version\': pd.__version__})\n    print(\'Saved\', stem, flush=True)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
