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
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2)
        f.write('\n')

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key, env in [('input-root', 'SCRNA_PUBLICATION_INPUT_ROOT'), ('reference-root', 'SCRNA_PUBLICATION_REFERENCE_ROOT'), ('script-root', 'SCRIPT_ROOT'), ('output-root', 'SCRNA_PUBLICATION_FIGURE_ROOT')]:
        p.add_argument('--' + key, default=os.getenv(env), required=not os.getenv(env))
    args = p.parse_args()
    inp, ref, src, out = [Path(getattr(args, k)).expanduser().resolve() for k in ('input_root', 'reference_root', 'script_root', 'output_root')]
    for key, path in [('input_root', inp), ('reference_root', ref), ('script_root', src), ('output_root', out)]:
        print(key, path, flush=True)
    reference_files = [inp / 'metadata/PREPARATION_COMPLETE.json', ref / 'metadata/RENDER_MANIFEST.json', ref / 'metadata/feature_colour_limits.csv', ref / 'metadata/coordinate_translation_ledger.csv', ref / 'scripts/draw_features_unlettered_v5.py', src / Path(__file__).name]
    for path in reference_files:
        if not path.is_file():
            raise FileNotFoundError('Missing path: ' + str(path))
    if (src / Path(__file__).name).resolve() != Path(__file__).resolve():
        raise ValueError('Resolved script root does not contain this entrypoint')
    figdir, meta = (out / 'figures', out / 'metadata')
    figdir.mkdir(exist_ok=True, parents=True)
    meta.mkdir(exist_ok=True, parents=True)
    if any(figdir.iterdir()) or any(meta.iterdir()):
        raise FileExistsError('Refusing to overwrite existing figure or metadata outputs')
    os.environ.setdefault('MPLCONFIGDIR', str(meta / 'matplotlib_cache'))
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.font_manager import fontManager
    original = json.loads((ref / 'metadata/RENDER_MANIFEST.json').read_text())
    original_source = ref / 'scripts/draw_features_unlettered_v5.py'
    assert sha(original_source) == original['script_sha256'], 'Original source changed'
    fontfile = Path('/System/Library/Fonts/Helvetica.ttc')
    if fontfile.is_file():
        fontManager.addfont(str(fontfile))
    font_sizes = {'gene_heading': 11.05, 'axis_label': 9.1, 'colorbar_tick': 9.1, 'cohort_label': 10.14, 'existing_expression_note': 9.36}
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'], 'font.size': 10.4, 'axes.labelsize': 10.4, 'axes.titlesize': 11.05, 'xtick.labelsize': 9.1, 'ytick.labelsize': 9.1, 'legend.fontsize': 9.49, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.975, 'axes.edgecolor': '#202020', 'text.color': '#202020', 'axes.labelcolor': '#202020', 'savefig.facecolor': 'white', 'figure.facecolor': 'white', 'svg.hashsalt': 'scrna-style-20261105-v1'})
    dslist = ['GSE289422', 'GSE255566']
    cohort = {'GSE289422': 'CAD cohort', 'GSE255566': 'T2D cohort'}
    disease = {'GSE289422': 'CAD', 'GSE255566': 'T2D'}
    cohortcolor = {'GSE289422': '#B65F64', 'GSE255566': '#5683A8'}
    genes = ['FCRL6', 'GZMA', 'GZMH', 'KLRD1']
    dat, sourcehash, featurehash = ({}, {}, {})
    for ds in dslist:
        path = inp / f'results/{ds}/coordinates.csv.gz'
        if not path.is_file():
            raise FileNotFoundError(path)
        sourcehash[str(path)] = sha(path)
        assert sourcehash[str(path)] == original['coordinate_input_sha256'][str(path)], 'Frozen input changed'
        d = pd.read_csv(path)
        assert d.cell_id.is_unique and d.groupby('sample').size().le(1000).all()
        assert set(d.group) == {'Control', disease[ds]}
        dat[ds] = d
        featurehash[ds] = {'cell_id_order': hashlib.sha256('\n'.join(d.cell_id.astype(str)).encode()).hexdigest(), 'UMAP_coordinates_float64': hashlib.sha256(d[['UMAP1', 'UMAP2']].to_numpy(dtype='<f8').tobytes()).hexdigest(), 'expression_float64': {gene: hashlib.sha256(d[gene + '_logCP10k'].to_numpy(dtype='<f8').tobytes()).hexdigest() for gene in genes}}
    limits = pd.read_csv(ref / 'metadata/feature_colour_limits.csv')
    caps = {}
    for ds in dslist:
        for gene in genes:
            row = limits.loc[(limits.dataset == ds) & (limits.gene == gene)]
            assert len(row) == 1
            values = dat[ds][gene + '_logCP10k'].to_numpy()
            positive = values[values > 0]
            cap = float(np.quantile(positive, 0.99)) if len(positive) else 1.0
            assert float(row.iloc[0].vmin) == 0
            assert np.isclose(float(row.iloc[0].vmax), cap, rtol=1e-15, atol=0)
            caps[ds, gene] = cap
    cmap_stops = [0, 0.15, 0.45, 0.72, 1]
    cmap_colors = ['#ECEDEF', '#F1D8BD', '#D78D67', '#C55432', '#7651A6']
    cmap = mcolors.LinearSegmentedColormap.from_list('YY_gene', list(zip(cmap_stops, cmap_colors)), N=256)
    cols = ['UMAP1', 'UMAP2']
    centers = {ds: ((dat[ds][cols].max() + dat[ds][cols].min()) / 2).to_numpy() for ds in dslist}
    span = max((float(np.ptp(dat[ds][col])) for ds in dslist for col in cols)) * 1.12
    local = {ds: dat[ds][cols].to_numpy() - centers[ds] for ds in dslist}
    ledger, axes_records, text_records = ([], [], [])
    fw, fh = (7.8, 5.12)
    fig = plt.figure(figsize=(fw, fh))
    left, pw, gx = (0.075, 0.204, 0.032)
    ph = pw * fw / fh
    rowbottom = [0.573, 0.128]
    feature_axes, bar_axes = ([], [])
    for row, ds in enumerate(dslist):
        d, xy = (dat[ds], local[ds])
        text_records.append(fig.text(0.022, rowbottom[row] + ph / 2, cohort[ds] + '\n' + ds, fontsize=font_sizes['cohort_label'], fontweight='bold', color=cohortcolor[ds], rotation=90, ha='center', va='center'))
        for col, gene in enumerate(genes):
            x0, y0 = (left + col * (pw + gx), rowbottom[row])
            ax = fig.add_axes([x0, y0, pw, ph])
            feature_axes.append((row, col, ax))
            ax.set_aspect('equal', adjustable='box')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_xlim(-span / 2, span / 2)
            ax.set_ylim(-span / 2, span / 2)
            values = d[gene + '_logCP10k'].to_numpy()
            ix = np.argsort(values, kind='stable')
            sc = ax.scatter(xy[ix, 0], xy[ix, 1], c=values[ix], s=0.68, alpha=0.9, linewidths=0, rasterized=True, cmap=cmap, vmin=0, vmax=caps[ds, gene], zorder=1)
            assert np.array_equal(sc.get_offsets(), xy[ix])
            assert np.array_equal(sc.get_array(), values[ix])
            ax.set_title(gene, fontsize=font_sizes['gene_heading'], fontstyle='italic', fontweight='bold', pad=4)
            ax.set_xlabel('UMAP 1', fontsize=font_sizes['axis_label'], fontweight='bold', labelpad=0)
            ax.set_ylabel('UMAP 2' if col == 0 else '', fontsize=font_sizes['axis_label'], fontweight='bold', labelpad=0)
            cax = fig.add_axes([x0 + 0.018, y0 - 0.047, pw - 0.036, 0.012])
            bar_axes.append(cax)
            cb = fig.colorbar(sc, cax=cax, orientation='horizontal', ticks=[0, caps[ds, gene] / 2, caps[ds, gene]])
            cb.outline.set_linewidth(0.35)
            cb.ax.tick_params(labelsize=font_sizes['colorbar_tick'], length=1.5, width=0.35, pad=1.5)
            cb.set_ticklabels([f'{v:.2g}' for v in [0, caps[ds, gene] / 2, caps[ds, gene]]])
            ledger.append({'algorithm': 'UMAP', 'view': f'FeaturePlot_{gene}', 'dataset': ds, 'group': 'all', 'source_x': cols[0], 'source_y': cols[1], 'dx': float(-centers[ds][0]), 'dy': float(-centers[ds][1]), 'scale_x': 1, 'scale_y': 1, 'rotation_degrees': 0, 'display_cells': len(d), 'slot': col, 'slot_span': span, 'gap': 0})
    text_records.append(fig.text(0.53, 0.02, 'Expression: log1p(CP10k); each cohort–gene scale is independent', ha='center', fontsize=font_sizes['existing_expression_note']))
    fig.canvas.draw()
    for row in [0, 1]:
        boxes = [ax.get_position() for r, c, ax in feature_axes if r == row]
        assert np.allclose([b.y0 for b in boxes], boxes[0].y0)
        assert np.allclose([b.height for b in boxes], boxes[0].height)
    for col in range(4):
        boxes = [ax.get_position() for r, c, ax in feature_axes if c == col]
        assert np.allclose([b.x0 for b in boxes], boxes[0].x0)
        assert np.allclose([b.width for b in boxes], boxes[0].width)
    reference_ledger = pd.read_csv(ref / 'metadata/coordinate_translation_ledger.csv')
    reference_ledger = reference_ledger.loc[reference_ledger.algorithm == 'UMAP'].reset_index(drop=True)
    pd.testing.assert_frame_equal(pd.DataFrame(ledger), reference_ledger, check_dtype=False, check_exact=False, rtol=1e-14)
    renderer = fig.canvas.get_renderer()
    all_text = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    for item in all_text:
        b = item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        assert min(b.x0, b.y0) >= 0 and max(b.x1, b.y1) <= 1, 'Text clipped: ' + item.get_text()
        text_records.append({'text': item.get_text(), 'font_size_pt': item.get_fontsize(), 'bounds_normalized': [b.x0, b.y0, b.width, b.height]})
    text_records = [r for r in text_records if isinstance(r, dict)]
    for ax in fig.axes:
        b = ax.get_position()
        axes_records.append([b.x0, b.y0, b.width, b.height])
    files = []
    for ext, dpi in [('.pdf', 600), ('.svg', 600), ('_600dpi.png', 600)]:
        path = figdir / ('Figure_features_UMAP' + ext)
        if path.exists():
            raise FileExistsError(path)
        fig.savefig(path, dpi=dpi)
        files.append(path)
        print('Saved', path, flush=True)
    plt.close(fig)
    for path, digest in sourcehash.items():
        assert sha(path) == digest
    with (meta / 'feature_colour_limits.csv').open('x') as f:
        f.write((ref / 'metadata/feature_colour_limits.csv').read_text())
    with (meta / 'coordinate_translation_ledger.csv').open('x') as f:
        pd.DataFrame(ledger).to_csv(f, index=False)
    write_json(meta / 'RENDER_MANIFEST.json', {'source_figure': str(ref / 'figures/Figure_features_UMAP_600dpi.png'), 'source_figure_sha256': sha(ref / 'figures/Figure_features_UMAP_600dpi.png'), 'source_entrypoint': str(original_source), 'source_entrypoint_sha256': sha(original_source), 'script': str(Path(__file__).resolve()), 'script_sha256': sha(__file__), 'coordinate_input_sha256': sourcehash, 'data_array_sha256': featurehash, 'display_cells': {ds: len(d) for ds, d in dat.items()}, 'fixed_genes': genes, 'font': 'Helvetica', 'font_sizes_pt': font_sizes, 'font_increase_percent': 30, 'cohort_label_colors': cohortcolor, 'expression_colormap': {'name': 'YY_gene', 'positions': cmap_stops, 'colors': cmap_colors, 'N': 256}, 'expression_colormap_unchanged': True, 'expression_limits_unchanged': True, 'expression_values_unchanged': True, 'coordinate_values_unchanged': True, 'only_original_coordinate_translation': True, 'point_size': 0.68, 'point_alpha': 0.9, 'stable_expression_point_order_unchanged': True, 'no_new_titles_letters_or_notes': True, 'raw_h5ad_not_read': True, 'no_fitting_or_new_embeddings': True, 'figure_size_inches': [fw, fh], 'figure_width_mm': fw * 25.4, 'layout_change': 'All text +30%; cohort colors fixed; row/column whitespace adjusted only', 'equal_square_panels': True, 'row_column_edges_aligned': True, 'axes_normalized_bounds': axes_records, 'text_audit': text_records, 'no_clipped_text': True, 'raster_dpi': 600, 'editable_vector_text': True, 'software_versions': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'matplotlib': matplotlib.__version__}, 'files': {str(f): {'sha256': sha(f), 'bytes': f.stat().st_size} for f in files}, 'source_inputs_unchanged': True})

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Typography-only restyle of frozen core4 UMAP FeaturePlots; no new analysis."""\nimport argparse\nimport hashlib\nimport json\nimport os\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\n\ndef sha(path):\n    h = hashlib.sha256()\n    with Path(path).open(\'rb\') as f:\n        for block in iter(lambda: f.read(1048576), b\'\'):\n            h.update(block)\n    return h.hexdigest()\n\ndef write_json(path, obj):\n    with Path(path).open(\'x\') as f:\n        json.dump(obj, f, indent=2)\n        f.write(\'\\n\')\n\ndef main():\n    p = argparse.ArgumentParser(description=__doc__)\n    for key, env in [(\'input-root\', \'SCRNA_PUBLICATION_INPUT_ROOT\'), (\'reference-root\', \'SCRNA_PUBLICATION_REFERENCE_ROOT\'), (\'script-root\', \'SCRIPT_ROOT\'), (\'output-root\', \'SCRNA_PUBLICATION_FIGURE_ROOT\')]:\n        p.add_argument(\'--\' + key, default=os.getenv(env), required=not os.getenv(env))\n    args = p.parse_args()\n    inp, ref, src, out = [Path(getattr(args, k)).expanduser().resolve() for k in (\'input_root\', \'reference_root\', \'script_root\', \'output_root\')]\n    for key, path in [(\'input_root\', inp), (\'reference_root\', ref), (\'script_root\', src), (\'output_root\', out)]:\n        print(key, path, flush=True)\n    reference_files = [inp / \'metadata/PREPARATION_COMPLETE.json\', ref / \'metadata/RENDER_MANIFEST.json\', ref / \'metadata/feature_colour_limits.csv\', ref / \'metadata/coordinate_translation_ledger.csv\', ref / \'scripts/draw_features_unlettered_v5.py\', src / Path(__file__).name]\n    for path in reference_files:\n        if not path.is_file():\n            raise FileNotFoundError(\'Missing path: \' + str(path))\n    if (src / Path(__file__).name).resolve() != Path(__file__).resolve():\n        raise ValueError(\'Resolved script root does not contain this entrypoint\')\n    figdir, meta = (out / \'figures\', out / \'metadata\')\n    figdir.mkdir(exist_ok=True, parents=True)\n    meta.mkdir(exist_ok=True, parents=True)\n    if any(figdir.iterdir()) or any(meta.iterdir()):\n        raise FileExistsError(\'Refusing to overwrite existing figure or metadata outputs\')\n    os.environ.setdefault(\'MPLCONFIGDIR\', str(meta / \'matplotlib_cache\'))\n    import numpy as np\n    import pandas as pd\n    import matplotlib\n    matplotlib.use(\'Agg\')\n    import matplotlib.pyplot as plt\n    import matplotlib.colors as mcolors\n    from matplotlib.font_manager import fontManager\n    original = json.loads((ref / \'metadata/RENDER_MANIFEST.json\').read_text())\n    original_source = ref / \'scripts/draw_features_unlettered_v5.py\'\n    assert sha(original_source) == original[\'script_sha256\'], \'Original source changed\'\n    fontfile = Path(\'/System/Library/Fonts/Helvetica.ttc\')\n    if fontfile.is_file():\n        fontManager.addfont(str(fontfile))\n    font_sizes = {\'gene_heading\': 11.05, \'axis_label\': 9.1, \'colorbar_tick\': 9.1, \'cohort_label\': 10.14, \'existing_expression_note\': 9.36}\n    plt.rcParams.update({\'font.family\': \'sans-serif\', \'font.sans-serif\': [\'Helvetica\', \'Arial\', \'DejaVu Sans\'], \'font.size\': 10.4, \'axes.labelsize\': 10.4, \'axes.titlesize\': 11.05, \'xtick.labelsize\': 9.1, \'ytick.labelsize\': 9.1, \'legend.fontsize\': 9.49, \'pdf.fonttype\': 42, \'ps.fonttype\': 42, \'svg.fonttype\': \'none\', \'axes.linewidth\': 0.975, \'axes.edgecolor\': \'#202020\', \'text.color\': \'#202020\', \'axes.labelcolor\': \'#202020\', \'savefig.facecolor\': \'white\', \'figure.facecolor\': \'white\', \'svg.hashsalt\': \'scrna-style-20261105-v1\'})\n    dslist = [\'GSE289422\', \'GSE255566\']\n    cohort = {\'GSE289422\': \'CAD cohort\', \'GSE255566\': \'T2D cohort\'}\n    disease = {\'GSE289422\': \'CAD\', \'GSE255566\': \'T2D\'}\n    cohortcolor = {\'GSE289422\': \'#B65F64\', \'GSE255566\': \'#5683A8\'}\n    genes = [\'FCRL6\', \'GZMA\', \'GZMH\', \'KLRD1\']\n    dat, sourcehash, featurehash = ({}, {}, {})\n    for ds in dslist:\n        path = inp / f\'results/{ds}/coordinates.csv.gz\'\n        if not path.is_file():\n            raise FileNotFoundError(path)\n        sourcehash[str(path)] = sha(path)\n        assert sourcehash[str(path)] == original[\'coordinate_input_sha256\'][str(path)], \'Frozen input changed\'\n        d = pd.read_csv(path)\n        assert d.cell_id.is_unique and d.groupby(\'sample\').size().le(1000).all()\n        assert set(d.group) == {\'Control\', disease[ds]}\n        dat[ds] = d\n        featurehash[ds] = {\'cell_id_order\': hashlib.sha256(\'\\n\'.join(d.cell_id.astype(str)).encode()).hexdigest(), \'UMAP_coordinates_float64\': hashlib.sha256(d[[\'UMAP1\', \'UMAP2\']].to_numpy(dtype=\'<f8\').tobytes()).hexdigest(), \'expression_float64\': {gene: hashlib.sha256(d[gene + \'_logCP10k\'].to_numpy(dtype=\'<f8\').tobytes()).hexdigest() for gene in genes}}\n    limits = pd.read_csv(ref / \'metadata/feature_colour_limits.csv\')\n    caps = {}\n    for ds in dslist:\n        for gene in genes:\n            row = limits.loc[(limits.dataset == ds) & (limits.gene == gene)]\n            assert len(row) == 1\n            values = dat[ds][gene + \'_logCP10k\'].to_numpy()\n            positive = values[values > 0]\n            cap = float(np.quantile(positive, 0.99)) if len(positive) else 1.0\n            assert float(row.iloc[0].vmin) == 0\n            assert np.isclose(float(row.iloc[0].vmax), cap, rtol=1e-15, atol=0)\n            caps[ds, gene] = cap\n    cmap_stops = [0, 0.15, 0.45, 0.72, 1]\n    cmap_colors = [\'#ECEDEF\', \'#F1D8BD\', \'#D78D67\', \'#C55432\', \'#7651A6\']\n    cmap = mcolors.LinearSegmentedColormap.from_list(\'YY_gene\', list(zip(cmap_stops, cmap_colors)), N=256)\n    cols = [\'UMAP1\', \'UMAP2\']\n    centers = {ds: ((dat[ds][cols].max() + dat[ds][cols].min()) / 2).to_numpy() for ds in dslist}\n    span = max((float(np.ptp(dat[ds][col])) for ds in dslist for col in cols)) * 1.12\n    local = {ds: dat[ds][cols].to_numpy() - centers[ds] for ds in dslist}\n    ledger, axes_records, text_records = ([], [], [])\n    fw, fh = (7.8, 5.12)\n    fig = plt.figure(figsize=(fw, fh))\n    left, pw, gx = (0.075, 0.204, 0.032)\n    ph = pw * fw / fh\n    rowbottom = [0.573, 0.128]\n    feature_axes, bar_axes = ([], [])\n    for row, ds in enumerate(dslist):\n        d, xy = (dat[ds], local[ds])\n        text_records.append(fig.text(0.022, rowbottom[row] + ph / 2, cohort[ds] + \'\\n\' + ds, fontsize=font_sizes[\'cohort_label\'], fontweight=\'bold\', color=cohortcolor[ds], rotation=90, ha=\'center\', va=\'center\'))\n        for col, gene in enumerate(genes):\n            x0, y0 = (left + col * (pw + gx), rowbottom[row])\n            ax = fig.add_axes([x0, y0, pw, ph])\n            feature_axes.append((row, col, ax))\n            ax.set_aspect(\'equal\', adjustable=\'box\')\n            ax.set_xticks([])\n            ax.set_yticks([])\n            for spine in ax.spines.values():\n                spine.set_visible(False)\n            ax.set_xlim(-span / 2, span / 2)\n            ax.set_ylim(-span / 2, span / 2)\n            values = d[gene + \'_logCP10k\'].to_numpy()\n            ix = np.argsort(values, kind=\'stable\')\n            sc = ax.scatter(xy[ix, 0], xy[ix, 1], c=values[ix], s=0.68, alpha=0.9, linewidths=0, rasterized=True, cmap=cmap, vmin=0, vmax=caps[ds, gene], zorder=1)\n            assert np.array_equal(sc.get_offsets(), xy[ix])\n            assert np.array_equal(sc.get_array(), values[ix])\n            ax.set_title(gene, fontsize=font_sizes[\'gene_heading\'], fontstyle=\'italic\', fontweight=\'bold\', pad=4)\n            ax.set_xlabel(\'UMAP 1\', fontsize=font_sizes[\'axis_label\'], fontweight=\'bold\', labelpad=0)\n            ax.set_ylabel(\'UMAP 2\' if col == 0 else \'\', fontsize=font_sizes[\'axis_label\'], fontweight=\'bold\', labelpad=0)\n            cax = fig.add_axes([x0 + 0.018, y0 - 0.047, pw - 0.036, 0.012])\n            bar_axes.append(cax)\n            cb = fig.colorbar(sc, cax=cax, orientation=\'horizontal\', ticks=[0, caps[ds, gene] / 2, caps[ds, gene]])\n            cb.outline.set_linewidth(0.35)\n            cb.ax.tick_params(labelsize=font_sizes[\'colorbar_tick\'], length=1.5, width=0.35, pad=1.5)\n            cb.set_ticklabels([f\'{v:.2g}\' for v in [0, caps[ds, gene] / 2, caps[ds, gene]]])\n            ledger.append({\'algorithm\': \'UMAP\', \'view\': f\'FeaturePlot_{gene}\', \'dataset\': ds, \'group\': \'all\', \'source_x\': cols[0], \'source_y\': cols[1], \'dx\': float(-centers[ds][0]), \'dy\': float(-centers[ds][1]), \'scale_x\': 1, \'scale_y\': 1, \'rotation_degrees\': 0, \'display_cells\': len(d), \'slot\': col, \'slot_span\': span, \'gap\': 0})\n    text_records.append(fig.text(0.53, 0.02, \'Expression: log1p(CP10k); each cohort–gene scale is independent\', ha=\'center\', fontsize=font_sizes[\'existing_expression_note\']))\n    fig.canvas.draw()\n    for row in [0, 1]:\n        boxes = [ax.get_position() for r, c, ax in feature_axes if r == row]\n        assert np.allclose([b.y0 for b in boxes], boxes[0].y0)\n        assert np.allclose([b.height for b in boxes], boxes[0].height)\n    for col in range(4):\n        boxes = [ax.get_position() for r, c, ax in feature_axes if c == col]\n        assert np.allclose([b.x0 for b in boxes], boxes[0].x0)\n        assert np.allclose([b.width for b in boxes], boxes[0].width)\n    reference_ledger = pd.read_csv(ref / \'metadata/coordinate_translation_ledger.csv\')\n    reference_ledger = reference_ledger.loc[reference_ledger.algorithm == \'UMAP\'].reset_index(drop=True)\n    pd.testing.assert_frame_equal(pd.DataFrame(ledger), reference_ledger, check_dtype=False, check_exact=False, rtol=1e-14)\n    renderer = fig.canvas.get_renderer()\n    all_text = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]\n    for item in all_text:\n        b = item.get_window_extent(renderer).transformed(fig.transFigure.inverted())\n        assert min(b.x0, b.y0) >= 0 and max(b.x1, b.y1) <= 1, \'Text clipped: \' + item.get_text()\n        text_records.append({\'text\': item.get_text(), \'font_size_pt\': item.get_fontsize(), \'bounds_normalized\': [b.x0, b.y0, b.width, b.height]})\n    text_records = [r for r in text_records if isinstance(r, dict)]\n    for ax in fig.axes:\n        b = ax.get_position()\n        axes_records.append([b.x0, b.y0, b.width, b.height])\n    files = []\n    for ext, dpi in [(\'.pdf\', 600), (\'.svg\', 600), (\'_600dpi.png\', 600)]:\n        path = figdir / (\'Figure_features_UMAP\' + ext)\n        if path.exists():\n            raise FileExistsError(path)\n        fig.savefig(path, dpi=dpi)\n        files.append(path)\n        print(\'Saved\', path, flush=True)\n    plt.close(fig)\n    for path, digest in sourcehash.items():\n        assert sha(path) == digest\n    with (meta / \'feature_colour_limits.csv\').open(\'x\') as f:\n        f.write((ref / \'metadata/feature_colour_limits.csv\').read_text())\n    with (meta / \'coordinate_translation_ledger.csv\').open(\'x\') as f:\n        pd.DataFrame(ledger).to_csv(f, index=False)\n    write_json(meta / \'RENDER_MANIFEST.json\', {\'source_figure\': str(ref / \'figures/Figure_features_UMAP_600dpi.png\'), \'source_figure_sha256\': sha(ref / \'figures/Figure_features_UMAP_600dpi.png\'), \'source_entrypoint\': str(original_source), \'source_entrypoint_sha256\': sha(original_source), \'script\': str(Path(__file__).resolve()), \'script_sha256\': sha(__file__), \'coordinate_input_sha256\': sourcehash, \'data_array_sha256\': featurehash, \'display_cells\': {ds: len(d) for ds, d in dat.items()}, \'fixed_genes\': genes, \'font\': \'Helvetica\', \'font_sizes_pt\': font_sizes, \'font_increase_percent\': 30, \'cohort_label_colors\': cohortcolor, \'expression_colormap\': {\'name\': \'YY_gene\', \'positions\': cmap_stops, \'colors\': cmap_colors, \'N\': 256}, \'expression_colormap_unchanged\': True, \'expression_limits_unchanged\': True, \'expression_values_unchanged\': True, \'coordinate_values_unchanged\': True, \'only_original_coordinate_translation\': True, \'point_size\': 0.68, \'point_alpha\': 0.9, \'stable_expression_point_order_unchanged\': True, \'no_new_titles_letters_or_notes\': True, \'raw_h5ad_not_read\': True, \'no_fitting_or_new_embeddings\': True, \'figure_size_inches\': [fw, fh], \'figure_width_mm\': fw * 25.4, \'layout_change\': \'All text +30%; cohort colors fixed; row/column whitespace adjusted only\', \'equal_square_panels\': True, \'row_column_edges_aligned\': True, \'axes_normalized_bounds\': axes_records, \'text_audit\': text_records, \'no_clipped_text\': True, \'raster_dpi\': 600, \'editable_vector_text\': True, \'software_versions\': {\'python\': sys.version, \'numpy\': np.__version__, \'pandas\': pd.__version__, \'matplotlib\': matplotlib.__version__}, \'files\': {str(f): {\'sha256\': sha(f), \'bytes\': f.stat().st_size} for f in files}, \'source_inputs_unchanged\': True})\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
