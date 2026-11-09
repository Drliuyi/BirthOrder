"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import importlib.util
import inspect
import json
import os
from pathlib import Path
import sys
FONT_MULTIPLIER = 1.3
STAR_MULTIPLIER = 1.5
COHORT_COLORS = {'GSE289422': '#B65F64', 'GSE255566': '#5683A8'}

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def replace_one(source, before, after):
    if source.count(before) != 1:
        raise ValueError(f'Frozen drawing block changed: {before[:100]!r}')
    return source.replace(before, after)

def baseline_source(module):
    """Apply the precise original render_box_stars_unlettered_v3.py changes."""
    source = inspect.getsource(module.donor_plot)
    before = '                ax.plot([x - .34, x - .34, x + .34, x + .34],\n                        [y - .14, y, y, y - .14], color="#3C3C3C", lw=.55, clip_on=False)\n                ax.text(x, y + .15, annotation.annotation, ha="center", va="bottom", fontsize=6.25)'
    after = '                sig = stars(float(annotation.P))\n                if sig:\n                    ax.plot([x - .34, x - .34, x + .34, x + .34],\n                            [y - .14, y, y, y - .14], color="#3C3C3C", lw=.55, clip_on=False)\n                    ax.text(x, y + .15, sig, ha="center", va="bottom", fontsize=8.5, fontweight="bold")'
    source = replace_one(source, before, after)
    panel = '\n'.join(['        ax.text(-.086, 1.075, chr(ord("A") + rownum), transform=ax.transAxes,', '                fontsize=11, fontweight="bold", ha="left", va="bottom")'])
    return replace_one(source, panel, '')

def styled_source(source):
    source = replace_one(source, '                    color = GROUP_COLORS[group]', '                    color = COHORT_COLORS[ds]\n                    is_control = group == "Control"')
    source = replace_one(source, 'matplotlib.colors.to_rgba(color, .16)', 'matplotlib.colors.to_rgba(color, .06 if is_control else .28)')
    source = replace_one(source, '                    ax.scatter(xs, values, s=9, color=color, alpha=.68, linewidth=0, zorder=3)', '                    ax.scatter(xs, values, s=9, facecolors="none" if is_control else color,\n                               edgecolors=color, alpha=.80, linewidth=.60 if is_control else 0, zorder=3)')
    source = replace_one(source, 'fontsize=8.5, fontweight="bold")', 'fontsize=12.75, fontweight="bold")')
    source = replace_one(source, 'fontsize=8.5, fontweight="bold", ha="left", va="bottom")', 'fontsize=11.05, fontweight="bold", ha="left", va="bottom")')
    source = replace_one(source, 'fontsize=6.65, fontweight="bold", va="center", ha="center", clip_on=False)', 'fontsize=8.645, fontweight="bold", va="center", ha="center", clip_on=False)')
    old_legend = '    fig.legend(handles=[Patch(facecolor=matplotlib.colors.to_rgba(BLUE, .16), edgecolor=BLUE, label="Control"),\n                        Patch(facecolor=matplotlib.colors.to_rgba(ORANGE, .16), edgecolor=ORANGE, label="Disease")],'
    new_legend = '    fig.legend(handles=[Patch(facecolor="none", edgecolor="#686868", label="Control"),\n                        Patch(facecolor="#686868", edgecolor="#686868", label="Disease")],'
    source = replace_one(source, old_legend, new_legend)
    return replace_one(source, 'fontsize=7.4, handlelength=1.6, columnspacing=2.0)', 'fontsize=9.62, handlelength=1.6, columnspacing=2.0)')

def capture_geometry(fig, axes):
    """Capture every donor, box/whisker/bracket coordinate and annotation."""
    fig.canvas.draw()
    return {'canvas_inches': fig.get_size_inches().tolist(), 'axes_bounds': [list(ax.get_position().bounds) for ax in axes], 'limits': [[list(ax.get_xlim()), list(ax.get_ylim())] for ax in axes], 'point_offsets': [[artist.get_offsets().tolist() for artist in ax.collections] for ax in axes], 'lines': [[{'x': artist.get_xdata().tolist(), 'y': artist.get_ydata().tolist()} for artist in ax.lines] for ax in axes], 'box_vertices': [[patch.get_path().vertices.tolist() for patch in ax.patches[:16]] for ax in axes], 'stars': [[{'text': text.get_text(), 'position': list(text.get_position())} for text in ax.texts if text.get_text() in ('*', '**', '***')] for ax in axes], 'ticks': [{'x': ax.get_xticks().tolist(), 'y': ax.get_yticks().tolist()} for ax in axes], 'labels': [{'x': [text.get_text() for text in ax.get_xticklabels()], 'y': [text.get_text() for text in ax.get_yticklabels()], 'ylabel': ax.get_ylabel(), 'text': [text.get_text() for text in ax.texts]} for ax in axes]}

def check_visible_text(fig, axes):
    """Check all nonempty labels/stars remain inside the unchanged canvas."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    boxes = []
    texts = list(fig.legends[0].get_texts())
    for ax in axes:
        texts += [ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels(), *ax.texts]
    for text in texts:
        if not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < canvas.x0 - 0.1 or box.x1 > canvas.x1 + 0.1 or box.y0 < canvas.y0 - 0.1 or (box.y1 > canvas.y1 + 0.1):
            raise ValueError(f'Text is clipped by canvas: {text.get_text()}')
        boxes.append({'text': text.get_text(), 'fontsize_pt': text.get_fontsize(), 'color': text.get_color(), 'bbox_pixels': list(box.bounds)})
    return boxes

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('input-root', 'script-root', 'output-root', 'reference-root'):
        default = os.environ.get('SCRNA_BOX_STYLE_' + key.replace('-', '_').upper())
        parser.add_argument('--' + key, default=default, required=default is None)
    args = parser.parse_args()
    roots = {key: Path(value).expanduser().resolve() for key, value in vars(args).items()}
    print(json.dumps({key: str(value) for key, value in roots.items()}), flush=True)
    inp, scripts, out, reference = [roots[key] for key in ('input_root', 'script_root', 'output_root', 'reference_root')]
    source_data = inp / 'donor_expression_plot_data.csv'
    source_p = inp / 'frozen_P_annotations.csv'
    renderer_path = scripts / 'render_donor_effects.py'
    unlettered_path = scripts / 'render_box_stars_unlettered_v3.py'
    frozen_manifest = reference / 'MANIFEST.json'
    required = [source_data, source_p, renderer_path, unlettered_path, frozen_manifest, reference / 'Fig_donor_expression_core4_600dpi.png']
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f'Missing frozen input: {path}')
    if out.exists():
        raise FileExistsError(f'Refusing existing output root: {out}')
    if out.is_relative_to(inp) or out.is_relative_to(scripts) or out.is_relative_to(reference):
        raise ValueError('Output root must not be inside frozen input/source/reference directories')
    original_manifest = json.loads(frozen_manifest.read_text())
    input_hashes = {str(path): sha(path) for path in required}
    for path in (source_data, source_p, renderer_path):
        assert input_hashes[str(path)] == original_manifest['input_sha256'][str(path)]
    assert input_hashes[str(unlettered_path)] == original_manifest['script_sha256']
    png_record = next((row for row in original_manifest['outputs'] if row['path'].endswith('_600dpi.png')))
    assert input_hashes[str(required[-1])] == png_record['sha256']
    out.mkdir(parents=True)
    os.environ.setdefault('MPLCONFIGDIR', str(out / 'matplotlib_cache'))
    spec = importlib.util.spec_from_file_location('frozen_box_renderer', renderer_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    import matplotlib.pyplot as plt
    import pandas as pd
    data = pd.read_csv(source_data)
    annotations = pd.read_csv(source_p)
    assert len(data) == 400 and len(annotations) == 16
    saved_data = data.copy(deep=True)
    saved_annotations = annotations.copy(deep=True)
    source = baseline_source(module)
    baseline = {}

    def baseline_save(fig, axes, destination, stem, expected_mm):
        baseline.update(capture_geometry(fig, axes))
        plt.close(fig)
        return {'axes_bounds': baseline['axes_bounds']}
    namespace = dict(module.__dict__)
    namespace['save_figure'] = baseline_save
    exec(compile(source, str(unlettered_path), 'exec'), namespace)
    font = module.configure_style()
    namespace['donor_plot'](data, annotations, out)
    audit = {}

    def checked_save(fig, axes, destination, stem, expected_mm):
        geometry = capture_geometry(fig, axes)
        assert geometry == baseline, 'Coordinates, values, labels or frozen annotations changed'
        assert sum((len(points) for row in geometry['point_offsets'] for points in row)) == 400
        assert sum((len(row) for row in geometry['stars'])) == 5
        assert len(axes[0].collections) == len(axes[1].collections) == 16
        text_audit = check_visible_text(fig, axes)
        for ax in axes:
            for collection, ds in zip(ax.collections, [dataset for dataset in module.DATASETS for _ in range(8)]):
                expected = module.matplotlib.colors.to_rgba(COHORT_COLORS[ds], 0.8)
                assert module.np.allclose(collection.get_edgecolors()[0], expected)
            assert [patch.get_facecolor() for patch in ax.patches[-2:]] == [module.matplotlib.colors.to_rgba(COHORT_COLORS[ds]) for ds in module.DATASETS]
        audit.update(geometry=geometry, text=text_audit, all_400_coordinates_equal=True, all_box_whisker_bracket_coordinates_equal=True, all_16_P_q_equal=True, original_5_star_positions_equal=True, all_visible_text_inside_canvas=True)
        return module.save_figure(fig, axes, destination, stem, expected_mm)
    module.configure_style()
    for key in ('font.size', 'axes.labelsize', 'axes.titlesize', 'xtick.labelsize', 'ytick.labelsize'):
        plt.rcParams[key] *= FONT_MULTIPLIER
    for key in ('text.color', 'axes.labelcolor', 'xtick.color', 'ytick.color'):
        plt.rcParams[key] = '#000000'
    namespace = dict(module.__dict__)
    namespace.update(COHORT_COLORS=COHORT_COLORS, save_figure=checked_save)
    exec(compile(styled_source(source), str(Path(__file__).resolve()), 'exec'), namespace)
    alignment = namespace['donor_plot'](data, annotations, out)
    alignment['group_colors'] = {'CAD Control': COHORT_COLORS['GSE289422'], 'CAD Disease': COHORT_COLORS['GSE289422'], 'T2D Control': COHORT_COLORS['GSE255566'], 'T2D Disease': COHORT_COLORS['GSE255566']}
    pd.testing.assert_frame_equal(data, saved_data, check_exact=True)
    pd.testing.assert_frame_equal(annotations, saved_annotations, check_exact=True)
    for path in required:
        assert sha(path) == input_hashes[str(path)], f'Frozen input changed: {path}'
    annotations.to_csv(out / 'FROZEN_P_ANNOTATIONS_UNCHANGED.csv', index=False, mode='x')
    with (out / 'STYLE_AND_GEOMETRY_AUDIT.json').open('x') as handle:
        json.dump(audit, handle, indent=2)
    files = [path for path in sorted(out.iterdir()) if path.is_file()]
    manifest = {'scope': 'Style only: fonts +30%, stars +50%; CAD red, T2D blue in both control and disease', 'roots': {key: str(value) for key, value in roots.items()}, 'input_sha256': input_hashes, 'script_sha256': sha(Path(__file__).resolve()), 'frozen_source_verified_against_original_manifest': True, 'font': font, 'font_multiplier': FONT_MULTIPLIER, 'star_multiplier': STAR_MULTIPLIER, 'colors': COHORT_COLORS, 'control': 'cohort hue, hollow donor circles and 6% box fill', 'disease': 'cohort hue, solid donor circles and 28% box fill', 'legend': 'neutral hollow Control / filled Disease keys', 'alignment': alignment, 'all_400_expression_points_and_coordinates_preserved': True, 'all_16_original_P_q_preserved': True, 'original_significance_brackets_and_stars_preserved': True, 'panel_lettering': False, 'added_title': False, 'added_note': False, 'new_tests': 0, 'model_refits': 0, 'versions': {name: importlib.metadata.version(name) for name in ('matplotlib', 'numpy', 'pandas')}, 'outputs': [{'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in files]}
    with (out / 'MANIFEST.json').open('x') as handle:
        json.dump(manifest, handle, indent=2)
    print(json.dumps({'status': 'COMPLETE', 'output_root': str(out), 'frozen_points': 400, 'frozen_stars': 5, 'new_tests': 0, 'model_refits': 0}), flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Style-only rendering of the frozen 400-point core4 donor expression boxplot.\n\nUses the exact unlettered-v3 drawing function as the baseline. No expression,\nnormalization, testing, adjustment or model fit is performed. A new output root\nis mandatory. The previous source/data/figure hashes are checked before drawing.\n"""\nfrom __future__ import annotations\nimport argparse\nimport hashlib\nimport importlib.metadata\nimport importlib.util\nimport inspect\nimport json\nimport os\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\nFONT_MULTIPLIER = 1.3\nSTAR_MULTIPLIER = 1.5\nCOHORT_COLORS = {\'GSE289422\': \'#B65F64\', \'GSE255566\': \'#5683A8\'}\n\ndef sha(path):\n    digest = hashlib.sha256()\n    with Path(path).open(\'rb\') as handle:\n        for chunk in iter(lambda: handle.read(1024 * 1024), b\'\'):\n            digest.update(chunk)\n    return digest.hexdigest()\n\ndef replace_one(source, before, after):\n    if source.count(before) != 1:\n        raise ValueError(f\'Frozen drawing block changed: {before[:100]!r}\')\n    return source.replace(before, after)\n\ndef baseline_source(module):\n    """Apply the precise original render_box_stars_unlettered_v3.py changes."""\n    source = inspect.getsource(module.donor_plot)\n    before = \'                ax.plot([x - .34, x - .34, x + .34, x + .34],\\n                        [y - .14, y, y, y - .14], color="#3C3C3C", lw=.55, clip_on=False)\\n                ax.text(x, y + .15, annotation.annotation, ha="center", va="bottom", fontsize=6.25)\'\n    after = \'                sig = stars(float(annotation.P))\\n                if sig:\\n                    ax.plot([x - .34, x - .34, x + .34, x + .34],\\n                            [y - .14, y, y, y - .14], color="#3C3C3C", lw=.55, clip_on=False)\\n                    ax.text(x, y + .15, sig, ha="center", va="bottom", fontsize=8.5, fontweight="bold")\'\n    source = replace_one(source, before, after)\n    panel = \'\\n\'.join([\'        ax.text(-.086, 1.075, chr(ord("A") + rownum), transform=ax.transAxes,\', \'                fontsize=11, fontweight="bold", ha="left", va="bottom")\'])\n    return replace_one(source, panel, \'\')\n\ndef styled_source(source):\n    source = replace_one(source, \'                    color = GROUP_COLORS[group]\', \'                    color = COHORT_COLORS[ds]\\n                    is_control = group == "Control"\')\n    source = replace_one(source, \'matplotlib.colors.to_rgba(color, .16)\', \'matplotlib.colors.to_rgba(color, .06 if is_control else .28)\')\n    source = replace_one(source, \'                    ax.scatter(xs, values, s=9, color=color, alpha=.68, linewidth=0, zorder=3)\', \'                    ax.scatter(xs, values, s=9, facecolors="none" if is_control else color,\\n                               edgecolors=color, alpha=.80, linewidth=.60 if is_control else 0, zorder=3)\')\n    source = replace_one(source, \'fontsize=8.5, fontweight="bold")\', \'fontsize=12.75, fontweight="bold")\')\n    source = replace_one(source, \'fontsize=8.5, fontweight="bold", ha="left", va="bottom")\', \'fontsize=11.05, fontweight="bold", ha="left", va="bottom")\')\n    source = replace_one(source, \'fontsize=6.65, fontweight="bold", va="center", ha="center", clip_on=False)\', \'fontsize=8.645, fontweight="bold", va="center", ha="center", clip_on=False)\')\n    old_legend = \'    fig.legend(handles=[Patch(facecolor=matplotlib.colors.to_rgba(BLUE, .16), edgecolor=BLUE, label="Control"),\\n                        Patch(facecolor=matplotlib.colors.to_rgba(ORANGE, .16), edgecolor=ORANGE, label="Disease")],\'\n    new_legend = \'    fig.legend(handles=[Patch(facecolor="none", edgecolor="#686868", label="Control"),\\n                        Patch(facecolor="#686868", edgecolor="#686868", label="Disease")],\'\n    source = replace_one(source, old_legend, new_legend)\n    return replace_one(source, \'fontsize=7.4, handlelength=1.6, columnspacing=2.0)\', \'fontsize=9.62, handlelength=1.6, columnspacing=2.0)\')\n\ndef capture_geometry(fig, axes):\n    """Capture every donor, box/whisker/bracket coordinate and annotation."""\n    fig.canvas.draw()\n    return {\'canvas_inches\': fig.get_size_inches().tolist(), \'axes_bounds\': [list(ax.get_position().bounds) for ax in axes], \'limits\': [[list(ax.get_xlim()), list(ax.get_ylim())] for ax in axes], \'point_offsets\': [[artist.get_offsets().tolist() for artist in ax.collections] for ax in axes], \'lines\': [[{\'x\': artist.get_xdata().tolist(), \'y\': artist.get_ydata().tolist()} for artist in ax.lines] for ax in axes], \'box_vertices\': [[patch.get_path().vertices.tolist() for patch in ax.patches[:16]] for ax in axes], \'stars\': [[{\'text\': text.get_text(), \'position\': list(text.get_position())} for text in ax.texts if text.get_text() in (\'*\', \'**\', \'***\')] for ax in axes], \'ticks\': [{\'x\': ax.get_xticks().tolist(), \'y\': ax.get_yticks().tolist()} for ax in axes], \'labels\': [{\'x\': [text.get_text() for text in ax.get_xticklabels()], \'y\': [text.get_text() for text in ax.get_yticklabels()], \'ylabel\': ax.get_ylabel(), \'text\': [text.get_text() for text in ax.texts]} for ax in axes]}\n\ndef check_visible_text(fig, axes):\n    """Check all nonempty labels/stars remain inside the unchanged canvas."""\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    canvas = fig.bbox\n    boxes = []\n    texts = list(fig.legends[0].get_texts())\n    for ax in axes:\n        texts += [ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels(), *ax.texts]\n    for text in texts:\n        if not text.get_text():\n            continue\n        box = text.get_window_extent(renderer)\n        if box.x0 < canvas.x0 - 0.1 or box.x1 > canvas.x1 + 0.1 or box.y0 < canvas.y0 - 0.1 or (box.y1 > canvas.y1 + 0.1):\n            raise ValueError(f\'Text is clipped by canvas: {text.get_text()}\')\n        boxes.append({\'text\': text.get_text(), \'fontsize_pt\': text.get_fontsize(), \'color\': text.get_color(), \'bbox_pixels\': list(box.bounds)})\n    return boxes\n\ndef main():\n    parser = argparse.ArgumentParser(description=__doc__)\n    for key in (\'input-root\', \'script-root\', \'output-root\', \'reference-root\'):\n        default = os.environ.get(\'SCRNA_BOX_STYLE_\' + key.replace(\'-\', \'_\').upper())\n        parser.add_argument(\'--\' + key, default=default, required=default is None)\n    args = parser.parse_args()\n    roots = {key: Path(value).expanduser().resolve() for key, value in vars(args).items()}\n    print(json.dumps({key: str(value) for key, value in roots.items()}), flush=True)\n    inp, scripts, out, reference = [roots[key] for key in (\'input_root\', \'script_root\', \'output_root\', \'reference_root\')]\n    source_data = inp / \'donor_expression_plot_data.csv\'\n    source_p = inp / \'frozen_P_annotations.csv\'\n    renderer_path = scripts / \'render_donor_effects.py\'\n    unlettered_path = scripts / \'render_box_stars_unlettered_v3.py\'\n    frozen_manifest = reference / \'MANIFEST.json\'\n    required = [source_data, source_p, renderer_path, unlettered_path, frozen_manifest, reference / \'Fig_donor_expression_core4_600dpi.png\']\n    for path in required:\n        if not path.is_file():\n            raise FileNotFoundError(f\'Missing frozen input: {path}\')\n    if out.exists():\n        raise FileExistsError(f\'Refusing existing output root: {out}\')\n    if out.is_relative_to(inp) or out.is_relative_to(scripts) or out.is_relative_to(reference):\n        raise ValueError(\'Output root must not be inside frozen input/source/reference directories\')\n    original_manifest = json.loads(frozen_manifest.read_text())\n    input_hashes = {str(path): sha(path) for path in required}\n    for path in (source_data, source_p, renderer_path):\n        assert input_hashes[str(path)] == original_manifest[\'input_sha256\'][str(path)]\n    assert input_hashes[str(unlettered_path)] == original_manifest[\'script_sha256\']\n    png_record = next((row for row in original_manifest[\'outputs\'] if row[\'path\'].endswith(\'_600dpi.png\')))\n    assert input_hashes[str(required[-1])] == png_record[\'sha256\']\n    out.mkdir(parents=True)\n    os.environ.setdefault(\'MPLCONFIGDIR\', str(out / \'matplotlib_cache\'))\n    spec = importlib.util.spec_from_file_location(\'frozen_box_renderer\', renderer_path)\n    module = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(module)\n    import matplotlib.pyplot as plt\n    import pandas as pd\n    data = pd.read_csv(source_data)\n    annotations = pd.read_csv(source_p)\n    assert len(data) == 400 and len(annotations) == 16\n    saved_data = data.copy(deep=True)\n    saved_annotations = annotations.copy(deep=True)\n    source = baseline_source(module)\n    baseline = {}\n\n    def baseline_save(fig, axes, destination, stem, expected_mm):\n        baseline.update(capture_geometry(fig, axes))\n        plt.close(fig)\n        return {\'axes_bounds\': baseline[\'axes_bounds\']}\n    namespace = dict(module.__dict__)\n    namespace[\'save_figure\'] = baseline_save\n    exec(compile(source, str(unlettered_path), \'exec\'), namespace)\n    font = module.configure_style()\n    namespace[\'donor_plot\'](data, annotations, out)\n    audit = {}\n\n    def checked_save(fig, axes, destination, stem, expected_mm):\n        geometry = capture_geometry(fig, axes)\n        assert geometry == baseline, \'Coordinates, values, labels or frozen annotations changed\'\n        assert sum((len(points) for row in geometry[\'point_offsets\'] for points in row)) == 400\n        assert sum((len(row) for row in geometry[\'stars\'])) == 5\n        assert len(axes[0].collections) == len(axes[1].collections) == 16\n        text_audit = check_visible_text(fig, axes)\n        for ax in axes:\n            for collection, ds in zip(ax.collections, [dataset for dataset in module.DATASETS for _ in range(8)]):\n                expected = module.matplotlib.colors.to_rgba(COHORT_COLORS[ds], 0.8)\n                assert module.np.allclose(collection.get_edgecolors()[0], expected)\n            assert [patch.get_facecolor() for patch in ax.patches[-2:]] == [module.matplotlib.colors.to_rgba(COHORT_COLORS[ds]) for ds in module.DATASETS]\n        audit.update(geometry=geometry, text=text_audit, all_400_coordinates_equal=True, all_box_whisker_bracket_coordinates_equal=True, all_16_P_q_equal=True, original_5_star_positions_equal=True, all_visible_text_inside_canvas=True)\n        return module.save_figure(fig, axes, destination, stem, expected_mm)\n    module.configure_style()\n    for key in (\'font.size\', \'axes.labelsize\', \'axes.titlesize\', \'xtick.labelsize\', \'ytick.labelsize\'):\n        plt.rcParams[key] *= FONT_MULTIPLIER\n    for key in (\'text.color\', \'axes.labelcolor\', \'xtick.color\', \'ytick.color\'):\n        plt.rcParams[key] = \'#000000\'\n    namespace = dict(module.__dict__)\n    namespace.update(COHORT_COLORS=COHORT_COLORS, save_figure=checked_save)\n    exec(compile(styled_source(source), str(Path(__file__).resolve()), \'exec\'), namespace)\n    alignment = namespace[\'donor_plot\'](data, annotations, out)\n    alignment[\'group_colors\'] = {\'CAD Control\': COHORT_COLORS[\'GSE289422\'], \'CAD Disease\': COHORT_COLORS[\'GSE289422\'], \'T2D Control\': COHORT_COLORS[\'GSE255566\'], \'T2D Disease\': COHORT_COLORS[\'GSE255566\']}\n    pd.testing.assert_frame_equal(data, saved_data, check_exact=True)\n    pd.testing.assert_frame_equal(annotations, saved_annotations, check_exact=True)\n    for path in required:\n        assert sha(path) == input_hashes[str(path)], f\'Frozen input changed: {path}\'\n    annotations.to_csv(out / \'FROZEN_P_ANNOTATIONS_UNCHANGED.csv\', index=False, mode=\'x\')\n    with (out / \'STYLE_AND_GEOMETRY_AUDIT.json\').open(\'x\') as handle:\n        json.dump(audit, handle, indent=2)\n    files = [path for path in sorted(out.iterdir()) if path.is_file()]\n    manifest = {\'scope\': \'Style only: fonts +30%, stars +50%; CAD red, T2D blue in both control and disease\', \'roots\': {key: str(value) for key, value in roots.items()}, \'input_sha256\': input_hashes, \'script_sha256\': sha(Path(__file__).resolve()), \'frozen_source_verified_against_original_manifest\': True, \'font\': font, \'font_multiplier\': FONT_MULTIPLIER, \'star_multiplier\': STAR_MULTIPLIER, \'colors\': COHORT_COLORS, \'control\': \'cohort hue, hollow donor circles and 6% box fill\', \'disease\': \'cohort hue, solid donor circles and 28% box fill\', \'legend\': \'neutral hollow Control / filled Disease keys\', \'alignment\': alignment, \'all_400_expression_points_and_coordinates_preserved\': True, \'all_16_original_P_q_preserved\': True, \'original_significance_brackets_and_stars_preserved\': True, \'panel_lettering\': False, \'added_title\': False, \'added_note\': False, \'new_tests\': 0, \'model_refits\': 0, \'versions\': {name: importlib.metadata.version(name) for name in (\'matplotlib\', \'numpy\', \'pandas\')}, \'outputs\': [{\'path\': str(path), \'bytes\': path.stat().st_size, \'sha256\': sha(path)} for path in files]}\n    with (out / \'MANIFEST.json\').open(\'x\') as handle:\n        json.dump(manifest, handle, indent=2)\n    print(json.dumps({\'status\': \'COMPLETE\', \'output_root\': str(out), \'frozen_points\': 400, \'frozen_stars\': 5, \'new_tests\': 0, \'model_refits\': 0}), flush=True)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
