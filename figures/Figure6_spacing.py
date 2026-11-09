"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    source = Path(__file__).resolve().parents[1] / 'umap/draw_umap_expression_style_v1.py'
    if not source.is_file():
        raise FileNotFoundError(f'Missing unchanged v1 renderer: {source}')
    source_hash = digest(source)
    spec = importlib.util.spec_from_file_location('frozen_umap_style_v1', source)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    base.__file__ = str(Path(__file__).resolve())
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    import numpy as np
    old_add_artist, old_text, old_savefig = (Figure.add_artist, Figure.text, Figure.savefig)
    delta_fraction = 0.03
    captured = {'cohort_band_count': 0, 'cohort_name_count': 0}

    def add_artist(self, artist, *args, **kwargs):
        if isinstance(artist, Line2D) and np.array_equal(artist.get_ydata(), [0.181, 0.181]):
            artist.set_ydata(np.asarray(artist.get_ydata()) + delta_fraction)
            captured['cohort_band_count'] += 1
        return old_add_artist(self, artist, *args, **kwargs)

    def text(self, x, y, s, *args, **kwargs):
        if s in ['CAD cohort · GSE289422', 'T2D cohort · GSE255566']:
            assert y == 0.16
            y += delta_fraction
            captured['cohort_name_count'] += 1
        return old_text(self, x, y, s, *args, **kwargs)

    def savefig(self, fname, *args, **kwargs):
        if 'output_root' not in captured:
            self.canvas.draw()
            renderer = self.canvas.get_renderer()
            to_mm = 25.4 / self.dpi
            cohort = [t for t in self.texts if 'cohort · GSE' in t.get_text()]
            legend_text = [t for legend in self.legends for t in legend.get_texts()]
            cohort_bottom = min((t.get_window_extent(renderer).y0 for t in cohort))
            legend_top = max((t.get_window_extent(renderer).y1 for t in legend_text))
            bar_top = (0.181 + delta_fraction) * self.bbox.height
            xlabel_bottom = self.axes[0].xaxis.label.get_window_extent(renderer).y0
            captured.update({'output_root': str(Path(fname).resolve().parents[1]), 'upward_shift_mm': delta_fraction * self.get_figheight() * 25.4, 'cohort_name_to_celltype_legend_clear_gap_mm': (cohort_bottom - legend_top) * to_mm, 'cohort_band_to_UMAP_xlabel_clear_gap_mm': (xlabel_bottom - bar_top) * to_mm, 'new_cohort_band_y_figure': 0.181 + delta_fraction, 'new_cohort_name_y_figure': 0.16 + delta_fraction, 'legend_bbox_to_anchor': list(self.legends[0].get_bbox_to_anchor().bounds), 'axes_bounds': list(self.axes[0].get_position().bounds)})
            assert captured['cohort_name_to_celltype_legend_clear_gap_mm'] > 2
            assert captured['cohort_band_to_UMAP_xlabel_clear_gap_mm'] > 1
        return old_savefig(self, fname, *args, **kwargs)
    Figure.add_artist, Figure.text, Figure.savefig = (add_artist, text, savefig)
    try:
        base.main()
    finally:
        Figure.add_artist, Figure.text, Figure.savefig = (old_add_artist, old_text, old_savefig)
    assert captured['cohort_band_count'] == 2
    assert captured['cohort_name_count'] == 2
    assert digest(source) == source_hash
    out = Path(captured.pop('output_root'))
    manifest = json.loads((out / 'metadata/RENDER_MANIFEST.json').read_text())
    audit = {'renderer': str(Path(__file__).resolve()), 'renderer_sha256': digest(__file__), 'unchanged_v1_source_dependency': str(source), 'v1_source_sha256': source_hash, 'change': 'Only both cohort-colored bands and their cohort names move upward together.', 'preserved': ['cell-type legend position', 'UMAP frame and axes labels', 'all font sizes', 'all colors', 'cell point coordinates and translations', 'cell IDs/samples/counts', 'original files'], 'geometry_and_spacing': captured, 'coordinate_checks': manifest['coordinate_checks'], 'source_coordinate_sha256': {path: value for path, value in manifest['input_sha256'].items() if path.endswith('coordinates.csv.gz')}, 'export_paths': list(manifest['files'])}
    with (out / 'metadata/SPACING_AUDIT.json').open('x') as handle:
        json.dump(audit, handle, indent=2)
        handle.write('\n')
    print(json.dumps(captured, indent=2), flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Reproduce UMAP style v1 with only the two cohort bands/names 2.58 mm higher."""\nimport hashlib\nimport importlib.util\nimport json\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\n\ndef digest(path):\n    return hashlib.sha256(Path(path).read_bytes()).hexdigest()\n\ndef main():\n    source = Path(__file__).resolve().parents[1] / \'umap/draw_umap_expression_style_v1.py\'\n    if not source.is_file():\n        raise FileNotFoundError(f\'Missing unchanged v1 renderer: {source}\')\n    source_hash = digest(source)\n    spec = importlib.util.spec_from_file_location(\'frozen_umap_style_v1\', source)\n    base = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(base)\n    base.__file__ = str(Path(__file__).resolve())\n    from matplotlib.figure import Figure\n    from matplotlib.lines import Line2D\n    import numpy as np\n    old_add_artist, old_text, old_savefig = (Figure.add_artist, Figure.text, Figure.savefig)\n    delta_fraction = 0.03\n    captured = {\'cohort_band_count\': 0, \'cohort_name_count\': 0}\n\n    def add_artist(self, artist, *args, **kwargs):\n        if isinstance(artist, Line2D) and np.array_equal(artist.get_ydata(), [0.181, 0.181]):\n            artist.set_ydata(np.asarray(artist.get_ydata()) + delta_fraction)\n            captured[\'cohort_band_count\'] += 1\n        return old_add_artist(self, artist, *args, **kwargs)\n\n    def text(self, x, y, s, *args, **kwargs):\n        if s in [\'CAD cohort · GSE289422\', \'T2D cohort · GSE255566\']:\n            assert y == 0.16\n            y += delta_fraction\n            captured[\'cohort_name_count\'] += 1\n        return old_text(self, x, y, s, *args, **kwargs)\n\n    def savefig(self, fname, *args, **kwargs):\n        if \'output_root\' not in captured:\n            self.canvas.draw()\n            renderer = self.canvas.get_renderer()\n            to_mm = 25.4 / self.dpi\n            cohort = [t for t in self.texts if \'cohort · GSE\' in t.get_text()]\n            legend_text = [t for legend in self.legends for t in legend.get_texts()]\n            cohort_bottom = min((t.get_window_extent(renderer).y0 for t in cohort))\n            legend_top = max((t.get_window_extent(renderer).y1 for t in legend_text))\n            bar_top = (0.181 + delta_fraction) * self.bbox.height\n            xlabel_bottom = self.axes[0].xaxis.label.get_window_extent(renderer).y0\n            captured.update({\'output_root\': str(Path(fname).resolve().parents[1]), \'upward_shift_mm\': delta_fraction * self.get_figheight() * 25.4, \'cohort_name_to_celltype_legend_clear_gap_mm\': (cohort_bottom - legend_top) * to_mm, \'cohort_band_to_UMAP_xlabel_clear_gap_mm\': (xlabel_bottom - bar_top) * to_mm, \'new_cohort_band_y_figure\': 0.181 + delta_fraction, \'new_cohort_name_y_figure\': 0.16 + delta_fraction, \'legend_bbox_to_anchor\': list(self.legends[0].get_bbox_to_anchor().bounds), \'axes_bounds\': list(self.axes[0].get_position().bounds)})\n            assert captured[\'cohort_name_to_celltype_legend_clear_gap_mm\'] > 2\n            assert captured[\'cohort_band_to_UMAP_xlabel_clear_gap_mm\'] > 1\n        return old_savefig(self, fname, *args, **kwargs)\n    Figure.add_artist, Figure.text, Figure.savefig = (add_artist, text, savefig)\n    try:\n        base.main()\n    finally:\n        Figure.add_artist, Figure.text, Figure.savefig = (old_add_artist, old_text, old_savefig)\n    assert captured[\'cohort_band_count\'] == 2\n    assert captured[\'cohort_name_count\'] == 2\n    assert digest(source) == source_hash\n    out = Path(captured.pop(\'output_root\'))\n    manifest = json.loads((out / \'metadata/RENDER_MANIFEST.json\').read_text())\n    audit = {\'renderer\': str(Path(__file__).resolve()), \'renderer_sha256\': digest(__file__), \'unchanged_v1_source_dependency\': str(source), \'v1_source_sha256\': source_hash, \'change\': \'Only both cohort-colored bands and their cohort names move upward together.\', \'preserved\': [\'cell-type legend position\', \'UMAP frame and axes labels\', \'all font sizes\', \'all colors\', \'cell point coordinates and translations\', \'cell IDs/samples/counts\', \'original files\'], \'geometry_and_spacing\': captured, \'coordinate_checks\': manifest[\'coordinate_checks\'], \'source_coordinate_sha256\': {path: value for path, value in manifest[\'input_sha256\'].items() if path.endswith(\'coordinates.csv.gz\')}, \'export_paths\': list(manifest[\'files\'])}\n    with (out / \'metadata/SPACING_AUDIT.json\').open(\'x\') as handle:\n        json.dump(audit, handle, indent=2)\n        handle.write(\'\\n\')\n    print(json.dumps(captured, indent=2), flush=True)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
