"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import sys

def main():
    script = Path(__file__).resolve()
    project_default = Path(os.environ.get('PROJECT_ROOT', str(script.parents[4])))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=project_default)
    parser.add_argument('--helper-dir', type=Path, default=Path(os.environ.get('FIGURE2_HELPER_DIR', str(script.parent))))
    parser.add_argument('--input', type=Path, default=Path(os.environ.get('FIGURE2_C_CURVES_CSV', str(project_default / 'reports/figure2_c_curve_20261103T151741/nature_v5/PREDICTION_CURVES.csv'))))
    parser.add_argument('--output', type=Path, default=Path(os.environ.get('FIGURE2_C_OUTPUT_ROOT', str(project_default / 'reports/figure2_typography_20261105_v1/C'))))
    args = parser.parse_args()
    project, helpers, source, output = [p.expanduser().resolve() for p in (args.project_root, args.helper_dir, args.input, args.output)]
    print(f'Project root: {project}\nInput aggregate: {source}\nHelper root: {helpers}\nOutput root: {output}')
    helper = helpers / 'frozen_curve_helpers.py'
    if not helper.is_file() or not source.is_file():
        raise FileNotFoundError(f'Required aggregate/helper missing: {source}; {helper}')
    if output.exists():
        raise FileExistsError(f'New isolated output directory required: {output}')
    spec = importlib.util.spec_from_file_location('frozen_curve_helpers', helper)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    base.FONT_PT = 7.0
    fields, rows, groups, source_hash, common_n, checks = base.read_curves(source)
    fonts = base.configure_fonts()
    fig = base.plt.figure(figsize=base.MAIN_SIZE, facecolor='white')
    ax = fig.add_axes((0.225, 0.295, 0.735, 0.64))
    base.style_axes(ax, base.y_limits(groups))
    ax.grid(False, axis='both', which='both')
    lines = base.draw_curves(ax, groups)
    ax.set_xlabel('Relative birth position', fontsize=7.5, fontweight='bold', labelpad=5)
    ax.set_ylabel('Standardized probability difference', fontsize=7.5, fontweight='bold', labelpad=5)
    handles = [base.Line2D([], [], label=base.NAMES[c], linewidth=1.0, **base.STYLES[c]) for c in ('DM', 'T2D')]
    handles.append(base.Line2D([], [], color='none', linestyle='none', label=' '))
    handles.extend((base.Line2D([], [], label=base.NAMES[c], linewidth=1.0, **base.STYLES[c]) for c in ('IHD', 'Angina', 'CAS')))
    legend = fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.525, 0.025), ncol=2, fontsize=7.0, frameon=False, borderpad=0, handlelength=2.0, handletextpad=0.45, columnspacing=0.9, labelspacing=0.5)
    visible = base.check_text_bounds(fig)
    renderer = fig.canvas.get_renderer()
    legend_box = legend.get_window_extent(renderer)
    xlabel_box = ax.xaxis.label.get_window_extent(renderer)
    if legend_box.y1 + 6 > xlabel_box.y0:
        raise AssertionError(f'Legend/x-label gap too small: {xlabel_box.y0 - legend_box.y1}')
    legend_texts = [t for t in legend.get_texts() if t.get_text().strip()]
    for i, a in enumerate(legend_texts):
        for b in legend_texts[i + 1:]:
            if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                raise AssertionError('Legend text overlaps')
    if any((line.get_visible() for line in ax.get_xgridlines() + ax.get_ygridlines())):
        raise AssertionError('Only the saved zero reference line should remain')
    if fig.texts or ax.get_title():
        raise AssertionError('No titles, subtitles, panel letters, or notes permitted')
    coordinates = base.check_lines(lines, groups)
    prior = project / 'reports/figure2_c_curve_20261103T151741/nature_v5/Panel_C_plot_audit.json'
    if not prior.is_file():
        raise FileNotFoundError(f'Frozen curve audit missing: {prior}')
    prior_hash = base.sha256(prior)
    previous = json.loads(prior.read_text())
    expected_coordinates = previous['main']['QA']['disease_line_coordinates']
    if coordinates != expected_coordinates:
        raise AssertionError('Curve coordinates differ from the original Nature panel')
    output.mkdir(parents=True, exist_ok=False)
    saved = base.save_formats(fig, output, 'Panel_C')
    base.plt.close(fig)
    with (output / 'Panel_C_plot_data.csv').open('x', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (output / 'Panel_C_plot_data.csv').open(newline='', encoding='utf-8') as f:
        if list(csv.DictReader(f)) != rows:
            raise AssertionError('Saved source cells changed')
    if base.sha256(source) != source_hash or base.sha256(prior) != prior_hash:
        raise AssertionError('An original input was modified')
    (output / 'code').mkdir()
    for path in (script, helper):
        with (output / 'code' / path.name).open('xb') as f:
            f.write(path.read_bytes())
    audit = {'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS', 'scope': 'Typography only; exact saved five-disease curves and intervals', 'source': str(source), 'source_sha256': source_hash, 'source_unchanged': True, 'previous_audit': str(prior), 'previous_audit_sha256': prior_hash, 'common_N_from_source': common_n, 'rows': len(rows), 'models_fitted': 0, 'predictions_computed': 0, 'intervals_computed': 0, 'font_family': 'Arial', 'fonts': fonts, 'old_font_pt': 6.1, 'body_and_legend_font_pt': 7.0, 'axis_title_font_pt': 7.5, 'axis_title_weight': 'bold', 'figure_inches': list(base.MAIN_SIZE), 'plot_bounds': list(ax.get_position().bounds), 'legend_columns': {'left': [base.NAMES[c] for c in ('DM', 'T2D')], 'right': [base.NAMES[c] for c in ('IHD', 'Angina', 'CAS')]}, 'legend_xlabel_gap_pixels': xlabel_box.y0 - legend_box.y1, 'no_clipped_or_overlapping_labels': True, 'no_new_title_or_notes': True, 'source_validation': checks, 'exact_original_coordinates': coordinates, 'CI_alpha': 0.085, 'styles': base.STYLES, 'outputs': saved, 'visible_text': visible, 'script_sha256': base.sha256(script), 'helper_sha256': base.sha256(helper), 'python': sys.executable}
    with (output / 'Panel_C_audit.json').open('x', encoding='utf-8') as f:
        json.dump(audit, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps({'status': 'PASS', 'output': str(output), 'body_font_pt': 7.0, 'bold_axis_title_pt': 7.5, 'source_rows': len(rows)}))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Typography-only revision of the five exact saved standardized curves.\n\nNo fitted model, prediction, confidence interval, or plotted coordinate changes.\nThe frozen helper reads only aggregate data and checks all 105 saved rows.\n"""\nfrom __future__ import annotations\nimport argparse\nimport csv\nfrom datetime import datetime, timezone\nimport importlib.util\nimport json\nimport os\nfrom pathlib import Path\nimport sys\n\ndef main():\n    script = Path(__file__).resolve()\n    project_default = Path(os.environ.get(\'PROJECT_ROOT\', str(script.parents[4])))\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--project-root\', type=Path, default=project_default)\n    parser.add_argument(\'--helper-dir\', type=Path, default=Path(os.environ.get(\'FIGURE2_HELPER_DIR\', str(script.parent))))\n    parser.add_argument(\'--input\', type=Path, default=Path(os.environ.get(\'FIGURE2_C_CURVES_CSV\', str(project_default / \'reports/figure2_c_curve_20261103T151741/nature_v5/PREDICTION_CURVES.csv\'))))\n    parser.add_argument(\'--output\', type=Path, default=Path(os.environ.get(\'FIGURE2_C_OUTPUT_ROOT\', str(project_default / \'reports/figure2_typography_20261105_v1/C\'))))\n    args = parser.parse_args()\n    project, helpers, source, output = [p.expanduser().resolve() for p in (args.project_root, args.helper_dir, args.input, args.output)]\n    print(f\'Project root: {project}\\nInput aggregate: {source}\\nHelper root: {helpers}\\nOutput root: {output}\')\n    helper = helpers / \'frozen_curve_helpers.py\'\n    if not helper.is_file() or not source.is_file():\n        raise FileNotFoundError(f\'Required aggregate/helper missing: {source}; {helper}\')\n    if output.exists():\n        raise FileExistsError(f\'New isolated output directory required: {output}\')\n    spec = importlib.util.spec_from_file_location(\'frozen_curve_helpers\', helper)\n    base = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(base)\n    base.FONT_PT = 7.0\n    fields, rows, groups, source_hash, common_n, checks = base.read_curves(source)\n    fonts = base.configure_fonts()\n    fig = base.plt.figure(figsize=base.MAIN_SIZE, facecolor=\'white\')\n    ax = fig.add_axes((0.225, 0.295, 0.735, 0.64))\n    base.style_axes(ax, base.y_limits(groups))\n    ax.grid(False, axis=\'both\', which=\'both\')\n    lines = base.draw_curves(ax, groups)\n    ax.set_xlabel(\'Relative birth position\', fontsize=7.5, fontweight=\'bold\', labelpad=5)\n    ax.set_ylabel(\'Standardized probability difference\', fontsize=7.5, fontweight=\'bold\', labelpad=5)\n    handles = [base.Line2D([], [], label=base.NAMES[c], linewidth=1.0, **base.STYLES[c]) for c in (\'DM\', \'T2D\')]\n    handles.append(base.Line2D([], [], color=\'none\', linestyle=\'none\', label=\' \'))\n    handles.extend((base.Line2D([], [], label=base.NAMES[c], linewidth=1.0, **base.STYLES[c]) for c in (\'IHD\', \'Angina\', \'CAS\')))\n    legend = fig.legend(handles=handles, loc=\'lower center\', bbox_to_anchor=(0.525, 0.025), ncol=2, fontsize=7.0, frameon=False, borderpad=0, handlelength=2.0, handletextpad=0.45, columnspacing=0.9, labelspacing=0.5)\n    visible = base.check_text_bounds(fig)\n    renderer = fig.canvas.get_renderer()\n    legend_box = legend.get_window_extent(renderer)\n    xlabel_box = ax.xaxis.label.get_window_extent(renderer)\n    if legend_box.y1 + 6 > xlabel_box.y0:\n        raise AssertionError(f\'Legend/x-label gap too small: {xlabel_box.y0 - legend_box.y1}\')\n    legend_texts = [t for t in legend.get_texts() if t.get_text().strip()]\n    for i, a in enumerate(legend_texts):\n        for b in legend_texts[i + 1:]:\n            if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):\n                raise AssertionError(\'Legend text overlaps\')\n    if any((line.get_visible() for line in ax.get_xgridlines() + ax.get_ygridlines())):\n        raise AssertionError(\'Only the saved zero reference line should remain\')\n    if fig.texts or ax.get_title():\n        raise AssertionError(\'No titles, subtitles, panel letters, or notes permitted\')\n    coordinates = base.check_lines(lines, groups)\n    prior = project / \'reports/figure2_c_curve_20261103T151741/nature_v5/Panel_C_plot_audit.json\'\n    if not prior.is_file():\n        raise FileNotFoundError(f\'Frozen curve audit missing: {prior}\')\n    prior_hash = base.sha256(prior)\n    previous = json.loads(prior.read_text())\n    expected_coordinates = previous[\'main\'][\'QA\'][\'disease_line_coordinates\']\n    if coordinates != expected_coordinates:\n        raise AssertionError(\'Curve coordinates differ from the original Nature panel\')\n    output.mkdir(parents=True, exist_ok=False)\n    saved = base.save_formats(fig, output, \'Panel_C\')\n    base.plt.close(fig)\n    with (output / \'Panel_C_plot_data.csv\').open(\'x\', encoding=\'utf-8\', newline=\'\') as f:\n        writer = csv.DictWriter(f, fieldnames=fields)\n        writer.writeheader()\n        writer.writerows(rows)\n    with (output / \'Panel_C_plot_data.csv\').open(newline=\'\', encoding=\'utf-8\') as f:\n        if list(csv.DictReader(f)) != rows:\n            raise AssertionError(\'Saved source cells changed\')\n    if base.sha256(source) != source_hash or base.sha256(prior) != prior_hash:\n        raise AssertionError(\'An original input was modified\')\n    (output / \'code\').mkdir()\n    for path in (script, helper):\n        with (output / \'code\' / path.name).open(\'xb\') as f:\n            f.write(path.read_bytes())\n    audit = {\'created_utc\': datetime.now(timezone.utc).isoformat(), \'status\': \'PASS\', \'scope\': \'Typography only; exact saved five-disease curves and intervals\', \'source\': str(source), \'source_sha256\': source_hash, \'source_unchanged\': True, \'previous_audit\': str(prior), \'previous_audit_sha256\': prior_hash, \'common_N_from_source\': common_n, \'rows\': len(rows), \'models_fitted\': 0, \'predictions_computed\': 0, \'intervals_computed\': 0, \'font_family\': \'Arial\', \'fonts\': fonts, \'old_font_pt\': 6.1, \'body_and_legend_font_pt\': 7.0, \'axis_title_font_pt\': 7.5, \'axis_title_weight\': \'bold\', \'figure_inches\': list(base.MAIN_SIZE), \'plot_bounds\': list(ax.get_position().bounds), \'legend_columns\': {\'left\': [base.NAMES[c] for c in (\'DM\', \'T2D\')], \'right\': [base.NAMES[c] for c in (\'IHD\', \'Angina\', \'CAS\')]}, \'legend_xlabel_gap_pixels\': xlabel_box.y0 - legend_box.y1, \'no_clipped_or_overlapping_labels\': True, \'no_new_title_or_notes\': True, \'source_validation\': checks, \'exact_original_coordinates\': coordinates, \'CI_alpha\': 0.085, \'styles\': base.STYLES, \'outputs\': saved, \'visible_text\': visible, \'script_sha256\': base.sha256(script), \'helper_sha256\': base.sha256(helper), \'python\': sys.executable}\n    with (output / \'Panel_C_audit.json\').open(\'x\', encoding=\'utf-8\') as f:\n        json.dump(audit, f, ensure_ascii=False, indent=2)\n        f.write(\'\\n\')\n    print(json.dumps({\'status\': \'PASS\', \'output\': str(output), \'body_font_pt\': 7.0, \'bold_axis_title_pt\': 7.5, \'source_rows\': len(rows)}))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
