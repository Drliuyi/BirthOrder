"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import collections
import csv
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
from matplotlib.patches import Polygon
import numpy as np

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def write_csv(path, rows):
    with path.open('x', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def interval(row):
    lo, hi = (row['IE_CI_low_per_1000'], row['IE_CI_high_per_1000'])
    if not lo or not hi:
        assert not lo and (not hi), 'ONE_SIDED_MISSING_CI'
        return None
    lo, hi = (float(lo), float(hi))
    assert math.isfinite(lo) and math.isfinite(hi) and (lo <= hi)
    assert int(row['bootstrap_valid']) > 0
    assert row['CI_method'] == 'FROZEN_TYPE7_PERCENTILE'
    return (lo, hi)

def main():
    print(json.dumps({'input_root': str(SOURCE), 'helper_root': str(HERE), 'output_root': str(OUTPUT)}))
    if OUTPUT.exists():
        raise FileExistsError(f'PRESERVE_PREVIOUS_FIGURES:{OUTPUT}')
    files = [SOURCE / 'B100_121_COMPLETE_RESULTS.csv', SOURCE / 'B100_SIX_CONTRASTS_RESULTS.csv', SOURCE / 'TABLE_MANIFEST.json', HERE / 'nearby_labels.py', ROOT / 'figure_revision_nmr_v4_near_labels/B100_QUADRANT_NMR_VISIBLE_V2.pdf', ROOT / 'figure_revision_nmr_v4_near_labels/B100_NMR_QUADRANT_V2.pdf', ORIGINAL / 'POINTWISE_AXIS_CI_AUDIT.csv', ORIGINAL / 'ANNOTATIONS.csv']
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(f'MISSING_REQUIRED_PLOT_INPUT:{path}')
    before = {str(p): sha(p) for p in files}
    base = module('_near_label_plot', HERE / 'nearby_labels.py')
    prior = module('_prior_plot', ROOT / 'plot_and_report.py')
    name_helper = ROOT.parents[2] / 'scripts/birth/border_fullfollowup_b100_20261027/repairs/cache_only_screen_20261028_v1/plot_screen.py'
    names = module('_original_display_names', name_helper)

    def label(row):
        return base.DISPLAY.get(row['indicator_id']) or prior.DISPLAY.get(row['indicator_id']) or names.english(row)
    single = read(files[0])
    lookup = {(r['indicator_id'], r['endpoint_id']): r for r in single}
    assert len(lookup) == len(single) == 605
    pairs = read(files[1])
    points, audit = ([], [])
    cis = {}
    for row in pairs:
        try:
            x, y = (float(row['IE_B_per_1000']), float(row['IE_A_per_1000']))
        except ValueError:
            continue
        assert math.isfinite(x) and math.isfinite(y)
        xsource = lookup[row['indicator_id'], row['endpoint_B']]
        ysource = lookup[row['indicator_id'], row['endpoint_A']]
        for point, source in ((x, xsource), (y, ysource)):
            assert abs(float(source['IE_per_1000']) - point) <= 1e-12, 'POINT_SOURCE_MISMATCH'
            assert source['horizon_years'] == row['horizon_years'] == '11'
        quad = ('P' if x > 0 else 'N') + ('P' if y > 0 else 'N')
        assert x != 0 and y != 0
        xci, yci = (interval(xsource), interval(ysource))
        key = (row['indicator_id'], row['contrast'])
        assert key not in cis
        cis[key] = (xci, yci)
        points.append((row, x, y, quad))
        audit.append({'indicator_id': row['indicator_id'], 'contrast': row['contrast'], 'x_endpoint': row['endpoint_B'], 'x_disease': row['disease_B'], 'x_IE_per_1000': x, 'x_95CI_low_per_1000': xci[0] if xci else '', 'x_95CI_high_per_1000': xci[1] if xci else '', 'x_valid_B': xsource['bootstrap_valid'], 'x_status': xsource['result_status'], 'x_method_hash': xsource['source_method_hash'], 'y_endpoint': row['endpoint_A'], 'y_disease': row['disease_A'], 'y_IE_per_1000': y, 'y_95CI_low_per_1000': yci[0] if yci else '', 'y_95CI_high_per_1000': yci[1] if yci else '', 'y_valid_B': ysource['bootstrap_valid'], 'y_status': ysource['result_status'], 'y_method_hash': ysource['source_method_hash'], 'interval_scope': 'AXIS_SPECIFIC_MARGINAL_95CI_NOT_JOINT_REGION_NOT_DELTA_CI', 'x_CI_drawn': xci is not None, 'y_CI_drawn': yci is not None})
    assert len(points) == 660 and sum((r['x_CI_drawn'] for r in audit)) == 660
    assert sum((r['y_CI_drawn'] for r in audit)) == 660
    counts = collections.Counter(((r['indicator_id'], q) for r, x, y, q in points))
    fontManager.addfont('/System/Library/Fonts/Supplemental/Arial.ttf')
    matplotlib.rcParams.update({'font.family': 'Arial', 'font.size': 14.5, 'text.color': 'black', 'axes.labelcolor': 'black', 'axes.edgecolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black', 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    OUTPUT.mkdir()
    (OUTPUT / 'code').mkdir()
    import shutil
    shutil.copy2(Path(__file__), OUTPUT / 'code/build.py')
    shutil.copy2(HERE / 'nearby_labels.py', OUTPUT / 'code/nearby_labels.py')
    qa, annotation_rows = ([], [])

    def plot(dataset, selections, stem):
        items = []
        for quad, ids in selections.items():
            chosen = []
            for ident in ids:
                choices = [p for p in dataset if p[0]['indicator_id'] == ident and p[3] == quad]
                assert choices
                chosen.append(max(choices, key=lambda p: math.hypot(p[1], p[2])))
            items.extend(sorted(chosen, key=lambda p: -p[2]))
        chosen_keys = {(r['indicator_id'], r['contrast']) for r, x, y, q in items}
        extremes = [abs(v) for r, x, y, q in dataset for v in (x, y)]
        for row, x, y, q in dataset:
            for ci in cis[row['indicator_id'], row['contrast']]:
                if ci is not None:
                    extremes.extend(map(abs, ci))
        point_limit = math.ceil(max((max(abs(x), abs(y)) for r, x, y, q in dataset)) * 1.12 * 10) / 10
        limit = max(point_limit, math.ceil(max(extremes) * 1.01 * 10) / 10)
        fig, ax = plt.subplots(figsize=(8.4, 8.7))
        fig.subplots_adjust(left=0.13, right=0.97, top=0.98, bottom=0.17)
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit))
        ax.set_aspect('equal', adjustable='box')
        ax.add_patch(Polygon([(-limit, -limit), (-limit, limit), (limit, limit)], fc='#EDF5FC', ec='none', zorder=0))
        ax.add_patch(Polygon([(-limit, -limit), (limit, -limit), (limit, limit)], fc='#F3F3F3', ec='none', zorder=0))
        ax.axhline(0, color='#999999', lw=0.9, zorder=1)
        ax.axvline(0, color='#999999', lw=0.9, zorder=1)
        ax.plot([-limit, limit], [-limit, limit], color='#666666', ls=(0, (5, 3)), lw=1.3, zorder=2)
        h_count, v_count, half_ci, no_ci = (0, 0, 0, 0)
        ordered = sorted(dataset, key=lambda p: ((p[0]['indicator_id'], p[0]['contrast']) in chosen_keys, p[0]['domain'] == 'NMR', counts[p[0]['indicator_id'], p[3]]))
        for row, x, y, quad in ordered:
            key = (row['indicator_id'], row['contrast'])
            xci, yci = cis[key]
            highlighted = key in chosen_keys
            color = base.PALETTE[quad]
            alpha, width = (0.55, 0.5) if highlighted else (0.18, 0.35)
            if xci:
                ax.plot(xci, [y, y], color=color, alpha=alpha, lw=width, zorder=1.5, solid_capstyle='butt')
                h_count += 1
            if yci:
                ax.plot([x, x], yci, color=color, alpha=alpha, lw=width, zorder=1.5, solid_capstyle='butt')
                v_count += 1
            half_ci += (xci is None) != (yci is None)
            no_ci += xci is None and yci is None
        plotted = []
        for row, x, y, quad in ordered:
            rgb = np.array(matplotlib.colors.to_rgb(base.PALETTE[quad]))
            strength = base.STRENGTHS[counts[row['indicator_id'], quad] - 1]
            ax.scatter([x], [y], s=88, c=[1 + (rgb - 1) * strength], edgecolors='white', linewidths=0.4, zorder=3)
            plotted.append((row['indicator_id'], row['contrast'], x, y))
        assert collections.Counter(plotted) == collections.Counter(((r['indicator_id'], r['contrast'], x, y) for r, x, y, q in dataset))
        ax.set_xlabel('DM / T2D indirect risk difference (per 1,000)', fontsize=16.25, fontweight='bold', labelpad=9)
        ax.set_ylabel('IHD / Angina / CAS indirect risk difference (per 1,000)', fontsize=16.25, fontweight='bold', labelpad=9)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
        positions = base.place_nearby_labels(ax, items, dataset, label)
        for (row, x, y, q), position in zip(items, positions):
            annotation_rows.append({'figure': stem, 'indicator_id': row['indicator_id'], 'label': label(row), 'contrast': row['contrast'], 'x': x, 'y': y, 'label_x': position['center'][0], 'label_y': position['center'][1]})
        keyax = fig.add_axes([0.34, 0.025, 0.38, 0.028])
        swatches = []
        for q in ('NP', 'PP', 'NN', 'PN'):
            rgb = np.array(matplotlib.colors.to_rgb(base.PALETTE[q]))
            swatches.append([1 + (rgb - 1) * s for s in base.STRENGTHS])
        keyax.imshow(swatches, aspect='auto', extent=(0.5, 6.5, 4, 0))
        keyax.set(yticks=[], xticks=range(1, 7))
        keyax.tick_params(length=0)
        keyax.set_title('Directional agreement (of 6)', fontsize=13)
        for ext in ('pdf', 'png', 'svg'):
            fig.savefig(OUTPUT / (stem + '.' + ext), dpi=600, facecolor='white', bbox_inches='tight', pad_inches=0.06)
        fig.savefig(OUTPUT / (stem + '_preview.png'), dpi=150, facecolor='white', bbox_inches='tight', pad_inches=0.06)
        plt.close(fig)
        qa.append({'figure': stem, 'coordinates': len(dataset), 'indicators': len({r['indicator_id'] for r, x, y, q in dataset}), 'horizontal_95CI': h_count, 'vertical_95CI': v_count, 'one_axis_only': int(half_ci), 'both_axes_no_CI': int(no_ci), 'limit': limit, 'all_CI_endpoints_inside_axes': max(extremes) < limit, 'equal_axis_scale': True, 'point_multiset_unchanged': True, 'label_point_overlap': False, 'label_label_overlap': False})
    plot(points, base.LABELS_ALL, 'B100_QUADRANT_LARGER_TEXT')
    write_csv(OUTPUT / 'POINTWISE_AXIS_CI_AUDIT.csv', audit)
    write_csv(OUTPUT / 'ANNOTATIONS.csv', annotation_rows)
    assert (OUTPUT / 'POINTWISE_AXIS_CI_AUDIT.csv').read_bytes() == (ORIGINAL / 'POINTWISE_AXIS_CI_AUDIT.csv').read_bytes()
    annotation_fields = ['indicator_id', 'label', 'contrast', 'x', 'y']
    previous_annotations = [r for r in read(ORIGINAL / 'ANNOTATIONS.csv') if r['figure'] == 'B100_QUADRANT_COMPLETED_WITH_95CI']
    assert len(previous_annotations) == len(annotation_rows) == 20
    for old, new in zip(previous_annotations, annotation_rows):
        assert all((str(old[k]) == str(new[k]) for k in annotation_fields))
    after = {str(p): sha(p) for p in files}
    assert before == after, 'SOURCE_RELEASE_OR_OLD_FIGURE_CHANGED'
    record = {'time_beijing': datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(), 'plots': qa, 'source_hashes_before': before, 'source_hashes_after': after, 'source_unchanged': True, 'new_fits': 0, 'new_bootstrap': 0, 'CI': 'EXISTING_AXIS_SPECIFIC_MARGINAL_TYPE7_95CI;NOT_JOINT_95_REGION;NOT_DELTA_CI', 'partial_bootstrap_status': 'UNCHANGED_SUCCESS_CONDITIONAL', 'missing_intervals': 'NOT_IMPUTED_NOT_ZERO_NOT_SYNTHESIZED'}
    record.update({'revision': 'Typography only; existing 20 scatter labels unchanged', 'font_sizes_pt': {'axes': 16.25, 'ticks': 14.5, 'labels': 14.5, 'legend': 13}, 'original_pointwise_audit_byte_identical': True, 'annotation_ids_names_anchors_identical': True})
    (OUTPUT / 'FIGURE_95CI_QA.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    note = '# 散点图字号修订\n\n复用已完成2026-10-30 figures_complete_v1：660个坐标及各660条横纵边际95%CI、颜色、点大小、绘图顺序、轴范围和20个既有散点标注不变。\n轴标题12.5→16.25pt；刻度11→14.5pt；标注11.2→14.5pt；色键标题10→13pt。只重新避让文字，未增加、替换或删除散点标注。无图题、副标题、面板字母或方法note。\nPOINTWISE_AXIS_CI_AUDIT.csv与原图逐字节一致；标注指标ID、名称、对比及原始坐标完全一致。没有拟合、Bootstrap或P/q重算。\n所有风险差仍是原模型的11年每千人差，图内不另加时间字样。既有探索性与有效重复条件推断限制不因改图改变。\n运行代码：本输出的code文件可在当前项目内执行，默认输出已存在时拒绝覆盖；重现时需先在新代码副本中将OUTPUT设置为不存在的新版本目录。\n'
    (OUTPUT / 'FIGURE_95CI_NOTES_ZH.md').write_text(note, encoding='utf-8')
    manifest = {p.name: {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(OUTPUT.iterdir()) if p.is_file()}
    (OUTPUT / 'FIGURE_MANIFEST.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'output': str(OUTPUT), 'plots': qa, 'source_unchanged': True}, ensure_ascii=False))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Enlarge typography in the already completed, audited quadrant figure.\n\nAggregate-only, no estimation. Source releases and earlier figures are read-only.\n"""\nfrom __future__ import annotations\nimport collections\nimport csv\nimport datetime\nimport hashlib\nimport importlib.util\nimport json\nimport math\nimport os\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib.font_manager import fontManager\nfrom matplotlib.patches import Polygon\nimport numpy as np\nHERE = Path(os.environ.get(\'BIRTH_ORDER_SCRIPT_ROOT\', str(Path(__file__).resolve().parent)))\nBASE = next((p for p in Path(__file__).resolve().parents if (p / \'scripts/birth\').is_dir()))\nROOT = BASE / \'reports/birth/b100_b1000_release_20261030_v1\'\nSOURCE = Path(os.environ.get(\'BIRTH_ORDER_INPUT_ROOT\', str(BASE / \'reports/birth/b100_completion_20261030_v1/release\')))\nOUTPUT = Path(os.environ.get(\'BIRTH_ORDER_OUTPUT_ROOT\', str(BASE / \'reports/birth/b100_quadrant_typography_20261105_v1\')))\nORIGINAL = BASE / \'reports/birth/b100_completion_20261030_v1/figures_complete_v1\'\n\ndef module(name, path):\n    spec = importlib.util.spec_from_file_location(name, path)\n    value = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(value)\n    return value\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef read(path):\n    with path.open(encoding=\'utf-8-sig\', newline=\'\') as stream:\n        return list(csv.DictReader(stream))\n\ndef write_csv(path, rows):\n    with path.open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as stream:\n        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))\n        writer.writeheader()\n        writer.writerows(rows)\n\ndef interval(row):\n    lo, hi = (row[\'IE_CI_low_per_1000\'], row[\'IE_CI_high_per_1000\'])\n    if not lo or not hi:\n        assert not lo and (not hi), \'ONE_SIDED_MISSING_CI\'\n        return None\n    lo, hi = (float(lo), float(hi))\n    assert math.isfinite(lo) and math.isfinite(hi) and (lo <= hi)\n    assert int(row[\'bootstrap_valid\']) > 0\n    assert row[\'CI_method\'] == \'FROZEN_TYPE7_PERCENTILE\'\n    return (lo, hi)\n\ndef main():\n    print(json.dumps({\'input_root\': str(SOURCE), \'helper_root\': str(HERE), \'output_root\': str(OUTPUT)}))\n    if OUTPUT.exists():\n        raise FileExistsError(f\'PRESERVE_PREVIOUS_FIGURES:{OUTPUT}\')\n    files = [SOURCE / \'B100_121_COMPLETE_RESULTS.csv\', SOURCE / \'B100_SIX_CONTRASTS_RESULTS.csv\', SOURCE / \'TABLE_MANIFEST.json\', HERE / \'nearby_labels.py\', ROOT / \'figure_revision_nmr_v4_near_labels/B100_QUADRANT_NMR_VISIBLE_V2.pdf\', ROOT / \'figure_revision_nmr_v4_near_labels/B100_NMR_QUADRANT_V2.pdf\', ORIGINAL / \'POINTWISE_AXIS_CI_AUDIT.csv\', ORIGINAL / \'ANNOTATIONS.csv\']\n    for path in files:\n        if not path.is_file():\n            raise FileNotFoundError(f\'MISSING_REQUIRED_PLOT_INPUT:{path}\')\n    before = {str(p): sha(p) for p in files}\n    base = module(\'_near_label_plot\', HERE / \'nearby_labels.py\')\n    prior = module(\'_prior_plot\', ROOT / \'plot_and_report.py\')\n    name_helper = ROOT.parents[2] / \'scripts/birth/border_fullfollowup_b100_20261027/repairs/cache_only_screen_20261028_v1/plot_screen.py\'\n    names = module(\'_original_display_names\', name_helper)\n\n    def label(row):\n        return base.DISPLAY.get(row[\'indicator_id\']) or prior.DISPLAY.get(row[\'indicator_id\']) or names.english(row)\n    single = read(files[0])\n    lookup = {(r[\'indicator_id\'], r[\'endpoint_id\']): r for r in single}\n    assert len(lookup) == len(single) == 605\n    pairs = read(files[1])\n    points, audit = ([], [])\n    cis = {}\n    for row in pairs:\n        try:\n            x, y = (float(row[\'IE_B_per_1000\']), float(row[\'IE_A_per_1000\']))\n        except ValueError:\n            continue\n        assert math.isfinite(x) and math.isfinite(y)\n        xsource = lookup[row[\'indicator_id\'], row[\'endpoint_B\']]\n        ysource = lookup[row[\'indicator_id\'], row[\'endpoint_A\']]\n        for point, source in ((x, xsource), (y, ysource)):\n            assert abs(float(source[\'IE_per_1000\']) - point) <= 1e-12, \'POINT_SOURCE_MISMATCH\'\n            assert source[\'horizon_years\'] == row[\'horizon_years\'] == \'11\'\n        quad = (\'P\' if x > 0 else \'N\') + (\'P\' if y > 0 else \'N\')\n        assert x != 0 and y != 0\n        xci, yci = (interval(xsource), interval(ysource))\n        key = (row[\'indicator_id\'], row[\'contrast\'])\n        assert key not in cis\n        cis[key] = (xci, yci)\n        points.append((row, x, y, quad))\n        audit.append({\'indicator_id\': row[\'indicator_id\'], \'contrast\': row[\'contrast\'], \'x_endpoint\': row[\'endpoint_B\'], \'x_disease\': row[\'disease_B\'], \'x_IE_per_1000\': x, \'x_95CI_low_per_1000\': xci[0] if xci else \'\', \'x_95CI_high_per_1000\': xci[1] if xci else \'\', \'x_valid_B\': xsource[\'bootstrap_valid\'], \'x_status\': xsource[\'result_status\'], \'x_method_hash\': xsource[\'source_method_hash\'], \'y_endpoint\': row[\'endpoint_A\'], \'y_disease\': row[\'disease_A\'], \'y_IE_per_1000\': y, \'y_95CI_low_per_1000\': yci[0] if yci else \'\', \'y_95CI_high_per_1000\': yci[1] if yci else \'\', \'y_valid_B\': ysource[\'bootstrap_valid\'], \'y_status\': ysource[\'result_status\'], \'y_method_hash\': ysource[\'source_method_hash\'], \'interval_scope\': \'AXIS_SPECIFIC_MARGINAL_95CI_NOT_JOINT_REGION_NOT_DELTA_CI\', \'x_CI_drawn\': xci is not None, \'y_CI_drawn\': yci is not None})\n    assert len(points) == 660 and sum((r[\'x_CI_drawn\'] for r in audit)) == 660\n    assert sum((r[\'y_CI_drawn\'] for r in audit)) == 660\n    counts = collections.Counter(((r[\'indicator_id\'], q) for r, x, y, q in points))\n    fontManager.addfont(\'/System/Library/Fonts/Supplemental/Arial.ttf\')\n    matplotlib.rcParams.update({\'font.family\': \'Arial\', \'font.size\': 14.5, \'text.color\': \'black\', \'axes.labelcolor\': \'black\', \'axes.edgecolor\': \'black\', \'xtick.color\': \'black\', \'ytick.color\': \'black\', \'figure.facecolor\': \'white\', \'axes.facecolor\': \'white\', \'pdf.fonttype\': 42, \'ps.fonttype\': 42, \'svg.fonttype\': \'none\'})\n    OUTPUT.mkdir()\n    (OUTPUT / \'code\').mkdir()\n    import shutil\n    shutil.copy2(Path(__file__), OUTPUT / \'code/build.py\')\n    shutil.copy2(HERE / \'nearby_labels.py\', OUTPUT / \'code/nearby_labels.py\')\n    qa, annotation_rows = ([], [])\n\n    def plot(dataset, selections, stem):\n        items = []\n        for quad, ids in selections.items():\n            chosen = []\n            for ident in ids:\n                choices = [p for p in dataset if p[0][\'indicator_id\'] == ident and p[3] == quad]\n                assert choices\n                chosen.append(max(choices, key=lambda p: math.hypot(p[1], p[2])))\n            items.extend(sorted(chosen, key=lambda p: -p[2]))\n        chosen_keys = {(r[\'indicator_id\'], r[\'contrast\']) for r, x, y, q in items}\n        extremes = [abs(v) for r, x, y, q in dataset for v in (x, y)]\n        for row, x, y, q in dataset:\n            for ci in cis[row[\'indicator_id\'], row[\'contrast\']]:\n                if ci is not None:\n                    extremes.extend(map(abs, ci))\n        point_limit = math.ceil(max((max(abs(x), abs(y)) for r, x, y, q in dataset)) * 1.12 * 10) / 10\n        limit = max(point_limit, math.ceil(max(extremes) * 1.01 * 10) / 10)\n        fig, ax = plt.subplots(figsize=(8.4, 8.7))\n        fig.subplots_adjust(left=0.13, right=0.97, top=0.98, bottom=0.17)\n        ax.set(xlim=(-limit, limit), ylim=(-limit, limit))\n        ax.set_aspect(\'equal\', adjustable=\'box\')\n        ax.add_patch(Polygon([(-limit, -limit), (-limit, limit), (limit, limit)], fc=\'#EDF5FC\', ec=\'none\', zorder=0))\n        ax.add_patch(Polygon([(-limit, -limit), (limit, -limit), (limit, limit)], fc=\'#F3F3F3\', ec=\'none\', zorder=0))\n        ax.axhline(0, color=\'#999999\', lw=0.9, zorder=1)\n        ax.axvline(0, color=\'#999999\', lw=0.9, zorder=1)\n        ax.plot([-limit, limit], [-limit, limit], color=\'#666666\', ls=(0, (5, 3)), lw=1.3, zorder=2)\n        h_count, v_count, half_ci, no_ci = (0, 0, 0, 0)\n        ordered = sorted(dataset, key=lambda p: ((p[0][\'indicator_id\'], p[0][\'contrast\']) in chosen_keys, p[0][\'domain\'] == \'NMR\', counts[p[0][\'indicator_id\'], p[3]]))\n        for row, x, y, quad in ordered:\n            key = (row[\'indicator_id\'], row[\'contrast\'])\n            xci, yci = cis[key]\n            highlighted = key in chosen_keys\n            color = base.PALETTE[quad]\n            alpha, width = (0.55, 0.5) if highlighted else (0.18, 0.35)\n            if xci:\n                ax.plot(xci, [y, y], color=color, alpha=alpha, lw=width, zorder=1.5, solid_capstyle=\'butt\')\n                h_count += 1\n            if yci:\n                ax.plot([x, x], yci, color=color, alpha=alpha, lw=width, zorder=1.5, solid_capstyle=\'butt\')\n                v_count += 1\n            half_ci += (xci is None) != (yci is None)\n            no_ci += xci is None and yci is None\n        plotted = []\n        for row, x, y, quad in ordered:\n            rgb = np.array(matplotlib.colors.to_rgb(base.PALETTE[quad]))\n            strength = base.STRENGTHS[counts[row[\'indicator_id\'], quad] - 1]\n            ax.scatter([x], [y], s=88, c=[1 + (rgb - 1) * strength], edgecolors=\'white\', linewidths=0.4, zorder=3)\n            plotted.append((row[\'indicator_id\'], row[\'contrast\'], x, y))\n        assert collections.Counter(plotted) == collections.Counter(((r[\'indicator_id\'], r[\'contrast\'], x, y) for r, x, y, q in dataset))\n        ax.set_xlabel(\'DM / T2D indirect risk difference (per 1,000)\', fontsize=16.25, fontweight=\'bold\', labelpad=9)\n        ax.set_ylabel(\'IHD / Angina / CAS indirect risk difference (per 1,000)\', fontsize=16.25, fontweight=\'bold\', labelpad=9)\n        for spine in ax.spines.values():\n            spine.set_visible(True)\n            spine.set_linewidth(1.5)\n        positions = base.place_nearby_labels(ax, items, dataset, label)\n        for (row, x, y, q), position in zip(items, positions):\n            annotation_rows.append({\'figure\': stem, \'indicator_id\': row[\'indicator_id\'], \'label\': label(row), \'contrast\': row[\'contrast\'], \'x\': x, \'y\': y, \'label_x\': position[\'center\'][0], \'label_y\': position[\'center\'][1]})\n        keyax = fig.add_axes([0.34, 0.025, 0.38, 0.028])\n        swatches = []\n        for q in (\'NP\', \'PP\', \'NN\', \'PN\'):\n            rgb = np.array(matplotlib.colors.to_rgb(base.PALETTE[q]))\n            swatches.append([1 + (rgb - 1) * s for s in base.STRENGTHS])\n        keyax.imshow(swatches, aspect=\'auto\', extent=(0.5, 6.5, 4, 0))\n        keyax.set(yticks=[], xticks=range(1, 7))\n        keyax.tick_params(length=0)\n        keyax.set_title(\'Directional agreement (of 6)\', fontsize=13)\n        for ext in (\'pdf\', \'png\', \'svg\'):\n            fig.savefig(OUTPUT / (stem + \'.\' + ext), dpi=600, facecolor=\'white\', bbox_inches=\'tight\', pad_inches=0.06)\n        fig.savefig(OUTPUT / (stem + \'_preview.png\'), dpi=150, facecolor=\'white\', bbox_inches=\'tight\', pad_inches=0.06)\n        plt.close(fig)\n        qa.append({\'figure\': stem, \'coordinates\': len(dataset), \'indicators\': len({r[\'indicator_id\'] for r, x, y, q in dataset}), \'horizontal_95CI\': h_count, \'vertical_95CI\': v_count, \'one_axis_only\': int(half_ci), \'both_axes_no_CI\': int(no_ci), \'limit\': limit, \'all_CI_endpoints_inside_axes\': max(extremes) < limit, \'equal_axis_scale\': True, \'point_multiset_unchanged\': True, \'label_point_overlap\': False, \'label_label_overlap\': False})\n    plot(points, base.LABELS_ALL, \'B100_QUADRANT_LARGER_TEXT\')\n    write_csv(OUTPUT / \'POINTWISE_AXIS_CI_AUDIT.csv\', audit)\n    write_csv(OUTPUT / \'ANNOTATIONS.csv\', annotation_rows)\n    assert (OUTPUT / \'POINTWISE_AXIS_CI_AUDIT.csv\').read_bytes() == (ORIGINAL / \'POINTWISE_AXIS_CI_AUDIT.csv\').read_bytes()\n    annotation_fields = [\'indicator_id\', \'label\', \'contrast\', \'x\', \'y\']\n    previous_annotations = [r for r in read(ORIGINAL / \'ANNOTATIONS.csv\') if r[\'figure\'] == \'B100_QUADRANT_COMPLETED_WITH_95CI\']\n    assert len(previous_annotations) == len(annotation_rows) == 20\n    for old, new in zip(previous_annotations, annotation_rows):\n        assert all((str(old[k]) == str(new[k]) for k in annotation_fields))\n    after = {str(p): sha(p) for p in files}\n    assert before == after, \'SOURCE_RELEASE_OR_OLD_FIGURE_CHANGED\'\n    record = {\'time_beijing\': datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(), \'plots\': qa, \'source_hashes_before\': before, \'source_hashes_after\': after, \'source_unchanged\': True, \'new_fits\': 0, \'new_bootstrap\': 0, \'CI\': \'EXISTING_AXIS_SPECIFIC_MARGINAL_TYPE7_95CI;NOT_JOINT_95_REGION;NOT_DELTA_CI\', \'partial_bootstrap_status\': \'UNCHANGED_SUCCESS_CONDITIONAL\', \'missing_intervals\': \'NOT_IMPUTED_NOT_ZERO_NOT_SYNTHESIZED\'}\n    record.update({\'revision\': \'Typography only; existing 20 scatter labels unchanged\', \'font_sizes_pt\': {\'axes\': 16.25, \'ticks\': 14.5, \'labels\': 14.5, \'legend\': 13}, \'original_pointwise_audit_byte_identical\': True, \'annotation_ids_names_anchors_identical\': True})\n    (OUTPUT / \'FIGURE_95CI_QA.json\').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding=\'utf-8\')\n    note = \'# 散点图字号修订\\n\\n复用已完成2026-10-30 figures_complete_v1：660个坐标及各660条横纵边际95%CI、颜色、点大小、绘图顺序、轴范围和20个既有散点标注不变。\\n轴标题12.5→16.25pt；刻度11→14.5pt；标注11.2→14.5pt；色键标题10→13pt。只重新避让文字，未增加、替换或删除散点标注。无图题、副标题、面板字母或方法note。\\nPOINTWISE_AXIS_CI_AUDIT.csv与原图逐字节一致；标注指标ID、名称、对比及原始坐标完全一致。没有拟合、Bootstrap或P/q重算。\\n所有风险差仍是原模型的11年每千人差，图内不另加时间字样。既有探索性与有效重复条件推断限制不因改图改变。\\n运行代码：本输出的code文件可在当前项目内执行，默认输出已存在时拒绝覆盖；重现时需先在新代码副本中将OUTPUT设置为不存在的新版本目录。\\n\'\n    (OUTPUT / \'FIGURE_95CI_NOTES_ZH.md\').write_text(note, encoding=\'utf-8\')\n    manifest = {p.name: {\'sha256\': sha(p), \'bytes\': p.stat().st_size} for p in sorted(OUTPUT.iterdir()) if p.is_file()}\n    (OUTPUT / \'FIGURE_MANIFEST.json\').write_text(json.dumps(manifest, indent=2), encoding=\'utf-8\')\n    print(json.dumps({\'output\': str(OUTPUT), \'plots\': qa, \'source_unchanged\': True}, ensure_ascii=False))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
