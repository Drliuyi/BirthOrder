"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
ORDER = ('DM', 'T2D', 'IHD', 'Angina', 'CAS')
WRAPPED = {'DM': 'Diabetes\nmellitus', 'T2D': 'Type 2 diabetes', 'IHD': 'Ischemic heart\ndisease', 'Angina': 'Angina pectoris', 'CAS': 'Coronary\natherosclerosis'}
SIZE = (3.9, 3.3)
BOUNDS = ((0.015, 0.16, 0.305, 0.79), (0.345, 0.16, 0.225, 0.79), (0.59, 0.16, 0.245, 0.79), (0.855, 0.16, 0.13, 0.79))
FONTS = {'default': 8.05, 'header': 7.8, 'disease': 8.3, 'model': 7.5, 'N': 7.1, 'effect': 7.35, 'q': 7.1, 'ticks': 7.5, 'xlabel': 7.6}
ORIGINAL_FONTS = {'default': 7, 'header': 6.8, 'disease': 7.2, 'model': 6.5, 'N': 6.2, 'effect': 6.4, 'q': 6.2, 'ticks': 6.5, 'xlabel': 6.6}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def readcsv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def qtext(q):
    q = float(q)
    e = math.floor(math.log10(q))
    s = f'{q / 10 ** e:.2f}' + f'\\times10^{{{e}}}'
    return '$' + ('\\mathbf{' + s + '}' if q < 0.05 else s) + '$'

def settings():
    for weight in ('normal', 'bold'):
        font_manager.findfont(font_manager.FontProperties(family='Arial', weight=weight), fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': FONTS['default'], 'axes.unicode_minus': False, 'savefig.facecolor': 'white', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular', 'mathtext.fallback': None, 'axes.linewidth': 0.6, 'text.color': '#222222', 'axes.labelcolor': '#222222', 'xtick.color': '#333333', 'ytick.color': '#333333'})

def forest_axes(measure, second):
    fig = plt.figure(figsize=SIZE, facecolor='white')
    axes = [fig.add_axes(p) for p in BOUNDS]
    left, forest, estimate, qax = axes
    for ax in axes:
        ax.set_ylim(0.2, 17.8)
        ax.set_yticks([])
    for ax in (left, estimate, qax):
        ax.set_xlim(0, 1)
        ax.axis('off')
    for label, x, align in [('Disease', 0.015, 'left'), (second, 0.015 + 0.305 * 0.86, 'center'), (measure + ' (95% CI)', 0.59, 'left'), ('BH q', 0.855, 'left')]:
        fig.text(x, 0.965, label, fontsize=FONTS['header'], fontweight='bold', ha=align)
    for i, d in enumerate(ORDER):
        left.text(0, 16 - 3.5 * i, WRAPPED[d], va='center', fontsize=FONTS['disease'], linespacing=1.05, color='black')
        if i < 4:
            for ax in axes:
                ax.axhline(14.25 - 3.5 * i, lw=0.45, color='#E6E6E6', zorder=0)
    forest.set_xscale('log')
    forest.set_xlim(0.8, 1.3)
    forest.xaxis.set_major_locator(ticker.FixedLocator([0.8, 1, 1.2]))
    forest.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:g}'))
    forest.xaxis.set_minor_locator(ticker.NullLocator())
    forest.tick_params(axis='x', length=2.5, width=0.6, pad=3, labelsize=FONTS['ticks'])
    forest.spines[['top', 'left', 'right']].set_visible(False)
    forest.spines['bottom'].set_color('#777777')
    forest.axvline(1, color='#888888', ls=(0, (3, 3)), lw=0.65, zorder=0)
    forest.set_xlabel('Relative birth order (0→1)\n' + ('Odds ratio' if measure == 'OR' else 'Hazard ratio'), fontsize=FONTS['xlabel'], labelpad=4)
    return (fig, axes)

def interval(ax, x, lo, hi, y, color, alpha=1, size=3.8):
    ax.plot([lo, hi], [y, y], color=color, lw=0.85, alpha=alpha, zorder=2)
    for v in (lo, hi):
        ax.plot([v, v], [y - 0.11, y + 0.11], color=color, lw=0.7, alpha=alpha, zorder=2)
    ax.plot(x, y, 'o', ms=size, color=color, alpha=alpha, mec='white', mew=0.45, zorder=3)

def layout_check(fig, axes):
    fig.set_dpi(600)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    left, forest, estimate, qax = axes
    names = [t.get_window_extent(renderer) for t in left.texts if t.get_text() in WRAPPED.values()]
    column2 = [t.get_window_extent(renderer) for t in left.texts if t.get_text() not in WRAPPED.values()]
    disease_gap = min((b.x0 for b in column2)) - max((b.x1 for b in names))
    numeric_gap = min((t.get_window_extent(renderer).x0 for t in qax.texts)) - max((t.get_window_extent(renderer).x1 for t in estimate.texts))
    assert disease_gap / fig.dpi * 72 > 2, 'DISEASE_COLUMN2_OVERLAP'
    assert numeric_gap / fig.dpi * 72 > 2, 'EFFECT_Q_OVERLAP'
    texts = [*fig.texts]
    for ax in axes:
        texts.extend(ax.texts)
    texts.extend(forest.get_xticklabels())
    texts.append(forest.xaxis.label)
    for t in texts:
        box = t.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= fig.bbox.x1) and (box.y1 <= fig.bbox.y1), 'TEXT_OUTSIDE_CANVAS: ' + t.get_text()
    estimates = [t.get_window_extent(renderer) for t in estimate.texts]
    assert all((a.y0 > b.y1 for a, b in zip(estimates, estimates[1:]))), 'EFFECT_ROWS_OVERLAP'
    qboxes = [t.get_window_extent(renderer) for t in qax.texts]
    assert all((a.y0 > b.y1 for a, b in zip(qboxes, qboxes[1:]))), 'Q_ROWS_OVERLAP'
    return {'disease_column2_gap_pt': disease_gap / fig.dpi * 72, 'numeric_column_gap_pt': numeric_gap / fig.dpi * 72, 'all_text_within_canvas': True, 'effect_and_q_rows_nonoverlapping': True, 'axes_bounds': [list(ax.get_position().bounds) for ax in axes], 'x_limits': list(forest.get_xlim()), 'y_limits': list(forest.get_ylim())}

def render(letter, rows):
    assert len(rows) == (15 if letter == 'A' else 5)
    fig, axes = forest_axes('OR' if letter == 'A' else 'HR', 'Model' if letter == 'A' else 'N')
    left, forest, estimate, qax = axes
    for i, d in enumerate(ORDER):
        disease_rows = [r for r in rows if r['disease'] == d]
        assert len(disease_rows) == (3 if letter == 'A' else 1)
        for j, r in enumerate(disease_rows):
            if letter == 'A':
                assert r['model'] == ('M1_common', 'M2_common', 'M3_common')[j]
                y = 17 - 3.5 * i - j
                alpha, size = ((0.55, 0.75, 1)[j], 3.6 if j == 2 else 2.8)
                left.text(0.86, y, r['model'][:2], ha='center', va='center', fontsize=FONTS['model'], fontweight='bold' if j == 2 else 'normal')
            else:
                y, alpha, size = (16 - 3.5 * i, 1, 3.8)
                left.text(0.86, y, f"{int(r['N']):,}", ha='center', va='center', fontsize=FONTS['N'])
            interval(forest, *[float(r[k]) for k in ('effect', 'CI_low', 'CI_high')], y, COLORS[d], alpha, size)
            estimate.text(0, y, f"{float(r['effect']):.3f} ({float(r['CI_low']):.3f}–{float(r['CI_high']):.3f})", va='center', fontsize=FONTS['effect'])
            qax.text(0, y, qtext(r['BH_q']), va='center', fontsize=FONTS['q']).set_math_fontfamily('custom')
    return (fig, layout_check(fig, axes))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=os.getenv('BIRTH_PROJECT_ROOT'))
    parser.add_argument('--out', type=Path, default=os.getenv('BIRTH_FIGURE_OUT'))
    args = parser.parse_args()
    if args.project_root is None or args.out is None:
        parser.error('Explicit project root and output root required')
    root = args.project_root.resolve(strict=True)
    out = args.out.resolve()
    source_dir = root / 'reports/figure2_formal_20261103_v1/final_png/plot_data'
    sources = {letter: source_dir / f'Panel_{letter}_data.csv' for letter in 'AB'}
    for source in sources.values():
        if not source.is_file():
            parser.error('Missing frozen plot table: ' + str(source))
    print(f'INPUT_ROOT={source_dir}\nSCRIPT_ROOT={Path(__file__).parent}\nOUTPUT_ROOT={out}', flush=True)
    out.mkdir(parents=True, exist_ok=True)
    data = out / 'plot_data'
    data.mkdir(exist_ok=True)
    targets = [out / f'Panel_{letter}.{suffix}' for letter in 'AB' for suffix in ('png', 'svg')]
    targets.extend((out / f'Panel_{letter}_audit.json' for letter in 'AB'))
    targets.extend((data / f'Panel_{letter}_data.csv' for letter in 'AB'))
    assert not any((p.exists() for p in targets)), 'REFUSE_EXISTING_OUTPUT'
    settings()
    for letter, source in sources.items():
        before = sha(source)
        rows = readcsv(source)
        fig, audit = render(letter, rows)
        with (out / f'Panel_{letter}.png').open('xb') as f:
            fig.savefig(f, format='png', dpi=600)
        with (out / f'Panel_{letter}.svg').open('xb') as f:
            fig.savefig(f, format='svg')
        plt.close(fig)
        shutil.copyfile(source, data / source.name)
        assert sha(data / source.name) == before == sha(source)
        audit.update(panel=letter, figure_inches=list(SIZE), pixels=[2340, 1980], png_dpi=600, rows=len(rows), title=False, panel_letter=False, extra_annotations=False, font='Arial', disease_order=list(ORDER), disease_names_black=True, font_points=FONTS, original_font_points=ORIGINAL_FONTS, font_increase_pct={k: (FONTS[k] / v - 1) * 100 for k, v in ORIGINAL_FONTS.items()}, disease_label_wrap_change='Diabetes mellitus wrapped on two lines', plot_data_sha256=before, original_numbers_unchanged=True, inferential_status_flags_unchanged=True, model_refits_by_plot_code=0, original_marker_and_interval_styles_unchanged=True, source=str(source))
        with (out / f'Panel_{letter}_audit.json').open('x') as f:
            json.dump(audit, f, indent=2, ensure_ascii=False)
        print(json.dumps({'panel': letter, 'PNG': str(out / f'Panel_{letter}.png'), **audit}), flush=True)
    code = out / 'code'
    code.mkdir(exist_ok=True)
    target = code / Path(__file__).name
    assert not target.exists(), 'REFUSE_EXISTING_CODE_COPY'
    shutil.copyfile(Path(__file__), target)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Typography-only Figure 2 A/B revision from frozen final-panel plot tables."""\nfrom __future__ import annotations\nimport argparse\nimport csv\nimport hashlib\nimport json\nimport math\nimport os\nfrom pathlib import Path\nimport shutil\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib import font_manager, ticker\nORDER = (\'DM\', \'T2D\', \'IHD\', \'Angina\', \'CAS\')\nWRAPPED = {\'DM\': \'Diabetes\\nmellitus\', \'T2D\': \'Type 2 diabetes\', \'IHD\': \'Ischemic heart\\ndisease\', \'Angina\': \'Angina pectoris\', \'CAS\': \'Coronary\\natherosclerosis\'}\nBLUE, RED = (\'#5683A8\', \'#B65F64\')\nCOLORS = {d: BLUE if d in (\'DM\', \'T2D\') else RED for d in ORDER}\nSIZE = (3.9, 3.3)\nBOUNDS = ((0.015, 0.16, 0.305, 0.79), (0.345, 0.16, 0.225, 0.79), (0.59, 0.16, 0.245, 0.79), (0.855, 0.16, 0.13, 0.79))\nFONTS = {\'default\': 8.05, \'header\': 7.8, \'disease\': 8.3, \'model\': 7.5, \'N\': 7.1, \'effect\': 7.35, \'q\': 7.1, \'ticks\': 7.5, \'xlabel\': 7.6}\nORIGINAL_FONTS = {\'default\': 7, \'header\': 6.8, \'disease\': 7.2, \'model\': 6.5, \'N\': 6.2, \'effect\': 6.4, \'q\': 6.2, \'ticks\': 6.5, \'xlabel\': 6.6}\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef readcsv(path):\n    with path.open(encoding=\'utf-8-sig\', newline=\'\') as f:\n        return list(csv.DictReader(f))\n\ndef qtext(q):\n    q = float(q)\n    e = math.floor(math.log10(q))\n    s = f\'{q / 10 ** e:.2f}\' + f\'\\\\times10^{{{e}}}\'\n    return \'$\' + (\'\\\\mathbf{\' + s + \'}\' if q < 0.05 else s) + \'$\'\n\ndef settings():\n    for weight in (\'normal\', \'bold\'):\n        font_manager.findfont(font_manager.FontProperties(family=\'Arial\', weight=weight), fallback_to_default=False)\n    plt.rcParams.update({\'font.family\': \'Arial\', \'font.size\': FONTS[\'default\'], \'axes.unicode_minus\': False, \'savefig.facecolor\': \'white\', \'mathtext.fontset\': \'custom\', \'mathtext.rm\': \'Arial\', \'mathtext.it\': \'Arial:italic\', \'mathtext.bf\': \'Arial:bold\', \'mathtext.default\': \'regular\', \'mathtext.fallback\': None, \'axes.linewidth\': 0.6, \'text.color\': \'#222222\', \'axes.labelcolor\': \'#222222\', \'xtick.color\': \'#333333\', \'ytick.color\': \'#333333\'})\n\ndef forest_axes(measure, second):\n    fig = plt.figure(figsize=SIZE, facecolor=\'white\')\n    axes = [fig.add_axes(p) for p in BOUNDS]\n    left, forest, estimate, qax = axes\n    for ax in axes:\n        ax.set_ylim(0.2, 17.8)\n        ax.set_yticks([])\n    for ax in (left, estimate, qax):\n        ax.set_xlim(0, 1)\n        ax.axis(\'off\')\n    for label, x, align in [(\'Disease\', 0.015, \'left\'), (second, 0.015 + 0.305 * 0.86, \'center\'), (measure + \' (95% CI)\', 0.59, \'left\'), (\'BH q\', 0.855, \'left\')]:\n        fig.text(x, 0.965, label, fontsize=FONTS[\'header\'], fontweight=\'bold\', ha=align)\n    for i, d in enumerate(ORDER):\n        left.text(0, 16 - 3.5 * i, WRAPPED[d], va=\'center\', fontsize=FONTS[\'disease\'], linespacing=1.05, color=\'black\')\n        if i < 4:\n            for ax in axes:\n                ax.axhline(14.25 - 3.5 * i, lw=0.45, color=\'#E6E6E6\', zorder=0)\n    forest.set_xscale(\'log\')\n    forest.set_xlim(0.8, 1.3)\n    forest.xaxis.set_major_locator(ticker.FixedLocator([0.8, 1, 1.2]))\n    forest.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f\'{x:g}\'))\n    forest.xaxis.set_minor_locator(ticker.NullLocator())\n    forest.tick_params(axis=\'x\', length=2.5, width=0.6, pad=3, labelsize=FONTS[\'ticks\'])\n    forest.spines[[\'top\', \'left\', \'right\']].set_visible(False)\n    forest.spines[\'bottom\'].set_color(\'#777777\')\n    forest.axvline(1, color=\'#888888\', ls=(0, (3, 3)), lw=0.65, zorder=0)\n    forest.set_xlabel(\'Relative birth order (0→1)\\n\' + (\'Odds ratio\' if measure == \'OR\' else \'Hazard ratio\'), fontsize=FONTS[\'xlabel\'], labelpad=4)\n    return (fig, axes)\n\ndef interval(ax, x, lo, hi, y, color, alpha=1, size=3.8):\n    ax.plot([lo, hi], [y, y], color=color, lw=0.85, alpha=alpha, zorder=2)\n    for v in (lo, hi):\n        ax.plot([v, v], [y - 0.11, y + 0.11], color=color, lw=0.7, alpha=alpha, zorder=2)\n    ax.plot(x, y, \'o\', ms=size, color=color, alpha=alpha, mec=\'white\', mew=0.45, zorder=3)\n\ndef layout_check(fig, axes):\n    fig.set_dpi(600)\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    left, forest, estimate, qax = axes\n    names = [t.get_window_extent(renderer) for t in left.texts if t.get_text() in WRAPPED.values()]\n    column2 = [t.get_window_extent(renderer) for t in left.texts if t.get_text() not in WRAPPED.values()]\n    disease_gap = min((b.x0 for b in column2)) - max((b.x1 for b in names))\n    numeric_gap = min((t.get_window_extent(renderer).x0 for t in qax.texts)) - max((t.get_window_extent(renderer).x1 for t in estimate.texts))\n    assert disease_gap / fig.dpi * 72 > 2, \'DISEASE_COLUMN2_OVERLAP\'\n    assert numeric_gap / fig.dpi * 72 > 2, \'EFFECT_Q_OVERLAP\'\n    texts = [*fig.texts]\n    for ax in axes:\n        texts.extend(ax.texts)\n    texts.extend(forest.get_xticklabels())\n    texts.append(forest.xaxis.label)\n    for t in texts:\n        box = t.get_window_extent(renderer)\n        assert box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= fig.bbox.x1) and (box.y1 <= fig.bbox.y1), \'TEXT_OUTSIDE_CANVAS: \' + t.get_text()\n    estimates = [t.get_window_extent(renderer) for t in estimate.texts]\n    assert all((a.y0 > b.y1 for a, b in zip(estimates, estimates[1:]))), \'EFFECT_ROWS_OVERLAP\'\n    qboxes = [t.get_window_extent(renderer) for t in qax.texts]\n    assert all((a.y0 > b.y1 for a, b in zip(qboxes, qboxes[1:]))), \'Q_ROWS_OVERLAP\'\n    return {\'disease_column2_gap_pt\': disease_gap / fig.dpi * 72, \'numeric_column_gap_pt\': numeric_gap / fig.dpi * 72, \'all_text_within_canvas\': True, \'effect_and_q_rows_nonoverlapping\': True, \'axes_bounds\': [list(ax.get_position().bounds) for ax in axes], \'x_limits\': list(forest.get_xlim()), \'y_limits\': list(forest.get_ylim())}\n\ndef render(letter, rows):\n    assert len(rows) == (15 if letter == \'A\' else 5)\n    fig, axes = forest_axes(\'OR\' if letter == \'A\' else \'HR\', \'Model\' if letter == \'A\' else \'N\')\n    left, forest, estimate, qax = axes\n    for i, d in enumerate(ORDER):\n        disease_rows = [r for r in rows if r[\'disease\'] == d]\n        assert len(disease_rows) == (3 if letter == \'A\' else 1)\n        for j, r in enumerate(disease_rows):\n            if letter == \'A\':\n                assert r[\'model\'] == (\'M1_common\', \'M2_common\', \'M3_common\')[j]\n                y = 17 - 3.5 * i - j\n                alpha, size = ((0.55, 0.75, 1)[j], 3.6 if j == 2 else 2.8)\n                left.text(0.86, y, r[\'model\'][:2], ha=\'center\', va=\'center\', fontsize=FONTS[\'model\'], fontweight=\'bold\' if j == 2 else \'normal\')\n            else:\n                y, alpha, size = (16 - 3.5 * i, 1, 3.8)\n                left.text(0.86, y, f"{int(r[\'N\']):,}", ha=\'center\', va=\'center\', fontsize=FONTS[\'N\'])\n            interval(forest, *[float(r[k]) for k in (\'effect\', \'CI_low\', \'CI_high\')], y, COLORS[d], alpha, size)\n            estimate.text(0, y, f"{float(r[\'effect\']):.3f} ({float(r[\'CI_low\']):.3f}–{float(r[\'CI_high\']):.3f})", va=\'center\', fontsize=FONTS[\'effect\'])\n            qax.text(0, y, qtext(r[\'BH_q\']), va=\'center\', fontsize=FONTS[\'q\']).set_math_fontfamily(\'custom\')\n    return (fig, layout_check(fig, axes))\n\ndef main():\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--project-root\', type=Path, default=os.getenv(\'BIRTH_PROJECT_ROOT\'))\n    parser.add_argument(\'--out\', type=Path, default=os.getenv(\'BIRTH_FIGURE_OUT\'))\n    args = parser.parse_args()\n    if args.project_root is None or args.out is None:\n        parser.error(\'Explicit project root and output root required\')\n    root = args.project_root.resolve(strict=True)\n    out = args.out.resolve()\n    source_dir = root / \'reports/figure2_formal_20261103_v1/final_png/plot_data\'\n    sources = {letter: source_dir / f\'Panel_{letter}_data.csv\' for letter in \'AB\'}\n    for source in sources.values():\n        if not source.is_file():\n            parser.error(\'Missing frozen plot table: \' + str(source))\n    print(f\'INPUT_ROOT={source_dir}\\nSCRIPT_ROOT={Path(__file__).parent}\\nOUTPUT_ROOT={out}\', flush=True)\n    out.mkdir(parents=True, exist_ok=True)\n    data = out / \'plot_data\'\n    data.mkdir(exist_ok=True)\n    targets = [out / f\'Panel_{letter}.{suffix}\' for letter in \'AB\' for suffix in (\'png\', \'svg\')]\n    targets.extend((out / f\'Panel_{letter}_audit.json\' for letter in \'AB\'))\n    targets.extend((data / f\'Panel_{letter}_data.csv\' for letter in \'AB\'))\n    assert not any((p.exists() for p in targets)), \'REFUSE_EXISTING_OUTPUT\'\n    settings()\n    for letter, source in sources.items():\n        before = sha(source)\n        rows = readcsv(source)\n        fig, audit = render(letter, rows)\n        with (out / f\'Panel_{letter}.png\').open(\'xb\') as f:\n            fig.savefig(f, format=\'png\', dpi=600)\n        with (out / f\'Panel_{letter}.svg\').open(\'xb\') as f:\n            fig.savefig(f, format=\'svg\')\n        plt.close(fig)\n        shutil.copyfile(source, data / source.name)\n        assert sha(data / source.name) == before == sha(source)\n        audit.update(panel=letter, figure_inches=list(SIZE), pixels=[2340, 1980], png_dpi=600, rows=len(rows), title=False, panel_letter=False, extra_annotations=False, font=\'Arial\', disease_order=list(ORDER), disease_names_black=True, font_points=FONTS, original_font_points=ORIGINAL_FONTS, font_increase_pct={k: (FONTS[k] / v - 1) * 100 for k, v in ORIGINAL_FONTS.items()}, disease_label_wrap_change=\'Diabetes mellitus wrapped on two lines\', plot_data_sha256=before, original_numbers_unchanged=True, inferential_status_flags_unchanged=True, model_refits_by_plot_code=0, original_marker_and_interval_styles_unchanged=True, source=str(source))\n        with (out / f\'Panel_{letter}_audit.json\').open(\'x\') as f:\n            json.dump(audit, f, indent=2, ensure_ascii=False)\n        print(json.dumps({\'panel\': letter, \'PNG\': str(out / f\'Panel_{letter}.png\'), **audit}), flush=True)\n    code = out / \'code\'\n    code.mkdir(exist_ok=True)\n    target = code / Path(__file__).name\n    assert not target.exists(), \'REFUSE_EXISTING_CODE_COPY\'\n    shutil.copyfile(Path(__file__), target)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
