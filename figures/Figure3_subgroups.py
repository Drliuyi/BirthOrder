"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import Rectangle
from PIL import Image
RUN = 'figure2_typography_20261105_v1'
ORDER = ('DM', 'T2D', 'IHD', 'Angina', 'CAS')
NAMES = ('Diabetes\nmellitus', 'Type 2 diabetes', 'Ischemic heart\ndisease', 'Angina pectoris', 'Coronary\natherosclerosis')
COLS = ('Female', 'Male', '<50', '50–59', '≥60', 'S=2', 'S=3', 'S=4', 'S=5–8')
GROUPS = ((0, 2, 'Sex', '#8D86A0'), (2, 5, 'Baseline age', '#B8A081'), (5, 9, 'Sibship size', '#7F9D94'))
SIZE = (4.7, 3.04)
DPI = 600
FONT_PT = {'disease': 7.9, 'OR': 7.6, 'star': 8.7, 'top_subgroup': 7.8, 'colorbar_tick': 6.7, 'colorbar_title': 6.8}
TEXT_OFFSET = 0.18
STRIP = 0.18
REFERENCE = '/var/folders/w0/b8bjmt_14_g22kl1yk_8wxz80000gn/T/codex-clipboard-7d5b67bc-e987-45de-aa01-798e427ae9cf.png'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def root_path():
    if os.environ.get('PROJECT_ROOT'):
        return Path(os.environ['PROJECT_ROOT']).expanduser().resolve()
    for p in Path(__file__).resolve().parents:
        if (p / 'reports/figure2_cd_revision_20261103_v1/clean_final_v7/Panel_D_data.csv').is_file():
            return p
    raise FileNotFoundError('Set PROJECT_ROOT to the existing Birth Order workspace.')

def fonts():
    regular = Path(os.environ.get('ARIAL_FONT', '/System/Library/Fonts/Supplemental/Arial.ttf'))
    bold = Path(os.environ.get('ARIAL_BOLD_FONT', '/System/Library/Fonts/Supplemental/Arial Bold.ttf'))
    for p in (regular, bold):
        if not p.is_file():
            raise FileNotFoundError(p)
        font_manager.fontManager.addfont(str(p))
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 6.0, 'text.color': 'black', 'axes.unicode_minus': True, 'figure.dpi': DPI, 'savefig.dpi': DPI, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    return str(regular)

def luminance(color):
    rgb = to_rgb(color)
    rgb = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb]
    return sum((a * b for a, b in zip(rgb, (0.2126, 0.7152, 0.0722))))

def render(rows, old_audit):
    fig = plt.figure(figsize=SIZE, facecolor='white')
    labels = fig.add_axes((LABEL_LEFT, BODY_BOTTOM, 0.17, BODY_HEIGHT))
    heat = fig.add_axes((HEAT_LEFT, BODY_BOTTOM, HEAT_WIDTH, BODY_HEIGHT))
    for ax in (labels, heat):
        ax.set_ylim(-0.5, 4.5)
        ax.axis('off')
    labels.set_xlim(0, 1)
    heat.set_xlim(0, 9)
    cmap = LinearSegmentedColormap.from_list('OR_unchanged', [BLUE, '#FFFFFF', RED], N=257)
    limit = float(old_audit['visual_QA']['symmetric_log_OR_limit'])
    norm = Normalize(-limit, limit)
    texts = []
    values = []
    pairs = []
    patches = []
    name_texts = []
    stars = []
    for i, name in enumerate(NAMES):
        y = 4 - i
        t = labels.text(0, y, name, ha='left', va='center', fontsize=FONT_PT['disease'], fontweight='normal', linespacing=1.15)
        texts.append(t)
        name_texts.append(t)
        heat.add_patch(Rectangle((-STRIP, y - 0.5), STRIP, 1, facecolor=BLUE if i < 2 else RED, edgecolor='none', clip_on=False))
    prior = old_audit['visual_QA']['cell_visual_encodings']
    fills = []
    for r in rows:
        i, j = (int(r['row_order']) - 1, int(r['column_order']) - 1)
        y = 4 - i
        x = j
        color = cmap(norm(float(r['log_OR_for_fill'])))
        old = next((v for v in prior if v['disease'] == r['disease'] and v['column'] == r['column_label']))
        assert all((abs(a - b) < 1e-12 for a, b in zip(color, old['fill_RGBA']))), 'MATRIX_COLOR_CHANGED'
        patch = Rectangle((x, y - 0.5), 1, 1, facecolor=color, edgecolor='white', linewidth=0.35)
        heat.add_patch(patch)
        patches.append(patch)
        number = heat.text(x + 0.5, y - TEXT_OFFSET, r['OR_display'], ha='center', va='center', fontsize=FONT_PT['OR'], color='#171717')
        texts.append(number)
        values.append(number)
        if r['significance_symbol'] == '*':
            star = heat.text(x + 0.5, y + TEXT_OFFSET, '*', ha='center', va='center', fontsize=FONT_PT['star'], color='black')
            texts.append(star)
            stars.append(star)
            pairs.append((number, star, patch))
        fills.append({'disease': r['disease'], 'column': r['column_label'], 'RGBA': list(color), 'OR_text': r['OR_display'], 'symbol': r['significance_symbol']})
    top_strips = []
    for lo, hi, name, color in GROUPS:
        p = Rectangle((lo, 4.5), hi - lo, STRIP, facecolor=color, edgecolor='none', clip_on=False)
        heat.add_patch(p)
        top_strips.append(p)
    column_texts = []
    for j, label in enumerate(COLS):
        t = heat.text(j + 0.5, 1.095, label, transform=heat.get_xaxis_transform(), ha='center', va='center', fontsize=FONT_PT['top_subgroup'], clip_on=False)
        texts.append(t)
        column_texts.append(t)
    cax = fig.add_axes((148 / 385, 29 / 249, 123 / 385, 6 / 249))
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
    cb.set_ticks([math.log(0.8), 0, math.log(1.25)], labels=['0.80', '1.00', '1.25'])
    cb.ax.tick_params(labelsize=FONT_PT['colorbar_tick'], length=2, width=0.4, pad=2)
    cb.outline.set_visible(False)
    cb.set_label('Odds ratio', fontsize=FONT_PT['colorbar_title'], labelpad=2)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    texts.extend(cb.ax.get_xticklabels())
    texts.append(cb.ax.xaxis.label)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in texts:
        box = text.get_window_extent(renderer)
        assert 2 < box.x0 < box.x1 < fig.bbox.width - 2 and 2 < box.y0 < box.y1 < fig.bbox.height - 2, ('CLIPPED_TEXT', text.get_text())
    for left, right in zip(column_texts, column_texts[1:]):
        assert left.get_window_extent(renderer).x1 + 5 < right.get_window_extent(renderer).x0, 'COLUMN_LABEL_OVERLAP'
    for number, star, patch in pairs:
        nb, sb, pb = [a.get_window_extent(renderer) for a in (number, star, patch)]
        assert nb.y1 + 5 < sb.y0 and pb.y0 + 2 < nb.y0 < sb.y1 < pb.y1 - 2, 'CELL_TEXT_OVERLAP'
    stripe_left = heat.transData.transform((-STRIP, 0))[0]
    gap = stripe_left - max((t.get_window_extent(renderer).x1 for t in name_texts))
    assert gap > 8, 'DISEASE_LABEL_OVERLAP'
    origin = heat.transData.transform((0, 0))
    one = heat.transData.transform((1, 1))
    cell_w, cell_h = one - origin
    assert abs(cell_w - cell_h) < 1e-08, 'NONSQUARE_CELL'
    left_width = cell_w * STRIP
    top_height = top_strips[0].get_window_extent(renderer).height
    assert abs(left_width - top_height) < 1e-08, 'STRIP_THICKNESS_NOT_EQUAL'
    assert len(values) == 45 and len(stars) == 34 and (len(fig.legends) == 0)
    assert not heat.lines and all((p.get_linestyle() in ('solid', '-') for p in heat.patches))
    assert tuple(cmap(norm(0)))[:3] == (1.0, 1.0, 1.0)
    old_top = ('#E7E4EE', '#F4ECD9', '#E2EBE5')
    contrasts = []
    for (_, _, name, color), previous in zip(GROUPS, old_top):
        new_contrast = 1.05 / (luminance(color) + 0.05)
        old_contrast = 1.05 / (luminance(previous) + 0.05)
        assert new_contrast > old_contrast * 2
        contrasts.append({'group': name, 'new_hex': color, 'old_hex': previous, 'alpha': 1.0, 'new_white_contrast': new_contrast, 'old_white_contrast': old_contrast})
    qa = {'text_clipping': 'PASS', 'cell_text_nonoverlap': 'PASS', 'all45_OR_texts_and34_stars_preserved': True, 'matrix_fill_colors_unchanged': True, 'symmetric_log_OR_limits': [-limit, limit], 'cell_side_pixels': [float(cell_w), float(cell_h)], 'top_strip_height_pixels': float(top_height), 'left_strip_width_pixels': float(left_width), 'strip_thickness_equal': 'PASS', 'all_cells_and_strips_touch': True, 'top_group_contrasts': contrasts, 'full_disease_names_black_regular': True, 'disease_label_to_strip_min_gap_pixels': float(gap), 'disease_label_to_strip_min_gap_pt': float(gap / DPI * 72), 'disease_font_pt': FONT_PT['disease'], 'font_sizes_pt': FONT_PT, 'cell_text_center_offsets': {'OR': -TEXT_OFFSET, 'star': TEXT_OFFSET}, 'reference_layout': 'User-adjusted 384x249 screenshot; compact label-to-grid spacing reproduced', 'figure_inches': list(SIZE), 'heat_axes_fraction': list(heat.get_position().bounds), 'no_title_subtitle_panel_letter_legend_note_dagger_or_dashed_box': True, 'cell_encodings': fills}
    return (fig, qa)

def main():
    project = root_path()
    parser = argparse.ArgumentParser(description=__doc__)
    default_source = project / 'reports/figure2_cd_revision_20261103_v1/clean_final_v7'
    parser.add_argument('--source-dir', type=Path, default=Path(os.environ.get('FIGURE2_D_SOURCE_ROOT', default_source)))
    parser.add_argument('--output', type=Path, default=Path(os.environ.get('FIGURE2_D_OUTPUT_ROOT', project / 'reports' / RUN / 'D')))
    args = parser.parse_args()
    source = args.source_dir.resolve()
    out = args.output.resolve()
    print(f'PROJECT_ROOT={project}\nSOURCE_ROOT={source}\nSCRIPT_ROOT={Path(__file__).parent}\nOUTPUT_ROOT={out}')
    data = source / 'Panel_D_data.csv'
    audit_path = source / 'Panel_D_audit.json'
    if not data.is_file() or not audit_path.is_file():
        raise FileNotFoundError('Missing exact prior D plot data or audit.')
    old = json.loads(audit_path.read_text())
    before = sha(data)
    assert before == old['output_sha256']['Panel_D_data.csv'], 'PRIOR_D_DATA_HASH_CHANGED'
    with data.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 45 and len({(r['disease'], r['column_order']) for r in rows}) == 45
    assert {r['disease'] for r in rows} == set(ORDER) and sum((r['significance_symbol'] == '*' for r in rows)) == 34
    for i, d in enumerate(ORDER):
        rr = [r for r in rows if r['disease'] == d]
        assert [int(r['column_order']) for r in rr] == list(range(1, 10)) and all((int(r['row_order']) == i + 1 for r in rr))
    targets = [out / n for n in ('Panel_D.png', 'Panel_D.svg', 'Panel_D_data.csv', 'Panel_D_audit.json', 'README.md', 'code/render_compact_d.py')]
    if any((p.exists() for p in targets)):
        raise FileExistsError('Choose a fresh output; previous files are never replaced.')
    font = fonts()
    fig, qa = render(rows, old)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'svg'):
        with (out / f'Panel_D.{ext}').open('xb') as handle:
            fig.savefig(handle, format=ext, dpi=DPI, facecolor='white')
    plt.close(fig)
    with Image.open(out / 'Panel_D.png') as im:
        assert im.size == (2820, 1824)
    assert '<text' in (out / 'Panel_D.svg').read_text() and '<image' not in (out / 'Panel_D.svg').read_text()
    with (out / 'Panel_D_data.csv').open('xb') as handle:
        handle.write(data.read_bytes())
    assert sha(out / 'Panel_D_data.csv') == before == sha(data), 'SOURCE_DATA_CHANGED'
    (out / 'code').mkdir(exist_ok=True)
    with (out / 'code/render_compact_d.py').open('xb') as handle:
        handle.write(Path(__file__).read_bytes())
    with (out / 'README.md').open('x', encoding='utf-8') as handle:
        handle.write('# Compact D larger-typography revision\n\nTypography-only revision of the authenticated muted-colour heatmap. Top subgroup labels increase from 5.9 to 7.8 pt (+32%); disease labels 6.9 to 7.9 pt; cell OR 6.5 to 7.6 pt; stars 7.4 to 8.7 pt; colorbar ticks 5.6 to 6.7 pt; colorbar title 5.7 to 6.8 pt. The first disease name wraps on two lines; the disease-label left edge shifts by 1/385 of the canvas to preserve compact readable spacing. Cell text separates vertically to keep each star above its OR.\n\nThe 4.7 x 3.04-inch canvas, heatmap geometry, exact 45 OR/q data rows, 34 original stars, all blue/red cell RGBA fills and muted category-strip colours are preserved. No fitting, q recalculation, or remote work.\n\nTop-strip height exactly equals left disease-strip width. Every matrix cell remains square and touching its neighbours. Main output is 600-dpi PNG; editable SVG is also saved. The original continuous odds-ratio colorbar is retained; no panel letter, title, subtitle, classification legend, notes, dagger, or dashed box is added.\n\nStars retain the previous definition: valid original unadjusted 95% CI excludes 1; they do not denote BH significance. The pending CAS/<50 estimate retains its original source status and no star.\n\nReproduce into a new directory:\n\n```bash\nPROJECT_ROOT=' + str(project) + ' ' + str(project / '.codex_tmp/border_figures_20261016_venv/bin/python') + ' -B ' + str(out / 'code/render_compact_d.py') + ' --source-dir ' + str(source) + ' --output /absolute/path/to/NEW_D_DIRECTORY\n```\n')
    result = {'created_UTC': datetime.now(timezone.utc).isoformat(), 'source_data_sha256': before, 'source_data_unchanged': True, 'prior_audit_sha256': sha(audit_path), 'reference_image': REFERENCE, 'reference_image_sha256': sha(Path(REFERENCE)) if Path(REFERENCE).is_file() else None, 'model_fits': 0, 'new_P_q_calculations': 0, 'remote_operations': 0, 'fonts': font, 'DPI': DPI, 'PNG_pixels': [2820, 1824], 'visual_QA': qa, 'output_sha256': {p.name: sha(p) for p in targets if p.exists()}}
    with (out / 'Panel_D_audit.json').open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    print(json.dumps({'status': 'PASS', 'output': str(out), 'PNG_pixels': [2820, 1824], 'cell_pixels': qa['cell_side_pixels'], 'strip_thickness_pixels': qa['top_strip_height_pixels']}))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Enlarge typography in the authenticated compact D heatmap.\n\nAll heatmap estimates, fills, stars and inference fields come verbatim from the\nprevious authenticated plot data. No statistical computation or remote access.\n"""\nfrom __future__ import annotations\nimport argparse\nimport csv\nfrom datetime import datetime, timezone\nimport hashlib\nimport json\nimport math\nimport os\nfrom pathlib import Path\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib import font_manager\nfrom matplotlib.cm import ScalarMappable\nfrom matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb\nfrom matplotlib.patches import Rectangle\nfrom PIL import Image\nRUN = \'figure2_typography_20261105_v1\'\nORDER = (\'DM\', \'T2D\', \'IHD\', \'Angina\', \'CAS\')\nNAMES = (\'Diabetes\\nmellitus\', \'Type 2 diabetes\', \'Ischemic heart\\ndisease\', \'Angina pectoris\', \'Coronary\\natherosclerosis\')\nCOLS = (\'Female\', \'Male\', \'<50\', \'50–59\', \'≥60\', \'S=2\', \'S=3\', \'S=4\', \'S=5–8\')\nBLUE, RED = (\'#5683A8\', \'#B65F64\')\nGROUPS = ((0, 2, \'Sex\', \'#8D86A0\'), (2, 5, \'Baseline age\', \'#B8A081\'), (5, 9, \'Sibship size\', \'#7F9D94\'))\nSIZE = (4.7, 3.04)\nDPI = 600\nFONT_PT = {\'disease\': 7.9, \'OR\': 7.6, \'star\': 8.7, \'top_subgroup\': 7.8, \'colorbar_tick\': 6.7, \'colorbar_title\': 6.8}\nTEXT_OFFSET = 0.18\nLABEL_LEFT = 7 / 385\nSTRIP = 0.18\nHEAT_LEFT = 78 / 385\nHEAT_WIDTH = 287 / 385\nBODY_BOTTOM = 61 / 249\nBODY_HEIGHT = 5 * SIZE[0] * HEAT_WIDTH / 9 / SIZE[1]\nREFERENCE = \'/var/folders/w0/b8bjmt_14_g22kl1yk_8wxz80000gn/T/codex-clipboard-7d5b67bc-e987-45de-aa01-798e427ae9cf.png\'\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef root_path():\n    if os.environ.get(\'PROJECT_ROOT\'):\n        return Path(os.environ[\'PROJECT_ROOT\']).expanduser().resolve()\n    for p in Path(__file__).resolve().parents:\n        if (p / \'reports/figure2_cd_revision_20261103_v1/clean_final_v7/Panel_D_data.csv\').is_file():\n            return p\n    raise FileNotFoundError(\'Set PROJECT_ROOT to the existing Birth Order workspace.\')\n\ndef fonts():\n    regular = Path(os.environ.get(\'ARIAL_FONT\', \'/System/Library/Fonts/Supplemental/Arial.ttf\'))\n    bold = Path(os.environ.get(\'ARIAL_BOLD_FONT\', \'/System/Library/Fonts/Supplemental/Arial Bold.ttf\'))\n    for p in (regular, bold):\n        if not p.is_file():\n            raise FileNotFoundError(p)\n        font_manager.fontManager.addfont(str(p))\n    plt.rcParams.update({\'font.family\': \'Arial\', \'font.size\': 6.0, \'text.color\': \'black\', \'axes.unicode_minus\': True, \'figure.dpi\': DPI, \'savefig.dpi\': DPI, \'svg.fonttype\': \'none\', \'pdf.fonttype\': 42})\n    return str(regular)\n\ndef luminance(color):\n    rgb = to_rgb(color)\n    rgb = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb]\n    return sum((a * b for a, b in zip(rgb, (0.2126, 0.7152, 0.0722))))\n\ndef render(rows, old_audit):\n    fig = plt.figure(figsize=SIZE, facecolor=\'white\')\n    labels = fig.add_axes((LABEL_LEFT, BODY_BOTTOM, 0.17, BODY_HEIGHT))\n    heat = fig.add_axes((HEAT_LEFT, BODY_BOTTOM, HEAT_WIDTH, BODY_HEIGHT))\n    for ax in (labels, heat):\n        ax.set_ylim(-0.5, 4.5)\n        ax.axis(\'off\')\n    labels.set_xlim(0, 1)\n    heat.set_xlim(0, 9)\n    cmap = LinearSegmentedColormap.from_list(\'OR_unchanged\', [BLUE, \'#FFFFFF\', RED], N=257)\n    limit = float(old_audit[\'visual_QA\'][\'symmetric_log_OR_limit\'])\n    norm = Normalize(-limit, limit)\n    texts = []\n    values = []\n    pairs = []\n    patches = []\n    name_texts = []\n    stars = []\n    for i, name in enumerate(NAMES):\n        y = 4 - i\n        t = labels.text(0, y, name, ha=\'left\', va=\'center\', fontsize=FONT_PT[\'disease\'], fontweight=\'normal\', linespacing=1.15)\n        texts.append(t)\n        name_texts.append(t)\n        heat.add_patch(Rectangle((-STRIP, y - 0.5), STRIP, 1, facecolor=BLUE if i < 2 else RED, edgecolor=\'none\', clip_on=False))\n    prior = old_audit[\'visual_QA\'][\'cell_visual_encodings\']\n    fills = []\n    for r in rows:\n        i, j = (int(r[\'row_order\']) - 1, int(r[\'column_order\']) - 1)\n        y = 4 - i\n        x = j\n        color = cmap(norm(float(r[\'log_OR_for_fill\'])))\n        old = next((v for v in prior if v[\'disease\'] == r[\'disease\'] and v[\'column\'] == r[\'column_label\']))\n        assert all((abs(a - b) < 1e-12 for a, b in zip(color, old[\'fill_RGBA\']))), \'MATRIX_COLOR_CHANGED\'\n        patch = Rectangle((x, y - 0.5), 1, 1, facecolor=color, edgecolor=\'white\', linewidth=0.35)\n        heat.add_patch(patch)\n        patches.append(patch)\n        number = heat.text(x + 0.5, y - TEXT_OFFSET, r[\'OR_display\'], ha=\'center\', va=\'center\', fontsize=FONT_PT[\'OR\'], color=\'#171717\')\n        texts.append(number)\n        values.append(number)\n        if r[\'significance_symbol\'] == \'*\':\n            star = heat.text(x + 0.5, y + TEXT_OFFSET, \'*\', ha=\'center\', va=\'center\', fontsize=FONT_PT[\'star\'], color=\'black\')\n            texts.append(star)\n            stars.append(star)\n            pairs.append((number, star, patch))\n        fills.append({\'disease\': r[\'disease\'], \'column\': r[\'column_label\'], \'RGBA\': list(color), \'OR_text\': r[\'OR_display\'], \'symbol\': r[\'significance_symbol\']})\n    top_strips = []\n    for lo, hi, name, color in GROUPS:\n        p = Rectangle((lo, 4.5), hi - lo, STRIP, facecolor=color, edgecolor=\'none\', clip_on=False)\n        heat.add_patch(p)\n        top_strips.append(p)\n    column_texts = []\n    for j, label in enumerate(COLS):\n        t = heat.text(j + 0.5, 1.095, label, transform=heat.get_xaxis_transform(), ha=\'center\', va=\'center\', fontsize=FONT_PT[\'top_subgroup\'], clip_on=False)\n        texts.append(t)\n        column_texts.append(t)\n    cax = fig.add_axes((148 / 385, 29 / 249, 123 / 385, 6 / 249))\n    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation=\'horizontal\')\n    cb.set_ticks([math.log(0.8), 0, math.log(1.25)], labels=[\'0.80\', \'1.00\', \'1.25\'])\n    cb.ax.tick_params(labelsize=FONT_PT[\'colorbar_tick\'], length=2, width=0.4, pad=2)\n    cb.outline.set_visible(False)\n    cb.set_label(\'Odds ratio\', fontsize=FONT_PT[\'colorbar_title\'], labelpad=2)\n    if cb.solids is not None:\n        cb.solids.set_rasterized(False)\n    texts.extend(cb.ax.get_xticklabels())\n    texts.append(cb.ax.xaxis.label)\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    for text in texts:\n        box = text.get_window_extent(renderer)\n        assert 2 < box.x0 < box.x1 < fig.bbox.width - 2 and 2 < box.y0 < box.y1 < fig.bbox.height - 2, (\'CLIPPED_TEXT\', text.get_text())\n    for left, right in zip(column_texts, column_texts[1:]):\n        assert left.get_window_extent(renderer).x1 + 5 < right.get_window_extent(renderer).x0, \'COLUMN_LABEL_OVERLAP\'\n    for number, star, patch in pairs:\n        nb, sb, pb = [a.get_window_extent(renderer) for a in (number, star, patch)]\n        assert nb.y1 + 5 < sb.y0 and pb.y0 + 2 < nb.y0 < sb.y1 < pb.y1 - 2, \'CELL_TEXT_OVERLAP\'\n    stripe_left = heat.transData.transform((-STRIP, 0))[0]\n    gap = stripe_left - max((t.get_window_extent(renderer).x1 for t in name_texts))\n    assert gap > 8, \'DISEASE_LABEL_OVERLAP\'\n    origin = heat.transData.transform((0, 0))\n    one = heat.transData.transform((1, 1))\n    cell_w, cell_h = one - origin\n    assert abs(cell_w - cell_h) < 1e-08, \'NONSQUARE_CELL\'\n    left_width = cell_w * STRIP\n    top_height = top_strips[0].get_window_extent(renderer).height\n    assert abs(left_width - top_height) < 1e-08, \'STRIP_THICKNESS_NOT_EQUAL\'\n    assert len(values) == 45 and len(stars) == 34 and (len(fig.legends) == 0)\n    assert not heat.lines and all((p.get_linestyle() in (\'solid\', \'-\') for p in heat.patches))\n    assert tuple(cmap(norm(0)))[:3] == (1.0, 1.0, 1.0)\n    old_top = (\'#E7E4EE\', \'#F4ECD9\', \'#E2EBE5\')\n    contrasts = []\n    for (_, _, name, color), previous in zip(GROUPS, old_top):\n        new_contrast = 1.05 / (luminance(color) + 0.05)\n        old_contrast = 1.05 / (luminance(previous) + 0.05)\n        assert new_contrast > old_contrast * 2\n        contrasts.append({\'group\': name, \'new_hex\': color, \'old_hex\': previous, \'alpha\': 1.0, \'new_white_contrast\': new_contrast, \'old_white_contrast\': old_contrast})\n    qa = {\'text_clipping\': \'PASS\', \'cell_text_nonoverlap\': \'PASS\', \'all45_OR_texts_and34_stars_preserved\': True, \'matrix_fill_colors_unchanged\': True, \'symmetric_log_OR_limits\': [-limit, limit], \'cell_side_pixels\': [float(cell_w), float(cell_h)], \'top_strip_height_pixels\': float(top_height), \'left_strip_width_pixels\': float(left_width), \'strip_thickness_equal\': \'PASS\', \'all_cells_and_strips_touch\': True, \'top_group_contrasts\': contrasts, \'full_disease_names_black_regular\': True, \'disease_label_to_strip_min_gap_pixels\': float(gap), \'disease_label_to_strip_min_gap_pt\': float(gap / DPI * 72), \'disease_font_pt\': FONT_PT[\'disease\'], \'font_sizes_pt\': FONT_PT, \'cell_text_center_offsets\': {\'OR\': -TEXT_OFFSET, \'star\': TEXT_OFFSET}, \'reference_layout\': \'User-adjusted 384x249 screenshot; compact label-to-grid spacing reproduced\', \'figure_inches\': list(SIZE), \'heat_axes_fraction\': list(heat.get_position().bounds), \'no_title_subtitle_panel_letter_legend_note_dagger_or_dashed_box\': True, \'cell_encodings\': fills}\n    return (fig, qa)\n\ndef main():\n    project = root_path()\n    parser = argparse.ArgumentParser(description=__doc__)\n    default_source = project / \'reports/figure2_cd_revision_20261103_v1/clean_final_v7\'\n    parser.add_argument(\'--source-dir\', type=Path, default=Path(os.environ.get(\'FIGURE2_D_SOURCE_ROOT\', default_source)))\n    parser.add_argument(\'--output\', type=Path, default=Path(os.environ.get(\'FIGURE2_D_OUTPUT_ROOT\', project / \'reports\' / RUN / \'D\')))\n    args = parser.parse_args()\n    source = args.source_dir.resolve()\n    out = args.output.resolve()\n    print(f\'PROJECT_ROOT={project}\\nSOURCE_ROOT={source}\\nSCRIPT_ROOT={Path(__file__).parent}\\nOUTPUT_ROOT={out}\')\n    data = source / \'Panel_D_data.csv\'\n    audit_path = source / \'Panel_D_audit.json\'\n    if not data.is_file() or not audit_path.is_file():\n        raise FileNotFoundError(\'Missing exact prior D plot data or audit.\')\n    old = json.loads(audit_path.read_text())\n    before = sha(data)\n    assert before == old[\'output_sha256\'][\'Panel_D_data.csv\'], \'PRIOR_D_DATA_HASH_CHANGED\'\n    with data.open(encoding=\'utf-8-sig\', newline=\'\') as handle:\n        rows = list(csv.DictReader(handle))\n    assert len(rows) == 45 and len({(r[\'disease\'], r[\'column_order\']) for r in rows}) == 45\n    assert {r[\'disease\'] for r in rows} == set(ORDER) and sum((r[\'significance_symbol\'] == \'*\' for r in rows)) == 34\n    for i, d in enumerate(ORDER):\n        rr = [r for r in rows if r[\'disease\'] == d]\n        assert [int(r[\'column_order\']) for r in rr] == list(range(1, 10)) and all((int(r[\'row_order\']) == i + 1 for r in rr))\n    targets = [out / n for n in (\'Panel_D.png\', \'Panel_D.svg\', \'Panel_D_data.csv\', \'Panel_D_audit.json\', \'README.md\', \'code/render_compact_d.py\')]\n    if any((p.exists() for p in targets)):\n        raise FileExistsError(\'Choose a fresh output; previous files are never replaced.\')\n    font = fonts()\n    fig, qa = render(rows, old)\n    out.mkdir(parents=True, exist_ok=True)\n    for ext in (\'png\', \'svg\'):\n        with (out / f\'Panel_D.{ext}\').open(\'xb\') as handle:\n            fig.savefig(handle, format=ext, dpi=DPI, facecolor=\'white\')\n    plt.close(fig)\n    with Image.open(out / \'Panel_D.png\') as im:\n        assert im.size == (2820, 1824)\n    assert \'<text\' in (out / \'Panel_D.svg\').read_text() and \'<image\' not in (out / \'Panel_D.svg\').read_text()\n    with (out / \'Panel_D_data.csv\').open(\'xb\') as handle:\n        handle.write(data.read_bytes())\n    assert sha(out / \'Panel_D_data.csv\') == before == sha(data), \'SOURCE_DATA_CHANGED\'\n    (out / \'code\').mkdir(exist_ok=True)\n    with (out / \'code/render_compact_d.py\').open(\'xb\') as handle:\n        handle.write(Path(__file__).read_bytes())\n    with (out / \'README.md\').open(\'x\', encoding=\'utf-8\') as handle:\n        handle.write(\'# Compact D larger-typography revision\\n\\nTypography-only revision of the authenticated muted-colour heatmap. Top subgroup labels increase from 5.9 to 7.8 pt (+32%); disease labels 6.9 to 7.9 pt; cell OR 6.5 to 7.6 pt; stars 7.4 to 8.7 pt; colorbar ticks 5.6 to 6.7 pt; colorbar title 5.7 to 6.8 pt. The first disease name wraps on two lines; the disease-label left edge shifts by 1/385 of the canvas to preserve compact readable spacing. Cell text separates vertically to keep each star above its OR.\\n\\nThe 4.7 x 3.04-inch canvas, heatmap geometry, exact 45 OR/q data rows, 34 original stars, all blue/red cell RGBA fills and muted category-strip colours are preserved. No fitting, q recalculation, or remote work.\\n\\nTop-strip height exactly equals left disease-strip width. Every matrix cell remains square and touching its neighbours. Main output is 600-dpi PNG; editable SVG is also saved. The original continuous odds-ratio colorbar is retained; no panel letter, title, subtitle, classification legend, notes, dagger, or dashed box is added.\\n\\nStars retain the previous definition: valid original unadjusted 95% CI excludes 1; they do not denote BH significance. The pending CAS/<50 estimate retains its original source status and no star.\\n\\nReproduce into a new directory:\\n\\n```bash\\nPROJECT_ROOT=\' + str(project) + \' \' + str(project / \'.codex_tmp/border_figures_20261016_venv/bin/python\') + \' -B \' + str(out / \'code/render_compact_d.py\') + \' --source-dir \' + str(source) + \' --output /absolute/path/to/NEW_D_DIRECTORY\\n```\\n\')\n    result = {\'created_UTC\': datetime.now(timezone.utc).isoformat(), \'source_data_sha256\': before, \'source_data_unchanged\': True, \'prior_audit_sha256\': sha(audit_path), \'reference_image\': REFERENCE, \'reference_image_sha256\': sha(Path(REFERENCE)) if Path(REFERENCE).is_file() else None, \'model_fits\': 0, \'new_P_q_calculations\': 0, \'remote_operations\': 0, \'fonts\': font, \'DPI\': DPI, \'PNG_pixels\': [2820, 1824], \'visual_QA\': qa, \'output_sha256\': {p.name: sha(p) for p in targets if p.exists()}}\n    with (out / \'Panel_D_audit.json\').open(\'x\', encoding=\'utf-8\') as handle:\n        json.dump(result, handle, indent=2, ensure_ascii=False)\n        handle.write(\'\\n\')\n    print(json.dumps({\'status\': \'PASS\', \'output\': str(out), \'PNG_pixels\': [2820, 1824], \'cell_pixels\': qa[\'cell_side_pixels\'], \'strip_thickness_pixels\': qa[\'top_strip_height_pixels\']}))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
