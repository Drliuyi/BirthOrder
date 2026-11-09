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
from pathlib import Path
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
from matplotlib.patches import Polygon
import matplotlib.patheffects as pe
import numpy as np
PALETTE = {'NP': '#009E73', 'PP': '#D55E00', 'NN': '#0072B2', 'PN': '#E69F00'}
STRENGTHS = (0.2, 0.36, 0.52, 0.68, 0.84, 1.0)
LABELS_ALL = {'NP': ['olink__fcrl6', 'clinical__glucose', 'nmr__nmr_23404', 'nmr__nmr_23523', 'nmr__nmr_23506'], 'PP': ['qualification', 'smoking', 'age_education', 'income', 'nmr__nmr_23451'], 'NN': ['birth_weight', 'bmi', 'olink__gzma', 'nmr__nmr_23613', 'nmr__nmr_23406'], 'PN': ['employment_employed', 'participation_religious', 'employment_voluntary', 'confiding', 'employment_student']}
LABELS_NMR = {'NP': ['nmr__nmr_23404', 'nmr__nmr_23523', 'nmr__nmr_23506', 'nmr__nmr_23526', 'nmr__nmr_23584'], 'PP': ['nmr__nmr_23451'], 'NN': ['nmr__nmr_23613', 'nmr__nmr_23406', 'nmr__nmr_23440', 'nmr__nmr_23608', 'nmr__nmr_23648'], 'PN': []}
DISPLAY = {'nmr__nmr_23404': 'NMR LDL-C', 'nmr__nmr_23406': 'NMR HDL-C', 'nmr__nmr_23440': 'NMR ApoA1', 'nmr__nmr_23451': 'Omega-3 FA%', 'nmr__nmr_23506': 'M-VLDL CE', 'nmr__nmr_23523': 'IDL-P', 'nmr__nmr_23526': 'IDL-C', 'nmr__nmr_23584': 'Very large VLDL PL%', 'nmr__nmr_23608': 'XS-VLDL TG%', 'nmr__nmr_23613': 'IDL TG%', 'nmr__nmr_23648': 'Small HDL TG%', 'employment_employed': 'Employed / self-employed', 'employment_student': 'Student', 'employment_voluntary': 'Unpaid / voluntary work', 'participation_religious': 'Religious group', 'confiding': 'Confiding frequency'}

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_csv(name, rows):
    with (OUTPUT / name).open('x', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def place_nearby_labels(ax, items, dataset, label):
    """Bounded display-space placement; avoid every point and every other label.

    The fixed RNG is for visual layout only, never analysis or Bootstrap.
    Candidates are scored by point-to-text-edge distance, not fixed plot slots.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    extent = ax.get_window_extent(renderer)
    cloud = ax.transData.transform([(x, y) for _, x, y, _ in dataset])
    anchors = ax.transData.transform([(x, y) for _, x, y, _ in items])
    radius = math.sqrt(88) * fig.dpi / 72 / 2 + 2.4
    candidates, distances, centers = ([], [], [])
    texts = []
    for idx, (row, x, y, quad) in enumerate(items):
        text = ax.text(x, y, label(row), ha='center', va='center', fontsize=14.5, color='black', zorder=6)
        bbox = text.get_window_extent(renderer)
        width, height = (bbox.width + 6, bbox.height + 6)
        boxes, costs, cent = ([], [], [])
        px, py = anchors[idx]
        for gap in (9, 14, 20, 28, 38, 50, 66, 86, 110, 140, 175, 215):
            for angle in np.linspace(0, 2 * math.pi, 40, endpoint=False):
                cx = px + math.cos(angle) * (width / 2 + gap)
                cy = py + math.sin(angle) * (height / 2 + gap)
                box = np.array([cx - width / 2, cy - height / 2, cx + width / 2, cy + height / 2])
                if box[0] < extent.x0 + 4 or box[2] > extent.x1 - 4 or box[1] < extent.y0 + 4 or (box[3] > extent.y1 - 4):
                    continue
                dx = np.maximum(np.maximum(box[0] - cloud[:, 0], cloud[:, 0] - box[2]), 0)
                dy = np.maximum(np.maximum(box[1] - cloud[:, 1], cloud[:, 1] - box[3]), 0)
                if np.any(dx * dx + dy * dy < radius * radius):
                    continue
                edge = np.clip(anchors[idx], box[:2], box[2:])
                near = float(np.linalg.norm(edge - anchors[idx]))
                cost = near + 0.035 * math.hypot(cx - px, cy - py)
                data_center = ax.transData.inverted().transform([cx, cy])
                cost += 8 * ((data_center[0] * x < 0) + (data_center[1] * y < 0))
                boxes.append(box)
                costs.append(cost)
                cent.append((cx, cy))
        assert boxes, f"NO_VALID_NEAR_LABEL_POSITION:{row['indicator_id']}"
        order = np.argsort(costs)
        candidates.append(np.array(boxes)[order])
        distances.append(np.array(costs)[order])
        centers.append(np.array(cent)[order])
        texts.append(text)

    def collision(boxes, others):
        if not len(others):
            return np.zeros(len(boxes), dtype=bool)
        return np.any((boxes[:, None, 0] < others[None, :, 2]) & (boxes[:, None, 2] > others[None, :, 0]) & (boxes[:, None, 1] < others[None, :, 3]) & (boxes[:, None, 3] > others[None, :, 1]), axis=1)
    rng = np.random.default_rng(20261030)
    density = [int(np.sum(np.linalg.norm(cloud - p, axis=1) < 75)) for p in anchors]
    orders = [np.argsort([-v for v in density]), np.arange(len(items))]
    orders += [rng.permutation(len(items)) for _ in range(30)]
    best, best_cost = (None, math.inf)
    for order in orders:
        choice = {}
        for idx in order:
            others = np.array([candidates[j][k] for j, k in choice.items()])
            legal = np.flatnonzero(~collision(candidates[idx], others))
            if not len(legal):
                break
            choice[idx] = int(legal[0])
        if len(choice) != len(items):
            continue
        for _ in range(3):
            for idx in rng.permutation(len(items)):
                others = np.array([candidates[j][k] for j, k in choice.items() if j != idx])
                legal = np.flatnonzero(~collision(candidates[idx], others))
                choice[idx] = int(legal[0])
        cost = sum((distances[j][k] for j, k in choice.items()))
        if cost < best_cost:
            best, best_cost = (choice, cost)
    assert best is not None, 'LABEL_LAYOUT_UNRESOLVED_NO_OVERLAP_FREE_ASSIGNMENT'
    placements = []
    for idx, text in enumerate(texts):
        box = candidates[idx][best[idx]]
        center = ax.transData.inverted().transform(centers[idx][best[idx]])
        text.set_position(center)
        bg = '#EDF5FC' if center[0] < center[1] else '#F3F3F3'
        text.set_path_effects([pe.withStroke(linewidth=2.1, foreground=bg)])
        edge_display = np.clip(anchors[idx], box[:2], box[2:])
        edge = ax.transData.inverted().transform(edge_display)
        gap = float(np.linalg.norm(edge_display - anchors[idx]))
        if gap > radius + 2:
            ax.annotate('', xy=items[idx][1:3], xytext=edge, arrowprops={'arrowstyle': '-', 'color': '#777777', 'lw': 0.6, 'shrinkB': 5}, zorder=2.5)
        placements.append({'text': text, 'center': center, 'edge_gap_pixels': gap, 'box': box})
    fig.canvas.draw()
    actual_boxes = [p['text'].get_window_extent(fig.canvas.get_renderer()).expanded(1.02, 1.15) for p in placements]
    assert not any((a.overlaps(b) for i, a in enumerate(actual_boxes) for b in actual_boxes[i + 1:])), 'LABEL_OVERLAP'
    for box in actual_boxes:
        dx = np.maximum(np.maximum(box.x0 - cloud[:, 0], cloud[:, 0] - box.x1), 0)
        dy = np.maximum(np.maximum(box.y0 - cloud[:, 1], cloud[:, 1] - box.y1), 0)
        assert not np.any(dx * dx + dy * dy < (radius - 1.7) ** 2), 'LABEL_POINT_OVERLAP'
    return placements

def main():
    if OUTPUT.exists():
        raise FileExistsError(f'PRESERVE_EXISTING_FIGURE_REVISION: {OUTPUT}')
    tracked = [SOURCE / 'B100_SIX_CONTRASTS_RESULTS.csv', SOURCE / 'B100_121_COMPLETE_RESULTS.csv', SOURCE / 'B100_QUADRANT_UPDATED.pdf', SOURCE / 'B100_QUADRANT_UPDATED.png', SOURCE / 'RELEASE_MANIFEST.json', REGISTRY, OLD_PLOT]
    source_before = {str(path): sha(path) for path in tracked}
    rows = read(tracked[0])
    registry = {r['indicator_id']: r for r in read(REGISTRY)}
    assert len(rows) == 726 and len({r['indicator_id'] for r in rows}) == 121
    points = []
    for row in rows:
        try:
            x, y = (float(row['IE_B_per_1000']), float(row['IE_A_per_1000']))
        except (TypeError, ValueError):
            continue
        if math.isfinite(x) and math.isfinite(y):
            quad = ('P' if x > 0 else 'N') + ('P' if y > 0 else 'N')
            assert x != 0 and y != 0, 'ZERO_COORDINATE_REQUIRES_EXPLICIT_AXIS_CATEGORY'
            points.append((row, x, y, quad))
    assert len(points) == 660
    assert len({r['indicator_id'] for r, _, _, _ in points}) == 110
    nmr = [p for p in points if p[0]['domain'] == 'NMR']
    assert len(nmr) == 318 and len({p[0]['indicator_id'] for p in nmr}) == 53
    assert len({r['indicator_id'] for r in rows if r['domain'] == 'NMR'}) == 54
    counts = collections.Counter(((r['indicator_id'], q) for r, x, y, q in points))
    spec = importlib.util.spec_from_file_location('_prior_plot_names', OLD_PLOT)
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    helper = ROOT.parents[2] / 'scripts/birth/border_fullfollowup_b100_20261027/repairs/cache_only_screen_20261028_v1/plot_screen.py'
    spec = importlib.util.spec_from_file_location('_frozen_names', helper)
    names = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(names)

    def label(row):
        ident = row['indicator_id']
        return DISPLAY.get(ident) or prior.DISPLAY.get(ident) or names.english(row)
    fontManager.addfont('/System/Library/Fonts/Supplemental/Arial.ttf')
    matplotlib.rcParams.update({'font.family': 'Arial', 'font.size': 11, 'text.color': 'black', 'axes.labelcolor': 'black', 'axes.edgecolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black', 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    OUTPUT.mkdir()
    annotations = []
    plots = []

    def plot(dataset, selections, stem):
        limit = math.ceil(max((max(abs(x), abs(y)) for _, x, y, _ in dataset)) * 1.12 * 10) / 10
        fig, ax = plt.subplots(figsize=(8.4, 8.7))
        fig.subplots_adjust(left=0.13, right=0.97, top=0.98, bottom=0.17)
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit))
        ax.set_aspect('equal', adjustable='box')
        ax.add_patch(Polygon([(-limit, -limit), (-limit, limit), (limit, limit)], fc='#EDF5FC', ec='none', zorder=0))
        ax.add_patch(Polygon([(-limit, -limit), (limit, -limit), (limit, limit)], fc='#F3F3F3', ec='none', zorder=0))
        ax.axhline(0, color='#999999', lw=0.9, zorder=1)
        ax.axvline(0, color='#999999', lw=0.9, zorder=1)
        ax.plot([-limit, limit], [-limit, limit], color='#666666', ls=(0, (5, 3)), lw=1.3, zorder=2)
        ordered = sorted(dataset, key=lambda p: (p[0]['domain'] == 'NMR', counts[p[0]['indicator_id'], p[3]]))
        plotted_keys = []
        for row, x, y, quad in ordered:
            strength = STRENGTHS[counts[row['indicator_id'], quad] - 1]
            rgb = np.array(matplotlib.colors.to_rgb(PALETTE[quad]))
            color = 1 + (rgb - 1) * strength
            ax.scatter([x], [y], s=88, c=[color], edgecolors='white', linewidths=0.4, zorder=3)
            plotted_keys.append((row['indicator_id'], row['contrast'], x, y))
        assert collections.Counter(plotted_keys) == collections.Counter(((r['indicator_id'], r['contrast'], x, y) for r, x, y, q in dataset))
        items = []
        for quad, ids in selections.items():
            chosen = []
            for ident in ids:
                candidates = [p for p in dataset if p[0]['indicator_id'] == ident and p[3] == quad]
                assert candidates, (ident, quad)
                chosen.append(max(candidates, key=lambda p: math.hypot(p[1], p[2])))
            chosen.sort(key=lambda p: -p[2])
            items.extend(chosen)
        ax.set_xlabel('DM / T2D indirect risk difference (per 1,000)', fontsize=12.5, fontweight='bold', labelpad=9)
        ax.set_ylabel('IHD / Angina / CAS indirect risk difference (per 1,000)', fontsize=12.5, fontweight='bold', labelpad=9)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
        placements = place_nearby_labels(ax, items, dataset, label)
        for (row, x, y, quad), placement in zip(items, placements):
            annotations.append({'figure': stem, 'indicator_id': row['indicator_id'], 'label': label(row), 'domain': row['domain'], 'quadrant': quad, 'contrast': row['contrast'], 'x_per_1000': x, 'y_per_1000': y, 'label_x': placement['center'][0], 'label_y': placement['center'][1], 'label_edge_gap_pixels_at_100dpi': placement['edge_gap_pixels'], 'registry_label_zh': registry[row['indicator_id']]['label_zh'], 'registry_original_label': registry[row['indicator_id']]['original_label'], 'selection': 'DISPLAY_ONLY_PRESPECIFIED_LABELS_NOT_STATISTICAL_FILTER'})
        key = fig.add_axes([0.34, 0.025, 0.38, 0.028])
        swatches = []
        for q in ('NP', 'PP', 'NN', 'PN'):
            rgb = np.array(matplotlib.colors.to_rgb(PALETTE[q]))
            swatches.append([1 + (rgb - 1) * strength for strength in STRENGTHS])
        key.imshow(swatches, aspect='auto', extent=(0.5, 6.5, 4, 0))
        key.set(yticks=[], xticks=range(1, 7))
        key.tick_params(length=0)
        key.set_title('Directional agreement (of 6)', fontsize=10)
        for ext in ('pdf', 'png', 'svg'):
            target = OUTPUT / f'{stem}.{ext}'
            assert not target.exists()
            fig.savefig(target, dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.06)
        plt.close(fig)
        plots.append({'stem': stem, 'coordinates': len(dataset), 'indicators': len({p[0]['indicator_id'] for p in dataset}), 'NMR_coordinates': sum((p[0]['domain'] == 'NMR' for p in dataset)), 'NMR_indicators': len({p[0]['indicator_id'] for p in dataset if p[0]['domain'] == 'NMR'}), 'NMR_labels': sum((r['figure'] == stem and r['domain'] == 'NMR' for r in annotations)), 'axis_equal_scale': True, 'border_linewidth': 1.5, 'diagonal_linewidth': 1.3, 'coordinate_multiset_exact': True, 'label_overlap': False, 'label_point_overlap': False, 'label_placement': 'LOCAL_DISTANCE_MINIMIZED_NO_FIXED_SLOTS', 'median_label_gap_pixels': float(np.median([p['edge_gap_pixels'] for p in placements])), 'max_label_gap_pixels': max((p['edge_gap_pixels'] for p in placements))})
    plot(points, LABELS_ALL, 'B100_QUADRANT_NMR_VISIBLE_V2')
    plot(nmr, LABELS_NMR, 'B100_NMR_QUADRANT_V2')
    write_csv('ANNOTATION_COORDINATES.csv', annotations)
    nmr_points = []
    for row, x, y, quad in nmr:
        nmr_points.append(dict(row, x_per_1000=x, y_per_1000=y, quadrant=quad))
    write_csv('NMR_PLOTTED_COORDINATES.csv', nmr_points)
    coverage = []
    for ident in sorted((r['indicator_id'] for r in registry.values() if r['category'] == 'NMR')):
        group = [r for r in rows if r['indicator_id'] == ident]
        legal = [p for p in nmr if p[0]['indicator_id'] == ident]
        coverage.append({'indicator_id': ident, 'label_zh': registry[ident]['label_zh'], 'registry_original_label': registry[ident]['original_label'], 'planned_contrasts': 6, 'plotted_contrasts': len(legal), 'statuses': ';'.join(sorted({r['result_status'] for r in group})), 'reason_if_missing': ';'.join(sorted({r['failure_reason'] for r in group if r['failure_reason']}))})
    write_csv('NMR_COVERAGE_AUDIT.csv', coverage)
    source_after = {str(path): sha(path) for path in tracked}
    assert source_before == source_after, 'SOURCE_RELEASE_CHANGED'
    qa = {'time_beijing': datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(), 'plots': plots, 'NMR_registered': 54, 'NMR_plotted': 53, 'NMR_unplotted': [r['indicator_id'] for r in coverage if r['plotted_contrasts'] == 0], 'model_fits': 0, 'bootstrap_runs': 0, 'source_hashes_before': source_before, 'source_hashes_after': source_after, 'source_unchanged': True, 'background': {'x_lt_y': '#EDF5FC', 'x_gt_y': '#F3F3F3'}, 'method_or_inference_changes': False, 'name_correction': {'nmr__nmr_23523': {'old_display': 'IDL TG proportion', 'new_display': 'IDL-P', 'source_registry_label': registry['nmr__nmr_23523']['original_label']}}}
    (OUTPUT / 'FIGURE_REVISION_QA.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding='utf-8')
    note = '# B100四象限图：NMR显示修订V2\n\n只修改绘图；原release、旧图和模型结果均完整保留。全部坐标来自全随访B100结果，11年风险，每千人风险差。\n\n- 原图实际已绘110项、660坐标，其中NMR 53项、318坐标。原按每象限离原点最远的5项标名，NMR名称未被选中；原数据没有遗漏这53项。\n- 主图仍绘全部110项、660坐标；每象限5个名称是展示层代表标签，不按P值筛选。NMR置于点图上层避免被其他领域遮盖。另附只显示53项NMR全部318坐标的专图。\n- 54项NMR中23569（中型HDL胆固醇酯）仍没有合法坐标，保留NA并列入缺失清单，不当作0或无效应。\n- 23523经冻结注册表核实为IDL颗粒浓度，显示IDL-P；旧绘图别名IDL TG proportion有误。真正IDL甘油三酯比例为23613，显示IDL TG%。修正只涉及显示名称，不改指标ID或估计值。\n- 四边框线宽从1增至1.5；y=x虚线线宽从0.65增至1.3并加深。虚线左上侧x<y为浅蓝，右下侧x>y为浅灰；两侧不是显著性分组。\n- 标签改为对应点附近的局部避让，不再按象限边缘排成一列。对全部实际点和所有文字框检查无覆盖，必要时保留短引线；固定可视化随机数只用于排版，不接触分析或Bootstrap随机数。\n- 四象限点颜色保持原绿/橙红/蓝/黄，深浅继续表示同一指标6组比较落在相应象限的次数，不是显著性支持次数。无效应坐标抖动或重新缩放；横纵轴单位和比例相同。\n- 六组比较是IHD、Angina、CAS各自与DM、T2D比较；轴是间接风险差IE，不是a/b系数。不同结局目标人群可能不同，不能把所有点视为同一人群或独立证据。\n- 有警告恢复点仍仅为探索性点值，无Bootstrap支持者不获得CI/P/q。未重算模型、A或Bootstrap，未更改任何推断状态。\n\n`NMR_COVERAGE_AUDIT.csv`列全54项；`NMR_PLOTTED_COORDINATES.csv`列318个实际点；`ANNOTATION_COORDINATES.csv`记录每个标签对应的真实疾病比较和坐标。\n'
    (OUTPUT / 'FIGURE_REVISION_ZH.md').write_text(note, encoding='utf-8')
    manifest = {path.name: {'sha256': sha(path), 'bytes': path.stat().st_size} for path in sorted(OUTPUT.iterdir()) if path.is_file()}
    (OUTPUT / 'FIGURE_MANIFEST.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'output': str(OUTPUT), 'plots': plots, 'source_unchanged': True}, ensure_ascii=False))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Additive, aggregate-only figure revision; never refit or modify source releases."""\nfrom __future__ import annotations\nimport collections\nimport csv\nimport datetime\nimport hashlib\nimport importlib.util\nimport json\nimport math\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib.font_manager import fontManager\nfrom matplotlib.patches import Polygon\nimport matplotlib.patheffects as pe\nimport numpy as np\nROOT = Path(__file__).resolve().parent\nSOURCE = ROOT / \'final_release\'\nOUTPUT = ROOT / \'figure_revision_nmr_v4_near_labels\'\nREGISTRY = ROOT.parents[2] / \'scripts/birth/border_fullfollowup_b100_20261027/input/APPROVED_121_NO_CMR.csv\'\nOLD_PLOT = ROOT / \'plot_and_report.py\'\nPALETTE = {\'NP\': \'#009E73\', \'PP\': \'#D55E00\', \'NN\': \'#0072B2\', \'PN\': \'#E69F00\'}\nSTRENGTHS = (0.2, 0.36, 0.52, 0.68, 0.84, 1.0)\nLABELS_ALL = {\'NP\': [\'olink__fcrl6\', \'clinical__glucose\', \'nmr__nmr_23404\', \'nmr__nmr_23523\', \'nmr__nmr_23506\'], \'PP\': [\'qualification\', \'smoking\', \'age_education\', \'income\', \'nmr__nmr_23451\'], \'NN\': [\'birth_weight\', \'bmi\', \'olink__gzma\', \'nmr__nmr_23613\', \'nmr__nmr_23406\'], \'PN\': [\'employment_employed\', \'participation_religious\', \'employment_voluntary\', \'confiding\', \'employment_student\']}\nLABELS_NMR = {\'NP\': [\'nmr__nmr_23404\', \'nmr__nmr_23523\', \'nmr__nmr_23506\', \'nmr__nmr_23526\', \'nmr__nmr_23584\'], \'PP\': [\'nmr__nmr_23451\'], \'NN\': [\'nmr__nmr_23613\', \'nmr__nmr_23406\', \'nmr__nmr_23440\', \'nmr__nmr_23608\', \'nmr__nmr_23648\'], \'PN\': []}\nDISPLAY = {\'nmr__nmr_23404\': \'NMR LDL-C\', \'nmr__nmr_23406\': \'NMR HDL-C\', \'nmr__nmr_23440\': \'NMR ApoA1\', \'nmr__nmr_23451\': \'Omega-3 FA%\', \'nmr__nmr_23506\': \'M-VLDL CE\', \'nmr__nmr_23523\': \'IDL-P\', \'nmr__nmr_23526\': \'IDL-C\', \'nmr__nmr_23584\': \'Very large VLDL PL%\', \'nmr__nmr_23608\': \'XS-VLDL TG%\', \'nmr__nmr_23613\': \'IDL TG%\', \'nmr__nmr_23648\': \'Small HDL TG%\', \'employment_employed\': \'Employed / self-employed\', \'employment_student\': \'Student\', \'employment_voluntary\': \'Unpaid / voluntary work\', \'participation_religious\': \'Religious group\', \'confiding\': \'Confiding frequency\'}\n\ndef read(path):\n    with path.open(encoding=\'utf-8-sig\', newline=\'\') as handle:\n        return list(csv.DictReader(handle))\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef write_csv(name, rows):\n    with (OUTPUT / name).open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as handle:\n        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))\n        writer.writeheader()\n        writer.writerows(rows)\n\ndef place_nearby_labels(ax, items, dataset, label):\n    """Bounded display-space placement; avoid every point and every other label.\n\n    The fixed RNG is for visual layout only, never analysis or Bootstrap.\n    Candidates are scored by point-to-text-edge distance, not fixed plot slots.\n    """\n    fig = ax.figure\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    extent = ax.get_window_extent(renderer)\n    cloud = ax.transData.transform([(x, y) for _, x, y, _ in dataset])\n    anchors = ax.transData.transform([(x, y) for _, x, y, _ in items])\n    radius = math.sqrt(88) * fig.dpi / 72 / 2 + 2.4\n    candidates, distances, centers = ([], [], [])\n    texts = []\n    for idx, (row, x, y, quad) in enumerate(items):\n        text = ax.text(x, y, label(row), ha=\'center\', va=\'center\', fontsize=14.5, color=\'black\', zorder=6)\n        bbox = text.get_window_extent(renderer)\n        width, height = (bbox.width + 6, bbox.height + 6)\n        boxes, costs, cent = ([], [], [])\n        px, py = anchors[idx]\n        for gap in (9, 14, 20, 28, 38, 50, 66, 86, 110, 140, 175, 215):\n            for angle in np.linspace(0, 2 * math.pi, 40, endpoint=False):\n                cx = px + math.cos(angle) * (width / 2 + gap)\n                cy = py + math.sin(angle) * (height / 2 + gap)\n                box = np.array([cx - width / 2, cy - height / 2, cx + width / 2, cy + height / 2])\n                if box[0] < extent.x0 + 4 or box[2] > extent.x1 - 4 or box[1] < extent.y0 + 4 or (box[3] > extent.y1 - 4):\n                    continue\n                dx = np.maximum(np.maximum(box[0] - cloud[:, 0], cloud[:, 0] - box[2]), 0)\n                dy = np.maximum(np.maximum(box[1] - cloud[:, 1], cloud[:, 1] - box[3]), 0)\n                if np.any(dx * dx + dy * dy < radius * radius):\n                    continue\n                edge = np.clip(anchors[idx], box[:2], box[2:])\n                near = float(np.linalg.norm(edge - anchors[idx]))\n                cost = near + 0.035 * math.hypot(cx - px, cy - py)\n                data_center = ax.transData.inverted().transform([cx, cy])\n                cost += 8 * ((data_center[0] * x < 0) + (data_center[1] * y < 0))\n                boxes.append(box)\n                costs.append(cost)\n                cent.append((cx, cy))\n        assert boxes, f"NO_VALID_NEAR_LABEL_POSITION:{row[\'indicator_id\']}"\n        order = np.argsort(costs)\n        candidates.append(np.array(boxes)[order])\n        distances.append(np.array(costs)[order])\n        centers.append(np.array(cent)[order])\n        texts.append(text)\n\n    def collision(boxes, others):\n        if not len(others):\n            return np.zeros(len(boxes), dtype=bool)\n        return np.any((boxes[:, None, 0] < others[None, :, 2]) & (boxes[:, None, 2] > others[None, :, 0]) & (boxes[:, None, 1] < others[None, :, 3]) & (boxes[:, None, 3] > others[None, :, 1]), axis=1)\n    rng = np.random.default_rng(20261030)\n    density = [int(np.sum(np.linalg.norm(cloud - p, axis=1) < 75)) for p in anchors]\n    orders = [np.argsort([-v for v in density]), np.arange(len(items))]\n    orders += [rng.permutation(len(items)) for _ in range(30)]\n    best, best_cost = (None, math.inf)\n    for order in orders:\n        choice = {}\n        for idx in order:\n            others = np.array([candidates[j][k] for j, k in choice.items()])\n            legal = np.flatnonzero(~collision(candidates[idx], others))\n            if not len(legal):\n                break\n            choice[idx] = int(legal[0])\n        if len(choice) != len(items):\n            continue\n        for _ in range(3):\n            for idx in rng.permutation(len(items)):\n                others = np.array([candidates[j][k] for j, k in choice.items() if j != idx])\n                legal = np.flatnonzero(~collision(candidates[idx], others))\n                choice[idx] = int(legal[0])\n        cost = sum((distances[j][k] for j, k in choice.items()))\n        if cost < best_cost:\n            best, best_cost = (choice, cost)\n    assert best is not None, \'LABEL_LAYOUT_UNRESOLVED_NO_OVERLAP_FREE_ASSIGNMENT\'\n    placements = []\n    for idx, text in enumerate(texts):\n        box = candidates[idx][best[idx]]\n        center = ax.transData.inverted().transform(centers[idx][best[idx]])\n        text.set_position(center)\n        bg = \'#EDF5FC\' if center[0] < center[1] else \'#F3F3F3\'\n        text.set_path_effects([pe.withStroke(linewidth=2.1, foreground=bg)])\n        edge_display = np.clip(anchors[idx], box[:2], box[2:])\n        edge = ax.transData.inverted().transform(edge_display)\n        gap = float(np.linalg.norm(edge_display - anchors[idx]))\n        if gap > radius + 2:\n            ax.annotate(\'\', xy=items[idx][1:3], xytext=edge, arrowprops={\'arrowstyle\': \'-\', \'color\': \'#777777\', \'lw\': 0.6, \'shrinkB\': 5}, zorder=2.5)\n        placements.append({\'text\': text, \'center\': center, \'edge_gap_pixels\': gap, \'box\': box})\n    fig.canvas.draw()\n    actual_boxes = [p[\'text\'].get_window_extent(fig.canvas.get_renderer()).expanded(1.02, 1.15) for p in placements]\n    assert not any((a.overlaps(b) for i, a in enumerate(actual_boxes) for b in actual_boxes[i + 1:])), \'LABEL_OVERLAP\'\n    for box in actual_boxes:\n        dx = np.maximum(np.maximum(box.x0 - cloud[:, 0], cloud[:, 0] - box.x1), 0)\n        dy = np.maximum(np.maximum(box.y0 - cloud[:, 1], cloud[:, 1] - box.y1), 0)\n        assert not np.any(dx * dx + dy * dy < (radius - 1.7) ** 2), \'LABEL_POINT_OVERLAP\'\n    return placements\n\ndef main():\n    if OUTPUT.exists():\n        raise FileExistsError(f\'PRESERVE_EXISTING_FIGURE_REVISION: {OUTPUT}\')\n    tracked = [SOURCE / \'B100_SIX_CONTRASTS_RESULTS.csv\', SOURCE / \'B100_121_COMPLETE_RESULTS.csv\', SOURCE / \'B100_QUADRANT_UPDATED.pdf\', SOURCE / \'B100_QUADRANT_UPDATED.png\', SOURCE / \'RELEASE_MANIFEST.json\', REGISTRY, OLD_PLOT]\n    source_before = {str(path): sha(path) for path in tracked}\n    rows = read(tracked[0])\n    registry = {r[\'indicator_id\']: r for r in read(REGISTRY)}\n    assert len(rows) == 726 and len({r[\'indicator_id\'] for r in rows}) == 121\n    points = []\n    for row in rows:\n        try:\n            x, y = (float(row[\'IE_B_per_1000\']), float(row[\'IE_A_per_1000\']))\n        except (TypeError, ValueError):\n            continue\n        if math.isfinite(x) and math.isfinite(y):\n            quad = (\'P\' if x > 0 else \'N\') + (\'P\' if y > 0 else \'N\')\n            assert x != 0 and y != 0, \'ZERO_COORDINATE_REQUIRES_EXPLICIT_AXIS_CATEGORY\'\n            points.append((row, x, y, quad))\n    assert len(points) == 660\n    assert len({r[\'indicator_id\'] for r, _, _, _ in points}) == 110\n    nmr = [p for p in points if p[0][\'domain\'] == \'NMR\']\n    assert len(nmr) == 318 and len({p[0][\'indicator_id\'] for p in nmr}) == 53\n    assert len({r[\'indicator_id\'] for r in rows if r[\'domain\'] == \'NMR\'}) == 54\n    counts = collections.Counter(((r[\'indicator_id\'], q) for r, x, y, q in points))\n    spec = importlib.util.spec_from_file_location(\'_prior_plot_names\', OLD_PLOT)\n    prior = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(prior)\n    helper = ROOT.parents[2] / \'scripts/birth/border_fullfollowup_b100_20261027/repairs/cache_only_screen_20261028_v1/plot_screen.py\'\n    spec = importlib.util.spec_from_file_location(\'_frozen_names\', helper)\n    names = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(names)\n\n    def label(row):\n        ident = row[\'indicator_id\']\n        return DISPLAY.get(ident) or prior.DISPLAY.get(ident) or names.english(row)\n    fontManager.addfont(\'/System/Library/Fonts/Supplemental/Arial.ttf\')\n    matplotlib.rcParams.update({\'font.family\': \'Arial\', \'font.size\': 11, \'text.color\': \'black\', \'axes.labelcolor\': \'black\', \'axes.edgecolor\': \'black\', \'xtick.color\': \'black\', \'ytick.color\': \'black\', \'figure.facecolor\': \'white\', \'axes.facecolor\': \'white\', \'pdf.fonttype\': 42, \'ps.fonttype\': 42, \'svg.fonttype\': \'none\'})\n    OUTPUT.mkdir()\n    annotations = []\n    plots = []\n\n    def plot(dataset, selections, stem):\n        limit = math.ceil(max((max(abs(x), abs(y)) for _, x, y, _ in dataset)) * 1.12 * 10) / 10\n        fig, ax = plt.subplots(figsize=(8.4, 8.7))\n        fig.subplots_adjust(left=0.13, right=0.97, top=0.98, bottom=0.17)\n        ax.set(xlim=(-limit, limit), ylim=(-limit, limit))\n        ax.set_aspect(\'equal\', adjustable=\'box\')\n        ax.add_patch(Polygon([(-limit, -limit), (-limit, limit), (limit, limit)], fc=\'#EDF5FC\', ec=\'none\', zorder=0))\n        ax.add_patch(Polygon([(-limit, -limit), (limit, -limit), (limit, limit)], fc=\'#F3F3F3\', ec=\'none\', zorder=0))\n        ax.axhline(0, color=\'#999999\', lw=0.9, zorder=1)\n        ax.axvline(0, color=\'#999999\', lw=0.9, zorder=1)\n        ax.plot([-limit, limit], [-limit, limit], color=\'#666666\', ls=(0, (5, 3)), lw=1.3, zorder=2)\n        ordered = sorted(dataset, key=lambda p: (p[0][\'domain\'] == \'NMR\', counts[p[0][\'indicator_id\'], p[3]]))\n        plotted_keys = []\n        for row, x, y, quad in ordered:\n            strength = STRENGTHS[counts[row[\'indicator_id\'], quad] - 1]\n            rgb = np.array(matplotlib.colors.to_rgb(PALETTE[quad]))\n            color = 1 + (rgb - 1) * strength\n            ax.scatter([x], [y], s=88, c=[color], edgecolors=\'white\', linewidths=0.4, zorder=3)\n            plotted_keys.append((row[\'indicator_id\'], row[\'contrast\'], x, y))\n        assert collections.Counter(plotted_keys) == collections.Counter(((r[\'indicator_id\'], r[\'contrast\'], x, y) for r, x, y, q in dataset))\n        items = []\n        for quad, ids in selections.items():\n            chosen = []\n            for ident in ids:\n                candidates = [p for p in dataset if p[0][\'indicator_id\'] == ident and p[3] == quad]\n                assert candidates, (ident, quad)\n                chosen.append(max(candidates, key=lambda p: math.hypot(p[1], p[2])))\n            chosen.sort(key=lambda p: -p[2])\n            items.extend(chosen)\n        ax.set_xlabel(\'DM / T2D indirect risk difference (per 1,000)\', fontsize=12.5, fontweight=\'bold\', labelpad=9)\n        ax.set_ylabel(\'IHD / Angina / CAS indirect risk difference (per 1,000)\', fontsize=12.5, fontweight=\'bold\', labelpad=9)\n        for spine in ax.spines.values():\n            spine.set_visible(True)\n            spine.set_linewidth(1.5)\n        placements = place_nearby_labels(ax, items, dataset, label)\n        for (row, x, y, quad), placement in zip(items, placements):\n            annotations.append({\'figure\': stem, \'indicator_id\': row[\'indicator_id\'], \'label\': label(row), \'domain\': row[\'domain\'], \'quadrant\': quad, \'contrast\': row[\'contrast\'], \'x_per_1000\': x, \'y_per_1000\': y, \'label_x\': placement[\'center\'][0], \'label_y\': placement[\'center\'][1], \'label_edge_gap_pixels_at_100dpi\': placement[\'edge_gap_pixels\'], \'registry_label_zh\': registry[row[\'indicator_id\']][\'label_zh\'], \'registry_original_label\': registry[row[\'indicator_id\']][\'original_label\'], \'selection\': \'DISPLAY_ONLY_PRESPECIFIED_LABELS_NOT_STATISTICAL_FILTER\'})\n        key = fig.add_axes([0.34, 0.025, 0.38, 0.028])\n        swatches = []\n        for q in (\'NP\', \'PP\', \'NN\', \'PN\'):\n            rgb = np.array(matplotlib.colors.to_rgb(PALETTE[q]))\n            swatches.append([1 + (rgb - 1) * strength for strength in STRENGTHS])\n        key.imshow(swatches, aspect=\'auto\', extent=(0.5, 6.5, 4, 0))\n        key.set(yticks=[], xticks=range(1, 7))\n        key.tick_params(length=0)\n        key.set_title(\'Directional agreement (of 6)\', fontsize=10)\n        for ext in (\'pdf\', \'png\', \'svg\'):\n            target = OUTPUT / f\'{stem}.{ext}\'\n            assert not target.exists()\n            fig.savefig(target, dpi=300, facecolor=\'white\', bbox_inches=\'tight\', pad_inches=0.06)\n        plt.close(fig)\n        plots.append({\'stem\': stem, \'coordinates\': len(dataset), \'indicators\': len({p[0][\'indicator_id\'] for p in dataset}), \'NMR_coordinates\': sum((p[0][\'domain\'] == \'NMR\' for p in dataset)), \'NMR_indicators\': len({p[0][\'indicator_id\'] for p in dataset if p[0][\'domain\'] == \'NMR\'}), \'NMR_labels\': sum((r[\'figure\'] == stem and r[\'domain\'] == \'NMR\' for r in annotations)), \'axis_equal_scale\': True, \'border_linewidth\': 1.5, \'diagonal_linewidth\': 1.3, \'coordinate_multiset_exact\': True, \'label_overlap\': False, \'label_point_overlap\': False, \'label_placement\': \'LOCAL_DISTANCE_MINIMIZED_NO_FIXED_SLOTS\', \'median_label_gap_pixels\': float(np.median([p[\'edge_gap_pixels\'] for p in placements])), \'max_label_gap_pixels\': max((p[\'edge_gap_pixels\'] for p in placements))})\n    plot(points, LABELS_ALL, \'B100_QUADRANT_NMR_VISIBLE_V2\')\n    plot(nmr, LABELS_NMR, \'B100_NMR_QUADRANT_V2\')\n    write_csv(\'ANNOTATION_COORDINATES.csv\', annotations)\n    nmr_points = []\n    for row, x, y, quad in nmr:\n        nmr_points.append(dict(row, x_per_1000=x, y_per_1000=y, quadrant=quad))\n    write_csv(\'NMR_PLOTTED_COORDINATES.csv\', nmr_points)\n    coverage = []\n    for ident in sorted((r[\'indicator_id\'] for r in registry.values() if r[\'category\'] == \'NMR\')):\n        group = [r for r in rows if r[\'indicator_id\'] == ident]\n        legal = [p for p in nmr if p[0][\'indicator_id\'] == ident]\n        coverage.append({\'indicator_id\': ident, \'label_zh\': registry[ident][\'label_zh\'], \'registry_original_label\': registry[ident][\'original_label\'], \'planned_contrasts\': 6, \'plotted_contrasts\': len(legal), \'statuses\': \';\'.join(sorted({r[\'result_status\'] for r in group})), \'reason_if_missing\': \';\'.join(sorted({r[\'failure_reason\'] for r in group if r[\'failure_reason\']}))})\n    write_csv(\'NMR_COVERAGE_AUDIT.csv\', coverage)\n    source_after = {str(path): sha(path) for path in tracked}\n    assert source_before == source_after, \'SOURCE_RELEASE_CHANGED\'\n    qa = {\'time_beijing\': datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(), \'plots\': plots, \'NMR_registered\': 54, \'NMR_plotted\': 53, \'NMR_unplotted\': [r[\'indicator_id\'] for r in coverage if r[\'plotted_contrasts\'] == 0], \'model_fits\': 0, \'bootstrap_runs\': 0, \'source_hashes_before\': source_before, \'source_hashes_after\': source_after, \'source_unchanged\': True, \'background\': {\'x_lt_y\': \'#EDF5FC\', \'x_gt_y\': \'#F3F3F3\'}, \'method_or_inference_changes\': False, \'name_correction\': {\'nmr__nmr_23523\': {\'old_display\': \'IDL TG proportion\', \'new_display\': \'IDL-P\', \'source_registry_label\': registry[\'nmr__nmr_23523\'][\'original_label\']}}}\n    (OUTPUT / \'FIGURE_REVISION_QA.json\').write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding=\'utf-8\')\n    note = \'# B100四象限图：NMR显示修订V2\\n\\n只修改绘图；原release、旧图和模型结果均完整保留。全部坐标来自全随访B100结果，11年风险，每千人风险差。\\n\\n- 原图实际已绘110项、660坐标，其中NMR 53项、318坐标。原按每象限离原点最远的5项标名，NMR名称未被选中；原数据没有遗漏这53项。\\n- 主图仍绘全部110项、660坐标；每象限5个名称是展示层代表标签，不按P值筛选。NMR置于点图上层避免被其他领域遮盖。另附只显示53项NMR全部318坐标的专图。\\n- 54项NMR中23569（中型HDL胆固醇酯）仍没有合法坐标，保留NA并列入缺失清单，不当作0或无效应。\\n- 23523经冻结注册表核实为IDL颗粒浓度，显示IDL-P；旧绘图别名IDL TG proportion有误。真正IDL甘油三酯比例为23613，显示IDL TG%。修正只涉及显示名称，不改指标ID或估计值。\\n- 四边框线宽从1增至1.5；y=x虚线线宽从0.65增至1.3并加深。虚线左上侧x<y为浅蓝，右下侧x>y为浅灰；两侧不是显著性分组。\\n- 标签改为对应点附近的局部避让，不再按象限边缘排成一列。对全部实际点和所有文字框检查无覆盖，必要时保留短引线；固定可视化随机数只用于排版，不接触分析或Bootstrap随机数。\\n- 四象限点颜色保持原绿/橙红/蓝/黄，深浅继续表示同一指标6组比较落在相应象限的次数，不是显著性支持次数。无效应坐标抖动或重新缩放；横纵轴单位和比例相同。\\n- 六组比较是IHD、Angina、CAS各自与DM、T2D比较；轴是间接风险差IE，不是a/b系数。不同结局目标人群可能不同，不能把所有点视为同一人群或独立证据。\\n- 有警告恢复点仍仅为探索性点值，无Bootstrap支持者不获得CI/P/q。未重算模型、A或Bootstrap，未更改任何推断状态。\\n\\n`NMR_COVERAGE_AUDIT.csv`列全54项；`NMR_PLOTTED_COORDINATES.csv`列318个实际点；`ANNOTATION_COORDINATES.csv`记录每个标签对应的真实疾病比较和坐标。\\n\'\n    (OUTPUT / \'FIGURE_REVISION_ZH.md\').write_text(note, encoding=\'utf-8\')\n    manifest = {path.name: {\'sha256\': sha(path), \'bytes\': path.stat().st_size} for path in sorted(OUTPUT.iterdir()) if path.is_file()}\n    (OUTPUT / \'FIGURE_MANIFEST.json\').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding=\'utf-8\')\n    print(json.dumps({\'output\': str(OUTPUT), \'plots\': plots, \'source_unchanged\': True}, ensure_ascii=False))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
