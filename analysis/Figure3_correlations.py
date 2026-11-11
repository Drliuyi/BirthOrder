"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import collections
import csv
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import sys
import zipfile
from zoneinfo import ZoneInfo
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.font_manager import fontManager
from matplotlib.patches import Rectangle, Patch
import numpy as np
DOM = {'危险因素': 'Clinical / socio-behavioural', '社会、行为及生命历程': 'Clinical / socio-behavioural', 'NMR': 'NMR', '蛋白组': 'Olink'}
COLORS = {'Clinical / socio-behavioural': '#7FA3B6', 'NMR': '#8CAB9F', 'Olink': '#A899B7'}

def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, rows):
    assert rows
    with p.open('x', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys((k for r in rows for k in r))))
        w.writeheader()
        w.writerows(rows)

def js(p, x):
    with p.open('x', encoding='utf-8') as f:
        json.dump(x, f, ensure_ascii=False, indent=2, allow_nan=False)

def text(p, s):
    with p.open('x', encoding='utf-8') as f:
        f.write(s)

def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None

def matrixwrite(path, ids, matrix):
    write(path, [{'indicator_id': i, **{j: matrix[i][j] for j in ids}} for i in ids])

def hierarchy(ids, rho):
    """Exact deterministic agglomerative complete linkage; no scipy required."""
    ids = list(ids)
    n = len(ids)
    members = {k: (i,) for k, i in enumerate(ids)}
    active = list(range(n))
    merges = []
    for step in range(n - 1):
        candidates = []
        for a, b in itertools.combinations(active, 2):
            distance = max((1 - abs(rho[i][j]) for i in members[a] for j in members[b]))
            candidates.append((distance, tuple(sorted(members[a] + members[b])), a, b))
        dist, _, a, b = min(candidates)
        if members[b] < members[a]:
            a, b = (b, a)
        node = n + step
        members[node] = tuple(sorted(members[a] + members[b]))
        merges.append((a, b, max(0.0, dist), len(members[node])))
        active = [v for v in active if v not in (a, b)] + [node]

    def leaves(node):
        if node < n:
            return [ids[node]]
        a, b, _, _ = merges[node - n]
        return leaves(a) + leaves(b)

    def cut(node):
        if node < n:
            return [(ids[node],)]
        a, b, h, _ = merges[node - n]
        if h <= 0.1 + 1e-14:
            return [members[node]]
        return cut(a) + cut(b)
    root = active[0]
    return {'ids': ids, 'merges': merges, 'order': leaves(root), 'groups': cut(root), 'members': members}

def tree(ax, h, orientation):
    n = len(h['ids'])
    position = {i: k for k, i in enumerate(h['order'])}
    xs = {k: position[i] for k, i in enumerate(h['ids'])}
    ys = {k: 0.0 for k in range(n)}
    for k, (a, b, height, _) in enumerate(h['merges']):
        xx = [xs[a], xs[a], xs[b], xs[b]]
        yy = [ys[a], height, height, ys[b]]
        ax.plot(xx if orientation == 'top' else yy, yy if orientation == 'top' else xx, color='#414141', lw=0.65)
        xs[n + k] = (xs[a] + xs[b]) / 2
        ys[n + k] = height
    limit = max([x[2] for x in h['merges']] + [0.15]) * 1.04
    if orientation == 'top':
        ax.set(xlim=(-0.5, n - 0.5), ylim=(0, limit))
        ax.axhline(0.1, color='#B75871', lw=0.65, ls=(0, (3, 2)))
    else:
        ax.set(xlim=(limit, 0), ylim=(n - 0.5, -0.5))
        ax.axvline(0.1, color='#B75871', lw=0.65, ls=(0, (3, 2)))
    ax.axis('off')

def heatmap(out, stem, ids, rho, labels, domains, clustered=False, annotate=False):
    h = hierarchy(ids, rho)
    order = h['order'] if clustered else list(ids)
    n = len(order)
    cell = 0.235 if n > 35 else 0.275 if n > 20 else 0.41 if n > 6 else 0.78
    side = n * cell
    dendro = 0.9 if clustered else 0.0
    left_labels = 2.5 if n > 20 else 2.12
    bottom = 2.46 if n > 20 else 2.16 if n > 6 else 1.23
    if n <= 6:
        left_labels = 0.88
    top = 0.62
    right_labels = 0 if not clustered else left_labels
    left = 0.3 + dendro if clustered else left_labels
    right = 0.6 + right_labels
    fw = left + side + right
    fh = bottom + side + dendro + top
    fig = plt.figure(figsize=(fw, fh), dpi=100)
    ax = fig.add_axes([left / fw, bottom / fh, side / fw, side / fh])
    mat = np.array([[rho[i][j] for j in order] for i in order])
    im = ax.imshow(mat, cmap=CMAP, norm=NORM, interpolation='nearest', origin='upper')
    ticks = np.arange(n)
    names = [labels[i] for i in order]
    ax.set_xticks(ticks, labels=names, rotation=90, fontsize=8.0 if n > 20 else 9)
    ax.set_yticks(ticks, labels=names, fontsize=8.0 if n > 20 else 9)
    ax.tick_params(length=0, pad=4)
    if clustered:
        ax.yaxis.tick_right()
    for spine in ax.spines.values():
        spine.set(color='#404040', linewidth=0.75)
    bar = fig.add_axes([left / fw, (bottom + side + 0.025) / fh, side / fw, 0.065 / fh])
    bar.imshow(np.array([[list(COLORS).index(domains[i]) for i in order]]), cmap=ListedColormap(list(COLORS.values())), vmin=0, vmax=2, aspect='auto', interpolation='nearest')
    bar.axis('off')
    if clustered:
        ta = fig.add_axes([left / fw, (bottom + side + 0.13) / fh, side / fw, dendro / fh])
        tree(ta, h, 'top')
        la = fig.add_axes([0.2 / fw, bottom / fh, dendro / fw, side / fh])
        tree(la, h, 'left')
        for group in h['groups']:
            if len(group) > 1:
                positions = [order.index(i) for i in group]
                start = min(positions)
                assert sorted(positions) == list(range(start, start + len(group)))
                ax.add_patch(Rectangle((start - 0.5, start - 0.5), len(group), len(group), fc='none', ec='#444444', lw=0.75))
    else:
        for k in range(1, n):
            if domains[order[k]] != domains[order[k - 1]]:
                ax.axhline(k - 0.5, color='#666666', lw=0.8)
                ax.axvline(k - 0.5, color='#666666', lw=0.8)
    if annotate:
        for y in range(n):
            for x in range(n):
                ax.text(x, y, f'{mat[y, x]:.2f}', ha='center', va='center', fontsize=9, color='white' if abs(mat[y, x]) > 0.65 else 'black')
    cbax = fig.add_axes([left / fw, 0.37 / fh, min(side, 4.5) / fw, 0.13 / fh])
    cb = fig.colorbar(im, cax=cbax, orientation='horizontal', ticks=[-1, -0.5, 0, 0.5, 1])
    cb.ax.tick_params(labelsize=8, length=2)
    cb.outline.set_linewidth(0.5)
    cb.set_label('Spearman ρ', fontsize=9, labelpad=3)
    title = {'CORR_HEATMAP_ALL_44': 'Selected mediators (n = 44)', 'CORR_HEATMAP_ALL_44_CLUSTERED': 'Selected mediators (n = 44)', 'CORR_HEATMAP_CLINICAL_SOCIAL': 'Clinical and behavioural mediators (n = 11)', 'CORR_HEATMAP_NMR': 'NMR metabolomic mediators (n = 29)', 'CORR_HEATMAP_PROTEIN': 'Olink proteomic mediators (n = 4)'}[stem]
    fig.text(left / fw, (fh - 0.24) / fh, title, ha='left', fontsize=11, color='black')
    if n == 44:
        fig.legend(handles=[Patch(facecolor=c, label=k) for k, c in COLORS.items()], loc='upper right', bbox_to_anchor=(0.995, 0.998), ncol=3, fontsize=7.2, frameon=False, handlelength=0.8, columnspacing=0.9)
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    for a in fig.axes:
        if not a.axison:
            continue
        for q in [*a.get_xticklabels(), *a.get_yticklabels()]:
            if not q.get_visible() or not q.get_text():
                continue
            box = q.get_window_extent(ren)
            fb = fig.bbox
            assert fb.contains(box.x0, box.y0) and fb.contains(box.x1, box.y1), ('TEXT_CLIPPED', stem, q.get_text())
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(out / (stem + '.' + ext), dpi=320, facecolor='white')
    fig.savefig(out / (stem + '_preview.png'), dpi=110, facecolor='white')
    plt.close(fig)
    return {'stem': stem, 'n': n, 'clustered': clustered, 'order': order, 'size_inches': [fw, fh], 'text_bounds_passed': True}

def main(args):
    src = args.master
    corrpath = args.correlation
    out = args.out
    if out is None:
        v = 1
        while (BASE / f'reports/birth/b100_correlation44_20261030_v{v}').exists():
            v += 1
        out = BASE / f'reports/birth/b100_correlation44_20261030_v{v}'
    assert not out.exists(), 'ADDITIVE_OUTPUT_ONLY'
    required = ['B100_121_IHD_T2D_MASTER.csv', 'B100_121_ALL_DISEASES_MASTER.csv', 'B100_121_COMPLETE_RESULTS.csv', 'B100_SIX_CONTRASTS_RESULTS.csv', 'PRIORITY_NAMED_CURRENT_RECORDS.csv', 'REVIEW_REPORT_ZH.md']
    input_files = [src / f for f in required] + [src / 'FINAL_MANIFEST.json', src / 'SUMMARY_AND_PROVENANCE.json', corrpath, src / 'evidence/NMR_BIOCHEMICAL_CLASSIFICATION.csv', BASE / 'scripts/birth/border_fullfollowup_b100_20261027/input/APPROVED_121_NO_CMR.csv', BASE / 'reports/birth/figures_b100_pm_4domains_20261030_v2/LABEL_DICTIONARY.csv']
    for p in input_files:
        assert p.is_file(), str(p)
    before = {str(p): sha(p) for p in input_files}
    manifest = json.loads((src / 'FINAL_MANIFEST.json').read_text())
    for name in required:
        assert sha(src / name) == manifest[name]['sha256']
    masters = read(src / required[0])
    wide = read(src / required[1])
    single = read(src / required[2])
    pairs = read(src / required[3])
    priority = read(src / required[4])
    report_text = (src / required[5]).read_text()
    assert len(masters) == len(wide) == 121 and len(single) == 605 and (len(pairs) == 726) and (len(priority) == 12) and (len(report_text) > 1000)
    registry = read(input_files[-2])
    assert {r['indicator_id'] for r in registry} == {r['indicator_id'] for r in masters}
    s = {(r['indicator_id'], r['endpoint_id']): r for r in single}
    p = {r['indicator_id']: r for r in pairs if r['endpoint_A'] == '411' and r['endpoint_B'] == '250.2'}
    assert len(s) == 605 and len(p) == 121
    audit = []
    selected = []
    for m in masters:
        i = m['indicator_id']
        a = s[i, '250.2']
        b = s[i, '411']
        q = [num(a['q_BH_global605_exploratory']), num(b['q_BH_global605_exploratory'])]
        delta = num(p[i]['Delta_IE'])
        dq = num(p[i]['q_BH_within_contrast121_exploratory'])
        pm = [num(a['PM_percent']), num(b['PM_percent'])]
        c1 = any((v is not None and v < 0.05 for v in q))
        c2 = delta is not None and delta > 0
        c3 = dq is not None and dq < 0.05
        c4 = any((v is not None and v > 5 for v in pm))
        first = c1 and c2 and c3
        final = first and c4
        assert first == (m['all_three_conditions'] == 'TRUE')
        audit.append({**m, 'condition1_IE_q605_lt005': c1, 'condition2_Delta_IE_positive': c2, 'condition3_Delta_IE_q121_lt005': c3, 'condition4_T2D_or_IHD_PM_gt5_not_absolute': c4, 'selected_44': final})
        if final:
            selected.append(m)
    assert sum((r['all_three_conditions'] == 'TRUE' for r in masters)) == 64 and len(selected) == 44
    assert collections.Counter((r['domain'] for r in selected)) == {'危险因素': 11, 'NMR': 29, '蛋白组': 4}
    ids = [r['indicator_id'] for r in selected]
    m = {r['indicator_id']: r for r in masters}
    domains = {i: DOM[m[i]['domain']] for i in ids}
    dictionary = {r['indicator_id']: r['label_en'] for r in read(input_files[-1])}
    labels = {i: dictionary.get(i, s[i, '411']['label_en']) for i in ids}
    labels.update({'smoking': 'Smoking status', 'clinical__hdl': 'Clinical HDL-C', 'clinical__glucose': 'Glucose', 'clinical__hba1c': 'HbA1c', 'clinical__tg': 'Triglycerides'})
    evidence = json.loads(corrpath.read_text())
    assert evidence['ids'] == ids and evidence['models_fit'] == evidence['bootstrap_added'] == evidence['participant_rows_exported'] == 0
    for path, digest in evidence['provenance_sha256'].items():
        assert sha(Path(path)) == digest
    rho = evidence['rho']
    nn = evidence['n_pair']
    mdata = evidence['metadata']
    mat = np.array([[rho[i][j] for j in ids] for i in ids])
    nmat = np.array([[nn[i][j] for j in ids] for i in ids])
    assert np.all(np.isfinite(mat)) and np.all(np.abs(mat) <= 1 + 1e-14) and np.allclose(mat, mat.T, rtol=0, atol=1e-14)
    assert np.allclose(np.diag(mat), 1, rtol=0, atol=1e-14) and np.array_equal(nmat, nmat.T) and np.all(nmat > 2)
    assert all((nn[i][j] <= min(nn[i][i], nn[j][j]) for i in ids for j in ids))
    assert all((nn[i][i] == mdata[i]['n_finite_applicable'] for i in ids))
    bio = {r['indicator_id']: r for r in read(src / 'evidence/NMR_BIOCHEMICAL_CLASSIFICATION.csv')}
    raw = hierarchy(ids, rho)
    raw_group = {i: f'RAW_{k:02d}' for k, g in enumerate(raw['groups'], 1) for i in g}
    by_domain = {d: [i for i in ids if domains[i] == d] for d in COLORS}
    nmr = by_domain['NMR']
    display_groups = []
    for kind in sorted({bio[i]['quantity_kind'] for i in nmr}):
        subset = [i for i in nmr if bio[i]['quantity_kind'] == kind]
        display_groups.extend(((g, 'NMR', kind) for g in hierarchy(subset, rho)['groups']))
    for i in by_domain['Clinical / socio-behavioural']:
        display_groups.append(((i,), 'Clinical / socio-behavioural', 'DISTINCT_CLINICAL_CONSTRUCT'))
    for i in by_domain['Olink']:
        display_groups.append(((i,), 'Olink', 'PROTEINS_RETAINED_SEPARATELY'))
    display_groups.sort(key=lambda x: (list(COLORS).index(x[1]), min((ids.index(i) for i in x[0]))))
    assignments = []
    representatives = []
    for k, (members, domain, construct) in enumerate(display_groups, 1):
        group = f'CORR_{k:02d}'
        central = {i: sum((abs(rho[i][j]) for j in members if j != i)) / (len(members) - 1) if len(members) > 1 else 1.0 for i in members}
        rep = min(members, key=lambda i: (-central[i], -nn[i][i], i))
        minrho = min([abs(rho[i][j]) for i, j in itertools.combinations(members, 2)] + [1.0])
        negative = any((rho[i][j] < 0 for i, j in itertools.combinations(members, 2)))
        note = 'Display candidate only; no result pooling or deletion.'
        if negative:
            note += ' Contains inverse-correlated phenotypes; retain direction/definitions, not interchangeable effects.'
        if domain == 'Olink':
            note += ' Keep each protein irrespective of moderate correlations.'
        if rep == 'smoking':
            note += ' Nominal codes; rank correlation is code-order dependent and not redundancy evidence.'
        for i in members:
            assignments.append({'indicator_id': i, 'indicator_label': labels[i], 'domain': domain, 'cluster_id': group, 'cluster_size': len(members), 'representative_candidate': rep, 'is_representative_candidate': i == rep, 'construct_constraint': construct, 'mean_abs_rho_within_cluster': central[i], 'minimum_abs_rho_within_cluster': minrho, 'rho_to_representative': rho[i][rep], 'orientation_to_representative': 'same' if rho[i][rep] >= 0 else 'inverse', 'raw_hierarchical_cluster_id': raw_group[i], 'n_finite_applicable': nn[i][i], 'note': note})
        representatives.append({'cluster_id': group, 'domain': domain, 'construct_constraint': construct, 'cluster_size': len(members), 'representative_candidate': rep, 'representative_label': labels[rep], 'members': ';'.join(members), 'member_labels': '; '.join((labels[i] for i in members)), 'mean_abs_rho_of_candidate': central[rep], 'candidate_n_finite_applicable': nn[rep][rep], 'minimum_abs_rho_within_cluster': minrho, 'contains_inverse_correlations': negative, 'selection_reason': 'Same broad domain/measurement construct; maximum mean absolute rho (correlation medoid), then finite applicable coverage, then stable ID; no IE/PM/P/q ranking', 'note': note})
    assert len(assignments) == 44 and sum((r['domain'] == 'NMR' for r in representatives)) == 13
    out.mkdir(parents=True)
    (out / 'source').mkdir()
    (out / 'evidence').mkdir()
    (out / 'code').mkdir()
    for name in required:
        shutil.copy2(src / name, out / 'source' / name)
    for file in [corrpath, src / 'evidence/NMR_BIOCHEMICAL_CLASSIFICATION.csv', src / 'evidence/NMR_CORRELATION_SUMMARY.json']:
        shutil.copy2(file, out / 'evidence' / file.name)
    for file in [Path(__file__), HERE / 'extract44.py']:
        shutil.copy2(file, out / 'code' / file.name)
    write(out / 'SELECTION_AUDIT_121.csv', audit)
    write(out / 'SELECTED_44_CURRENT.csv', [r for r in audit if r['selected_44']])
    matrixwrite(out / 'CORR_MATRIX_44.csv', ids, rho)
    matrixwrite(out / 'NPAIR_MATRIX_44.csv', ids, nn)
    matrixwrite(out / 'CORR_MATRIX_44_CLUSTERED.csv', raw['order'], rho)
    matrixwrite(out / 'NPAIR_MATRIX_44_CLUSTERED.csv', raw['order'], nn)
    write(out / 'CLUSTER_ASSIGNMENT_44.csv', assignments)
    write(out / 'REPRESENTATIVE_CANDIDATE_TABLE.csv', representatives)
    write(out / 'LINKAGE_44.csv', [{'node_id': 44 + k, 'left_node': a, 'right_node': b, 'distance_1_minus_abs_rho': d, 'n_members': count} for k, (a, b, d, count) in enumerate(raw['merges'])])
    write(out / 'LEAF_ORDER_44.csv', [{'leaf_position': k + 1, 'indicator_id': i, 'indicator_label': labels[i], 'domain': domains[i]} for k, i in enumerate(raw['order'])])
    pairs_out = []
    for i, j in itertools.combinations(ids, 2):
        pairs_out.append({'indicator_i': i, 'label_i': labels[i], 'domain_i': domains[i], 'indicator_j': j, 'label_j': labels[j], 'domain_j': domains[j], 'rho': rho[i][j], 'abs_rho': abs(rho[i][j]), 'n_pair': nn[i][j], 'calculation_source': evidence['origin'][i][j], 'contains_nominal_smoking': i == 'smoking' or j == 'smoking'})
    write(out / 'CORR_PAIRS_ALL_946.csv', pairs_out)
    high = [r for r in pairs_out if r['abs_rho'] >= 0.9]
    write(out / 'HIGH_CORRELATION_PAIRS.csv', high)
    write(out / 'RAW_CLUSTER_CUT_090.csv', [{'indicator_id': i, 'raw_cluster_id': raw_group[i], 'domain': domains[i], 'indicator_label': labels[i]} for i in raw['order']])
    write(out / 'VARIABLE_VALUES_AND_CODING.csv', [{'indicator_id': i, 'label': labels[i], 'domain': domains[i], **mdata[i]} for i in ids])
    for domain, suffix in [('Clinical / socio-behavioural', 'CLINICAL_SOCIAL'), ('NMR', 'NMR'), ('Olink', 'PROTEIN')]:
        matrixwrite(out / f'CORR_MATRIX_{suffix}.csv', by_domain[domain], rho)
        matrixwrite(out / f'NPAIR_MATRIX_{suffix}.csv', by_domain[domain], nn)
    fontManager.addfont('/System/Library/Fonts/Supplemental/Arial.ttf')
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'text.color': 'black', 'axes.labelcolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    qa = []
    for stem, subset, cl, annotate in [('CORR_HEATMAP_ALL_44', ids, False, False), ('CORR_HEATMAP_ALL_44_CLUSTERED', ids, True, False), ('CORR_HEATMAP_CLINICAL_SOCIAL', by_domain['Clinical / socio-behavioural'], True, False), ('CORR_HEATMAP_NMR', nmr, True, False), ('CORR_HEATMAP_PROTEIN', by_domain['Olink'], True, True)]:
        qa.append(heatmap(out, stem, subset, rho, labels, domains, cl, annotate))
        print('RENDERED ' + stem, flush=True)
    ps = [r for r in pairs_out if r['domain_i'] == r['domain_j'] == 'Olink']
    clin = [r for r in high if r['domain_i'] == r['domain_j'] == 'Clinical / socio-behavioural']
    nmr_groups = [r for r in representatives if r['domain'] == 'NMR']
    fullraw_nmr = hierarchy(nmr, rho)
    report = f"# Indicator correlations and representative candidates\n\nSpearman correlations use pairwise complete observations from the fixed analysis inputs. Pair-specific counts are reported separately. Hierarchical clustering uses 1 - |rho| and complete linkage, with 0.90 as the high-correlation reference. Representatives are chosen within domains and measurement constructs by correlation centrality, data completeness and stable identifier; effects and P/q values are not selection criteria. Smoking codes are descriptive ranks, not an ordinal model assumption.\n\nParent risk-set N: {evidence['N_parent']:,}. Numeric NMR clusters: {len(fullraw_nmr['groups'])}.\n\n|Group|Representative|Members|Minimum absolute correlation|Inverse correlations|\n|---|---|---:|---:|---|\n"
    for r in nmr_groups:
        report += f"\nClinical high-correlation pairs: {len(clin)}. Protein pairwise rho range: {min(r['rho'] for r in ps):.4f} to {max(r['rho'] for r in ps):.4f}. All indicator records remain in the supplementary outputs; representative groups are for display, not pooled statistical effects.\n"
    report += f"\nClinical high-correlation pairs: {len(clin)}. Protein pairwise rho range: {min(r['rho'] for r in ps):.4f} to {max(r['rho'] for r in ps):.4f}. All indicator records remain in the supplementary outputs; representative groups are for display, not pooled statistical effects.\n"
    text(out / 'FIGURE_CORRELATION_REPORT_ZH.md', report)
    text(out / 'REPRODUCE_COMMANDS.md', f'''# 已实际运行命令\n\n仅从本次小型封存相关证据重画（不需读取参与者数据）：\n\n```bash\n"{sys.executable}" "{Path(__file__).resolve()}" --correlation "{corrpath}"\n```\n\n默认自动选择下一个空版本目录，不覆盖旧文件。仅如需重新核查个体相关性：\n\n```bash\n"{sys.executable}" "{HERE / 'extract44.py'}" --out "{HERE / 'CORRELATION44_READONLY_recheck.json'}"\n```\n\n后一个命令只有在目标JSON不存在时才运行，真实参与者始终留Win/G；不运行模型，不增加Bootstrap。\n''')
    assert all((sha(Path(p)) == v for p, v in before.items()))
    for name in required:
        assert sha(out / 'source' / name) == sha(src / name)
    js(out / 'FINAL_QA.json', {'passed': True, 'registered': 121, 'three_condition': 64, 'PM_gt5_selected': 44, 'domains': dict(collections.Counter(domains.values())), 'unique_pairs': 946, 'matrix_shape': [44, 44], 'symmetric_and_valid': True, 'n_pair_bound_checks_passed': True, 'nmr_numeric_clusters': len(fullraw_nmr['groups']), 'nmr_construct_reviewed_groups': 13, 'representative_candidates': len(representatives), 'global_raw_clusters': len(raw['groups']), 'clinical_high_correlation_pairs': len(clin), 'figures': qa, 'source_hashes': before, 'models_fit': 0, 'Bootstrap_added': 0, 'P_q_recalculated': False, 'sankey_or_PM_figures_generated': 0, 'participant_rows_exported': 0, 'all_source_files_unchanged': True})
    js(out / 'FINAL_MANIFEST.json', {str(p.relative_to(out)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(out.rglob('*')) if p.is_file()})
    with zipfile.ZipFile(out.with_suffix('.zip'), 'x', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file():
                z.write(p, p.relative_to(out))
    print(json.dumps({'out': str(out), 'selected': 44, 'NMR_groups': 13, 'total_candidate_groups': len(representatives), 'models_fit': 0, 'sankey_generated': 0}, ensure_ascii=False), flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Selected 44 real-value correlation arrays only. No mediation/P/q calculation."""\nfrom __future__ import annotations\nimport argparse\nimport collections\nimport csv\nfrom datetime import datetime\nimport hashlib\nimport itertools\nimport json\nimport math\nfrom pathlib import Path\nimport shutil\nimport sys\nimport zipfile\nfrom zoneinfo import ZoneInfo\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib.colors import ListedColormap, BoundaryNorm\nfrom matplotlib.font_manager import fontManager\nfrom matplotlib.patches import Rectangle, Patch\nimport numpy as np\nBASE = next((p for p in Path(__file__).resolve().parents if (p / \'scripts/birth\').is_dir() and (p / \'reports/birth\').is_dir()))\nHERE = Path(__file__).resolve().parent\nDOM = {\'危险因素\': \'Clinical / socio-behavioural\', \'社会、行为及生命历程\': \'Clinical / socio-behavioural\', \'NMR\': \'NMR\', \'蛋白组\': \'Olink\'}\nCOLORS = {\'Clinical / socio-behavioural\': \'#7FA3B6\', \'NMR\': \'#8CAB9F\', \'Olink\': \'#A899B7\'}\nCMAP = ListedColormap([\'#1B6B78\', \'#388894\', \'#68A5AF\', \'#A0C5CB\', \'#D1E4E7\', \'#FFFFFF\', \'#F3D5DC\', \'#EAA8B7\', \'#D97C91\', \'#B75871\', \'#8E3652\'])\nNORM = BoundaryNorm(np.linspace(-1, 1, 12), CMAP.N, clip=True)\n\ndef read(p):\n    with p.open(encoding=\'utf-8-sig\', newline=\'\') as f:\n        return list(csv.DictReader(f))\n\ndef sha(p):\n    return hashlib.sha256(p.read_bytes()).hexdigest()\n\ndef write(p, rows):\n    assert rows\n    with p.open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as f:\n        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys((k for r in rows for k in r))))\n        w.writeheader()\n        w.writerows(rows)\n\ndef js(p, x):\n    with p.open(\'x\', encoding=\'utf-8\') as f:\n        json.dump(x, f, ensure_ascii=False, indent=2, allow_nan=False)\n\ndef text(p, s):\n    with p.open(\'x\', encoding=\'utf-8\') as f:\n        f.write(s)\n\ndef num(x):\n    try:\n        v = float(x)\n        return v if math.isfinite(v) else None\n    except (ValueError, TypeError):\n        return None\n\ndef matrixwrite(path, ids, matrix):\n    write(path, [{\'indicator_id\': i, **{j: matrix[i][j] for j in ids}} for i in ids])\n\ndef hierarchy(ids, rho):\n    """Exact deterministic agglomerative complete linkage; no scipy required."""\n    ids = list(ids)\n    n = len(ids)\n    members = {k: (i,) for k, i in enumerate(ids)}\n    active = list(range(n))\n    merges = []\n    for step in range(n - 1):\n        candidates = []\n        for a, b in itertools.combinations(active, 2):\n            distance = max((1 - abs(rho[i][j]) for i in members[a] for j in members[b]))\n            candidates.append((distance, tuple(sorted(members[a] + members[b])), a, b))\n        dist, _, a, b = min(candidates)\n        if members[b] < members[a]:\n            a, b = (b, a)\n        node = n + step\n        members[node] = tuple(sorted(members[a] + members[b]))\n        merges.append((a, b, max(0.0, dist), len(members[node])))\n        active = [v for v in active if v not in (a, b)] + [node]\n\n    def leaves(node):\n        if node < n:\n            return [ids[node]]\n        a, b, _, _ = merges[node - n]\n        return leaves(a) + leaves(b)\n\n    def cut(node):\n        if node < n:\n            return [(ids[node],)]\n        a, b, h, _ = merges[node - n]\n        if h <= 0.1 + 1e-14:\n            return [members[node]]\n        return cut(a) + cut(b)\n    root = active[0]\n    return {\'ids\': ids, \'merges\': merges, \'order\': leaves(root), \'groups\': cut(root), \'members\': members}\n\ndef tree(ax, h, orientation):\n    n = len(h[\'ids\'])\n    position = {i: k for k, i in enumerate(h[\'order\'])}\n    xs = {k: position[i] for k, i in enumerate(h[\'ids\'])}\n    ys = {k: 0.0 for k in range(n)}\n    for k, (a, b, height, _) in enumerate(h[\'merges\']):\n        xx = [xs[a], xs[a], xs[b], xs[b]]\n        yy = [ys[a], height, height, ys[b]]\n        ax.plot(xx if orientation == \'top\' else yy, yy if orientation == \'top\' else xx, color=\'#414141\', lw=0.65)\n        xs[n + k] = (xs[a] + xs[b]) / 2\n        ys[n + k] = height\n    limit = max([x[2] for x in h[\'merges\']] + [0.15]) * 1.04\n    if orientation == \'top\':\n        ax.set(xlim=(-0.5, n - 0.5), ylim=(0, limit))\n        ax.axhline(0.1, color=\'#B75871\', lw=0.65, ls=(0, (3, 2)))\n    else:\n        ax.set(xlim=(limit, 0), ylim=(n - 0.5, -0.5))\n        ax.axvline(0.1, color=\'#B75871\', lw=0.65, ls=(0, (3, 2)))\n    ax.axis(\'off\')\n\ndef heatmap(out, stem, ids, rho, labels, domains, clustered=False, annotate=False):\n    h = hierarchy(ids, rho)\n    order = h[\'order\'] if clustered else list(ids)\n    n = len(order)\n    cell = 0.235 if n > 35 else 0.275 if n > 20 else 0.41 if n > 6 else 0.78\n    side = n * cell\n    dendro = 0.9 if clustered else 0.0\n    left_labels = 2.5 if n > 20 else 2.12\n    bottom = 2.46 if n > 20 else 2.16 if n > 6 else 1.23\n    if n <= 6:\n        left_labels = 0.88\n    top = 0.62\n    right_labels = 0 if not clustered else left_labels\n    left = 0.3 + dendro if clustered else left_labels\n    right = 0.6 + right_labels\n    fw = left + side + right\n    fh = bottom + side + dendro + top\n    fig = plt.figure(figsize=(fw, fh), dpi=100)\n    ax = fig.add_axes([left / fw, bottom / fh, side / fw, side / fh])\n    mat = np.array([[rho[i][j] for j in order] for i in order])\n    im = ax.imshow(mat, cmap=CMAP, norm=NORM, interpolation=\'nearest\', origin=\'upper\')\n    ticks = np.arange(n)\n    names = [labels[i] for i in order]\n    ax.set_xticks(ticks, labels=names, rotation=90, fontsize=8.0 if n > 20 else 9)\n    ax.set_yticks(ticks, labels=names, fontsize=8.0 if n > 20 else 9)\n    ax.tick_params(length=0, pad=4)\n    if clustered:\n        ax.yaxis.tick_right()\n    for spine in ax.spines.values():\n        spine.set(color=\'#404040\', linewidth=0.75)\n    bar = fig.add_axes([left / fw, (bottom + side + 0.025) / fh, side / fw, 0.065 / fh])\n    bar.imshow(np.array([[list(COLORS).index(domains[i]) for i in order]]), cmap=ListedColormap(list(COLORS.values())), vmin=0, vmax=2, aspect=\'auto\', interpolation=\'nearest\')\n    bar.axis(\'off\')\n    if clustered:\n        ta = fig.add_axes([left / fw, (bottom + side + 0.13) / fh, side / fw, dendro / fh])\n        tree(ta, h, \'top\')\n        la = fig.add_axes([0.2 / fw, bottom / fh, dendro / fw, side / fh])\n        tree(la, h, \'left\')\n        for group in h[\'groups\']:\n            if len(group) > 1:\n                positions = [order.index(i) for i in group]\n                start = min(positions)\n                assert sorted(positions) == list(range(start, start + len(group)))\n                ax.add_patch(Rectangle((start - 0.5, start - 0.5), len(group), len(group), fc=\'none\', ec=\'#444444\', lw=0.75))\n    else:\n        for k in range(1, n):\n            if domains[order[k]] != domains[order[k - 1]]:\n                ax.axhline(k - 0.5, color=\'#666666\', lw=0.8)\n                ax.axvline(k - 0.5, color=\'#666666\', lw=0.8)\n    if annotate:\n        for y in range(n):\n            for x in range(n):\n                ax.text(x, y, f\'{mat[y, x]:.2f}\', ha=\'center\', va=\'center\', fontsize=9, color=\'white\' if abs(mat[y, x]) > 0.65 else \'black\')\n    cbax = fig.add_axes([left / fw, 0.37 / fh, min(side, 4.5) / fw, 0.13 / fh])\n    cb = fig.colorbar(im, cax=cbax, orientation=\'horizontal\', ticks=[-1, -0.5, 0, 0.5, 1])\n    cb.ax.tick_params(labelsize=8, length=2)\n    cb.outline.set_linewidth(0.5)\n    cb.set_label(\'Spearman ρ\', fontsize=9, labelpad=3)\n    title = {\'CORR_HEATMAP_ALL_44\': \'Selected mediators (n = 44)\', \'CORR_HEATMAP_ALL_44_CLUSTERED\': \'Selected mediators (n = 44)\', \'CORR_HEATMAP_CLINICAL_SOCIAL\': \'Clinical and behavioural mediators (n = 11)\', \'CORR_HEATMAP_NMR\': \'NMR metabolomic mediators (n = 29)\', \'CORR_HEATMAP_PROTEIN\': \'Olink proteomic mediators (n = 4)\'}[stem]\n    fig.text(left / fw, (fh - 0.24) / fh, title, ha=\'left\', fontsize=11, color=\'black\')\n    if n == 44:\n        fig.legend(handles=[Patch(facecolor=c, label=k) for k, c in COLORS.items()], loc=\'upper right\', bbox_to_anchor=(0.995, 0.998), ncol=3, fontsize=7.2, frameon=False, handlelength=0.8, columnspacing=0.9)\n    fig.canvas.draw()\n    ren = fig.canvas.get_renderer()\n    for a in fig.axes:\n        if not a.axison:\n            continue\n        for q in [*a.get_xticklabels(), *a.get_yticklabels()]:\n            if not q.get_visible() or not q.get_text():\n                continue\n            box = q.get_window_extent(ren)\n            fb = fig.bbox\n            assert fb.contains(box.x0, box.y0) and fb.contains(box.x1, box.y1), (\'TEXT_CLIPPED\', stem, q.get_text())\n    for ext in [\'pdf\', \'svg\', \'png\']:\n        fig.savefig(out / (stem + \'.\' + ext), dpi=320, facecolor=\'white\')\n    fig.savefig(out / (stem + \'_preview.png\'), dpi=110, facecolor=\'white\')\n    plt.close(fig)\n    return {\'stem\': stem, \'n\': n, \'clustered\': clustered, \'order\': order, \'size_inches\': [fw, fh], \'text_bounds_passed\': True}\n\ndef main(args):\n    src = args.master\n    corrpath = args.correlation\n    out = args.out\n    if out is None:\n        v = 1\n        while (BASE / f\'reports/birth/b100_correlation44_20261030_v{v}\').exists():\n            v += 1\n        out = BASE / f\'reports/birth/b100_correlation44_20261030_v{v}\'\n    assert not out.exists(), \'ADDITIVE_OUTPUT_ONLY\'\n    required = [\'B100_121_IHD_T2D_MASTER.csv\', \'B100_121_ALL_DISEASES_MASTER.csv\', \'B100_121_COMPLETE_RESULTS.csv\', \'B100_SIX_CONTRASTS_RESULTS.csv\', \'PRIORITY_NAMED_CURRENT_RECORDS.csv\', \'REVIEW_REPORT_ZH.md\']\n    input_files = [src / f for f in required] + [src / \'FINAL_MANIFEST.json\', src / \'SUMMARY_AND_PROVENANCE.json\', corrpath, src / \'evidence/NMR_BIOCHEMICAL_CLASSIFICATION.csv\', BASE / \'scripts/birth/border_fullfollowup_b100_20261027/input/APPROVED_121_NO_CMR.csv\', BASE / \'reports/birth/figures_b100_pm_4domains_20261030_v2/LABEL_DICTIONARY.csv\']\n    for p in input_files:\n        assert p.is_file(), str(p)\n    before = {str(p): sha(p) for p in input_files}\n    manifest = json.loads((src / \'FINAL_MANIFEST.json\').read_text())\n    for name in required:\n        assert sha(src / name) == manifest[name][\'sha256\']\n    masters = read(src / required[0])\n    wide = read(src / required[1])\n    single = read(src / required[2])\n    pairs = read(src / required[3])\n    priority = read(src / required[4])\n    report_text = (src / required[5]).read_text()\n    assert len(masters) == len(wide) == 121 and len(single) == 605 and (len(pairs) == 726) and (len(priority) == 12) and (len(report_text) > 1000)\n    registry = read(input_files[-2])\n    assert {r[\'indicator_id\'] for r in registry} == {r[\'indicator_id\'] for r in masters}\n    s = {(r[\'indicator_id\'], r[\'endpoint_id\']): r for r in single}\n    p = {r[\'indicator_id\']: r for r in pairs if r[\'endpoint_A\'] == \'411\' and r[\'endpoint_B\'] == \'250.2\'}\n    assert len(s) == 605 and len(p) == 121\n    audit = []\n    selected = []\n    for m in masters:\n        i = m[\'indicator_id\']\n        a = s[i, \'250.2\']\n        b = s[i, \'411\']\n        q = [num(a[\'q_BH_global605_exploratory\']), num(b[\'q_BH_global605_exploratory\'])]\n        delta = num(p[i][\'Delta_IE\'])\n        dq = num(p[i][\'q_BH_within_contrast121_exploratory\'])\n        pm = [num(a[\'PM_percent\']), num(b[\'PM_percent\'])]\n        c1 = any((v is not None and v < 0.05 for v in q))\n        c2 = delta is not None and delta > 0\n        c3 = dq is not None and dq < 0.05\n        c4 = any((v is not None and v > 5 for v in pm))\n        first = c1 and c2 and c3\n        final = first and c4\n        assert first == (m[\'all_three_conditions\'] == \'TRUE\')\n        audit.append({**m, \'condition1_IE_q605_lt005\': c1, \'condition2_Delta_IE_positive\': c2, \'condition3_Delta_IE_q121_lt005\': c3, \'condition4_T2D_or_IHD_PM_gt5_not_absolute\': c4, \'selected_44\': final})\n        if final:\n            selected.append(m)\n    assert sum((r[\'all_three_conditions\'] == \'TRUE\' for r in masters)) == 64 and len(selected) == 44\n    assert collections.Counter((r[\'domain\'] for r in selected)) == {\'危险因素\': 11, \'NMR\': 29, \'蛋白组\': 4}\n    ids = [r[\'indicator_id\'] for r in selected]\n    m = {r[\'indicator_id\']: r for r in masters}\n    domains = {i: DOM[m[i][\'domain\']] for i in ids}\n    dictionary = {r[\'indicator_id\']: r[\'label_en\'] for r in read(input_files[-1])}\n    labels = {i: dictionary.get(i, s[i, \'411\'][\'label_en\']) for i in ids}\n    labels.update({\'smoking\': \'Smoking status\', \'clinical__hdl\': \'Clinical HDL-C\', \'clinical__glucose\': \'Glucose\', \'clinical__hba1c\': \'HbA1c\', \'clinical__tg\': \'Triglycerides\'})\n    evidence = json.loads(corrpath.read_text())\n    assert evidence[\'ids\'] == ids and evidence[\'models_fit\'] == evidence[\'bootstrap_added\'] == evidence[\'participant_rows_exported\'] == 0\n    for path, digest in evidence[\'provenance_sha256\'].items():\n        assert sha(Path(path)) == digest\n    rho = evidence[\'rho\']\n    nn = evidence[\'n_pair\']\n    mdata = evidence[\'metadata\']\n    mat = np.array([[rho[i][j] for j in ids] for i in ids])\n    nmat = np.array([[nn[i][j] for j in ids] for i in ids])\n    assert np.all(np.isfinite(mat)) and np.all(np.abs(mat) <= 1 + 1e-14) and np.allclose(mat, mat.T, rtol=0, atol=1e-14)\n    assert np.allclose(np.diag(mat), 1, rtol=0, atol=1e-14) and np.array_equal(nmat, nmat.T) and np.all(nmat > 2)\n    assert all((nn[i][j] <= min(nn[i][i], nn[j][j]) for i in ids for j in ids))\n    assert all((nn[i][i] == mdata[i][\'n_finite_applicable\'] for i in ids))\n    bio = {r[\'indicator_id\']: r for r in read(src / \'evidence/NMR_BIOCHEMICAL_CLASSIFICATION.csv\')}\n    raw = hierarchy(ids, rho)\n    raw_group = {i: f\'RAW_{k:02d}\' for k, g in enumerate(raw[\'groups\'], 1) for i in g}\n    by_domain = {d: [i for i in ids if domains[i] == d] for d in COLORS}\n    nmr = by_domain[\'NMR\']\n    display_groups = []\n    for kind in sorted({bio[i][\'quantity_kind\'] for i in nmr}):\n        subset = [i for i in nmr if bio[i][\'quantity_kind\'] == kind]\n        display_groups.extend(((g, \'NMR\', kind) for g in hierarchy(subset, rho)[\'groups\']))\n    for i in by_domain[\'Clinical / socio-behavioural\']:\n        display_groups.append(((i,), \'Clinical / socio-behavioural\', \'DISTINCT_CLINICAL_CONSTRUCT\'))\n    for i in by_domain[\'Olink\']:\n        display_groups.append(((i,), \'Olink\', \'PROTEINS_RETAINED_SEPARATELY\'))\n    display_groups.sort(key=lambda x: (list(COLORS).index(x[1]), min((ids.index(i) for i in x[0]))))\n    assignments = []\n    representatives = []\n    for k, (members, domain, construct) in enumerate(display_groups, 1):\n        group = f\'CORR_{k:02d}\'\n        central = {i: sum((abs(rho[i][j]) for j in members if j != i)) / (len(members) - 1) if len(members) > 1 else 1.0 for i in members}\n        rep = min(members, key=lambda i: (-central[i], -nn[i][i], i))\n        minrho = min([abs(rho[i][j]) for i, j in itertools.combinations(members, 2)] + [1.0])\n        negative = any((rho[i][j] < 0 for i, j in itertools.combinations(members, 2)))\n        note = \'Display candidate only; no result pooling or deletion.\'\n        if negative:\n            note += \' Contains inverse-correlated phenotypes; retain direction/definitions, not interchangeable effects.\'\n        if domain == \'Olink\':\n            note += \' Keep each protein irrespective of moderate correlations.\'\n        if rep == \'smoking\':\n            note += \' Nominal codes; rank correlation is code-order dependent and not redundancy evidence.\'\n        for i in members:\n            assignments.append({\'indicator_id\': i, \'indicator_label\': labels[i], \'domain\': domain, \'cluster_id\': group, \'cluster_size\': len(members), \'representative_candidate\': rep, \'is_representative_candidate\': i == rep, \'construct_constraint\': construct, \'mean_abs_rho_within_cluster\': central[i], \'minimum_abs_rho_within_cluster\': minrho, \'rho_to_representative\': rho[i][rep], \'orientation_to_representative\': \'same\' if rho[i][rep] >= 0 else \'inverse\', \'raw_hierarchical_cluster_id\': raw_group[i], \'n_finite_applicable\': nn[i][i], \'note\': note})\n        representatives.append({\'cluster_id\': group, \'domain\': domain, \'construct_constraint\': construct, \'cluster_size\': len(members), \'representative_candidate\': rep, \'representative_label\': labels[rep], \'members\': \';\'.join(members), \'member_labels\': \'; \'.join((labels[i] for i in members)), \'mean_abs_rho_of_candidate\': central[rep], \'candidate_n_finite_applicable\': nn[rep][rep], \'minimum_abs_rho_within_cluster\': minrho, \'contains_inverse_correlations\': negative, \'selection_reason\': \'Same broad domain/measurement construct; maximum mean absolute rho (correlation medoid), then finite applicable coverage, then stable ID; no IE/PM/P/q ranking\', \'note\': note})\n    assert len(assignments) == 44 and sum((r[\'domain\'] == \'NMR\' for r in representatives)) == 13\n    out.mkdir(parents=True)\n    (out / \'source\').mkdir()\n    (out / \'evidence\').mkdir()\n    (out / \'code\').mkdir()\n    for name in required:\n        shutil.copy2(src / name, out / \'source\' / name)\n    for file in [corrpath, src / \'evidence/NMR_BIOCHEMICAL_CLASSIFICATION.csv\', src / \'evidence/NMR_CORRELATION_SUMMARY.json\']:\n        shutil.copy2(file, out / \'evidence\' / file.name)\n    for file in [Path(__file__), HERE / \'extract44.py\']:\n        shutil.copy2(file, out / \'code\' / file.name)\n    write(out / \'SELECTION_AUDIT_121.csv\', audit)\n    write(out / \'SELECTED_44_CURRENT.csv\', [r for r in audit if r[\'selected_44\']])\n    matrixwrite(out / \'CORR_MATRIX_44.csv\', ids, rho)\n    matrixwrite(out / \'NPAIR_MATRIX_44.csv\', ids, nn)\n    matrixwrite(out / \'CORR_MATRIX_44_CLUSTERED.csv\', raw[\'order\'], rho)\n    matrixwrite(out / \'NPAIR_MATRIX_44_CLUSTERED.csv\', raw[\'order\'], nn)\n    write(out / \'CLUSTER_ASSIGNMENT_44.csv\', assignments)\n    write(out / \'REPRESENTATIVE_CANDIDATE_TABLE.csv\', representatives)\n    write(out / \'LINKAGE_44.csv\', [{\'node_id\': 44 + k, \'left_node\': a, \'right_node\': b, \'distance_1_minus_abs_rho\': d, \'n_members\': count} for k, (a, b, d, count) in enumerate(raw[\'merges\'])])\n    write(out / \'LEAF_ORDER_44.csv\', [{\'leaf_position\': k + 1, \'indicator_id\': i, \'indicator_label\': labels[i], \'domain\': domains[i]} for k, i in enumerate(raw[\'order\'])])\n    pairs_out = []\n    for i, j in itertools.combinations(ids, 2):\n        pairs_out.append({\'indicator_i\': i, \'label_i\': labels[i], \'domain_i\': domains[i], \'indicator_j\': j, \'label_j\': labels[j], \'domain_j\': domains[j], \'rho\': rho[i][j], \'abs_rho\': abs(rho[i][j]), \'n_pair\': nn[i][j], \'calculation_source\': evidence[\'origin\'][i][j], \'contains_nominal_smoking\': i == \'smoking\' or j == \'smoking\'})\n    write(out / \'CORR_PAIRS_ALL_946.csv\', pairs_out)\n    high = [r for r in pairs_out if r[\'abs_rho\'] >= 0.9]\n    write(out / \'HIGH_CORRELATION_PAIRS.csv\', high)\n    write(out / \'RAW_CLUSTER_CUT_090.csv\', [{\'indicator_id\': i, \'raw_cluster_id\': raw_group[i], \'domain\': domains[i], \'indicator_label\': labels[i]} for i in raw[\'order\']])\n    write(out / \'VARIABLE_VALUES_AND_CODING.csv\', [{\'indicator_id\': i, \'label\': labels[i], \'domain\': domains[i], **mdata[i]} for i in ids])\n    for domain, suffix in [(\'Clinical / socio-behavioural\', \'CLINICAL_SOCIAL\'), (\'NMR\', \'NMR\'), (\'Olink\', \'PROTEIN\')]:\n        matrixwrite(out / f\'CORR_MATRIX_{suffix}.csv\', by_domain[domain], rho)\n        matrixwrite(out / f\'NPAIR_MATRIX_{suffix}.csv\', by_domain[domain], nn)\n    fontManager.addfont(\'/System/Library/Fonts/Supplemental/Arial.ttf\')\n    plt.rcParams.update({\'font.family\': \'Arial\', \'font.size\': 9, \'text.color\': \'black\', \'axes.labelcolor\': \'black\', \'xtick.color\': \'black\', \'ytick.color\': \'black\', \'pdf.fonttype\': 42, \'svg.fonttype\': \'none\'})\n    qa = []\n    for stem, subset, cl, annotate in [(\'CORR_HEATMAP_ALL_44\', ids, False, False), (\'CORR_HEATMAP_ALL_44_CLUSTERED\', ids, True, False), (\'CORR_HEATMAP_CLINICAL_SOCIAL\', by_domain[\'Clinical / socio-behavioural\'], True, False), (\'CORR_HEATMAP_NMR\', nmr, True, False), (\'CORR_HEATMAP_PROTEIN\', by_domain[\'Olink\'], True, True)]:\n        qa.append(heatmap(out, stem, subset, rho, labels, domains, cl, annotate))\n        print(\'RENDERED \' + stem, flush=True)\n    ps = [r for r in pairs_out if r[\'domain_i\'] == r[\'domain_j\'] == \'Olink\']\n    clin = [r for r in high if r[\'domain_i\'] == r[\'domain_j\'] == \'Clinical / socio-behavioural\']\n    nmr_groups = [r for r in representatives if r[\'domain\'] == \'NMR\']\n    fullraw_nmr = hierarchy(nmr, rho)\n    report = f"# Indicator correlations and representative candidates\\n\\nSpearman correlations use pairwise complete observations from the fixed analysis inputs. Pair-specific counts are reported separately. Hierarchical clustering uses 1 - |rho| and complete linkage, with 0.90 as the high-correlation reference. Representatives are chosen within domains and measurement constructs by correlation centrality, data completeness and stable identifier; effects and P/q values are not selection criteria. Smoking codes are descriptive ranks, not an ordinal model assumption.\\n\\nParent risk-set N: {evidence[\'N_parent\']:,}. Numeric NMR clusters: {len(fullraw_nmr[\'groups\'])}.\\n\\n|Group|Representative|Members|Minimum absolute correlation|Inverse correlations|\\n|---|---|---:|---:|---|\\n"\n    for r in nmr_groups:\n        report += f"\\nClinical high-correlation pairs: {len(clin)}. Protein pairwise rho range: {min(r[\'rho\'] for r in ps):.4f} to {max(r[\'rho\'] for r in ps):.4f}. All indicator records remain in the supplementary outputs; representative groups are for display, not pooled statistical effects.\\n"\n    report += f"\\nClinical high-correlation pairs: {len(clin)}. Protein pairwise rho range: {min(r[\'rho\'] for r in ps):.4f} to {max(r[\'rho\'] for r in ps):.4f}. All indicator records remain in the supplementary outputs; representative groups are for display, not pooled statistical effects.\\n"\n    text(out / \'FIGURE_CORRELATION_REPORT_ZH.md\', report)\n    text(out / \'REPRODUCE_COMMANDS.md\', f\'\'\'# 已实际运行命令\\n\\n仅从本次小型封存相关证据重画（不需读取参与者数据）：\\n\\n```bash\\n"{sys.executable}" "{Path(__file__).resolve()}" --correlation "{corrpath}"\\n```\\n\\n默认自动选择下一个空版本目录，不覆盖旧文件。仅如需重新核查个体相关性：\\n\\n```bash\\n"{sys.executable}" "{HERE / \'extract44.py\'}" --out "{HERE / \'CORRELATION44_READONLY_recheck.json\'}"\\n```\\n\\n后一个命令只有在目标JSON不存在时才运行，真实参与者始终留Win/G；不运行模型，不增加Bootstrap。\\n\'\'\')\n    assert all((sha(Path(p)) == v for p, v in before.items()))\n    for name in required:\n        assert sha(out / \'source\' / name) == sha(src / name)\n    js(out / \'FINAL_QA.json\', {\'passed\': True, \'registered\': 121, \'three_condition\': 64, \'PM_gt5_selected\': 44, \'domains\': dict(collections.Counter(domains.values())), \'unique_pairs\': 946, \'matrix_shape\': [44, 44], \'symmetric_and_valid\': True, \'n_pair_bound_checks_passed\': True, \'nmr_numeric_clusters\': len(fullraw_nmr[\'groups\']), \'nmr_construct_reviewed_groups\': 13, \'representative_candidates\': len(representatives), \'global_raw_clusters\': len(raw[\'groups\']), \'clinical_high_correlation_pairs\': len(clin), \'figures\': qa, \'source_hashes\': before, \'models_fit\': 0, \'Bootstrap_added\': 0, \'P_q_recalculated\': False, \'sankey_or_PM_figures_generated\': 0, \'participant_rows_exported\': 0, \'all_source_files_unchanged\': True})\n    js(out / \'FINAL_MANIFEST.json\', {str(p.relative_to(out)): {\'sha256\': sha(p), \'bytes\': p.stat().st_size} for p in sorted(out.rglob(\'*\')) if p.is_file()})\n    with zipfile.ZipFile(out.with_suffix(\'.zip\'), \'x\', zipfile.ZIP_DEFLATED) as z:\n        for p in sorted(out.rglob(\'*\')):\n            if p.is_file():\n                z.write(p, p.relative_to(out))\n    print(json.dumps({\'out\': str(out), \'selected\': 44, \'NMR_groups\': 13, \'total_candidate_groups\': len(representatives), \'models_fit\': 0, \'sankey_generated\': 0}, ensure_ascii=False), flush=True)\nif __name__ == \'__main__\':\n    p = argparse.ArgumentParser(description=__doc__)\n    p.add_argument(\'--master\', type=Path, default=BASE / \'reports/birth/b100_master121_audit_20261030_v1\')\n    p.add_argument(\'--correlation\', type=Path, default=HERE / \'CORRELATION44_READONLY.json\')\n    p.add_argument(\'--out\', type=Path)\n    main(p.parse_args())', '<private-stage>', 'exec'), context)
