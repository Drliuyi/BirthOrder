"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import ast
import collections
import csv
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import shutil
import sys
from zoneinfo import ZoneInfo
EPS = {'250': 'DM', '250.2': 'T2D', '411': 'IHD', '411.3': 'Angina', '411.4': 'CAS'}
METHOD = birth_private_setting('PRIVATE_7AA11EFF31702D72C96A44823457AB674B93D8CC8A1AD7965B72932F60DA')
SOURCE = birth_private_setting('PRIVATE_1095A231CC6F1187B50096D4E8F4F675170A4FA973A9BB301E2CDD931030')
VERSIONS = {'B100_ORIGINAL_ESTIMATOR_WARNING_AWARE_RELEASE_V1', 'B100_MISSING_INFERENCE_COMPLETION_20261030_V1_ORIGINAL_ESTIMATOR'}
FOCUS = ['qualification', 'income', 'age_education', 'birth_weight', 'insomnia', 'townsend', 'smoking', 'commute_car', 'commute_public_transport', 'commute_walking', 'packyears', 'maternal_smoking']
UNKNOWN = 'NOT_EVALUABLE_CURRENT_RECORD'
UNVERIFIED = 'CURRENT_VERSION_NOT_VERIFIED'
FLAG_LABELS = {'T2D_IE_q605_lt005': 'T2D单病IE：当前q605<0.05', 'IHD_IE_q605_lt005': 'IHD单病IE：当前q605<0.05', 'at_least_one_IE_q605_lt005': '条件1：两病至少一病IE q605<0.05', 'Delta_IE_positive': '条件2：IHD−T2D的ΔIE>0', 'Delta_IE_q121_lt005': '条件3：配对ΔIE当前q121<0.05', 'all_three_conditions': '上述三条件同时满足（仅展示筛选，不是结果是否存在）'}

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write(path, rows):
    assert rows
    fields = list(dict.fromkeys((k for r in rows for k in r)))
    with path.open('x', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def jwrite(path, obj):
    with path.open('x', encoding='utf-8') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)

def txt(path, value):
    with path.open('x', encoding='utf-8') as f:
        f.write(value)

def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None

def state(value, predicate):
    n = number(value)
    return UNKNOWN if n is None else 'TRUE' if predicate(n) else 'FALSE'

def either(a, b):
    if 'TRUE' in [a, b]:
        return 'TRUE'
    return 'FALSE' if a == b == 'FALSE' else UNKNOWN

def all_flags(flags):
    if any((v not in ['TRUE', 'FALSE'] for v in flags)):
        return UNKNOWN
    return 'TRUE' if all((v == 'TRUE' for v in flags)) else 'FALSE'

def percent_set(s):
    return re.sub('[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?', lambda m: format(float(m.group()) * 100, '.12g'), s)

def close(a, b):
    return abs(float(a) - float(b)) <= 1e-12 + 1e-09 * max(abs(float(a)), abs(float(b)))

def pretty(v, digits=4):
    n = number(v)
    return 'NA' if n is None else f'{n:.{digits}g}'

def ci(row, prefix):
    return f"{pretty(row.get(prefix + '_IE_per_1000'))} [{pretty(row.get(prefix + '_IE_CI_low_per_1000'))}, {pretty(row.get(prefix + '_IE_CI_high_per_1000'))}]"

def run(args):
    source = args.source
    out = args.out
    if out is None:
        v = 1
        while (BASE / f'reports/birth/b100_master121_audit_20261030_v{v}').exists():
            v += 1
        out = BASE / f'reports/birth/b100_master121_audit_20261030_v{v}'
    assert not out.exists(), 'PRESERVE_EXISTING_OUTPUT:' + str(out)
    project = BASE / 'scripts/birth/border_fullfollowup_b100_20261027'
    registry_path = project / 'input/APPROVED_121_NO_CMR.csv'
    matrix_path = project / 'input/TASK_MATRIX_605.csv'
    single_path = source / 'B100_121_COMPLETE_RESULTS.csv'
    pair_path = source / 'B100_SIX_CONTRASTS_RESULTS.csv'
    release_root = source.parent
    table_manifest = source / 'TABLE_MANIFEST.json'
    receipt = release_root / 'PUBLICATION_RECEIPT_v2.json'
    transfer = release_root / 'PUBLICATION_TRANSFER_MANIFEST_v2.json'
    qa_path = source / 'INDEPENDENT_AGGREGATE_QA.json'
    corr_root = args.correlation
    corr_names = ['NMR_CORRELATION_MATRIX.csv', 'NMR_CORRELATION_PAIR_N.csv', 'NMR_CORRELATION_VARIABLE_AUDIT.csv', 'NMR_BIOCHEMICAL_CLASSIFICATION.csv', 'NMR_CORRELATION_SUMMARY.json', 'VALIDATION_RECEIPT.json', 'CORRELATION_EXECUTION_RECEIPT.json']
    cluster_source = BASE / 'scripts/birth/b100_pm_4domain_sankey_20261030_v1/build.py'
    inputs = [registry_path, matrix_path, single_path, pair_path, table_manifest, receipt, transfer, qa_path, cluster_source] + [corr_root / n for n in corr_names]
    for p in inputs:
        assert p.is_file(), 'MISSING_INPUT:' + str(p)
    hashes = {str(p): sha(p) for p in inputs}
    single_digest, pair_digest = (hashes[str(single_path)], hashes[str(pair_path)])
    manifest = json.loads(table_manifest.read_text())
    transferred = json.loads(transfer.read_text())
    publication = json.loads(receipt.read_text())
    for p in [single_path, pair_path]:
        assert sha(p) == manifest[p.name]['sha256'] == transferred[p.name]['sha256']
        assert p.stat().st_size == manifest[p.name]['bytes'] == transferred[p.name]['bytes']
    assert publication['complete']['registered_indicators'] == 121
    assert publication['complete']['delivery_complete']
    assert json.loads(qa_path.read_text())['passed']
    registry = read(registry_path)
    reg = {r['indicator_id']: r for r in registry}
    matrix = read(matrix_path)
    singles = read(single_path)
    pairs = read(pair_path)
    s = {(r['indicator_id'], r['endpoint_id']): r for r in singles}
    p = {(r['indicator_id'], r['contrast']): r for r in pairs}
    assert len(registry) == len(reg) == 121
    expected_single = set(itertools.product(reg, EPS))
    expected_pair = set(itertools.product(reg, CONTRASTS))
    assert len(singles) == len(s) == 605 and set(s) == expected_single
    assert len(pairs) == len(p) == 726 and set(p) == expected_pair
    assert {(r['indicator_id'], r['endpoint_id']) for r in matrix} == expected_single
    assert all((r['source_method_hash'] == METHOD and r['source_hash'] == SOURCE and (r['version'] in VERSIONS) and (r['horizon_years'] == '11') and (r['model_followup'] == 'FULL_FOLLOWUP') and (r['bootstrap_planned'] == '100') for r in singles))
    for r in singles:
        assert r['domain'] == reg[r['indicator_id']]['category']
        if number(r['PM_percent']) is not None:
            assert close(r['PM_percent'], 100 * float(r['IE']) / float(r['TE']))
    for r in pairs:
        assert r['horizon_years'] == '11' and r['bootstrap_planned'] == '100'
        a, b = (s[r['indicator_id'], r['endpoint_A']], s[r['indicator_id'], r['endpoint_B']])
        if number(r['Delta_IE']) is not None:
            assert close(r['IE_A_per_1000'], a['IE_per_1000'])
            assert close(r['IE_B_per_1000'], b['IE_per_1000'])
            assert close(r['Delta_IE'], float(a['IE']) - float(b['IE']))
            assert r['pairing_evidence'] in {'ACTUAL_SHARED_PARENT_DRAW_HASH_VERIFIED_FOR_EACH_VALID_REPEAT', 'EXACT_OLD_ACTUAL_PARENT_DRAW_HASH_PER_ID_VERIFIED;NEW_A_ORIGIN_RETAINED'}
    export_version = 'B100_MASTER121_AUDIT_20261030_V1_CURRENT_COMPLETED_RELEASE'
    master = []
    filters = []
    records = []
    for rank, entry in enumerate(registry, 1):
        i = entry['indicator_id']
        t, h, d = (s[i, '250.2'], s[i, '411'], p[i, '411_minus_250.2'])
        flags = {'T2D_IE_q605_lt005': state(t['q_BH_global605_exploratory'], lambda x: x < 0.05), 'IHD_IE_q605_lt005': state(h['q_BH_global605_exploratory'], lambda x: x < 0.05), 'Delta_IE_positive': state(d['Delta_IE'], lambda x: x > 0), 'Delta_IE_q121_lt005': state(d['q_BH_within_contrast121_exploratory'], lambda x: x < 0.05)}
        flags['at_least_one_IE_q605_lt005'] = either(flags['T2D_IE_q605_lt005'], flags['IHD_IE_q605_lt005'])
        flags['all_three_conditions'] = all_flags([flags['at_least_one_IE_q605_lt005'], flags['Delta_IE_positive'], flags['Delta_IE_q121_lt005']])
        reason = '；'.join((FLAG_LABELS[k] + '=FALSE' for k in ['at_least_one_IE_q605_lt005', 'Delta_IE_positive', 'Delta_IE_q121_lt005'] if flags[k] == 'FALSE'))
        if flags['all_three_conditions'] == UNKNOWN:
            reason = '已核对当前版本，但无可估计数值，不能判成不显著或方向不符；' + d['failure_reason']
        elif not reason:
            reason = '三条件均TRUE；仅符合既定展示筛选，不表示独立确认或PM有界'
        record_status = '已核对当前版本：有估计' if number(d['Delta_IE']) is not None else '已核对当前版本：不可估，原因保留'
        row = {'registry_order': rank, 'indicator_id': i, 'indicator_name_zh': entry['label_zh'], 'domain': entry['category'], 'indicator_type': entry['type'], 'current_record_status': record_status, 'current_version_verified': 'TRUE', 'export_version': export_version, 'method_hash_verified': METHOD, 'source_hash_verified': SOURCE, 'single_source_sha256': single_digest, 'paired_source_sha256': pair_digest, 'T2D_source_row_number': singles.index(t) + 2, 'IHD_source_row_number': singles.index(h) + 2, 'paired_source_row_number': pairs.index(d) + 2, 'three_condition_reason_zh': reason, **flags}
        for prefix, r in [('T2D', t), ('IHD', h), ('paired', d)]:
            row.update({prefix + '_' + k: v for k, v in r.items() if k not in ['indicator_id', 'domain']})
        for prefix, r in [('T2D', t), ('IHD', h)]:
            row[prefix + '_PM_Fieller_set_percent'] = percent_set(r['PM_Fieller_set'])
            row[prefix + '_PM_Fieller_original_set_unit'] = 'ratio; PM_percent is already percent'
        master.append(row)
        raw_values = {'T2D_IE_q605_lt005': t['q_BH_global605_exploratory'], 'IHD_IE_q605_lt005': h['q_BH_global605_exploratory'], 'at_least_one_IE_q605_lt005': json.dumps([t['q_BH_global605_exploratory'], h['q_BH_global605_exploratory']]), 'Delta_IE_positive': d['Delta_IE'], 'Delta_IE_q121_lt005': d['q_BH_within_contrast121_exploratory'], 'all_three_conditions': json.dumps([flags[k] for k in ['at_least_one_IE_q605_lt005', 'Delta_IE_positive', 'Delta_IE_q121_lt005']])}
        for k, label in FLAG_LABELS.items():
            filters.append({'indicator_id': i, 'indicator_name_zh': entry['label_zh'], 'domain': entry['category'], 'condition': k, 'condition_zh': label, 'actual_current_value': raw_values[k], 'condition_state': flags[k], 'current_record_status': record_status, 'source_single_sha256': single_digest, 'source_paired_sha256': pair_digest, 'reason_if_not_evaluable': d['failure_reason'] if flags[k] == UNKNOWN else ''})
        for ep in EPS:
            r = s[i, ep]
            records.append({'indicator_id': i, 'endpoint_id': ep, 'source_row_number': singles.index(r) + 2, 'record_found': 'TRUE', 'version_verified': 'TRUE', 'version': r['version'], 'source_method_hash': r['source_method_hash'], 'source_hash': r['source_hash'], 'source_file_sha256': single_digest, 'result_status': r['result_status'], 'failure_reason': r['failure_reason']})
    by_id = {r['indicator_id']: r for r in master}
    all_disease = []
    for row in master:
        i = row['indicator_id']
        wide = {k: row[k] for k in ['registry_order', 'indicator_id', 'indicator_name_zh', 'domain', 'indicator_type', 'current_record_status', 'export_version', 'method_hash_verified']}
        for ep, name in EPS.items():
            wide.update({name + '_' + k: v for k, v in s[i, ep].items() if k not in ['indicator_id', 'domain']})
        for contrast in CONTRASTS:
            wide.update({contrast + '_' + k: v for k, v in p[i, contrast].items() if k not in ['indicator_id', 'domain']})
        all_disease.append(wide)
    corr_summary = json.loads((corr_root / 'NMR_CORRELATION_SUMMARY.json').read_text())
    validation = json.loads((corr_root / 'VALIDATION_RECEIPT.json').read_text())
    assert corr_summary['registry_sha256'] == sha(registry_path)
    assert corr_summary['schema'] == 'NMR_REAL_PARTICIPANT_SPEARMAN_V1'
    assert corr_summary['nmr_count'] == 54 and validation['status'] == 'PASS'
    assert corr_summary['models_run'] == 0 and corr_summary['participant_rows_exported'] == 0
    rho = {r['indicator_id']: {k: float(v) for k, v in r.items() if k != 'indicator_id'} for r in read(corr_root / 'NMR_CORRELATION_MATRIX.csv')}
    ns = {r['indicator_id']: {k: int(v) for k, v in r.items() if k != 'indicator_id'} for r in read(corr_root / 'NMR_CORRELATION_PAIR_N.csv')}
    bio = {r['indicator_id']: r for r in read(corr_root / 'NMR_BIOCHEMICAL_CLASSIFICATION.csv')}
    cov = {r['indicator_id']: int(r['N_finite']) for r in read(corr_root / 'NMR_CORRELATION_VARIABLE_AUDIT.csv')}
    nmr = [r['indicator_id'] for r in master if r['domain'] == 'NMR']
    assert len(nmr) == 54 and set(nmr) == set(rho) == set(ns) == set(bio) == set(cov)
    for a, b in itertools.product(nmr, repeat=2):
        assert math.isfinite(rho[a][b]) and -1 <= rho[a][b] <= 1
        assert abs(rho[a][b] - rho[b][a]) < 1e-12
        assert ns[a][b] == ns[b][a] and 0 < ns[a][b] <= min(cov[a], cov[b])
    tree = ast.parse(cluster_source.read_text())
    fn = next((n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'cluster_nmr'))
    snippet = ast.get_source_segment(cluster_source.read_text(), fn)
    env = {'collections': collections, 'itertools': itertools, 'math': math}
    exec(compile(snippet, str(cluster_source), 'exec'), env)
    order = {r['indicator_id']: r['registry_order'] for r in master}
    eligible_nmr = [i for i in nmr if by_id[i]['all_three_conditions'] == 'TRUE']
    proposals, merges = ([], [])
    group_counts = {}
    for label, ids in [('ALL_REGISTERED_NMR54', nmr), ('THREE_CONDITION_ELIGIBLE_NMR', eligible_nmr)]:
        clusters, log = env['cluster_nmr'](ids, bio, rho, cov, order)
        group_counts[label] = {'input_markers': len(ids), 'proposed_display_groups': len(clusters), 'proposed_display_reduction': len(ids) - len(clusters), 'non_singleton_groups': sum((len(cl['members']) > 1 for cl in clusters))}
        for cl in clusters:
            for i in cl['members']:
                proposals.append({'assessment_scope': label, 'cluster_id': cl['cluster_id'], 'indicator_id': i, 'indicator_name_zh': reg[i]['label_zh'], 'quantity_kind': cl['quantity_kind'], 'members': ';'.join(cl['members']), 'representative_proposal_only': cl['representative'], 'mean_within_cluster_rho': cl['mean_rho'][i], 'min_within_cluster_rho': cl['minimum_within_rho'], 'marker_finite_N': cov[i], 'proposal_only_not_applied': True, 'IE_PM_never_pooled': True})
        merges.extend(({'assessment_scope': label, **r, 'proposal_only_not_applied': True} for r in log))
    pairwise = [{'indicator_A': a, 'indicator_B': b, 'rho': rho[a][b], 'actual_pair_N': ns[a][b], 'same_quantity_kind': bio[a]['quantity_kind'] == bio[b]['quantity_kind'], 'positive_rho_ge090': rho[a][b] >= 0.9, 'both_three_condition_TRUE': a in eligible_nmr and b in eligible_nmr, 'correlation_source': 'Existing real participant-level descriptive audit, not result correlations'} for a, b in itertools.combinations(nmr, 2)]
    out.mkdir(parents=True)
    (out / 'source').mkdir()
    (out / 'evidence').mkdir()
    (out / 'code').mkdir()
    for pth in [single_path, pair_path]:
        shutil.copy2(pth, out / pth.name)
        assert sha(out / pth.name) == sha(pth)
    for pth in [registry_path, matrix_path, table_manifest, receipt, transfer, qa_path]:
        shutil.copy2(pth, out / 'source' / pth.name)
    for name in corr_names:
        shutil.copy2(corr_root / name, out / 'evidence' / name)
    shutil.copy2(Path(__file__), out / 'code/export.py')
    txt(out / 'code/reused_cluster_nmr.py.txt', snippet)
    write(out / 'B100_121_IHD_T2D_MASTER.csv', master)
    write(out / 'B100_121_ALL_DISEASES_MASTER.csv', all_disease)
    write(out / 'B100_121_FILTER_EXPLANATIONS.csv', filters)
    write(out / 'CURRENT_VERSION_605_RECORD_AUDIT.csv', records)
    priority = [by_id[i] for i in FOCUS]
    write(out / 'PRIORITY_NAMED_CURRENT_RECORDS.csv', priority)
    write(out / 'NMR_DISPLAY_REDUNDANCY_PROPOSAL_ONLY.csv', proposals)
    if merges:
        write(out / 'NMR_CLUSTER_STEPS_PROPOSAL_ONLY.csv', merges)
    write(out / 'NMR54_ACTUAL_CORRELATION_PAIR_N.csv', pairwise)
    counts = {'registered': len(master), 'paired_estimable': sum((number(r['paired_Delta_IE']) is not None for r in master)), 'unverified_current_records': 0, 'single_cells': len(singles), 'paired_rows': len(pairs), 'single_key_duplicates': len(singles) - len(s), 'pair_key_duplicates': len(pairs) - len(p)}
    for key in FLAG_LABELS:
        counts[key] = dict(collections.Counter((r[key] for r in master)))
    domains = []
    for domain in ['危险因素', '社会、行为及生命历程', 'NMR', '蛋白组']:
        rows = [r for r in master if r['domain'] == domain]
        domains.append({'domain': domain, 'registered': len(rows), 'paired_estimable': sum((number(r['paired_Delta_IE']) is not None for r in rows)), 'at_least_one_IE_q605_lt005': sum((r['at_least_one_IE_q605_lt005'] == 'TRUE' for r in rows)), 'IE_supported_and_Delta_positive': sum((r['at_least_one_IE_q605_lt005'] == r['Delta_IE_positive'] == 'TRUE' for r in rows)), 'all_three_TRUE': sum((r['all_three_conditions'] == 'TRUE' for r in rows))})
    write(out / 'COUNTS_BY_DOMAIN.csv', domains)
    summary = {'created_beijing': datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(), 'source_local_release': str(source), 'source_remote_release_receipt': publication['remote'], 'source_method_hash': METHOD, 'source_source_hash': SOURCE, 'source_result_versions': dict(collections.Counter((r['version'] for r in singles))), 'universe_origin': str(registry_path), 'never_reverse_inferred_from_77_or_64': True, 'counts': counts, 'domains': domains, 'NMR_display_redundancy_assessment_only': group_counts, 'models_fit': 0, 'bootstrap_run': 0, 'P_q_recalculated': False, 'PM_figures_generated': 0, 'actual_merges_applied': 0, 'old_result_fallback_used': False, 'source_sha256': hashes, 'old_results_unchanged': True}
    jwrite(out / 'SUMMARY_AND_PROVENANCE.json', summary)
    lines = ['# 当前完整B100：121项主表、筛选和合并评估（未出图）', '', '仅从当前完整605单病表和726配对表读取结果，以原冻结121项稳定ID注册表构建全集。没有读取77/64筛后名单，没有用旧缓存填充结果。', '', f'已核对605/605单病、726/726配对记录：主表121项，其中110项有IHD−T2D估计、11项有当前不可估记录；尚未核对的当前记录0项。重复主键0项，不需要为清理重复结果而合并。', '', '两种行版本标签分别为原有效结果565格及同估计器补全40格，均属于当前完整发布；统计method/source哈希一致。不是把旧的无蛋白区间版本混入。', '', '## 先看多少项', '', '|领域|完整注册|可估配对|至少一病IE q605<0.05|再要求ΔIE>0|三条件全部满足|', '|---|---:|---:|---:|---:|---:|']
    for r in domains:
        lines.append('|' + '|'.join((str(r[k]) for k in ['domain', 'registered', 'paired_estimable', 'at_least_one_IE_q605_lt005', 'IE_supported_and_Delta_positive', 'all_three_TRUE'])) + '|')
    lines += ['', '共96项至少一病IE q605<0.05；其中77项ΔIE>0；加入配对q121<0.05后为64项。这些数字是筛选层级，不是把其余指标从结果全集删除。', '', '## 点名指标当前记录', '', 'IE和ΔIE均为每千人11年风险差；下表q为各自当前固定家族。', '', '|指标（稳定ID）|T2D IE [95%CI]|IHD IE [95%CI]|ΔIE [95%CI]|T2D / IHD q605|ΔIE P / q121|有效B T2D/IHD/配对|条件1/2/3|', '|---|---|---|---|---|---|---|---|']
    for r in priority:
        delta = f"{pretty(r['paired_Delta_IE_per_1000'])} [{pretty(r['paired_Delta_IE_CI_low_per_1000'])}, {pretty(r['paired_Delta_IE_CI_high_per_1000'])}]"
        lines.append(f"|{r['indicator_name_zh']} (`{r['indicator_id']}`)|{ci(r, 'T2D')}|{ci(r, 'IHD')}|{delta}|{pretty(r['T2D_q_BH_global605_exploratory'])} / {pretty(r['IHD_q_BH_global605_exploratory'])}|{pretty(r['paired_P_Delta_IE'])} / {pretty(r['paired_q_BH_within_contrast121_exploratory'])}|{r['T2D_bootstrap_valid']}/{r['IHD_bootstrap_valid']}/{r['paired_paired_valid']}|{r['at_least_one_IE_q605_lt005']} / {r['Delta_IE_positive']} / {r['Delta_IE_q121_lt005']}|")
    lines += ['', '点名记录的完整单病P、两种q、PM、Fieller比例集合及百分数集合、版本和原因均在PRIORITY_NAMED_CURRENT_RECORDS.csv；此处没有用表格缩略值作筛选。', '', '学历、收入、教育结束年龄、出生体重、失眠和Townsend均有当前估计，不应说成“没有结果”。学历/收入/教育结束年龄的ΔIE点值为负；出生体重/失眠/Townsend为正，但本次配对q未达0.05。吸烟与汽车通勤满足三条件；公共交通与步行通勤未达到配对q标准。吸烟包年是当前明确不可估，不是阴性。', '', '## 是否需要合并', '', '1. **结果记录不需要合并**：单病605个唯一键、配对726个唯一键，121个稳定指标ID；未发现同版本重复行。不能把不同疾病、不同样本的IE或PM相加，也不因高相关就删除完整结果。', '2. **NMR可以讨论展示去冗余**：复用实际参与者层Spearman矩阵及每对N，仅作描述性评估。沿用已有“同测量构念、正rho≥0.90、完全链接、簇内平均相关最高作展示代表”的规则，未按IE/PM/P/q挑代表。浓度和组成比例不混簇。', '']
    for label, c in group_counts.items():
        lines.append(f"- {label}：{c['input_markers']}项可形成{c['proposed_display_groups']}个展示组（{c['non_singleton_groups']}个非单项组），主图可少展示{c['proposed_display_reduction']}个高度重叠条目。仅是建议，未实际合并。")
    lines += ['', '相关矩阵来自冻结总体风险集的基线NMR、两两完整记录、未加权Spearman；不是每个结局的选择加权相关，也不能证明独立中介通路。每对实际N见NMR54_ACTUAL_CORRELATION_PAIR_N.csv。其他领域缺少本轮已核对的对应相关矩阵，不凭相似名称或效应大小合并；蛋白不按样本量小而删除或虚增N。', '', '3. **若要统一图例/重复名称，只合并展示映射，不合并统计效应**。121项主表、605单病、726配对总表及固定校正家族始终保留。用户尚未确认展示选择，本轮不生成四张PM图。', '', '## 判据与版本说明', '', '- TRUE/FALSE只用于已有可核对数值的判据；NOT_EVALUABLE_CURRENT_RECORD表示当前记录已经核对但估计不可得，不等于不显著/方向不符；若未找到同版本记录，应使用CURRENT_VERSION_NOT_VERIFIED（当前版本尚未核对），本次为0。', '- 条件1采用单病IE的q_BH_global605_exploratory；条件2采用IHD−T2D ΔIE点值>0；条件3采用配对q_BH_within_contrast121_exploratory。其他q列全部保留，不替换或重算。', '- P沿用冻结Type-7百分位CI反演、1/有效B网格；原有P和所有q原字符串复制。两病及配对区间、实际有效B、失败数、warning-aware和缺失排序敏感性全部保留；部分B区间条件于有效估计，不冒充完整100次。', '- PM_percent已是百分数；Fieller原集合为比例尺度，新_percent字段只做端点×100的单位换算，保留无界与不连通集合。无界PM不因移除图形符号而变有界。', '- 所有当前IHD−T2D可估比较均保存实际父簇抽样证明；目标/权重兼容性原字段保留。未核对原始个体数据或重新计算模型，本次核对范围是已认证发布包及完整注册/任务矩阵。', '', '## 全部121项逐项说明', '', '|稳定ID|中文名称|条件1|条件2|条件3|三条件|说明|', '|---|---|---|---|---|---|---|']
    for r in master:
        lines.append('|' + '|'.join((str(r[k]).replace('|', '/') for k in ['indicator_id', 'indicator_name_zh', 'at_least_one_IE_q605_lt005', 'Delta_IE_positive', 'Delta_IE_q121_lt005', 'all_three_conditions', 'three_condition_reason_zh'])) + '|')
    lines += ['', '## 文件', '', '- B100_121_IHD_T2D_MASTER.csv：121项，完整保留两病及主要配对原字段。', '- B100_121_ALL_DISEASES_MASTER.csv：121项宽表，完整展开五病及六组比较。', '- B100_121_COMPLETE_RESULTS.csv / B100_SIX_CONTRASTS_RESULTS.csv：605 / 726行当前原表的逐字节副本。', '- B100_121_FILTER_EXPLANATIONS.csv：121×6=726条判据记录，每条有真/假/不可判状态及实际取值。', '- NMR_DISPLAY_REDUNDANCY_PROPOSAL_ONLY.csv：仅展示去冗余建议；未改动任何结果。', '- SUMMARY_AND_PROVENANCE.json / FINAL_QA.json：版本、SHA、计数与未改变P/q证明。', '', '没有重拟合、抽样、重算P/q、生成PM图、改变候选或自动派发任务。']
    txt(out / 'REVIEW_REPORT_ZH.md', '\n'.join(lines) + '\n')
    txt(out / 'RUN_COMMAND.txt', f'"{sys.executable}" "{Path(__file__).resolve()}"\n')
    exported = read(out / 'B100_121_IHD_T2D_MASTER.csv')
    for row in exported:
        i = row['indicator_id']
        for prefix, original in [('T2D', s[i, '250.2']), ('IHD', s[i, '411']), ('paired', p[i, '411_minus_250.2'])]:
            for k, v in original.items():
                if k not in ['indicator_id', 'domain']:
                    assert row[prefix + '_' + k] == v, ('CHANGED_ORIGINAL_VALUE', i, prefix, k)
    assert len(exported) == len({r['indicator_id'] for r in exported}) == 121
    assert all((sha(Path(path)) == digest for path, digest in hashes.items()))
    jwrite(out / 'FINAL_QA.json', {'passed': True, 'registered_ids_from_original_registry': 121, 'master_rows': 121, 'single_rows': 605, 'paired_rows': 726, 'filter_condition_rows': len(filters), 'raw_copies_byte_identical': True, 'all_two_disease_and_pair_original_columns_string_identical': True, 'all_five_disease_and_six_contrast_wide_columns_present': True, 'no_old_cache_fill': True, 'no_filtered_list_as_universe': True, 'all_current_records_verified': True, 'missing_estimate_not_misclassified_as_FALSE': True, 'source_hashes_unchanged': True, 'models_fit': 0, 'P_q_recalculated': False, 'PM_figures_generated': 0, 'actual_merges_applied': 0})
    jwrite(out / 'FINAL_MANIFEST.json', {str(f.relative_to(out)): {'bytes': f.stat().st_size, 'sha256': sha(f)} for f in sorted(out.rglob('*')) if f.is_file()})
    print(json.dumps({'output': str(out), 'counts': counts, 'domains': domains, 'NMR': group_counts, 'models_fit': 0, 'P_q_recalculated': False, 'PM_figures_generated': 0}, ensure_ascii=False))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Export the registered 121 from the current complete B100 release, not a shortlist.\n\nAggregate-only joins, Boolean selection audit and descriptive NMR display-redundancy\nassessment. No estimation, P/q calculation, model access, plotting or task dispatch.\n"""\nfrom __future__ import annotations\nimport argparse\nimport ast\nimport collections\nimport csv\nfrom datetime import datetime\nimport hashlib\nimport itertools\nimport json\nimport math\nfrom pathlib import Path\nimport re\nimport shutil\nimport sys\nfrom zoneinfo import ZoneInfo\nBASE = next((p for p in Path(__file__).resolve().parents if (p / \'scripts/birth\').is_dir() and (p / \'reports/birth\').is_dir()))\nEPS = {\'250\': \'DM\', \'250.2\': \'T2D\', \'411\': \'IHD\', \'411.3\': \'Angina\', \'411.4\': \'CAS\'}\nCONTRASTS = [f\'{a}_minus_{b}\' for a in [\'411\', \'411.3\', \'411.4\'] for b in [\'250\', \'250.2\']]\nMETHOD = birth_private_setting(\'PRIVATE_7AA11EFF31702D72C96A44823457AB674B93D8CC8A1AD7965B72932F60DA\')\nSOURCE = birth_private_setting(\'PRIVATE_1095A231CC6F1187B50096D4E8F4F675170A4FA973A9BB301E2CDD931030\')\nVERSIONS = {\'B100_ORIGINAL_ESTIMATOR_WARNING_AWARE_RELEASE_V1\', \'B100_MISSING_INFERENCE_COMPLETION_20261030_V1_ORIGINAL_ESTIMATOR\'}\nFOCUS = [\'qualification\', \'income\', \'age_education\', \'birth_weight\', \'insomnia\', \'townsend\', \'smoking\', \'commute_car\', \'commute_public_transport\', \'commute_walking\', \'packyears\', \'maternal_smoking\']\nUNKNOWN = \'NOT_EVALUABLE_CURRENT_RECORD\'\nUNVERIFIED = \'CURRENT_VERSION_NOT_VERIFIED\'\nFLAG_LABELS = {\'T2D_IE_q605_lt005\': \'T2D单病IE：当前q605<0.05\', \'IHD_IE_q605_lt005\': \'IHD单病IE：当前q605<0.05\', \'at_least_one_IE_q605_lt005\': \'条件1：两病至少一病IE q605<0.05\', \'Delta_IE_positive\': \'条件2：IHD−T2D的ΔIE>0\', \'Delta_IE_q121_lt005\': \'条件3：配对ΔIE当前q121<0.05\', \'all_three_conditions\': \'上述三条件同时满足（仅展示筛选，不是结果是否存在）\'}\n\ndef read(path):\n    with path.open(encoding=\'utf-8-sig\', newline=\'\') as f:\n        return list(csv.DictReader(f))\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef write(path, rows):\n    assert rows\n    fields = list(dict.fromkeys((k for r in rows for k in r)))\n    with path.open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as f:\n        w = csv.DictWriter(f, fieldnames=fields)\n        w.writeheader()\n        w.writerows(rows)\n\ndef jwrite(path, obj):\n    with path.open(\'x\', encoding=\'utf-8\') as f:\n        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)\n\ndef txt(path, value):\n    with path.open(\'x\', encoding=\'utf-8\') as f:\n        f.write(value)\n\ndef number(value):\n    try:\n        x = float(value)\n        return x if math.isfinite(x) else None\n    except (ValueError, TypeError):\n        return None\n\ndef state(value, predicate):\n    n = number(value)\n    return UNKNOWN if n is None else \'TRUE\' if predicate(n) else \'FALSE\'\n\ndef either(a, b):\n    if \'TRUE\' in [a, b]:\n        return \'TRUE\'\n    return \'FALSE\' if a == b == \'FALSE\' else UNKNOWN\n\ndef all_flags(flags):\n    if any((v not in [\'TRUE\', \'FALSE\'] for v in flags)):\n        return UNKNOWN\n    return \'TRUE\' if all((v == \'TRUE\' for v in flags)) else \'FALSE\'\n\ndef percent_set(s):\n    return re.sub(\'[-+]?(?:\\\\d+\\\\.?\\\\d*|\\\\.\\\\d+)(?:[eE][-+]?\\\\d+)?\', lambda m: format(float(m.group()) * 100, \'.12g\'), s)\n\ndef close(a, b):\n    return abs(float(a) - float(b)) <= 1e-12 + 1e-09 * max(abs(float(a)), abs(float(b)))\n\ndef pretty(v, digits=4):\n    n = number(v)\n    return \'NA\' if n is None else f\'{n:.{digits}g}\'\n\ndef ci(row, prefix):\n    return f"{pretty(row.get(prefix + \'_IE_per_1000\'))} [{pretty(row.get(prefix + \'_IE_CI_low_per_1000\'))}, {pretty(row.get(prefix + \'_IE_CI_high_per_1000\'))}]"\n\ndef run(args):\n    source = args.source\n    out = args.out\n    if out is None:\n        v = 1\n        while (BASE / f\'reports/birth/b100_master121_audit_20261030_v{v}\').exists():\n            v += 1\n        out = BASE / f\'reports/birth/b100_master121_audit_20261030_v{v}\'\n    assert not out.exists(), \'PRESERVE_EXISTING_OUTPUT:\' + str(out)\n    project = BASE / \'scripts/birth/border_fullfollowup_b100_20261027\'\n    registry_path = project / \'input/APPROVED_121_NO_CMR.csv\'\n    matrix_path = project / \'input/TASK_MATRIX_605.csv\'\n    single_path = source / \'B100_121_COMPLETE_RESULTS.csv\'\n    pair_path = source / \'B100_SIX_CONTRASTS_RESULTS.csv\'\n    release_root = source.parent\n    table_manifest = source / \'TABLE_MANIFEST.json\'\n    receipt = release_root / \'PUBLICATION_RECEIPT_v2.json\'\n    transfer = release_root / \'PUBLICATION_TRANSFER_MANIFEST_v2.json\'\n    qa_path = source / \'INDEPENDENT_AGGREGATE_QA.json\'\n    corr_root = args.correlation\n    corr_names = [\'NMR_CORRELATION_MATRIX.csv\', \'NMR_CORRELATION_PAIR_N.csv\', \'NMR_CORRELATION_VARIABLE_AUDIT.csv\', \'NMR_BIOCHEMICAL_CLASSIFICATION.csv\', \'NMR_CORRELATION_SUMMARY.json\', \'VALIDATION_RECEIPT.json\', \'CORRELATION_EXECUTION_RECEIPT.json\']\n    cluster_source = BASE / \'scripts/birth/b100_pm_4domain_sankey_20261030_v1/build.py\'\n    inputs = [registry_path, matrix_path, single_path, pair_path, table_manifest, receipt, transfer, qa_path, cluster_source] + [corr_root / n for n in corr_names]\n    for p in inputs:\n        assert p.is_file(), \'MISSING_INPUT:\' + str(p)\n    hashes = {str(p): sha(p) for p in inputs}\n    single_digest, pair_digest = (hashes[str(single_path)], hashes[str(pair_path)])\n    manifest = json.loads(table_manifest.read_text())\n    transferred = json.loads(transfer.read_text())\n    publication = json.loads(receipt.read_text())\n    for p in [single_path, pair_path]:\n        assert sha(p) == manifest[p.name][\'sha256\'] == transferred[p.name][\'sha256\']\n        assert p.stat().st_size == manifest[p.name][\'bytes\'] == transferred[p.name][\'bytes\']\n    assert publication[\'complete\'][\'registered_indicators\'] == 121\n    assert publication[\'complete\'][\'delivery_complete\']\n    assert json.loads(qa_path.read_text())[\'passed\']\n    registry = read(registry_path)\n    reg = {r[\'indicator_id\']: r for r in registry}\n    matrix = read(matrix_path)\n    singles = read(single_path)\n    pairs = read(pair_path)\n    s = {(r[\'indicator_id\'], r[\'endpoint_id\']): r for r in singles}\n    p = {(r[\'indicator_id\'], r[\'contrast\']): r for r in pairs}\n    assert len(registry) == len(reg) == 121\n    expected_single = set(itertools.product(reg, EPS))\n    expected_pair = set(itertools.product(reg, CONTRASTS))\n    assert len(singles) == len(s) == 605 and set(s) == expected_single\n    assert len(pairs) == len(p) == 726 and set(p) == expected_pair\n    assert {(r[\'indicator_id\'], r[\'endpoint_id\']) for r in matrix} == expected_single\n    assert all((r[\'source_method_hash\'] == METHOD and r[\'source_hash\'] == SOURCE and (r[\'version\'] in VERSIONS) and (r[\'horizon_years\'] == \'11\') and (r[\'model_followup\'] == \'FULL_FOLLOWUP\') and (r[\'bootstrap_planned\'] == \'100\') for r in singles))\n    for r in singles:\n        assert r[\'domain\'] == reg[r[\'indicator_id\']][\'category\']\n        if number(r[\'PM_percent\']) is not None:\n            assert close(r[\'PM_percent\'], 100 * float(r[\'IE\']) / float(r[\'TE\']))\n    for r in pairs:\n        assert r[\'horizon_years\'] == \'11\' and r[\'bootstrap_planned\'] == \'100\'\n        a, b = (s[r[\'indicator_id\'], r[\'endpoint_A\']], s[r[\'indicator_id\'], r[\'endpoint_B\']])\n        if number(r[\'Delta_IE\']) is not None:\n            assert close(r[\'IE_A_per_1000\'], a[\'IE_per_1000\'])\n            assert close(r[\'IE_B_per_1000\'], b[\'IE_per_1000\'])\n            assert close(r[\'Delta_IE\'], float(a[\'IE\']) - float(b[\'IE\']))\n            assert r[\'pairing_evidence\'] in {\'ACTUAL_SHARED_PARENT_DRAW_HASH_VERIFIED_FOR_EACH_VALID_REPEAT\', \'EXACT_OLD_ACTUAL_PARENT_DRAW_HASH_PER_ID_VERIFIED;NEW_A_ORIGIN_RETAINED\'}\n    export_version = \'B100_MASTER121_AUDIT_20261030_V1_CURRENT_COMPLETED_RELEASE\'\n    master = []\n    filters = []\n    records = []\n    for rank, entry in enumerate(registry, 1):\n        i = entry[\'indicator_id\']\n        t, h, d = (s[i, \'250.2\'], s[i, \'411\'], p[i, \'411_minus_250.2\'])\n        flags = {\'T2D_IE_q605_lt005\': state(t[\'q_BH_global605_exploratory\'], lambda x: x < 0.05), \'IHD_IE_q605_lt005\': state(h[\'q_BH_global605_exploratory\'], lambda x: x < 0.05), \'Delta_IE_positive\': state(d[\'Delta_IE\'], lambda x: x > 0), \'Delta_IE_q121_lt005\': state(d[\'q_BH_within_contrast121_exploratory\'], lambda x: x < 0.05)}\n        flags[\'at_least_one_IE_q605_lt005\'] = either(flags[\'T2D_IE_q605_lt005\'], flags[\'IHD_IE_q605_lt005\'])\n        flags[\'all_three_conditions\'] = all_flags([flags[\'at_least_one_IE_q605_lt005\'], flags[\'Delta_IE_positive\'], flags[\'Delta_IE_q121_lt005\']])\n        reason = \'；\'.join((FLAG_LABELS[k] + \'=FALSE\' for k in [\'at_least_one_IE_q605_lt005\', \'Delta_IE_positive\', \'Delta_IE_q121_lt005\'] if flags[k] == \'FALSE\'))\n        if flags[\'all_three_conditions\'] == UNKNOWN:\n            reason = \'已核对当前版本，但无可估计数值，不能判成不显著或方向不符；\' + d[\'failure_reason\']\n        elif not reason:\n            reason = \'三条件均TRUE；仅符合既定展示筛选，不表示独立确认或PM有界\'\n        record_status = \'已核对当前版本：有估计\' if number(d[\'Delta_IE\']) is not None else \'已核对当前版本：不可估，原因保留\'\n        row = {\'registry_order\': rank, \'indicator_id\': i, \'indicator_name_zh\': entry[\'label_zh\'], \'domain\': entry[\'category\'], \'indicator_type\': entry[\'type\'], \'current_record_status\': record_status, \'current_version_verified\': \'TRUE\', \'export_version\': export_version, \'method_hash_verified\': METHOD, \'source_hash_verified\': SOURCE, \'single_source_sha256\': single_digest, \'paired_source_sha256\': pair_digest, \'T2D_source_row_number\': singles.index(t) + 2, \'IHD_source_row_number\': singles.index(h) + 2, \'paired_source_row_number\': pairs.index(d) + 2, \'three_condition_reason_zh\': reason, **flags}\n        for prefix, r in [(\'T2D\', t), (\'IHD\', h), (\'paired\', d)]:\n            row.update({prefix + \'_\' + k: v for k, v in r.items() if k not in [\'indicator_id\', \'domain\']})\n        for prefix, r in [(\'T2D\', t), (\'IHD\', h)]:\n            row[prefix + \'_PM_Fieller_set_percent\'] = percent_set(r[\'PM_Fieller_set\'])\n            row[prefix + \'_PM_Fieller_original_set_unit\'] = \'ratio; PM_percent is already percent\'\n        master.append(row)\n        raw_values = {\'T2D_IE_q605_lt005\': t[\'q_BH_global605_exploratory\'], \'IHD_IE_q605_lt005\': h[\'q_BH_global605_exploratory\'], \'at_least_one_IE_q605_lt005\': json.dumps([t[\'q_BH_global605_exploratory\'], h[\'q_BH_global605_exploratory\']]), \'Delta_IE_positive\': d[\'Delta_IE\'], \'Delta_IE_q121_lt005\': d[\'q_BH_within_contrast121_exploratory\'], \'all_three_conditions\': json.dumps([flags[k] for k in [\'at_least_one_IE_q605_lt005\', \'Delta_IE_positive\', \'Delta_IE_q121_lt005\']])}\n        for k, label in FLAG_LABELS.items():\n            filters.append({\'indicator_id\': i, \'indicator_name_zh\': entry[\'label_zh\'], \'domain\': entry[\'category\'], \'condition\': k, \'condition_zh\': label, \'actual_current_value\': raw_values[k], \'condition_state\': flags[k], \'current_record_status\': record_status, \'source_single_sha256\': single_digest, \'source_paired_sha256\': pair_digest, \'reason_if_not_evaluable\': d[\'failure_reason\'] if flags[k] == UNKNOWN else \'\'})\n        for ep in EPS:\n            r = s[i, ep]\n            records.append({\'indicator_id\': i, \'endpoint_id\': ep, \'source_row_number\': singles.index(r) + 2, \'record_found\': \'TRUE\', \'version_verified\': \'TRUE\', \'version\': r[\'version\'], \'source_method_hash\': r[\'source_method_hash\'], \'source_hash\': r[\'source_hash\'], \'source_file_sha256\': single_digest, \'result_status\': r[\'result_status\'], \'failure_reason\': r[\'failure_reason\']})\n    by_id = {r[\'indicator_id\']: r for r in master}\n    all_disease = []\n    for row in master:\n        i = row[\'indicator_id\']\n        wide = {k: row[k] for k in [\'registry_order\', \'indicator_id\', \'indicator_name_zh\', \'domain\', \'indicator_type\', \'current_record_status\', \'export_version\', \'method_hash_verified\']}\n        for ep, name in EPS.items():\n            wide.update({name + \'_\' + k: v for k, v in s[i, ep].items() if k not in [\'indicator_id\', \'domain\']})\n        for contrast in CONTRASTS:\n            wide.update({contrast + \'_\' + k: v for k, v in p[i, contrast].items() if k not in [\'indicator_id\', \'domain\']})\n        all_disease.append(wide)\n    corr_summary = json.loads((corr_root / \'NMR_CORRELATION_SUMMARY.json\').read_text())\n    validation = json.loads((corr_root / \'VALIDATION_RECEIPT.json\').read_text())\n    assert corr_summary[\'registry_sha256\'] == sha(registry_path)\n    assert corr_summary[\'schema\'] == \'NMR_REAL_PARTICIPANT_SPEARMAN_V1\'\n    assert corr_summary[\'nmr_count\'] == 54 and validation[\'status\'] == \'PASS\'\n    assert corr_summary[\'models_run\'] == 0 and corr_summary[\'participant_rows_exported\'] == 0\n    rho = {r[\'indicator_id\']: {k: float(v) for k, v in r.items() if k != \'indicator_id\'} for r in read(corr_root / \'NMR_CORRELATION_MATRIX.csv\')}\n    ns = {r[\'indicator_id\']: {k: int(v) for k, v in r.items() if k != \'indicator_id\'} for r in read(corr_root / \'NMR_CORRELATION_PAIR_N.csv\')}\n    bio = {r[\'indicator_id\']: r for r in read(corr_root / \'NMR_BIOCHEMICAL_CLASSIFICATION.csv\')}\n    cov = {r[\'indicator_id\']: int(r[\'N_finite\']) for r in read(corr_root / \'NMR_CORRELATION_VARIABLE_AUDIT.csv\')}\n    nmr = [r[\'indicator_id\'] for r in master if r[\'domain\'] == \'NMR\']\n    assert len(nmr) == 54 and set(nmr) == set(rho) == set(ns) == set(bio) == set(cov)\n    for a, b in itertools.product(nmr, repeat=2):\n        assert math.isfinite(rho[a][b]) and -1 <= rho[a][b] <= 1\n        assert abs(rho[a][b] - rho[b][a]) < 1e-12\n        assert ns[a][b] == ns[b][a] and 0 < ns[a][b] <= min(cov[a], cov[b])\n    tree = ast.parse(cluster_source.read_text())\n    fn = next((n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == \'cluster_nmr\'))\n    snippet = ast.get_source_segment(cluster_source.read_text(), fn)\n    env = {\'collections\': collections, \'itertools\': itertools, \'math\': math}\n    exec(compile(snippet, str(cluster_source), \'exec\'), env)\n    order = {r[\'indicator_id\']: r[\'registry_order\'] for r in master}\n    eligible_nmr = [i for i in nmr if by_id[i][\'all_three_conditions\'] == \'TRUE\']\n    proposals, merges = ([], [])\n    group_counts = {}\n    for label, ids in [(\'ALL_REGISTERED_NMR54\', nmr), (\'THREE_CONDITION_ELIGIBLE_NMR\', eligible_nmr)]:\n        clusters, log = env[\'cluster_nmr\'](ids, bio, rho, cov, order)\n        group_counts[label] = {\'input_markers\': len(ids), \'proposed_display_groups\': len(clusters), \'proposed_display_reduction\': len(ids) - len(clusters), \'non_singleton_groups\': sum((len(cl[\'members\']) > 1 for cl in clusters))}\n        for cl in clusters:\n            for i in cl[\'members\']:\n                proposals.append({\'assessment_scope\': label, \'cluster_id\': cl[\'cluster_id\'], \'indicator_id\': i, \'indicator_name_zh\': reg[i][\'label_zh\'], \'quantity_kind\': cl[\'quantity_kind\'], \'members\': \';\'.join(cl[\'members\']), \'representative_proposal_only\': cl[\'representative\'], \'mean_within_cluster_rho\': cl[\'mean_rho\'][i], \'min_within_cluster_rho\': cl[\'minimum_within_rho\'], \'marker_finite_N\': cov[i], \'proposal_only_not_applied\': True, \'IE_PM_never_pooled\': True})\n        merges.extend(({\'assessment_scope\': label, **r, \'proposal_only_not_applied\': True} for r in log))\n    pairwise = [{\'indicator_A\': a, \'indicator_B\': b, \'rho\': rho[a][b], \'actual_pair_N\': ns[a][b], \'same_quantity_kind\': bio[a][\'quantity_kind\'] == bio[b][\'quantity_kind\'], \'positive_rho_ge090\': rho[a][b] >= 0.9, \'both_three_condition_TRUE\': a in eligible_nmr and b in eligible_nmr, \'correlation_source\': \'Existing real participant-level descriptive audit, not result correlations\'} for a, b in itertools.combinations(nmr, 2)]\n    out.mkdir(parents=True)\n    (out / \'source\').mkdir()\n    (out / \'evidence\').mkdir()\n    (out / \'code\').mkdir()\n    for pth in [single_path, pair_path]:\n        shutil.copy2(pth, out / pth.name)\n        assert sha(out / pth.name) == sha(pth)\n    for pth in [registry_path, matrix_path, table_manifest, receipt, transfer, qa_path]:\n        shutil.copy2(pth, out / \'source\' / pth.name)\n    for name in corr_names:\n        shutil.copy2(corr_root / name, out / \'evidence\' / name)\n    shutil.copy2(Path(__file__), out / \'code/export.py\')\n    txt(out / \'code/reused_cluster_nmr.py.txt\', snippet)\n    write(out / \'B100_121_IHD_T2D_MASTER.csv\', master)\n    write(out / \'B100_121_ALL_DISEASES_MASTER.csv\', all_disease)\n    write(out / \'B100_121_FILTER_EXPLANATIONS.csv\', filters)\n    write(out / \'CURRENT_VERSION_605_RECORD_AUDIT.csv\', records)\n    priority = [by_id[i] for i in FOCUS]\n    write(out / \'PRIORITY_NAMED_CURRENT_RECORDS.csv\', priority)\n    write(out / \'NMR_DISPLAY_REDUNDANCY_PROPOSAL_ONLY.csv\', proposals)\n    if merges:\n        write(out / \'NMR_CLUSTER_STEPS_PROPOSAL_ONLY.csv\', merges)\n    write(out / \'NMR54_ACTUAL_CORRELATION_PAIR_N.csv\', pairwise)\n    counts = {\'registered\': len(master), \'paired_estimable\': sum((number(r[\'paired_Delta_IE\']) is not None for r in master)), \'unverified_current_records\': 0, \'single_cells\': len(singles), \'paired_rows\': len(pairs), \'single_key_duplicates\': len(singles) - len(s), \'pair_key_duplicates\': len(pairs) - len(p)}\n    for key in FLAG_LABELS:\n        counts[key] = dict(collections.Counter((r[key] for r in master)))\n    domains = []\n    for domain in [\'危险因素\', \'社会、行为及生命历程\', \'NMR\', \'蛋白组\']:\n        rows = [r for r in master if r[\'domain\'] == domain]\n        domains.append({\'domain\': domain, \'registered\': len(rows), \'paired_estimable\': sum((number(r[\'paired_Delta_IE\']) is not None for r in rows)), \'at_least_one_IE_q605_lt005\': sum((r[\'at_least_one_IE_q605_lt005\'] == \'TRUE\' for r in rows)), \'IE_supported_and_Delta_positive\': sum((r[\'at_least_one_IE_q605_lt005\'] == r[\'Delta_IE_positive\'] == \'TRUE\' for r in rows)), \'all_three_TRUE\': sum((r[\'all_three_conditions\'] == \'TRUE\' for r in rows))})\n    write(out / \'COUNTS_BY_DOMAIN.csv\', domains)\n    summary = {\'created_beijing\': datetime.now(ZoneInfo(\'Asia/Shanghai\')).isoformat(), \'source_local_release\': str(source), \'source_remote_release_receipt\': publication[\'remote\'], \'source_method_hash\': METHOD, \'source_source_hash\': SOURCE, \'source_result_versions\': dict(collections.Counter((r[\'version\'] for r in singles))), \'universe_origin\': str(registry_path), \'never_reverse_inferred_from_77_or_64\': True, \'counts\': counts, \'domains\': domains, \'NMR_display_redundancy_assessment_only\': group_counts, \'models_fit\': 0, \'bootstrap_run\': 0, \'P_q_recalculated\': False, \'PM_figures_generated\': 0, \'actual_merges_applied\': 0, \'old_result_fallback_used\': False, \'source_sha256\': hashes, \'old_results_unchanged\': True}\n    jwrite(out / \'SUMMARY_AND_PROVENANCE.json\', summary)\n    lines = [\'# 当前完整B100：121项主表、筛选和合并评估（未出图）\', \'\', \'仅从当前完整605单病表和726配对表读取结果，以原冻结121项稳定ID注册表构建全集。没有读取77/64筛后名单，没有用旧缓存填充结果。\', \'\', f\'已核对605/605单病、726/726配对记录：主表121项，其中110项有IHD−T2D估计、11项有当前不可估记录；尚未核对的当前记录0项。重复主键0项，不需要为清理重复结果而合并。\', \'\', \'两种行版本标签分别为原有效结果565格及同估计器补全40格，均属于当前完整发布；统计method/source哈希一致。不是把旧的无蛋白区间版本混入。\', \'\', \'## 先看多少项\', \'\', \'|领域|完整注册|可估配对|至少一病IE q605<0.05|再要求ΔIE>0|三条件全部满足|\', \'|---|---:|---:|---:|---:|---:|\']\n    for r in domains:\n        lines.append(\'|\' + \'|\'.join((str(r[k]) for k in [\'domain\', \'registered\', \'paired_estimable\', \'at_least_one_IE_q605_lt005\', \'IE_supported_and_Delta_positive\', \'all_three_TRUE\'])) + \'|\')\n    lines += [\'\', \'共96项至少一病IE q605<0.05；其中77项ΔIE>0；加入配对q121<0.05后为64项。这些数字是筛选层级，不是把其余指标从结果全集删除。\', \'\', \'## 点名指标当前记录\', \'\', \'IE和ΔIE均为每千人11年风险差；下表q为各自当前固定家族。\', \'\', \'|指标（稳定ID）|T2D IE [95%CI]|IHD IE [95%CI]|ΔIE [95%CI]|T2D / IHD q605|ΔIE P / q121|有效B T2D/IHD/配对|条件1/2/3|\', \'|---|---|---|---|---|---|---|---|\']\n    for r in priority:\n        delta = f"{pretty(r[\'paired_Delta_IE_per_1000\'])} [{pretty(r[\'paired_Delta_IE_CI_low_per_1000\'])}, {pretty(r[\'paired_Delta_IE_CI_high_per_1000\'])}]"\n        lines.append(f"|{r[\'indicator_name_zh\']} (`{r[\'indicator_id\']}`)|{ci(r, \'T2D\')}|{ci(r, \'IHD\')}|{delta}|{pretty(r[\'T2D_q_BH_global605_exploratory\'])} / {pretty(r[\'IHD_q_BH_global605_exploratory\'])}|{pretty(r[\'paired_P_Delta_IE\'])} / {pretty(r[\'paired_q_BH_within_contrast121_exploratory\'])}|{r[\'T2D_bootstrap_valid\']}/{r[\'IHD_bootstrap_valid\']}/{r[\'paired_paired_valid\']}|{r[\'at_least_one_IE_q605_lt005\']} / {r[\'Delta_IE_positive\']} / {r[\'Delta_IE_q121_lt005\']}|")\n    lines += [\'\', \'点名记录的完整单病P、两种q、PM、Fieller比例集合及百分数集合、版本和原因均在PRIORITY_NAMED_CURRENT_RECORDS.csv；此处没有用表格缩略值作筛选。\', \'\', \'学历、收入、教育结束年龄、出生体重、失眠和Townsend均有当前估计，不应说成“没有结果”。学历/收入/教育结束年龄的ΔIE点值为负；出生体重/失眠/Townsend为正，但本次配对q未达0.05。吸烟与汽车通勤满足三条件；公共交通与步行通勤未达到配对q标准。吸烟包年是当前明确不可估，不是阴性。\', \'\', \'## 是否需要合并\', \'\', \'1. **结果记录不需要合并**：单病605个唯一键、配对726个唯一键，121个稳定指标ID；未发现同版本重复行。不能把不同疾病、不同样本的IE或PM相加，也不因高相关就删除完整结果。\', \'2. **NMR可以讨论展示去冗余**：复用实际参与者层Spearman矩阵及每对N，仅作描述性评估。沿用已有“同测量构念、正rho≥0.90、完全链接、簇内平均相关最高作展示代表”的规则，未按IE/PM/P/q挑代表。浓度和组成比例不混簇。\', \'\']\n    for label, c in group_counts.items():\n        lines.append(f"- {label}：{c[\'input_markers\']}项可形成{c[\'proposed_display_groups\']}个展示组（{c[\'non_singleton_groups\']}个非单项组），主图可少展示{c[\'proposed_display_reduction\']}个高度重叠条目。仅是建议，未实际合并。")\n    lines += [\'\', \'相关矩阵来自冻结总体风险集的基线NMR、两两完整记录、未加权Spearman；不是每个结局的选择加权相关，也不能证明独立中介通路。每对实际N见NMR54_ACTUAL_CORRELATION_PAIR_N.csv。其他领域缺少本轮已核对的对应相关矩阵，不凭相似名称或效应大小合并；蛋白不按样本量小而删除或虚增N。\', \'\', \'3. **若要统一图例/重复名称，只合并展示映射，不合并统计效应**。121项主表、605单病、726配对总表及固定校正家族始终保留。用户尚未确认展示选择，本轮不生成四张PM图。\', \'\', \'## 判据与版本说明\', \'\', \'- TRUE/FALSE只用于已有可核对数值的判据；NOT_EVALUABLE_CURRENT_RECORD表示当前记录已经核对但估计不可得，不等于不显著/方向不符；若未找到同版本记录，应使用CURRENT_VERSION_NOT_VERIFIED（当前版本尚未核对），本次为0。\', \'- 条件1采用单病IE的q_BH_global605_exploratory；条件2采用IHD−T2D ΔIE点值>0；条件3采用配对q_BH_within_contrast121_exploratory。其他q列全部保留，不替换或重算。\', \'- P沿用冻结Type-7百分位CI反演、1/有效B网格；原有P和所有q原字符串复制。两病及配对区间、实际有效B、失败数、warning-aware和缺失排序敏感性全部保留；部分B区间条件于有效估计，不冒充完整100次。\', \'- PM_percent已是百分数；Fieller原集合为比例尺度，新_percent字段只做端点×100的单位换算，保留无界与不连通集合。无界PM不因移除图形符号而变有界。\', \'- 所有当前IHD−T2D可估比较均保存实际父簇抽样证明；目标/权重兼容性原字段保留。未核对原始个体数据或重新计算模型，本次核对范围是已认证发布包及完整注册/任务矩阵。\', \'\', \'## 全部121项逐项说明\', \'\', \'|稳定ID|中文名称|条件1|条件2|条件3|三条件|说明|\', \'|---|---|---|---|---|---|---|\']\n    for r in master:\n        lines.append(\'|\' + \'|\'.join((str(r[k]).replace(\'|\', \'/\') for k in [\'indicator_id\', \'indicator_name_zh\', \'at_least_one_IE_q605_lt005\', \'Delta_IE_positive\', \'Delta_IE_q121_lt005\', \'all_three_conditions\', \'three_condition_reason_zh\'])) + \'|\')\n    lines += [\'\', \'## 文件\', \'\', \'- B100_121_IHD_T2D_MASTER.csv：121项，完整保留两病及主要配对原字段。\', \'- B100_121_ALL_DISEASES_MASTER.csv：121项宽表，完整展开五病及六组比较。\', \'- B100_121_COMPLETE_RESULTS.csv / B100_SIX_CONTRASTS_RESULTS.csv：605 / 726行当前原表的逐字节副本。\', \'- B100_121_FILTER_EXPLANATIONS.csv：121×6=726条判据记录，每条有真/假/不可判状态及实际取值。\', \'- NMR_DISPLAY_REDUNDANCY_PROPOSAL_ONLY.csv：仅展示去冗余建议；未改动任何结果。\', \'- SUMMARY_AND_PROVENANCE.json / FINAL_QA.json：版本、SHA、计数与未改变P/q证明。\', \'\', \'没有重拟合、抽样、重算P/q、生成PM图、改变候选或自动派发任务。\']\n    txt(out / \'REVIEW_REPORT_ZH.md\', \'\\n\'.join(lines) + \'\\n\')\n    txt(out / \'RUN_COMMAND.txt\', f\'"{sys.executable}" "{Path(__file__).resolve()}"\\n\')\n    exported = read(out / \'B100_121_IHD_T2D_MASTER.csv\')\n    for row in exported:\n        i = row[\'indicator_id\']\n        for prefix, original in [(\'T2D\', s[i, \'250.2\']), (\'IHD\', s[i, \'411\']), (\'paired\', p[i, \'411_minus_250.2\'])]:\n            for k, v in original.items():\n                if k not in [\'indicator_id\', \'domain\']:\n                    assert row[prefix + \'_\' + k] == v, (\'CHANGED_ORIGINAL_VALUE\', i, prefix, k)\n    assert len(exported) == len({r[\'indicator_id\'] for r in exported}) == 121\n    assert all((sha(Path(path)) == digest for path, digest in hashes.items()))\n    jwrite(out / \'FINAL_QA.json\', {\'passed\': True, \'registered_ids_from_original_registry\': 121, \'master_rows\': 121, \'single_rows\': 605, \'paired_rows\': 726, \'filter_condition_rows\': len(filters), \'raw_copies_byte_identical\': True, \'all_two_disease_and_pair_original_columns_string_identical\': True, \'all_five_disease_and_six_contrast_wide_columns_present\': True, \'no_old_cache_fill\': True, \'no_filtered_list_as_universe\': True, \'all_current_records_verified\': True, \'missing_estimate_not_misclassified_as_FALSE\': True, \'source_hashes_unchanged\': True, \'models_fit\': 0, \'P_q_recalculated\': False, \'PM_figures_generated\': 0, \'actual_merges_applied\': 0})\n    jwrite(out / \'FINAL_MANIFEST.json\', {str(f.relative_to(out)): {\'bytes\': f.stat().st_size, \'sha256\': sha(f)} for f in sorted(out.rglob(\'*\')) if f.is_file()})\n    print(json.dumps({\'output\': str(out), \'counts\': counts, \'domains\': domains, \'NMR\': group_counts, \'models_fit\': 0, \'P_q_recalculated\': False, \'PM_figures_generated\': 0}, ensure_ascii=False))\nif __name__ == \'__main__\':\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--source\', type=Path, default=BASE / \'reports/birth/b100_completion_20261030_v1/release\')\n    parser.add_argument(\'--correlation\', type=Path, default=BASE / \'scripts/birth/border_fullfollowup_b1000_16_20261028_v4/input/candidate_evidence\')\n    parser.add_argument(\'--out\', type=Path)\n    run(parser.parse_args())', '<private-stage>', 'exec'), context)
