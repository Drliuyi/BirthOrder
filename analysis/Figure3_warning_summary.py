"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import os
import sys
from pathlib import Path
IDS = ('olink__fcrl6', 'olink__gzma', 'olink__gzmh', 'olink__klrd1')
EPS = ('250', '250.2', '411', '411.3', '411.4')
GOOD = ('VALID', 'WARNING_QUALIFIED_CONDITIONAL_ESTIMATE')
REVIEW_REL = 'reports/birth/b1000_add9_acceptance_20261031_v1/review_20261031_v1'

def require(condition, message):
    if not condition:
        raise ValueError(message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def jsread(path, expected=None):
    path = Path(path)
    require(path.is_file() and path.stat().st_size <= 32 * 1024 ** 2, 'SMALL_AGGREGATE_SOURCE_REQUIRED:' + str(path))
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    require(expected is None or digest == expected, 'AGGREGATE_HASH_CHANGED:' + str(path))
    return (json.loads(content), digest)

def obj(value):
    return value if isinstance(value, dict) else {}

def rows(value):
    return value if isinstance(value, list) else list(value.values()) if isinstance(value, dict) else []

def numeric(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None

def same(a, b, tolerance=1e-12):
    a, b = (numeric(a), numeric(b))
    return a is not None and b is not None and (abs(a - b) <= tolerance + tolerance * max(abs(a), abs(b)))

def identity(row, expected):
    require(all((row.get(k) == v for k, v in expected.items())), 'SEALED_IDENTITY_CHANGED:' + str(expected))

def sealed_result(directory, expected, expected_seal=None):
    directory = Path(directory)
    seal, sealsha = jsread(directory / 'COMPLETE.json', expected_seal)
    require(seal.get('complete') is True, 'INCOMPLETE_SEAL:' + str(directory))
    identity(seal, expected)
    artifacts = seal.get('artifacts_sha256', {})
    require('RESULT.json' in artifacts, 'SEALED_RESULT_NOT_INDEXED:' + str(directory))
    for relative in artifacts:
        require(not Path(relative).is_absolute() and '..' not in Path(relative).parts, 'UNSAFE_SEALED_RELATIVE_PATH')
    result, resultsha = jsread(directory / 'RESULT.json', artifacts['RESULT.json'])
    identity(result, expected)
    return (result, seal, sealsha, resultsha)

def effects_check(row):
    require(all((numeric(row.get(k)) is not None for k in ('R00', 'R10', 'R11', 'TE', 'DE', 'IE'))), 'NONFINITE_EFFECT_OR_RISK')
    require(all((0 <= row[k] <= 1 for k in ('R00', 'R10', 'R11'))), 'MEAN_RISK_OUTSIDE_UNIT_INTERVAL')
    require(same(row['TE'], row['R11'] - row['R00']) and same(row['DE'], row['R10'] - row['R00']) and same(row['IE'], row['R11'] - row['R10']) and same(row['TE'], row['DE'] + row['IE']), 'SAVED_EFFECT_IDENTITY_FAILED')
    require(row['TE'] == 0 and numeric(row.get('PM')) is None or (row['TE'] != 0 and same(row.get('PM'), row['IE'] / row['TE'])), 'PM_NOT_SAME_MODEL_RATIO')

def probability_check(row):
    predictions = obj(row.get('predictions'))
    require(set(predictions) == {'observed', 'R00', 'R10', 'R11'}, 'FOUR_PREDICTION_AUDITS_REQUIRED')
    for scenario, pred in predictions.items():
        for family in ('original', 'warm'):
            audit = pred[family]
            require(audit.get('legal') is True and audit.get('nonfinite') == 0 and (audit.get('outside_unit_interval') == 0), 'INDIVIDUAL_PROBABILITY_LEGALITY_FAILED')
            require(numeric(audit.get('max_abs_closure')) is not None and audit['max_abs_closure'] < 1e-08, 'INDIVIDUAL_PROBABILITY_CLOSURE_FAILED')
            require(all((0 <= v <= 1 for v in obj(audit.get('state_min')).values())) and all((0 <= v <= 1 for v in obj(audit.get('state_max')).values())), 'INDIVIDUAL_PROBABILITY_RANGE_FAILED')
        if scenario != 'observed':
            selected = row.get(scenario, obj(row.get('original_effects')).get(scenario))
            require(same(selected, pred['original_weighted_disease_risk']), 'ORIGINAL_FAMILY_RISK_NOT_SELECTED')
    require(row.get('native_validation', {}).get('original', {}).get('passed') is True and row.get('native_validation', {}).get('warm', {}).get('passed') is True, 'NATIVE_VALIDATION_FAILED')
    diagnostics = obj(row.get('coefficient_diagnostics'))
    require(set(diagnostics) == {'disease', 'death'}, 'BOTH_CAUSE_DIAGNOSTICS_REQUIRED')
    for diagnostic in diagnostics.values():
        stable = diagnostic.get('core_stability_passed', diagnostic.get('core_stability_descriptive_rule_passed'))
        require(stable is True and diagnostic.get('core_positive_information_both') is True, 'CORE_INFORMATION_OR_STABILITY_FAILED')
    if 'individual_probability_1e9_is_acceptance_gate' in row:
        require(row['individual_probability_1e9_is_acceptance_gate'] is False, 'UNAPPROVED_INDIVIDUAL_CIF_ACCEPTANCE_GATE')

def draw_proof(row, group=None):
    if group is not None:
        return (group['draw_hash'], group['draw_multiplicity_hash'], group['parent_frame_hash'], group['parent_cluster_membership_hash'])
    draw = row['_draw']
    return (draw['draw_hash'], draw['multiplicity_hash'], draw['parent_frame_hash'], draw['parent_cluster_membership_hash'])

def pair_proof(a, b):
    require(a['indicator_id'] == b['indicator_id'] and a['bootstrap_id'] == b['bootstrap_id'], 'PAIR_NOT_SAME_INDICATOR_AND_ACTUAL_ID')
    require(a['_actual_pair_proof'] == b['_actual_pair_proof'], 'ACTUAL_DRAW_MULTIPLICITY_PARENT_PAIR_CONFLICT')
    require(all((a[k] == b[k] for k in ('target_hash', 'sample_hash', 'selection_weight_hash'))), 'PROTEIN_COMMON_SAMPLE_OR_WEIGHT_PAIR_CONFLICT')

def finish_presence(manifest, recovery_root):
    workers = [recovery_root / ('worker_%02d' % worker) / 'COMPLETE.json' for worker in range(1, 21)]
    groups = [recovery_root / 'groups' / group['group_id'] / 'COMPLETE.json' for group in manifest['groups']]
    cells = [recovery_root / 'groups' / group['group_id'] / ('endpoint_' + endpoint['endpoint_id'].replace('.', '_')) / 'COMPLETE.json' for group in manifest['groups'] for endpoint in group['endpoints']]
    found = {name: sum((p.is_file() for p in paths)) for name, paths in (('workers', workers), ('groups', groups), ('cells', cells))}
    expected = dict(workers=20, groups=4000, cells=13958)
    require({'workers': len(workers), 'groups': len(groups), 'cells': len(cells)} == expected, 'FROZEN_RECOVERY_MANIFEST_SCOPE_CHANGED')
    return dict(complete=found == expected, found=found, expected=expected, missing_examples=[str(p) for p in workers + groups + cells if not p.is_file()][:20])

def aggregate(manifest_path, recovery_root, source_code_root, phase1_root):
    """Run on WinPC only; return small aggregate tables, never individual data."""
    import sqlite3
    from concurrent.futures import ThreadPoolExecutor
    manifest_path, recovery_root = (Path(manifest_path), Path(recovery_root))
    manifest, manifestsha = jsread(manifest_path)
    require(manifest['schema'] == 'FOUR_PROTEIN_WARNING_PHASE2_V1' and manifest['worker_count'] == 20, 'FROZEN_RECOVERY_MANIFEST_REQUIRED')
    presence = finish_presence(manifest, recovery_root)
    if not presence['complete']:
        return dict(schema='FOUR_PROTEIN_RECOVERY_AGGREGATE_V1', complete=False, status='PARTIAL_RECOVERY_NO_FINAL_INFERENCE_OR_RELEASE', presence=presence)
    source_code_root = Path(source_code_root)
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(source_code_root / 'tools'))
    os.environ['FF_HOST'] = 'winpc'
    import summary_add9 as kernels
    from common import require_v4_row, v4_decision
    decision = v4_decision(verify_files=False)
    frozen = {k: manifest[k] for k in ('method_hash', 'source_hash', 'input_hash')}
    require(frozen == dict(method_hash=decision['new_method_hash'], source_hash=decision['new_source_hash'], input_hash=decision['input_lock_sha256']), 'ORIGINAL_FROZEN_DECISION_CHANGED')
    expected_base = dict(manifest_sha256=manifestsha, recovery_policy_sha256=manifest['recovery_policy_sha256'], **frozen)
    group_by_key = {(g['indicator_id'], int(g['bootstrap_id'])): g for g in manifest['groups']}
    worker_results = []
    for worker in range(1, 21):
        result, _, _, _ = sealed_result(recovery_root / ('worker_%02d' % worker), dict(expected_base, worker_id=worker))
        require(result.get('new_draws') == 0 and result.get('new_A_models') == 0 and (result.get('new_selection_models') == 0) and (result.get('old_VALID_overwritten') is False), 'WORKER_REUSE_OR_SAFETY_BOUNDARY_VIOLATED')
        worker_results.append(result)
    runtime_maps = [r.get('recovery_runtime_sha256.1') for r in worker_results]
    runtime_digests = {r.get('recovery_runtime_sha256') for r in worker_results}
    require(all((isinstance(m, dict) and set(m) == {'phase2_worker.R', 'phase2_policy.R', 'phase1_helpers.R'} for m in runtime_maps)), 'EXACT_SEALED_RUNTIME_PER_FILE_MAP_REQUIRED')
    runtime_signatures = {json.dumps(m, sort_keys=True) for m in runtime_maps}
    frozen_signatures = {json.dumps(obj(r.get('frozen_source_sha256')), sort_keys=True) for r in worker_results}
    require(len(runtime_signatures) == 1 and len(frozen_signatures) == 1 and (len(runtime_digests) == 1) and all((isinstance(d, str) and len(d) == 64 for d in runtime_digests)) and all((m['phase2_policy.R'] == manifest['recovery_policy_sha256'] for m in runtime_maps)), 'CROSS_HOST_OR_WORKER_SCIENTIFIC_SOURCE_CHANGED')
    db = sqlite3.connect('file:' + manifest['database_path'] + '?mode=ro', uri=True, timeout=30)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA query_only=ON')
    db.execute('BEGIN')
    source_counts = {str(row[0]): row[1] for row in db.execute('SELECT state,count(*) FROM cells WHERE bid>0 GROUP BY state')}
    require(source_counts == {'VALID': 30917, 'FAILED_STAT': 6083, 'BLOCKED_POINT_NOT_RUN': 8000}, 'ORIGINAL_ADD9_COUNTS_CHANGED')
    sqlids = ','.join(('?' for _ in IDS))
    records, old_rows, original_failures = ({}, {}, {})
    digest, preserved = (hashlib.sha256(), 0)
    for cell in db.execute('SELECT indicator,bid,endpoint,state,row_json FROM cells WHERE indicator IN (' + sqlids + ') ORDER BY indicator,bid,endpoint', IDS):
        key = (cell['indicator'], int(cell['bid']), str(cell['endpoint']))
        row = json.loads(cell['row_json'])
        identity(row, dict(indicator_id=key[0], bootstrap_id=key[1], endpoint_id=key[2], **frozen))
        require(row['status'] == cell['state'] and key not in old_rows, 'ORIGINAL_STATE_OR_DUPLICATE_KEY_CHANGED')
        require_v4_row(row, decision)
        old_rows[key] = row
        if row['status'] == 'VALID':
            if key[1] > 0:
                preserved += 1
                digest.update(json.dumps([*key, cell['row_json']], ensure_ascii=False, separators=(',', ':')).encode())
                digest.update(b'\n')
            effects_check(row)
            row['_actual_pair_proof'] = draw_proof(row)
            if key[1] > 0:
                group = group_by_key[key[:2]]
                require(row['_actual_pair_proof'] == draw_proof(row, group), 'ORIGINAL_VALID_ACTUAL_DRAW_DIFFERS_FROM_GROUP')
                identity(row, {k: group[k] for k in ('target_hash', 'sample_hash', 'selection_weight_hash')})
            row['_review_origin'] = 'IMMUTABLE_CURRENT_V4_ORIGINAL_VALID'
            records[key] = row
        else:
            original_failures[key] = row
    db.rollback()
    db.close()
    require(len(old_rows) == 20020 and preserved == 5917 and (digest.hexdigest() == manifest['preserved_valid_row_json_sha256']), 'IMMUTABLE_5917_VALID_ROWS_DIGEST_CHANGED')
    singular_keys = {(r['indicator_id'], int(r['bootstrap_id']), str(r['endpoint_id'])) for r in manifest['known_singular_cases']}
    require(len(singular_keys) == 125, 'KNOWN_SINGULAR_LEDGER_SCOPE_CHANGED')
    failures = [dict(indicator_id=k[0], bootstrap_id=k[1], endpoint_id=k[2], status=old_rows[k]['status'], reason=old_rows[k]['reason'], origin='KNOWN_ORIGINAL_SINGULAR_RETAINED_NO_REFIT') for k in sorted(singular_keys)]
    adjudication, adjudicationsha = jsread(manifest['phase1_adjudication_path'], manifest['phase1_adjudication_sha256'])
    require(adjudication['schema'] == 'FOUR_PROTEIN_PHASE1_POSTPROCESS_ADJUDICATION_V1', 'PHASE1_ADJUDICATION_CHANGED')
    accepted = {entry['case_id']: entry for entry in adjudication['entries']}
    b0_accepted = 0
    for protein in IDS:
        for endpoint in EPS[:2]:
            case_id = '%s_B0000_E%s' % (protein, endpoint.replace('.', '_'))
            entry = accepted[case_id]
            require(entry['point_recovery_accepted'] is True and entry['postrun_task_A_draw_model_source_hashes_match'] is True, 'B0_RECOVERY_NOT_ADJUDICATED')
            require(str(entry['result_path']).startswith(str(phase1_root) + '/'), 'B0_RESULT_OUTSIDE_ADJUDICATED_PHASE1_ROOT')
            result, _ = jsread(entry['result_path'], entry['result_sha256'])
            seal, _ = jsread(Path(entry['result_path']).parent / 'COMPLETE.json', entry['complete_sha256'])
            require(seal['complete'] and seal['result_sha256'] == entry['result_sha256'], 'PHASE1_B0_RESULT_SEAL_CHANGED')
            identity(result, dict(indicator_id=protein, bootstrap_id=0, endpoint_id=endpoint, **frozen))
            probability_check(result)
            key = (protein, 0, endpoint)
            original = old_rows[key]
            require(original['status'] == 'BLOCKED_POINT', 'B0_NOT_ORIGINAL_BLOCKED_POINT')
            require(all((result[k] == original[k] for k in ('target_hash', 'sample_hash', 'selection_weight_hash'))), 'B0_SAMPLE_OR_A_WEIGHT_CHANGED')
            require(draw_proof(original)[0] == result['actual_draw_hash'] and draw_proof(original)[1] == result['draw_multiplicity_hash'] and (original['_draw']['actual_draw_private_sha256'] == result['draw_private_sha256']), 'B0_ACTUAL_DRAW_OR_FREQUENCY_CHANGED')
            row = dict(original)
            row.update(result['original_effects'])
            row.update(status='WARNING_QUALIFIED_POINT_CONDITIONAL_RECOVERY', recovery_policy_sha256=manifest['recovery_policy_sha256'], _review_origin='ADJUDICATED_PHASE1_B0_ORIGINAL_NATIVE_FAMILY', _actual_pair_proof=draw_proof(original))
            effects_check(row)
            records[key] = row
            b0_accepted += 1

    def reduce_group(group):
        gid, protein, bid = (group['group_id'], group['indicator_id'], int(group['bootstrap_id']))
        expected = dict(expected_base, worker_id=int(group['assigned_worker']), group_id=gid, indicator_id=protein, bootstrap_id=bid)
        directory = recovery_root / 'groups' / gid
        gr, gs, _, _ = sealed_result(directory, expected)
        identity(gr, dict(A_fits=0, selection_fits=0, draws_generated=0, production_modified=False))
        entries = obj(gr['cases'])
        require(set(entries) == {e['endpoint_id'] for e in group['endpoints']}, 'RECOVERY_GROUP_ENDPOINT_SCOPE_CHANGED')
        valid, failed, stats = ([], [], collections.Counter())
        for item in group['endpoints']:
            endpoint = item['endpoint_id']
            key = (protein, bid, endpoint)
            require(key not in singular_keys and old_rows[key]['status'] != 'VALID', 'SELECTIVE_RECOVERY_OVERLAPS_PROTECTED_CELL')
            entry = entries[endpoint]
            relative = 'endpoint_' + endpoint.replace('.', '_')
            require(entry['relative_directory'] == relative and gs['artifacts_sha256'][relative + '/COMPLETE.json'] == entry['complete_sha256'], 'GROUP_CELL_SEAL_INDEX_CONFLICT')
            row, seal, _, _ = sealed_result(directory / relative, dict(expected, endpoint_id=endpoint), entry['complete_sha256'])
            require(row['status'] == seal['status'] == entry['status'], 'RECOVERY_CELL_STATUS_SEAL_CONFLICT')
            identity(row, dict(original_status=item['original_state'], original_reason=item['reason'], A_fits=0, selection_fits=0, draws_generated=0, original_VALID_outputs_modified=False, production_queue_touched=False))
            require_v4_row(row, decision)
            stats['sealed_cells'] += 1
            stats['Y_cause_fits_attempted'] += int(row.get('Y_cause_fits_attempted', 0))
            stats['source_' + item['original_state']] += 1
            stats['source_' + item['original_state'] + '_Y_fit_attempted'] += int(row.get('Y_cause_fits_attempted', 0) > 0)
            stats['phase1_reused'] += int(row.get('phase1_reused_without_Y_refit') is True)
            if row['status'] in GOOD:
                identity(row, {k: group[k] for k in ('target_hash', 'sample_hash', 'selection_weight_hash')})
                identity(row, dict(actual_draw_hash=group['draw_hash'], draw_private_sha256=group['draw_rds_sha256'], A_descriptor_sha256=group['a_descriptors_sha256']))
                effects_check(row)
                probability_check(row)
                require(all((same(row[k], row['original_effects'][k]) for k in ('R00', 'R10', 'R11', 'TE', 'DE', 'IE'))), 'TOP_LEVEL_EFFECTS_NOT_ORIGINAL_FAMILY')
                row['_actual_pair_proof'] = draw_proof(row, group)
                row['_review_origin'] = 'SEPARATE_SELECTIVE_WARNING_RECOVERY_POLICY'
                row['_original_state'] = item['original_state']
                valid.append((key, row))
                stats['accepted_' + item['original_state']] += 1
                stats['accepted_status_' + row['status']] += 1
                stats['individual_1e_minus9_diagnostic_false'] += int(row.get('individual_probability_1e9_diagnostic_passed') is False)
            else:
                require(row['status'] in ('FAILED_STAT', 'FAILED_TECH', 'BLOCKED_POINT'), 'UNREVIEWED_RECOVERY_TERMINAL_STATE')
                failed.append(dict(indicator_id=protein, bootstrap_id=bid, endpoint_id=endpoint, status=row['status'], reason=row.get('reason', ''), origin='NEW_RECOVERY_FAILURE_NO_FALLBACK', original_state=item['original_state'], group_id=gid, result_path=str(directory / relative / 'RESULT.json')))
                stats['failed_' + item['original_state']] += 1
        return (valid, failed, stats)
    totals = collections.Counter()
    with ThreadPoolExecutor(max_workers=8) as pool:
        for valid, failed, stats in pool.map(reduce_group, manifest['groups']):
            for key, row in valid:
                require(key not in records, 'RECOVERY_OVERWRITES_OR_DUPLICATES_VALID_CELL')
                records[key] = row
            failures.extend(failed)
            totals.update(stats)
    require(totals['sealed_cells'] == 13958 and totals['source_FAILED_STAT'] == 5958 and (totals['source_BLOCKED_POINT_NOT_RUN'] == 8000), 'FINAL_RECOVERY_CELL_COUNTS_CONFLICT')
    require(totals['phase1_reused'] == 3 and b0_accepted == 8, 'PHASE1_REUSE_COUNTS_CHANGED')
    require(len(records) + len(failures) == 20020, 'FINAL_PROTEIN_SCOPE_NOT_EXHAUSTIVE')
    parents = {record['_actual_pair_proof'][2:] for record in records.values()}
    require(len(parents) == 1, 'FROZEN_PARENT_FRAME_OR_CLUSTER_CHANGED')
    singles, contrasts = ([], [])
    failure_by_series = collections.Counter(((f['indicator_id'], f['endpoint_id'], f['status']) for f in failures))
    for protein in IDS:
        for endpoint in EPS:
            point = records.get((protein, 0, endpoint))
            draws = [records[protein, bid, endpoint] for bid in range(1, 1001) if (protein, bid, endpoint) in records]
            row = dict(indicator_id=protein, indicator_name=protein.split('__')[1].upper(), endpoint_id=endpoint, disease=NAMES[endpoint], B0_status=point['status'] if point else 'NOT_ESTIMABLE', bootstrap_valid=len(draws), bootstrap_planned=1000, failed_stat=failure_by_series[protein, endpoint, 'FAILED_STAT'], failed_tech=failure_by_series[protein, endpoint, 'FAILED_TECH'], blocked_not_run=0, q=None, result_status='COMPLETE_NUMBERED_BOOTSTRAP_CONDITIONAL_WARNING_POLICY' if len(draws) == 1000 else 'PARTIAL_NUMBERED_BOOTSTRAP_CONDITIONAL_WARNING_POLICY', origin='ORIGINAL_VALID_PLUS_SELECTIVE_WARNING_RECOVERY_NO_REPLACEMENT', source_run=str(recovery_root), horizon_years=11, full_followup=True, finite_joint_MLE_certified=False, **frozen)
            require(point is not None, 'ADJUDICATED_B0_MISSING')
            row.update(N_target=point.get('N_target'), N_complete=point.get('N_complete'), analysis_n=point.get('N_complete'), original_valid_replicates=sum((d['_review_origin'] == 'IMMUTABLE_CURRENT_V4_ORIGINAL_VALID' for d in draws)), recovered_replicates=sum((d['_review_origin'] != 'IMMUTABLE_CURRENT_V4_ORIGINAL_VALID' for d in draws)), **kernels.estimates(point, draws))
            row['inference_caveat'] = 'SEPARATELY_AUTHORIZED_NUISANCE_WARNING_CONDITIONAL_RECOVERY;KNOWN_AND_NEW_FAILURES_RETAINED;NOT_RANDOM_MISSING;NO_FINITE_JOINT_MLE_CERTIFICATE;NO_NEW_BH'
            singles.append(row)
        for a, b in PAIRS:
            pa, pb = (records[protein, 0, a], records[protein, 0, b])
            pair_proof(pa, pb)
            matched = []
            for bid in range(1, 1001):
                ka, kb = ((protein, bid, a), (protein, bid, b))
                if ka in records and kb in records:
                    pair_proof(records[ka], records[kb])
                    matched.append({k: records[ka][k] - records[kb][k] for k in ('TE', 'DE', 'IE')})
            row = dict(indicator_id=protein, indicator_name=protein.split('__')[1].upper(), contrast=a + '_minus_' + b, disease_A=NAMES[a], disease_B=NAMES[b], paired_valid=len(matched), paired_missing=1000 - len(matched), bootstrap_planned=1000, actual_draw_frequency_parent_proof_verified=True, same_actual_complete_sample_and_weights_verified=True, q=None, result_status='COMPLETE_PAIRED_CONDITIONAL_WARNING_POLICY' if len(matched) == 1000 else 'PARTIAL_PAIRED_CONDITIONAL_WARNING_POLICY', origin='ORIGINAL_VALID_PLUS_SELECTIVE_WARNING_RECOVERY_NO_REPLACEMENT', source_run=str(recovery_root), IE_A_per_1000=1000 * pa['IE'], IE_B_per_1000=1000 * pb['IE'], horizon_years=11, test_family='ADD9_DESCRIPTIVE_UNADJUSTED_NO_NEW_BH_FAMILY', **frozen)
            for effect in ('TE', 'DE', 'IE'):
                row['Delta_' + effect] = pa[effect] - pb[effect]
                row['Delta_' + effect + '_per_1000'] = 1000 * row['Delta_' + effect]
                for side, probability in (('low', 0.025), ('high', 0.975)):
                    value = kernels.q7([d[effect] for d in matched], probability) if len(matched) >= 2 else None
                    row['Delta_' + effect + '_CI_' + side] = value
                    row['Delta_' + effect + '_CI_' + side + '_per_1000'] = 1000 * value if value is not None else None
            if len(matched) >= 2:
                p = kernels.invert([d['IE'] for d in matched], row['Delta_IE'])
                row.update(P_Delta_IE=p['P_raw'], P_method=p['P_method'], P_resolution=p['P_grid_step'])
            row['inference_caveat'] = 'MATCHED_ACTUAL_ID_NOT_SEED_ONLY;CONDITIONAL_WARNING_POLICY;ACTUAL_VALID_DENOMINATOR;NO_NEW_BH'
            contrasts.append(row)
    require(sha(manifest_path) == manifestsha, 'MANIFEST_CHANGED_DURING_AGGREGATION')
    summary = dict(schema='FOUR_PROTEIN_RECOVERY_FINAL_SUMMARY_V1', complete=True, presence=presence, original_source_bootstrap_counts=source_counts, original_valid_protein_cells_preserved=5917, preserved_valid_row_json_sha256=digest.hexdigest(), digest_matches_preparation=True, original_failed_coronary_cells=6083, selected_nonsingular_failed_coronary_cells=5958, known_singular_cells_not_refitted=125, original_DM_T2D_blocked_unrun_cells=8000, formerly_failed_coronary_cells_recovered=totals['accepted_FAILED_STAT'], formerly_failed_coronary_cells_still_unresolved=125 + totals['failed_FAILED_STAT'], formerly_blocked_DM_T2D_cells_recovered=totals['accepted_BLOCKED_POINT_NOT_RUN'], formerly_blocked_DM_T2D_cells_newly_failed=totals['failed_BLOCKED_POINT_NOT_RUN'], formerly_blocked_DM_T2D_cells_with_Y_fit_attempted=totals['source_BLOCKED_POINT_NOT_RUN_Y_fit_attempted'], B0_newly_adjudicated=8, B0_original_coronary_valid=12, totals=dict(totals), final_valid_protein_bootstrap_cells=len(records) - 20, final_failure_cells=len(failures), paired_coverages=[dict(indicator_id=r['indicator_id'], contrast=r['contrast'], paired_valid=r['paired_valid']) for r in contrasts], final_all_series_have_1000_valid_replicates=all((r['bootstrap_valid'] == 1000 for r in singles)), finite_joint_MLE_certified=False, additional_BH_computed=False, q_for_recovered_rows=None, models_fit_by_summarizer=0, new_draws=0, production_db_writes=0, all_13958_result_hashes_verified=True, all_4000_group_seals_and_20_worker_seals_verified=True, private_RDS_or_model_blobs_read_or_transferred=False, private_blob_hash_audit='EXPECTED_HASHES_RETAINED_IN_WORKER_SEALS_NOT_BULK_REHASHED_DURING_AGGREGATION', manifest_path=str(manifest_path), manifest_sha256=manifestsha, recovery_policy_sha256=manifest['recovery_policy_sha256'], phase1_adjudication_sha256=adjudicationsha, recovery_runtime_source_sha256=runtime_maps[0], recovery_runtime_vector_digest_sha256=next(iter(runtime_digests)), worker_runtime_map_field='recovery_runtime_sha256.1', frozen_scientific_source_sha256=obj(worker_results[0].get('frozen_source_sha256')), inference_kernels_source_sha256=sha(source_code_root / 'tools/summary_add9.py'), quantile_Fieller_kernel_source_sha256=sha(source_code_root / 'inference/methods.py'), **frozen)
    return dict(schema='FOUR_PROTEIN_RECOVERY_AGGREGATE_V1', complete=True, summary=summary, protein_single20=singles, protein_six_contrasts24=contrasts, failures=sorted(failures, key=lambda f: (f['indicator_id'], f['bootstrap_id'], f['endpoint_id'])))

def csvread(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def csvwrite(path, data):
    fields = list(dict.fromkeys((k for row in data for k in row)))
    require(bool(fields), 'EMPTY_TABLE_REQUIRES_EXPLICIT_SCHEMA')
    with Path(path).open('x', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(data)

def display_tables(aggregate_result, review):
    old_primary = csvread(review / 'DISPLAY18_PRIMARY_REVIEW.csv')
    old_single = csvread(review / 'DISPLAY18_SINGLE90_REVIEW.csv')
    replacements = {r['indicator_id']: r for r in aggregate_result['protein_six_contrasts24'] if r['contrast'] == '411_minus_250.2'}
    single_replacements = {(r['indicator_id'], r['endpoint_id']): r for r in aggregate_result['protein_single20']}
    primary, single = ([], [])
    for old in old_primary:
        if old['indicator_id'] in IDS:
            new = dict(replacements[old['indicator_id']])
            new.update(P=new.get('P_Delta_IE'), q=None, q_family='NO_NEW_BH_FAMILY_AUTHORIZED', current_count_point_pair_hash_verified=True, display_review_source='NEW_SELECTIVE_WARNING_RECOVERY')
        else:
            new = dict(old, display_review_source='UNCHANGED_ACCEPTED_OLD9_OR_NEW5')
            require(all((new[k] == v for k, v in old.items())), 'NONPROTEIN_PRIMARY_ROW_CHANGED')
        primary.append(new)
    for old in old_single:
        if old['indicator_id'] in IDS:
            new = dict(single_replacements[old['indicator_id'], old['endpoint_id']], display_review_source='NEW_SELECTIVE_WARNING_RECOVERY')
        else:
            new = dict(old, display_review_source='UNCHANGED_ACCEPTED_OLD9_OR_NEW5')
            require(all((new[k] == v for k, v in old.items())), 'NONPROTEIN_SINGLE_ROW_CHANGED')
        single.append(new)
    require(len(primary) == 18 and len(single) == 90 and (len(replacements) == 4) and (len(single_replacements) == 20) and (len({r['indicator_id'] for r in primary}) == 18) and (len({(r['indicator_id'], r['endpoint_id']) for r in single}) == 90), 'EXACT_UNIFIED_DISPLAY18_SCOPE_REQUIRED')
    return (primary, single)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    remote = sub.add_parser('remote-main')
    remote.add_argument('--manifest', type=Path, required=True)
    remote.add_argument('--recovery-root', type=Path, required=True)
    remote.add_argument('--source-code-root', type=Path, default=Path(birth_private_setting('PRIVATE_BORDER_FULLFOLLOWUP_B1000_ADD9_20261030_V1')))
    remote.add_argument('--phase1-root', type=Path, default=Path(birth_private_setting('PRIVATE_RUN_20261031_PHASE1_V1')))
    local = sub.add_parser('collect')
    local.add_argument('--remote-manifest', required=True)
    local.add_argument('--recovery-root', required=True)
    local.add_argument('--out', type=Path, required=True)
    local.add_argument('--review-dir', type=Path, default=HERE.parents[2] / REVIEW_REL)
    local.add_argument('--remote-code-root', default=birth_private_setting('PRIVATE_BORDER_FOURPROTEIN_WARNING_RECOVERY_20261031_V1'))
    args = parser.parse_args()
    if args.action == 'remote-main':
        value = aggregate(args.manifest, args.recovery_root, args.source_code_root, args.phase1_root)
        print(json.dumps(value, ensure_ascii=False, separators=(',', ':'), allow_nan=False))
        return
    require(not args.out.exists(), 'PRESERVE_EXISTING_LOCAL_RELEASE')
    from remote_control import remote as call_remote
    code = 'import sys,json\nsys.dont_write_bytecode=True\nsys.path.insert(0,' + repr(args.remote_code_root) + ')\nimport summarize_recovery_schema_v2 as s\nz=s.aggregate(' + repr(args.remote_manifest) + ',' + repr(args.recovery_root) + ",'/mnt/d/scripts/border/border_fullfollowup_b1000_add9_20261030_v1','/mnt/g/UKB_results/birth/border_fourprotein_warning_recovery_20261031_v1/run_20261031_phase1_v1')\nprint(json.dumps(z,ensure_ascii=False,separators=(',',':'),allow_nan=False))"
    value = call_remote(code, timeout=1800)
    require(value.get('complete') is True, 'RECOVERY_NOT_COMPLETE_NO_RELEASE:' + json.dumps(value, ensure_ascii=False))
    primary, single = display_tables(value, args.review_dir)
    args.out.mkdir(parents=True, exist_ok=False)
    tables = {'PROTEIN_SINGLE20_REVIEW.csv': value['protein_single20'], 'PROTEIN_SIX_CONTRASTS24_REVIEW.csv': value['protein_six_contrasts24'], 'PROTEIN_IHD_T2D_REVIEW.csv': [r for r in value['protein_six_contrasts24'] if r['contrast'] == '411_minus_250.2'], 'ALL_PROTEIN_FAILURE_RECORDS.csv': value['failures'], 'DISPLAY18_PRIMARY_REVIEW.csv': primary, 'DISPLAY18_SINGLE90_REVIEW.csv': single}
    for name, data in tables.items():
        csvwrite(args.out / name, data)
    summary = value['summary']
    summary['unchanged_nonprotein_display_indicators'] = 14
    summary['original_display_sources'] = [dict(path=str((args.review_dir / name).resolve()), sha256=sha(args.review_dir / name)) for name in ('DISPLAY18_PRIMARY_REVIEW.csv', 'DISPLAY18_SINGLE90_REVIEW.csv')]
    with (args.out / 'RECOVERY_SUMMARY.json').open('x', encoding='utf-8') as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    release = dict(schema='FOUR_PROTEIN_RECOVERY_LOCAL_RELEASE_MANIFEST_V1', complete=True, files_sha256={name: sha(args.out / name) for name in [*tables, 'RECOVERY_SUMMARY.json']}, source_code_sha256=sha(Path(__file__)), production_modified=False, models_fit=0, new_BH=False)
    with (args.out / 'RELEASE_MANIFEST.json').open('x', encoding='utf-8') as stream:
        json.dump(release, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.out.resolve()), table_rows={name: len(data) for name, data in tables.items()}, final_valid_protein_bootstrap_cells=summary['final_valid_protein_bootstrap_cells'], final_failure_cells=summary['final_failure_cells'])))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Seal-checked aggregate-only recovery reduction; no fits, draws or DB writes.\n\nDeploy this code-only file through the reviewed transport, then use `collect`\non the Mac. Remote reduction reads aggregate JSON/SQLite rows, never private\nRDS/model/participant arrays. Existing production and releases remain intact.\n"""\nfrom __future__ import annotations\nimport argparse\nimport collections\nimport csv\nimport hashlib\nimport io\nimport json\nimport math\nimport os\nimport sys\nfrom pathlib import Path\nIDS = (\'olink__fcrl6\', \'olink__gzma\', \'olink__gzmh\', \'olink__klrd1\')\nEPS = (\'250\', \'250.2\', \'411\', \'411.3\', \'411.4\')\nNAMES = dict(zip(EPS, (\'DM\', \'T2D\', \'IHD\', \'Angina\', \'CAS\')))\nPAIRS = tuple(((a, b) for a in EPS[2:] for b in EPS[:2]))\nGOOD = (\'VALID\', \'WARNING_QUALIFIED_CONDITIONAL_ESTIMATE\')\nHERE = Path(__file__).resolve().parent\nREVIEW_REL = \'reports/birth/b1000_add9_acceptance_20261031_v1/review_20261031_v1\'\n\ndef require(condition, message):\n    if not condition:\n        raise ValueError(message)\n\ndef sha(path):\n    h = hashlib.sha256()\n    with Path(path).open(\'rb\') as stream:\n        for block in iter(lambda: stream.read(1024 * 1024), b\'\'):\n            h.update(block)\n    return h.hexdigest()\n\ndef jsread(path, expected=None):\n    path = Path(path)\n    require(path.is_file() and path.stat().st_size <= 32 * 1024 ** 2, \'SMALL_AGGREGATE_SOURCE_REQUIRED:\' + str(path))\n    content = path.read_bytes()\n    digest = hashlib.sha256(content).hexdigest()\n    require(expected is None or digest == expected, \'AGGREGATE_HASH_CHANGED:\' + str(path))\n    return (json.loads(content), digest)\n\ndef obj(value):\n    return value if isinstance(value, dict) else {}\n\ndef rows(value):\n    return value if isinstance(value, list) else list(value.values()) if isinstance(value, dict) else []\n\ndef numeric(value):\n    if value is None or isinstance(value, bool):\n        return None\n    try:\n        value = float(value)\n        return value if math.isfinite(value) else None\n    except (TypeError, ValueError):\n        return None\n\ndef same(a, b, tolerance=1e-12):\n    a, b = (numeric(a), numeric(b))\n    return a is not None and b is not None and (abs(a - b) <= tolerance + tolerance * max(abs(a), abs(b)))\n\ndef identity(row, expected):\n    require(all((row.get(k) == v for k, v in expected.items())), \'SEALED_IDENTITY_CHANGED:\' + str(expected))\n\ndef sealed_result(directory, expected, expected_seal=None):\n    directory = Path(directory)\n    seal, sealsha = jsread(directory / \'COMPLETE.json\', expected_seal)\n    require(seal.get(\'complete\') is True, \'INCOMPLETE_SEAL:\' + str(directory))\n    identity(seal, expected)\n    artifacts = seal.get(\'artifacts_sha256\', {})\n    require(\'RESULT.json\' in artifacts, \'SEALED_RESULT_NOT_INDEXED:\' + str(directory))\n    for relative in artifacts:\n        require(not Path(relative).is_absolute() and \'..\' not in Path(relative).parts, \'UNSAFE_SEALED_RELATIVE_PATH\')\n    result, resultsha = jsread(directory / \'RESULT.json\', artifacts[\'RESULT.json\'])\n    identity(result, expected)\n    return (result, seal, sealsha, resultsha)\n\ndef effects_check(row):\n    require(all((numeric(row.get(k)) is not None for k in (\'R00\', \'R10\', \'R11\', \'TE\', \'DE\', \'IE\'))), \'NONFINITE_EFFECT_OR_RISK\')\n    require(all((0 <= row[k] <= 1 for k in (\'R00\', \'R10\', \'R11\'))), \'MEAN_RISK_OUTSIDE_UNIT_INTERVAL\')\n    require(same(row[\'TE\'], row[\'R11\'] - row[\'R00\']) and same(row[\'DE\'], row[\'R10\'] - row[\'R00\']) and same(row[\'IE\'], row[\'R11\'] - row[\'R10\']) and same(row[\'TE\'], row[\'DE\'] + row[\'IE\']), \'SAVED_EFFECT_IDENTITY_FAILED\')\n    require(row[\'TE\'] == 0 and numeric(row.get(\'PM\')) is None or (row[\'TE\'] != 0 and same(row.get(\'PM\'), row[\'IE\'] / row[\'TE\'])), \'PM_NOT_SAME_MODEL_RATIO\')\n\ndef probability_check(row):\n    predictions = obj(row.get(\'predictions\'))\n    require(set(predictions) == {\'observed\', \'R00\', \'R10\', \'R11\'}, \'FOUR_PREDICTION_AUDITS_REQUIRED\')\n    for scenario, pred in predictions.items():\n        for family in (\'original\', \'warm\'):\n            audit = pred[family]\n            require(audit.get(\'legal\') is True and audit.get(\'nonfinite\') == 0 and (audit.get(\'outside_unit_interval\') == 0), \'INDIVIDUAL_PROBABILITY_LEGALITY_FAILED\')\n            require(numeric(audit.get(\'max_abs_closure\')) is not None and audit[\'max_abs_closure\'] < 1e-08, \'INDIVIDUAL_PROBABILITY_CLOSURE_FAILED\')\n            require(all((0 <= v <= 1 for v in obj(audit.get(\'state_min\')).values())) and all((0 <= v <= 1 for v in obj(audit.get(\'state_max\')).values())), \'INDIVIDUAL_PROBABILITY_RANGE_FAILED\')\n        if scenario != \'observed\':\n            selected = row.get(scenario, obj(row.get(\'original_effects\')).get(scenario))\n            require(same(selected, pred[\'original_weighted_disease_risk\']), \'ORIGINAL_FAMILY_RISK_NOT_SELECTED\')\n    require(row.get(\'native_validation\', {}).get(\'original\', {}).get(\'passed\') is True and row.get(\'native_validation\', {}).get(\'warm\', {}).get(\'passed\') is True, \'NATIVE_VALIDATION_FAILED\')\n    diagnostics = obj(row.get(\'coefficient_diagnostics\'))\n    require(set(diagnostics) == {\'disease\', \'death\'}, \'BOTH_CAUSE_DIAGNOSTICS_REQUIRED\')\n    for diagnostic in diagnostics.values():\n        stable = diagnostic.get(\'core_stability_passed\', diagnostic.get(\'core_stability_descriptive_rule_passed\'))\n        require(stable is True and diagnostic.get(\'core_positive_information_both\') is True, \'CORE_INFORMATION_OR_STABILITY_FAILED\')\n    if \'individual_probability_1e9_is_acceptance_gate\' in row:\n        require(row[\'individual_probability_1e9_is_acceptance_gate\'] is False, \'UNAPPROVED_INDIVIDUAL_CIF_ACCEPTANCE_GATE\')\n\ndef draw_proof(row, group=None):\n    if group is not None:\n        return (group[\'draw_hash\'], group[\'draw_multiplicity_hash\'], group[\'parent_frame_hash\'], group[\'parent_cluster_membership_hash\'])\n    draw = row[\'_draw\']\n    return (draw[\'draw_hash\'], draw[\'multiplicity_hash\'], draw[\'parent_frame_hash\'], draw[\'parent_cluster_membership_hash\'])\n\ndef pair_proof(a, b):\n    require(a[\'indicator_id\'] == b[\'indicator_id\'] and a[\'bootstrap_id\'] == b[\'bootstrap_id\'], \'PAIR_NOT_SAME_INDICATOR_AND_ACTUAL_ID\')\n    require(a[\'_actual_pair_proof\'] == b[\'_actual_pair_proof\'], \'ACTUAL_DRAW_MULTIPLICITY_PARENT_PAIR_CONFLICT\')\n    require(all((a[k] == b[k] for k in (\'target_hash\', \'sample_hash\', \'selection_weight_hash\'))), \'PROTEIN_COMMON_SAMPLE_OR_WEIGHT_PAIR_CONFLICT\')\n\ndef finish_presence(manifest, recovery_root):\n    workers = [recovery_root / (\'worker_%02d\' % worker) / \'COMPLETE.json\' for worker in range(1, 21)]\n    groups = [recovery_root / \'groups\' / group[\'group_id\'] / \'COMPLETE.json\' for group in manifest[\'groups\']]\n    cells = [recovery_root / \'groups\' / group[\'group_id\'] / (\'endpoint_\' + endpoint[\'endpoint_id\'].replace(\'.\', \'_\')) / \'COMPLETE.json\' for group in manifest[\'groups\'] for endpoint in group[\'endpoints\']]\n    found = {name: sum((p.is_file() for p in paths)) for name, paths in ((\'workers\', workers), (\'groups\', groups), (\'cells\', cells))}\n    expected = dict(workers=20, groups=4000, cells=13958)\n    require({\'workers\': len(workers), \'groups\': len(groups), \'cells\': len(cells)} == expected, \'FROZEN_RECOVERY_MANIFEST_SCOPE_CHANGED\')\n    return dict(complete=found == expected, found=found, expected=expected, missing_examples=[str(p) for p in workers + groups + cells if not p.is_file()][:20])\n\ndef aggregate(manifest_path, recovery_root, source_code_root, phase1_root):\n    """Run on WinPC only; return small aggregate tables, never individual data."""\n    import sqlite3\n    from concurrent.futures import ThreadPoolExecutor\n    manifest_path, recovery_root = (Path(manifest_path), Path(recovery_root))\n    manifest, manifestsha = jsread(manifest_path)\n    require(manifest[\'schema\'] == \'FOUR_PROTEIN_WARNING_PHASE2_V1\' and manifest[\'worker_count\'] == 20, \'FROZEN_RECOVERY_MANIFEST_REQUIRED\')\n    presence = finish_presence(manifest, recovery_root)\n    if not presence[\'complete\']:\n        return dict(schema=\'FOUR_PROTEIN_RECOVERY_AGGREGATE_V1\', complete=False, status=\'PARTIAL_RECOVERY_NO_FINAL_INFERENCE_OR_RELEASE\', presence=presence)\n    source_code_root = Path(source_code_root)\n    sys.dont_write_bytecode = True\n    sys.path.insert(0, str(source_code_root / \'tools\'))\n    os.environ[\'FF_HOST\'] = \'winpc\'\n    import summary_add9 as kernels\n    from common import require_v4_row, v4_decision\n    decision = v4_decision(verify_files=False)\n    frozen = {k: manifest[k] for k in (\'method_hash\', \'source_hash\', \'input_hash\')}\n    require(frozen == dict(method_hash=decision[\'new_method_hash\'], source_hash=decision[\'new_source_hash\'], input_hash=decision[\'input_lock_sha256\']), \'ORIGINAL_FROZEN_DECISION_CHANGED\')\n    expected_base = dict(manifest_sha256=manifestsha, recovery_policy_sha256=manifest[\'recovery_policy_sha256\'], **frozen)\n    group_by_key = {(g[\'indicator_id\'], int(g[\'bootstrap_id\'])): g for g in manifest[\'groups\']}\n    worker_results = []\n    for worker in range(1, 21):\n        result, _, _, _ = sealed_result(recovery_root / (\'worker_%02d\' % worker), dict(expected_base, worker_id=worker))\n        require(result.get(\'new_draws\') == 0 and result.get(\'new_A_models\') == 0 and (result.get(\'new_selection_models\') == 0) and (result.get(\'old_VALID_overwritten\') is False), \'WORKER_REUSE_OR_SAFETY_BOUNDARY_VIOLATED\')\n        worker_results.append(result)\n    runtime_maps = [r.get(\'recovery_runtime_sha256.1\') for r in worker_results]\n    runtime_digests = {r.get(\'recovery_runtime_sha256\') for r in worker_results}\n    require(all((isinstance(m, dict) and set(m) == {\'phase2_worker.R\', \'phase2_policy.R\', \'phase1_helpers.R\'} for m in runtime_maps)), \'EXACT_SEALED_RUNTIME_PER_FILE_MAP_REQUIRED\')\n    runtime_signatures = {json.dumps(m, sort_keys=True) for m in runtime_maps}\n    frozen_signatures = {json.dumps(obj(r.get(\'frozen_source_sha256\')), sort_keys=True) for r in worker_results}\n    require(len(runtime_signatures) == 1 and len(frozen_signatures) == 1 and (len(runtime_digests) == 1) and all((isinstance(d, str) and len(d) == 64 for d in runtime_digests)) and all((m[\'phase2_policy.R\'] == manifest[\'recovery_policy_sha256\'] for m in runtime_maps)), \'CROSS_HOST_OR_WORKER_SCIENTIFIC_SOURCE_CHANGED\')\n    db = sqlite3.connect(\'file:\' + manifest[\'database_path\'] + \'?mode=ro\', uri=True, timeout=30)\n    db.row_factory = sqlite3.Row\n    db.execute(\'PRAGMA query_only=ON\')\n    db.execute(\'BEGIN\')\n    source_counts = {str(row[0]): row[1] for row in db.execute(\'SELECT state,count(*) FROM cells WHERE bid>0 GROUP BY state\')}\n    require(source_counts == {\'VALID\': 30917, \'FAILED_STAT\': 6083, \'BLOCKED_POINT_NOT_RUN\': 8000}, \'ORIGINAL_ADD9_COUNTS_CHANGED\')\n    sqlids = \',\'.join((\'?\' for _ in IDS))\n    records, old_rows, original_failures = ({}, {}, {})\n    digest, preserved = (hashlib.sha256(), 0)\n    for cell in db.execute(\'SELECT indicator,bid,endpoint,state,row_json FROM cells WHERE indicator IN (\' + sqlids + \') ORDER BY indicator,bid,endpoint\', IDS):\n        key = (cell[\'indicator\'], int(cell[\'bid\']), str(cell[\'endpoint\']))\n        row = json.loads(cell[\'row_json\'])\n        identity(row, dict(indicator_id=key[0], bootstrap_id=key[1], endpoint_id=key[2], **frozen))\n        require(row[\'status\'] == cell[\'state\'] and key not in old_rows, \'ORIGINAL_STATE_OR_DUPLICATE_KEY_CHANGED\')\n        require_v4_row(row, decision)\n        old_rows[key] = row\n        if row[\'status\'] == \'VALID\':\n            if key[1] > 0:\n                preserved += 1\n                digest.update(json.dumps([*key, cell[\'row_json\']], ensure_ascii=False, separators=(\',\', \':\')).encode())\n                digest.update(b\'\\n\')\n            effects_check(row)\n            row[\'_actual_pair_proof\'] = draw_proof(row)\n            if key[1] > 0:\n                group = group_by_key[key[:2]]\n                require(row[\'_actual_pair_proof\'] == draw_proof(row, group), \'ORIGINAL_VALID_ACTUAL_DRAW_DIFFERS_FROM_GROUP\')\n                identity(row, {k: group[k] for k in (\'target_hash\', \'sample_hash\', \'selection_weight_hash\')})\n            row[\'_review_origin\'] = \'IMMUTABLE_CURRENT_V4_ORIGINAL_VALID\'\n            records[key] = row\n        else:\n            original_failures[key] = row\n    db.rollback()\n    db.close()\n    require(len(old_rows) == 20020 and preserved == 5917 and (digest.hexdigest() == manifest[\'preserved_valid_row_json_sha256\']), \'IMMUTABLE_5917_VALID_ROWS_DIGEST_CHANGED\')\n    singular_keys = {(r[\'indicator_id\'], int(r[\'bootstrap_id\']), str(r[\'endpoint_id\'])) for r in manifest[\'known_singular_cases\']}\n    require(len(singular_keys) == 125, \'KNOWN_SINGULAR_LEDGER_SCOPE_CHANGED\')\n    failures = [dict(indicator_id=k[0], bootstrap_id=k[1], endpoint_id=k[2], status=old_rows[k][\'status\'], reason=old_rows[k][\'reason\'], origin=\'KNOWN_ORIGINAL_SINGULAR_RETAINED_NO_REFIT\') for k in sorted(singular_keys)]\n    adjudication, adjudicationsha = jsread(manifest[\'phase1_adjudication_path\'], manifest[\'phase1_adjudication_sha256\'])\n    require(adjudication[\'schema\'] == \'FOUR_PROTEIN_PHASE1_POSTPROCESS_ADJUDICATION_V1\', \'PHASE1_ADJUDICATION_CHANGED\')\n    accepted = {entry[\'case_id\']: entry for entry in adjudication[\'entries\']}\n    b0_accepted = 0\n    for protein in IDS:\n        for endpoint in EPS[:2]:\n            case_id = \'%s_B0000_E%s\' % (protein, endpoint.replace(\'.\', \'_\'))\n            entry = accepted[case_id]\n            require(entry[\'point_recovery_accepted\'] is True and entry[\'postrun_task_A_draw_model_source_hashes_match\'] is True, \'B0_RECOVERY_NOT_ADJUDICATED\')\n            require(str(entry[\'result_path\']).startswith(str(phase1_root) + \'/\'), \'B0_RESULT_OUTSIDE_ADJUDICATED_PHASE1_ROOT\')\n            result, _ = jsread(entry[\'result_path\'], entry[\'result_sha256\'])\n            seal, _ = jsread(Path(entry[\'result_path\']).parent / \'COMPLETE.json\', entry[\'complete_sha256\'])\n            require(seal[\'complete\'] and seal[\'result_sha256\'] == entry[\'result_sha256\'], \'PHASE1_B0_RESULT_SEAL_CHANGED\')\n            identity(result, dict(indicator_id=protein, bootstrap_id=0, endpoint_id=endpoint, **frozen))\n            probability_check(result)\n            key = (protein, 0, endpoint)\n            original = old_rows[key]\n            require(original[\'status\'] == \'BLOCKED_POINT\', \'B0_NOT_ORIGINAL_BLOCKED_POINT\')\n            require(all((result[k] == original[k] for k in (\'target_hash\', \'sample_hash\', \'selection_weight_hash\'))), \'B0_SAMPLE_OR_A_WEIGHT_CHANGED\')\n            require(draw_proof(original)[0] == result[\'actual_draw_hash\'] and draw_proof(original)[1] == result[\'draw_multiplicity_hash\'] and (original[\'_draw\'][\'actual_draw_private_sha256\'] == result[\'draw_private_sha256\']), \'B0_ACTUAL_DRAW_OR_FREQUENCY_CHANGED\')\n            row = dict(original)\n            row.update(result[\'original_effects\'])\n            row.update(status=\'WARNING_QUALIFIED_POINT_CONDITIONAL_RECOVERY\', recovery_policy_sha256=manifest[\'recovery_policy_sha256\'], _review_origin=\'ADJUDICATED_PHASE1_B0_ORIGINAL_NATIVE_FAMILY\', _actual_pair_proof=draw_proof(original))\n            effects_check(row)\n            records[key] = row\n            b0_accepted += 1\n\n    def reduce_group(group):\n        gid, protein, bid = (group[\'group_id\'], group[\'indicator_id\'], int(group[\'bootstrap_id\']))\n        expected = dict(expected_base, worker_id=int(group[\'assigned_worker\']), group_id=gid, indicator_id=protein, bootstrap_id=bid)\n        directory = recovery_root / \'groups\' / gid\n        gr, gs, _, _ = sealed_result(directory, expected)\n        identity(gr, dict(A_fits=0, selection_fits=0, draws_generated=0, production_modified=False))\n        entries = obj(gr[\'cases\'])\n        require(set(entries) == {e[\'endpoint_id\'] for e in group[\'endpoints\']}, \'RECOVERY_GROUP_ENDPOINT_SCOPE_CHANGED\')\n        valid, failed, stats = ([], [], collections.Counter())\n        for item in group[\'endpoints\']:\n            endpoint = item[\'endpoint_id\']\n            key = (protein, bid, endpoint)\n            require(key not in singular_keys and old_rows[key][\'status\'] != \'VALID\', \'SELECTIVE_RECOVERY_OVERLAPS_PROTECTED_CELL\')\n            entry = entries[endpoint]\n            relative = \'endpoint_\' + endpoint.replace(\'.\', \'_\')\n            require(entry[\'relative_directory\'] == relative and gs[\'artifacts_sha256\'][relative + \'/COMPLETE.json\'] == entry[\'complete_sha256\'], \'GROUP_CELL_SEAL_INDEX_CONFLICT\')\n            row, seal, _, _ = sealed_result(directory / relative, dict(expected, endpoint_id=endpoint), entry[\'complete_sha256\'])\n            require(row[\'status\'] == seal[\'status\'] == entry[\'status\'], \'RECOVERY_CELL_STATUS_SEAL_CONFLICT\')\n            identity(row, dict(original_status=item[\'original_state\'], original_reason=item[\'reason\'], A_fits=0, selection_fits=0, draws_generated=0, original_VALID_outputs_modified=False, production_queue_touched=False))\n            require_v4_row(row, decision)\n            stats[\'sealed_cells\'] += 1\n            stats[\'Y_cause_fits_attempted\'] += int(row.get(\'Y_cause_fits_attempted\', 0))\n            stats[\'source_\' + item[\'original_state\']] += 1\n            stats[\'source_\' + item[\'original_state\'] + \'_Y_fit_attempted\'] += int(row.get(\'Y_cause_fits_attempted\', 0) > 0)\n            stats[\'phase1_reused\'] += int(row.get(\'phase1_reused_without_Y_refit\') is True)\n            if row[\'status\'] in GOOD:\n                identity(row, {k: group[k] for k in (\'target_hash\', \'sample_hash\', \'selection_weight_hash\')})\n                identity(row, dict(actual_draw_hash=group[\'draw_hash\'], draw_private_sha256=group[\'draw_rds_sha256\'], A_descriptor_sha256=group[\'a_descriptors_sha256\']))\n                effects_check(row)\n                probability_check(row)\n                require(all((same(row[k], row[\'original_effects\'][k]) for k in (\'R00\', \'R10\', \'R11\', \'TE\', \'DE\', \'IE\'))), \'TOP_LEVEL_EFFECTS_NOT_ORIGINAL_FAMILY\')\n                row[\'_actual_pair_proof\'] = draw_proof(row, group)\n                row[\'_review_origin\'] = \'SEPARATE_SELECTIVE_WARNING_RECOVERY_POLICY\'\n                row[\'_original_state\'] = item[\'original_state\']\n                valid.append((key, row))\n                stats[\'accepted_\' + item[\'original_state\']] += 1\n                stats[\'accepted_status_\' + row[\'status\']] += 1\n                stats[\'individual_1e_minus9_diagnostic_false\'] += int(row.get(\'individual_probability_1e9_diagnostic_passed\') is False)\n            else:\n                require(row[\'status\'] in (\'FAILED_STAT\', \'FAILED_TECH\', \'BLOCKED_POINT\'), \'UNREVIEWED_RECOVERY_TERMINAL_STATE\')\n                failed.append(dict(indicator_id=protein, bootstrap_id=bid, endpoint_id=endpoint, status=row[\'status\'], reason=row.get(\'reason\', \'\'), origin=\'NEW_RECOVERY_FAILURE_NO_FALLBACK\', original_state=item[\'original_state\'], group_id=gid, result_path=str(directory / relative / \'RESULT.json\')))\n                stats[\'failed_\' + item[\'original_state\']] += 1\n        return (valid, failed, stats)\n    totals = collections.Counter()\n    with ThreadPoolExecutor(max_workers=8) as pool:\n        for valid, failed, stats in pool.map(reduce_group, manifest[\'groups\']):\n            for key, row in valid:\n                require(key not in records, \'RECOVERY_OVERWRITES_OR_DUPLICATES_VALID_CELL\')\n                records[key] = row\n            failures.extend(failed)\n            totals.update(stats)\n    require(totals[\'sealed_cells\'] == 13958 and totals[\'source_FAILED_STAT\'] == 5958 and (totals[\'source_BLOCKED_POINT_NOT_RUN\'] == 8000), \'FINAL_RECOVERY_CELL_COUNTS_CONFLICT\')\n    require(totals[\'phase1_reused\'] == 3 and b0_accepted == 8, \'PHASE1_REUSE_COUNTS_CHANGED\')\n    require(len(records) + len(failures) == 20020, \'FINAL_PROTEIN_SCOPE_NOT_EXHAUSTIVE\')\n    parents = {record[\'_actual_pair_proof\'][2:] for record in records.values()}\n    require(len(parents) == 1, \'FROZEN_PARENT_FRAME_OR_CLUSTER_CHANGED\')\n    singles, contrasts = ([], [])\n    failure_by_series = collections.Counter(((f[\'indicator_id\'], f[\'endpoint_id\'], f[\'status\']) for f in failures))\n    for protein in IDS:\n        for endpoint in EPS:\n            point = records.get((protein, 0, endpoint))\n            draws = [records[protein, bid, endpoint] for bid in range(1, 1001) if (protein, bid, endpoint) in records]\n            row = dict(indicator_id=protein, indicator_name=protein.split(\'__\')[1].upper(), endpoint_id=endpoint, disease=NAMES[endpoint], B0_status=point[\'status\'] if point else \'NOT_ESTIMABLE\', bootstrap_valid=len(draws), bootstrap_planned=1000, failed_stat=failure_by_series[protein, endpoint, \'FAILED_STAT\'], failed_tech=failure_by_series[protein, endpoint, \'FAILED_TECH\'], blocked_not_run=0, q=None, result_status=\'COMPLETE_NUMBERED_BOOTSTRAP_CONDITIONAL_WARNING_POLICY\' if len(draws) == 1000 else \'PARTIAL_NUMBERED_BOOTSTRAP_CONDITIONAL_WARNING_POLICY\', origin=\'ORIGINAL_VALID_PLUS_SELECTIVE_WARNING_RECOVERY_NO_REPLACEMENT\', source_run=str(recovery_root), horizon_years=11, full_followup=True, finite_joint_MLE_certified=False, **frozen)\n            require(point is not None, \'ADJUDICATED_B0_MISSING\')\n            row.update(N_target=point.get(\'N_target\'), N_complete=point.get(\'N_complete\'), analysis_n=point.get(\'N_complete\'), original_valid_replicates=sum((d[\'_review_origin\'] == \'IMMUTABLE_CURRENT_V4_ORIGINAL_VALID\' for d in draws)), recovered_replicates=sum((d[\'_review_origin\'] != \'IMMUTABLE_CURRENT_V4_ORIGINAL_VALID\' for d in draws)), **kernels.estimates(point, draws))\n            row[\'inference_caveat\'] = \'SEPARATELY_AUTHORIZED_NUISANCE_WARNING_CONDITIONAL_RECOVERY;KNOWN_AND_NEW_FAILURES_RETAINED;NOT_RANDOM_MISSING;NO_FINITE_JOINT_MLE_CERTIFICATE;NO_NEW_BH\'\n            singles.append(row)\n        for a, b in PAIRS:\n            pa, pb = (records[protein, 0, a], records[protein, 0, b])\n            pair_proof(pa, pb)\n            matched = []\n            for bid in range(1, 1001):\n                ka, kb = ((protein, bid, a), (protein, bid, b))\n                if ka in records and kb in records:\n                    pair_proof(records[ka], records[kb])\n                    matched.append({k: records[ka][k] - records[kb][k] for k in (\'TE\', \'DE\', \'IE\')})\n            row = dict(indicator_id=protein, indicator_name=protein.split(\'__\')[1].upper(), contrast=a + \'_minus_\' + b, disease_A=NAMES[a], disease_B=NAMES[b], paired_valid=len(matched), paired_missing=1000 - len(matched), bootstrap_planned=1000, actual_draw_frequency_parent_proof_verified=True, same_actual_complete_sample_and_weights_verified=True, q=None, result_status=\'COMPLETE_PAIRED_CONDITIONAL_WARNING_POLICY\' if len(matched) == 1000 else \'PARTIAL_PAIRED_CONDITIONAL_WARNING_POLICY\', origin=\'ORIGINAL_VALID_PLUS_SELECTIVE_WARNING_RECOVERY_NO_REPLACEMENT\', source_run=str(recovery_root), IE_A_per_1000=1000 * pa[\'IE\'], IE_B_per_1000=1000 * pb[\'IE\'], horizon_years=11, test_family=\'ADD9_DESCRIPTIVE_UNADJUSTED_NO_NEW_BH_FAMILY\', **frozen)\n            for effect in (\'TE\', \'DE\', \'IE\'):\n                row[\'Delta_\' + effect] = pa[effect] - pb[effect]\n                row[\'Delta_\' + effect + \'_per_1000\'] = 1000 * row[\'Delta_\' + effect]\n                for side, probability in ((\'low\', 0.025), (\'high\', 0.975)):\n                    value = kernels.q7([d[effect] for d in matched], probability) if len(matched) >= 2 else None\n                    row[\'Delta_\' + effect + \'_CI_\' + side] = value\n                    row[\'Delta_\' + effect + \'_CI_\' + side + \'_per_1000\'] = 1000 * value if value is not None else None\n            if len(matched) >= 2:\n                p = kernels.invert([d[\'IE\'] for d in matched], row[\'Delta_IE\'])\n                row.update(P_Delta_IE=p[\'P_raw\'], P_method=p[\'P_method\'], P_resolution=p[\'P_grid_step\'])\n            row[\'inference_caveat\'] = \'MATCHED_ACTUAL_ID_NOT_SEED_ONLY;CONDITIONAL_WARNING_POLICY;ACTUAL_VALID_DENOMINATOR;NO_NEW_BH\'\n            contrasts.append(row)\n    require(sha(manifest_path) == manifestsha, \'MANIFEST_CHANGED_DURING_AGGREGATION\')\n    summary = dict(schema=\'FOUR_PROTEIN_RECOVERY_FINAL_SUMMARY_V1\', complete=True, presence=presence, original_source_bootstrap_counts=source_counts, original_valid_protein_cells_preserved=5917, preserved_valid_row_json_sha256=digest.hexdigest(), digest_matches_preparation=True, original_failed_coronary_cells=6083, selected_nonsingular_failed_coronary_cells=5958, known_singular_cells_not_refitted=125, original_DM_T2D_blocked_unrun_cells=8000, formerly_failed_coronary_cells_recovered=totals[\'accepted_FAILED_STAT\'], formerly_failed_coronary_cells_still_unresolved=125 + totals[\'failed_FAILED_STAT\'], formerly_blocked_DM_T2D_cells_recovered=totals[\'accepted_BLOCKED_POINT_NOT_RUN\'], formerly_blocked_DM_T2D_cells_newly_failed=totals[\'failed_BLOCKED_POINT_NOT_RUN\'], formerly_blocked_DM_T2D_cells_with_Y_fit_attempted=totals[\'source_BLOCKED_POINT_NOT_RUN_Y_fit_attempted\'], B0_newly_adjudicated=8, B0_original_coronary_valid=12, totals=dict(totals), final_valid_protein_bootstrap_cells=len(records) - 20, final_failure_cells=len(failures), paired_coverages=[dict(indicator_id=r[\'indicator_id\'], contrast=r[\'contrast\'], paired_valid=r[\'paired_valid\']) for r in contrasts], final_all_series_have_1000_valid_replicates=all((r[\'bootstrap_valid\'] == 1000 for r in singles)), finite_joint_MLE_certified=False, additional_BH_computed=False, q_for_recovered_rows=None, models_fit_by_summarizer=0, new_draws=0, production_db_writes=0, all_13958_result_hashes_verified=True, all_4000_group_seals_and_20_worker_seals_verified=True, private_RDS_or_model_blobs_read_or_transferred=False, private_blob_hash_audit=\'EXPECTED_HASHES_RETAINED_IN_WORKER_SEALS_NOT_BULK_REHASHED_DURING_AGGREGATION\', manifest_path=str(manifest_path), manifest_sha256=manifestsha, recovery_policy_sha256=manifest[\'recovery_policy_sha256\'], phase1_adjudication_sha256=adjudicationsha, recovery_runtime_source_sha256=runtime_maps[0], recovery_runtime_vector_digest_sha256=next(iter(runtime_digests)), worker_runtime_map_field=\'recovery_runtime_sha256.1\', frozen_scientific_source_sha256=obj(worker_results[0].get(\'frozen_source_sha256\')), inference_kernels_source_sha256=sha(source_code_root / \'tools/summary_add9.py\'), quantile_Fieller_kernel_source_sha256=sha(source_code_root / \'inference/methods.py\'), **frozen)\n    return dict(schema=\'FOUR_PROTEIN_RECOVERY_AGGREGATE_V1\', complete=True, summary=summary, protein_single20=singles, protein_six_contrasts24=contrasts, failures=sorted(failures, key=lambda f: (f[\'indicator_id\'], f[\'bootstrap_id\'], f[\'endpoint_id\'])))\n\ndef csvread(path):\n    with Path(path).open(encoding=\'utf-8-sig\', newline=\'\') as stream:\n        return list(csv.DictReader(stream))\n\ndef csvwrite(path, data):\n    fields = list(dict.fromkeys((k for row in data for k in row)))\n    require(bool(fields), \'EMPTY_TABLE_REQUIRES_EXPLICIT_SCHEMA\')\n    with Path(path).open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as stream:\n        writer = csv.DictWriter(stream, fields)\n        writer.writeheader()\n        writer.writerows(data)\n\ndef display_tables(aggregate_result, review):\n    old_primary = csvread(review / \'DISPLAY18_PRIMARY_REVIEW.csv\')\n    old_single = csvread(review / \'DISPLAY18_SINGLE90_REVIEW.csv\')\n    replacements = {r[\'indicator_id\']: r for r in aggregate_result[\'protein_six_contrasts24\'] if r[\'contrast\'] == \'411_minus_250.2\'}\n    single_replacements = {(r[\'indicator_id\'], r[\'endpoint_id\']): r for r in aggregate_result[\'protein_single20\']}\n    primary, single = ([], [])\n    for old in old_primary:\n        if old[\'indicator_id\'] in IDS:\n            new = dict(replacements[old[\'indicator_id\']])\n            new.update(P=new.get(\'P_Delta_IE\'), q=None, q_family=\'NO_NEW_BH_FAMILY_AUTHORIZED\', current_count_point_pair_hash_verified=True, display_review_source=\'NEW_SELECTIVE_WARNING_RECOVERY\')\n        else:\n            new = dict(old, display_review_source=\'UNCHANGED_ACCEPTED_OLD9_OR_NEW5\')\n            require(all((new[k] == v for k, v in old.items())), \'NONPROTEIN_PRIMARY_ROW_CHANGED\')\n        primary.append(new)\n    for old in old_single:\n        if old[\'indicator_id\'] in IDS:\n            new = dict(single_replacements[old[\'indicator_id\'], old[\'endpoint_id\']], display_review_source=\'NEW_SELECTIVE_WARNING_RECOVERY\')\n        else:\n            new = dict(old, display_review_source=\'UNCHANGED_ACCEPTED_OLD9_OR_NEW5\')\n            require(all((new[k] == v for k, v in old.items())), \'NONPROTEIN_SINGLE_ROW_CHANGED\')\n        single.append(new)\n    require(len(primary) == 18 and len(single) == 90 and (len(replacements) == 4) and (len(single_replacements) == 20) and (len({r[\'indicator_id\'] for r in primary}) == 18) and (len({(r[\'indicator_id\'], r[\'endpoint_id\']) for r in single}) == 90), \'EXACT_UNIFIED_DISPLAY18_SCOPE_REQUIRED\')\n    return (primary, single)\n\ndef main():\n    parser = argparse.ArgumentParser(description=__doc__)\n    sub = parser.add_subparsers(dest=\'action\', required=True)\n    remote = sub.add_parser(\'remote-main\')\n    remote.add_argument(\'--manifest\', type=Path, required=True)\n    remote.add_argument(\'--recovery-root\', type=Path, required=True)\n    remote.add_argument(\'--source-code-root\', type=Path, default=Path(birth_private_setting(\'PRIVATE_BORDER_FULLFOLLOWUP_B1000_ADD9_20261030_V1\')))\n    remote.add_argument(\'--phase1-root\', type=Path, default=Path(birth_private_setting(\'PRIVATE_RUN_20261031_PHASE1_V1\')))\n    local = sub.add_parser(\'collect\')\n    local.add_argument(\'--remote-manifest\', required=True)\n    local.add_argument(\'--recovery-root\', required=True)\n    local.add_argument(\'--out\', type=Path, required=True)\n    local.add_argument(\'--review-dir\', type=Path, default=HERE.parents[2] / REVIEW_REL)\n    local.add_argument(\'--remote-code-root\', default=birth_private_setting(\'PRIVATE_BORDER_FOURPROTEIN_WARNING_RECOVERY_20261031_V1\'))\n    args = parser.parse_args()\n    if args.action == \'remote-main\':\n        value = aggregate(args.manifest, args.recovery_root, args.source_code_root, args.phase1_root)\n        print(json.dumps(value, ensure_ascii=False, separators=(\',\', \':\'), allow_nan=False))\n        return\n    require(not args.out.exists(), \'PRESERVE_EXISTING_LOCAL_RELEASE\')\n    from remote_control import remote as call_remote\n    code = \'import sys,json\\nsys.dont_write_bytecode=True\\nsys.path.insert(0,\' + repr(args.remote_code_root) + \')\\nimport summarize_recovery_schema_v2 as s\\nz=s.aggregate(\' + repr(args.remote_manifest) + \',\' + repr(args.recovery_root) + ",\'/mnt/d/scripts/border/border_fullfollowup_b1000_add9_20261030_v1\',\'/mnt/g/UKB_results/birth/border_fourprotein_warning_recovery_20261031_v1/run_20261031_phase1_v1\')\\nprint(json.dumps(z,ensure_ascii=False,separators=(\',\',\':\'),allow_nan=False))"\n    value = call_remote(code, timeout=1800)\n    require(value.get(\'complete\') is True, \'RECOVERY_NOT_COMPLETE_NO_RELEASE:\' + json.dumps(value, ensure_ascii=False))\n    primary, single = display_tables(value, args.review_dir)\n    args.out.mkdir(parents=True, exist_ok=False)\n    tables = {\'PROTEIN_SINGLE20_REVIEW.csv\': value[\'protein_single20\'], \'PROTEIN_SIX_CONTRASTS24_REVIEW.csv\': value[\'protein_six_contrasts24\'], \'PROTEIN_IHD_T2D_REVIEW.csv\': [r for r in value[\'protein_six_contrasts24\'] if r[\'contrast\'] == \'411_minus_250.2\'], \'ALL_PROTEIN_FAILURE_RECORDS.csv\': value[\'failures\'], \'DISPLAY18_PRIMARY_REVIEW.csv\': primary, \'DISPLAY18_SINGLE90_REVIEW.csv\': single}\n    for name, data in tables.items():\n        csvwrite(args.out / name, data)\n    summary = value[\'summary\']\n    summary[\'unchanged_nonprotein_display_indicators\'] = 14\n    summary[\'original_display_sources\'] = [dict(path=str((args.review_dir / name).resolve()), sha256=sha(args.review_dir / name)) for name in (\'DISPLAY18_PRIMARY_REVIEW.csv\', \'DISPLAY18_SINGLE90_REVIEW.csv\')]\n    with (args.out / \'RECOVERY_SUMMARY.json\').open(\'x\', encoding=\'utf-8\') as stream:\n        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)\n        stream.write(\'\\n\')\n    release = dict(schema=\'FOUR_PROTEIN_RECOVERY_LOCAL_RELEASE_MANIFEST_V1\', complete=True, files_sha256={name: sha(args.out / name) for name in [*tables, \'RECOVERY_SUMMARY.json\']}, source_code_sha256=sha(Path(__file__)), production_modified=False, models_fit=0, new_BH=False)\n    with (args.out / \'RELEASE_MANIFEST.json\').open(\'x\', encoding=\'utf-8\') as stream:\n        json.dump(release, stream, ensure_ascii=False, indent=2, allow_nan=False)\n        stream.write(\'\\n\')\n    print(json.dumps(dict(output=str(args.out.resolve()), table_rows={name: len(data) for name, data in tables.items()}, final_valid_protein_bootstrap_cells=summary[\'final_valid_protein_bootstrap_cells\'], final_failure_cells=summary[\'final_failure_cells\'])))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
