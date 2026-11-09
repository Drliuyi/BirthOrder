"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
MODEL = 'S44_BASE11'
DOMAINS = ['SMOKING', 'BODY_COMPOSITION', 'CLINICAL_LIPIDS', 'GLYCAEMIA']
ORDERS = ('forward', 'reverse')
ENDPOINTS = ('IHD', 'T2D', 'Delta')
INDEX_SCHEMA = 'GLYCAEMIA_SENSITIVITY_BASE11_B100_AGGREGATE_INDEX_V1'
TERMINAL = {'SUCCESS', 'STATISTICAL_FAILURE', 'STRUCTURAL_FAILURE', 'TECHNICAL_FAILURE', 'SYSTEM_FAILURE', 'INPUT_PREPARATION_FAILED'}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_json(path):
    path = Path(path)
    require(path.is_file() and path.stat().st_size < 40 * 1024 * 1024, 'SMALL_AGGREGATE_JSON_REQUIRED:' + str(path))
    with path.open(encoding='utf-8') as handle:
        return json.load(handle)

def finite(value):
    if value is None or isinstance(value, bool) or value == '':
        return None
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return number if math.isfinite(number) else None

def pinned(path, expected):
    path = Path(path).resolve()
    require(path.is_file() and sha(path) == expected, 'PINNED_SOURCE_CHANGED_OR_MISSING:' + str(path))
    return path

def module(paths, expected, name):
    path = next((Path(p).resolve() for p in paths if p is not None and Path(p).is_file()), None)
    require(path is not None, 'MISSING_FROZEN_MODULE:' + name)
    pinned(path, expected)
    spec = importlib.util.spec_from_file_location(name, path)
    item = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(item)
    return (item, {'path': str(path), 'sha256': expected})

def value(vectors, b, order, kind, domain, endpoint):
    vector = vectors[b][order]
    if endpoint == 'Delta':
        return vector[kind, domain, 'IHD'] - vector[kind, domain, 'T2D']
    return vector[kind, domain, endpoint]

def validate_aggregate(result, b, prior, expected_binding=None):
    """Transport adapter only; frozen validator owns scientific vector masks."""
    require(result.get('status') == 'SUCCESS', 'SAVED_SUCCESS_REQUIRED')
    require(result.get('model_id') == MODEL and result.get('b') == b, 'ACTUAL_SAME_MODEL_AND_BOOTSTRAP_ID_REQUIRED')
    require(result.get('participant_rows_in_reports', 0) == 0, 'AGGREGATE_ONLY_RESULT_REQUIRED')
    if b:
        require(result.get('schema') == RESULT_SCHEMA, 'B100_REPLICATE_SCHEMA_REQUIRED')
        require(type(result.get('seed')) is int and result['seed'] == 20261103 + 1009 * b, 'FROZEN_NUMBERED_REPLICATE_SEED_REQUIRED')
    else:
        require(result.get('schema') == 'GLYCAEMIA_SENSITIVITY_BASE11_B0_RESULT_V1', 'PINNED_NEW_SENSITIVITY_B0_SCHEMA_REQUIRED')
    binding = result.get('binding', {})
    require(isinstance(binding, dict) and binding.get('bootstrap_id') == b, 'ACTUAL_BOOTSTRAP_BINDING_REQUIRED')
    for key in ('sample_hash', 'population_hash', 'input_hash', 'model_spec_hash', 'rng_spec_hash'):
        require(isinstance(binding.get(key), str) and len(binding[key]) == 64, 'ACTUAL_SCIENTIFIC_BINDING_HASH_REQUIRED:' + key)
    for key, expected in (expected_binding or {}).items():
        require(binding.get(key) == expected, 'INDEX_FROZEN_BINDING_MISMATCH:' + key)
    audit = result.get('audit', {})
    kernel = audit.get('kernel', {})
    require(kernel.get('n_mc', kernel.get('actual_k')) == 50, 'ACTUAL_COMPLETE_K50_REQUIRED')
    require(audit.get('fresh_all_A') is True and audit.get('fresh_four_Y') is True, 'ACTUAL_FRESH_A_AND_FOUR_Y_REQUIRED')
    require(audit.get('old_A_Y_or_prediction_cache_reuse', False) is False, 'OLD_A_Y_OR_PREDICTION_CACHE_FORBIDDEN')
    require(isinstance(result.get('scenarios'), list) and result['scenarios'], 'SAVED_JOINT_PRODUCT_RISK_SCENARIOS_REQUIRED')
    adapted = dict(result, schema='NEW121_PRODUCT_RESULT_V1', panel_id=MODEL, bootstrap_id=b, scientific_binding=binding, audit={**audit, 'n_mc': 50})
    parsed, valid, errors = prior.validate_result(adapted, {'bootstrap_id': b}, {'panel_id': MODEL, 'domain_ids': DOMAINS, 'scientific_binding': binding})
    paired_checks = []
    for row in result.get('paired', []):
        order, kind = (row.get('allocation_order'), row.get('component_type'))
        domain = row.get('domain_id') if kind == 'IIE' else None
        observed = finite(row.get('Delta'))
        ihd, t2d = (parsed.get(order, {}).get((kind, domain, 'IHD')), parsed.get(order, {}).get((kind, domain, 'T2D')))
        expected = ihd - t2d if ihd is not None and t2d is not None else None
        paired_checks.append(dict(bootstrap_id=b, allocation_order=order, component_type=kind, domain_id=domain, saved_Delta=observed, derived_same_ID_Delta=expected, agrees_within_1e_10=observed is not None and expected is not None and (abs(observed - expected) <= 1e-10), validity_role='AUDIT_ONLY_FROZEN_VECTOR_SUBTRACTION_DEFINES_INFERENCE'))
    return (parsed, valid, errors, paired_checks)

def task_path(task, index_path):
    name = task.get('result_path')
    if not name:
        return None
    path = Path(name)
    return (index_path.parent / path).resolve() if not path.is_absolute() else path.resolve()

def load_inputs(index_path, b0_path, prior):
    index_path = Path(index_path).resolve()
    index = read_json(index_path)
    require(index.get('schema') == INDEX_SCHEMA, 'FIXED100_AGGREGATE_INDEX_SCHEMA_REQUIRED')
    tasks = index.get('tasks', [])
    require(isinstance(tasks, list) and len(tasks) == 100, 'EXACT100_NUMBERED_TASK_ROWS_REQUIRED')
    indexed = {}
    for task in tasks:
        b = task.get('b')
        require(type(b) is int and 1 <= b <= 100 and (b not in indexed), 'UNIQUE_FIXED_B1_TO_B100_REQUIRED')
        require(task.get('model_id', MODEL) == MODEL, 'BASE11_ONLY_INDEX_REQUIRED')
        require(task.get('status') in TERMINAL, 'ALL_FIXED_TASKS_MUST_BE_TERMINAL_BEFORE_RELEASE')
        indexed[b] = task
    require(set(indexed) == set(range(1, 101)), 'EXACT_FIXED_B1_TO_B100_COVERAGE_REQUIRED')
    if index.get('b0'):
        require(index['b0'].get('sha256') == B0_SHA, 'INDEX_NEW_B0_PIN_MISMATCH')
    pinned(b0_path, B0_SHA)
    b0 = read_json(b0_path)
    vectors, ledger, risks, pair_checks, diagnostics = ({}, [], [], [], [])
    success_results = {}
    for b in range(101):
        task = {'b': 0, 'status': 'SUCCESS', 'result_path': str(b0_path), 'sha256': B0_SHA} if b == 0 else indexed[b]
        path = Path(b0_path).resolve() if b == 0 else task_path(task, index_path)
        row = dict(bootstrap_id=b, status=task['status'], model_id=MODEL, failure_reason=task.get('failure_reason', task.get('reason')), source_path=str(path) if path else None, source_sha256=task.get('sha256'), forward_valid=False, reverse_valid=False, actual_k=None, transport_or_binding_error=None, order_validation_error=None, recorded_reverse_failure=False, reverse_failure_reason=None, order_success_json=json.dumps(task.get('order_success'), ensure_ascii=False), attempts_json=json.dumps(task.get('attempts'), ensure_ascii=False))
        result = {}
        try:
            if path is not None:
                require(task.get('sha256') and path.is_file() and (sha(path) == task['sha256']), 'ACTUAL_AGGREGATE_RESULT_SHA_MISMATCH')
                result = read_json(path)
                row['saved_result_status'] = result.get('status')
            if task['status'] == 'SUCCESS':
                require(path is not None, 'SUCCESS_CANNOT_HAVE_MISSING_RESULT')
                expected = {**index.get('required_result_binding', {}), **task.get('binding', {})} if b else {}
                if b:
                    for key, receipt in (('input_hash', index.get('input_sha256')), ('sample_hash', task.get('sample_hash'))):
                        if receipt is not None:
                            require(key not in expected or expected[key] == receipt, 'CONFLICTING_INDEX_BINDING_RECEIPTS:' + key)
                            expected[key] = receipt
                parsed, valid, errors, checks = validate_aggregate(result, b, prior, expected)
                if b and task.get('N_fit') is not None:
                    require(result.get('N', {}).get('N_fit') == task['N_fit'], 'INDEX_ACTUAL_FIT_POPULATION_COUNT_MISMATCH')
                vectors[b] = parsed
                pair_checks += checks
                row.update(forward_valid=valid['forward'], reverse_valid=valid['reverse'], actual_k=50, order_validation_error=json.dumps(errors, ensure_ascii=False) if any(errors.values()) else None)
                audit, kernel = (result.get('audit', {}), result.get('audit', {}).get('kernel', {}))
                order_success = result.get('order_success', task.get('order_success', {}))
                reverse_failure = kernel.get('reverse_failure')
                row.update(recorded_reverse_failure=bool(kernel.get('reverse_failure_recorded') or (isinstance(order_success, dict) and order_success.get('reverse') is False)), reverse_failure_reason=json.dumps(reverse_failure, ensure_ascii=False) if reverse_failure else None)
                row['order_success_json'] = json.dumps(order_success, ensure_ascii=False)
                success_results[b] = result
                risks += [dict(bootstrap_id=b, **risk) for risk in result.get('scenarios', [])]
                for diagnostic in result.get('mc_se', kernel.get('mc_se', [])):
                    diagnostics.append(dict(bootstrap_id=b, **diagnostic))
            elif path is None:
                require(row['failure_reason'], 'FAILED_TASK_WITHOUT_RESULT_REQUIRES_EXPLICIT_REASON')
        except (ValueError, KeyError, TypeError, OSError) as error:
            row['transport_or_binding_error'] = str(error)
            row['forward_valid'], row['reverse_valid'] = (False, False)
        ledger.append(row)
    require(ledger[0]['forward_valid'], 'PINNED_NEW_B0_FORWARD_NOT_VALID_NO_RELEASE')
    return (index, ledger, vectors, risks, pair_checks, diagnostics, success_results)

def normalize_effect(row):
    out = dict(row, model_id=MODEL, component_type=row['metric'], unit='probability')
    for key in ('point_B0', 'CI_type7_low', 'CI_type7_high'):
        out[key + '_per_1000'] = row[key]
        out[key] = row[key] / 1000 if row[key] is not None else None
    out['valid_ids'] = row['actual_ids']
    if out.get('test_role') == 'DESCRIPTIVE_ONLY':
        out['P_status'] = 'NOT_TESTED_UNDER_FROZEN_RULES_DESCRIPTIVE_CI_ONLY'
    else:
        out['P_status'] = 'FROZEN_ACTUAL_B_GRID_CI_INVERSION' if out.get('P') is not None else 'TOO_FEW_VALID_B_OR_MISSING_B0'
    out['BH_family'] = MODEL + '_FORWARD_DELTA_4_DOMAINS_PLUS_DEP' if out.get('BH_family_n') else None
    return out

def normalize_ratio(row):
    return dict(row, model_id=MODEL, P=None, q=None, test_role='DESCRIPTIVE_RATIO_ONLY_NO_NEW_TEST_FAMILY')

def review_rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return [row for row in csv.DictReader(handle) if row.get('model_id') == MODEL]

def frozen_original(reference_path, review_root):
    pinned(reference_path, REF_SHA)
    reference = read_json(reference_path)
    source_receipts = [{'path': str(Path(reference_path).resolve()), 'sha256': REF_SHA}]
    tables = {}
    for name, expected in REVIEW_SHAS.items():
        path = pinned(Path(review_root) / name, expected)
        tables[name] = review_rows(path)
        source_receipts.append({'path': str(path), 'sha256': expected})
    effects, ratios = ({}, {})
    for row in tables['FIVE_MODEL_OVERALL_EFFECTS.csv'] + tables['FIVE_MODEL_DOMAIN_EFFECTS.csv']:
        if row.get('allocation_order') != 'forward' or row.get('metric') not in ('TE', 'DE', 'IE_joint', 'IIE', 'DEP'):
            continue
        key = (row['metric'], row.get('domain_id') or None, row['endpoint'])
        require(key not in effects, 'ORIGINAL_REVIEW_DUPLICATE_EFFECT')
        effects[key] = row
        if row['metric'] in ('IIE', 'DEP'):
            ratio_name = ('domain_' if row['metric'] == 'IIE' else 'DEP_') + ('Psep' if row['endpoint'] == 'Delta' else 'PM')
            ratios[ratio_name, key[1], key[2]] = dict(point_B0=row.get('ratio_point_B0_percent'), CI_type7_low=row.get('ratio_CI_type7_low_percent'), CI_type7_high=row.get('ratio_CI_type7_high_percent'), effect_valid_N=row.get('effect_valid_N'), ratio_valid_N=row.get('ratio_valid_N'), denominator_unstable=row.get('denominator_unstable'), interval_status=row.get('ratio_interval_status'), source_path=row.get('source_path'))
    for row in tables['FIVE_MODEL_RATIOS.csv']:
        ratios[row['metric'], row.get('domain_id') or None, row['endpoint']] = dict(point_B0=row.get('point_percent'), CI_type7_low=row.get('CI_low_percent'), CI_type7_high=row.get('CI_high_percent'), ratio_valid_N=row.get('valid_N'), denominator_unstable=row.get('denominator_unstable'), source_path=row.get('source_path'))
    require(len(effects) == 24, 'EXACT_ORIGINAL_BASE11_NINE_TOTAL_FIFTEEN_DOMAIN_ROWS_REQUIRED')
    for (kind, domain, endpoint), row in effects.items():
        expected = reference['domain_IIE'][domain][endpoint] if kind == 'IIE' else reference['DEP'][endpoint] if kind == 'DEP' else reference['overall'][endpoint][kind]
        require(math.isclose(finite(row['point_B0']), expected, abs_tol=1e-15, rel_tol=1e-12), 'ORIGINAL_REFERENCE_AND_REVIEW_POINT_MISMATCH')
    return (effects, ratios, source_receipts)

def change(new, old):
    new, old = (finite(new), finite(old))
    difference = new - old if new is not None and old is not None else None
    difference = difference if difference is not None and math.isfinite(difference) else None
    relative = 100 * difference / old if difference is not None and old != 0 else None
    return (difference, relative if relative is not None and math.isfinite(relative) else None)

def compare_original(effects, ratios, original_effects, original_ratios):
    rows = []
    for row in effects:
        if row['allocation_order'] != 'forward':
            continue
        key = (row['metric'], row.get('domain_id'), row['endpoint'])
        old = original_effects.get(key, {})
        delta, relative = change(row['point_B0'], old.get('point_B0'))
        rows.append(dict(model_id=MODEL, allocation_order='forward', metric=row['metric'], domain_id=row.get('domain_id'), endpoint=row['endpoint'], unit='probability', original_point_B0=finite(old.get('point_B0')), new_point_B0=row['point_B0'], original_CI_type7_low=finite(old.get('CI_type7_low')), original_CI_type7_high=finite(old.get('CI_type7_high')), new_CI_type7_low=row['CI_type7_low'], new_CI_type7_high=row['CI_type7_high'], original_P=finite(old.get('P')), original_q=finite(old.get('q')), new_P=row.get('P'), new_q=row.get('q'), original_BH_family=old.get('family'), new_BH_family=row.get('BH_family'), original_valid_B=finite(old.get('effect_valid_N')), new_valid_B=row['effect_valid_N'], change_probability=delta, change_per_1000=1000 * delta if delta is not None else None, relative_change_percent=relative, relative_change_defined=relative is not None, signed_relative_formula='100*(new_B0-original_B0)/original_B0_NO_ABS', original_interval_source='HASH_PINNED_SAVED_REVIEW' if old else 'NOT_AVAILABLE', comparison_inference='DESCRIPTIVE_POINT_CHANGE_NO_CROSS_COHORT_CI_OR_TEST'))
    for row in ratios:
        if row['allocation_order'] != 'forward':
            continue
        key = (row['metric'], row.get('domain_id'), row['endpoint'])
        old = original_ratios.get(key, {})
        delta, relative = change(row['point_B0'], old.get('point_B0'))
        rows.append(dict(model_id=MODEL, allocation_order='forward', metric=row['metric'], domain_id=row.get('domain_id'), endpoint=row['endpoint'], unit='percent', original_point_B0=finite(old.get('point_B0')), new_point_B0=row['point_B0'], original_CI_type7_low=finite(old.get('CI_type7_low')), original_CI_type7_high=finite(old.get('CI_type7_high')), new_CI_type7_low=row['CI_type7_low'], new_CI_type7_high=row['CI_type7_high'], original_valid_B=finite(old.get('ratio_valid_N')), new_valid_B=row['ratio_valid_N'], original_denominator_unstable=old.get('denominator_unstable'), new_denominator_unstable=row['denominator_unstable'], change_percentage_points=delta, relative_change_percent=relative, relative_change_defined=relative is not None, signed_relative_formula='100*(new_B0-original_B0)/original_B0_NO_ABS', original_interval_source='HASH_PINNED_SAVED_REVIEW' if old else 'NOT_PUBLISHED_IN_FROZEN_REVIEW', comparison_inference='DESCRIPTIVE_POINT_CHANGE_NO_CROSS_COHORT_CI_OR_TEST'))
    return rows

def inference_details(rows, vectors, methods):
    details = []
    for row in rows:
        if row.get('test_role') not in ('PRIMARY_JOINT_DELTA_SEPARATE', 'FORWARD_DELTA_FIXED_G_PLUS_1_BH'):
            continue
        ids = [int(b) for b in row['actual_ids'].split(',') if b]
        values = [value(vectors, b, row['allocation_order'], row['metric'], row.get('domain_id'), row['endpoint']) for b in ids]
        info = methods.invert(values, row['point_B0']) if len(values) >= 2 and row['point_B0'] is not None else {}
        details.append(dict(model_id=MODEL, allocation_order=row['allocation_order'], component_type=row['metric'], domain_id=row.get('domain_id'), endpoint=row['endpoint'], point_B0=row['point_B0'], q=row.get('q'), BH_family=row.get('BH_family'), BH_family_n=row.get('BH_family_n'), actual_ids=row['actual_ids'], **info))
    return details

def count_status(index, ledger):
    bootstrap = [row for row in ledger if row['bootstrap_id'] > 0]
    counts = Counter((row['status'] for row in bootstrap))
    if 'attempt_failures' in index:
        indexed_attempts = list(index.get('attempt_failures') or [])
        attempt_scope = 'INDEX_ATTEMPT_FAILURES_ONLY_TASK_ATTEMPTS_NOT_ADDED'
    else:
        indexed_attempts = [attempt for task in index['tasks'] for attempt in task.get('attempts') or []]
        attempt_scope = 'TASK_INDEXED_ATTEMPTS_ONLY_NO_FAILURES_INFERRED'

    def exit_code(attempt):
        status = attempt.get('status', '')
        if isinstance(status, str) and status.startswith('EXIT_'):
            try:
                return int(status[5:])
            except ValueError:
                return None
        return None
    failed_attempts = [attempt for attempt in indexed_attempts if attempt.get('status') in TERMINAL - {'SUCCESS'} or exit_code(attempt) not in (None, 0)]
    return dict(terminal_status_counts={name: counts.get(name, 0) for name in sorted(TERMINAL)}, all_fixed100_terminal=True, terminal_task_N=len(bootstrap), statistical_terminal_failure_N=counts.get('STATISTICAL_FAILURE', 0), technical_terminal_failure_N=counts.get('TECHNICAL_FAILURE', 0), input_preparation_terminal_failure_N=counts.get('INPUT_PREPARATION_FAILED', 0), structural_terminal_failure_N=counts.get('STRUCTURAL_FAILURE', 0), system_terminal_failure_N=counts.get('SYSTEM_FAILURE', 0), aggregate_transport_or_binding_rejection_N=sum((bool(row['transport_or_binding_error']) for row in bootstrap)), successful_task_forward_invalid_N=sum((row['status'] == 'SUCCESS' and (not row['forward_valid']) for row in bootstrap)), successful_task_reverse_invalid_N=sum((row['status'] == 'SUCCESS' and (not row['reverse_valid']) for row in bootstrap)), recorded_reverse_allocation_failure_N=sum((row['recorded_reverse_failure'] for row in bootstrap)), reverse_only_failure_forward_retained_N=sum((row['forward_valid'] and (not row['reverse_valid']) for row in bootstrap)), attempt_failure_count_scope=attempt_scope, indexed_failed_attempts_N=len(failed_attempts), indexed_nonzero_process_exit_attempts_N=sum((exit_code(row) not in (None, 0) for row in indexed_attempts)), attempt_statistical_or_technical_class_inferred_from_exit_code=False, indexed_failed_attempt_status_counts=dict(Counter((row.get('status') for row in failed_attempts))))

def write_csv(path, rows):
    fields = list(dict.fromkeys((key for row in rows for key in row)))
    if not fields:
        fields, rows = (['status'], [{'status': 'NO_SAVED_AGGREGATE_ROWS'}])
    with Path(path).open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def write_json(path, value_):
    with Path(path).open('x', encoding='utf-8') as handle:
        json.dump(value_, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')

def release(args):
    workspace = Path(args.workspace_root).resolve()
    old_root = workspace / 'scripts/birth/border_selected44_stratified_20261103_v1'
    prior, prior_receipt = module([args.prior_reducer] if args.prior_reducer else [ROOT / 'vendor/prior121/tools/reduce.py', old_root / 'vendor/prior121/tools/reduce.py'], PRIOR_SHA, '_gs100_mature_reducer')
    methods, methods_receipt = module([args.methods] if args.methods else [ROOT / 'vendor/prior_inference/methods.py', old_root / 'vendor/prior_inference/methods.py'], METHOD_SHA, '_gs100_frozen_methods')
    b0 = Path(args.b0_result or workspace / 'reports/birth/GLYCAEMIA_SENSITIVITY_BASE11_B0_20261106_v1/remote_aggregate_evidence/RESULT.json').resolve()
    original_ref = Path(args.original_reference or workspace / 'scripts/birth/border_glycaemia_sensitivity_BASE11_B0_20261106_v1/tests/ORIGINAL_BASE11_B0_REFERENCE.json').resolve()
    review_root = Path(args.original_review_root or workspace / 'reports/birth/SELECTED44_ALL_RESULTS_REVIEW_20261105_v1/review_20261105T004052Z').resolve()
    output = Path(args.out).resolve()
    require(not output.exists(), 'PRESERVE_EXISTING_REPORT_CHOOSE_FRESH_VERSIONED_OUTPUT_DIRECTORY')
    index, ledger, vectors, risks, pair_checks, mc, results = load_inputs(args.index, b0, prior)
    old_effects, old_ratios, original_sources = frozen_original(original_ref, review_root)
    index_binding = index.get('binding') or {key: index[key] for key in ('input_sha256', 'config_sha256', 'draw_manifest_sha256', 'bootstrap_unit') if key in index}
    pseudo_manifest = {'panel_id': MODEL, 'domain_ids': DOMAINS, 'scientific_binding': index_binding}
    totals, domain_rows, raw, sensitivity, inference_audit = prior.summarize(pseudo_manifest, ledger, vectors, methods, 1e-06)
    effects = [normalize_effect(row) for row in totals if row['unit'] != 'percent']
    domains = [normalize_effect(row) for row in domain_rows if row['unit'] != 'percent']
    ratios = [normalize_ratio(row) for row in totals + domain_rows if row['unit'] == 'percent']
    require(inference_audit['BH_family_n'] == 5, 'FROZEN_EXACT5_BH_FAMILY_REQUIRED')
    original_comparison = compare_original(effects + domains, ratios, old_effects, old_ratios)
    detailed = inference_details(effects + domains, vectors, methods)
    population = []
    for b, result in results.items():
        n = result.get('N', {})
        population.append(dict(bootstrap_id=b, model_id=MODEL, source='ACTUAL_SAVED_REPLICATE', N_target=n.get('N_target'), N_fit=n.get('N_fit', n.get('N_joint_fit')), N_standardization=n.get('N_standardization'), N_excluded=n.get('N_excluded'), aggregate_population_and_events_json=json.dumps(n, ensure_ascii=False)))
    tables = {'RESULTS.csv': effects, 'DOMAIN.csv': domains, 'RATIOS.csv': ratios, 'ORIGINAL_BASE11_COMPARISON.csv': original_comparison, 'BOOTSTRAP_STATUS.csv': ledger, 'RAW_EFFECT_AND_RATIO_VECTORS.csv': raw, 'RAW_AGGREGATE_RISKS.csv': risks, 'SAVED_PAIRED_VECTOR_CHECKS.csv': pair_checks, 'FORWARD_REVERSE_SENSITIVITY.csv': sensitivity, 'FROZEN_INFERENCE_DETAILS.csv': detailed, 'POPULATION_FLOW.csv': population, 'MC_DIAGNOSTICS.csv': mc}
    counts = count_status(index, ledger)
    inference_audit['legacy_terminal_enum_metadata'] = {key: inference_audit[key] for key in ('fixed_ids_processed_to_terminal', 'all_fixed_ids_terminal')}
    inference_audit.update(fixed_ids_processed_to_terminal=counts['terminal_task_N'], all_fixed_ids_terminal=counts['all_fixed100_terminal'], terminal_enum_extension='INPUT_PREPARATION_FAILED_TERMINAL_NO_EFFECT_VECTOR')
    audit = dict(schema='GLYCAEMIA_SENSITIVITY_BASE11_B100_SUMMARY_V1', model_id=MODEL, source_index={'path': str(Path(args.index).resolve()), 'sha256': sha(args.index)}, source_new_B0={'path': str(b0), 'sha256': B0_SHA}, source_original_aggregate_review=original_sources, frozen_methods=methods_receipt, mature_reducer=prior_receipt, script_sha256=sha(__file__), fixed_B=100, inner_K=50, frozen_inference=inference_audit, failure_accounting=counts, bootstrap_point_estimator='NEW_SENSITIVITY_B0_NOT_BOOTSTRAP_MEAN', ratios_clipped=False, denominator_abs_applied=False, descriptive_endpoint_P_missing_per_freeze=True, new_test_families=0, original_CI_P_q_recomputed=False, original_comparison_inference='DESCRIPTIVE_CHANGE_ONLY_NO_PAIRED_CROSS_COHORT_CI_OR_P', publication_gate='ACTUAL_ALL_FIXED100_TERMINAL_NOT_SUCCESS_PERCENTAGE', participant_reads=False, participant_rows=0, model_reads=False, remote_access=False, new_model_fits=0, new_integrations=0, numerical_pilot_tests=0)
    output.mkdir(parents=True, exist_ok=False)
    for name, rows in tables.items():
        write_csv(output / name, rows)
    rules_path = ROOT / 'docs/FROZEN_INFERENCE_RULES.md'
    with (output / 'METHODS.md').open('x', encoding='utf-8') as handle:
        handle.write(rules_path.read_text(encoding='utf-8'))
    fwd, rev = (inference_audit['forward_effect_valid_N'], inference_audit['reverse_effect_valid_N'])
    report = f"# BASE11 HbA1c敏感性：固定B100聚合结果\n\n本目录仅对实际完成的固定B1–100终态记录归约。点估计固定取已完成的新敏感性B0，不是bootstrap均值；原BASE11结果未改动。\n\n固定100编号全部终态；正序两病共同完整有效B={fwd}，反序有效B={rev}。统计终态失败={counts['statistical_terminal_failure_N']}；技术终态失败={counts['technical_terminal_failure_N']}；输入准备终态失败={counts['input_preparation_terminal_failure_N']}；结构终态失败={counts['structural_terminal_failure_N']}；系统终态失败={counts['system_terminal_failure_N']}；聚合传输/绑定校验拒绝={counts['aggregate_transport_or_binding_rejection_N']}；仅反序失败而正序保留={counts['reverse_only_failure_forward_retained_N']}。\n\nRESULTS.csv列示IHD、T2D及同编号配对Delta的TE/DE/IE_joint；DOMAIN.csv列示四域IIE和DEP。RATIOS.csv保留全部逐编号PM/Psep与分母稳定性。正反序各自使用自身完整编号；FORWARD_REVERSE_SENSITIVITY.csv仅用交集比较。\n\nCI严格沿用Type-7百分位。P仅按冻结规则用于正序Delta joint IE及四域+DEP；四域+DEP固定五项BH，joint IE单列。IHD/T2D、TE/DE、比例、反序P/q为空是冻结规则，并非漏算。P是实际有效B的离散CI反演筛查值，不是精确尾概率；FROZEN_INFERENCE_DETAILS.csv公开实际B、网格步长、搜索括号与限制。\n\n比例分母保留符号，负值或超过100%不截断。零/非有限只排除该比例，接近零、变号或CI跨零均标记；分母不稳定时有限比例分位数不代表可靠解释率区间。\n\nORIGINAL_BASE11_COMPARISON.csv包含原结果已冻结CI/P/q及新B0的绝对/有符号相对变化；比例绝对变化单位为百分点。两队列变化仅描述，不创造配对变化CI/P。未在冻结旧报告发布的指标以NOT_PUBLISHED标明，不回算旧区间。\n\n详见METHODS.md、SUMMARY_AUDIT.json及逐编号BOOTSTRAP_STATUS.csv。本聚合没有参与者级读取、模型拟合或新增积分。\n"
    with (output / 'REPORT_ZH.md').open('x', encoding='utf-8') as handle:
        handle.write(report)
    audit['outputs'] = {name: {'sha256': sha(output / name), 'scientific_rows': len(rows)} for name, rows in tables.items()}
    audit['outputs'].update({name: {'sha256': sha(output / name)} for name in ('METHODS.md', 'REPORT_ZH.md')})
    write_json(output / 'SUMMARY_AUDIT.json', audit)
    write_json(output / 'OUTPUT_MANIFEST.json', dict(schema='AGGREGATE_EXCLUSIVE_OUTPUT_MANIFEST_V1', files={path.name: sha(path) for path in sorted(output.iterdir()) if path.is_file()}, old_results_unchanged_by_aggregate_tool=True))
    return {'output_root': str(output), 'forward_valid_B': fwd, 'reverse_valid_B': rev, 'fixed100_all_terminal': True, 'failure_accounting': counts}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', default=os.environ.get('AGGREGATE_INDEX'), required=not os.environ.get('AGGREGATE_INDEX'))
    parser.add_argument('--out', '--report-root', dest='out', default=os.environ.get('REPORT_OUTPUT_ROOT'), required=not os.environ.get('REPORT_OUTPUT_ROOT'))
    parser.add_argument('--workspace-root', default=os.environ.get('WORKSPACE_ROOT', str(WORKSPACE)))
    parser.add_argument('--b0-result', '--b0', dest='b0_result', default=os.environ.get('NEW_B0_RESULT'))
    parser.add_argument('--original-reference', '--reference', dest='original_reference', default=os.environ.get('ORIGINAL_B0_REFERENCE'))
    parser.add_argument('--original-review-root', default=os.environ.get('ORIGINAL_REVIEW_ROOT'))
    parser.add_argument('--methods', type=Path)
    parser.add_argument('--prior-reducer', type=Path)
    args = parser.parse_args()
    print(json.dumps({'aggregate_input_root': str(Path(args.index).resolve().parent), 'script_root': str(ROOT), 'output_root': str(Path(args.out).resolve()), 'workspace_root': str(Path(args.workspace_root).resolve())}))
    print(json.dumps(release(args), ensure_ascii=False, allow_nan=False))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Release aggregate-only BASE11 HbA1c B100 using hash-pinned frozen inference.\n\nNo participant/model/RDS reads, fitting, integration, bootstrap generation or\nremote access. Only release against all 100 actual terminal numbered tasks.\nNo numerical pilot/self-test is embedded or run on import.\n"""\nfrom __future__ import annotations\nimport argparse\nfrom collections import Counter\nimport csv\nimport hashlib\nimport importlib.util\nimport json\nimport math\nimport os\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\nROOT = Path(__file__).resolve().parents[1]\nWORKSPACE = ROOT.parents[2]\nMODEL = \'S44_BASE11\'\nDOMAINS = [\'SMOKING\', \'BODY_COMPOSITION\', \'CLINICAL_LIPIDS\', \'GLYCAEMIA\']\nORDERS = (\'forward\', \'reverse\')\nENDPOINTS = (\'IHD\', \'T2D\', \'Delta\')\nINDEX_SCHEMA = \'GLYCAEMIA_SENSITIVITY_BASE11_B100_AGGREGATE_INDEX_V1\'\nRESULT_SCHEMA = \'GLYCAEMIA_SENSITIVITY_BASE11_B100_REPLICATE_V1\'\nTERMINAL = {\'SUCCESS\', \'STATISTICAL_FAILURE\', \'STRUCTURAL_FAILURE\', \'TECHNICAL_FAILURE\', \'SYSTEM_FAILURE\', \'INPUT_PREPARATION_FAILED\'}\nMETHOD_SHA = birth_private_setting(\'PRIVATE_A6390305CC75B51962884894E004082FF02D67BBA2B0E9F0946893A19DD0\')\nPRIOR_SHA = birth_private_setting(\'PRIVATE_0C218D2EEE0898966EEC1274832A638D68C18AFFCDC340F8B59C54DA9A9E\')\nB0_SHA = birth_private_setting(\'PRIVATE_CC4FD25EF3C2D5D2CFF6E6A9297BE332329632017BF4BD218407CE9CA400\')\nREF_SHA = birth_private_setting(\'PRIVATE_E49F4321E09E1154F2218C1293253C1A2D924713126F129A6D364B32E3FD\')\nREVIEW_SHAS = {\'FIVE_MODEL_OVERALL_EFFECTS.csv\': birth_private_setting(\'PRIVATE_574A418D92A89BFC368791C2BEF6227F8B0CBCE771BE6D009AA94D00F2C8\'), \'FIVE_MODEL_DOMAIN_EFFECTS.csv\': birth_private_setting(\'PRIVATE_83A642F02CFA4A074AF87B8B4B956011E0AA4FB0796AF16803ACE0CB2B7F\'), \'FIVE_MODEL_RATIOS.csv\': birth_private_setting(\'PRIVATE_6042736BE9FFE75C95ACF92DB05EE8C393E337BA16222377C83D2B52055D\')}\n\ndef require(ok, message):\n    if not ok:\n        raise ValueError(message)\n\ndef sha(path):\n    return hashlib.sha256(Path(path).read_bytes()).hexdigest()\n\ndef read_json(path):\n    path = Path(path)\n    require(path.is_file() and path.stat().st_size < 40 * 1024 * 1024, \'SMALL_AGGREGATE_JSON_REQUIRED:\' + str(path))\n    with path.open(encoding=\'utf-8\') as handle:\n        return json.load(handle)\n\ndef finite(value):\n    if value is None or isinstance(value, bool) or value == \'\':\n        return None\n    try:\n        number = float(value)\n    except (ValueError, TypeError, OverflowError):\n        return None\n    return number if math.isfinite(number) else None\n\ndef pinned(path, expected):\n    path = Path(path).resolve()\n    require(path.is_file() and sha(path) == expected, \'PINNED_SOURCE_CHANGED_OR_MISSING:\' + str(path))\n    return path\n\ndef module(paths, expected, name):\n    path = next((Path(p).resolve() for p in paths if p is not None and Path(p).is_file()), None)\n    require(path is not None, \'MISSING_FROZEN_MODULE:\' + name)\n    pinned(path, expected)\n    spec = importlib.util.spec_from_file_location(name, path)\n    item = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(item)\n    return (item, {\'path\': str(path), \'sha256\': expected})\n\ndef value(vectors, b, order, kind, domain, endpoint):\n    vector = vectors[b][order]\n    if endpoint == \'Delta\':\n        return vector[kind, domain, \'IHD\'] - vector[kind, domain, \'T2D\']\n    return vector[kind, domain, endpoint]\n\ndef validate_aggregate(result, b, prior, expected_binding=None):\n    """Transport adapter only; frozen validator owns scientific vector masks."""\n    require(result.get(\'status\') == \'SUCCESS\', \'SAVED_SUCCESS_REQUIRED\')\n    require(result.get(\'model_id\') == MODEL and result.get(\'b\') == b, \'ACTUAL_SAME_MODEL_AND_BOOTSTRAP_ID_REQUIRED\')\n    require(result.get(\'participant_rows_in_reports\', 0) == 0, \'AGGREGATE_ONLY_RESULT_REQUIRED\')\n    if b:\n        require(result.get(\'schema\') == RESULT_SCHEMA, \'B100_REPLICATE_SCHEMA_REQUIRED\')\n        require(type(result.get(\'seed\')) is int and result[\'seed\'] == 20261103 + 1009 * b, \'FROZEN_NUMBERED_REPLICATE_SEED_REQUIRED\')\n    else:\n        require(result.get(\'schema\') == \'GLYCAEMIA_SENSITIVITY_BASE11_B0_RESULT_V1\', \'PINNED_NEW_SENSITIVITY_B0_SCHEMA_REQUIRED\')\n    binding = result.get(\'binding\', {})\n    require(isinstance(binding, dict) and binding.get(\'bootstrap_id\') == b, \'ACTUAL_BOOTSTRAP_BINDING_REQUIRED\')\n    for key in (\'sample_hash\', \'population_hash\', \'input_hash\', \'model_spec_hash\', \'rng_spec_hash\'):\n        require(isinstance(binding.get(key), str) and len(binding[key]) == 64, \'ACTUAL_SCIENTIFIC_BINDING_HASH_REQUIRED:\' + key)\n    for key, expected in (expected_binding or {}).items():\n        require(binding.get(key) == expected, \'INDEX_FROZEN_BINDING_MISMATCH:\' + key)\n    audit = result.get(\'audit\', {})\n    kernel = audit.get(\'kernel\', {})\n    require(kernel.get(\'n_mc\', kernel.get(\'actual_k\')) == 50, \'ACTUAL_COMPLETE_K50_REQUIRED\')\n    require(audit.get(\'fresh_all_A\') is True and audit.get(\'fresh_four_Y\') is True, \'ACTUAL_FRESH_A_AND_FOUR_Y_REQUIRED\')\n    require(audit.get(\'old_A_Y_or_prediction_cache_reuse\', False) is False, \'OLD_A_Y_OR_PREDICTION_CACHE_FORBIDDEN\')\n    require(isinstance(result.get(\'scenarios\'), list) and result[\'scenarios\'], \'SAVED_JOINT_PRODUCT_RISK_SCENARIOS_REQUIRED\')\n    adapted = dict(result, schema=\'NEW121_PRODUCT_RESULT_V1\', panel_id=MODEL, bootstrap_id=b, scientific_binding=binding, audit={**audit, \'n_mc\': 50})\n    parsed, valid, errors = prior.validate_result(adapted, {\'bootstrap_id\': b}, {\'panel_id\': MODEL, \'domain_ids\': DOMAINS, \'scientific_binding\': binding})\n    paired_checks = []\n    for row in result.get(\'paired\', []):\n        order, kind = (row.get(\'allocation_order\'), row.get(\'component_type\'))\n        domain = row.get(\'domain_id\') if kind == \'IIE\' else None\n        observed = finite(row.get(\'Delta\'))\n        ihd, t2d = (parsed.get(order, {}).get((kind, domain, \'IHD\')), parsed.get(order, {}).get((kind, domain, \'T2D\')))\n        expected = ihd - t2d if ihd is not None and t2d is not None else None\n        paired_checks.append(dict(bootstrap_id=b, allocation_order=order, component_type=kind, domain_id=domain, saved_Delta=observed, derived_same_ID_Delta=expected, agrees_within_1e_10=observed is not None and expected is not None and (abs(observed - expected) <= 1e-10), validity_role=\'AUDIT_ONLY_FROZEN_VECTOR_SUBTRACTION_DEFINES_INFERENCE\'))\n    return (parsed, valid, errors, paired_checks)\n\ndef task_path(task, index_path):\n    name = task.get(\'result_path\')\n    if not name:\n        return None\n    path = Path(name)\n    return (index_path.parent / path).resolve() if not path.is_absolute() else path.resolve()\n\ndef load_inputs(index_path, b0_path, prior):\n    index_path = Path(index_path).resolve()\n    index = read_json(index_path)\n    require(index.get(\'schema\') == INDEX_SCHEMA, \'FIXED100_AGGREGATE_INDEX_SCHEMA_REQUIRED\')\n    tasks = index.get(\'tasks\', [])\n    require(isinstance(tasks, list) and len(tasks) == 100, \'EXACT100_NUMBERED_TASK_ROWS_REQUIRED\')\n    indexed = {}\n    for task in tasks:\n        b = task.get(\'b\')\n        require(type(b) is int and 1 <= b <= 100 and (b not in indexed), \'UNIQUE_FIXED_B1_TO_B100_REQUIRED\')\n        require(task.get(\'model_id\', MODEL) == MODEL, \'BASE11_ONLY_INDEX_REQUIRED\')\n        require(task.get(\'status\') in TERMINAL, \'ALL_FIXED_TASKS_MUST_BE_TERMINAL_BEFORE_RELEASE\')\n        indexed[b] = task\n    require(set(indexed) == set(range(1, 101)), \'EXACT_FIXED_B1_TO_B100_COVERAGE_REQUIRED\')\n    if index.get(\'b0\'):\n        require(index[\'b0\'].get(\'sha256\') == B0_SHA, \'INDEX_NEW_B0_PIN_MISMATCH\')\n    pinned(b0_path, B0_SHA)\n    b0 = read_json(b0_path)\n    vectors, ledger, risks, pair_checks, diagnostics = ({}, [], [], [], [])\n    success_results = {}\n    for b in range(101):\n        task = {\'b\': 0, \'status\': \'SUCCESS\', \'result_path\': str(b0_path), \'sha256\': B0_SHA} if b == 0 else indexed[b]\n        path = Path(b0_path).resolve() if b == 0 else task_path(task, index_path)\n        row = dict(bootstrap_id=b, status=task[\'status\'], model_id=MODEL, failure_reason=task.get(\'failure_reason\', task.get(\'reason\')), source_path=str(path) if path else None, source_sha256=task.get(\'sha256\'), forward_valid=False, reverse_valid=False, actual_k=None, transport_or_binding_error=None, order_validation_error=None, recorded_reverse_failure=False, reverse_failure_reason=None, order_success_json=json.dumps(task.get(\'order_success\'), ensure_ascii=False), attempts_json=json.dumps(task.get(\'attempts\'), ensure_ascii=False))\n        result = {}\n        try:\n            if path is not None:\n                require(task.get(\'sha256\') and path.is_file() and (sha(path) == task[\'sha256\']), \'ACTUAL_AGGREGATE_RESULT_SHA_MISMATCH\')\n                result = read_json(path)\n                row[\'saved_result_status\'] = result.get(\'status\')\n            if task[\'status\'] == \'SUCCESS\':\n                require(path is not None, \'SUCCESS_CANNOT_HAVE_MISSING_RESULT\')\n                expected = {**index.get(\'required_result_binding\', {}), **task.get(\'binding\', {})} if b else {}\n                if b:\n                    for key, receipt in ((\'input_hash\', index.get(\'input_sha256\')), (\'sample_hash\', task.get(\'sample_hash\'))):\n                        if receipt is not None:\n                            require(key not in expected or expected[key] == receipt, \'CONFLICTING_INDEX_BINDING_RECEIPTS:\' + key)\n                            expected[key] = receipt\n                parsed, valid, errors, checks = validate_aggregate(result, b, prior, expected)\n                if b and task.get(\'N_fit\') is not None:\n                    require(result.get(\'N\', {}).get(\'N_fit\') == task[\'N_fit\'], \'INDEX_ACTUAL_FIT_POPULATION_COUNT_MISMATCH\')\n                vectors[b] = parsed\n                pair_checks += checks\n                row.update(forward_valid=valid[\'forward\'], reverse_valid=valid[\'reverse\'], actual_k=50, order_validation_error=json.dumps(errors, ensure_ascii=False) if any(errors.values()) else None)\n                audit, kernel = (result.get(\'audit\', {}), result.get(\'audit\', {}).get(\'kernel\', {}))\n                order_success = result.get(\'order_success\', task.get(\'order_success\', {}))\n                reverse_failure = kernel.get(\'reverse_failure\')\n                row.update(recorded_reverse_failure=bool(kernel.get(\'reverse_failure_recorded\') or (isinstance(order_success, dict) and order_success.get(\'reverse\') is False)), reverse_failure_reason=json.dumps(reverse_failure, ensure_ascii=False) if reverse_failure else None)\n                row[\'order_success_json\'] = json.dumps(order_success, ensure_ascii=False)\n                success_results[b] = result\n                risks += [dict(bootstrap_id=b, **risk) for risk in result.get(\'scenarios\', [])]\n                for diagnostic in result.get(\'mc_se\', kernel.get(\'mc_se\', [])):\n                    diagnostics.append(dict(bootstrap_id=b, **diagnostic))\n            elif path is None:\n                require(row[\'failure_reason\'], \'FAILED_TASK_WITHOUT_RESULT_REQUIRES_EXPLICIT_REASON\')\n        except (ValueError, KeyError, TypeError, OSError) as error:\n            row[\'transport_or_binding_error\'] = str(error)\n            row[\'forward_valid\'], row[\'reverse_valid\'] = (False, False)\n        ledger.append(row)\n    require(ledger[0][\'forward_valid\'], \'PINNED_NEW_B0_FORWARD_NOT_VALID_NO_RELEASE\')\n    return (index, ledger, vectors, risks, pair_checks, diagnostics, success_results)\n\ndef normalize_effect(row):\n    out = dict(row, model_id=MODEL, component_type=row[\'metric\'], unit=\'probability\')\n    for key in (\'point_B0\', \'CI_type7_low\', \'CI_type7_high\'):\n        out[key + \'_per_1000\'] = row[key]\n        out[key] = row[key] / 1000 if row[key] is not None else None\n    out[\'valid_ids\'] = row[\'actual_ids\']\n    if out.get(\'test_role\') == \'DESCRIPTIVE_ONLY\':\n        out[\'P_status\'] = \'NOT_TESTED_UNDER_FROZEN_RULES_DESCRIPTIVE_CI_ONLY\'\n    else:\n        out[\'P_status\'] = \'FROZEN_ACTUAL_B_GRID_CI_INVERSION\' if out.get(\'P\') is not None else \'TOO_FEW_VALID_B_OR_MISSING_B0\'\n    out[\'BH_family\'] = MODEL + \'_FORWARD_DELTA_4_DOMAINS_PLUS_DEP\' if out.get(\'BH_family_n\') else None\n    return out\n\ndef normalize_ratio(row):\n    return dict(row, model_id=MODEL, P=None, q=None, test_role=\'DESCRIPTIVE_RATIO_ONLY_NO_NEW_TEST_FAMILY\')\n\ndef review_rows(path):\n    with Path(path).open(encoding=\'utf-8-sig\', newline=\'\') as handle:\n        return [row for row in csv.DictReader(handle) if row.get(\'model_id\') == MODEL]\n\ndef frozen_original(reference_path, review_root):\n    pinned(reference_path, REF_SHA)\n    reference = read_json(reference_path)\n    source_receipts = [{\'path\': str(Path(reference_path).resolve()), \'sha256\': REF_SHA}]\n    tables = {}\n    for name, expected in REVIEW_SHAS.items():\n        path = pinned(Path(review_root) / name, expected)\n        tables[name] = review_rows(path)\n        source_receipts.append({\'path\': str(path), \'sha256\': expected})\n    effects, ratios = ({}, {})\n    for row in tables[\'FIVE_MODEL_OVERALL_EFFECTS.csv\'] + tables[\'FIVE_MODEL_DOMAIN_EFFECTS.csv\']:\n        if row.get(\'allocation_order\') != \'forward\' or row.get(\'metric\') not in (\'TE\', \'DE\', \'IE_joint\', \'IIE\', \'DEP\'):\n            continue\n        key = (row[\'metric\'], row.get(\'domain_id\') or None, row[\'endpoint\'])\n        require(key not in effects, \'ORIGINAL_REVIEW_DUPLICATE_EFFECT\')\n        effects[key] = row\n        if row[\'metric\'] in (\'IIE\', \'DEP\'):\n            ratio_name = (\'domain_\' if row[\'metric\'] == \'IIE\' else \'DEP_\') + (\'Psep\' if row[\'endpoint\'] == \'Delta\' else \'PM\')\n            ratios[ratio_name, key[1], key[2]] = dict(point_B0=row.get(\'ratio_point_B0_percent\'), CI_type7_low=row.get(\'ratio_CI_type7_low_percent\'), CI_type7_high=row.get(\'ratio_CI_type7_high_percent\'), effect_valid_N=row.get(\'effect_valid_N\'), ratio_valid_N=row.get(\'ratio_valid_N\'), denominator_unstable=row.get(\'denominator_unstable\'), interval_status=row.get(\'ratio_interval_status\'), source_path=row.get(\'source_path\'))\n    for row in tables[\'FIVE_MODEL_RATIOS.csv\']:\n        ratios[row[\'metric\'], row.get(\'domain_id\') or None, row[\'endpoint\']] = dict(point_B0=row.get(\'point_percent\'), CI_type7_low=row.get(\'CI_low_percent\'), CI_type7_high=row.get(\'CI_high_percent\'), ratio_valid_N=row.get(\'valid_N\'), denominator_unstable=row.get(\'denominator_unstable\'), source_path=row.get(\'source_path\'))\n    require(len(effects) == 24, \'EXACT_ORIGINAL_BASE11_NINE_TOTAL_FIFTEEN_DOMAIN_ROWS_REQUIRED\')\n    for (kind, domain, endpoint), row in effects.items():\n        expected = reference[\'domain_IIE\'][domain][endpoint] if kind == \'IIE\' else reference[\'DEP\'][endpoint] if kind == \'DEP\' else reference[\'overall\'][endpoint][kind]\n        require(math.isclose(finite(row[\'point_B0\']), expected, abs_tol=1e-15, rel_tol=1e-12), \'ORIGINAL_REFERENCE_AND_REVIEW_POINT_MISMATCH\')\n    return (effects, ratios, source_receipts)\n\ndef change(new, old):\n    new, old = (finite(new), finite(old))\n    difference = new - old if new is not None and old is not None else None\n    difference = difference if difference is not None and math.isfinite(difference) else None\n    relative = 100 * difference / old if difference is not None and old != 0 else None\n    return (difference, relative if relative is not None and math.isfinite(relative) else None)\n\ndef compare_original(effects, ratios, original_effects, original_ratios):\n    rows = []\n    for row in effects:\n        if row[\'allocation_order\'] != \'forward\':\n            continue\n        key = (row[\'metric\'], row.get(\'domain_id\'), row[\'endpoint\'])\n        old = original_effects.get(key, {})\n        delta, relative = change(row[\'point_B0\'], old.get(\'point_B0\'))\n        rows.append(dict(model_id=MODEL, allocation_order=\'forward\', metric=row[\'metric\'], domain_id=row.get(\'domain_id\'), endpoint=row[\'endpoint\'], unit=\'probability\', original_point_B0=finite(old.get(\'point_B0\')), new_point_B0=row[\'point_B0\'], original_CI_type7_low=finite(old.get(\'CI_type7_low\')), original_CI_type7_high=finite(old.get(\'CI_type7_high\')), new_CI_type7_low=row[\'CI_type7_low\'], new_CI_type7_high=row[\'CI_type7_high\'], original_P=finite(old.get(\'P\')), original_q=finite(old.get(\'q\')), new_P=row.get(\'P\'), new_q=row.get(\'q\'), original_BH_family=old.get(\'family\'), new_BH_family=row.get(\'BH_family\'), original_valid_B=finite(old.get(\'effect_valid_N\')), new_valid_B=row[\'effect_valid_N\'], change_probability=delta, change_per_1000=1000 * delta if delta is not None else None, relative_change_percent=relative, relative_change_defined=relative is not None, signed_relative_formula=\'100*(new_B0-original_B0)/original_B0_NO_ABS\', original_interval_source=\'HASH_PINNED_SAVED_REVIEW\' if old else \'NOT_AVAILABLE\', comparison_inference=\'DESCRIPTIVE_POINT_CHANGE_NO_CROSS_COHORT_CI_OR_TEST\'))\n    for row in ratios:\n        if row[\'allocation_order\'] != \'forward\':\n            continue\n        key = (row[\'metric\'], row.get(\'domain_id\'), row[\'endpoint\'])\n        old = original_ratios.get(key, {})\n        delta, relative = change(row[\'point_B0\'], old.get(\'point_B0\'))\n        rows.append(dict(model_id=MODEL, allocation_order=\'forward\', metric=row[\'metric\'], domain_id=row.get(\'domain_id\'), endpoint=row[\'endpoint\'], unit=\'percent\', original_point_B0=finite(old.get(\'point_B0\')), new_point_B0=row[\'point_B0\'], original_CI_type7_low=finite(old.get(\'CI_type7_low\')), original_CI_type7_high=finite(old.get(\'CI_type7_high\')), new_CI_type7_low=row[\'CI_type7_low\'], new_CI_type7_high=row[\'CI_type7_high\'], original_valid_B=finite(old.get(\'ratio_valid_N\')), new_valid_B=row[\'ratio_valid_N\'], original_denominator_unstable=old.get(\'denominator_unstable\'), new_denominator_unstable=row[\'denominator_unstable\'], change_percentage_points=delta, relative_change_percent=relative, relative_change_defined=relative is not None, signed_relative_formula=\'100*(new_B0-original_B0)/original_B0_NO_ABS\', original_interval_source=\'HASH_PINNED_SAVED_REVIEW\' if old else \'NOT_PUBLISHED_IN_FROZEN_REVIEW\', comparison_inference=\'DESCRIPTIVE_POINT_CHANGE_NO_CROSS_COHORT_CI_OR_TEST\'))\n    return rows\n\ndef inference_details(rows, vectors, methods):\n    details = []\n    for row in rows:\n        if row.get(\'test_role\') not in (\'PRIMARY_JOINT_DELTA_SEPARATE\', \'FORWARD_DELTA_FIXED_G_PLUS_1_BH\'):\n            continue\n        ids = [int(b) for b in row[\'actual_ids\'].split(\',\') if b]\n        values = [value(vectors, b, row[\'allocation_order\'], row[\'metric\'], row.get(\'domain_id\'), row[\'endpoint\']) for b in ids]\n        info = methods.invert(values, row[\'point_B0\']) if len(values) >= 2 and row[\'point_B0\'] is not None else {}\n        details.append(dict(model_id=MODEL, allocation_order=row[\'allocation_order\'], component_type=row[\'metric\'], domain_id=row.get(\'domain_id\'), endpoint=row[\'endpoint\'], point_B0=row[\'point_B0\'], q=row.get(\'q\'), BH_family=row.get(\'BH_family\'), BH_family_n=row.get(\'BH_family_n\'), actual_ids=row[\'actual_ids\'], **info))\n    return details\n\ndef count_status(index, ledger):\n    bootstrap = [row for row in ledger if row[\'bootstrap_id\'] > 0]\n    counts = Counter((row[\'status\'] for row in bootstrap))\n    if \'attempt_failures\' in index:\n        indexed_attempts = list(index.get(\'attempt_failures\') or [])\n        attempt_scope = \'INDEX_ATTEMPT_FAILURES_ONLY_TASK_ATTEMPTS_NOT_ADDED\'\n    else:\n        indexed_attempts = [attempt for task in index[\'tasks\'] for attempt in task.get(\'attempts\') or []]\n        attempt_scope = \'TASK_INDEXED_ATTEMPTS_ONLY_NO_FAILURES_INFERRED\'\n\n    def exit_code(attempt):\n        status = attempt.get(\'status\', \'\')\n        if isinstance(status, str) and status.startswith(\'EXIT_\'):\n            try:\n                return int(status[5:])\n            except ValueError:\n                return None\n        return None\n    failed_attempts = [attempt for attempt in indexed_attempts if attempt.get(\'status\') in TERMINAL - {\'SUCCESS\'} or exit_code(attempt) not in (None, 0)]\n    return dict(terminal_status_counts={name: counts.get(name, 0) for name in sorted(TERMINAL)}, all_fixed100_terminal=True, terminal_task_N=len(bootstrap), statistical_terminal_failure_N=counts.get(\'STATISTICAL_FAILURE\', 0), technical_terminal_failure_N=counts.get(\'TECHNICAL_FAILURE\', 0), input_preparation_terminal_failure_N=counts.get(\'INPUT_PREPARATION_FAILED\', 0), structural_terminal_failure_N=counts.get(\'STRUCTURAL_FAILURE\', 0), system_terminal_failure_N=counts.get(\'SYSTEM_FAILURE\', 0), aggregate_transport_or_binding_rejection_N=sum((bool(row[\'transport_or_binding_error\']) for row in bootstrap)), successful_task_forward_invalid_N=sum((row[\'status\'] == \'SUCCESS\' and (not row[\'forward_valid\']) for row in bootstrap)), successful_task_reverse_invalid_N=sum((row[\'status\'] == \'SUCCESS\' and (not row[\'reverse_valid\']) for row in bootstrap)), recorded_reverse_allocation_failure_N=sum((row[\'recorded_reverse_failure\'] for row in bootstrap)), reverse_only_failure_forward_retained_N=sum((row[\'forward_valid\'] and (not row[\'reverse_valid\']) for row in bootstrap)), attempt_failure_count_scope=attempt_scope, indexed_failed_attempts_N=len(failed_attempts), indexed_nonzero_process_exit_attempts_N=sum((exit_code(row) not in (None, 0) for row in indexed_attempts)), attempt_statistical_or_technical_class_inferred_from_exit_code=False, indexed_failed_attempt_status_counts=dict(Counter((row.get(\'status\') for row in failed_attempts))))\n\ndef write_csv(path, rows):\n    fields = list(dict.fromkeys((key for row in rows for key in row)))\n    if not fields:\n        fields, rows = ([\'status\'], [{\'status\': \'NO_SAVED_AGGREGATE_ROWS\'}])\n    with Path(path).open(\'x\', encoding=\'utf-8\', newline=\'\') as handle:\n        writer = csv.DictWriter(handle, fieldnames=fields)\n        writer.writeheader()\n        writer.writerows(rows)\n\ndef write_json(path, value_):\n    with Path(path).open(\'x\', encoding=\'utf-8\') as handle:\n        json.dump(value_, handle, ensure_ascii=False, indent=2, allow_nan=False)\n        handle.write(\'\\n\')\n\ndef release(args):\n    workspace = Path(args.workspace_root).resolve()\n    old_root = workspace / \'scripts/birth/border_selected44_stratified_20261103_v1\'\n    prior, prior_receipt = module([args.prior_reducer] if args.prior_reducer else [ROOT / \'vendor/prior121/tools/reduce.py\', old_root / \'vendor/prior121/tools/reduce.py\'], PRIOR_SHA, \'_gs100_mature_reducer\')\n    methods, methods_receipt = module([args.methods] if args.methods else [ROOT / \'vendor/prior_inference/methods.py\', old_root / \'vendor/prior_inference/methods.py\'], METHOD_SHA, \'_gs100_frozen_methods\')\n    b0 = Path(args.b0_result or workspace / \'reports/birth/GLYCAEMIA_SENSITIVITY_BASE11_B0_20261106_v1/remote_aggregate_evidence/RESULT.json\').resolve()\n    original_ref = Path(args.original_reference or workspace / \'scripts/birth/border_glycaemia_sensitivity_BASE11_B0_20261106_v1/tests/ORIGINAL_BASE11_B0_REFERENCE.json\').resolve()\n    review_root = Path(args.original_review_root or workspace / \'reports/birth/SELECTED44_ALL_RESULTS_REVIEW_20261105_v1/review_20261105T004052Z\').resolve()\n    output = Path(args.out).resolve()\n    require(not output.exists(), \'PRESERVE_EXISTING_REPORT_CHOOSE_FRESH_VERSIONED_OUTPUT_DIRECTORY\')\n    index, ledger, vectors, risks, pair_checks, mc, results = load_inputs(args.index, b0, prior)\n    old_effects, old_ratios, original_sources = frozen_original(original_ref, review_root)\n    index_binding = index.get(\'binding\') or {key: index[key] for key in (\'input_sha256\', \'config_sha256\', \'draw_manifest_sha256\', \'bootstrap_unit\') if key in index}\n    pseudo_manifest = {\'panel_id\': MODEL, \'domain_ids\': DOMAINS, \'scientific_binding\': index_binding}\n    totals, domain_rows, raw, sensitivity, inference_audit = prior.summarize(pseudo_manifest, ledger, vectors, methods, 1e-06)\n    effects = [normalize_effect(row) for row in totals if row[\'unit\'] != \'percent\']\n    domains = [normalize_effect(row) for row in domain_rows if row[\'unit\'] != \'percent\']\n    ratios = [normalize_ratio(row) for row in totals + domain_rows if row[\'unit\'] == \'percent\']\n    require(inference_audit[\'BH_family_n\'] == 5, \'FROZEN_EXACT5_BH_FAMILY_REQUIRED\')\n    original_comparison = compare_original(effects + domains, ratios, old_effects, old_ratios)\n    detailed = inference_details(effects + domains, vectors, methods)\n    population = []\n    for b, result in results.items():\n        n = result.get(\'N\', {})\n        population.append(dict(bootstrap_id=b, model_id=MODEL, source=\'ACTUAL_SAVED_REPLICATE\', N_target=n.get(\'N_target\'), N_fit=n.get(\'N_fit\', n.get(\'N_joint_fit\')), N_standardization=n.get(\'N_standardization\'), N_excluded=n.get(\'N_excluded\'), aggregate_population_and_events_json=json.dumps(n, ensure_ascii=False)))\n    tables = {\'RESULTS.csv\': effects, \'DOMAIN.csv\': domains, \'RATIOS.csv\': ratios, \'ORIGINAL_BASE11_COMPARISON.csv\': original_comparison, \'BOOTSTRAP_STATUS.csv\': ledger, \'RAW_EFFECT_AND_RATIO_VECTORS.csv\': raw, \'RAW_AGGREGATE_RISKS.csv\': risks, \'SAVED_PAIRED_VECTOR_CHECKS.csv\': pair_checks, \'FORWARD_REVERSE_SENSITIVITY.csv\': sensitivity, \'FROZEN_INFERENCE_DETAILS.csv\': detailed, \'POPULATION_FLOW.csv\': population, \'MC_DIAGNOSTICS.csv\': mc}\n    counts = count_status(index, ledger)\n    inference_audit[\'legacy_terminal_enum_metadata\'] = {key: inference_audit[key] for key in (\'fixed_ids_processed_to_terminal\', \'all_fixed_ids_terminal\')}\n    inference_audit.update(fixed_ids_processed_to_terminal=counts[\'terminal_task_N\'], all_fixed_ids_terminal=counts[\'all_fixed100_terminal\'], terminal_enum_extension=\'INPUT_PREPARATION_FAILED_TERMINAL_NO_EFFECT_VECTOR\')\n    audit = dict(schema=\'GLYCAEMIA_SENSITIVITY_BASE11_B100_SUMMARY_V1\', model_id=MODEL, source_index={\'path\': str(Path(args.index).resolve()), \'sha256\': sha(args.index)}, source_new_B0={\'path\': str(b0), \'sha256\': B0_SHA}, source_original_aggregate_review=original_sources, frozen_methods=methods_receipt, mature_reducer=prior_receipt, script_sha256=sha(__file__), fixed_B=100, inner_K=50, frozen_inference=inference_audit, failure_accounting=counts, bootstrap_point_estimator=\'NEW_SENSITIVITY_B0_NOT_BOOTSTRAP_MEAN\', ratios_clipped=False, denominator_abs_applied=False, descriptive_endpoint_P_missing_per_freeze=True, new_test_families=0, original_CI_P_q_recomputed=False, original_comparison_inference=\'DESCRIPTIVE_CHANGE_ONLY_NO_PAIRED_CROSS_COHORT_CI_OR_P\', publication_gate=\'ACTUAL_ALL_FIXED100_TERMINAL_NOT_SUCCESS_PERCENTAGE\', participant_reads=False, participant_rows=0, model_reads=False, remote_access=False, new_model_fits=0, new_integrations=0, numerical_pilot_tests=0)\n    output.mkdir(parents=True, exist_ok=False)\n    for name, rows in tables.items():\n        write_csv(output / name, rows)\n    rules_path = ROOT / \'docs/FROZEN_INFERENCE_RULES.md\'\n    with (output / \'METHODS.md\').open(\'x\', encoding=\'utf-8\') as handle:\n        handle.write(rules_path.read_text(encoding=\'utf-8\'))\n    fwd, rev = (inference_audit[\'forward_effect_valid_N\'], inference_audit[\'reverse_effect_valid_N\'])\n    report = f"# BASE11 HbA1c敏感性：固定B100聚合结果\\n\\n本目录仅对实际完成的固定B1–100终态记录归约。点估计固定取已完成的新敏感性B0，不是bootstrap均值；原BASE11结果未改动。\\n\\n固定100编号全部终态；正序两病共同完整有效B={fwd}，反序有效B={rev}。统计终态失败={counts[\'statistical_terminal_failure_N\']}；技术终态失败={counts[\'technical_terminal_failure_N\']}；输入准备终态失败={counts[\'input_preparation_terminal_failure_N\']}；结构终态失败={counts[\'structural_terminal_failure_N\']}；系统终态失败={counts[\'system_terminal_failure_N\']}；聚合传输/绑定校验拒绝={counts[\'aggregate_transport_or_binding_rejection_N\']}；仅反序失败而正序保留={counts[\'reverse_only_failure_forward_retained_N\']}。\\n\\nRESULTS.csv列示IHD、T2D及同编号配对Delta的TE/DE/IE_joint；DOMAIN.csv列示四域IIE和DEP。RATIOS.csv保留全部逐编号PM/Psep与分母稳定性。正反序各自使用自身完整编号；FORWARD_REVERSE_SENSITIVITY.csv仅用交集比较。\\n\\nCI严格沿用Type-7百分位。P仅按冻结规则用于正序Delta joint IE及四域+DEP；四域+DEP固定五项BH，joint IE单列。IHD/T2D、TE/DE、比例、反序P/q为空是冻结规则，并非漏算。P是实际有效B的离散CI反演筛查值，不是精确尾概率；FROZEN_INFERENCE_DETAILS.csv公开实际B、网格步长、搜索括号与限制。\\n\\n比例分母保留符号，负值或超过100%不截断。零/非有限只排除该比例，接近零、变号或CI跨零均标记；分母不稳定时有限比例分位数不代表可靠解释率区间。\\n\\nORIGINAL_BASE11_COMPARISON.csv包含原结果已冻结CI/P/q及新B0的绝对/有符号相对变化；比例绝对变化单位为百分点。两队列变化仅描述，不创造配对变化CI/P。未在冻结旧报告发布的指标以NOT_PUBLISHED标明，不回算旧区间。\\n\\n详见METHODS.md、SUMMARY_AUDIT.json及逐编号BOOTSTRAP_STATUS.csv。本聚合没有参与者级读取、模型拟合或新增积分。\\n"\n    with (output / \'REPORT_ZH.md\').open(\'x\', encoding=\'utf-8\') as handle:\n        handle.write(report)\n    audit[\'outputs\'] = {name: {\'sha256\': sha(output / name), \'scientific_rows\': len(rows)} for name, rows in tables.items()}\n    audit[\'outputs\'].update({name: {\'sha256\': sha(output / name)} for name in (\'METHODS.md\', \'REPORT_ZH.md\')})\n    write_json(output / \'SUMMARY_AUDIT.json\', audit)\n    write_json(output / \'OUTPUT_MANIFEST.json\', dict(schema=\'AGGREGATE_EXCLUSIVE_OUTPUT_MANIFEST_V1\', files={path.name: sha(path) for path in sorted(output.iterdir()) if path.is_file()}, old_results_unchanged_by_aggregate_tool=True))\n    return {\'output_root\': str(output), \'forward_valid_B\': fwd, \'reverse_valid_B\': rev, \'fixed100_all_terminal\': True, \'failure_accounting\': counts}\n\ndef main():\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--index\', default=os.environ.get(\'AGGREGATE_INDEX\'), required=not os.environ.get(\'AGGREGATE_INDEX\'))\n    parser.add_argument(\'--out\', \'--report-root\', dest=\'out\', default=os.environ.get(\'REPORT_OUTPUT_ROOT\'), required=not os.environ.get(\'REPORT_OUTPUT_ROOT\'))\n    parser.add_argument(\'--workspace-root\', default=os.environ.get(\'WORKSPACE_ROOT\', str(WORKSPACE)))\n    parser.add_argument(\'--b0-result\', \'--b0\', dest=\'b0_result\', default=os.environ.get(\'NEW_B0_RESULT\'))\n    parser.add_argument(\'--original-reference\', \'--reference\', dest=\'original_reference\', default=os.environ.get(\'ORIGINAL_B0_REFERENCE\'))\n    parser.add_argument(\'--original-review-root\', default=os.environ.get(\'ORIGINAL_REVIEW_ROOT\'))\n    parser.add_argument(\'--methods\', type=Path)\n    parser.add_argument(\'--prior-reducer\', type=Path)\n    args = parser.parse_args()\n    print(json.dumps({\'aggregate_input_root\': str(Path(args.index).resolve().parent), \'script_root\': str(ROOT), \'output_root\': str(Path(args.out).resolve()), \'workspace_root\': str(Path(args.workspace_root).resolve())}))\n    print(json.dumps(release(args), ensure_ascii=False, allow_nan=False))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
