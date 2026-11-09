"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse, concurrent.futures, csv, datetime, fcntl, hashlib, json, os, socket, subprocess, sys, time, uuid
from pathlib import Path
MODELS = ('S44_PROT_REF11', 'S44_PROT_PLUS15')

def read(p):
    return json.loads(Path(p).read_bytes())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def canonical(v):
    return json.dumps(v, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)

def digest(v):
    return hashlib.sha256(canonical(v).encode()).hexdigest()

def file_record(p):
    return dict(path=str(p), sha256=sha(p))

def check_record(v):
    p = Path(v['path'])
    assert p.is_file() and (not p.is_symlink()) and (sha(p) == v['sha256']), 'EXISTING_RESULT_BYTES_CHANGED'
    return p

def publish(p, v):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        assert not p.is_symlink() and read(p) == v, 'PRESERVE_EXISTING_DIFFERENT_METADATA'
        return
    with p.open('x') as f:
        json.dump(v, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.flush()
        os.fsync(f.fileno())

def text_file(p, text):
    with Path(p).open('x') as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())

def csv_file(p, rows):
    with Path(p).open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
        f.flush()
        os.fsync(f.fileno())

def publication(parent, prefix, key, writer):
    """Reuse an exact complete publication, otherwise publish a fresh directory.

    COMPLETE is fsynced inside an unpublished, unique directory. The final
    same-filesystem directory rename publishes all files and its complete
    manifest together. Interrupted staging directories are retained, not reused.
    """
    parent = Path(parent)
    parent.mkdir(parents=True, exist_ok=True)
    assert not parent.is_symlink()
    for folder in sorted(parent.glob(prefix + '*')):
        marker = folder / 'COMPLETE.json'
        if not marker.is_file():
            continue
        assert not folder.is_symlink() and (not marker.is_symlink()), 'PUBLICATION_SYMLINK_REFUSED'
        completed = read(marker)
        if completed.get('publication_key') != key:
            continue
        assert completed.get('schema') == 'SELECTED44_FINISH_V2_COMPLETE_V1' and completed.get('complete') is True
        for item in completed['files']:
            name = item['path']
            assert Path(name).name == name, 'PUBLICATION_RELATIVE_FILE_REQUIRED'
            check_record(dict(path=str(folder / name), sha256=item['sha256']))
        return folder
    token = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ') + '_' + uuid.uuid4().hex
    final = parent / (prefix + token)
    staging = parent / ('.publishing_' + token)
    assert not final.exists()
    staging.mkdir(exist_ok=False)
    writer(staging, final)
    files = [dict(path=p.name, sha256=sha(p)) for p in sorted(staging.iterdir()) if p.is_file()]
    assert files and (not (staging / 'COMPLETE.json').exists())
    publish(staging / 'COMPLETE.json', dict(schema='SELECTED44_FINISH_V2_COMPLETE_V1', complete=True, publication_key=key, files=files))
    assert not final.exists()
    os.rename(staging, final)
    return final

def failure_handoff(run, variant, reason, details, sensitivity=None):
    details = dict(details, reason=str(reason), status='NOT_ACCEPTED', variant=variant, sensitivity_publication=str(sensitivity) if sensitivity else None, new_draws=0, old_outputs_unchanged=True)

    def write(folder, final):
        publish(folder / 'HANDOFF_STATUS.json', details)
        text_file(folder / 'RELEASE_REPORT_ZH.md', '# SELECTED44 蛋白新版：尚未完成交接\n\n' + '主汇总失败、未完整或敏感性未完整，未授予计算验收。有效结果、旧结果及所有attempt均保留；未重算成功结果。\n\n' + '原因：' + str(reason) + '\n\n' + canonical(details) + '\n')
    folder = publication(run / variant / 'release_v2', 'FAILED_', digest(details), write)
    print(canonical(dict(ok=False, exit_code=2, failure_handoff=str(folder), reason=str(reason))), flush=True)
    return 2

def task_request(cfg, idx, variant, model, setting):
    return dict(schema='SELECTED44_FINISH_V2_B0_REQUEST_V1', model_id=model, b=0, setting=setting, variant=variant, source_manifest=idx['source_manifest'], source_adoption=idx.get('source_adoption'), all_y_gate=idx['all_y_gate'], main_B0_result=idx['jobs']['MC'][model + ':0']['result'], runtime_sha256=sha(CODE / 'ops/runtime.json'), finisher_sha256=sha(Path(__file__)), worker_sha256=sha(CODE / 'R/worker.R'), new_bootstrap_draws=0)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', choices=('primary', 'body_ridge'), required=True)
    parser.add_argument('--wait', action='store_true')
    args = parser.parse_args()
    cfg = read(CODE / 'ops/runtime.json')
    state = Path(cfg['state_root'])
    run = Path(cfg['run_root'])
    assert socket.gethostname() == cfg['hostname'] and os.uname().sysname == 'Darwin'
    assert Path(cfg['code_root']).resolve() == CODE and 'protein_ridge_b100_20261104_v1' in str(state) and ('protein_ridge_b100_20261104_v1' in str(run))
    with (state / 'protein_ridge_finish.lock').open('a+b') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        index_path = state / 'controller' / args.variant / 'INDEX.json'
        while True:
            idx = read(index_path)
            if all((j['state'] == 'SUCCESS' for j in idx['jobs']['MC'].values())):
                break
            if idx['stage'] in ('FAILED_AWAIT_REVIEW', 'STOPPED'):
                raise RuntimeError('MAIN_REPAIR_NOT_COMPLETE:' + idx['stage'])
            if not args.wait:
                raise RuntimeError('MAIN_MC_NOT_YET_COMPLETE')
            time.sleep(20)
        pin = read(check_record(idx['source_manifest']))
        assert all((sha(CODE / p) == s for p, s in pin['files'].items()))
        assert read(check_record(idx['all_y_gate']))['passed']
        target = run / args.variant / 'B0_PENALTY_SENSITIVITY_V2'
        target.mkdir(parents=True, exist_ok=True)
        assert not target.is_symlink(), 'VERSIONED_SENSITIVITY_TARGET_SYMLINK_REFUSED'
        settings = ['ethnicity_low', 'ethnicity_high']
        if args.variant == 'body_ridge':
            settings += ['body_low', 'body_high']
        pending = [(m, s) for m in MODELS for s in settings]

        def task(pair):
            model, setting = pair
            parent = target / 'jobs' / model / setting
            parent.mkdir(parents=True, exist_ok=True)
            assert not parent.is_symlink(), 'SENSITIVITY_JOB_SYMLINK_REFUSED'
            request = task_request(cfg, idx, args.variant, model, setting)
            request_sha = digest(request)
            request_path = parent / 'REQUEST.json'
            publish(request_path, request)
            receipt = parent / 'COMPLETED.json'
            if receipt.exists():
                old = read(receipt)
                assert old['request_sha256'] == request_sha, 'EXACT_SENSITIVITY_REQUEST_REQUIRED'
                assert old['model_id'] == model and old['setting'] == setting and (old['b'] == 0) and (old['status'] == 'SUCCESS')
                check_record(dict(path=old['result_path'], sha256=old['result_sha256']))
                assert len(old['phase_exits']) == 2, 'BOTH_SAVED_PHASE_EXITS_REQUIRED'
                for proof in old['phase_exits']:
                    exit_record = read(check_record(proof))
                    check_record(exit_record['result'])
                    assert exit_record['exit_code'] == 0 and exit_record['request_sha256'] == request_sha
                return old
            claim = parent / 'ATTEMPT.json'
            if claim.exists():
                saved = read(claim)
                assert saved['request_sha256'] == request_sha, 'EXISTING_ATTEMPT_REQUEST_CHANGED'
                token = saved['token']
                attempt = parent / ('attempt_' + token)
                assert attempt.is_dir() and (not attempt.is_symlink()), 'EXISTING_ATTEMPT_MISSING_OR_SYMLINK'
            else:
                token = uuid.uuid4().hex
                attempt = parent / ('attempt_' + token)
                attempt.mkdir()
                publish(claim, dict(token=token, model_id=model, b=0, setting=setting, request_sha256=request_sha))
            env = dict(os.environ, **cfg['environment'], PRY_B0_SENSITIVITY=setting)
            phase_exits = []
            for phase in ('Y', 'MC'):
                out = attempt / phase
                log = attempt / (phase + '.log')
                if phase == 'MC':
                    yr = read(attempt / 'Y/result.json')
                    yp = attempt / 'Y' / yr['Y']['path']
                    assert yr['status'] == 'SUCCESS' and sha(yp) == yr['Y']['sha256']
                    env.update(PRY_B0_SENSITIVITY_Y=str(yp), PRY_B0_SENSITIVITY_Y_SHA256=yr['Y']['sha256'])
                assert all((sha(CODE / p) == s for p, s in pin['files'].items()))
                argv = [cfg['rscript'], '--vanilla', str(CODE / 'R/worker.R'), str(CODE), phase, model, '0', args.variant, str(out)]
                phase_request = attempt / (phase + '_REQUEST.json')
                exit_path = attempt / (phase + '_EXIT.json')
                intended = dict(request_sha256=request_sha, phase=phase, argv=argv, out=str(out), sensitivity_environment={k: v for k, v in env.items() if k.startswith('PRY_B0_SENSITIVITY')})
                if exit_path.exists():
                    observed = read(exit_path)
                    assert read(phase_request) == intended
                    assert observed['phase_request_sha256'] == sha(phase_request) and observed['request_sha256'] == request_sha
                    if observed.get('result'):
                        check_record(observed['result'])
                else:
                    assert not phase_request.exists(), 'EXISTING_PHASE_WITHOUT_EXIT_REVIEW_REQUIRED_NO_RETRY'
                    assert not out.exists() and (not log.exists()), 'EXISTING_PHASE_OUTPUT_PRESERVED'
                    publish(phase_request, intended)
                    started = time.time()
                    with log.open('xb') as f:
                        r = subprocess.run(argv, stdin=subprocess.DEVNULL, stdout=f, stderr=subprocess.STDOUT, env=env, cwd=CODE)
                    result_path = out / 'result.json'
                    observed = dict(exit_code=r.returncode, request_sha256=request_sha, phase_request_sha256=sha(phase_request), started_unix=started, finished_unix=time.time(), result=file_record(result_path) if result_path.is_file() else None)
                    publish(exit_path, observed)
                phase_exits.append(file_record(exit_path))
                result = read(out / 'result.json') if observed.get('result') else {}
                if observed['exit_code'] or result.get('status') != 'SUCCESS':
                    failed = dict(model_id=model, b=0, setting=setting, status='FAILED', phase=phase, reason=result.get('reason', 'PROCESS_FAILED_OR_RESULT_MISSING'), attempt=str(attempt), request_sha256=request_sha, phase_exits=phase_exits)
                    publish(parent / 'FAILED.json', failed)
                    return failed
            completed = dict(model_id=model, b=0, setting=setting, status='SUCCESS', result_path=str(attempt / 'MC/result.json'), result_sha256=sha(attempt / 'MC/result.json'), request_sha256=request_sha, phase_exits=phase_exits, new_bootstrap_draws=0)
            publish(receipt, completed)
            return completed
        slots = min(4, idx['calibration']['MC']['max_workers'])
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=slots) as pool:
            futures = [(pair, pool.submit(task, pair)) for pair in pending]
            for pair, future in futures:
                try:
                    results.append(future.result())
                except Exception as error:
                    results.append(dict(model_id=pair[0], b=0, setting=pair[1], status='FAILED', reason='SAVED_ATTEMPT_PRESERVED_REVIEW_REQUIRED:' + str(error), new_bootstrap_draws=0))
        rows = []
        failures = []

        def keyed(z):
            return {(r['allocation_order'], r['endpoint'], r['component_type'], r.get('domain_id')): r['effect'] for r in z['components']}
        for record in results:
            if record['status'] != 'SUCCESS':
                failures.append(record)
                continue
            model = record['model_id']
            point = read(idx['jobs']['MC'][model + ':0']['result']['path'])
            actual = read(record['result_path'])
            assert actual['binding'] == point['binding'] and actual['pair_binding'] == point['pair_binding']
            base = keyed(point)
            sens = keyed(actual)
            assert set(base) == set(sens)
            for order, e, kind, domain in base:
                rows.append(dict(model_id=model, setting=record['setting'], allocation_order=order, endpoint=e, component_type=kind, domain_id=domain, main_per_1000=1000 * base[order, e, kind, domain], sensitivity_per_1000=1000 * sens[order, e, kind, domain], change_per_1000=1000 * (sens[order, e, kind, domain] - base[order, e, kind, domain]), main_PM_percent=100 * base[order, e, kind, domain] / base[order, e, 'TE', None] if base[order, e, 'TE', None] else None, sensitivity_PM_percent=100 * sens[order, e, kind, domain] / sens[order, e, 'TE', None] if sens[order, e, 'TE', None] else None))
            for order, kind, domain in {(o, k, d) for o, e, k, d in base}:
                a = base[order, 'IHD', kind, domain] - base[order, 'T2D', kind, domain]
                v = sens[order, 'IHD', kind, domain] - sens[order, 'T2D', kind, domain]
                den = sens[order, 'IHD', 'TE', None] - sens[order, 'T2D', 'TE', None]
                bden = base[order, 'IHD', 'TE', None] - base[order, 'T2D', 'TE', None]
                rows.append(dict(model_id=model, setting=record['setting'], allocation_order=order, endpoint='Delta', component_type=kind, domain_id=domain, main_per_1000=1000 * a, sensitivity_per_1000=1000 * v, change_per_1000=1000 * (v - a), main_PM_percent=100 * a / bden if bden else None, sensitivity_PM_percent=100 * v / den if den else None))
        evidence = dict(passed=not failures, settings=settings, models=MODELS, results=results, failures=failures, source_manifest_sha256=idx['source_manifest']['sha256'], fixed_bootstrap_draws_added=0, interpretation='Prespecified B0 model-dependence analysis, not parameter selection')
        sensitivity_key = digest(evidence)

        def write_sensitivity(folder, final):
            if rows:
                csv_file(folder / 'B0_PENALTY_SENSITIVITY.csv', rows)
            publish(folder / 'SENSITIVITY_ACCEPTANCE.json', dict(evidence, comparison_csv=str(final / 'B0_PENALTY_SENSITIVITY.csv') if rows else None))
            text_file(folder / 'REPORT_ZH.md', '# 原样本惩罚敏感性\n\n' + f"完成{sum((r['status'] == 'SUCCESS' for r in results))}/{len(results)}个B0设置；未新增Bootstrap。\n\n" + '主规格始终df=k/2；附加设置为k/3与2k/3，一次只改变一块。A、抽样、人群和MC随机流不变。数值差异见配套表；不按IE方向或显著性选惩罚。\n')
        sensitivity = publication(target / 'publications', 'PUBLISH_', sensitivity_key, write_sensitivity)
        csv_path = sensitivity / 'B0_PENALTY_SENSITIVITY.csv'
        while True:
            latest = read(index_path)
            summary = latest.get('summary', {})
            if summary.get('state') in ('SUCCESS', 'FAILED_REVIEW'):
                break
            if latest.get('stage') in ('FAILED_AWAIT_REVIEW', 'STOPPED', 'SUMMARY_FAILED_REVIEW'):
                return failure_handoff(run, args.variant, '主汇总尚未完成且控制器已停止', dict(controller_stage=latest.get('stage'), summary_state=summary.get('state')), sensitivity)
            time.sleep(5)
        main_dir = None
        if summary.get('folder'):
            log = Path(summary['folder']) / 'summary.log'
            if log.is_file():
                for line in reversed(log.read_text().splitlines()):
                    try:
                        doc = json.loads(line)
                    except ValueError:
                        continue
                    if isinstance(doc, dict) and doc.get('summary_dir'):
                        main_dir = Path(doc['summary_dir'])
                        break
        failure_details = dict(summary_state=summary.get('state'), summary_dir=str(main_dir) if main_dir else None, sensitivity_failures=failures)
        if summary.get('state') != 'SUCCESS':
            return failure_handoff(run, args.variant, '主汇总执行失败；未读取正式结果表', failure_details, sensitivity)
        if main_dir is None:
            return failure_handoff(run, args.variant, '主汇总未提供完整发布目录', failure_details, sensitivity)
        try:
            assert main_dir.absolute().is_relative_to(run.absolute()) and (not main_dir.is_symlink()), 'SUMMARY_PATH_OUTSIDE_REPAIR'
            main_audit = read(main_dir / 'SUMMARY_AUDIT.json')
            assert main_audit.get('computation_passed') is True and main_audit.get('full_100_paired_valid') is True, 'MAIN_202_MC_OR_100_PAIRED_VALID_NOT_ACCEPTED'
            assert main_audit.get('old_B0_joint_comparison_complete') is True, 'OLD_B0_COMPARISON_INCOMPLETE'
            assert main_audit.get('source_manifest') == idx['source_manifest'], 'MAIN_SUMMARY_SOURCE_ADOPTION_MISMATCH'
            assert main_audit.get('BH_family_n') == {MODELS[0]: 5, MODELS[1]: 6} and main_audit.get('joint_Delta_IE_separate') is True, 'MAIN_FROZEN_INFERENCE_CONTRACT_MISMATCH'
            if failures:
                raise RuntimeError('B0_SENSITIVITY_INCOMPLETE_OR_FAILED_RESULTS_PRESERVED')
            for name in ('RESULTS', 'DOMAIN', 'RATIOS', 'OLD_B0_COMPARISON'):
                assert (main_dir / (name + '.csv')).is_file() and (not (main_dir / (name + '.csv')).is_symlink()), 'MISSING_MAIN_TABLE:' + name

            def table(name):
                with (main_dir / (name + '.csv')).open(newline='') as f:
                    return list(csv.DictReader(f))
            effects = table('RESULTS')
            domain = table('DOMAIN')
            ratios = table('RATIOS')
            assert effects and domain and ratios, 'EMPTY_MAIN_RESULT_TABLE'
        except Exception as error:
            return failure_handoff(run, args.variant, '主汇总或必要附加结果不完整：' + str(error), failure_details, sensitivity)
        complete = True
        accepted = dict(schema='SELECTED44_PROTEIN_RIDGE_FINAL_COMPUTATION_ACCEPTANCE_V1', computation_accepted=True, full_new_method_B0_and_B100=True, full_100_true_paired_draws=True, primary_method=args.variant, body_response_penalty_used=args.variant == 'body_ridge', sensitivity_completed=True, sensitivity_evidence=str(sensitivity / 'SENSITIVITY_ACCEPTANCE.json'), sensitivity_evidence_sha256=sha(sensitivity / 'SENSITIVITY_ACCEPTANCE.json'), main_summary=str(main_dir), main_audit_sha256=sha(main_dir / 'SUMMARY_AUDIT.json'), new_draws=0, B101_started=False, B1000_started=False, note='Computational acceptance is not proof of causal identification, unbiasedness, coverage or selection-wide FDR')

        def fmt(x):
            if x in (None, ''):
                return 'NA'
            try:
                return f'{float(x):.6g}'
            except (ValueError, TypeError):
                return str(x)
        text = ['# SELECTED44 蛋白联合中介：新版固定B100最终交付', '', '计算验收：' + ('通过。' if complete else '尚未通过；有效产物均保留，见验收JSON。'), '', '范围仅两个蛋白队列模型，原样本9,174人，IHD/T2D、全随访拟合、11年风险、K50及原固定家系B1–100。没有重抽样或启动B101/B1000。', '', '四个原因别分支统一中心分层；原族裔编码使用原生ridge，主规则为各次拟合df=k/2、Efron、选择权重及robust=FALSE。', '使用统一体成分响应岭惩罚后备，不能视为纯数值修复。' if args.variant == 'body_ridge' else '未启用体成分响应惩罚；该块仅做等价缩放和全列可逆QR。', '兼容A及实际索引复用，旧Y与旧成功重复不混入新区间。', '', '## 联合效应', '', '单位为每1,000人的11年绝对风险差。区间来自实际共同有效编号的Type7百分位。', '', '| 模型 | 结局 | TE | IE | IE 95%CI | 有效重复 |', '|---|---|---:|---:|---|---:|']
        for m in MODELS:
            for e in ('IHD', 'T2D', 'Delta'):
                selected = {r['component_type']: r for r in effects if r['model_id'] == m and r['allocation_order'] == 'forward' and (r['endpoint'] == e)}
                if 'IE_joint' not in selected:
                    continue
                ie = selected['IE_joint']
                te = selected.get('TE', {})
                text.append('| ' + ' | '.join((m, e, fmt(te.get('point_B0_per_1000')), fmt(ie.get('point_B0_per_1000')), '[' + fmt(ie.get('CI_type7_low_per_1000')) + ', ' + fmt(ie.get('CI_type7_high_per_1000')) + ']', fmt(ie.get('effect_valid_N')))) + ' |')
        text += ['', '## 正式蛋白域贡献', '', '这来自PLUS15的正式域分解，不是PLUS−REF联合IE差。', '', '| 结局 | 蛋白域IIE/1,000 | 95%CI | P | q（仅前向Δ，BH6） |', '|---|---:|---|---:|---:|']
        for r in domain:
            if r['model_id'] == MODELS[1] and r['allocation_order'] == 'forward' and (r['component_type'] == 'IIE') and (r['domain_id'] == 'CANDIDATE_PROTEINS'):
                text.append('| ' + ' | '.join((r['endpoint'], fmt(r.get('point_B0_per_1000')), '[' + fmt(r.get('CI_type7_low_per_1000')) + ', ' + fmt(r.get('CI_type7_high_per_1000')) + ']', fmt(r.get('P')), fmt(r.get('q')))) + ' |')
        text += ['', '## 惩罚敏感性及限制', '', f'B0预指定附加设置完成{len(results) - len(failures)}/{len(results)}；同一A、样本、权重及MC随机流，一次只改一个惩罚块。', '实际逐项变化见 ' + str(csv_path) + '。敏感性不带未经计算的Bootstrap区间，也不用于挑选更显著的设置。']
        relevant = [r for r in rows if r['allocation_order'] == 'forward' and (r['component_type'] == 'IE_joint' or r.get('domain_id') == 'CANDIDATE_PROTEINS')]
        flips = [r for r in relevant if r['main_per_1000'] * r['sensitivity_per_1000'] < 0]
        text.append('上述关键前向点估计出现符号改变的敏感性记录数：' + str(len(flips)) + '；最大绝对变化/1,000：' + fmt(max((abs(r['change_per_1000']) for r in relevant), default=None)) + '。小的近零点估计翻号不自动等于重要结论改变，需结合主区间及绝对量解释。')
        unstable = sum((str(r.get('denominator_unstable', '')).lower() == 'true' for r in ratios))
        text += ['', f'比值分母不稳定记录数：{unstable}。PM/Psep不裁剪；近零或跨零TE及Fieller无界/不连通集合均在原表保留。', 'P、q沿用冻结的B100 Type7区间反演及模型内域差值+DEP家族（REF BH5、PLUS BH6）；联合Δ单列，不套用BH18。B100的尾部与P/q精度有限，不是独立确认或筛选全流程FDR保证。', 'ridge的df约定不是最优参数或无偏保证；计算通过并不证明所有中介可交换性、模型正确性或竞争风险因果假设。', '', '## 文件', '', str(main_dir), str(csv_path)]

        def write_release(folder, final):
            publish(folder / 'RELEASE_ACCEPTANCE.json', accepted)
            text_file(folder / 'RELEASE_REPORT_ZH.md', '\n'.join(text + ['', str(final / 'RELEASE_ACCEPTANCE.json')]) + '\n')
            publish(folder / 'OUTPUT_MANIFEST.json', dict(files=[dict(path=str(final / p.name), sha256=sha(p)) for p in folder.iterdir() if p.is_file()], main_summary=str(main_dir), sensitivity=str(sensitivity), actual_results_only=True))
        final = publication(run / args.variant / 'release_v2', 'RELEASE_', digest(accepted), write_release)
        print(canonical(dict(ok=True, computation_accepted=True, release=str(final), sensitivity=str(sensitivity))))
        return 0

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Restart-safe Mini-only B0 penalty sensitivity and final delivery; no new draws.\n\nWaits for this repair\'s 202 main MC attempts only, then runs the preauthorized\nB0 settings, one block at a time, with unchanged worker/A/streams. Never retries\nstatistical failures, changes a model, or touches other production queues.\n"""\nimport argparse, concurrent.futures, csv, datetime, fcntl, hashlib, json, os, socket, subprocess, sys, time, uuid\nfrom pathlib import Path\nCODE = Path(__file__).resolve().parents[1]\nMODELS = (\'S44_PROT_REF11\', \'S44_PROT_PLUS15\')\n\ndef read(p):\n    return json.loads(Path(p).read_bytes())\n\ndef sha(p):\n    return hashlib.sha256(Path(p).read_bytes()).hexdigest()\n\ndef canonical(v):\n    return json.dumps(v, ensure_ascii=False, sort_keys=True, separators=(\',\', \':\'), allow_nan=False)\n\ndef digest(v):\n    return hashlib.sha256(canonical(v).encode()).hexdigest()\n\ndef file_record(p):\n    return dict(path=str(p), sha256=sha(p))\n\ndef check_record(v):\n    p = Path(v[\'path\'])\n    assert p.is_file() and (not p.is_symlink()) and (sha(p) == v[\'sha256\']), \'EXISTING_RESULT_BYTES_CHANGED\'\n    return p\n\ndef publish(p, v):\n    p = Path(p)\n    p.parent.mkdir(parents=True, exist_ok=True)\n    if p.exists():\n        assert not p.is_symlink() and read(p) == v, \'PRESERVE_EXISTING_DIFFERENT_METADATA\'\n        return\n    with p.open(\'x\') as f:\n        json.dump(v, f, ensure_ascii=False, indent=2, allow_nan=False)\n        f.flush()\n        os.fsync(f.fileno())\n\ndef text_file(p, text):\n    with Path(p).open(\'x\') as f:\n        f.write(text)\n        f.flush()\n        os.fsync(f.fileno())\n\ndef csv_file(p, rows):\n    with Path(p).open(\'x\', newline=\'\') as f:\n        w = csv.DictWriter(f, fieldnames=list(rows[0]))\n        w.writeheader()\n        w.writerows(rows)\n        f.flush()\n        os.fsync(f.fileno())\n\ndef publication(parent, prefix, key, writer):\n    """Reuse an exact complete publication, otherwise publish a fresh directory.\n\n    COMPLETE is fsynced inside an unpublished, unique directory. The final\n    same-filesystem directory rename publishes all files and its complete\n    manifest together. Interrupted staging directories are retained, not reused.\n    """\n    parent = Path(parent)\n    parent.mkdir(parents=True, exist_ok=True)\n    assert not parent.is_symlink()\n    for folder in sorted(parent.glob(prefix + \'*\')):\n        marker = folder / \'COMPLETE.json\'\n        if not marker.is_file():\n            continue\n        assert not folder.is_symlink() and (not marker.is_symlink()), \'PUBLICATION_SYMLINK_REFUSED\'\n        completed = read(marker)\n        if completed.get(\'publication_key\') != key:\n            continue\n        assert completed.get(\'schema\') == \'SELECTED44_FINISH_V2_COMPLETE_V1\' and completed.get(\'complete\') is True\n        for item in completed[\'files\']:\n            name = item[\'path\']\n            assert Path(name).name == name, \'PUBLICATION_RELATIVE_FILE_REQUIRED\'\n            check_record(dict(path=str(folder / name), sha256=item[\'sha256\']))\n        return folder\n    token = datetime.datetime.now(datetime.timezone.utc).strftime(\'%Y%m%dT%H%M%S%fZ\') + \'_\' + uuid.uuid4().hex\n    final = parent / (prefix + token)\n    staging = parent / (\'.publishing_\' + token)\n    assert not final.exists()\n    staging.mkdir(exist_ok=False)\n    writer(staging, final)\n    files = [dict(path=p.name, sha256=sha(p)) for p in sorted(staging.iterdir()) if p.is_file()]\n    assert files and (not (staging / \'COMPLETE.json\').exists())\n    publish(staging / \'COMPLETE.json\', dict(schema=\'SELECTED44_FINISH_V2_COMPLETE_V1\', complete=True, publication_key=key, files=files))\n    assert not final.exists()\n    os.rename(staging, final)\n    return final\n\ndef failure_handoff(run, variant, reason, details, sensitivity=None):\n    details = dict(details, reason=str(reason), status=\'NOT_ACCEPTED\', variant=variant, sensitivity_publication=str(sensitivity) if sensitivity else None, new_draws=0, old_outputs_unchanged=True)\n\n    def write(folder, final):\n        publish(folder / \'HANDOFF_STATUS.json\', details)\n        text_file(folder / \'RELEASE_REPORT_ZH.md\', \'# SELECTED44 蛋白新版：尚未完成交接\\n\\n\' + \'主汇总失败、未完整或敏感性未完整，未授予计算验收。有效结果、旧结果及所有attempt均保留；未重算成功结果。\\n\\n\' + \'原因：\' + str(reason) + \'\\n\\n\' + canonical(details) + \'\\n\')\n    folder = publication(run / variant / \'release_v2\', \'FAILED_\', digest(details), write)\n    print(canonical(dict(ok=False, exit_code=2, failure_handoff=str(folder), reason=str(reason))), flush=True)\n    return 2\n\ndef task_request(cfg, idx, variant, model, setting):\n    return dict(schema=\'SELECTED44_FINISH_V2_B0_REQUEST_V1\', model_id=model, b=0, setting=setting, variant=variant, source_manifest=idx[\'source_manifest\'], source_adoption=idx.get(\'source_adoption\'), all_y_gate=idx[\'all_y_gate\'], main_B0_result=idx[\'jobs\'][\'MC\'][model + \':0\'][\'result\'], runtime_sha256=sha(CODE / \'ops/runtime.json\'), finisher_sha256=sha(Path(__file__)), worker_sha256=sha(CODE / \'R/worker.R\'), new_bootstrap_draws=0)\n\ndef main():\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--variant\', choices=(\'primary\', \'body_ridge\'), required=True)\n    parser.add_argument(\'--wait\', action=\'store_true\')\n    args = parser.parse_args()\n    cfg = read(CODE / \'ops/runtime.json\')\n    state = Path(cfg[\'state_root\'])\n    run = Path(cfg[\'run_root\'])\n    assert socket.gethostname() == cfg[\'hostname\'] and os.uname().sysname == \'Darwin\'\n    assert Path(cfg[\'code_root\']).resolve() == CODE and \'protein_ridge_b100_20261104_v1\' in str(state) and (\'protein_ridge_b100_20261104_v1\' in str(run))\n    with (state / \'protein_ridge_finish.lock\').open(\'a+b\') as lock:\n        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)\n        index_path = state / \'controller\' / args.variant / \'INDEX.json\'\n        while True:\n            idx = read(index_path)\n            if all((j[\'state\'] == \'SUCCESS\' for j in idx[\'jobs\'][\'MC\'].values())):\n                break\n            if idx[\'stage\'] in (\'FAILED_AWAIT_REVIEW\', \'STOPPED\'):\n                raise RuntimeError(\'MAIN_REPAIR_NOT_COMPLETE:\' + idx[\'stage\'])\n            if not args.wait:\n                raise RuntimeError(\'MAIN_MC_NOT_YET_COMPLETE\')\n            time.sleep(20)\n        pin = read(check_record(idx[\'source_manifest\']))\n        assert all((sha(CODE / p) == s for p, s in pin[\'files\'].items()))\n        assert read(check_record(idx[\'all_y_gate\']))[\'passed\']\n        target = run / args.variant / \'B0_PENALTY_SENSITIVITY_V2\'\n        target.mkdir(parents=True, exist_ok=True)\n        assert not target.is_symlink(), \'VERSIONED_SENSITIVITY_TARGET_SYMLINK_REFUSED\'\n        settings = [\'ethnicity_low\', \'ethnicity_high\']\n        if args.variant == \'body_ridge\':\n            settings += [\'body_low\', \'body_high\']\n        pending = [(m, s) for m in MODELS for s in settings]\n\n        def task(pair):\n            model, setting = pair\n            parent = target / \'jobs\' / model / setting\n            parent.mkdir(parents=True, exist_ok=True)\n            assert not parent.is_symlink(), \'SENSITIVITY_JOB_SYMLINK_REFUSED\'\n            request = task_request(cfg, idx, args.variant, model, setting)\n            request_sha = digest(request)\n            request_path = parent / \'REQUEST.json\'\n            publish(request_path, request)\n            receipt = parent / \'COMPLETED.json\'\n            if receipt.exists():\n                old = read(receipt)\n                assert old[\'request_sha256\'] == request_sha, \'EXACT_SENSITIVITY_REQUEST_REQUIRED\'\n                assert old[\'model_id\'] == model and old[\'setting\'] == setting and (old[\'b\'] == 0) and (old[\'status\'] == \'SUCCESS\')\n                check_record(dict(path=old[\'result_path\'], sha256=old[\'result_sha256\']))\n                assert len(old[\'phase_exits\']) == 2, \'BOTH_SAVED_PHASE_EXITS_REQUIRED\'\n                for proof in old[\'phase_exits\']:\n                    exit_record = read(check_record(proof))\n                    check_record(exit_record[\'result\'])\n                    assert exit_record[\'exit_code\'] == 0 and exit_record[\'request_sha256\'] == request_sha\n                return old\n            claim = parent / \'ATTEMPT.json\'\n            if claim.exists():\n                saved = read(claim)\n                assert saved[\'request_sha256\'] == request_sha, \'EXISTING_ATTEMPT_REQUEST_CHANGED\'\n                token = saved[\'token\']\n                attempt = parent / (\'attempt_\' + token)\n                assert attempt.is_dir() and (not attempt.is_symlink()), \'EXISTING_ATTEMPT_MISSING_OR_SYMLINK\'\n            else:\n                token = uuid.uuid4().hex\n                attempt = parent / (\'attempt_\' + token)\n                attempt.mkdir()\n                publish(claim, dict(token=token, model_id=model, b=0, setting=setting, request_sha256=request_sha))\n            env = dict(os.environ, **cfg[\'environment\'], PRY_B0_SENSITIVITY=setting)\n            phase_exits = []\n            for phase in (\'Y\', \'MC\'):\n                out = attempt / phase\n                log = attempt / (phase + \'.log\')\n                if phase == \'MC\':\n                    yr = read(attempt / \'Y/result.json\')\n                    yp = attempt / \'Y\' / yr[\'Y\'][\'path\']\n                    assert yr[\'status\'] == \'SUCCESS\' and sha(yp) == yr[\'Y\'][\'sha256\']\n                    env.update(PRY_B0_SENSITIVITY_Y=str(yp), PRY_B0_SENSITIVITY_Y_SHA256=yr[\'Y\'][\'sha256\'])\n                assert all((sha(CODE / p) == s for p, s in pin[\'files\'].items()))\n                argv = [cfg[\'rscript\'], \'--vanilla\', str(CODE / \'R/worker.R\'), str(CODE), phase, model, \'0\', args.variant, str(out)]\n                phase_request = attempt / (phase + \'_REQUEST.json\')\n                exit_path = attempt / (phase + \'_EXIT.json\')\n                intended = dict(request_sha256=request_sha, phase=phase, argv=argv, out=str(out), sensitivity_environment={k: v for k, v in env.items() if k.startswith(\'PRY_B0_SENSITIVITY\')})\n                if exit_path.exists():\n                    observed = read(exit_path)\n                    assert read(phase_request) == intended\n                    assert observed[\'phase_request_sha256\'] == sha(phase_request) and observed[\'request_sha256\'] == request_sha\n                    if observed.get(\'result\'):\n                        check_record(observed[\'result\'])\n                else:\n                    assert not phase_request.exists(), \'EXISTING_PHASE_WITHOUT_EXIT_REVIEW_REQUIRED_NO_RETRY\'\n                    assert not out.exists() and (not log.exists()), \'EXISTING_PHASE_OUTPUT_PRESERVED\'\n                    publish(phase_request, intended)\n                    started = time.time()\n                    with log.open(\'xb\') as f:\n                        r = subprocess.run(argv, stdin=subprocess.DEVNULL, stdout=f, stderr=subprocess.STDOUT, env=env, cwd=CODE)\n                    result_path = out / \'result.json\'\n                    observed = dict(exit_code=r.returncode, request_sha256=request_sha, phase_request_sha256=sha(phase_request), started_unix=started, finished_unix=time.time(), result=file_record(result_path) if result_path.is_file() else None)\n                    publish(exit_path, observed)\n                phase_exits.append(file_record(exit_path))\n                result = read(out / \'result.json\') if observed.get(\'result\') else {}\n                if observed[\'exit_code\'] or result.get(\'status\') != \'SUCCESS\':\n                    failed = dict(model_id=model, b=0, setting=setting, status=\'FAILED\', phase=phase, reason=result.get(\'reason\', \'PROCESS_FAILED_OR_RESULT_MISSING\'), attempt=str(attempt), request_sha256=request_sha, phase_exits=phase_exits)\n                    publish(parent / \'FAILED.json\', failed)\n                    return failed\n            completed = dict(model_id=model, b=0, setting=setting, status=\'SUCCESS\', result_path=str(attempt / \'MC/result.json\'), result_sha256=sha(attempt / \'MC/result.json\'), request_sha256=request_sha, phase_exits=phase_exits, new_bootstrap_draws=0)\n            publish(receipt, completed)\n            return completed\n        slots = min(4, idx[\'calibration\'][\'MC\'][\'max_workers\'])\n        results = []\n        with concurrent.futures.ThreadPoolExecutor(max_workers=slots) as pool:\n            futures = [(pair, pool.submit(task, pair)) for pair in pending]\n            for pair, future in futures:\n                try:\n                    results.append(future.result())\n                except Exception as error:\n                    results.append(dict(model_id=pair[0], b=0, setting=pair[1], status=\'FAILED\', reason=\'SAVED_ATTEMPT_PRESERVED_REVIEW_REQUIRED:\' + str(error), new_bootstrap_draws=0))\n        rows = []\n        failures = []\n\n        def keyed(z):\n            return {(r[\'allocation_order\'], r[\'endpoint\'], r[\'component_type\'], r.get(\'domain_id\')): r[\'effect\'] for r in z[\'components\']}\n        for record in results:\n            if record[\'status\'] != \'SUCCESS\':\n                failures.append(record)\n                continue\n            model = record[\'model_id\']\n            point = read(idx[\'jobs\'][\'MC\'][model + \':0\'][\'result\'][\'path\'])\n            actual = read(record[\'result_path\'])\n            assert actual[\'binding\'] == point[\'binding\'] and actual[\'pair_binding\'] == point[\'pair_binding\']\n            base = keyed(point)\n            sens = keyed(actual)\n            assert set(base) == set(sens)\n            for order, e, kind, domain in base:\n                rows.append(dict(model_id=model, setting=record[\'setting\'], allocation_order=order, endpoint=e, component_type=kind, domain_id=domain, main_per_1000=1000 * base[order, e, kind, domain], sensitivity_per_1000=1000 * sens[order, e, kind, domain], change_per_1000=1000 * (sens[order, e, kind, domain] - base[order, e, kind, domain]), main_PM_percent=100 * base[order, e, kind, domain] / base[order, e, \'TE\', None] if base[order, e, \'TE\', None] else None, sensitivity_PM_percent=100 * sens[order, e, kind, domain] / sens[order, e, \'TE\', None] if sens[order, e, \'TE\', None] else None))\n            for order, kind, domain in {(o, k, d) for o, e, k, d in base}:\n                a = base[order, \'IHD\', kind, domain] - base[order, \'T2D\', kind, domain]\n                v = sens[order, \'IHD\', kind, domain] - sens[order, \'T2D\', kind, domain]\n                den = sens[order, \'IHD\', \'TE\', None] - sens[order, \'T2D\', \'TE\', None]\n                bden = base[order, \'IHD\', \'TE\', None] - base[order, \'T2D\', \'TE\', None]\n                rows.append(dict(model_id=model, setting=record[\'setting\'], allocation_order=order, endpoint=\'Delta\', component_type=kind, domain_id=domain, main_per_1000=1000 * a, sensitivity_per_1000=1000 * v, change_per_1000=1000 * (v - a), main_PM_percent=100 * a / bden if bden else None, sensitivity_PM_percent=100 * v / den if den else None))\n        evidence = dict(passed=not failures, settings=settings, models=MODELS, results=results, failures=failures, source_manifest_sha256=idx[\'source_manifest\'][\'sha256\'], fixed_bootstrap_draws_added=0, interpretation=\'Prespecified B0 model-dependence analysis, not parameter selection\')\n        sensitivity_key = digest(evidence)\n\n        def write_sensitivity(folder, final):\n            if rows:\n                csv_file(folder / \'B0_PENALTY_SENSITIVITY.csv\', rows)\n            publish(folder / \'SENSITIVITY_ACCEPTANCE.json\', dict(evidence, comparison_csv=str(final / \'B0_PENALTY_SENSITIVITY.csv\') if rows else None))\n            text_file(folder / \'REPORT_ZH.md\', \'# 原样本惩罚敏感性\\n\\n\' + f"完成{sum((r[\'status\'] == \'SUCCESS\' for r in results))}/{len(results)}个B0设置；未新增Bootstrap。\\n\\n" + \'主规格始终df=k/2；附加设置为k/3与2k/3，一次只改变一块。A、抽样、人群和MC随机流不变。数值差异见配套表；不按IE方向或显著性选惩罚。\\n\')\n        sensitivity = publication(target / \'publications\', \'PUBLISH_\', sensitivity_key, write_sensitivity)\n        csv_path = sensitivity / \'B0_PENALTY_SENSITIVITY.csv\'\n        while True:\n            latest = read(index_path)\n            summary = latest.get(\'summary\', {})\n            if summary.get(\'state\') in (\'SUCCESS\', \'FAILED_REVIEW\'):\n                break\n            if latest.get(\'stage\') in (\'FAILED_AWAIT_REVIEW\', \'STOPPED\', \'SUMMARY_FAILED_REVIEW\'):\n                return failure_handoff(run, args.variant, \'主汇总尚未完成且控制器已停止\', dict(controller_stage=latest.get(\'stage\'), summary_state=summary.get(\'state\')), sensitivity)\n            time.sleep(5)\n        main_dir = None\n        if summary.get(\'folder\'):\n            log = Path(summary[\'folder\']) / \'summary.log\'\n            if log.is_file():\n                for line in reversed(log.read_text().splitlines()):\n                    try:\n                        doc = json.loads(line)\n                    except ValueError:\n                        continue\n                    if isinstance(doc, dict) and doc.get(\'summary_dir\'):\n                        main_dir = Path(doc[\'summary_dir\'])\n                        break\n        failure_details = dict(summary_state=summary.get(\'state\'), summary_dir=str(main_dir) if main_dir else None, sensitivity_failures=failures)\n        if summary.get(\'state\') != \'SUCCESS\':\n            return failure_handoff(run, args.variant, \'主汇总执行失败；未读取正式结果表\', failure_details, sensitivity)\n        if main_dir is None:\n            return failure_handoff(run, args.variant, \'主汇总未提供完整发布目录\', failure_details, sensitivity)\n        try:\n            assert main_dir.absolute().is_relative_to(run.absolute()) and (not main_dir.is_symlink()), \'SUMMARY_PATH_OUTSIDE_REPAIR\'\n            main_audit = read(main_dir / \'SUMMARY_AUDIT.json\')\n            assert main_audit.get(\'computation_passed\') is True and main_audit.get(\'full_100_paired_valid\') is True, \'MAIN_202_MC_OR_100_PAIRED_VALID_NOT_ACCEPTED\'\n            assert main_audit.get(\'old_B0_joint_comparison_complete\') is True, \'OLD_B0_COMPARISON_INCOMPLETE\'\n            assert main_audit.get(\'source_manifest\') == idx[\'source_manifest\'], \'MAIN_SUMMARY_SOURCE_ADOPTION_MISMATCH\'\n            assert main_audit.get(\'BH_family_n\') == {MODELS[0]: 5, MODELS[1]: 6} and main_audit.get(\'joint_Delta_IE_separate\') is True, \'MAIN_FROZEN_INFERENCE_CONTRACT_MISMATCH\'\n            if failures:\n                raise RuntimeError(\'B0_SENSITIVITY_INCOMPLETE_OR_FAILED_RESULTS_PRESERVED\')\n            for name in (\'RESULTS\', \'DOMAIN\', \'RATIOS\', \'OLD_B0_COMPARISON\'):\n                assert (main_dir / (name + \'.csv\')).is_file() and (not (main_dir / (name + \'.csv\')).is_symlink()), \'MISSING_MAIN_TABLE:\' + name\n\n            def table(name):\n                with (main_dir / (name + \'.csv\')).open(newline=\'\') as f:\n                    return list(csv.DictReader(f))\n            effects = table(\'RESULTS\')\n            domain = table(\'DOMAIN\')\n            ratios = table(\'RATIOS\')\n            assert effects and domain and ratios, \'EMPTY_MAIN_RESULT_TABLE\'\n        except Exception as error:\n            return failure_handoff(run, args.variant, \'主汇总或必要附加结果不完整：\' + str(error), failure_details, sensitivity)\n        complete = True\n        accepted = dict(schema=\'SELECTED44_PROTEIN_RIDGE_FINAL_COMPUTATION_ACCEPTANCE_V1\', computation_accepted=True, full_new_method_B0_and_B100=True, full_100_true_paired_draws=True, primary_method=args.variant, body_response_penalty_used=args.variant == \'body_ridge\', sensitivity_completed=True, sensitivity_evidence=str(sensitivity / \'SENSITIVITY_ACCEPTANCE.json\'), sensitivity_evidence_sha256=sha(sensitivity / \'SENSITIVITY_ACCEPTANCE.json\'), main_summary=str(main_dir), main_audit_sha256=sha(main_dir / \'SUMMARY_AUDIT.json\'), new_draws=0, B101_started=False, B1000_started=False, note=\'Computational acceptance is not proof of causal identification, unbiasedness, coverage or selection-wide FDR\')\n\n        def fmt(x):\n            if x in (None, \'\'):\n                return \'NA\'\n            try:\n                return f\'{float(x):.6g}\'\n            except (ValueError, TypeError):\n                return str(x)\n        text = [\'# SELECTED44 蛋白联合中介：新版固定B100最终交付\', \'\', \'计算验收：\' + (\'通过。\' if complete else \'尚未通过；有效产物均保留，见验收JSON。\'), \'\', \'范围仅两个蛋白队列模型，原样本9,174人，IHD/T2D、全随访拟合、11年风险、K50及原固定家系B1–100。没有重抽样或启动B101/B1000。\', \'\', \'四个原因别分支统一中心分层；原族裔编码使用原生ridge，主规则为各次拟合df=k/2、Efron、选择权重及robust=FALSE。\', \'使用统一体成分响应岭惩罚后备，不能视为纯数值修复。\' if args.variant == \'body_ridge\' else \'未启用体成分响应惩罚；该块仅做等价缩放和全列可逆QR。\', \'兼容A及实际索引复用，旧Y与旧成功重复不混入新区间。\', \'\', \'## 联合效应\', \'\', \'单位为每1,000人的11年绝对风险差。区间来自实际共同有效编号的Type7百分位。\', \'\', \'| 模型 | 结局 | TE | IE | IE 95%CI | 有效重复 |\', \'|---|---|---:|---:|---|---:|\']\n        for m in MODELS:\n            for e in (\'IHD\', \'T2D\', \'Delta\'):\n                selected = {r[\'component_type\']: r for r in effects if r[\'model_id\'] == m and r[\'allocation_order\'] == \'forward\' and (r[\'endpoint\'] == e)}\n                if \'IE_joint\' not in selected:\n                    continue\n                ie = selected[\'IE_joint\']\n                te = selected.get(\'TE\', {})\n                text.append(\'| \' + \' | \'.join((m, e, fmt(te.get(\'point_B0_per_1000\')), fmt(ie.get(\'point_B0_per_1000\')), \'[\' + fmt(ie.get(\'CI_type7_low_per_1000\')) + \', \' + fmt(ie.get(\'CI_type7_high_per_1000\')) + \']\', fmt(ie.get(\'effect_valid_N\')))) + \' |\')\n        text += [\'\', \'## 正式蛋白域贡献\', \'\', \'这来自PLUS15的正式域分解，不是PLUS−REF联合IE差。\', \'\', \'| 结局 | 蛋白域IIE/1,000 | 95%CI | P | q（仅前向Δ，BH6） |\', \'|---|---:|---|---:|---:|\']\n        for r in domain:\n            if r[\'model_id\'] == MODELS[1] and r[\'allocation_order\'] == \'forward\' and (r[\'component_type\'] == \'IIE\') and (r[\'domain_id\'] == \'CANDIDATE_PROTEINS\'):\n                text.append(\'| \' + \' | \'.join((r[\'endpoint\'], fmt(r.get(\'point_B0_per_1000\')), \'[\' + fmt(r.get(\'CI_type7_low_per_1000\')) + \', \' + fmt(r.get(\'CI_type7_high_per_1000\')) + \']\', fmt(r.get(\'P\')), fmt(r.get(\'q\')))) + \' |\')\n        text += [\'\', \'## 惩罚敏感性及限制\', \'\', f\'B0预指定附加设置完成{len(results) - len(failures)}/{len(results)}；同一A、样本、权重及MC随机流，一次只改一个惩罚块。\', \'实际逐项变化见 \' + str(csv_path) + \'。敏感性不带未经计算的Bootstrap区间，也不用于挑选更显著的设置。\']\n        relevant = [r for r in rows if r[\'allocation_order\'] == \'forward\' and (r[\'component_type\'] == \'IE_joint\' or r.get(\'domain_id\') == \'CANDIDATE_PROTEINS\')]\n        flips = [r for r in relevant if r[\'main_per_1000\'] * r[\'sensitivity_per_1000\'] < 0]\n        text.append(\'上述关键前向点估计出现符号改变的敏感性记录数：\' + str(len(flips)) + \'；最大绝对变化/1,000：\' + fmt(max((abs(r[\'change_per_1000\']) for r in relevant), default=None)) + \'。小的近零点估计翻号不自动等于重要结论改变，需结合主区间及绝对量解释。\')\n        unstable = sum((str(r.get(\'denominator_unstable\', \'\')).lower() == \'true\' for r in ratios))\n        text += [\'\', f\'比值分母不稳定记录数：{unstable}。PM/Psep不裁剪；近零或跨零TE及Fieller无界/不连通集合均在原表保留。\', \'P、q沿用冻结的B100 Type7区间反演及模型内域差值+DEP家族（REF BH5、PLUS BH6）；联合Δ单列，不套用BH18。B100的尾部与P/q精度有限，不是独立确认或筛选全流程FDR保证。\', \'ridge的df约定不是最优参数或无偏保证；计算通过并不证明所有中介可交换性、模型正确性或竞争风险因果假设。\', \'\', \'## 文件\', \'\', str(main_dir), str(csv_path)]\n\n        def write_release(folder, final):\n            publish(folder / \'RELEASE_ACCEPTANCE.json\', accepted)\n            text_file(folder / \'RELEASE_REPORT_ZH.md\', \'\\n\'.join(text + [\'\', str(final / \'RELEASE_ACCEPTANCE.json\')]) + \'\\n\')\n            publish(folder / \'OUTPUT_MANIFEST.json\', dict(files=[dict(path=str(final / p.name), sha256=sha(p)) for p in folder.iterdir() if p.is_file()], main_summary=str(main_dir), sensitivity=str(sensitivity), actual_results_only=True))\n        final = publication(run / args.variant / \'release_v2\', \'RELEASE_\', digest(accepted), write_release)\n        print(canonical(dict(ok=True, computation_accepted=True, release=str(final), sensitivity=str(sensitivity))))\n        return 0\nif __name__ == \'__main__\':\n    try:\n        sys.exit(main())\n    except Exception as error:\n        print(canonical(dict(ok=False, exit_code=2, reason=str(error), existing_outputs_preserved=True)), file=sys.stderr)\n        sys.exit(2)', '<private-stage>', 'exec'), context)
