"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
METHOD = birth_private_setting('PRIVATE_9F6C8582BD7A7EF1BE56E3FA93AB415C076C49D26C480D9C912159595545')
FAMILY = 'SELECTED18_IHD_MINUS_T2D_DELTA_IE_BH18_20261031_V1'
LABELS = {'bmi': 'BMI', 'clinical__glucose': '葡萄糖', 'clinical__hba1c': '糖化血红蛋白', 'clinical__hdl': 'HDL-C', 'clinical__tg': '甘油三酯', 'fat_free_mass': '去脂体重', 'fat_mass': '脂肪量', 'waist': '腰围', 'smoking': '吸烟状态', 'hip': '臀围', 'nmr__nmr_23596': '中型VLDL胆固醇酯比例', 'nmr__nmr_23611': 'IDL胆固醇酯比例', 'nmr__nmr_23618': '大型LDL甘油三酯比例', 'nmr__nmr_23639': '中型HDL磷脂比例', 'olink__fcrl6': 'FCRL6', 'olink__gzma': 'GZMA', 'olink__gzmh': 'GZMH', 'olink__klrd1': 'KLRD1'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def bh(p):
    order = sorted(range(len(p)), key=p.__getitem__)
    adjusted = [0.0] * len(p)
    running = 1.0
    for rank in range(len(p), 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * len(p) / rank)
        adjusted[i] = running
    return adjusted

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=os.getenv('B18_BH_INPUT', SOURCE))
    ap.add_argument('--output', type=Path, default=os.getenv('B18_BH_OUT', OUTPUT))
    args = ap.parse_args()
    source, out = (args.input.resolve(), args.output.resolve())
    print(json.dumps({'input': str(source), 'output': str(out)}, ensure_ascii=False))
    assert source.is_file(), f'MISSING_INPUT:{source}'
    before = sha(source)
    assert before == EXPECTED_SOURCE_SHA, 'SOURCE_VERSION_MISMATCH_NO_SUBSTITUTION'
    with source.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        fields, rows = (reader.fieldnames, list(reader))
    assert len(rows) == 18 and len({r['indicator_id'] for r in rows}) == 18
    assert {r['indicator_id'] for r in rows} == set(LABELS)
    p = []
    for r in rows:
        value = float(r['P'] or r['P_Delta_IE'])
        assert math.isfinite(value) and 0 <= value <= 1
        if r['P'] and r['P_Delta_IE']:
            assert r['P'] == r['P_Delta_IE'], 'CONFLICTING_EXISTING_P'
        assert r['method_hash'] == METHOD
        assert r['contrast'] in ('', '411_minus_250.2')
        assert r['horizon_years'] in ('', '11')
        assert r['current_count_point_pair_hash_verified'].lower() == 'true'
        assert int(r['bootstrap_planned']) == 1000
        assert 0 < int(r['paired_valid']) <= 1000
        assert all((math.isfinite(float(r[k])) for k in ('Delta_IE_per_1000', 'Delta_IE_CI_low_per_1000', 'Delta_IE_CI_high_per_1000')))
        p.append(value)
    q = bh(p)
    sorted_p = sorted(p)
    direct = [min(1.0, min((sorted_p[j] * 18 / (j + 1) for j in range(18) if sorted_p[j] >= value))) for value in p]
    assert all((math.isclose(a, b, rel_tol=1e-14, abs_tol=1e-16) for a, b in zip(q, direct)))
    assert all((value <= adjusted <= 1 for value, adjusted in zip(p, q)))
    ordered_q = [q[i] for i in sorted(range(18), key=p.__getitem__)]
    assert ordered_q == sorted(ordered_q)
    extra = ['label_zh_BH18', 'q_BH18', 'BH18_family_id', 'BH18_family_m', 'BH18_contrast', 'BH18_horizon_years', 'BH18_q_lt_0_05', 'BH18_inference_scope', 'BH18_adopted_beijing']
    assert not set(extra).intersection(fields)
    stamp = datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
    exports = []
    for r, adjusted in zip(rows, q):
        exports.append(dict(r, label_zh_BH18=LABELS[r['indicator_id']], q_BH18=repr(adjusted), BH18_family_id=FAMILY, BH18_family_m='18', BH18_contrast='411_minus_250.2', BH18_horizon_years='11', BH18_q_lt_0_05=str(adjusted < 0.05), BH18_inference_scope='POST_B100_SELECTED_PANEL_NOT_SELECTION_ADJUSTED_FDR;SOURCE_INFERENCE_STATUS_RETAINED', BH18_adopted_beijing=stamp))
    out.mkdir(parents=True, exist_ok=False)
    target = out / 'B1000_SELECTED18_IHD_T2D_BH18.csv'
    with target.open('x', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields + extra)
        writer.writeheader()
        writer.writerows(exports)
    with target.open(encoding='utf-8-sig', newline='') as f:
        reread = list(csv.DictReader(f))
    assert [{k: r[k] for k in fields} for r in reread] == rows
    assert sha(source) == before
    full = sum((int(r['paired_valid']) == 1000 for r in rows))
    lines = ['# 最终18项 IHD−T2D 配对ΔIE：BH18汇总', '', f'生成时间（北京时间）：{stamp}', '', '本次仅对最终18项已有的 IHD−T2D ΔIE P值做一次 BH 校正，固定 m=18；不是18×5或18×6家族。', 'BH18用于已选择的18项，不代表预先指定的全候选家族。原P、区间、q及q_family逐字段保留；新结果只在q_BH18列。', '未重拟合模型、未新增Bootstrap、未重算P值、未重新筛选或替换指标；旧16项及B100的家族和校正表均未覆盖。', '', f'18/18项 q_BH18<0.05；其中{full}项配对有效1000/1000，{18 - full}项为部分重复条件推断。', '', '|指标|ΔIE /1000 [95% CI]|原P|BH18 q|配对有效B/1000|', '|---|---|---|---|---|']
    for r, value, adjusted in zip(rows, p, q):
        effect = float(r['Delta_IE_per_1000'])
        low, high = (float(r['Delta_IE_CI_low_per_1000']), float(r['Delta_IE_CI_high_per_1000']))
        lines.append(f"|{LABELS[r['indicator_id']]}|{effect:.4f} [{low:.4f}, {high:.4f}]|{value:.8f}|{adjusted:.8f}|{r['paired_valid']}/1000|")
    lines += ['', '## 解释边界', '', '- ΔIE=IE(IHD)−IE(T2D)，全部随访拟合、11年风险报告；正值不代表IHD的IE必然为正。', '- 这是同一数据B100筛选后的重点面板精化，BH18只调整所选18项，不控制最初121项筛选流程的整体FDR，不是独立确认。', '- P沿用冻结的Type-7区间反演程序及实际有效n的离散网格，不是新算的精确尾概率；没有以零替代缺失重复。', '- 葡萄糖、HDL-C及四蛋白仍保留部分重复条件推断；蛋白保留warning-aware政策与其数值/缺失限制。校正不能消除这些限制。', '- 12项完整配对重复不意味着五疾病全部结果均完整。本表只处理IHD−T2D ΔIE。', '- NMR四项是相应成分比例，不能当作绝对浓度；本表不计算或合并IE、TE、PM。', '', '## 来源与复现', '', f'- 来源：`{source}`', f'- 来源SHA256：`{before}`', f'- 方法哈希：`{METHOD}`', f'- 输出：`{target}`', '- 核验：18项唯一名单、P字段一致、方法/当前配对审阅标记、BH两种实现一致、原字段回读逐值不变、源文件SHA不变。', '- 复现到新目录（现有目录不允许覆盖）：', '', '```bash', f'python3 "{Path(__file__).resolve()}" --input "{source}" --output "{out}_reproduce"', '```', '']
    report = out / 'REPORT_ZH.md'
    with report.open('x', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    manifest = {'schema': FAMILY, 'created_beijing': stamp, 'reviewer': 'Codex (OpenAI)', 'authorization': 'BH adjustment within the selected 18-indicator comparison', 'input': str(source), 'input_sha256': before, 'input_sha256_after': sha(source), 'family_m': 18, 'contrast': 'IHD_MINUS_T2D_DELTA_IE', 'horizon_years': 11, 'method_hash': METHOD, 'source_hashes_preserved': sorted({r['source_hash'] for r in rows}), 'status_counts_preserved': dict(Counter((r['result_status'] for r in rows))), 'all_source_columns_preserved': True, 'models_fit': 0, 'new_bootstraps': 0, 'p_recalculated': False, 'historical_q_overwritten': False, 'post_selection_independent_confirmation': False, 'controls_original_121_selection_FDR': False, 'complete_paired_1000': full, 'partial_paired_inference': 18 - full, 'q_lt_0_05': sum((value < 0.05 for value in q)), 'files_sha256': {p.name: sha(p) for p in (target, report)}, 'script': str(Path(__file__).resolve()), 'script_sha256': sha(Path(__file__))}
    with (out / 'METHOD_NOTE_AND_MANIFEST.json').open('x', encoding='utf-8') as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
        f.write('\n')
    for r in exports:
        print(json.dumps({k: r[k] for k in ('indicator_id', 'label_zh_BH18', 'P', 'q_BH18', 'paired_valid')}, ensure_ascii=False))
    print(json.dumps({'output': str(out), 'all18_preserved': True, 'q_lt_0_05': manifest['q_lt_0_05']}, ensure_ascii=False))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Additive BH18 readout; never refits models or changes a frozen release."""\nimport argparse\nimport csv\nimport hashlib\nimport json\nimport math\nimport os\nfrom collections import Counter\nfrom datetime import datetime\nfrom pathlib import Path\nfrom zoneinfo import ZoneInfo\nSOURCE = Path(birth_private_setting(\'PRIVATE_DISPLAY18_PRIMARY_REVIEW_CSV\'))\nOUTPUT = Path(birth_private_setting(\'PRIVATE_B1000_SELECTED18_BH18_20261031_V1\'))\nEXPECTED_SOURCE_SHA = birth_private_setting(\'PRIVATE_7CB3350EAB72CE68574C9B4D0AC9B76051A70A06D2F6DD1851C4C53D7DF2\')\nMETHOD = birth_private_setting(\'PRIVATE_9F6C8582BD7A7EF1BE56E3FA93AB415C076C49D26C480D9C912159595545\')\nFAMILY = \'SELECTED18_IHD_MINUS_T2D_DELTA_IE_BH18_20261031_V1\'\nLABELS = {\'bmi\': \'BMI\', \'clinical__glucose\': \'葡萄糖\', \'clinical__hba1c\': \'糖化血红蛋白\', \'clinical__hdl\': \'HDL-C\', \'clinical__tg\': \'甘油三酯\', \'fat_free_mass\': \'去脂体重\', \'fat_mass\': \'脂肪量\', \'waist\': \'腰围\', \'smoking\': \'吸烟状态\', \'hip\': \'臀围\', \'nmr__nmr_23596\': \'中型VLDL胆固醇酯比例\', \'nmr__nmr_23611\': \'IDL胆固醇酯比例\', \'nmr__nmr_23618\': \'大型LDL甘油三酯比例\', \'nmr__nmr_23639\': \'中型HDL磷脂比例\', \'olink__fcrl6\': \'FCRL6\', \'olink__gzma\': \'GZMA\', \'olink__gzmh\': \'GZMH\', \'olink__klrd1\': \'KLRD1\'}\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef bh(p):\n    order = sorted(range(len(p)), key=p.__getitem__)\n    adjusted = [0.0] * len(p)\n    running = 1.0\n    for rank in range(len(p), 0, -1):\n        i = order[rank - 1]\n        running = min(running, p[i] * len(p) / rank)\n        adjusted[i] = running\n    return adjusted\n\ndef main():\n    ap = argparse.ArgumentParser(description=__doc__)\n    ap.add_argument(\'--input\', type=Path, default=os.getenv(\'B18_BH_INPUT\', SOURCE))\n    ap.add_argument(\'--output\', type=Path, default=os.getenv(\'B18_BH_OUT\', OUTPUT))\n    args = ap.parse_args()\n    source, out = (args.input.resolve(), args.output.resolve())\n    print(json.dumps({\'input\': str(source), \'output\': str(out)}, ensure_ascii=False))\n    assert source.is_file(), f\'MISSING_INPUT:{source}\'\n    before = sha(source)\n    assert before == EXPECTED_SOURCE_SHA, \'SOURCE_VERSION_MISMATCH_NO_SUBSTITUTION\'\n    with source.open(encoding=\'utf-8-sig\', newline=\'\') as f:\n        reader = csv.DictReader(f)\n        fields, rows = (reader.fieldnames, list(reader))\n    assert len(rows) == 18 and len({r[\'indicator_id\'] for r in rows}) == 18\n    assert {r[\'indicator_id\'] for r in rows} == set(LABELS)\n    p = []\n    for r in rows:\n        value = float(r[\'P\'] or r[\'P_Delta_IE\'])\n        assert math.isfinite(value) and 0 <= value <= 1\n        if r[\'P\'] and r[\'P_Delta_IE\']:\n            assert r[\'P\'] == r[\'P_Delta_IE\'], \'CONFLICTING_EXISTING_P\'\n        assert r[\'method_hash\'] == METHOD\n        assert r[\'contrast\'] in (\'\', \'411_minus_250.2\')\n        assert r[\'horizon_years\'] in (\'\', \'11\')\n        assert r[\'current_count_point_pair_hash_verified\'].lower() == \'true\'\n        assert int(r[\'bootstrap_planned\']) == 1000\n        assert 0 < int(r[\'paired_valid\']) <= 1000\n        assert all((math.isfinite(float(r[k])) for k in (\'Delta_IE_per_1000\', \'Delta_IE_CI_low_per_1000\', \'Delta_IE_CI_high_per_1000\')))\n        p.append(value)\n    q = bh(p)\n    sorted_p = sorted(p)\n    direct = [min(1.0, min((sorted_p[j] * 18 / (j + 1) for j in range(18) if sorted_p[j] >= value))) for value in p]\n    assert all((math.isclose(a, b, rel_tol=1e-14, abs_tol=1e-16) for a, b in zip(q, direct)))\n    assert all((value <= adjusted <= 1 for value, adjusted in zip(p, q)))\n    ordered_q = [q[i] for i in sorted(range(18), key=p.__getitem__)]\n    assert ordered_q == sorted(ordered_q)\n    extra = [\'label_zh_BH18\', \'q_BH18\', \'BH18_family_id\', \'BH18_family_m\', \'BH18_contrast\', \'BH18_horizon_years\', \'BH18_q_lt_0_05\', \'BH18_inference_scope\', \'BH18_adopted_beijing\']\n    assert not set(extra).intersection(fields)\n    stamp = datetime.now(ZoneInfo(\'Asia/Shanghai\')).isoformat()\n    exports = []\n    for r, adjusted in zip(rows, q):\n        exports.append(dict(r, label_zh_BH18=LABELS[r[\'indicator_id\']], q_BH18=repr(adjusted), BH18_family_id=FAMILY, BH18_family_m=\'18\', BH18_contrast=\'411_minus_250.2\', BH18_horizon_years=\'11\', BH18_q_lt_0_05=str(adjusted < 0.05), BH18_inference_scope=\'POST_B100_SELECTED_PANEL_NOT_SELECTION_ADJUSTED_FDR;SOURCE_INFERENCE_STATUS_RETAINED\', BH18_adopted_beijing=stamp))\n    out.mkdir(parents=True, exist_ok=False)\n    target = out / \'B1000_SELECTED18_IHD_T2D_BH18.csv\'\n    with target.open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as f:\n        writer = csv.DictWriter(f, fieldnames=fields + extra)\n        writer.writeheader()\n        writer.writerows(exports)\n    with target.open(encoding=\'utf-8-sig\', newline=\'\') as f:\n        reread = list(csv.DictReader(f))\n    assert [{k: r[k] for k in fields} for r in reread] == rows\n    assert sha(source) == before\n    full = sum((int(r[\'paired_valid\']) == 1000 for r in rows))\n    lines = [\'# 最终18项 IHD−T2D 配对ΔIE：BH18汇总\', \'\', f\'生成时间（北京时间）：{stamp}\', \'\', \'本次仅对最终18项已有的 IHD−T2D ΔIE P值做一次 BH 校正，固定 m=18；不是18×5或18×6家族。\', \'BH18用于已选择的18项，不代表预先指定的全候选家族。原P、区间、q及q_family逐字段保留；新结果只在q_BH18列。\', \'未重拟合模型、未新增Bootstrap、未重算P值、未重新筛选或替换指标；旧16项及B100的家族和校正表均未覆盖。\', \'\', f\'18/18项 q_BH18<0.05；其中{full}项配对有效1000/1000，{18 - full}项为部分重复条件推断。\', \'\', \'|指标|ΔIE /1000 [95% CI]|原P|BH18 q|配对有效B/1000|\', \'|---|---|---|---|---|\']\n    for r, value, adjusted in zip(rows, p, q):\n        effect = float(r[\'Delta_IE_per_1000\'])\n        low, high = (float(r[\'Delta_IE_CI_low_per_1000\']), float(r[\'Delta_IE_CI_high_per_1000\']))\n        lines.append(f"|{LABELS[r[\'indicator_id\']]}|{effect:.4f} [{low:.4f}, {high:.4f}]|{value:.8f}|{adjusted:.8f}|{r[\'paired_valid\']}/1000|")\n    lines += [\'\', \'## 解释边界\', \'\', \'- ΔIE=IE(IHD)−IE(T2D)，全部随访拟合、11年风险报告；正值不代表IHD的IE必然为正。\', \'- 这是同一数据B100筛选后的重点面板精化，BH18只调整所选18项，不控制最初121项筛选流程的整体FDR，不是独立确认。\', \'- P沿用冻结的Type-7区间反演程序及实际有效n的离散网格，不是新算的精确尾概率；没有以零替代缺失重复。\', \'- 葡萄糖、HDL-C及四蛋白仍保留部分重复条件推断；蛋白保留warning-aware政策与其数值/缺失限制。校正不能消除这些限制。\', \'- 12项完整配对重复不意味着五疾病全部结果均完整。本表只处理IHD−T2D ΔIE。\', \'- NMR四项是相应成分比例，不能当作绝对浓度；本表不计算或合并IE、TE、PM。\', \'\', \'## 来源与复现\', \'\', f\'- 来源：`{source}`\', f\'- 来源SHA256：`{before}`\', f\'- 方法哈希：`{METHOD}`\', f\'- 输出：`{target}`\', \'- 核验：18项唯一名单、P字段一致、方法/当前配对审阅标记、BH两种实现一致、原字段回读逐值不变、源文件SHA不变。\', \'- 复现到新目录（现有目录不允许覆盖）：\', \'\', \'```bash\', f\'python3 "{Path(__file__).resolve()}" --input "{source}" --output "{out}_reproduce"\', \'```\', \'\']\n    report = out / \'REPORT_ZH.md\'\n    with report.open(\'x\', encoding=\'utf-8\') as f:\n        f.write(\'\\n\'.join(lines))\n    manifest = {\'schema\': FAMILY, \'created_beijing\': stamp, \'reviewer\': \'Codex (OpenAI)\', \'authorization\': \'BH adjustment within the selected 18-indicator comparison\', \'input\': str(source), \'input_sha256\': before, \'input_sha256_after\': sha(source), \'family_m\': 18, \'contrast\': \'IHD_MINUS_T2D_DELTA_IE\', \'horizon_years\': 11, \'method_hash\': METHOD, \'source_hashes_preserved\': sorted({r[\'source_hash\'] for r in rows}), \'status_counts_preserved\': dict(Counter((r[\'result_status\'] for r in rows))), \'all_source_columns_preserved\': True, \'models_fit\': 0, \'new_bootstraps\': 0, \'p_recalculated\': False, \'historical_q_overwritten\': False, \'post_selection_independent_confirmation\': False, \'controls_original_121_selection_FDR\': False, \'complete_paired_1000\': full, \'partial_paired_inference\': 18 - full, \'q_lt_0_05\': sum((value < 0.05 for value in q)), \'files_sha256\': {p.name: sha(p) for p in (target, report)}, \'script\': str(Path(__file__).resolve()), \'script_sha256\': sha(Path(__file__))}\n    with (out / \'METHOD_NOTE_AND_MANIFEST.json\').open(\'x\', encoding=\'utf-8\') as f:\n        json.dump(manifest, f, ensure_ascii=False, indent=2)\n        f.write(\'\\n\')\n    for r in exports:\n        print(json.dumps({k: r[k] for k in (\'indicator_id\', \'label_zh_BH18\', \'P\', \'q_BH18\', \'paired_valid\')}, ensure_ascii=False))\n    print(json.dumps({\'output\': str(out), \'all18_preserved\': True, \'q_lt_0_05\': manifest[\'q_lt_0_05\']}, ensure_ascii=False))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
