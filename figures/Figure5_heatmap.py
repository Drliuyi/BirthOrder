"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import types

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f))

def stars(q):
    q = float(q)
    return '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--no-significance-legend', action='store_true', help='Keep cell stars, omit their threshold explanation from the figure.')
    a = ap.parse_args()
    out = a.out.absolute()
    assert not out.exists(), 'OUTPUT_MUST_BE_NEW'
    previous = ROOT / 'final_heatmap_only_v1'
    source = ROOT / 'block_axes_v12/draw_figure2b_block_axes.py'
    frozen = previous / 'PLOTDATA_HEATMAP_EFFECTS_PERCENTAGES_PQ.csv'
    pins = {str(p): sha(p) for p in (source, frozen, previous / 'PLOTDATA_JOINT_RESULTS.csv', previous / 'PLOTDATA_JOINT_FOREST.csv')}
    data = read(frozen)
    assert len(data) == 42
    tiers = {s: sum((stars(r['q_BH']) == s for r in data)) for s in ('', '*', '**', '***')}
    assert tiers == {'': 28, '*': 14, '**': 0, '***': 0}
    code = source.read_text()
    edits = [("    SOURCE = ROOT / 'horizontal_v3'", "    SOURCE = ROOT / 'horizontal_v3'"), ("                    ha=ha, va='center', color=INK, zorder=5)", "                    ha=ha, va='center', color=INK, zorder=5)"), ('t = ax.text(x, y, s, fontsize=size, fontweight=', 't = ax.text(x, y, s, fontsize=size * .90, fontweight='), ("                text(x + CELL_W / 2 - (.6 if marked else 0), y + 3.1, fmt(val), 8.8, ha='center')\n", "                text(x + CELL_W / 2, y + (RH - .4) / 2, fmt(val), 9.6, ha='center')\n"), ("                text(x + CELL_W / 2, y + 7.4, '(' + r['proportion_display'] + ')', 8.0, ha='center')\n", ''), ("                    mark = text(x + CELL_W - 1.35, y + 3.1, '*', 12.8, True, ha='center')\n", "                    tier = '***' if q < .001 else '**' if q < .01 else '*'\n                    mark = text(x + CELL_W - 1.15, y + (RH - .4) / 2 - 1.1, tier, 11.3, True, ha='center')\n"), ('stars.append(dict(model_id=mid, domain_id=gid, outcome=endpoint, q_BH=q))', 'stars.append(dict(model_id=mid, domain_id=gid, outcome=endpoint, q_BH=q, symbol=tier))'), ("        text(gx + GRID_W, bar_y + 9.0, nominal_significance, 7.1, ha='right')\n", ''), ("    text(rx + 1.0, 83.8, '*', 11.5, True, ha='center')\n    text(rx + 4.0, 83.8, 'BH q < 0.05', 6.7)\n", "    for yy, symbol, threshold in ((83.8, '*', '0.05'), (88.4, '**', '0.01'), (93.0, '***', '0.001')):\n        text(rx + .5, yy, symbol, 10.5, True, ha='center')\n        text(rx + 4.0, yy, 'BH q < ' + threshold, 6.7)\n"), ('heatmap_asterisk_size_pt=12.8, previous_asterisk_size_pt=7.2,', 'heatmap_asterisk_size_pt=11.3 * .9, previous_asterisk_size_pt=12.8,'), ("joint_forest_significance_rule='EXPLICIT_ORIGINAL_P_THRESHOLD_NOT_BH_ASTERISK'", "joint_forest_significance_rule='SAVED_IN_DATA_NOT_DISPLAYED'"), ("significance_key_text='BH q < 0.05'", "significance_key_text='* BH q<0.05; ** BH q<0.01; *** BH q<0.001'"), ("(out / 'FIGURE_LEGEND_EN.md').write_text(", "(out / 'RENDERER_LEGEND_ORIGINAL.md').write_text("), ("(out / 'REPRODUCE.md').write_text(", "(out / 'REPRODUCE_ORIGINAL.md').write_text(")]
    for before, after in edits:
        if before == after:
            continue
        assert code.count(before) == 1, ('SOURCE_ANCHOR_DRIFT', before)
        code = code.replace(before, after)
    if a.no_significance_legend:
        significance_key = "    for yy, symbol, threshold in ((83.8, '*', '0.05'), (88.4, '**', '0.01'), (93.0, '***', '0.001')):\n        text(rx + .5, yy, symbol, 10.5, True, ha='center')\n        text(rx + 4.0, yy, 'BH q < ' + threshold, 6.7)\n"
        assert code.count(significance_key) == 1
        code = code.replace(significance_key, '')
        code = code.replace("significance_key_text='* BH q<0.05; ** BH q<0.01; *** BH q<0.001'", "significance_key_text='NOT_DISPLAYED_RESEARCHER_SUPPLIES_LEGEND'")
    code = code.replace('Figure2b_Heatmap.', 'Figure2b_Heatmap_Multistar.')
    code = code.replace('Figure2b_Heatmap_preview', 'Figure2b_Heatmap_Multistar_preview')
    code = code.replace('Figure2b_Heatmap_180mm_check', 'Figure2b_Heatmap_Multistar_180mm_check')
    module = types.ModuleType('figure2b_multistar_renderer')
    module.__file__ = str(source)
    exec(compile(code, str(source), 'exec'), module.__dict__)
    argv = sys.argv
    try:
        sys.argv = [str(source), '--out', str(out)]
        module.main()
    finally:
        sys.argv = argv
    for name in ('PLOTDATA_HEATMAP_EFFECTS_PERCENTAGES_PQ.csv', 'PLOTDATA_JOINT_RESULTS.csv', 'PLOTDATA_JOINT_FOREST.csv'):
        assert (out / name).read_bytes() == (previous / name).read_bytes(), name
    qc = json.loads((out / 'QC.json').read_text())
    assert len(qc['stars']) == 14
    assert [r['symbol'] for r in qc['stars']] == ['*'] * 14
    assert not qc['text_clipped'] and (not qc['graph_text_overlap'])
    assert pins == {p: sha(Path(p)) for p in pins}
    ledger = [dict(model_id=r['model_id'], domain_id=r['domain_id'], outcome=r['outcome'], point_per_1000=r['point_per_1000'], original_P=r['P'], original_q_BH=r['q_BH'], displayed_stars=stars(r['q_BH']), test_family=r['test_family']) for r in data]
    with (out / 'GRADED_STARS_AUDIT.csv').open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(ledger[0]))
        w.writeheader()
        w.writerows(ledger)
    (out / 'FIGURE_LEGEND_EN.md').write_text('Figure 2b | Joint mediation across three analysis populations. Original-sample domain IE contributions are shown as 11-year standardized absolute risk differences per 1,000. Asterisks use each original, unrounded BH q: * q<0.05, ** q<0.01, *** q<0.001. All 14 statistically supported cells in this frozen B100 table qualify only for a single star (minimum q=0.0125); none qualifies for a double or triple star. BH correction is separate within each model×effect-type family (m=5/6 including the unshown DEP), not whole-figure FDR. Single-disease tests are supplementary exploratory tests; ΔIE tests are the original paired comparisons. The three mediated-proportion bars and joint ΔIE forest estimates and 95% CIs are unchanged. No joint BH q is available, so no BH stars are assigned to those bars or forests. All original effects, CIs, P/q and result versions are unchanged. No models, bootstrap draws or P/q calculations were rerun.\n')
    (out / 'BUILD_REPORT_ZH.md').write_text('# Figure 2b：多级星号规则与小字号版\n\n仅修改绘图，全部字体较当前图缩小10%；数字9.6→8.64 pt，星号12.8→10.17 pt。热图尺寸、行列顺序、颜色映射、领域及队列色条、比例条和Joint ΔIE森林图完全保留。\n\n多级规则按原始未舍入BH q：*<0.05、**<0.01、***<0.001。当前42格为14格单星、0格双星、0格三星、28格无星；最小q=0.0125。不能为了出现多星而使用原P替代q或修改阈值。\n\n三份统计CSV与原final_heatmap_only_v1逐字节一致；无新增拟合、Bootstrap或P/q。PNG600dpi、矢量PDF、可编辑SVG已生成，文字无裁切/重叠。\n')
    (out / 'REPRODUCE.md').write_text('```bash\n' + f'/opt/homebrew/Caskroom/miniforge/base/bin/python3 "{Path(__file__).resolve()}"' + (' --no-significance-legend' if a.no_significance_legend else '') + f''' --out "{out.parent / (out.name + '_reproduced')}"\n```\n''')
    (out / 'REVISION_QC.json').write_text(json.dumps(dict(source_sha256=pins, star_counts=tiers, minimum_original_q=min((float(r['q_BH']) for r in data)), font_scale=0.9, significance_legend_displayed=not a.no_significance_legend, statistics_byte_identical=True, geometry_palette_order_unchanged=True, new_fits=0, new_bootstrap=0, new_P_q=0, text_clipped=qc['text_clipped'], text_overlap=qc['graph_text_overlap']), ensure_ascii=False, indent=2) + '\n')
    with (out / Path(__file__).name).open('xb') as f:
        f.write(Path(__file__).read_bytes())
    (out / 'RENDERED_SOURCE.py').write_text(code)
    (out / 'FINAL_MANIFEST.json').write_text(json.dumps(dict(files=[dict(name=p.name, sha256=sha(p)) for p in sorted(out.iterdir())]), indent=2) + '\n')
    print(json.dumps(dict(output=str(out), star_counts=tiers, font_scale=0.9, statistics_unchanged=True), ensure_ascii=False))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Rendering-only revision: graded BH-q notation and 10% smaller text."""\nimport argparse\nimport csv\nimport hashlib\nimport json\nfrom pathlib import Path\nimport sys\nimport types\nROOT = next((p for p in Path(__file__).resolve().parents if (p / \'final_heatmap_only_v1/PLOTDATA_HEATMAP_EFFECTS_PERCENTAGES_PQ.csv\').is_file()))\n\ndef sha(p):\n    return hashlib.sha256(p.read_bytes()).hexdigest()\n\ndef read(p):\n    with p.open(newline=\'\') as f:\n        return list(csv.DictReader(f))\n\ndef stars(q):\n    q = float(q)\n    return \'***\' if q < 0.001 else \'**\' if q < 0.01 else \'*\' if q < 0.05 else \'\'\n\ndef main():\n    ap = argparse.ArgumentParser(description=__doc__)\n    ap.add_argument(\'--out\', type=Path, required=True)\n    ap.add_argument(\'--no-significance-legend\', action=\'store_true\', help=\'Keep cell stars, omit their threshold explanation from the figure.\')\n    a = ap.parse_args()\n    out = a.out.absolute()\n    assert not out.exists(), \'OUTPUT_MUST_BE_NEW\'\n    previous = ROOT / \'final_heatmap_only_v1\'\n    source = ROOT / \'block_axes_v12/draw_figure2b_block_axes.py\'\n    frozen = previous / \'PLOTDATA_HEATMAP_EFFECTS_PERCENTAGES_PQ.csv\'\n    pins = {str(p): sha(p) for p in (source, frozen, previous / \'PLOTDATA_JOINT_RESULTS.csv\', previous / \'PLOTDATA_JOINT_FOREST.csv\')}\n    data = read(frozen)\n    assert len(data) == 42\n    tiers = {s: sum((stars(r[\'q_BH\']) == s for r in data)) for s in (\'\', \'*\', \'**\', \'***\')}\n    assert tiers == {\'\': 28, \'*\': 14, \'**\': 0, \'***\': 0}\n    code = source.read_text()\n    edits = [("    SOURCE = ROOT / \'horizontal_v3\'", "    SOURCE = ROOT / \'horizontal_v3\'"), ("                    ha=ha, va=\'center\', color=INK, zorder=5)", "                    ha=ha, va=\'center\', color=INK, zorder=5)"), (\'t = ax.text(x, y, s, fontsize=size, fontweight=\', \'t = ax.text(x, y, s, fontsize=size * .90, fontweight=\'), ("                text(x + CELL_W / 2 - (.6 if marked else 0), y + 3.1, fmt(val), 8.8, ha=\'center\')\\n", "                text(x + CELL_W / 2, y + (RH - .4) / 2, fmt(val), 9.6, ha=\'center\')\\n"), ("                text(x + CELL_W / 2, y + 7.4, \'(\' + r[\'proportion_display\'] + \')\', 8.0, ha=\'center\')\\n", \'\'), ("                    mark = text(x + CELL_W - 1.35, y + 3.1, \'*\', 12.8, True, ha=\'center\')\\n", "                    tier = \'***\' if q < .001 else \'**\' if q < .01 else \'*\'\\n                    mark = text(x + CELL_W - 1.15, y + (RH - .4) / 2 - 1.1, tier, 11.3, True, ha=\'center\')\\n"), (\'stars.append(dict(model_id=mid, domain_id=gid, outcome=endpoint, q_BH=q))\', \'stars.append(dict(model_id=mid, domain_id=gid, outcome=endpoint, q_BH=q, symbol=tier))\'), ("        text(gx + GRID_W, bar_y + 9.0, nominal_significance, 7.1, ha=\'right\')\\n", \'\'), ("    text(rx + 1.0, 83.8, \'*\', 11.5, True, ha=\'center\')\\n    text(rx + 4.0, 83.8, \'BH q < 0.05\', 6.7)\\n", "    for yy, symbol, threshold in ((83.8, \'*\', \'0.05\'), (88.4, \'**\', \'0.01\'), (93.0, \'***\', \'0.001\')):\\n        text(rx + .5, yy, symbol, 10.5, True, ha=\'center\')\\n        text(rx + 4.0, yy, \'BH q < \' + threshold, 6.7)\\n"), (\'heatmap_asterisk_size_pt=12.8, previous_asterisk_size_pt=7.2,\', \'heatmap_asterisk_size_pt=11.3 * .9, previous_asterisk_size_pt=12.8,\'), ("joint_forest_significance_rule=\'EXPLICIT_ORIGINAL_P_THRESHOLD_NOT_BH_ASTERISK\'", "joint_forest_significance_rule=\'SAVED_IN_DATA_NOT_DISPLAYED\'"), ("significance_key_text=\'BH q < 0.05\'", "significance_key_text=\'* BH q<0.05; ** BH q<0.01; *** BH q<0.001\'"), ("(out / \'FIGURE_LEGEND_EN.md\').write_text(", "(out / \'RENDERER_LEGEND_ORIGINAL.md\').write_text("), ("(out / \'REPRODUCE.md\').write_text(", "(out / \'REPRODUCE_ORIGINAL.md\').write_text(")]\n    for before, after in edits:\n        if before == after:\n            continue\n        assert code.count(before) == 1, (\'SOURCE_ANCHOR_DRIFT\', before)\n        code = code.replace(before, after)\n    if a.no_significance_legend:\n        significance_key = "    for yy, symbol, threshold in ((83.8, \'*\', \'0.05\'), (88.4, \'**\', \'0.01\'), (93.0, \'***\', \'0.001\')):\\n        text(rx + .5, yy, symbol, 10.5, True, ha=\'center\')\\n        text(rx + 4.0, yy, \'BH q < \' + threshold, 6.7)\\n"\n        assert code.count(significance_key) == 1\n        code = code.replace(significance_key, \'\')\n        code = code.replace("significance_key_text=\'* BH q<0.05; ** BH q<0.01; *** BH q<0.001\'", "significance_key_text=\'NOT_DISPLAYED_RESEARCHER_SUPPLIES_LEGEND\'")\n    code = code.replace(\'Figure2b_Heatmap.\', \'Figure2b_Heatmap_Multistar.\')\n    code = code.replace(\'Figure2b_Heatmap_preview\', \'Figure2b_Heatmap_Multistar_preview\')\n    code = code.replace(\'Figure2b_Heatmap_180mm_check\', \'Figure2b_Heatmap_Multistar_180mm_check\')\n    module = types.ModuleType(\'figure2b_multistar_renderer\')\n    module.__file__ = str(source)\n    exec(compile(code, str(source), \'exec\'), module.__dict__)\n    argv = sys.argv\n    try:\n        sys.argv = [str(source), \'--out\', str(out)]\n        module.main()\n    finally:\n        sys.argv = argv\n    for name in (\'PLOTDATA_HEATMAP_EFFECTS_PERCENTAGES_PQ.csv\', \'PLOTDATA_JOINT_RESULTS.csv\', \'PLOTDATA_JOINT_FOREST.csv\'):\n        assert (out / name).read_bytes() == (previous / name).read_bytes(), name\n    qc = json.loads((out / \'QC.json\').read_text())\n    assert len(qc[\'stars\']) == 14\n    assert [r[\'symbol\'] for r in qc[\'stars\']] == [\'*\'] * 14\n    assert not qc[\'text_clipped\'] and (not qc[\'graph_text_overlap\'])\n    assert pins == {p: sha(Path(p)) for p in pins}\n    ledger = [dict(model_id=r[\'model_id\'], domain_id=r[\'domain_id\'], outcome=r[\'outcome\'], point_per_1000=r[\'point_per_1000\'], original_P=r[\'P\'], original_q_BH=r[\'q_BH\'], displayed_stars=stars(r[\'q_BH\']), test_family=r[\'test_family\']) for r in data]\n    with (out / \'GRADED_STARS_AUDIT.csv\').open(\'x\', newline=\'\') as f:\n        w = csv.DictWriter(f, fieldnames=list(ledger[0]))\n        w.writeheader()\n        w.writerows(ledger)\n    (out / \'FIGURE_LEGEND_EN.md\').write_text(\'Figure 2b | Joint mediation across three analysis populations. Original-sample domain IE contributions are shown as 11-year standardized absolute risk differences per 1,000. Asterisks use each original, unrounded BH q: * q<0.05, ** q<0.01, *** q<0.001. All 14 statistically supported cells in this frozen B100 table qualify only for a single star (minimum q=0.0125); none qualifies for a double or triple star. BH correction is separate within each model×effect-type family (m=5/6 including the unshown DEP), not whole-figure FDR. Single-disease tests are supplementary exploratory tests; ΔIE tests are the original paired comparisons. The three mediated-proportion bars and joint ΔIE forest estimates and 95% CIs are unchanged. No joint BH q is available, so no BH stars are assigned to those bars or forests. All original effects, CIs, P/q and result versions are unchanged. No models, bootstrap draws or P/q calculations were rerun.\\n\')\n    (out / \'BUILD_REPORT_ZH.md\').write_text(\'# Figure 2b：多级星号规则与小字号版\\n\\n仅修改绘图，全部字体较当前图缩小10%；数字9.6→8.64 pt，星号12.8→10.17 pt。热图尺寸、行列顺序、颜色映射、领域及队列色条、比例条和Joint ΔIE森林图完全保留。\\n\\n多级规则按原始未舍入BH q：*<0.05、**<0.01、***<0.001。当前42格为14格单星、0格双星、0格三星、28格无星；最小q=0.0125。不能为了出现多星而使用原P替代q或修改阈值。\\n\\n三份统计CSV与原final_heatmap_only_v1逐字节一致；无新增拟合、Bootstrap或P/q。PNG600dpi、矢量PDF、可编辑SVG已生成，文字无裁切/重叠。\\n\')\n    (out / \'REPRODUCE.md\').write_text(\'```bash\\n\' + f\'/opt/homebrew/Caskroom/miniforge/base/bin/python3 "{Path(__file__).resolve()}"\' + (\' --no-significance-legend\' if a.no_significance_legend else \'\') + f\'\'\' --out "{out.parent / (out.name + \'_reproduced\')}"\\n```\\n\'\'\')\n    (out / \'REVISION_QC.json\').write_text(json.dumps(dict(source_sha256=pins, star_counts=tiers, minimum_original_q=min((float(r[\'q_BH\']) for r in data)), font_scale=0.9, significance_legend_displayed=not a.no_significance_legend, statistics_byte_identical=True, geometry_palette_order_unchanged=True, new_fits=0, new_bootstrap=0, new_P_q=0, text_clipped=qc[\'text_clipped\'], text_overlap=qc[\'graph_text_overlap\']), ensure_ascii=False, indent=2) + \'\\n\')\n    with (out / Path(__file__).name).open(\'xb\') as f:\n        f.write(Path(__file__).read_bytes())\n    (out / \'RENDERED_SOURCE.py\').write_text(code)\n    (out / \'FINAL_MANIFEST.json\').write_text(json.dumps(dict(files=[dict(name=p.name, sha256=sha(p)) for p in sorted(out.iterdir())]), indent=2) + \'\\n\')\n    print(json.dumps(dict(output=str(out), star_counts=tiers, font_scale=0.9, statistics_unchanged=True), ensure_ascii=False))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
