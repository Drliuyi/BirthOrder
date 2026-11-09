"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rows(path):
    return list(csv.DictReader(io.StringIO(path.read_text(encoding='utf-8-sig'))))

def main():
    project = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[3]))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(os.environ.get('PHEWAS_SOURCE_ROOT', project / 'reports/PheWAS_S8_relative_main_six_labels_20261030_v2')))
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'))
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError(f'Refuse existing output: {out}')
    print(f'INPUT_ROOT={source}\nSCRIPT_ROOT={Path(__file__).parent}\nOUTPUT_ROOT={out}', flush=True)
    required = [source / 'config.json', source / 'manuscript_style.json', source / 'code/plot.R']
    required += sorted((source / 'data/S8_relative').iterdir())
    if not all((path.is_file() for path in required)):
        raise FileNotFoundError('Missing frozen source or unexpected non-file plot input')
    bindings = [{'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in required]
    before = rows(source / 'data/S8_relative/disease_labels.csv')
    assert len(before) == 6
    assert [r['label'] for r in before if r['phecode'] == '411.4'] == ['Cor. atheroscl.']
    for folder in ['code', 'data/S8_relative', 'figures', 'logs', 'qa']:
        (out / folder).mkdir(parents=True, exist_ok=True)
    for path in required:
        target = out / path.relative_to(source)
        shutil.copy2(path, target)
    shutil.copy2(__file__, out / 'code/render.py')
    labels_path = out / 'data/S8_relative/disease_labels.csv'
    after = rows(labels_path)
    for row in after:
        if row['phecode'] == '411.4':
            row['label'] = 'CAS'
    labels_path.unlink()
    with labels_path.open('x', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(after[0]))
        writer.writeheader()
        writer.writerows(after)
    style_path = out / 'manuscript_style.json'
    style = json.loads(style_path.read_text())
    style['disease_short']['411.4'] = 'CAS'
    style_path.write_text(json.dumps(style, indent=2, ensure_ascii=False) + '\n')
    assert all(({k: v for k, v in a.items() if k != 'label'} == {k: v for k, v in b.items() if k != 'label'} for a, b in zip(before, after)))
    result = subprocess.run([args.rscript, str(out / 'code/plot.R'), str(out)], text=True, capture_output=True)
    (out / 'logs/PLOTTING_LOG.txt').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    stem = 'S8_relative_main_signed_PheWAS'
    old_svg = source / f'figures/S8_relative/{stem}.svg'
    new_svg = out / f'figures/S8_relative/{stem}.svg'
    old_xml, new_xml = (ET.parse(old_svg).getroot(), ET.parse(new_svg).getroot())
    ns = {'s': 'http://www.w3.org/2000/svg'}
    old_text = [''.join(node.itertext()) for node in old_xml.findall('.//s:text', ns)]
    new_text = [''.join(node.itertext()) for node in new_xml.findall('.//s:text', ns)]
    assert new_text == ['CAS' if value == 'Cor. atheroscl.' else value for value in old_text]

    def geometry(root):
        return [(node.tag, dict(node.attrib)) for node in root.iter() if node.tag.rsplit('}', 1)[-1] in {'circle', 'polyline', 'polygon', 'path', 'rect', 'line'}]
    assert geometry(old_xml) == geometry(new_xml), 'Unexpected non-text graphical change'
    for name in ['plot_points.csv', 'category_layout.csv', 'thresholds.json', 'source_branch_1582.csv']:
        assert sha(source / 'data/S8_relative' / name) == sha(out / 'data/S8_relative' / name)
    for binding in bindings:
        assert sha(Path(binding['path'])) == binding['sha256']
    audit = {'status': 'PASS', 'change': 'Cor. atheroscl. -> CAS', 'phecode': '411.4', 'full_name': 'Coronary atherosclerosis', 'all_other_labels_unchanged': True, 'non_text_svg_geometry_identical': True, 'statistical_inputs_byte_identical': True, 'models_fitted': 0, 'p_q_recomputed': False, 'source_preserved': True, 'source_files': bindings, 'outputs': [{'path': str(path.relative_to(out)), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in sorted((out / 'figures/S8_relative').iterdir()) if path.is_file()]}
    (out / 'qa/LABEL_REVISION_AUDIT.json').write_text(json.dumps(audit, indent=2) + '\n')
    (out / 'README.md').write_text(f"# Label-only revision\n\nCor. atheroscl. is now CAS (Coronary atherosclerosis; phecode 411.4). All other labels, points, thresholds, axes, colours and statistical inputs are unchanged. The original descriptive/pending-review status is retained. No models or corrections rerun.\n\nSource: {source}\n\nReproduce into a new directory:\n\n`python3 code/render.py --source '{source}' --out /absolute/path/to/NEW_DIRECTORY --rscript Rscript`\n")
    print(json.dumps({'status': 'PASS', 'output': str(out), 'models_fitted': 0}))

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Additive label-only revision; no modelling or numerical-result changes."""\nimport argparse\nimport csv\nimport hashlib\nimport io\nimport json\nimport os\nfrom pathlib import Path\nimport shutil\nimport subprocess\nimport sys\nimport xml.etree.ElementTree as ET\n\ndef sha(path):\n    return hashlib.sha256(path.read_bytes()).hexdigest()\n\ndef rows(path):\n    return list(csv.DictReader(io.StringIO(path.read_text(encoding=\'utf-8-sig\'))))\n\ndef main():\n    project = Path(os.environ.get(\'PROJECT_ROOT\', Path(__file__).resolve().parents[3]))\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--source\', type=Path, default=Path(os.environ.get(\'PHEWAS_SOURCE_ROOT\', project / \'reports/PheWAS_S8_relative_main_six_labels_20261030_v2\')))\n    parser.add_argument(\'--out\', type=Path, required=True)\n    parser.add_argument(\'--rscript\', default=os.environ.get(\'RSCRIPT\', \'Rscript\'))\n    args = parser.parse_args()\n    source = args.source.resolve(strict=True)\n    out = args.out.resolve()\n    if out.exists():\n        raise FileExistsError(f\'Refuse existing output: {out}\')\n    print(f\'INPUT_ROOT={source}\\nSCRIPT_ROOT={Path(__file__).parent}\\nOUTPUT_ROOT={out}\', flush=True)\n    required = [source / \'config.json\', source / \'manuscript_style.json\', source / \'code/plot.R\']\n    required += sorted((source / \'data/S8_relative\').iterdir())\n    if not all((path.is_file() for path in required)):\n        raise FileNotFoundError(\'Missing frozen source or unexpected non-file plot input\')\n    bindings = [{\'path\': str(path), \'bytes\': path.stat().st_size, \'sha256\': sha(path)} for path in required]\n    before = rows(source / \'data/S8_relative/disease_labels.csv\')\n    assert len(before) == 6\n    assert [r[\'label\'] for r in before if r[\'phecode\'] == \'411.4\'] == [\'Cor. atheroscl.\']\n    for folder in [\'code\', \'data/S8_relative\', \'figures\', \'logs\', \'qa\']:\n        (out / folder).mkdir(parents=True, exist_ok=True)\n    for path in required:\n        target = out / path.relative_to(source)\n        shutil.copy2(path, target)\n    shutil.copy2(__file__, out / \'code/render.py\')\n    labels_path = out / \'data/S8_relative/disease_labels.csv\'\n    after = rows(labels_path)\n    for row in after:\n        if row[\'phecode\'] == \'411.4\':\n            row[\'label\'] = \'CAS\'\n    labels_path.unlink()\n    with labels_path.open(\'x\', encoding=\'utf-8-sig\', newline=\'\') as handle:\n        writer = csv.DictWriter(handle, fieldnames=list(after[0]))\n        writer.writeheader()\n        writer.writerows(after)\n    style_path = out / \'manuscript_style.json\'\n    style = json.loads(style_path.read_text())\n    style[\'disease_short\'][\'411.4\'] = \'CAS\'\n    style_path.write_text(json.dumps(style, indent=2, ensure_ascii=False) + \'\\n\')\n    assert all(({k: v for k, v in a.items() if k != \'label\'} == {k: v for k, v in b.items() if k != \'label\'} for a, b in zip(before, after)))\n    result = subprocess.run([args.rscript, str(out / \'code/plot.R\'), str(out)], text=True, capture_output=True)\n    (out / \'logs/PLOTTING_LOG.txt\').write_text(result.stdout + result.stderr)\n    if result.returncode:\n        raise RuntimeError(result.stdout + result.stderr)\n    stem = \'S8_relative_main_signed_PheWAS\'\n    old_svg = source / f\'figures/S8_relative/{stem}.svg\'\n    new_svg = out / f\'figures/S8_relative/{stem}.svg\'\n    old_xml, new_xml = (ET.parse(old_svg).getroot(), ET.parse(new_svg).getroot())\n    ns = {\'s\': \'http://www.w3.org/2000/svg\'}\n    old_text = [\'\'.join(node.itertext()) for node in old_xml.findall(\'.//s:text\', ns)]\n    new_text = [\'\'.join(node.itertext()) for node in new_xml.findall(\'.//s:text\', ns)]\n    assert new_text == [\'CAS\' if value == \'Cor. atheroscl.\' else value for value in old_text]\n\n    def geometry(root):\n        return [(node.tag, dict(node.attrib)) for node in root.iter() if node.tag.rsplit(\'}\', 1)[-1] in {\'circle\', \'polyline\', \'polygon\', \'path\', \'rect\', \'line\'}]\n    assert geometry(old_xml) == geometry(new_xml), \'Unexpected non-text graphical change\'\n    for name in [\'plot_points.csv\', \'category_layout.csv\', \'thresholds.json\', \'source_branch_1582.csv\']:\n        assert sha(source / \'data/S8_relative\' / name) == sha(out / \'data/S8_relative\' / name)\n    for binding in bindings:\n        assert sha(Path(binding[\'path\'])) == binding[\'sha256\']\n    audit = {\'status\': \'PASS\', \'change\': \'Cor. atheroscl. -> CAS\', \'phecode\': \'411.4\', \'full_name\': \'Coronary atherosclerosis\', \'all_other_labels_unchanged\': True, \'non_text_svg_geometry_identical\': True, \'statistical_inputs_byte_identical\': True, \'models_fitted\': 0, \'p_q_recomputed\': False, \'source_preserved\': True, \'source_files\': bindings, \'outputs\': [{\'path\': str(path.relative_to(out)), \'bytes\': path.stat().st_size, \'sha256\': sha(path)} for path in sorted((out / \'figures/S8_relative\').iterdir()) if path.is_file()]}\n    (out / \'qa/LABEL_REVISION_AUDIT.json\').write_text(json.dumps(audit, indent=2) + \'\\n\')\n    (out / \'README.md\').write_text(f"# Label-only revision\\n\\nCor. atheroscl. is now CAS (Coronary atherosclerosis; phecode 411.4). All other labels, points, thresholds, axes, colours and statistical inputs are unchanged. The original descriptive/pending-review status is retained. No models or corrections rerun.\\n\\nSource: {source}\\n\\nReproduce into a new directory:\\n\\n`python3 code/render.py --source \'{source}\' --out /absolute/path/to/NEW_DIRECTORY --rscript Rscript`\\n")\n    print(json.dumps({\'status\': \'PASS\', \'output\': str(out), \'models_fitted\': 0}))\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
