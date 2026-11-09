"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key, env in [('input-root', 'SCRNA_EMBEDDING_INPUT_ROOT'), ('script-root', 'SCRNA_EMBEDDING_SCRIPT_ROOT'), ('output-root', 'SCRNA_PUBLICATION_OUTPUT_ROOT'), ('style-root', 'YY_STYLE_ROOT')]:
        p.add_argument('--' + key, default=os.getenv(env), required=not os.getenv(env))
    p.add_argument('--resume-completed-compute', action='store_true', help='Validate completed new maps and finish metadata after a verification-only interruption')
    args = p.parse_args()
    roots = {k: Path(v).expanduser().resolve() for k, v in vars(args).items() if k != 'resume_completed_compute'}
    for k, v in roots.items():
        print(k, v, flush=True)
    src, code, out, style = [roots[k] for k in ('input_root', 'script_root', 'output_root', 'style_root')]
    required = [src / 'input/GSE255566/subset_counts.h5ad', src / 'input/GSE289422/subset_counts.h5ad', src / 'results/GSE289422/coordinates.csv.gz', src / 'results/GSE289422/parameters.json', src / 'results/GSE289422/hvg.csv', src / 'input/GSE289422/all_retained_selection_ledger.csv.gz', src / 'input/GSE255566/selection_ledger.csv.gz', code / 'embedding_compute_20261102.py', code / 'embedding_draw_YY_20261102_visible_v3.R', style / 'transcriptome_figure.R']
    for f in required:
        if not f.is_file():
            raise FileNotFoundError(f'Missing frozen source: {f}')
    hashes = {str(f): {'sha256': sha(f), 'bytes': f.stat().st_size} for f in required}
    for folder in ('results', 'metadata', 'figures'):
        (out / folder).mkdir(parents=True, exist_ok=True)
    if (out / 'metadata/PREPARATION_COMPLETE.json').exists():
        raise FileExistsError('Refusing to overwrite completed preparation')
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scipy import sparse
    spec = importlib.util.spec_from_file_location('frozen_embedding_compute', code / 'embedding_compute_20261102.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cad = pd.read_csv(src / 'results/GSE289422/coordinates.csv.gz')
    cad_param = json.loads((src / 'results/GSE289422/parameters.json').read_text())
    assert len(cad) == 44000 and cad.groupby('sample').size().eq(1000).all()
    assert cad.groupby('group').size().to_dict() == {'CAD': 29000, 'Control': 15000}
    assert cad_param['hvg']['excluded_before_ranking'] == list(mod.FIXED_GENES)
    hvg = pd.read_csv(src / 'results/GSE289422/hvg.csv')
    assert not hvg.loc[hvg.highly_variable, 'gene'].isin(mod.FIXED_GENES).any()
    dest = out / 'results/GSE289422/coordinates.csv.gz'
    dest.parent.mkdir(exist_ok=True)
    if not args.resume_completed_compute:
        with dest.open('xb') as f:
            f.write((src / 'results/GSE289422/coordinates.csv.gz').read_bytes())
    assert sha(dest) == hashes[str(src / 'results/GSE289422/coordinates.csv.gz')]['sha256']
    a = ad.read_h5ad(src / 'input/GSE255566/subset_counts.h5ad')
    keep_libraries = {'Con1', 'Con2', 'Con3', 'Mod1', 'Mod2', 'Mod3'}
    excluded = a.obs.loc[~a.obs['sample'].isin(keep_libraries)].copy()
    assert set(excluded['sample'].astype(str)) == {'treat2', 'treat3'}
    assert set(excluded['group'].astype(str)) == {'Treatment'}
    if not args.resume_completed_compute:
        mod.write_csv_new(excluded.reset_index(drop=True), out / 'metadata/T2D_excluded_treatment_cells.csv.gz')
    a = a[a.obs['sample'].isin(keep_libraries)].copy()
    assert a.n_obs == 6000 and a.obs.groupby('sample', observed=True).size().eq(1000).all()
    assert set(a.obs['group'].astype(str)) == {'Control', 'T2D'}
    assert sparse.issparse(a.X) and np.equal(a.X.data, np.floor(a.X.data)).all()
    filtered = out / 'results/GSE255566_untreated_subset_counts.h5ad'
    if not args.resume_completed_compute:
        if filtered.exists():
            raise FileExistsError(filtered)
        a.write_h5ad(filtered, compression='gzip')
    else:
        restored = ad.read_h5ad(filtered)
        assert (a.X != restored.X).nnz == 0
        pd.testing.assert_frame_equal(a.obs, restored.obs)
        pd.testing.assert_frame_equal(a.var, restored.var)
    if filtered.stat().st_size >= 500000000:
        raise ValueError('File size limit exceeded')
    if not args.resume_completed_compute:
        t2dparam = mod.compute_embeddings(filtered, 'GSE255566', out / 'results/GSE255566', n_jobs=2)
    else:
        done = json.loads((out / 'results/GSE255566/COMPLETE.json').read_text())
        assert done['completed'] and done['input_sha256'] == sha(filtered)
        for name, info in done['outputs'].items():
            assert sha(out / 'results/GSE255566' / name) == info['sha256']
        t2dparam = json.loads((out / 'results/GSE255566/parameters.json').read_text())
    t2d = pd.read_csv(out / 'results/GSE255566/coordinates.csv.gz')
    for column in a.obs.columns:
        expected = a.obs[column].reset_index(drop=True)
        actual = t2d[column].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(expected):
            np.testing.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-12, equal_nan=True)
        else:
            assert actual.fillna('').astype(str).eq(expected.astype(object).fillna('').astype(str)).all(), column
    for ds, d in [('GSE289422', cad), ('GSE255566', t2d)]:
        counts = d.groupby(['sample', 'group', 'cell_type'], observed=True).size().reset_index(name='displayed_cells')
        counts.insert(0, 'dataset', ds)
        mod.write_csv_new(counts, out / f'metadata/{ds}_display_counts.csv')
    for f in required:
        if sha(f) != hashes[str(f)]['sha256']:
            raise ValueError(f'Source changed: {f}')
    mod.write_json_new(out / 'metadata/SOURCE_INPUT_HASHES.json', hashes)
    summary = {'cad_coordinate_reuse': str(src / 'results/GSE289422/coordinates.csv.gz'), 'cad_display_cells': 44000, 'cad_donors': {'Control': 15, 'CAD': 29}, 't2d_display_cells': 6000, 't2d_nominal_donors': {'Control': 3, 'T2D': 3}, 't2d_libraries': sorted(keep_libraries), 't2d_excluded_libraries': ['treat2', 'treat3'], 't2d_excluded_cells': len(excluded), 'sampling': 'Reuse original uniform without replacement <=1000 retained cells/library, sorted cell_id, seed 20261102. No additional sampling, expression selection, annotation, QC, or tests.', 'embedding': 'Independent within-cohort maps; no cross-cohort integration. T2D was recomputed after excluding treatment libraries.', 'candidate_genes_excluded_before_HVG': list(mod.FIXED_GENES), 'T2D_parameters': t2dparam, 'CAD_parameters': cad_param, 'source_inputs_unchanged': True, 'new_n_jobs': 2, 'completed': True}
    mod.write_json_new(out / 'metadata/PREPARATION_COMPLETE.json', summary)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Reuse frozen CAD maps and recompute untreated T2D maps without resampling.\n\nAll source inputs are read only. No QC, annotation, testing, or integration.\nUse the original bounded transcriptome embedding implementation and parameters.\n"""\nimport argparse\nimport hashlib\nimport importlib.util\nimport json\nimport os\nfrom pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\n\ndef sha(path):\n    h = hashlib.sha256()\n    with Path(path).open(\'rb\') as f:\n        for block in iter(lambda: f.read(1048576), b\'\'):\n            h.update(block)\n    return h.hexdigest()\n\ndef main():\n    p = argparse.ArgumentParser(description=__doc__)\n    for key, env in [(\'input-root\', \'SCRNA_EMBEDDING_INPUT_ROOT\'), (\'script-root\', \'SCRNA_EMBEDDING_SCRIPT_ROOT\'), (\'output-root\', \'SCRNA_PUBLICATION_OUTPUT_ROOT\'), (\'style-root\', \'YY_STYLE_ROOT\')]:\n        p.add_argument(\'--\' + key, default=os.getenv(env), required=not os.getenv(env))\n    p.add_argument(\'--resume-completed-compute\', action=\'store_true\', help=\'Validate completed new maps and finish metadata after a verification-only interruption\')\n    args = p.parse_args()\n    roots = {k: Path(v).expanduser().resolve() for k, v in vars(args).items() if k != \'resume_completed_compute\'}\n    for k, v in roots.items():\n        print(k, v, flush=True)\n    src, code, out, style = [roots[k] for k in (\'input_root\', \'script_root\', \'output_root\', \'style_root\')]\n    required = [src / \'input/GSE255566/subset_counts.h5ad\', src / \'input/GSE289422/subset_counts.h5ad\', src / \'results/GSE289422/coordinates.csv.gz\', src / \'results/GSE289422/parameters.json\', src / \'results/GSE289422/hvg.csv\', src / \'input/GSE289422/all_retained_selection_ledger.csv.gz\', src / \'input/GSE255566/selection_ledger.csv.gz\', code / \'embedding_compute_20261102.py\', code / \'embedding_draw_YY_20261102_visible_v3.R\', style / \'transcriptome_figure.R\']\n    for f in required:\n        if not f.is_file():\n            raise FileNotFoundError(f\'Missing frozen source: {f}\')\n    hashes = {str(f): {\'sha256\': sha(f), \'bytes\': f.stat().st_size} for f in required}\n    for folder in (\'results\', \'metadata\', \'figures\'):\n        (out / folder).mkdir(parents=True, exist_ok=True)\n    if (out / \'metadata/PREPARATION_COMPLETE.json\').exists():\n        raise FileExistsError(\'Refusing to overwrite completed preparation\')\n    import numpy as np\n    import pandas as pd\n    import anndata as ad\n    from scipy import sparse\n    spec = importlib.util.spec_from_file_location(\'frozen_embedding_compute\', code / \'embedding_compute_20261102.py\')\n    mod = importlib.util.module_from_spec(spec)\n    spec.loader.exec_module(mod)\n    cad = pd.read_csv(src / \'results/GSE289422/coordinates.csv.gz\')\n    cad_param = json.loads((src / \'results/GSE289422/parameters.json\').read_text())\n    assert len(cad) == 44000 and cad.groupby(\'sample\').size().eq(1000).all()\n    assert cad.groupby(\'group\').size().to_dict() == {\'CAD\': 29000, \'Control\': 15000}\n    assert cad_param[\'hvg\'][\'excluded_before_ranking\'] == list(mod.FIXED_GENES)\n    hvg = pd.read_csv(src / \'results/GSE289422/hvg.csv\')\n    assert not hvg.loc[hvg.highly_variable, \'gene\'].isin(mod.FIXED_GENES).any()\n    dest = out / \'results/GSE289422/coordinates.csv.gz\'\n    dest.parent.mkdir(exist_ok=True)\n    if not args.resume_completed_compute:\n        with dest.open(\'xb\') as f:\n            f.write((src / \'results/GSE289422/coordinates.csv.gz\').read_bytes())\n    assert sha(dest) == hashes[str(src / \'results/GSE289422/coordinates.csv.gz\')][\'sha256\']\n    a = ad.read_h5ad(src / \'input/GSE255566/subset_counts.h5ad\')\n    keep_libraries = {\'Con1\', \'Con2\', \'Con3\', \'Mod1\', \'Mod2\', \'Mod3\'}\n    excluded = a.obs.loc[~a.obs[\'sample\'].isin(keep_libraries)].copy()\n    assert set(excluded[\'sample\'].astype(str)) == {\'treat2\', \'treat3\'}\n    assert set(excluded[\'group\'].astype(str)) == {\'Treatment\'}\n    if not args.resume_completed_compute:\n        mod.write_csv_new(excluded.reset_index(drop=True), out / \'metadata/T2D_excluded_treatment_cells.csv.gz\')\n    a = a[a.obs[\'sample\'].isin(keep_libraries)].copy()\n    assert a.n_obs == 6000 and a.obs.groupby(\'sample\', observed=True).size().eq(1000).all()\n    assert set(a.obs[\'group\'].astype(str)) == {\'Control\', \'T2D\'}\n    assert sparse.issparse(a.X) and np.equal(a.X.data, np.floor(a.X.data)).all()\n    filtered = out / \'results/GSE255566_untreated_subset_counts.h5ad\'\n    if not args.resume_completed_compute:\n        if filtered.exists():\n            raise FileExistsError(filtered)\n        a.write_h5ad(filtered, compression=\'gzip\')\n    else:\n        restored = ad.read_h5ad(filtered)\n        assert (a.X != restored.X).nnz == 0\n        pd.testing.assert_frame_equal(a.obs, restored.obs)\n        pd.testing.assert_frame_equal(a.var, restored.var)\n    if filtered.stat().st_size >= 500000000:\n        raise ValueError(\'File size limit exceeded\')\n    if not args.resume_completed_compute:\n        t2dparam = mod.compute_embeddings(filtered, \'GSE255566\', out / \'results/GSE255566\', n_jobs=2)\n    else:\n        done = json.loads((out / \'results/GSE255566/COMPLETE.json\').read_text())\n        assert done[\'completed\'] and done[\'input_sha256\'] == sha(filtered)\n        for name, info in done[\'outputs\'].items():\n            assert sha(out / \'results/GSE255566\' / name) == info[\'sha256\']\n        t2dparam = json.loads((out / \'results/GSE255566/parameters.json\').read_text())\n    t2d = pd.read_csv(out / \'results/GSE255566/coordinates.csv.gz\')\n    for column in a.obs.columns:\n        expected = a.obs[column].reset_index(drop=True)\n        actual = t2d[column].reset_index(drop=True)\n        if pd.api.types.is_numeric_dtype(expected):\n            np.testing.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-12, equal_nan=True)\n        else:\n            assert actual.fillna(\'\').astype(str).eq(expected.astype(object).fillna(\'\').astype(str)).all(), column\n    for ds, d in [(\'GSE289422\', cad), (\'GSE255566\', t2d)]:\n        counts = d.groupby([\'sample\', \'group\', \'cell_type\'], observed=True).size().reset_index(name=\'displayed_cells\')\n        counts.insert(0, \'dataset\', ds)\n        mod.write_csv_new(counts, out / f\'metadata/{ds}_display_counts.csv\')\n    for f in required:\n        if sha(f) != hashes[str(f)][\'sha256\']:\n            raise ValueError(f\'Source changed: {f}\')\n    mod.write_json_new(out / \'metadata/SOURCE_INPUT_HASHES.json\', hashes)\n    summary = {\'cad_coordinate_reuse\': str(src / \'results/GSE289422/coordinates.csv.gz\'), \'cad_display_cells\': 44000, \'cad_donors\': {\'Control\': 15, \'CAD\': 29}, \'t2d_display_cells\': 6000, \'t2d_nominal_donors\': {\'Control\': 3, \'T2D\': 3}, \'t2d_libraries\': sorted(keep_libraries), \'t2d_excluded_libraries\': [\'treat2\', \'treat3\'], \'t2d_excluded_cells\': len(excluded), \'sampling\': \'Reuse original uniform without replacement <=1000 retained cells/library, sorted cell_id, seed 20261102. No additional sampling, expression selection, annotation, QC, or tests.\', \'embedding\': \'Independent within-cohort maps; no cross-cohort integration. T2D was recomputed after excluding treatment libraries.\', \'candidate_genes_excluded_before_HVG\': list(mod.FIXED_GENES), \'T2D_parameters\': t2dparam, \'CAD_parameters\': cad_param, \'source_inputs_unchanged\': True, \'new_n_jobs\': 2, \'completed\': True}\n    mod.write_json_new(out / \'metadata/PREPARATION_COMPLETE.json\', summary)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
