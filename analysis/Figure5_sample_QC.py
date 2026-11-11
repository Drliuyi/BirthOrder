"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
import gc
import gzip
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import re
import tarfile
import time
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
FIXED_GENES = ('FCRL6', 'GZMA', 'GZMH', 'KLRD1')
CELL_TYPES = ('NK', 'CD8_T', 'CD4_T', 'T_other', 'B', 'Plasma', 'Monocyte', 'DC', 'pDC', 'Platelet', 'Erythroid', 'Other/Unassigned')
SEED = 20261101
MARKERS = ('CD3D', 'CD3E', 'TRAC', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'CCR7', 'NCAM1', 'NCR1', 'FCGR3A', 'TYROBP', 'GNLY', 'PRF1', 'NKG7', 'KLRF1', 'MS4A1', 'CD79A', 'CD79B', 'CD19', 'JCHAIN', 'MZB1', 'SDC1', 'CD14', 'LYZ', 'S100A8', 'S100A9', 'FCER1A', 'CLEC10A', 'CLEC9A', 'GZMB', 'TCF4', 'IL3RA', 'PPBP', 'PF4', 'HBA1', 'HBB')

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def write_json_new(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, default=str)
        f.write('\n')

def write_csv_new(df, path, sep=','):
    path = Path(path)
    if path.suffix == '.gz':
        with path.open('xb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as gz:
                with io.TextIOWrapper(gz) as f:
                    df.to_csv(f, index=False, sep=sep)
    else:
        with path.open('x') as f:
            df.to_csv(f, index=False, sep=sep)

def write_mtx_new(matrix, path):
    with Path(path).open('xb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as gz:
            scipy.io.mmwrite(gz, matrix, field='integer', symmetry='general')

def broad_type(label):
    """Reference-label mapping fixed before seeing disease effects."""
    t = str(label).lower()
    if t in ('heterogeneous', 'unassigned', 'nan', ''):
        return 'Other/Unassigned'
    if 'mast' in t:
        return 'Other/Unassigned'
    if 'plasma' in t and 'plasmacytoid' not in t:
        return 'Plasma'
    if 'pdc' in t or 'plasmacytoid' in t:
        return 'pDC'
    if 'erythro' in t and 'macrophage' not in t:
        return 'Erythroid'
    if 'platelet' in t or 'megakaryo' in t or t == 'early mk':
        return 'Platelet'
    if 'nk' in t and 'nkt' not in t:
        return 'NK'
    if any((v in t for v in ('monocyt', 'macrophage', 'mono-mac', 'kupffer', 'hofbauer'))):
        return 'Monocyte'
    if re.search('(^|\\W)dc[123]?(\\W|$)', t) or 'dendritic' in t:
        return 'DC'
    if 'b cell' in t or '-b cell' in t:
        return 'B'
    if 'cd4' in t or 'helper t' in t or 'regulatory t' in t or ('treg' in t):
        return 'CD4_T'
    if 'cd8' in t or 'cytotoxic t' in t:
        return 'CD8_T'
    if re.search('(^|[ -])t cells?($|\\W)', t) or any((v in t for v in ('nkt', 'mait', 'thymocyte', 't(agonist)', 't lymphoid'))):
        return 'T_other'
    return 'Other/Unassigned'

def _dense_gene_values(counts, symbols, targets):
    """Sum duplicate symbols only for display; retain original features for PB."""
    return np.column_stack([np.asarray(counts[:, np.flatnonzero(symbols == g)].sum(axis=1)).ravel() if np.any(symbols == g) else np.full(counts.shape[0], np.nan) for g in targets])

def _concat_csv(inputs, output):
    """Bound memory by copying CSV records, one file at a time."""
    output = Path(output)
    with output.open('xb') as raw:
        if output.suffix == '.gz':
            holder = gzip.GzipFile(fileobj=raw, mode='wb', mtime=0)
        else:
            holder = raw
        try:
            for i, p in enumerate(inputs):
                opener = gzip.open if str(p).endswith('.gz') else open
                with opener(p, 'rb') as src:
                    header = src.readline()
                    if i == 0:
                        holder.write(header)
                    for b in iter(lambda: src.read(1024 * 1024), b''):
                        holder.write(b)
        finally:
            if holder is not raw:
                holder.close()

def process_sample(counts, genes, barcodes, *, sample, group, outdir, model_path, dataset, source=None, donor_status='nominal_sample'):
    """Process one cells×features CSR matrix; return completed JSON summary.

    genes: DataFrame with gene_id,gene in count-column order.
    source: optional donor_id,sex,age,batch plus immutable source provenance.
    Analyses retain QC-passed predicted singlets; failed automatic Scrublet
    estimation retains QC-passed cells with the failure explicitly recorded.
    """
    started = time.time()
    source = dict(source or {})
    sample = str(sample)
    if re.search('[^A-Za-z0-9_.-]', sample):
        raise ValueError('sample must be a safe single path component')
    root = Path(outdir)
    sdir = root / 'samples' / sample
    complete = sdir / 'COMPLETE_v1.json'
    if complete.exists():
        old = json.loads(complete.read_text())
        if old['source'] != source or old['group'] != group or old['dataset'] != dataset:
            raise ValueError(f'Completed source identity mismatch: {sdir}')
        print(f'{sample}: verified completed sample, reusing immutable checkpoint', flush=True)
        return old
    sdir.mkdir(parents=True, exist_ok=False)
    print(f'{sample}: starting {counts.shape}', flush=True)
    genes = genes.loc[:, ['gene_id', 'gene']].astype(str).reset_index(drop=True)
    counts = sp.csr_matrix(counts)
    counts.sum_duplicates()
    counts.eliminate_zeros()
    if counts.shape != (len(barcodes), len(genes)):
        raise ValueError('Count dimensions do not match barcodes/features')
    if len(set(map(str, barcodes))) != len(barcodes):
        raise ValueError('Duplicate source barcodes within a sample')
    valid = np.all(np.isfinite(counts.data)) and np.all(counts.data >= 0)
    integer = np.all(counts.data == np.floor(counts.data))
    if not (valid and integer):
        raise ValueError('Source matrix is not nonnegative finite integer counts')
    counts = counts.astype(np.int64, copy=False)
    symbols = genes.gene.to_numpy()
    total = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.diff(counts.indptr)
    is_mt = np.char.startswith(symbols.astype(str), 'MT-')
    mt = np.asarray(counts[:, is_mt].sum(axis=1)).ravel()
    mt_pct = np.divide(mt * 100.0, total, out=np.zeros(len(total)), where=total > 0)
    qc = (n_genes >= QC['min_genes']) & (n_genes <= QC['max_genes']) & (total >= QC['min_UMI']) & (total <= QC['max_UMI']) & (mt_pct <= QC['max_mito_percent'])
    meta = pd.DataFrame({'cell_id': [sample + '__' + str(b) for b in barcodes], 'barcode': list(map(str, barcodes)), 'dataset': dataset, 'sample': sample, 'group': group, 'donor_id': str(source.get('donor_id', sample)), 'donor_status': donor_status, 'sex': str(source.get('sex', 'unknown')), 'age': str(source.get('age', 'unknown')), 'batch': str(source.get('batch', 'unknown')), 'n_genes': n_genes, 'total_counts': total, 'pct_mt': mt_pct, 'qc_pass': qc, 'doublet_score': np.nan, 'predicted_doublet': False, 'retained': False, 'cell_type': 'Failed_QC', 'raw_CellTypist_label': 'not_annotated', 'majority_CellTypist_label': 'not_annotated', 'pred_confidence': np.nan, 'majority_label_probability': np.nan, 'over_clustering': 'not_annotated'})
    qix = np.flatnonzero(qc)
    scrub_info = {'expected_rate': 0.06, 'seed': SEED, 'status': 'not_run_no_qc_cells'}
    doublets = np.zeros(len(qix), dtype=bool)
    if len(qix):
        import scanpy as sc
        import anndata as ad
        scrub = ad.AnnData(counts[qix, :].copy())
        try:
            sc.pp.scrublet(scrub, expected_doublet_rate=0.06, random_state=SEED, n_prin_comps=30, verbose=True)
            scores = scrub.obs['doublet_score'].to_numpy()
            calls = scrub.obs.get('predicted_doublet')
            meta.loc[qix, 'doublet_score'] = scores
            threshold = scrub.uns.get('scrublet', {}).get('threshold')
            if calls is not None and threshold is not None and np.isfinite(threshold):
                doublets = np.asarray(calls, dtype=bool)
                scrub_info.update(status='automatic_threshold', threshold=float(threshold), predicted_doublets=int(doublets.sum()))
            else:
                scrub_info.update(status='automatic_threshold_failed_all_qc_retained', threshold=None, predicted_doublets=0)
        except Exception as e:
            if 'doublet_score' in scrub.obs:
                meta.loc[qix, 'doublet_score'] = scrub.obs['doublet_score'].to_numpy()
            scrub_info.update(status='not_estimable_all_qc_retained', error=repr(e), threshold=None, predicted_doublets=0)
        scrub_info['implementation'] = 'scanpy.pp.scrublet'
        del scrub
        gc.collect()
    meta.loc[qix, 'predicted_doublet'] = doublets
    meta.loc[qix, 'cell_type'] = 'Excluded_doublet'
    keep = qix[~doublets]
    meta.loc[keep, 'retained'] = True
    model_hash = sha256_file(model_path)
    model_info = {'path': str(Path(model_path).resolve()), 'sha256': model_hash, 'fixed_genes_excluded_from_annotation_and_overclustering': list(FIXED_GENES)}
    if len(keep):
        import scanpy as sc
        import anndata as ad
        import celltypist
        np.random.seed(SEED)
        model = celltypist.models.Model.load(str(model_path))
        excluded = np.isin(np.asarray(model.classifier.features), FIXED_GENES)
        model.classifier.coef_[:, excluded] = 0
        model_info['zeroed_model_features'] = np.asarray(model.classifier.features)[excluded].tolist()
        use = ~np.isin(symbols, FIXED_GENES)
        ann_counts = counts[keep, :][:, use].astype(np.float32)
        ann_symbols = symbols[use]
        codes, unique_symbols = pd.factorize(ann_symbols, sort=False)
        if len(unique_symbols) != len(ann_symbols):
            mapping = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(unique_symbols)))
            ann_counts = (ann_counts @ mapping).astype(np.float32).tocsr()
        annotation = ad.AnnData(ann_counts, obs=pd.DataFrame(index=meta.loc[keep, 'cell_id'].values), var=pd.DataFrame(index=pd.Index(unique_symbols)))
        sc.pp.normalize_total(annotation, target_sum=10000.0)
        sc.pp.log1p(annotation)
        pred = celltypist.annotate(annotation, model=model, majority_voting=True)
        labels = pred.predicted_labels
        rawlabel = labels['predicted_labels'].astype(str).to_numpy()
        majority = labels.get('majority_voting', labels['predicted_labels']).astype(str).to_numpy()
        probability = pred.probability_matrix
        conf = probability.max(axis=1).to_numpy()
        assigned_prob = np.array([probability.iloc[i].get(t, np.nan) for i, t in enumerate(majority)])
        meta.loc[keep, 'raw_CellTypist_label'] = rawlabel
        meta.loc[keep, 'majority_CellTypist_label'] = majority
        meta.loc[keep, 'pred_confidence'] = conf
        meta.loc[keep, 'majority_label_probability'] = assigned_prob
        meta.loc[keep, 'cell_type'] = [broad_type(t) for t in majority]
        if 'over_clustering' in labels:
            meta.loc[keep, 'over_clustering'] = labels['over_clustering'].astype(str).to_numpy()
        model_info['reference_labels_mapping'] = {str(t): broad_type(t) for t in model.cell_types}
        del pred, annotation, model, ann_counts
        gc.collect()
    fixed_raw = _dense_gene_values(counts, symbols, FIXED_GENES)
    fixed_log = np.log1p(np.divide(fixed_raw * 10000.0, total[:, None], out=np.zeros_like(fixed_raw, dtype=float), where=total[:, None] > 0))
    fx = meta[['cell_id', 'dataset', 'sample', 'group', 'donor_id', 'cell_type', 'qc_pass', 'predicted_doublet', 'retained']].copy()
    for j, g in enumerate(FIXED_GENES):
        fx[g + '_raw'] = fixed_raw[:, j]
        fx[g + '_logCP10k'] = fixed_log[:, j]
    write_csv_new(meta, sdir / 'cellmetadata_allinput_v1.csv.gz')
    write_csv_new(meta[qc], sdir / 'all_qc_cellmetadata_v1.csv.gz')
    write_csv_new(fx, sdir / 'fixed4_allcells_expression_v1.csv.gz')
    canonical_raw = _dense_gene_values(counts[keep, :], symbols, MARKERS) if len(keep) else np.empty((0, len(MARKERS)))
    canonical_log = np.log1p(np.divide(canonical_raw * 10000.0, total[keep, None], out=np.zeros_like(canonical_raw, dtype=float), where=total[keep, None] > 0))
    canonical_cells = meta.loc[keep, ['cell_id', 'dataset', 'sample', 'group', 'donor_id', 'cell_type', 'raw_CellTypist_label', 'majority_CellTypist_label', 'pred_confidence']].copy()
    for j, g in enumerate(MARKERS):
        canonical_cells[g + '_raw'] = canonical_raw[:, j]
        canonical_cells[g + '_logCP10k'] = canonical_log[:, j]
    write_csv_new(canonical_cells, sdir / 'canonical_markers_retainedcells_v1.csv.gz')
    ctype = meta.loc[keep, 'cell_type'].to_numpy()
    types = list(CELL_TYPES)
    base = {'dataset': dataset, 'sample': sample, 'group': group, 'donor_id': str(source.get('donor_id', sample)), 'sex': str(source.get('sex', 'unknown')), 'age': str(source.get('age', 'unknown')), 'batch': str(source.get('batch', 'unknown'))}
    pb_rows = []
    pb_meta = []
    summaries = []
    compositions = []
    marker_rows = []
    for t in types:
        mask = ctype == t
        ix = keep[mask]
        compositions.append(dict(base, cell_type=t, n_cells=int(mask.sum()), n_sample_retained=len(keep), proportion=float(mask.mean())))
        if not mask.any():
            for g in FIXED_GENES:
                summaries.append(dict(base, cell_type=t, gene=g, n_cells=0, measured=bool(np.any(symbols == g)), detected_cells=0, percent_detected=np.nan, mean_logCP10k=np.nan, mean_CP10k=np.nan, pseudobulk_count=0))
            continue
        pb_id = sample + '__' + t.replace('/', '_')
        pb_rows.append(sp.csr_matrix(counts[ix, :].sum(axis=0)))
        pb_meta.append(dict(base, pb_id=pb_id, cell_type=t, n_cells=int(mask.sum())))
        for j, g in enumerate(FIXED_GENES):
            summaries.append(dict(base, cell_type=t, gene=g, n_cells=int(mask.sum()), measured=bool(np.any(symbols == g)), detected_cells=int(np.sum(fixed_raw[ix, j] > 0)), percent_detected=float(np.mean(fixed_raw[ix, j] > 0) * 100), mean_logCP10k=float(np.mean(fixed_log[ix, j])), mean_CP10k=float(np.mean(np.expm1(fixed_log[ix, j]))), pseudobulk_count=int(np.nansum(fixed_raw[ix, j]))))
        for j, g in enumerate(MARKERS):
            marker_rows.append(dict(base, cell_type=t, gene=g, n_cells=int(mask.sum()), measured=bool(np.any(symbols == g)), percent_detected=float(np.mean(canonical_raw[mask, j] > 0) * 100), mean_logCP10k=float(np.mean(canonical_log[mask, j]))))
    if not pb_rows:
        raise ValueError(f'No retained cells in {sample}')
    pb = sp.vstack(pb_rows, format='csr').T.tocsr()
    write_mtx_new(pb, sdir / 'raw_pseudobulk_counts_v1.mtx.gz')
    write_csv_new(genes, sdir / 'raw_pseudobulk_genes_v1.tsv', sep='\t')
    write_csv_new(pd.DataFrame(pb_meta), sdir / 'raw_pseudobulk_samples_v1.tsv', sep='\t')
    write_csv_new(pd.DataFrame(summaries), sdir / 'sample_celltype_gene_summary_v1.csv')
    write_csv_new(pd.DataFrame(compositions), sdir / 'cell_composition_v1.csv')
    write_csv_new(pd.DataFrame(marker_rows), sdir / 'canonical_marker_means_v1.csv')
    summary = {'sample': sample, 'group': group, 'dataset': dataset, 'source': source, 'donor_status': donor_status, 'input_cells': counts.shape[0], 'input_genes': counts.shape[1], 'input_nnz': int(counts.nnz), 'raw_values_all_validated_nonnegative_integer': True, 'qc_pass': int(qc.sum()), 'qc_failed': int((~qc).sum()), 'retained': len(keep), 'doublets': scrub_info, 'qc': QC, 'annotation': model_info, 'cell_type_counts': meta.loc[keep, 'cell_type'].value_counts().to_dict(), 'fixed_gene_feature_rows': {g: np.flatnonzero(symbols == g).tolist() for g in FIXED_GENES}, 'elapsed_seconds': time.time() - started, 'software': {x: importlib.metadata.version(x) for x in ('numpy', 'scipy', 'pandas', 'scanpy', 'anndata', 'scrublet', 'celltypist', 'leidenalg')}}
    write_json_new(complete, summary)
    print(f"{sample}: COMPLETE retained {len(keep)}/{len(qix)} QC cells, {summary['cell_type_counts']}", flush=True)
    return summary

def finalize(outdir, *, dataset, source=None):
    """Merge only completed checkpoints, preserving one raw feature ordering."""
    root = Path(outdir)
    done = sorted(root.glob('samples/*/COMPLETE_v1.json'))
    if not done:
        raise ValueError('No completed sample checkpoints')
    if (root / 'FINAL_COMPLETE_v1.json').exists():
        old = json.loads((root / 'FINAL_COMPLETE_v1.json').read_text())
        if old['sample_names'] != [p.parent.name for p in done]:
            raise ValueError('Final checkpoint sample set mismatch')
        return old
    summaries = [json.loads(p.read_text()) for p in done]
    if any((s['dataset'] != dataset for s in summaries)):
        raise ValueError('Mixed datasets')
    dirs = [p.parent for p in done]
    genes = pd.read_csv(dirs[0] / 'raw_pseudobulk_genes_v1.tsv', sep='\t', dtype=str, keep_default_na=False)
    matrices = []
    pbmeta = []
    for d in dirs:
        g = pd.read_csv(d / 'raw_pseudobulk_genes_v1.tsv', sep='\t', dtype=str, keep_default_na=False)
        if not genes.equals(g):
            raise ValueError(f'Feature order mismatch: {d}')
        matrices.append(scipy.io.mmread(d / 'raw_pseudobulk_counts_v1.mtx.gz').tocsr())
        pbmeta.append(pd.read_csv(d / 'raw_pseudobulk_samples_v1.tsv', sep='\t', keep_default_na=False))
    pb = sp.hstack(matrices, format='csr')
    sm = pd.concat(pbmeta, ignore_index=True)
    if sm.pb_id.duplicated().any():
        raise ValueError('Duplicate pseudobulk IDs')
    write_mtx_new(pb, root / 'pseudobulk_counts.mtx.gz')
    write_csv_new(genes, root / 'pseudobulk_genes.tsv', sep='\t')
    write_csv_new(sm, root / 'pseudobulk_samples.tsv', sep='\t')
    for name in ('cellmetadata_allinput_v1.csv.gz', 'all_qc_cellmetadata_v1.csv.gz', 'fixed4_allcells_expression_v1.csv.gz', 'sample_celltype_gene_summary_v1.csv', 'cell_composition_v1.csv', 'canonical_marker_means_v1.csv', 'canonical_markers_retainedcells_v1.csv.gz'):
        _concat_csv([d / name for d in dirs], root / name)
    write_csv_new(pd.DataFrame([{'sample': s['sample'], 'group': s['group'], 'input_cells': s['input_cells'], 'qc_pass': s['qc_pass'], 'qc_failed': s['qc_failed'], 'retained': s['retained'], 'doublet_status': s['doublets']['status'], 'doublet_threshold': s['doublets'].get('threshold'), 'predicted_doublets': s['doublets'].get('predicted_doublets', 0)} for s in summaries]), root / 'sample_QC_summary_v1.csv')
    result = {'dataset': dataset, 'source': source or {}, 'sample_names': [p.parent.name for p in done], 'n_samples': len(done), 'input_cells': sum((s['input_cells'] for s in summaries)), 'qc_pass': sum((s['qc_pass'] for s in summaries)), 'retained': sum((s['retained'] for s in summaries)), 'pseudobulk_shape_genes_by_samples': list(pb.shape), 'counts_sha256': sha256_file(root / 'pseudobulk_counts.mtx.gz'), 'note': 'Source integer counts summed only over QC-passed retained singlets; treatment group remains separate.'}
    write_json_new(root / 'FINAL_COMPLETE_v1.json', result)
    return result

def run_t2d(mainroot, model_path):
    mainroot = Path(mainroot)
    archive = mainroot / 'data/GSE255566/GSE255566_RAW.tar'
    download = json.loads((mainroot / 'metadata/GSE255566_DOWNLOAD.json').read_text())
    actual = sha256_file(archive)
    expected = download.get('sha256')
    if expected and actual != expected:
        raise ValueError('Archive source-lock hash mismatch')
    outdir = mainroot / 'results/processed_v2/GSE255566'
    outdir.mkdir(parents=True, exist_ok=True)
    lock = {'archive': str(archive.resolve()), 'archive_sha256': actual, 'archive_bytes': archive.stat().st_size, 'model_sha256': sha256_file(model_path), 'script_sha256': sha256_file(__file__), 'protocol_sha256': sha256_file(mainroot / 'metadata/ANALYSIS_PROTOCOL_FROZEN.json')}
    if not (outdir / 'SOURCE_LOCK_v1.json').exists():
        write_json_new(outdir / 'SOURCE_LOCK_v1.json', lock)
    elif json.loads((outdir / 'SOURCE_LOCK_v1.json').read_text()) != lock:
        raise ValueError('Source or implementation lock differs from existing run')
    with tarfile.open(archive, 'r:') as tar:
        members = {Path(m.name).name: m for m in tar.getmembers() if m.isfile()}
        for i, label in enumerate(('Con1', 'Con2', 'Con3', 'Mod1', 'Mod2', 'Mod3', 'treat2', 'treat3')):
            gsm = 'GSM' + str(8075321 + i)
            prefix = gsm + '_' + label
            group = 'Control' if label.startswith('Con') else 'T2D' if label.startswith('Mod') else 'Treatment'
            ss = {'archive_sha256': actual, 'archive_member_prefix': prefix, 'donor_id': label if group != 'Treatment' else 'unknown', 'sex': 'unknown', 'age': 'unknown', 'batch': 'unknown', 'untreated_status': 'not_verified' if group == 'T2D' else 'not_applicable', 'treatment_pairing': 'unknown' if group == 'Treatment' else 'not_applicable'}
            cp = outdir / 'samples' / label / 'COMPLETE_v1.json'
            if cp.exists():
                old = json.loads(cp.read_text())
                if old['source'] != ss:
                    raise ValueError('Checkpoint source identity mismatch')
                continue
            with gzip.GzipFile(fileobj=tar.extractfile(members[prefix + '_features.tsv.gz'])) as f:
                genes = pd.read_csv(f, sep='\t', header=None, dtype=str, keep_default_na=False).iloc[:, :2]
                genes.columns = ['gene_id', 'gene']
            with gzip.GzipFile(fileobj=tar.extractfile(members[prefix + '_barcodes.tsv.gz'])) as f:
                barcodes = pd.read_csv(f, sep='\t', header=None, dtype=str)[0].tolist()
            with gzip.GzipFile(fileobj=tar.extractfile(members[prefix + '_matrix.mtx.gz'])) as f:
                counts = scipy.io.mmread(f).T.tocsr()
            process_sample(counts, genes, barcodes, sample=label, group=group, outdir=outdir, model_path=model_path, dataset='GSE255566', source=ss, donor_status='nominal_participant_paper_support' if group != 'Treatment' else 'identity_and_pairing_unknown')
            del counts
            gc.collect()
    print(json.dumps(finalize(outdir, dataset='GSE255566', source=lock), indent=2), flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Frozen, sample-wise PBMC QC/annotation and additive pseudobulk exports.\n\nAll functions preserve source counts. Output directories must be new, or contain\nan already completed sample with the same source identity. No source extraction,\nnormalized whole-transcriptome object, or all-gene differential testing is saved.\n"""\nfrom __future__ import annotations\nimport argparse\nimport gc\nimport gzip\nimport hashlib\nimport importlib.metadata\nimport io\nimport json\nimport os\nfrom pathlib import Path\nimport re\nimport tarfile\nimport time\nimport numpy as np\nimport pandas as pd\nimport scipy.io\nimport scipy.sparse as sp\nFIXED_GENES = (\'FCRL6\', \'GZMA\', \'GZMH\', \'KLRD1\')\nCELL_TYPES = (\'NK\', \'CD8_T\', \'CD4_T\', \'T_other\', \'B\', \'Plasma\', \'Monocyte\', \'DC\', \'pDC\', \'Platelet\', \'Erythroid\', \'Other/Unassigned\')\nSEED = 20261101\nMARKERS = (\'CD3D\', \'CD3E\', \'TRAC\', \'CD8A\', \'CD8B\', \'CD4\', \'IL7R\', \'CCR7\', \'NCAM1\', \'NCR1\', \'FCGR3A\', \'TYROBP\', \'GNLY\', \'PRF1\', \'NKG7\', \'KLRF1\', \'MS4A1\', \'CD79A\', \'CD79B\', \'CD19\', \'JCHAIN\', \'MZB1\', \'SDC1\', \'CD14\', \'LYZ\', \'S100A8\', \'S100A9\', \'FCER1A\', \'CLEC10A\', \'CLEC9A\', \'GZMB\', \'TCF4\', \'IL3RA\', \'PPBP\', \'PF4\', \'HBA1\', \'HBB\')\nQC = {\'min_genes\': 200, \'max_genes\': 6000, \'min_UMI\': 500, \'max_UMI\': 40000, \'max_mito_percent\': 20, \'doublet_expected_rate\': 0.06, \'seed\': SEED}\n\ndef sha256_file(path):\n    h = hashlib.sha256()\n    with open(path, \'rb\') as f:\n        for b in iter(lambda: f.read(1024 * 1024), b\'\'):\n            h.update(b)\n    return h.hexdigest()\n\ndef write_json_new(path, obj):\n    with Path(path).open(\'x\') as f:\n        json.dump(obj, f, indent=2, default=str)\n        f.write(\'\\n\')\n\ndef write_csv_new(df, path, sep=\',\'):\n    path = Path(path)\n    if path.suffix == \'.gz\':\n        with path.open(\'xb\') as raw:\n            with gzip.GzipFile(fileobj=raw, mode=\'wb\', mtime=0) as gz:\n                with io.TextIOWrapper(gz) as f:\n                    df.to_csv(f, index=False, sep=sep)\n    else:\n        with path.open(\'x\') as f:\n            df.to_csv(f, index=False, sep=sep)\n\ndef write_mtx_new(matrix, path):\n    with Path(path).open(\'xb\') as raw:\n        with gzip.GzipFile(fileobj=raw, mode=\'wb\', mtime=0) as gz:\n            scipy.io.mmwrite(gz, matrix, field=\'integer\', symmetry=\'general\')\n\ndef broad_type(label):\n    """Reference-label mapping fixed before seeing disease effects."""\n    t = str(label).lower()\n    if t in (\'heterogeneous\', \'unassigned\', \'nan\', \'\'):\n        return \'Other/Unassigned\'\n    if \'mast\' in t:\n        return \'Other/Unassigned\'\n    if \'plasma\' in t and \'plasmacytoid\' not in t:\n        return \'Plasma\'\n    if \'pdc\' in t or \'plasmacytoid\' in t:\n        return \'pDC\'\n    if \'erythro\' in t and \'macrophage\' not in t:\n        return \'Erythroid\'\n    if \'platelet\' in t or \'megakaryo\' in t or t == \'early mk\':\n        return \'Platelet\'\n    if \'nk\' in t and \'nkt\' not in t:\n        return \'NK\'\n    if any((v in t for v in (\'monocyt\', \'macrophage\', \'mono-mac\', \'kupffer\', \'hofbauer\'))):\n        return \'Monocyte\'\n    if re.search(\'(^|\\\\W)dc[123]?(\\\\W|$)\', t) or \'dendritic\' in t:\n        return \'DC\'\n    if \'b cell\' in t or \'-b cell\' in t:\n        return \'B\'\n    if \'cd4\' in t or \'helper t\' in t or \'regulatory t\' in t or (\'treg\' in t):\n        return \'CD4_T\'\n    if \'cd8\' in t or \'cytotoxic t\' in t:\n        return \'CD8_T\'\n    if re.search(\'(^|[ -])t cells?($|\\\\W)\', t) or any((v in t for v in (\'nkt\', \'mait\', \'thymocyte\', \'t(agonist)\', \'t lymphoid\'))):\n        return \'T_other\'\n    return \'Other/Unassigned\'\n\ndef _dense_gene_values(counts, symbols, targets):\n    """Sum duplicate symbols only for display; retain original features for PB."""\n    return np.column_stack([np.asarray(counts[:, np.flatnonzero(symbols == g)].sum(axis=1)).ravel() if np.any(symbols == g) else np.full(counts.shape[0], np.nan) for g in targets])\n\ndef _concat_csv(inputs, output):\n    """Bound memory by copying CSV records, one file at a time."""\n    output = Path(output)\n    with output.open(\'xb\') as raw:\n        if output.suffix == \'.gz\':\n            holder = gzip.GzipFile(fileobj=raw, mode=\'wb\', mtime=0)\n        else:\n            holder = raw\n        try:\n            for i, p in enumerate(inputs):\n                opener = gzip.open if str(p).endswith(\'.gz\') else open\n                with opener(p, \'rb\') as src:\n                    header = src.readline()\n                    if i == 0:\n                        holder.write(header)\n                    for b in iter(lambda: src.read(1024 * 1024), b\'\'):\n                        holder.write(b)\n        finally:\n            if holder is not raw:\n                holder.close()\n\ndef process_sample(counts, genes, barcodes, *, sample, group, outdir, model_path, dataset, source=None, donor_status=\'nominal_sample\'):\n    """Process one cells×features CSR matrix; return completed JSON summary.\n\n    genes: DataFrame with gene_id,gene in count-column order.\n    source: optional donor_id,sex,age,batch plus immutable source provenance.\n    Analyses retain QC-passed predicted singlets; failed automatic Scrublet\n    estimation retains QC-passed cells with the failure explicitly recorded.\n    """\n    started = time.time()\n    source = dict(source or {})\n    sample = str(sample)\n    if re.search(\'[^A-Za-z0-9_.-]\', sample):\n        raise ValueError(\'sample must be a safe single path component\')\n    root = Path(outdir)\n    sdir = root / \'samples\' / sample\n    complete = sdir / \'COMPLETE_v1.json\'\n    if complete.exists():\n        old = json.loads(complete.read_text())\n        if old[\'source\'] != source or old[\'group\'] != group or old[\'dataset\'] != dataset:\n            raise ValueError(f\'Completed source identity mismatch: {sdir}\')\n        print(f\'{sample}: verified completed sample, reusing immutable checkpoint\', flush=True)\n        return old\n    sdir.mkdir(parents=True, exist_ok=False)\n    print(f\'{sample}: starting {counts.shape}\', flush=True)\n    genes = genes.loc[:, [\'gene_id\', \'gene\']].astype(str).reset_index(drop=True)\n    counts = sp.csr_matrix(counts)\n    counts.sum_duplicates()\n    counts.eliminate_zeros()\n    if counts.shape != (len(barcodes), len(genes)):\n        raise ValueError(\'Count dimensions do not match barcodes/features\')\n    if len(set(map(str, barcodes))) != len(barcodes):\n        raise ValueError(\'Duplicate source barcodes within a sample\')\n    valid = np.all(np.isfinite(counts.data)) and np.all(counts.data >= 0)\n    integer = np.all(counts.data == np.floor(counts.data))\n    if not (valid and integer):\n        raise ValueError(\'Source matrix is not nonnegative finite integer counts\')\n    counts = counts.astype(np.int64, copy=False)\n    symbols = genes.gene.to_numpy()\n    total = np.asarray(counts.sum(axis=1)).ravel()\n    n_genes = np.diff(counts.indptr)\n    is_mt = np.char.startswith(symbols.astype(str), \'MT-\')\n    mt = np.asarray(counts[:, is_mt].sum(axis=1)).ravel()\n    mt_pct = np.divide(mt * 100.0, total, out=np.zeros(len(total)), where=total > 0)\n    qc = (n_genes >= QC[\'min_genes\']) & (n_genes <= QC[\'max_genes\']) & (total >= QC[\'min_UMI\']) & (total <= QC[\'max_UMI\']) & (mt_pct <= QC[\'max_mito_percent\'])\n    meta = pd.DataFrame({\'cell_id\': [sample + \'__\' + str(b) for b in barcodes], \'barcode\': list(map(str, barcodes)), \'dataset\': dataset, \'sample\': sample, \'group\': group, \'donor_id\': str(source.get(\'donor_id\', sample)), \'donor_status\': donor_status, \'sex\': str(source.get(\'sex\', \'unknown\')), \'age\': str(source.get(\'age\', \'unknown\')), \'batch\': str(source.get(\'batch\', \'unknown\')), \'n_genes\': n_genes, \'total_counts\': total, \'pct_mt\': mt_pct, \'qc_pass\': qc, \'doublet_score\': np.nan, \'predicted_doublet\': False, \'retained\': False, \'cell_type\': \'Failed_QC\', \'raw_CellTypist_label\': \'not_annotated\', \'majority_CellTypist_label\': \'not_annotated\', \'pred_confidence\': np.nan, \'majority_label_probability\': np.nan, \'over_clustering\': \'not_annotated\'})\n    qix = np.flatnonzero(qc)\n    scrub_info = {\'expected_rate\': 0.06, \'seed\': SEED, \'status\': \'not_run_no_qc_cells\'}\n    doublets = np.zeros(len(qix), dtype=bool)\n    if len(qix):\n        import scanpy as sc\n        import anndata as ad\n        scrub = ad.AnnData(counts[qix, :].copy())\n        try:\n            sc.pp.scrublet(scrub, expected_doublet_rate=0.06, random_state=SEED, n_prin_comps=30, verbose=True)\n            scores = scrub.obs[\'doublet_score\'].to_numpy()\n            calls = scrub.obs.get(\'predicted_doublet\')\n            meta.loc[qix, \'doublet_score\'] = scores\n            threshold = scrub.uns.get(\'scrublet\', {}).get(\'threshold\')\n            if calls is not None and threshold is not None and np.isfinite(threshold):\n                doublets = np.asarray(calls, dtype=bool)\n                scrub_info.update(status=\'automatic_threshold\', threshold=float(threshold), predicted_doublets=int(doublets.sum()))\n            else:\n                scrub_info.update(status=\'automatic_threshold_failed_all_qc_retained\', threshold=None, predicted_doublets=0)\n        except Exception as e:\n            if \'doublet_score\' in scrub.obs:\n                meta.loc[qix, \'doublet_score\'] = scrub.obs[\'doublet_score\'].to_numpy()\n            scrub_info.update(status=\'not_estimable_all_qc_retained\', error=repr(e), threshold=None, predicted_doublets=0)\n        scrub_info[\'implementation\'] = \'scanpy.pp.scrublet\'\n        del scrub\n        gc.collect()\n    meta.loc[qix, \'predicted_doublet\'] = doublets\n    meta.loc[qix, \'cell_type\'] = \'Excluded_doublet\'\n    keep = qix[~doublets]\n    meta.loc[keep, \'retained\'] = True\n    model_hash = sha256_file(model_path)\n    model_info = {\'path\': str(Path(model_path).resolve()), \'sha256\': model_hash, \'fixed_genes_excluded_from_annotation_and_overclustering\': list(FIXED_GENES)}\n    if len(keep):\n        import scanpy as sc\n        import anndata as ad\n        import celltypist\n        np.random.seed(SEED)\n        model = celltypist.models.Model.load(str(model_path))\n        excluded = np.isin(np.asarray(model.classifier.features), FIXED_GENES)\n        model.classifier.coef_[:, excluded] = 0\n        model_info[\'zeroed_model_features\'] = np.asarray(model.classifier.features)[excluded].tolist()\n        use = ~np.isin(symbols, FIXED_GENES)\n        ann_counts = counts[keep, :][:, use].astype(np.float32)\n        ann_symbols = symbols[use]\n        codes, unique_symbols = pd.factorize(ann_symbols, sort=False)\n        if len(unique_symbols) != len(ann_symbols):\n            mapping = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(unique_symbols)))\n            ann_counts = (ann_counts @ mapping).astype(np.float32).tocsr()\n        annotation = ad.AnnData(ann_counts, obs=pd.DataFrame(index=meta.loc[keep, \'cell_id\'].values), var=pd.DataFrame(index=pd.Index(unique_symbols)))\n        sc.pp.normalize_total(annotation, target_sum=10000.0)\n        sc.pp.log1p(annotation)\n        pred = celltypist.annotate(annotation, model=model, majority_voting=True)\n        labels = pred.predicted_labels\n        rawlabel = labels[\'predicted_labels\'].astype(str).to_numpy()\n        majority = labels.get(\'majority_voting\', labels[\'predicted_labels\']).astype(str).to_numpy()\n        probability = pred.probability_matrix\n        conf = probability.max(axis=1).to_numpy()\n        assigned_prob = np.array([probability.iloc[i].get(t, np.nan) for i, t in enumerate(majority)])\n        meta.loc[keep, \'raw_CellTypist_label\'] = rawlabel\n        meta.loc[keep, \'majority_CellTypist_label\'] = majority\n        meta.loc[keep, \'pred_confidence\'] = conf\n        meta.loc[keep, \'majority_label_probability\'] = assigned_prob\n        meta.loc[keep, \'cell_type\'] = [broad_type(t) for t in majority]\n        if \'over_clustering\' in labels:\n            meta.loc[keep, \'over_clustering\'] = labels[\'over_clustering\'].astype(str).to_numpy()\n        model_info[\'reference_labels_mapping\'] = {str(t): broad_type(t) for t in model.cell_types}\n        del pred, annotation, model, ann_counts\n        gc.collect()\n    fixed_raw = _dense_gene_values(counts, symbols, FIXED_GENES)\n    fixed_log = np.log1p(np.divide(fixed_raw * 10000.0, total[:, None], out=np.zeros_like(fixed_raw, dtype=float), where=total[:, None] > 0))\n    fx = meta[[\'cell_id\', \'dataset\', \'sample\', \'group\', \'donor_id\', \'cell_type\', \'qc_pass\', \'predicted_doublet\', \'retained\']].copy()\n    for j, g in enumerate(FIXED_GENES):\n        fx[g + \'_raw\'] = fixed_raw[:, j]\n        fx[g + \'_logCP10k\'] = fixed_log[:, j]\n    write_csv_new(meta, sdir / \'cellmetadata_allinput_v1.csv.gz\')\n    write_csv_new(meta[qc], sdir / \'all_qc_cellmetadata_v1.csv.gz\')\n    write_csv_new(fx, sdir / \'fixed4_allcells_expression_v1.csv.gz\')\n    canonical_raw = _dense_gene_values(counts[keep, :], symbols, MARKERS) if len(keep) else np.empty((0, len(MARKERS)))\n    canonical_log = np.log1p(np.divide(canonical_raw * 10000.0, total[keep, None], out=np.zeros_like(canonical_raw, dtype=float), where=total[keep, None] > 0))\n    canonical_cells = meta.loc[keep, [\'cell_id\', \'dataset\', \'sample\', \'group\', \'donor_id\', \'cell_type\', \'raw_CellTypist_label\', \'majority_CellTypist_label\', \'pred_confidence\']].copy()\n    for j, g in enumerate(MARKERS):\n        canonical_cells[g + \'_raw\'] = canonical_raw[:, j]\n        canonical_cells[g + \'_logCP10k\'] = canonical_log[:, j]\n    write_csv_new(canonical_cells, sdir / \'canonical_markers_retainedcells_v1.csv.gz\')\n    ctype = meta.loc[keep, \'cell_type\'].to_numpy()\n    types = list(CELL_TYPES)\n    base = {\'dataset\': dataset, \'sample\': sample, \'group\': group, \'donor_id\': str(source.get(\'donor_id\', sample)), \'sex\': str(source.get(\'sex\', \'unknown\')), \'age\': str(source.get(\'age\', \'unknown\')), \'batch\': str(source.get(\'batch\', \'unknown\'))}\n    pb_rows = []\n    pb_meta = []\n    summaries = []\n    compositions = []\n    marker_rows = []\n    for t in types:\n        mask = ctype == t\n        ix = keep[mask]\n        compositions.append(dict(base, cell_type=t, n_cells=int(mask.sum()), n_sample_retained=len(keep), proportion=float(mask.mean())))\n        if not mask.any():\n            for g in FIXED_GENES:\n                summaries.append(dict(base, cell_type=t, gene=g, n_cells=0, measured=bool(np.any(symbols == g)), detected_cells=0, percent_detected=np.nan, mean_logCP10k=np.nan, mean_CP10k=np.nan, pseudobulk_count=0))\n            continue\n        pb_id = sample + \'__\' + t.replace(\'/\', \'_\')\n        pb_rows.append(sp.csr_matrix(counts[ix, :].sum(axis=0)))\n        pb_meta.append(dict(base, pb_id=pb_id, cell_type=t, n_cells=int(mask.sum())))\n        for j, g in enumerate(FIXED_GENES):\n            summaries.append(dict(base, cell_type=t, gene=g, n_cells=int(mask.sum()), measured=bool(np.any(symbols == g)), detected_cells=int(np.sum(fixed_raw[ix, j] > 0)), percent_detected=float(np.mean(fixed_raw[ix, j] > 0) * 100), mean_logCP10k=float(np.mean(fixed_log[ix, j])), mean_CP10k=float(np.mean(np.expm1(fixed_log[ix, j]))), pseudobulk_count=int(np.nansum(fixed_raw[ix, j]))))\n        for j, g in enumerate(MARKERS):\n            marker_rows.append(dict(base, cell_type=t, gene=g, n_cells=int(mask.sum()), measured=bool(np.any(symbols == g)), percent_detected=float(np.mean(canonical_raw[mask, j] > 0) * 100), mean_logCP10k=float(np.mean(canonical_log[mask, j]))))\n    if not pb_rows:\n        raise ValueError(f\'No retained cells in {sample}\')\n    pb = sp.vstack(pb_rows, format=\'csr\').T.tocsr()\n    write_mtx_new(pb, sdir / \'raw_pseudobulk_counts_v1.mtx.gz\')\n    write_csv_new(genes, sdir / \'raw_pseudobulk_genes_v1.tsv\', sep=\'\\t\')\n    write_csv_new(pd.DataFrame(pb_meta), sdir / \'raw_pseudobulk_samples_v1.tsv\', sep=\'\\t\')\n    write_csv_new(pd.DataFrame(summaries), sdir / \'sample_celltype_gene_summary_v1.csv\')\n    write_csv_new(pd.DataFrame(compositions), sdir / \'cell_composition_v1.csv\')\n    write_csv_new(pd.DataFrame(marker_rows), sdir / \'canonical_marker_means_v1.csv\')\n    summary = {\'sample\': sample, \'group\': group, \'dataset\': dataset, \'source\': source, \'donor_status\': donor_status, \'input_cells\': counts.shape[0], \'input_genes\': counts.shape[1], \'input_nnz\': int(counts.nnz), \'raw_values_all_validated_nonnegative_integer\': True, \'qc_pass\': int(qc.sum()), \'qc_failed\': int((~qc).sum()), \'retained\': len(keep), \'doublets\': scrub_info, \'qc\': QC, \'annotation\': model_info, \'cell_type_counts\': meta.loc[keep, \'cell_type\'].value_counts().to_dict(), \'fixed_gene_feature_rows\': {g: np.flatnonzero(symbols == g).tolist() for g in FIXED_GENES}, \'elapsed_seconds\': time.time() - started, \'software\': {x: importlib.metadata.version(x) for x in (\'numpy\', \'scipy\', \'pandas\', \'scanpy\', \'anndata\', \'scrublet\', \'celltypist\', \'leidenalg\')}}\n    write_json_new(complete, summary)\n    print(f"{sample}: COMPLETE retained {len(keep)}/{len(qix)} QC cells, {summary[\'cell_type_counts\']}", flush=True)\n    return summary\n\ndef finalize(outdir, *, dataset, source=None):\n    """Merge only completed checkpoints, preserving one raw feature ordering."""\n    root = Path(outdir)\n    done = sorted(root.glob(\'samples/*/COMPLETE_v1.json\'))\n    if not done:\n        raise ValueError(\'No completed sample checkpoints\')\n    if (root / \'FINAL_COMPLETE_v1.json\').exists():\n        old = json.loads((root / \'FINAL_COMPLETE_v1.json\').read_text())\n        if old[\'sample_names\'] != [p.parent.name for p in done]:\n            raise ValueError(\'Final checkpoint sample set mismatch\')\n        return old\n    summaries = [json.loads(p.read_text()) for p in done]\n    if any((s[\'dataset\'] != dataset for s in summaries)):\n        raise ValueError(\'Mixed datasets\')\n    dirs = [p.parent for p in done]\n    genes = pd.read_csv(dirs[0] / \'raw_pseudobulk_genes_v1.tsv\', sep=\'\\t\', dtype=str, keep_default_na=False)\n    matrices = []\n    pbmeta = []\n    for d in dirs:\n        g = pd.read_csv(d / \'raw_pseudobulk_genes_v1.tsv\', sep=\'\\t\', dtype=str, keep_default_na=False)\n        if not genes.equals(g):\n            raise ValueError(f\'Feature order mismatch: {d}\')\n        matrices.append(scipy.io.mmread(d / \'raw_pseudobulk_counts_v1.mtx.gz\').tocsr())\n        pbmeta.append(pd.read_csv(d / \'raw_pseudobulk_samples_v1.tsv\', sep=\'\\t\', keep_default_na=False))\n    pb = sp.hstack(matrices, format=\'csr\')\n    sm = pd.concat(pbmeta, ignore_index=True)\n    if sm.pb_id.duplicated().any():\n        raise ValueError(\'Duplicate pseudobulk IDs\')\n    write_mtx_new(pb, root / \'pseudobulk_counts.mtx.gz\')\n    write_csv_new(genes, root / \'pseudobulk_genes.tsv\', sep=\'\\t\')\n    write_csv_new(sm, root / \'pseudobulk_samples.tsv\', sep=\'\\t\')\n    for name in (\'cellmetadata_allinput_v1.csv.gz\', \'all_qc_cellmetadata_v1.csv.gz\', \'fixed4_allcells_expression_v1.csv.gz\', \'sample_celltype_gene_summary_v1.csv\', \'cell_composition_v1.csv\', \'canonical_marker_means_v1.csv\', \'canonical_markers_retainedcells_v1.csv.gz\'):\n        _concat_csv([d / name for d in dirs], root / name)\n    write_csv_new(pd.DataFrame([{\'sample\': s[\'sample\'], \'group\': s[\'group\'], \'input_cells\': s[\'input_cells\'], \'qc_pass\': s[\'qc_pass\'], \'qc_failed\': s[\'qc_failed\'], \'retained\': s[\'retained\'], \'doublet_status\': s[\'doublets\'][\'status\'], \'doublet_threshold\': s[\'doublets\'].get(\'threshold\'), \'predicted_doublets\': s[\'doublets\'].get(\'predicted_doublets\', 0)} for s in summaries]), root / \'sample_QC_summary_v1.csv\')\n    result = {\'dataset\': dataset, \'source\': source or {}, \'sample_names\': [p.parent.name for p in done], \'n_samples\': len(done), \'input_cells\': sum((s[\'input_cells\'] for s in summaries)), \'qc_pass\': sum((s[\'qc_pass\'] for s in summaries)), \'retained\': sum((s[\'retained\'] for s in summaries)), \'pseudobulk_shape_genes_by_samples\': list(pb.shape), \'counts_sha256\': sha256_file(root / \'pseudobulk_counts.mtx.gz\'), \'note\': \'Source integer counts summed only over QC-passed retained singlets; treatment group remains separate.\'}\n    write_json_new(root / \'FINAL_COMPLETE_v1.json\', result)\n    return result\n\ndef run_t2d(mainroot, model_path):\n    mainroot = Path(mainroot)\n    archive = mainroot / \'data/GSE255566/GSE255566_RAW.tar\'\n    download = json.loads((mainroot / \'metadata/GSE255566_DOWNLOAD.json\').read_text())\n    actual = sha256_file(archive)\n    expected = download.get(\'sha256\')\n    if expected and actual != expected:\n        raise ValueError(\'Archive source-lock hash mismatch\')\n    outdir = mainroot / \'results/processed_v2/GSE255566\'\n    outdir.mkdir(parents=True, exist_ok=True)\n    lock = {\'archive\': str(archive.resolve()), \'archive_sha256\': actual, \'archive_bytes\': archive.stat().st_size, \'model_sha256\': sha256_file(model_path), \'script_sha256\': sha256_file(__file__), \'protocol_sha256\': sha256_file(mainroot / \'metadata/ANALYSIS_PROTOCOL_FROZEN.json\')}\n    if not (outdir / \'SOURCE_LOCK_v1.json\').exists():\n        write_json_new(outdir / \'SOURCE_LOCK_v1.json\', lock)\n    elif json.loads((outdir / \'SOURCE_LOCK_v1.json\').read_text()) != lock:\n        raise ValueError(\'Source or implementation lock differs from existing run\')\n    with tarfile.open(archive, \'r:\') as tar:\n        members = {Path(m.name).name: m for m in tar.getmembers() if m.isfile()}\n        for i, label in enumerate((\'Con1\', \'Con2\', \'Con3\', \'Mod1\', \'Mod2\', \'Mod3\', \'treat2\', \'treat3\')):\n            gsm = \'GSM\' + str(8075321 + i)\n            prefix = gsm + \'_\' + label\n            group = \'Control\' if label.startswith(\'Con\') else \'T2D\' if label.startswith(\'Mod\') else \'Treatment\'\n            ss = {\'archive_sha256\': actual, \'archive_member_prefix\': prefix, \'donor_id\': label if group != \'Treatment\' else \'unknown\', \'sex\': \'unknown\', \'age\': \'unknown\', \'batch\': \'unknown\', \'untreated_status\': \'not_verified\' if group == \'T2D\' else \'not_applicable\', \'treatment_pairing\': \'unknown\' if group == \'Treatment\' else \'not_applicable\'}\n            cp = outdir / \'samples\' / label / \'COMPLETE_v1.json\'\n            if cp.exists():\n                old = json.loads(cp.read_text())\n                if old[\'source\'] != ss:\n                    raise ValueError(\'Checkpoint source identity mismatch\')\n                continue\n            with gzip.GzipFile(fileobj=tar.extractfile(members[prefix + \'_features.tsv.gz\'])) as f:\n                genes = pd.read_csv(f, sep=\'\\t\', header=None, dtype=str, keep_default_na=False).iloc[:, :2]\n                genes.columns = [\'gene_id\', \'gene\']\n            with gzip.GzipFile(fileobj=tar.extractfile(members[prefix + \'_barcodes.tsv.gz\'])) as f:\n                barcodes = pd.read_csv(f, sep=\'\\t\', header=None, dtype=str)[0].tolist()\n            with gzip.GzipFile(fileobj=tar.extractfile(members[prefix + \'_matrix.mtx.gz\'])) as f:\n                counts = scipy.io.mmread(f).T.tocsr()\n            process_sample(counts, genes, barcodes, sample=label, group=group, outdir=outdir, model_path=model_path, dataset=\'GSE255566\', source=ss, donor_status=\'nominal_participant_paper_support\' if group != \'Treatment\' else \'identity_and_pairing_unknown\')\n            del counts\n            gc.collect()\n    print(json.dumps(finalize(outdir, dataset=\'GSE255566\', source=lock), indent=2), flush=True)\nif __name__ == \'__main__\':\n    ap = argparse.ArgumentParser()\n    ap.add_argument(\'--mainroot\', required=True)\n    ap.add_argument(\'--model\', required=True)\n    args = ap.parse_args()\n    run_t2d(args.mainroot, args.model)', '<private-stage>', 'exec'), context)
