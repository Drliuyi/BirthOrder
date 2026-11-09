"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gc
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import numpy as np
import pandas as pd
from scipy.io import mmread
FIXED = ('FCRL6', 'GZMA', 'GZMH', 'KLRD1')

def emit(event: str, **items) -> None:
    print(json.dumps({'time': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'event': event, **items}, ensure_ascii=False), flush=True)

def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while (block := handle.read(1024 * 1024)):
            digest.update(block)
    return digest.hexdigest()

def download_ram(url: str, *, expected_size: int | None=None, maximum_bytes: int=70000000) -> tuple[bytes, dict]:
    """Fetch a single bounded public source, checking length and hashing bytes."""
    if not url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/'):
        raise ValueError(f'Source URL outside the frozen public NCBI scope: {url}')
    failure = None
    for attempt in range(1, 6):
        try:
            if attempt <= 5:
                source_basename = url.rsplit('/', 1)[1]
                source_sample = source_basename.split('_', 1)[0]
                geo_download = 'https://www.ncbi.nlm.nih.gov/geo/download/?' + urlencode({'acc': source_sample, 'format': 'file', 'file': source_basename})
                candidates = [geo_download, url.replace('https://ftp.ncbi.nlm.nih.gov/', 'https://ftp.ncbi.nih.gov/'), url.replace('https://', 'ftp://', 1), url + '?download=1', url]
                fetch_url = candidates[attempt - 1]
                fetched = subprocess.run(['curl', '--http2', '--location', '--fail', '--silent', '--show-error', '--connect-timeout', '20', '--max-time', '600', '--max-filesize', str(maximum_bytes), fetch_url], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
                if fetched.returncode:
                    raise OSError(f'curl HTTP/2 failed ({fetched.returncode}): ' + fetched.stderr.decode(errors='replace')[:500])
                payload = fetched.stdout
                actual = len(payload)
                if not actual or actual > maximum_bytes:
                    raise OSError(f'curl source payload length invalid: {actual}')
                if expected_size is not None and actual != expected_size:
                    raise OSError(f'Frozen source length differs: {actual} != {expected_size}')
                return (payload, {'url': url, 'bytes': actual, 'sha256': hashlib.sha256(payload).hexdigest(), 'download_method': 'curl_stdout_RAM', 'actual_fetch_url': fetch_url})
            digest = hashlib.sha256()
            with urlopen(Request(url, headers={'User-Agent': 'birthorder-scrna-audit/1.0'}), timeout=90) as response:
                declared = response.headers.get('Content-Length')
                if declared and int(declared) > maximum_bytes:
                    raise ValueError(f'Source exceeds RAM byte ceiling: {declared}')
                buffer = io.BytesIO()
                while (block := response.read(1024 * 1024)):
                    buffer.write(block)
                    digest.update(block)
                    if buffer.tell() > maximum_bytes:
                        raise ValueError('Source exceeds RAM byte ceiling during read')
                actual = buffer.tell()
                if declared is not None and actual != int(declared):
                    raise OSError(f'Incomplete HTTP response: {actual} != {declared}')
                if expected_size is not None and actual != expected_size:
                    raise OSError(f'Frozen source length differs: {actual} != {expected_size}')
                payload = buffer.getvalue()
            return (payload, {'url': url, 'bytes': actual, 'sha256': digest.hexdigest(), 'download_method': 'urllib_RAM'})
        except Exception as error:
            failure = error
            emit('source_retry', url=url, attempt=attempt, error=repr(error))
            if attempt < 5:
                time.sleep(min(2 ** attempt, 15))
    raise RuntimeError(f'Five downloads failed: {url}') from failure

def write_new_or_verify(path: Path, value: dict) -> None:
    """An immutable checkpoint never replaces a prior receipt."""
    if path.exists():
        previous = json.loads(path.read_text())
        if previous != value:
            raise RuntimeError(f'Existing receipt differs; preserve and investigate: {path}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write('\n')

def load_api(path: Path):
    spec = importlib.util.spec_from_file_location('scrna_samplewise', path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Cannot import shared samplewise module: {path}')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

def download_raw_donor(row: dict) -> dict:
    """Fetch three independent files for exactly one donor concurrently."""
    with ThreadPoolExecutor(max_workers=3) as transport:
        futures = {'features': transport.submit(download_ram, row['feature_url'], maximum_bytes=2000000), 'barcodes': transport.submit(download_ram, row['barcode_url'], maximum_bytes=2000000), 'matrix': transport.submit(download_ram, row['matrix_url'], expected_size=int(row['matrix_bytes']))}
        return {kind: future.result() for kind, future in futures.items()}

def load_donor(row: dict, raw: dict) -> tuple[object, pd.DataFrame, list[str], dict]:
    features_payload, features_source = raw['features']
    features = pd.read_csv(io.BytesIO(features_payload), sep='\t', header=None, compression='gzip', dtype=str, keep_default_na=False)
    if features.shape[1] < 2:
        raise ValueError('Feature file must include gene identifiers and symbols')
    genes = features.iloc[:, :2].copy()
    genes.columns = ['gene_id', 'gene']
    if features.shape[1] >= 3 and (not (features.iloc[:, 2] == 'Gene Expression').all()):
        raise ValueError('Unexpected non-Gene-Expression features; do not silently drop rows')
    missing = sorted(set(FIXED) - set(genes['gene']))
    if missing:
        raise ValueError(f'Fixed quartet coverage gate failed before disease analysis: {missing}')
    barcodes_payload, barcode_source = raw['barcodes']
    barcodes = gzip.decompress(barcodes_payload).decode().splitlines()
    if len(set(barcodes)) != len(barcodes):
        raise ValueError('Duplicated barcodes within donor')
    matrix_payload, matrix_source = raw['matrix']
    with gzip.GzipFile(fileobj=io.BytesIO(matrix_payload), mode='rb') as stream:
        prefix = stream.readline().decode().strip()
        if prefix != '%%MatrixMarket matrix coordinate integer general':
            raise ValueError(f'Original integer MatrixMarket header required: {prefix}')
    with gzip.GzipFile(fileobj=io.BytesIO(matrix_payload), mode='rb') as stream:
        counts = mmread(stream).T.tocsr()
    if counts.shape != (len(barcodes), len(genes)):
        raise ValueError(f'Count dimensions disagree with feature/barcode files: {counts.shape}')
    if counts.data.size and (not np.isfinite(counts.data).all() or (counts.data < 0).any() or (not np.equal(counts.data, np.floor(counts.data)).all())):
        raise ValueError('Counts must be finite nonnegative integers')
    receipt = {'dataset': DATASET, 'sample_id': row['sample_id'], 'donor_id': row['donor_id'], 'group': row['group'], 'sex': row['sex'], 'files': {'features': features_source, 'barcodes': barcode_source, 'matrix': matrix_source}, 'matrix_header': prefix, 'n_input_cells': counts.shape[0], 'n_input_genes': counts.shape[1], 'n_nonzero': int(counts.nnz), 'fixed_genes_present': list(FIXED), 'full_counts_retained_on_disk': False}
    del matrix_payload, features_payload, barcodes_payload, features
    gc.collect()
    return (counts, genes, barcodes, receipt)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(os.environ.get('SCRNA_TASK_ROOT', Path(__file__).resolve().parents[1])))
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--helper', type=Path)
    parser.add_argument('--model', type=Path, default=Path(os.environ['CELLTYPIST_MODEL']) if os.environ.get('CELLTYPIST_MODEL') else None)
    args = parser.parse_args()
    root = args.root.resolve()
    input_root = Path(os.environ.get('CAD_INPUT_ROOT', root / 'cad_audit'))
    script_root = Path(os.environ.get('CAD_SCRIPT_ROOT', root / 'scripts'))
    output_root = Path(os.environ.get('CAD_OUTPUT_ROOT', root / 'results/processed_v2'))
    manifest = (args.manifest or input_root / 'GSE289422_analysis_manifest_v2.tsv').resolve()
    output = (args.output or output_root / DATASET).resolve()
    helper = (args.helper or script_root / '01_scrna_samplewise.py').resolve()
    model = args.model.resolve() if args.model else root / 'runtime/celltypist/data/models/Immune_All_Low.pkl'
    for required in (manifest, helper, model):
        if not required.is_file():
            raise FileNotFoundError(f'Required input path missing: {required}')
    if output == manifest.parent or output == helper.parent:
        raise ValueError('Output must be isolated from source metadata and shared helper directories')
    emit('resolved_roots', task_root=str(root), input_manifest=str(manifest), helper_root=str(helper.parent), output_root=str(output), model=str(model))
    rows = pd.read_csv(manifest, sep='\t', dtype=str, keep_default_na=False).to_dict('records')
    if len(rows) != 44 or len({r['sample_id'] for r in rows}) != 44 or len({r['donor_id'] for r in rows}) != 44:
        raise ValueError('Frozen design requires all 44 independent donors without duplication')
    groups = pd.Series([r['group'] for r in rows]).value_counts().to_dict()
    if groups != {'CAD': 29, 'Control': 15}:
        raise ValueError(f'Frozen disease/control design differs: {groups}')
    if {r['dataset'] for r in rows} != {DATASET}:
        raise ValueError('Manifest contains another dataset')
    output.mkdir(parents=True, exist_ok=True)
    design = {'dataset': DATASET, 'manifest_sha256': sha_file(manifest), 'helper_sha256': sha_file(helper), 'model_sha256': sha_file(model), 'expected_donors': 44, 'group_counts': groups, 'primary_disease_model': 'group plus public sex if full rank', 'missing_covariates': ['age', 'diabetes', 'statin', 'survival_subgroup'], 'count_input_mode': 'one donor in RAM; no full count matrix on disk', 'fixed_genes': list(FIXED), 'all_donors_selected_before_results': True}
    write_new_or_verify(output / 'frozen_streaming_design_v1.json', design)
    api = load_api(helper)

    def process_ready(row, raw):
        counts, genes, barcodes, receipt = load_donor(row, raw)
        raw.clear()
        receipt['manifest_sha256'] = design['manifest_sha256']
        receipt_path = output / 'source_receipts' / f"{row['sample_id']}.json"
        write_new_or_verify(receipt_path, receipt)
        source = {**receipt, 'donor_id': row['donor_id'], 'sex': row['sex'], 'batch': '', 'age': '', 'diabetes': '', 'statin': ''}
        result = api.process_sample(counts, genes, barcodes, sample=row['sample_id'], group=row['group'], outdir=output, model_path=model, dataset=DATASET, source=source, donor_status='independent_public_patient_id')
        del counts, genes, barcodes
        gc.collect()
        emit('donor_complete', sample=row['sample_id'], donor=row['donor_id'], result={k: result[k] for k in ('input_cells', 'qc_pass', 'retained', 'doublets', 'cell_type_counts')})
    pending = []
    for row in rows:
        checkpoint = output / 'samples' / row['sample_id'] / 'COMPLETE_v1.json'
        if not checkpoint.exists():
            pending.append(row)
            continue
        receipt_path = output / 'source_receipts' / f"{row['sample_id']}.json"
        receipt = json.loads(receipt_path.read_text())
        if receipt['manifest_sha256'] != design['manifest_sha256'] or receipt['sample_id'] != row['sample_id'] or receipt['donor_id'] != row['donor_id'] or (receipt['group'] != row['group']) or (receipt['sex'] != row['sex']) or (receipt['files']['matrix']['url'] != row['matrix_url']) or (receipt['files']['matrix']['bytes'] != int(row['matrix_bytes'])):
            raise RuntimeError(f'Checkpoint differs from frozen design: {checkpoint}')
        source = {**receipt, 'donor_id': row['donor_id'], 'sex': row['sex'], 'batch': '', 'age': '', 'diabetes': '', 'statin': ''}
        api.process_sample(None, None, None, sample=row['sample_id'], group=row['group'], outdir=output, model_path=model, dataset=DATASET, source=source, donor_status='independent_public_patient_id')
        emit('donor_reused', sample=row['sample_id'], donor=row['donor_id'])
    schedule = {'fixed_cohort_size': 44, 'maximum_concurrent_donor_downloads': 2, 'maximum_analysis_writers': 1, 'counts_saved_to_disk': False, 'processing_order': 'first completed transport; fixed cohort unchanged', 'failed_transport_policy': 'pending and retried after every fixed donor is attempted', 'scientific_QC_or_mapping_failure': 'stop and investigate'}
    write_new_or_verify(output / 'transport_schedule_v3.json', schedule)
    for transport_round in range(1, 4):
        if not pending:
            break
        todo = iter(pending)
        failed = []
        emit('transport_round', round=transport_round, pending=len(pending))
        with ThreadPoolExecutor(max_workers=2) as producer:
            active = {}

            def fill_queue():
                while len(active) < 2:
                    try:
                        row = next(todo)
                    except StopIteration:
                        break
                    active[producer.submit(download_raw_donor, row)] = row
                    emit('donor_download_start', sample=row['sample_id'], donor=row['donor_id'], round=transport_round)
            fill_queue()
            while active:
                future = next(as_completed(tuple(active)))
                row = active.pop(future)
                try:
                    raw = future.result()
                except Exception as error:
                    failed.append(row)
                    event = {'sample': row['sample_id'], 'donor': row['donor_id'], 'round': transport_round, 'error': repr(error), 'status': 'PENDING_transport_not_excluded'}
                    emit('donor_transport_pending', **event)
                    with (output / 'transport_pending_events_v3.jsonl').open('a') as handle:
                        handle.write(json.dumps(event) + '\n')
                    del future
                    gc.collect()
                    fill_queue()
                    continue
                del future
                process_ready(row, raw)
                del raw
                gc.collect()
                fill_queue()
        pending = failed
    if pending:
        raise RuntimeError('Full fixed cohort incomplete; no finalization or CAD fit permitted. Pending: ' + ', '.join((row['sample_id'] for row in pending)))
    completed = list(output.glob('samples/*/COMPLETE_v1.json'))
    if len(completed) != 44:
        raise RuntimeError(f'Final full-cohort gate failed: {len(completed)} != 44')
    summary = api.finalize(output, dataset=DATASET, source=design)
    emit('dataset_complete', summary=summary)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""GSE289422: process one donor at a time without staging count matrices.\n\nOriginal compressed count bytes exist only in RAM. Small source receipts,\nper-donor QC, fixed-gene expression and full-gene pseudobulk are retained.\nNo donor selection uses gene-expression effects or association P values.\n"""\nfrom __future__ import annotations\nimport argparse\nfrom concurrent.futures import ThreadPoolExecutor, as_completed\nimport gc\nimport gzip\nimport hashlib\nimport importlib.util\nimport io\nimport json\nimport os\nfrom pathlib import Path\nimport subprocess\nimport sys\nimport time\nfrom urllib.parse import urlencode\nfrom urllib.request import Request, urlopen\nfor variable in (\'OPENBLAS_NUM_THREADS\', \'OMP_NUM_THREADS\', \'MKL_NUM_THREADS\', \'NUMBA_NUM_THREADS\', \'VECLIB_MAXIMUM_THREADS\'):\n    os.environ.setdefault(variable, \'2\')\nimport numpy as np\nimport pandas as pd\nfrom scipy.io import mmread\nFIXED = (\'FCRL6\', \'GZMA\', \'GZMH\', \'KLRD1\')\nDATASET = \'GSE289422\'\n\ndef emit(event: str, **items) -> None:\n    print(json.dumps({\'time\': time.strftime(\'%Y-%m-%dT%H:%M:%S%z\'), \'event\': event, **items}, ensure_ascii=False), flush=True)\n\ndef sha_file(path: Path) -> str:\n    digest = hashlib.sha256()\n    with path.open(\'rb\') as handle:\n        while (block := handle.read(1024 * 1024)):\n            digest.update(block)\n    return digest.hexdigest()\n\ndef download_ram(url: str, *, expected_size: int | None=None, maximum_bytes: int=70000000) -> tuple[bytes, dict]:\n    """Fetch a single bounded public source, checking length and hashing bytes."""\n    if not url.startswith(\'https://ftp.ncbi.nlm.nih.gov/geo/\'):\n        raise ValueError(f\'Source URL outside the frozen public NCBI scope: {url}\')\n    failure = None\n    for attempt in range(1, 6):\n        try:\n            if attempt <= 5:\n                source_basename = url.rsplit(\'/\', 1)[1]\n                source_sample = source_basename.split(\'_\', 1)[0]\n                geo_download = \'https://www.ncbi.nlm.nih.gov/geo/download/?\' + urlencode({\'acc\': source_sample, \'format\': \'file\', \'file\': source_basename})\n                candidates = [geo_download, url.replace(\'https://ftp.ncbi.nlm.nih.gov/\', \'https://ftp.ncbi.nih.gov/\'), url.replace(\'https://\', \'ftp://\', 1), url + \'?download=1\', url]\n                fetch_url = candidates[attempt - 1]\n                fetched = subprocess.run([\'curl\', \'--http2\', \'--location\', \'--fail\', \'--silent\', \'--show-error\', \'--connect-timeout\', \'20\', \'--max-time\', \'600\', \'--max-filesize\', str(maximum_bytes), fetch_url], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)\n                if fetched.returncode:\n                    raise OSError(f\'curl HTTP/2 failed ({fetched.returncode}): \' + fetched.stderr.decode(errors=\'replace\')[:500])\n                payload = fetched.stdout\n                actual = len(payload)\n                if not actual or actual > maximum_bytes:\n                    raise OSError(f\'curl source payload length invalid: {actual}\')\n                if expected_size is not None and actual != expected_size:\n                    raise OSError(f\'Frozen source length differs: {actual} != {expected_size}\')\n                return (payload, {\'url\': url, \'bytes\': actual, \'sha256\': hashlib.sha256(payload).hexdigest(), \'download_method\': \'curl_stdout_RAM\', \'actual_fetch_url\': fetch_url})\n            digest = hashlib.sha256()\n            with urlopen(Request(url, headers={\'User-Agent\': \'birthorder-scrna-audit/1.0\'}), timeout=90) as response:\n                declared = response.headers.get(\'Content-Length\')\n                if declared and int(declared) > maximum_bytes:\n                    raise ValueError(f\'Source exceeds RAM byte ceiling: {declared}\')\n                buffer = io.BytesIO()\n                while (block := response.read(1024 * 1024)):\n                    buffer.write(block)\n                    digest.update(block)\n                    if buffer.tell() > maximum_bytes:\n                        raise ValueError(\'Source exceeds RAM byte ceiling during read\')\n                actual = buffer.tell()\n                if declared is not None and actual != int(declared):\n                    raise OSError(f\'Incomplete HTTP response: {actual} != {declared}\')\n                if expected_size is not None and actual != expected_size:\n                    raise OSError(f\'Frozen source length differs: {actual} != {expected_size}\')\n                payload = buffer.getvalue()\n            return (payload, {\'url\': url, \'bytes\': actual, \'sha256\': digest.hexdigest(), \'download_method\': \'urllib_RAM\'})\n        except Exception as error:\n            failure = error\n            emit(\'source_retry\', url=url, attempt=attempt, error=repr(error))\n            if attempt < 5:\n                time.sleep(min(2 ** attempt, 15))\n    raise RuntimeError(f\'Five downloads failed: {url}\') from failure\n\ndef write_new_or_verify(path: Path, value: dict) -> None:\n    """An immutable checkpoint never replaces a prior receipt."""\n    if path.exists():\n        previous = json.loads(path.read_text())\n        if previous != value:\n            raise RuntimeError(f\'Existing receipt differs; preserve and investigate: {path}\')\n        return\n    path.parent.mkdir(parents=True, exist_ok=True)\n    with path.open(\'x\') as handle:\n        json.dump(value, handle, indent=2, ensure_ascii=False)\n        handle.write(\'\\n\')\n\ndef load_api(path: Path):\n    spec = importlib.util.spec_from_file_location(\'scrna_samplewise\', path)\n    if spec is None or spec.loader is None:\n        raise RuntimeError(f\'Cannot import shared samplewise module: {path}\')\n    module = importlib.util.module_from_spec(spec)\n    sys.modules[spec.name] = module\n    spec.loader.exec_module(module)\n    return module\n\ndef download_raw_donor(row: dict) -> dict:\n    """Fetch three independent files for exactly one donor concurrently."""\n    with ThreadPoolExecutor(max_workers=3) as transport:\n        futures = {\'features\': transport.submit(download_ram, row[\'feature_url\'], maximum_bytes=2000000), \'barcodes\': transport.submit(download_ram, row[\'barcode_url\'], maximum_bytes=2000000), \'matrix\': transport.submit(download_ram, row[\'matrix_url\'], expected_size=int(row[\'matrix_bytes\']))}\n        return {kind: future.result() for kind, future in futures.items()}\n\ndef load_donor(row: dict, raw: dict) -> tuple[object, pd.DataFrame, list[str], dict]:\n    features_payload, features_source = raw[\'features\']\n    features = pd.read_csv(io.BytesIO(features_payload), sep=\'\\t\', header=None, compression=\'gzip\', dtype=str, keep_default_na=False)\n    if features.shape[1] < 2:\n        raise ValueError(\'Feature file must include gene identifiers and symbols\')\n    genes = features.iloc[:, :2].copy()\n    genes.columns = [\'gene_id\', \'gene\']\n    if features.shape[1] >= 3 and (not (features.iloc[:, 2] == \'Gene Expression\').all()):\n        raise ValueError(\'Unexpected non-Gene-Expression features; do not silently drop rows\')\n    missing = sorted(set(FIXED) - set(genes[\'gene\']))\n    if missing:\n        raise ValueError(f\'Fixed quartet coverage gate failed before disease analysis: {missing}\')\n    barcodes_payload, barcode_source = raw[\'barcodes\']\n    barcodes = gzip.decompress(barcodes_payload).decode().splitlines()\n    if len(set(barcodes)) != len(barcodes):\n        raise ValueError(\'Duplicated barcodes within donor\')\n    matrix_payload, matrix_source = raw[\'matrix\']\n    with gzip.GzipFile(fileobj=io.BytesIO(matrix_payload), mode=\'rb\') as stream:\n        prefix = stream.readline().decode().strip()\n        if prefix != \'%%MatrixMarket matrix coordinate integer general\':\n            raise ValueError(f\'Original integer MatrixMarket header required: {prefix}\')\n    with gzip.GzipFile(fileobj=io.BytesIO(matrix_payload), mode=\'rb\') as stream:\n        counts = mmread(stream).T.tocsr()\n    if counts.shape != (len(barcodes), len(genes)):\n        raise ValueError(f\'Count dimensions disagree with feature/barcode files: {counts.shape}\')\n    if counts.data.size and (not np.isfinite(counts.data).all() or (counts.data < 0).any() or (not np.equal(counts.data, np.floor(counts.data)).all())):\n        raise ValueError(\'Counts must be finite nonnegative integers\')\n    receipt = {\'dataset\': DATASET, \'sample_id\': row[\'sample_id\'], \'donor_id\': row[\'donor_id\'], \'group\': row[\'group\'], \'sex\': row[\'sex\'], \'files\': {\'features\': features_source, \'barcodes\': barcode_source, \'matrix\': matrix_source}, \'matrix_header\': prefix, \'n_input_cells\': counts.shape[0], \'n_input_genes\': counts.shape[1], \'n_nonzero\': int(counts.nnz), \'fixed_genes_present\': list(FIXED), \'full_counts_retained_on_disk\': False}\n    del matrix_payload, features_payload, barcodes_payload, features\n    gc.collect()\n    return (counts, genes, barcodes, receipt)\n\ndef main() -> None:\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--root\', type=Path, default=Path(os.environ.get(\'SCRNA_TASK_ROOT\', Path(__file__).resolve().parents[1])))\n    parser.add_argument(\'--manifest\', type=Path)\n    parser.add_argument(\'--output\', type=Path)\n    parser.add_argument(\'--helper\', type=Path)\n    parser.add_argument(\'--model\', type=Path, default=Path(os.environ[\'CELLTYPIST_MODEL\']) if os.environ.get(\'CELLTYPIST_MODEL\') else None)\n    args = parser.parse_args()\n    root = args.root.resolve()\n    input_root = Path(os.environ.get(\'CAD_INPUT_ROOT\', root / \'cad_audit\'))\n    script_root = Path(os.environ.get(\'CAD_SCRIPT_ROOT\', root / \'scripts\'))\n    output_root = Path(os.environ.get(\'CAD_OUTPUT_ROOT\', root / \'results/processed_v2\'))\n    manifest = (args.manifest or input_root / \'GSE289422_analysis_manifest_v2.tsv\').resolve()\n    output = (args.output or output_root / DATASET).resolve()\n    helper = (args.helper or script_root / \'01_scrna_samplewise.py\').resolve()\n    model = args.model.resolve() if args.model else root / \'runtime/celltypist/data/models/Immune_All_Low.pkl\'\n    for required in (manifest, helper, model):\n        if not required.is_file():\n            raise FileNotFoundError(f\'Required input path missing: {required}\')\n    if output == manifest.parent or output == helper.parent:\n        raise ValueError(\'Output must be isolated from source metadata and shared helper directories\')\n    emit(\'resolved_roots\', task_root=str(root), input_manifest=str(manifest), helper_root=str(helper.parent), output_root=str(output), model=str(model))\n    rows = pd.read_csv(manifest, sep=\'\\t\', dtype=str, keep_default_na=False).to_dict(\'records\')\n    if len(rows) != 44 or len({r[\'sample_id\'] for r in rows}) != 44 or len({r[\'donor_id\'] for r in rows}) != 44:\n        raise ValueError(\'Frozen design requires all 44 independent donors without duplication\')\n    groups = pd.Series([r[\'group\'] for r in rows]).value_counts().to_dict()\n    if groups != {\'CAD\': 29, \'Control\': 15}:\n        raise ValueError(f\'Frozen disease/control design differs: {groups}\')\n    if {r[\'dataset\'] for r in rows} != {DATASET}:\n        raise ValueError(\'Manifest contains another dataset\')\n    output.mkdir(parents=True, exist_ok=True)\n    design = {\'dataset\': DATASET, \'manifest_sha256\': sha_file(manifest), \'helper_sha256\': sha_file(helper), \'model_sha256\': sha_file(model), \'expected_donors\': 44, \'group_counts\': groups, \'primary_disease_model\': \'group plus public sex if full rank\', \'missing_covariates\': [\'age\', \'diabetes\', \'statin\', \'survival_subgroup\'], \'count_input_mode\': \'one donor in RAM; no full count matrix on disk\', \'fixed_genes\': list(FIXED), \'all_donors_selected_before_results\': True}\n    write_new_or_verify(output / \'frozen_streaming_design_v1.json\', design)\n    api = load_api(helper)\n\n    def process_ready(row, raw):\n        counts, genes, barcodes, receipt = load_donor(row, raw)\n        raw.clear()\n        receipt[\'manifest_sha256\'] = design[\'manifest_sha256\']\n        receipt_path = output / \'source_receipts\' / f"{row[\'sample_id\']}.json"\n        write_new_or_verify(receipt_path, receipt)\n        source = {**receipt, \'donor_id\': row[\'donor_id\'], \'sex\': row[\'sex\'], \'batch\': \'\', \'age\': \'\', \'diabetes\': \'\', \'statin\': \'\'}\n        result = api.process_sample(counts, genes, barcodes, sample=row[\'sample_id\'], group=row[\'group\'], outdir=output, model_path=model, dataset=DATASET, source=source, donor_status=\'independent_public_patient_id\')\n        del counts, genes, barcodes\n        gc.collect()\n        emit(\'donor_complete\', sample=row[\'sample_id\'], donor=row[\'donor_id\'], result={k: result[k] for k in (\'input_cells\', \'qc_pass\', \'retained\', \'doublets\', \'cell_type_counts\')})\n    pending = []\n    for row in rows:\n        checkpoint = output / \'samples\' / row[\'sample_id\'] / \'COMPLETE_v1.json\'\n        if not checkpoint.exists():\n            pending.append(row)\n            continue\n        receipt_path = output / \'source_receipts\' / f"{row[\'sample_id\']}.json"\n        receipt = json.loads(receipt_path.read_text())\n        if receipt[\'manifest_sha256\'] != design[\'manifest_sha256\'] or receipt[\'sample_id\'] != row[\'sample_id\'] or receipt[\'donor_id\'] != row[\'donor_id\'] or (receipt[\'group\'] != row[\'group\']) or (receipt[\'sex\'] != row[\'sex\']) or (receipt[\'files\'][\'matrix\'][\'url\'] != row[\'matrix_url\']) or (receipt[\'files\'][\'matrix\'][\'bytes\'] != int(row[\'matrix_bytes\'])):\n            raise RuntimeError(f\'Checkpoint differs from frozen design: {checkpoint}\')\n        source = {**receipt, \'donor_id\': row[\'donor_id\'], \'sex\': row[\'sex\'], \'batch\': \'\', \'age\': \'\', \'diabetes\': \'\', \'statin\': \'\'}\n        api.process_sample(None, None, None, sample=row[\'sample_id\'], group=row[\'group\'], outdir=output, model_path=model, dataset=DATASET, source=source, donor_status=\'independent_public_patient_id\')\n        emit(\'donor_reused\', sample=row[\'sample_id\'], donor=row[\'donor_id\'])\n    schedule = {\'fixed_cohort_size\': 44, \'maximum_concurrent_donor_downloads\': 2, \'maximum_analysis_writers\': 1, \'counts_saved_to_disk\': False, \'processing_order\': \'first completed transport; fixed cohort unchanged\', \'failed_transport_policy\': \'pending and retried after every fixed donor is attempted\', \'scientific_QC_or_mapping_failure\': \'stop and investigate\'}\n    write_new_or_verify(output / \'transport_schedule_v3.json\', schedule)\n    for transport_round in range(1, 4):\n        if not pending:\n            break\n        todo = iter(pending)\n        failed = []\n        emit(\'transport_round\', round=transport_round, pending=len(pending))\n        with ThreadPoolExecutor(max_workers=2) as producer:\n            active = {}\n\n            def fill_queue():\n                while len(active) < 2:\n                    try:\n                        row = next(todo)\n                    except StopIteration:\n                        break\n                    active[producer.submit(download_raw_donor, row)] = row\n                    emit(\'donor_download_start\', sample=row[\'sample_id\'], donor=row[\'donor_id\'], round=transport_round)\n            fill_queue()\n            while active:\n                future = next(as_completed(tuple(active)))\n                row = active.pop(future)\n                try:\n                    raw = future.result()\n                except Exception as error:\n                    failed.append(row)\n                    event = {\'sample\': row[\'sample_id\'], \'donor\': row[\'donor_id\'], \'round\': transport_round, \'error\': repr(error), \'status\': \'PENDING_transport_not_excluded\'}\n                    emit(\'donor_transport_pending\', **event)\n                    with (output / \'transport_pending_events_v3.jsonl\').open(\'a\') as handle:\n                        handle.write(json.dumps(event) + \'\\n\')\n                    del future\n                    gc.collect()\n                    fill_queue()\n                    continue\n                del future\n                process_ready(row, raw)\n                del raw\n                gc.collect()\n                fill_queue()\n        pending = failed\n    if pending:\n        raise RuntimeError(\'Full fixed cohort incomplete; no finalization or CAD fit permitted. Pending: \' + \', \'.join((row[\'sample_id\'] for row in pending)))\n    completed = list(output.glob(\'samples/*/COMPLETE_v1.json\'))\n    if len(completed) != 44:\n        raise RuntimeError(f\'Final full-cohort gate failed: {len(completed)} != 44\')\n    summary = api.finalize(output, dataset=DATASET, source=design)\n    emit(\'dataset_complete\', summary=summary)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
