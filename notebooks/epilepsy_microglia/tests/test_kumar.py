import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from epilepsy_microglia.download import fetch, sha256, verify
from epilepsy_microglia.kumar import clinical_metadata, matrix_digest, read_triplet

ROOT = Path(__file__).resolve().parents[1]


def test_completed_raw_is_never_replaced(tmp_path):
    raw, staging = tmp_path / 'raw', tmp_path / 'stage'
    raw.mkdir()
    staging.mkdir()
    path = raw / 'a.gz'
    path.write_bytes(gzip.compress(b'original'))
    before = path.read_bytes()
    with pytest.raises(ValueError, match='Inventory mismatch'):
        fetch(dict(filename='a.gz', url='https://invalid'), raw, staging,
              dict(sha256='wrong', size_bytes=path.stat().st_size))
    assert path.read_bytes() == before
    assert not list(staging.iterdir())


def test_gzip_corruption_rejected(tmp_path):
    path = tmp_path / 'bad.gz'
    path.write_bytes(gzip.compress(b'counts')[:-5])
    with pytest.raises((EOFError, OSError)):
        verify(path)


def test_download_resumes_then_preserves_completed_file(tmp_path, monkeypatch):
    raw, staging = tmp_path / 'raw', tmp_path / 'stage'
    raw.mkdir()
    staging.mkdir()
    payload = gzip.compress(b'count matrix')
    partial = staging / 'counts.mtx.gz'
    partial.write_bytes(payload[:8])
    calls = []

    def curl(command, check):
        assert check and command[command.index('--continue-at') + 1] == '-'
        assert partial.read_bytes() == payload[:8]
        calls.append(command)
        with partial.open('ab') as f:
            f.write(payload[8:])

    monkeypatch.setattr('epilepsy_microglia.download.subprocess.run', curl)
    source = dict(filename=partial.name, url='https://example.org/counts.mtx.gz', expected_size_bytes=len(payload))
    item = fetch(source, raw, staging)
    final = raw / partial.name
    before = (final.stat().st_ino, final.stat().st_mtime_ns)
    assert item['sha256'] == sha256(final) and not partial.exists()
    assert fetch(source, raw, staging, item) == item
    assert (final.stat().st_ino, final.stat().st_mtime_ns) == before
    assert len(calls) == 1


def test_download_allowlist_matches_sample_keys():
    spec = importlib.util.spec_from_file_location('kumar_download', ROOT / 'scripts/download/kumar.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sources = module.source_list(ROOT)
    samples = pd.read_csv(ROOT / 'data/metadata/kumar/source_samples.csv')
    assert {r['source_sample_id'] for r in sources} == set(samples.source_sample_id)
    for row in samples.itertuples():
        names = {r['filename'] for r in sources if r['source_sample_id'] == row.source_sample_id}
        assert names == {row.file_prefix + '_' + suffix for suffix in
                         ['matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz']}


def test_clinical_keys_and_missing_values():
    meta = ROOT / 'data/metadata/kumar'
    samples = pd.read_csv(meta / 'curated_samples.csv')
    geo = pd.DataFrame({'source_sample_id': samples.sample_id, 'sample_id': samples.assay_id,
                        'geo_brain_region': samples.brain_region_raw, 'files': [[]]*len(samples)})
    result = clinical_metadata(meta, geo)
    assert len(result) == 11 and result.donor_id.nunique() == 6
    assert result.loc[result.donor_id.isin(['P2', 'P3']), 'fcd_subtype'].eq('IIb').all()
    assert result.mutation_status.isna().all()
    assert result.epilepsy_duration_years.isna().all()
    assert not result.control_status.str.contains('healthy').any()
    geo.loc[0, 'source_sample_id'] = 'GSM6049643'
    with pytest.raises(ValueError, match='accessions differ'):
        clinical_metadata(meta, geo)


def test_wrong_sample_title_rejected():
    samples = pd.read_csv(ROOT / 'data/metadata/kumar/curated_samples.csv')
    geo = pd.DataFrame({'source_sample_id': samples.sample_id, 'sample_id': samples.assay_id,
                        'files': [[]]*len(samples)})
    geo.loc[0, 'sample_id'] = 'P6.C'
    with pytest.raises(ValueError, match='title'):
        clinical_metadata(ROOT / 'data/metadata/kumar', geo)


def test_mixed_triplet_retains_zero_cells_and_features(tmp_path):
    contents = {'x_features.tsv.gz': 'ENSG000001\tG\tGene Expression\nCD3\tCD3\tAntibody Capture\n',
                'x_barcodes.tsv.gz': 'a\nb\nc\n',
                'x_matrix.mtx.gz': '%%MatrixMarket matrix coordinate integer general\n2 3 2\n1 1 4\n2 2 9\n'}
    for name, value in contents.items():
        (tmp_path / name).write_bytes(gzip.compress(value.encode()))
    x, barcodes, features = read_triplet(tmp_path, list(contents))
    assert sparse.isspmatrix_csr(x) and x.dtype == np.int32
    assert x.shape == (3, 2) and x.sum() == 13
    assert barcodes.tolist() == ['a', 'b', 'c']
    assert features.feature_type.tolist() == ['Gene Expression', 'Antibody Capture']
    assert matrix_digest(x) == matrix_digest(x.astype(np.int64))
