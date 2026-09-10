"""Inspected Kumar GEO layout and explicit clinical joins (offline only)."""
import gzip
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from .metadata import join_metadata


def matrix_digest(matrix):
    matrix = matrix.tocsr(copy=True)
    matrix.sort_indices()
    h = hashlib.sha256()
    for a in (matrix.shape, matrix.indptr, matrix.indices, matrix.data):
        h.update(np.asarray(a, dtype='<i8').tobytes())
    return h.hexdigest()


def clinical_metadata(meta, geo):
    samples = pd.read_csv(meta / 'curated_samples.csv', index_col=0)
    tissues = pd.read_csv(meta / 'curated_tissues.csv', index_col=0)
    donors = pd.read_csv(meta / 'clinical_donors.csv')
    if len(samples) != 11 or samples.donor_id_source.nunique() != 6:
        raise ValueError('Expected 11 curated samples and 6 donors')
    if set(samples.sample_id) != set(geo.source_sample_id):
        raise ValueError('Curated and GEO accessions differ')
    mapping = geo.drop(columns='files').merge(samples, left_on='source_sample_id', right_on='sample_id',
                                             validate='one_to_one', suffixes=('', '_curated'))
    if not mapping.sample_id.eq(mapping.assay_id).all():
        raise ValueError('GEO title does not match curated assay key')
    mapping['tissue_id'] = 'Kumar2022_' + mapping.sample_id
    mapping = mapping.merge(tissues, on=['tissue_id', 'donor_id_global'], how='left',
                            validate='one_to_one', suffixes=('', '_tissue'))
    if mapping.histology_notes.isna().any():
        raise ValueError('Unmatched curated tissue key')
    mapping = join_metadata(mapping, donors, key='donor_id_source')
    if mapping.surgery_reported.isna().any():
        raise ValueError('Unmatched clinical donor key')
    fcd = mapping.donor_id_source.isin(['P2', 'P3'])
    if not mapping.loc[fcd, 'fcd_subtype'].eq('IIb').all() or mapping.loc[~fcd, 'fcd_subtype'].notna().any():
        raise ValueError('Clinical FCD mapping differs from Table S5')
    mapping['donor_id'] = mapping.donor_id_source
    mapping['diagnosis'] = mapping.diagnosis_raw
    mapping['pathology'] = mapping.histology_notes
    mapping['control_status'] = 'epilepsy_non_FCD'
    mapping.loc[fcd, 'control_status'] = 'FCD_lesion'
    # Unknown genetic testing is not a negative finding. Keep curation verbatim
    # in its original columns; standardized mutation_status is missing.
    mapping['mutation_status'] = np.nan
    return mapping.set_index('source_sample_id', drop=False)


def read_triplet(archive, files):
    def member(suffix):
        matches = [f for f in files if f.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f'Expected one {suffix}')
        return matches[0]
    def table(suffix):
        with (archive / member(suffix)).open('rb') as f, gzip.open(f, 'rt') as stream:
            return pd.read_csv(stream, sep='\t', header=None, dtype=str, keep_default_na=False)
    features = table('_features.tsv.gz')
    if features.shape[1] != 3:
        raise ValueError('Expected three feature columns')
    features.columns = ['source_gene_id', 'gene_symbol', 'feature_type']
    barcodes = pd.Index(table('_barcodes.tsv.gz')[0], name='source_cell_id')
    if not barcodes.is_unique or not features.source_gene_id.is_unique:
        raise ValueError('Duplicate source axes')
    with (archive / member('_matrix.mtx.gz')).open('rb') as f, gzip.open(f, 'rb') as stream:
        with threadpool_limits(limits=4):
            matrix = mmread(stream)
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError('Matrix/axis dimensions differ')
    if matrix.dtype.kind not in 'iu' or (matrix.data < 0).any() or matrix.data.max(initial=0) > np.iinfo(np.int32).max:
        raise ValueError('Invalid source counts')
    x = matrix.T.tocsr()
    if x.nnz != matrix.nnz:
        raise ValueError('Duplicate source coordinates')
    return x.astype(np.int32), barcodes, features
