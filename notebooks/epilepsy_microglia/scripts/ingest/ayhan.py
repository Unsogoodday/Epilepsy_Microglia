"""Offline, lossless Ayhan count ingestion. No normalization, filtering or integration.

CSV is genes x nuclei despite GEO's generic opposite-orientation description.
Library prefixes are joined exactly to GEO titles, never assigned by cell order.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts/download'))
from ayhan import sha256
from epilepsy_microglia.metadata import standardize_obs
from epilepsy_microglia.io import write_study
from epilepsy_microglia.validation import validate_adata

MATRIX = 'GSE160189_Hippo_Counts.csv.gz'


def read_geo(path):
    records = []
    with gzip.open(path, 'rt') as f:
        text = f.read()
    for block in text.split('^SAMPLE = ')[1:]:
        row = {'source_sample_id': block.splitlines()[0]}
        for line in block.splitlines():
            if line.startswith('!Sample_characteristics_ch1 = '):
                field = line.split(' = ', 1)[1]
                if ': ' not in field:
                    raise ValueError(f'Malformed GEO characteristic: {field}')
                key, val = field.split(': ', 1)
                key = 'geo_'+re.sub(r'\W+', '_', key.lower()).strip('_')
                if key in row and row[key] != val:
                    raise ValueError(f'Conflicting GEO characteristic in {row["source_sample_id"]}: {key}')
                row[key] = val
            elif line.startswith('!Sample_title = '):
                row['geo_title'] = line.split(' = ', 1)[1]
            elif line.startswith('!Sample_source_name_ch1 = '):
                row['geo_source_name_ch1'] = line.split(' = ', 1)[1]
                row['donor_id'] = row['geo_source_name_ch1']
        missing = {'source_sample_id', 'geo_title', 'geo_source_name_ch1', 'geo_tissue'} - set(row)
        if missing:
            raise ValueError(f'Incomplete GEO sample block: {sorted(missing)}')
        match = re.fullmatch(r'([AP]\d+)_(Donor\d+)_scRNA-seq', row['geo_title'])
        if not match or match[2] != row['donor_id']:
            raise ValueError('Unexpected GEO title/donor mapping')
        row['sample_id'] = match[1]
        row['brain_region'] = row['geo_tissue']
        row['hippocampal_region'] = {'A': 'anterior', 'P': 'posterior'}[match[1][0]]
        if not row['brain_region'].lower().startswith(row['hippocampal_region']):
            raise ValueError('GEO tissue/prefix conflict')
        records.append(row)
    frame = pd.DataFrame(records).set_index('sample_id', drop=False)
    if (len(frame) != 10 or not frame.index.is_unique or frame.source_sample_id.isna().any()
            or frame.source_sample_id.duplicated().any() or frame.donor_id.nunique() != 5):
        raise ValueError('Expected ten libraries from five donors')
    for _, group in frame.groupby('donor_id'):
        if set(group.hippocampal_region) != {'anterior', 'posterior'}:
            raise ValueError('Incomplete donor pair')
        for c in ['geo_age_yr', 'geo_sex', 'geo_race', 'geo_hemisphere', 'geo_epilepsy_duration_yr']:
            if group[c].nunique() != 1:
                raise ValueError(f'Conflicting donor metadata: {c}')
    return frame


def read_counts(path):
    genes, data, indices, indptr = [], [], [], [0]
    with gzip.open(path, 'rt') as f:
        header = next(csv.reader([next(f)]))
        cells = pd.Index(header[1:], name='source_cell_id')
        if header[0] != 'gene' or not cells.is_unique or (cells.str.len() == 0).any():
            raise ValueError('Unexpected CSV header or duplicate/empty cell IDs')
        for n, line in enumerate(f, 1):
            gene, sep, values = line.rstrip('\r\n').partition(',')
            if not sep or not gene or '"' in gene:
                raise ValueError(f'Unexpected gene record at row {n}')
            # Strict lexical check prevents fromstring silently accepting malformed text.
            if re.search(r'[^0-9,]', values) or ',,' in values or values.endswith(','):
                raise ValueError(f'Non-integer or missing count at row {n}')
            arr = np.fromstring(values, sep=',', dtype=np.int64)
            if arr.size != len(cells) or (arr < 0).any() or arr.max() > np.iinfo(np.int32).max:
                raise ValueError(f'Invalid count range/dimensions at row {n}')
            nz = np.flatnonzero(arr)
            genes.append(gene)
            data.append(arr[nz].astype(np.int32))
            indices.append(nz.astype(np.int32))
            indptr.append(indptr[-1]+len(nz))
            if n % 5000 == 0:
                print(f'{n:,} genes parsed; {indptr[-1]:,} nonzero counts', flush=True)
    if not pd.Index(genes).is_unique:
        raise ValueError('Duplicate source gene identifiers')
    gx = sparse.csr_matrix((np.concatenate(data), np.concatenate(indices), np.asarray(indptr)),
                           shape=(len(genes), len(cells)))
    return gx.T.tocsr(), cells, genes


def fingerprint(x):
    h = hashlib.sha256()
    for a in (np.asarray(x.shape), x.indptr, x.indices, x.data):
        h.update(np.asarray(a, dtype='<i8').tobytes())
    return h.hexdigest()


def run(root):
    raw, meta = root/'data/raw/ayhan', root/'data/metadata/ayhan'
    inventory_path = meta/'download_inventory.json'
    if not inventory_path.is_file():
        raise FileNotFoundError(
            f'{inventory_path} is missing. Complete the download first with: '
            'python -u ./download/ayhan.py')
    inventory = json.loads(inventory_path.read_text())
    required = {MATRIX, 'GSE160189_family.soft.gz'}
    if not required.issubset({r['filename'] for r in inventory}):
        raise ValueError('Download inventory is incomplete; run scripts/download/ayhan.py')
    for r in inventory:
        if sha256(raw/r['filename']) != r['sha256']:
            raise ValueError(f'Raw checksum mismatch: {r["filename"]}')
    samples = read_geo(raw/'GSE160189_family.soft.gz')
    geo_metadata_fields = list(samples.columns)
    x, cells, genes = read_counts(raw/MATRIX)
    prefixes = cells.str.split('_').str[0]
    unmatched_cells = int((~prefixes.isin(samples.index)).sum())
    unmatched_geo = sorted(set(samples.index) - set(prefixes))
    unexpected_prefixes = sorted(set(prefixes) - set(samples.index))
    if unmatched_cells or unexpected_prefixes:
        raise ValueError(f'Cell prefixes without exactly one GEO sample: {unexpected_prefixes}')
    if unmatched_geo:
        raise ValueError(f'GEO samples without cells: {unmatched_geo}')
    obs = samples.loc[prefixes].copy()
    obs.index = cells
    obs = standardize_obs(obs, 'ayhan')
    obs['assay'] = 'snRNA-seq'
    obs['platform'] = '10x Genomics Chromium; Illumina NovaSeq 6000'
    obs['source_accession'] = 'GSE160189'
    obs['source_file'] = 'ayhan/'+MATRIX
    obs['author_annotation_available'] = False
    for src, dst in [('geo_age_yr','age_years'), ('geo_epilepsy_duration_yr','epilepsy_duration_years'), ('geo_rin','rin')]:
        obs[dst] = pd.to_numeric(obs[src], errors='raise')
    obs['sex'] = obs.geo_sex
    obs['hemisphere'] = obs.geo_hemisphere
    obs['batch'] = obs.geo_batch
    if obs.sample_id.isna().any() or obs.source_sample_id.isna().any():
        raise ValueError('Sample identifiers are missing after GEO mapping')
    if samples.donor_id.notna().any() and obs.donor_id.isna().any():
        raise ValueError('Donor identifiers are missing after GEO mapping')
    obs.index = pd.Index('ayhan:'+cells, name='cell_id')
    for col in obs.select_dtypes('object'):
        obs[col] = pd.Categorical(obs[col])
    var = pd.DataFrame({'source_gene_id': genes, 'gene_symbol': genes}, index=pd.Index(genes, name='gene_name'))
    result = ad.AnnData(x, obs=obs, var=var)
    result.layers['counts'] = result.X
    result.uns['ingestion'] = dict(raw_counts_available=True, expression_type='raw_counts',
        source_files=[f'ayhan/{r["filename"]}' for r in inventory],
        metadata_sources=['GSE160189_family.soft.gz'],
        manifest_study_id='Ayhan2021', genome_build='hg19',
        matrix_provenance='GEO: Cell Ranger 3.0.2 count matrices; deposited gene-by-cell CSV',
        transformations='Transpose and lossless sparse int32 storage only; all deposited cells and genes retained',
        cell_mapping='Exact CSV cell prefix to GEO sample title; paired donors from GEO source_name',
        missing_metadata='Diagnosis, pathology, FCD subtype, mutation status and cell-type annotation are not supplied in GEO sample metadata')
    digest = fingerprint(x)
    samples['nuclei'] = pd.Series(prefixes.value_counts())
    samples.to_csv(meta/'sample_metadata.csv', index=False)
    print(f'Writing {result.shape}', flush=True)
    output = root/'data/processed/ayhan/ayhan.h5ad'
    if output.exists():
        print('Existing output: comparing all counts and metadata against freshly parsed sources', flush=True)
    else:
        output = write_study(result, study='ayhan', project_root=root)
    reopened = ad.read_h5ad(output)
    errors = validate_adata(reopened, study='ayhan')
    if errors or fingerprint(reopened.X) != digest or not reopened.obs_names.equals(result.obs_names) or not reopened.var_names.equals(result.var_names):
        raise ValueError(f'Round-trip validation failed: {errors}')
    pd.testing.assert_frame_equal(reopened.obs, result.obs)
    pd.testing.assert_frame_equal(reopened.var, result.var)
    summary = dict(shape=list(x.shape), nnz=x.nnz, total_counts=int(x.sum()),
                   count_csr_sha256=digest, sample_counts={str(k):int(v) for k,v in prefixes.value_counts().items()},
                   donors=int(obs.donor_id.nunique()), samples=int(obs.sample_id.nunique()),
                   path=str(output), count_dtype=str(reopened.X.dtype), count_format=reopened.X.format,
                   duplicate_cell_ids=0, duplicate_gene_ids=0, unmatched_cells=unmatched_cells,
                   unmatched_geo_samples=len(unmatched_geo), geo_metadata_fields=geo_metadata_fields,
                   excluded_cells=0, excluded_genes=0, source_checksums_verified=True,
                   validation='passed: full sparse count fingerprint, axes, schema and counts layer after reopening')
    (output.parent/'validation_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (output.parent/'dataset_report.md').write_text(
        '# Ayhan preparation report\n\n'
        f'Verified {reopened.n_obs:,} nuclei × {reopened.n_vars:,} genes; '
        f'{summary["samples"]} libraries and {summary["donors"]} donors.\n\n'
        'All deposited cells and genes retained. No duplicate identifiers, unmatched '
        'cell/library keys, or exclusions. X and counts are equal sparse int32 CSR matrices. '
        'Original cell IDs and gene symbols are preserved; gene IDs in this deposit are symbols. '
        'Metadata is joined by explicit library prefix to the GEO family SOFT sample record.\n\n'
        'GEO describes Cell Ranger 3.0.2 counts. Its format description incorrectly says '
        'tab-delimited cells-by-genes: the file is comma-delimited genes-by-cells. '
        '131,325 nuclei agrees with GEO; the paper reports 129,908 after author QC '
        '(1,417 fewer). No author QC was reapplied. No GEO diagnosis, cell-type annotations, FCD subtype, mutation status, '
        'or pathology fields are supplied, so these remain missing.\n\n'
        'See validation_summary.json for full count fingerprint and sample totals, '
        '../../metadata/ayhan/download_inventory.json for file sizes and SHA256, '
        'and ../../../scripts/AYHAN_README.md for source references and reproduction commands.\n')
    print(json.dumps(summary, indent=2))
    print(f'Ingestion complete: {output}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=ROOT)
    try:
        run(p.parse_args().project_root.resolve())
    except FileNotFoundError as error:
        p.error(str(error))
