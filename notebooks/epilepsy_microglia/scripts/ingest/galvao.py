"""Faithful Galvao GSE268807 RNA ingestion from inspected Cell Ranger ARC triplets.

All source barcodes and RNA features survive. Peaks stay in immutable raw files.
CELLxGENE is used for annotations and count cross-checks, not cell selection.
No network, normalization, QC filtering, or integration.
"""
from pathlib import Path
import argparse
import gc
import gzip
import hashlib
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from epilepsy_microglia.metadata import standardize_obs
from epilepsy_microglia.validation import validate_adata
from epilepsy_microglia.io import write_study

AUTHOR_FILE = '799e0777-8fe4-4d14-9183-01a7314923ec.h5ad'
# Verified by unique full barcode containment, GEO donor/histology/lobe and
# subsequently by equality of every shared gene count with CELLxGENE raw.X.
# Filename specimen, GEO source_name_ch1, CELLxGENE observation suffix.
SPECIMENS = [
    ('G120_D_FL', 'G120_D.3', '3'),
    ('G120_D_TL', 'G120_D.2', '2'),
    ('G120_F1_N', 'G120_N.1', '1'),
    ('G129_D', 'G129_D', '6'),
    ('G133_D_FL', 'G133_D.2', '4'),
    ('G133_N_FL', 'G133_N.2', '5'),
    ('G150_D', 'G150_D', '7'),
    ('G159_D', 'G159_D', '8'),
    ('G171_D', 'G171_D', '9'),
    ('G187_D', 'G187_D', '10'),
    ('G210_D', 'G210_D', '11'),
]
FEATURE_COLUMNS = ['source_gene_id', 'gene_symbol', 'feature_type',
                   'chromosome', 'source_start', 'source_end']


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def matrix_digest(matrix):
    """Exact CSR fingerprint independent of index/int storage width."""
    digest = hashlib.sha256()
    for array in (np.array(matrix.shape), matrix.indptr, matrix.indices, matrix.data):
        digest.update(np.asarray(array, dtype='<i8').tobytes())
    return digest.hexdigest()


def read_geo(path):
    samples = []
    with gzip.open(path, 'rt') as f:
        text = f.read()
    for block in text.split('^SAMPLE = ')[1:]:
        row = {'accession': block.splitlines()[0]}
        for line in block.splitlines():
            if line.startswith('!Sample_characteristics_ch1 = '):
                key, value = line.split(' = ', 1)[1].split(': ', 1)
                row[key] = value
            elif line.startswith('!Sample_source_name_ch1 = '):
                row['source_name'] = line.split(' = ', 1)[1]
            elif line.startswith('!Sample_platform_id = '):
                row['platform_id'] = line.split(' = ', 1)[1]
        samples.append(row)
    return samples


def run(root):
    raw = root / 'data/raw/galvao'
    meta = root / 'data/metadata/galvao'
    outdir = root / 'data/processed/galvao'
    destination = outdir / 'galvao.h5ad'
    if destination.exists():
        raise FileExistsError(destination)
    inventory = json.loads((meta / 'download_inventory.json').read_text())
    for item in inventory:
        file = raw / item['filename']
        if file.stat().st_size != item['size_bytes'] or sha256(file) != item['sha256']:
            raise ValueError(f'Raw original changed: {file}')
    geo = read_geo(raw / 'GSE268807_family.soft.gz')
    rna_geo = {r['source_name']: r for r in geo if r['library type'] == 'mRNA'}
    atac_geo = {r['source_name']: r for r in geo if r['library type'] == 'ATAC'}
    author = ad.read_h5ad(raw / AUTHOR_FILE, backed='r')
    author_obs = author.obs.copy()
    if not author_obs.index.is_unique:
        raise ValueError('Duplicate author cell IDs')
    suffixes = author_obs.index.str.rsplit('_', n=1).str[-1]
    if set(suffixes) != {s[2] for s in SPECIMENS}:
        raise ValueError('Unexpected author specimen suffixes')
    matrices, observations, audit = [], [], []
    shared_var = None
    feature_annotation = author.raw.var.copy()
    for file_specimen, sample, suffix in SPECIMENS:
        prefix = 'GSE268807_' + file_specimen
        features = pd.read_csv(raw / f'{prefix}_features.tsv.gz', sep='\t',
                               header=None, names=FEATURE_COLUMNS, dtype=str, keep_default_na=False)
        barcodes = pd.Index(pd.read_csv(raw / f'{prefix}_barcodes.tsv.gz',
                                      sep='\t', header=None, dtype=str, keep_default_na=False)[0])
        if not barcodes.is_unique or not features.source_gene_id.is_unique:
            raise ValueError('Duplicate source feature/barcode identifiers')
        gene_mask = features.feature_type.eq('Gene Expression').to_numpy()
        n_genes = int(gene_mask.sum())
        if n_genes != 36601 or not gene_mask[:n_genes].all() or gene_mask[n_genes:].any():
            raise ValueError('Source feature layout differs from inspection')
        if set(features.feature_type) != {'Gene Expression', 'Peaks'}:
            raise ValueError('Unexpected modality')
        genes = features.loc[gene_mask].copy()
        genes.index = pd.Index(genes.source_gene_id.to_numpy(), name="gene_id")
        if shared_var is None:
            shared_var = genes
        elif not genes.equals(shared_var):
            raise ValueError('Gene axes differ; refusing implicit feature intersection')
        record = rna_geo[sample]
        ids = np.flatnonzero(suffixes == suffix)
        annotations = author_obs.iloc[ids].copy()
        author_ids = annotations.index.copy()
        annotations.index = annotations.index.str.rsplit('_', n=1).str[0]
        if not annotations.index.is_unique or not annotations.index.isin(barcodes).all():
            raise ValueError(f'Invalid author barcode mapping for {sample}')
        for source, target in [('donor', 'donor_id'), ('histological assessment', 'histological_assessment'),
                               ('lateralization', 'lateralization')]:
            if not annotations[target].eq(record[source]).all():
                raise ValueError(f'GEO/CELLxGENE metadata conflict: {sample}, {source}')
        annotations['source_obs_id'] = author_ids.to_numpy()
        annotations = annotations.add_prefix('cellxgene_').reindex(barcodes)
        obs = standardize_obs(pd.DataFrame(index=barcodes), 'galvao')
        obs['sample_id'] = sample
        obs['donor_id'] = record['donor']
        obs['source_sample_id'] = record['accession']
        obs['paired_atac_accession'] = atac_geo[sample]['accession']
        obs['source_file_sample_id'] = file_specimen
        obs['source_accession'] = 'GSE268807'
        obs['source_file'] = f'galvao/{prefix}_matrix.mtx.gz'
        obs['diagnosis'] = record['disease']
        obs['pathology'] = record['histological assessment']
        subtype = {'FCD_typeIIa': 'IIa', 'FCD_typeIIb': 'IIb'}.get(record['histological assessment'])
        obs['fcd_type'] = pd.Categorical([subtype] * len(obs))
        obs['fcd_subtype'] = obs['fcd_type'].copy()
        obs['control_status'] = ('internal_histologically_normal' if record['histological assessment'] == 'Normal' else 'FCD_lesion')
        obs['brain_region'] = record['lobe']
        obs['assay'] = 'snRNA-seq'
        obs['platform'] = '10x Genomics Multiome ATAC + Gene Expression'
        for key, value in record.items():
            obs['geo_' + key.replace(' ', '_')] = value
        obs = pd.concat([obs, annotations], axis=1)
        obs['author_annotation_available'] = obs.cellxgene_source_obs_id.notna()
        obs.index = pd.Index([f'{sample}:{b}' for b in barcodes], name='cell_id')
        print(f'Reading {prefix}: {len(barcodes):,} nuclei', flush=True)
        with threadpool_limits(limits=4):
            full = mmread(raw / f'{prefix}_matrix.mtx.gz')
        if full.shape != (len(features), len(barcodes)) or full.dtype.kind not in 'iu':
            raise ValueError('Matrix dimensions/type differ from inspected axes')
        if (full.data < 0).any():
            raise ValueError('Negative source counts')
        keep = full.row < n_genes
        gene_entries = int(keep.sum())
        rna = sparse.coo_matrix((full.data[keep], (full.col[keep], full.row[keep])),
                                shape=(len(barcodes), n_genes)).tocsr()
        if rna.nnz != gene_entries:
            raise ValueError('Duplicate source coordinates; refusing implicit summation')
        if rna.data.max() > np.iinfo(np.int32).max:
            raise ValueError('Counts exceed int32 storage range')
        total_features, total_entries = full.shape[0], full.nnz
        del full, keep
        rna = rna.astype(np.int32)
        rna.sort_indices()
        # All author shared count values, not a sampled barcode-only check.
        source_rows = barcodes.get_indexer(author_ids.str.rsplit('_', n=1).str[0])
        source_cols = shared_var.index.get_indexer(author.raw.var_names)
        common = source_cols >= 0
        # The current CELLxGENE object has 18 IDs absent from GEO; compare exact
        # shared IDs only, without renaming any GEO genes.
        author_counts = author.raw.X[ids, :][:, common]
        mismatch = int((rna[source_rows][:, source_cols[common]] != author_counts).nnz)
        if mismatch:
            raise ValueError(f'{sample}: {mismatch} mismatched shared counts')
        del author_counts
        audit.append(dict(sample_id=sample, file_specimen=file_specimen,
                          donor_id=record['donor'], source_sample_id=record['accession'],
                          cellxgene_suffix=suffix, nuclei=len(barcodes), genes=n_genes,
                          peak_features=total_features - n_genes, source_entries=total_entries,
                          rna_nonzero_entries=rna.nnz, total_rna_counts=int(rna.sum()),
                          rna_csr_sha256=matrix_digest(rna), annotated_nuclei=len(ids),
                          shared_count_mismatches=mismatch, shared_genes_compared=int(common.sum()),
                          author_gene_ids_not_in_geo=author.raw.var_names[~common].tolist()))
        matrices.append(rna)
        observations.append(obs)
        print(f'  {rna.nnz:,} RNA entries; all shared counts match for {len(ids):,} annotated nuclei', flush=True)
        gc.collect()
    author.file.close()
    x = sparse.vstack(matrices, format='csr')
    del matrices, rna
    obs = pd.concat(observations)
    del observations
    # Preserve NA safely in AnnData 0.11 without stringifying missing values.
    for column in obs.select_dtypes(include=['object']).columns:
        obs[column] = pd.Categorical(obs[column])
    # Original feature columns come from GEO. Add only the author's gene
    # annotations, not their derived HVG statistics; full source remains in raw.
    feature_annotation = feature_annotation.add_prefix('cellxgene_').reindex(shared_var.index)
    var = pd.concat([shared_var, feature_annotation], axis=1)
    for column in var.select_dtypes(include=['object']).columns:
        var[column] = pd.Categorical(var[column])
    result = ad.AnnData(X=x, obs=obs, var=var)
    result.layers['counts'] = result.X
    result.uns['ingestion'] = dict(
        raw_counts_available=True, expression_type='raw_counts',
        source_files=[f'galvao/{r["filename"]}' for r in inventory],
        metadata_sources=['GSE268807_family.soft.gz', AUTHOR_FILE, 'cellxgene_collection.json'],
        manifest_study_id='Galvao2024',
        matrix_provenance='Cell Ranger ARC 2.0.2 integer gene-expression counts; all deposited barcodes and RNA features retained',
        cellxgene_usage='Annotation left join by verified specimen suffix and original barcode; every shared raw count matched',
        atac_provenance='Peaks present in original mixed-feature MTX triplets, preserved unchanged in raw; RNA-only output',
        transformations='None; explicit RNA feature selection, transpose, lossless int32 storage, specimen stacking and source gene-symbol indices with unique suffixes',
        missing_metadata='No mutation status supplied; fcd_type missing for histologically normal internal specimens',
    )
    outdir.mkdir(parents=True, exist_ok=True)
    (meta / 'ingestion_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    pd.DataFrame(audit).to_csv(meta / 'sample_mapping.csv', index=False)
    print(f'Writing {result.shape}', flush=True)
    output = write_study(result, study='galvao', project_root=root)
    del result, x, obs, var
    gc.collect()
    reopened = ad.read_h5ad(output)
    errors = validate_adata(reopened, study='galvao')
    if errors:
        raise ValueError('\n'.join(errors))
    offset = 0
    for row in audit:
        n = row['nuclei']
        if matrix_digest(reopened.X[offset:offset+n]) != row['rna_csr_sha256']:
            raise ValueError('Reopened counts differ from source RNA matrix')
        offset += n
    if offset != reopened.n_obs:
        raise ValueError('Source nuclei count not preserved')
    for item in inventory:
        if sha256(raw / item['filename']) != item['sha256']:
            raise ValueError('Raw source changed during ingestion')
    summary = dict(path=str(output), size_bytes=output.stat().st_size,
                   shape=list(reopened.shape), X_type=type(reopened.X).__name__,
                   X_dtype=str(reopened.X.dtype), X_nnz=reopened.X.nnz,
                   layers=list(reopened.layers), obs_columns=list(reopened.obs),
                   var_columns=list(reopened.var),
                   sample_counts={str(k):int(v) for k,v in reopened.obs.sample_id.value_counts().items()},
                   donor_counts={str(k):int(v) for k,v in reopened.obs.donor_id.value_counts().items()},
                   annotated_nuclei=int(reopened.obs.author_annotation_available.sum()),
                   values={c:{str(k):int(v) for k,v in reopened.obs[c].value_counts(dropna=False).items()}
                           for c in ['diagnosis','fcd_type','pathology','brain_region','control_status']},
                   validation='passed; reopened, source counts fingerprinted, raw SHA256 unchanged')
    (outdir / 'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    reopened.obs.head().to_csv(outdir / 'obs_head.tsv',sep='\t')
    print(json.dumps(summary,indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=ROOT)
    args = parser.parse_args()
    run(args.project_root.resolve())


if __name__ == '__main__':
    main()
