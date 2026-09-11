"""Offline Kumar ingestion: all deposited RNA cells/features, separate CITE counts."""
import argparse
import gc
import gzip
import json
import os
from pathlib import Path
import re
import sys
from tempfile import NamedTemporaryFile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from epilepsy_microglia.download import sha256, write_json
from epilepsy_microglia.io import write_study
from epilepsy_microglia.metadata import standardize_obs
from epilepsy_microglia.validation import validate_adata
from epilepsy_microglia.kumar import matrix_digest, read_triplet

SOFT = 'GSE201048_family.soft.gz'


def direct_geo_value(row, names, target):
    """Return an unambiguous direct GEO characteristic, without interpretation."""
    present = [(name, row['geo_' + name]) for name in names if 'geo_' + name in row]
    values = {value for _, value in present}
    if len(values) > 1:
        raise ValueError(f'Conflicting direct GEO fields for {target}: {present}')
    return next(iter(values)) if values else pd.NA


def parse_geo(path):
    """Return one source-driven metadata record for every GEO SAMPLE block."""
    records = []
    with gzip.open(path, 'rt') as stream:
        text = stream.read()
    for block in text.split('^SAMPLE = ')[1:]:
        row = {'source_sample_id': block.splitlines()[0]}
        for line in block.splitlines():
            if line.startswith('!Sample_title = '):
                row['geo_title'] = line.split(' = ', 1)[1]
            elif line.startswith('!Sample_source_name_ch1 = '):
                row['geo_source_name_ch1'] = line.split(' = ', 1)[1]
            elif line.startswith('!Sample_characteristics_ch1 = '):
                field = line.split(' = ', 1)[1]
                if ': ' not in field:
                    raise ValueError(f'Malformed GEO characteristic: {field}')
                key, value = field.split(': ', 1)
                key = 'geo_' + re.sub(r'\W+', '_', key.replace(' ', '_').lower()).strip('_')
                if key in row and row[key] != value:
                    raise ValueError(f'Conflicting GEO characteristic in {row["source_sample_id"]}: {key}')
                row[key] = value
        missing = {'source_sample_id', 'geo_title', 'geo_source_name_ch1'} - set(row)
        if missing:
            raise ValueError(f'Incomplete GEO sample block: {sorted(missing)}')
        row['sample_id'] = row['geo_title']
        row['brain_region'] = row['geo_source_name_ch1']
        row['donor_id'] = direct_geo_value(row, ('donor', 'donor_id', 'patient', 'patient_id'),
                                           'donor_id')
        row['diagnosis'] = direct_geo_value(row, ('diagnosis', 'disease'), 'diagnosis')
        records.append(row)
    frame = pd.DataFrame(records).set_index('source_sample_id', drop=False)
    if (len(frame) != 11 or not frame.index.is_unique or frame.index.isna().any()
            or frame.sample_id.isna().any() or frame.sample_id.duplicated().any()):
        raise ValueError('Expected 11 unique GEO accessions and sample titles')
    return frame


def categorical_strings(frame):
    for c in frame:
        if pd.api.types.is_object_dtype(frame[c]) or isinstance(frame[c].dtype, pd.StringDtype):
            frame[c] = pd.Categorical(frame[c])
    return frame


def write_cite(adata, path):
    """Publish new ancillary counts atomically; verify identical prior outputs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        prior = ad.read_h5ad(path)
        if (not prior.obs.equals(adata.obs) or not prior.var.equals(adata.var)
                or matrix_digest(prior.X) != matrix_digest(adata.X)
                or 'counts' not in prior.layers or (prior.layers['counts'] != prior.X).nnz):
            raise FileExistsError(f'Existing CITE output differs: {path}')
    else:
        with NamedTemporaryFile(dir=path.parent, suffix='.h5ad', delete=False) as f:
            temporary = Path(f.name)
        try:
            adata.write_h5ad(temporary)
            os.link(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    reopened = ad.read_h5ad(path)
    if (not reopened.obs_names.equals(adata.obs_names) or not reopened.var_names.equals(adata.var_names)
            or matrix_digest(reopened.X) != matrix_digest(adata.X)
            or (reopened.layers['counts'] != reopened.X).nnz):
        raise ValueError('CITE readback mismatch')


def run(root):
    raw, meta = root / 'data/raw/kumar', root / 'data/metadata/kumar'
    outdir = root / 'data/processed/kumar'
    destination = outdir / 'kumar.h5ad'
    if destination.exists():
        raise FileExistsError(destination)
    inventory = json.loads((meta / 'download_inventory.json').read_text())
    sources = pd.read_csv(meta / 'source_triplets.csv')
    required = set(sources.filename) | {SOFT}
    if required != {r['filename'] for r in inventory} or len(inventory) != len(required):
        raise ValueError('Incomplete download inventory')
    for item in inventory:
        path = raw / item['filename']
        if path.stat().st_size != item['size_bytes'] or sha256(path) != item['sha256']:
            raise ValueError(f'Raw source changed: {path}')
    geo = parse_geo(raw / SOFT)
    file_groups = sources.groupby('source_sample_id', sort=False).filename.agg(list)
    triplet_samples = set(file_groups.index)
    unmatched_triplets = sorted(triplet_samples - set(geo.index))
    unmatched_geo = sorted(set(geo.index) - triplet_samples)
    if unmatched_triplets:
        raise ValueError(f'Downloaded GSM triplets without GEO records: {unmatched_triplets}')
    if unmatched_geo:
        raise ValueError(f'GEO samples without expected triplets: {unmatched_geo}')
    if len(geo) != 11 or len(sources) != 33 or not sources.groupby('source_sample_id').size().eq(3).all():
        raise ValueError('Inconsistent 11-sample MTX triplet mapping')
    geo['files'] = geo.index.map(file_groups)
    geo.reset_index(drop=True).to_csv(meta / 'sample_mapping.csv', index=False)
    matrices, observations, audit = [], [], []
    shared_var = None
    for record in geo.to_dict('records'):
        gsm = record['source_sample_id']
        print(f'Reading {gsm} / {record["sample_id"]}', flush=True)
        full, barcodes, features = read_triplet(raw, record['files'])
        if set(features.feature_type) != {'Gene Expression', 'Antibody Capture'}:
            raise ValueError('Unexpected feature modality')
        mask = features.feature_type.eq('Gene Expression').to_numpy()
        genes = features.loc[mask].copy()
        genes.index = pd.Index(genes.source_gene_id, name='gene_id')
        if shared_var is None:
            shared_var = genes
        elif not genes.equals(shared_var):
            raise ValueError('RNA axes differ; refusing feature loss')
        rna = full[:, mask].tocsr()
        obs = standardize_obs(pd.DataFrame(index=barcodes), 'kumar')
        for key, value in record.items():
            if key != 'files':
                obs[key] = value
        if obs.sample_id.isna().any() or obs.source_sample_id.isna().any():
            raise ValueError(f'{gsm}: sample identifiers are missing after GEO mapping')
        if pd.notna(record['donor_id']) and obs.donor_id.isna().any():
            raise ValueError(f'{gsm}: donor identifier was not propagated from GEO')
        for key, value in record.items():
            if key != 'files' and pd.notna(value) and not obs[key].eq(value).all():
                raise ValueError(f'{gsm}: GEO field was not propagated to every cell: {key}')
        obs['source_file'] = 'kumar/' + next(f for f in record['files'] if f.endswith('_matrix.mtx.gz'))
        obs['source_accession'] = 'GSE201048'
        obs['assay'] = 'scRNA-seq'
        obs['platform'] = '10x Genomics Chromium 3 prime v2 + CITE-seq'
        obs['author_annotation_available'] = False
        annotation_audit = {'annotated_cells': 0}
        obs.index = pd.Index([f'{gsm}:{b}' for b in barcodes], name='cell_id')
        obs = categorical_strings(obs)
        antibody_var = features.loc[~mask].copy()
        antibody_var.index = pd.Index(antibody_var.source_gene_id, name='antibody_id')
        cite = ad.AnnData(X=full[:, ~mask].tocsr(), obs=obs.copy(), var=antibody_var)
        cite.layers['counts'] = cite.X
        cite.uns['ingestion'] = dict(expression_type='raw_counts', raw_counts_available=True,
                                    modality='Antibody Capture', source_file=str(obs.source_file.iloc[0]),
                                    transformations='Feature selection and transpose only; no normalization')
        write_cite(cite, outdir / 'cite_seq' / f'{gsm}.h5ad')
        audit.append(dict(sample_id=record['sample_id'], source_sample_id=gsm,
                          donor_id=None if pd.isna(record['donor_id']) else record['donor_id'],
                          cells=len(obs), genes=len(genes), rna_nnz=rna.nnz,
                          total_counts=int(rna.sum()), rna_sha256=matrix_digest(rna),
                          antibody_features=cite.n_vars, antibody_counts=int(cite.X.sum()),
                          antibody_sha256=matrix_digest(cite.X), **annotation_audit))
        matrices.append(rna)
        observations.append(obs)
        del full, cite
        gc.collect()
    x = sparse.vstack(matrices, format='csr')
    obs = categorical_strings(pd.concat(observations))
    result = ad.AnnData(X=x, obs=obs, var=shared_var)
    result.layers['counts'] = result.X
    result.uns['ingestion'] = dict(raw_counts_available=True, expression_type='raw_counts',
        manifest_study_id='Kumar2022', source_files=[f'kumar/{r["filename"]}' for r in inventory],
        metadata_sources=[SOFT, 'source_triplets.csv'],
        author_annotations='No CSV/TSV cell annotations supplied; unavailable under the MTX/tabular-only scope',
        matrix_provenance='Deposited Cell Ranger integer UMI counts from mixed RNA/Antibody Capture triplets',
        transformations='RNA feature selection, transpose, lossless int32 storage, sample stacking, unique source gene symbols; no QC or normalization',
        clinical_provenance='GEO family SOFT title, source_name_ch1 and Sample_characteristics_ch1 fields, broadcast by GSM triplet',
        missing_metadata='GEO does not supply donor/patient, disease/diagnosis, pathology, FCD subtype, mutation status, age or sex for these samples',
        cite_seq='Separate per-sample raw Antibody Capture h5ad files in cite_seq; panels differ',
        paper_post_qc_cells=85780)
    write_json(meta / 'ingestion_audit.json', audit)
    print(f'Writing RNA {result.shape}', flush=True)
    output = write_study(result, study='kumar', project_root=root)
    expected_obs = result.obs.copy()
    expected_var = result.var.copy()
    del result, matrices, observations, x, obs
    gc.collect()
    reopened = ad.read_h5ad(output)
    errors = validate_adata(reopened, study='kumar')
    if errors:
        raise ValueError('\n'.join(errors))
    if not reopened.obs.equals(expected_obs) or not reopened.var.equals(expected_var):
        raise ValueError('Reopened metadata differs')
    offset = 0
    for row in audit:
        end = offset + row['cells']
        if matrix_digest(reopened.X[offset:end]) != row['rna_sha256']:
            raise ValueError('Reopened counts differ from source')
        section = reopened.obs.iloc[offset:end]
        for key in ('sample_id', 'source_sample_id'):
            if not section[key].eq(row[key]).all():
                raise ValueError('Reopened sample/donor alignment differs')
        offset = end
    if (offset != reopened.n_obs or reopened.obs.sample_id.isna().any()
            or reopened.obs.source_sample_id.isna().any()
            or reopened.obs.sample_id.nunique() != 11 or reopened.obs.source_sample_id.nunique() != 11):
        raise ValueError('Cells or 11-sample mapping lost')
    for item in inventory:
        if sha256(raw / item['filename']) != item['sha256']:
            raise ValueError('Raw source changed during ingestion')
    summary = dict(validation='passed', shape=list(reopened.shape),
                   X_dtype=str(reopened.X.dtype), X_type=type(reopened.X).__name__, X_nnz=reopened.X.nnz,
                   total_counts=int(reopened.X.sum()),
                   sample_counts={str(k): int(v) for k,v in reopened.obs.sample_id.value_counts().items()},
                   samples=11, donor_count=(int(reopened.obs.donor_id.nunique())
                                            if reopened.obs.donor_id.notna().any() else None),
                   geo_metadata_fields=[c for c in geo.columns if c != 'files'],
                   unmatched_triplet_samples=len(unmatched_triplets), unmatched_geo_samples=len(unmatched_geo),
                   annotated_cells=int(reopened.obs.author_annotation_available.sum()),
                   raw_sha256_unchanged=True, all_source_counts_verified_after_reopening=True,
                   metadata_equal_after_reopening=True, cite_outputs_verified=11)
    write_json(outdir / 'validation_summary.json', summary)
    write_json(meta / 'validation_summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=ROOT)
    run(parser.parse_args().project_root.resolve())
