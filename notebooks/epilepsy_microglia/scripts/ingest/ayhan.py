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
import zipfile
import xml.etree.ElementTree as ET

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts/download'))
from ayhan import sha256, CLINICAL
from epilepsy_microglia.metadata import standardize_obs
from epilepsy_microglia.io import write_study
from epilepsy_microglia.validation import validate_adata

MATRIX = 'GSE160189_Hippo_Counts.csv.gz'


def read_clinical(path, samples):
    """Read the inspected Table S1 XLSX with stdlib; exclude its autopsy section.

    Retain the original field text, including the source's ambiguous PMI label.
    No spreadsheet engine, Excel application or manually curated table needed.
    """
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        strings = [''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si', ns)]
        sheet = ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
        rows = []
        for row in sheet.findall('s:sheetData/s:row', ns):
            values = {}
            for cell in row.findall('s:c', ns):
                v = cell.find('s:v', ns)
                if v is not None:
                    values[re.sub(r'\d', '', cell.attrib['r'])] = strings[int(v.text)] if cell.get('t') == 's' else v.text
            rows.append(values)
    header = rows[0]
    if header.get('A') != 'Library' or header.get('B') != 'Patient ID#':
        raise ValueError('Unexpected Table S1 layout')
    records = []
    for row in rows[1:]:
        if re.fullmatch(r'[AP]\d+', row.get('A', '')):
            records.append({header[k].strip(): v for k, v in row.items() if k in header})
    clinical = pd.DataFrame(records).set_index('Library')
    if not clinical.index.is_unique or set(clinical.index) != set(samples.index):
        raise ValueError('Clinical libraries do not match GEO')
    for sid, row in clinical.iterrows():
        geo = samples.loc[sid]
        for c, g in [('Patient ID#','donor_id'), ('Brain Region','brain_region'), ('SEX','geo_sex'),
                     ('RACE','geo_race'), ('Hemisphere','geo_hemisphere'), ('Library Batch','geo_batch')]:
            if row[c].strip().lower() != str(geo[g]).strip().lower():
                raise ValueError(f'Table S1/GEO conflict: {sid}, {c}')
        for c, g in [('AGE','geo_age_yr'), ('Epilepsy Duration','geo_epilepsy_duration_yr'), ('RIN','geo_rin')]:
            if float(row[c].split()[0]) != float(geo[g]):
                raise ValueError(f'Table S1/GEO numeric conflict: {sid}, {c}')
    clinical.columns = ['table_s1_'+re.sub(r'\W+', '_', c.lower()).strip('_') for c in clinical.columns]
    result = samples.join(clinical, validate='one_to_one')
    result['medications'] = result.table_s1_medications
    result['seizure_frequency_per_month'] = result.table_s1_seizure_frequency.str.split().str[0].astype(float)
    result['invasive_mapping'] = result.table_s1_invasive_mapping
    result['pathology'] = result.table_s1_neuropath_report_summary
    result['library_chemistry'] = result.table_s1_library_protocol
    # Table S1 says LEFT in its Hemisphere column but RIGHT in Donor4's report.
    # Preserve both and flag the contradiction rather than resolving it by guess.
    result['hemisphere_report_conflict'] = [
        ('RIGHT HIPPOCAMPUS' in r['pathology'].upper() and r['geo_hemisphere'] == 'left') or
        ('LEFT HIPPOCAMPUS' in r['pathology'].upper() and r['geo_hemisphere'] == 'right')
        for _, r in result.iterrows()]
    return result


def read_geo(path):
    records = []
    with gzip.open(path, 'rt') as f:
        text = f.read()
    for block in text.split('^SAMPLE = ')[1:]:
        row = {'source_sample_id': block.splitlines()[0]}
        for line in block.splitlines():
            if line.startswith('!Sample_characteristics_ch1 = '):
                key, val = line.split(' = ', 1)[1].split(': ', 1)
                row['geo_'+re.sub(r'\W+', '_', key.lower()).strip('_')] = val
            elif line.startswith('!Sample_title = '):
                row['geo_title'] = line.split(' = ', 1)[1]
            elif line.startswith('!Sample_source_name_ch1 = '):
                row['donor_id'] = line.split(' = ', 1)[1]
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
    if len(frame) != 10 or not frame.index.is_unique or frame.donor_id.nunique() != 5:
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
    inventory = json.loads((meta/'download_inventory.json').read_text())
    required = {MATRIX, CLINICAL, 'GSE160189_family.soft.gz', 'GSE160189_series_matrix.txt.gz'}
    if not required.issubset({r['filename'] for r in inventory}):
        raise ValueError('Download inventory is incomplete; run scripts/download/ayhan.py')
    for r in inventory:
        if sha256(raw/r['filename']) != r['sha256']:
            raise ValueError(f'Raw checksum mismatch: {r["filename"]}')
    samples = read_geo(raw/'GSE160189_family.soft.gz')
    samples = read_clinical(raw/CLINICAL, samples)
    x, cells, genes = read_counts(raw/MATRIX)
    prefixes = cells.str.split('_').str[0]
    if set(prefixes) != set(samples.index):
        raise ValueError('Cell library prefixes do not match GEO samples exactly')
    obs = samples.loc[prefixes].copy()
    obs.index = cells
    obs = standardize_obs(obs, 'ayhan')
    obs['diagnosis'] = 'temporal lobe epilepsy'
    obs['control_status'] = 'epilepsy'
    obs['assay'] = 'snRNA-seq'
    obs['platform'] = '10x Genomics Chromium; Illumina NovaSeq 6000'
    obs['source_accession'] = 'GSE160189'
    obs['source_file'] = 'ayhan/'+MATRIX
    obs['author_annotation_available'] = False
    obs['tissue_source'] = 'surgical resection'
    for src, dst in [('geo_age_yr','age_years'), ('geo_epilepsy_duration_yr','epilepsy_duration_years'), ('geo_rin','rin')]:
        obs[dst] = pd.to_numeric(obs[src], errors='raise')
    obs['sex'] = obs.geo_sex
    obs['hemisphere'] = obs.geo_hemisphere
    obs.index = pd.Index('ayhan:'+cells, name='cell_id')
    for col in obs.select_dtypes('object'):
        obs[col] = pd.Categorical(obs[col])
    var = pd.DataFrame({'source_gene_id': genes, 'gene_symbol': genes}, index=pd.Index(genes, name='gene_name'))
    result = ad.AnnData(x, obs=obs, var=var)
    result.layers['counts'] = result.X
    result.uns['ingestion'] = dict(raw_counts_available=True, expression_type='raw_counts',
        source_files=[f'ayhan/{r["filename"]}' for r in inventory],
        metadata_sources=['GSE160189_family.soft.gz', 'GSE160189_series_matrix.txt.gz', CLINICAL],
        manifest_study_id='Ayhan2021', genome_build='hg19',
        matrix_provenance='GEO: Cell Ranger 3.0.2 count matrices; deposited gene-by-cell CSV',
        transformations='Transpose and lossless sparse int32 storage only; all deposited cells and genes retained',
        cell_mapping='Exact CSV cell prefix to GEO sample title; paired donors from GEO source_name',
        missing_metadata='No GEO cell-type annotation; FCD subtype and mutation status unavailable',
        clinical_provenance='Table S1 surgical library rows only; source field text retained with table_s1 prefix; matched and cross-checked against GEO',
        metadata_caveats='Donor4 report says right while Hemisphere field and GEO say left; flagged. Table S1 A57 cell count is three larger than deposited CSV. PMI source label retained without interpreting as postmortem interval for surgical tissue.',
        cohort_caveat='Five surgical epilepsy donors, paired anterior/posterior hippocampus; no healthy-control libraries')
    digest = fingerprint(x)
    samples['nuclei'] = pd.Series(prefixes.value_counts())
    samples['table_s1_cell_count_difference'] = samples.nuclei - pd.to_numeric(samples.table_s1_cell_number)
    samples.to_csv(meta/'sample_metadata.csv', index=False)
    donor_cols = ['donor_id','geo_age_yr','geo_sex','geo_race','geo_hemisphere','geo_epilepsy_duration_yr',
                  'medications','seizure_frequency_per_month','invasive_mapping','pathology','hemisphere_report_conflict']
    samples[donor_cols].drop_duplicates().to_csv(meta/'clinical_metadata.csv', index=False)
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
                   duplicate_cell_ids=0, duplicate_gene_ids=0, unmatched_cells=0, unmatched_samples=0,
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
        'Metadata is joined by explicit library prefix and cross-checked with GEO and Table S1.\n\n'
        'GEO describes Cell Ranger 3.0.2 counts. Its format description incorrectly says '
        'tab-delimited cells-by-genes: the file is comma-delimited genes-by-cells. '
        '131,325 nuclei agrees with GEO; the paper reports 129,908 after author QC '
        '(1,417 fewer). No author QC was reapplied. Table S1 A57 exceeds the CSV by '
        'three cells. Donor4 hemisphere conflicts with its pathology narrative; both '
        'source values and a conflict flag are retained. No GEO cell-type annotations, '
        'FCD subtype, mutation status or healthy-control libraries are supplied. '
        'The paper links an external cell browser; its annotations are not part of this GEO deposit.\n\n'
        'See validation_summary.json for full count fingerprint and sample totals, '
        '../../metadata/ayhan/download_inventory.json for file sizes and SHA256, '
        'and ../../../scripts/AYHAN_README.md for source references and reproduction commands.\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=ROOT)
    run(p.parse_args().project_root.resolve())
