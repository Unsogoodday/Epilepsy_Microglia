# Epilepsy microglia ingestion workspace

Ayhan's reproducible download and ingestion commands, source interpretation,
and validation details are in [scripts/AYHAN_README.md](scripts/AYHAN_README.md).
Galvao and Ayhan are implemented; the pending-source statements below describe
the original restructuring and still apply to the other studies.

Standardize the output, not the input. This project lives at
`notebooks/epilepsy_microglia/` in the parent git repository. Run commands here.
No datasets have been downloaded or converted by this restructuring.

```text
config/datasets.yaml          seven-study registry and pipeline paths
config/dataset_manifest.csv   curated source instructions
data/raw/<study>/            immutable original downloads (created on download)
data/processed/<study>/      reports, scratch, then <study>.h5ad after inspection
data/metadata/               source-grounded mapping tables and provenance
scripts/download/            Galvao and Ayhan downloaders plus stage policy
scripts/ingest/<study>.py     explicit per-study ingestion entry points
scripts/validate/            shared validation CLI
src/epilepsy_microglia/       io.py, metadata.py, validation.py
notebooks/                   future exploratory notebooks
legacy/original/             preserved old project; reference only
tests/                       synthetic ingestion-contract checks
```

The studies are Galvao, Kumar, Ayhan, Pai, Thrupp, Pappalardo and Liu. Galvao
and Ayhan have inspected, executable download and ingestion workflows. Pending
entry points intentionally fail without writing anything until their source
formats have been inspected. There is no universal loader.

## Output contract

Each study writes `data/processed/<study>/<study>.h5ad`, cells/nuclei in rows,
features in columns. Preserve all source observations/features and expression;
no QC filtering, normalization, log transform, batch correction
or biological analysis. Document explicit axis decisions and within-study
sample assembly; never silently drop genes using an implicit inner join.

`obs` includes study, donor_id, sample_id, diagnosis, pathology, fcd_subtype,
control_status, mutation_status, brain_region, assay, platform, source_accession,
source_sample_id, source_cell_id and source_file. Unknown values remain actual
missing values, never inferred labels or strings such as "Unknown". Retain
original metadata columns. Only study, source_cell_id and source_file must be
populated; missing sample/donor information is documented in the report.
Use explicit source keys for joins before adding empty standardized columns.
Keep unique obs/var indices; if IDs repeat between samples, namespace explicitly
and preserve verbatim IDs in obs.source_cell_id and var.source_gene_id. Preserve
original feature names/symbols and feature types in additional var columns.

Use gene symbols as var_names. The shared writer converts Ensembl IDs using
source-grounded var.gene_symbol mappings and preserves var.source_gene_id.
Supply those mappings during each study ingestion; missing mappings fail rather
than inventing symbols. Duplicate symbols receive numeric suffixes (-1, -2,
etc.) without merging or dropping features. Existing symbol indices stay unchanged.

When raw counts are documented as available, retain them unchanged in X and
layers['counts']. Otherwise preserve supplied expression in X, omit counts,
and document why counts are unavailable. Do not reverse normalization or assume
AnnData.raw means raw counts. For multimodal sources explicitly identify RNA
features and preserve other modalities separately with documented source links.

Record uns['ingestion'] with raw_counts_available (bool), expression_type
('raw_counts' or an explicit source scale), source_files (relative raw paths),
and metadata_sources (source references, empty list when absent). Reports should
record filenames, sizes, hashes/integrity, assay evidence, count availability,
metadata keys/unmatched records, conversion plan, missing files and ambiguities.
Validation checks structure and numerical consistency; it cannot prove clinical
truth or count provenance. write_study validates and refuses to overwrite output.

## Local checks

Use a Python environment containing the dependencies in pyproject.toml. No old
GPU/analysis environment is required. From this directory:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/validate/validate_h5ad.py data/processed/kumar/kumar.h5ad --study kumar
python scripts/ingest/kumar.py --help
python scripts/download/ayhan.py --help
python scripts/ingest/ayhan.py --help
python tests/test_ayhan_ingest.py
```

The parent requirements/environment files remain historical and unchanged.
See legacy/README.md for recovered knowledge and unresolved source ambiguities.

The contract tests cover dense/sparse H5AD round trips, missing clinical
metadata, duplicate-key rejection, invalid counts, existing-output protection
and raw-directory symlink protection. The Ayhan tests add strict CSV parsing,
axis preservation, and resumable-download integrity checks.
