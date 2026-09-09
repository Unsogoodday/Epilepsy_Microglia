# Ayhan 2021 / GSE160189

Run from `notebooks/epilepsy_microglia` in a copy of this repository. Python 3.10
or 3.11 and the `curl` executable are required. No Codex, R, Excel, credentials,
or manually prepared metadata files are needed. The downloader uses the Python
standard library; ingestion uses the local `src/epilepsy_microglia` helpers.
Clone the repository and enter `notebooks/epilepsy_microglia`.
Allow approximately 12 GB RAM and 4 GB free disk space for one run (more for backups).

```bash
git clone git@github.com:Unsogoodday/Epilepsy_Microglia.git
cd Epilepsy_Microglia/notebooks/epilepsy_microglia
python3.10 -m venv .venv-ayhan
source .venv-ayhan/bin/activate
python -m pip install -r scripts/ayhan-requirements.txt
python scripts/download/ayhan.py
python scripts/ingest/ayhan.py
python tests/test_ayhan_ingest.py
```

On Ubuntu 22.04, system prerequisites can be installed with
`sudo apt-get install git curl python3.10 python3.10-venv`.

The downloader reads `config/dataset_manifest.csv` (the supplied manifest,
bundled for reproducibility). Override it with `--manifest /path/to/dataset_manifest.csv`.
It derives GEO URLs from that entry and discovers expression filenames in the
official supplementary listing. This format-specific reader rejects a changed accession.
The manifest does not enumerate clinical attachments: the inspected paper's
Table S1 publisher URL is an explicit supplemental source in the downloader.
There is no top-level `ingest.py` in this repository: Galvao's ingestion entry
point is `scripts/ingest/galvao.py`, and Ayhan follows it with `scripts/ingest/ayhan.py`.
Both scripts accept `--project-root /path/to/output-project`; use the same root
for both. This controls data locations independently of the working directory.
The corresponding Galvao commands are `python scripts/download/galvao.py` and
`python scripts/ingest/galvao.py`.

Downloads resume in `data/metadata/ayhan/download_staging`, validate full gzip
streams or Excel ZIP contents, and publish without replacing existing originals.
Rerunning the downloader verifies existing files and compares recorded SHA256
checksums. Invalid completed files are moved into staging with `.invalid-<timestamp>`
suffixes before replacement; no download directory is deleted. Invalid new downloads
are preserved in staging and fail before conversion; rerun to retry from zero.
Ingestion is offline. A rerun reconstructs the full matrix and metadata and compares
them against the existing h5ad, without replacing it. An interrupted write can leave
an unreferenced temporary h5ad, which is ignored on rerun. To rebuild a differing
output, archive the previous processed file first. No output is silently replaced.

## Sources and outputs

- [GEO GSE160189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160189):
  every file in its supplementary directory (currently one count CSV), family
  SOFT, and series matrix. Originals remain in `data/raw/ayhan/`.
- [Clinical Table S1](https://ars.els-cdn.com/content/image/1-s2.0-S0896627321003299-mmc2.xlsx)
  from [Ayhan et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8273123/): downloaded
  from the publisher; the PMC attachment currently returns a browser challenge.
- `data/metadata/ayhan/clinical_metadata.csv`: five donor records, demographics,
  epilepsy duration, medications, seizure frequency, invasive mapping,
  neuropathology report, and hemisphere conflict flag.
- `data/metadata/ayhan/sample_metadata.csv`: ten libraries, GEO characteristics,
  all Table S1 fields, RIN, chemistry, sequencing statistics, observed nuclei
  counts, and differences from Table S1.
- `data/metadata/ayhan/download_inventory.json`: original URLs, sizes, SHA256.
- A completed ingestion writes `data/processed/ayhan/ayhan.h5ad`. For the
  inspected source it contains 131,325 nuclei × 17,180 genes, with integer CSR
  counts in both `X` and `layers['counts']`. Original cell/gene identifiers and
  metadata provenance are retained. No new QC filtering, normalization,
  integration, or cell-type prediction occurs.
- `data/processed/ayhan/validation_summary.json`: counts, library totals, and
  full sparse matrix fingerprint verified after reopening the h5ad.

## Source-specific interpretation

The actual CSV is **genes × nuclei**, despite GEO's generic description stating
the opposite. Each original cell ID starts with its library (e.g. `P57_`), which
is joined exactly to a GEO title such as `P57_Donor2_scRNA-seq`. Donors are checked
to have paired anterior/posterior libraries. Table S1's ten surgical library
rows are joined and checked against GEO for donor, region, age, sex, race,
hemisphere, epilepsy duration, RIN, and batch. Its separate autopsy rows are not
snRNA libraries and never enter metadata or AnnData observations.

All five donors have temporal lobe epilepsy; none is a healthy control. Raw
counts here means deposited Cell Ranger-derived integer counts, not an
unfiltered droplet matrix. Retain all 131,325 deposited nuclei even though the
paper describes 129,908 nuclei after author QC (1,417 fewer); no author QC is
reapplied here. GEO provides no cell-level
author annotations; `author_annotation_available=False` specifically refers to
this GEO-based pipeline. FCD subtype and mutation status remain missing.

Two source discrepancies are retained explicitly:

1. Table S1 A57 reports 15,356 cells; the CSV has 15,353. The sample table records
   a difference of -3; no cells are fabricated or discarded.
2. Donor4's neuropathology report says RIGHT, while GEO and Table S1's Hemisphere
   column say left. `hemisphere` follows those structured source fields, the
   verbatim report is retained in `pathology`, and `hemisphere_report_conflict`
   is true for A76/P76.

Table S1 labels a surgical-tissue timing field `PMI` (0.2 hours). Preserve it as
`table_s1_pmi`; do not interpret it as evidence these were postmortem samples.
Original clinical text is retained in `table_s1_*` fields; derived numeric age,
duration, and monthly seizure frequency have explicit units in column names.

The deposited matrix is 171,687,706 compressed bytes, consistent with the GEO
listing's rounded 164M. GEO supplies no checksum in this listing: complete gzip
CRC/length checks and locally recorded SHA256 establish integrity and repeatability,
not independent publisher checksum authentication. Excel is checked using ZIP CRCs.
GEO incorrectly describes tab-delimited cells-by-genes input; the actual file is
comma-delimited genes-by-cells, which the parser verifies explicitly.

The pipeline does not claim that the output exists until ingestion completes and
the written H5AD passes its reopen, schema, axis, and full sparse-count checks.
