# Ayhan 2021 / GSE160189

Run from `notebooks/epilepsy_microglia` in a copy of this repository. Python 3.10
or 3.11 and the `curl` executable are required. No Codex, R, Excel, credentials,
or manually prepared metadata files are needed. The downloader uses the Python
standard library; ingestion uses the local `src/epilepsy_microglia` helpers.
Clone the repository and enter `notebooks/epilepsy_microglia`.
Allow approximately 12 GB RAM and 4 GB free disk space for one run (more for backups).

After the environment is installed, the complete pipeline is exactly these two
commands when run from this `scripts` directory:

```bash
python -u ./download/ayhan.py
python -u ./ingest/ayhan.py
```

The first command visibly resumes partial downloads, retries transient network
failures, and verifies the source files. The second prints the absolute path to
the completed `ayhan.h5ad`.

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
There is no top-level `ingest.py` in this repository: Galvao's ingestion entry
point is `scripts/ingest/galvao.py`, and Ayhan follows it with `scripts/ingest/ayhan.py`.
Both scripts accept `--project-root /path/to/output-project`; use the same root
for both. This controls data locations independently of the working directory.
The corresponding Galvao commands are `python scripts/download/galvao.py` and
`python scripts/ingest/galvao.py`.

Downloads resume in `data/metadata/ayhan/download_staging`, validate full gzip
streams, and publish without replacing existing originals.
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
- `data/metadata/ayhan/sample_metadata.csv`: ten libraries, GEO characteristics,
  and observed nuclei counts.
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
to have paired anterior/posterior libraries. The GEO sample records provide and
are checked for donor, region, age, sex, race, hemisphere, epilepsy duration,
RIN, and batch.

All five donors have temporal lobe epilepsy; none is a healthy control. Raw
counts here means deposited Cell Ranger-derived integer counts, not an
unfiltered droplet matrix. Retain all 131,325 deposited nuclei even though the
paper describes 129,908 nuclei after author QC (1,417 fewer); no author QC is
reapplied here. GEO provides no cell-level
author annotations; `author_annotation_available=False` specifically refers to
this GEO-based pipeline. FCD subtype and mutation status remain missing.

The deposited matrix is 171,687,706 compressed bytes, consistent with the GEO
listing's rounded 164M. GEO supplies no checksum in this listing: complete gzip
CRC/length checks and locally recorded SHA256 establish integrity and repeatability,
not independent publisher checksum authentication.
GEO incorrectly describes tab-delimited cells-by-genes input; the actual file is
comma-delimited genes-by-cells, which the parser verifies explicitly.

The pipeline does not claim that the output exists until ingestion completes and
the written H5AD passes its reopen, schema, axis, and full sparse-count checks.
