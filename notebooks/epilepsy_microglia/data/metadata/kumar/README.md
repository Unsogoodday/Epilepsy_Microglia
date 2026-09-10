# Kumar2022 / GSE201048 — MTX and tabular inputs only

Run from the `notebooks/epilepsy_microglia` project directory:

```bash
python -m pip install -e . -r scripts/kumar-requirements.txt
python scripts/download/kumar.py
python scripts/ingest/kumar.py
python scripts/validate/validate_h5ad.py data/processed/kumar/kumar.h5ad --study kumar
python -m pytest tests -q
```

Python >=3.10 and curl are sufficient; no R, R-related package, Codex, or Drive
connection is required. The downloader uses only the 33 explicitly listed MTX/TSV
URLs in `source_triplets.csv`. This user-requested scope supersedes the manifest's
older request for all supplementary files. No R object is downloaded, opened,
parsed, hashed, or required by either script. Previously completed raw files from
the earlier scope remain untouched and excluded from the inventory and ingestion.

Downloads resume in `download_staging`, validate complete gzip CRC/length streams,
check official byte sizes, and record SHA256 in `download_inventory.json` after
each success. Reruns verify completed files against that inventory without
overwriting them. Invalid partial files are preserved with a timestamp, allowing
a subsequent retry. Ingestion checks exact inventory coverage and every checksum
before reading the triplets; it performs no network access.

## Source metadata and conflict decisions

`source_samples.csv` records the official GEO accession/title/specimen mapping;
`source_triplets.csv` records the official sample supplementary URLs and byte sizes.
The mapping and sizes were checked against GEO before the format restriction;
only CSV records are inputs to the final workflow. Evidence:
[GSE201048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201048),
[GEO file list](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201048/suppl/filelist.txt).

The three `curated_*.csv` files are verbatim Kumar-only snapshots of the study,
sample, and tissue tabs in the
[curated workbook](https://docs.google.com/spreadsheets/d/1lD4crm_CdD3YSu-x_oGUvS52SpUbusD5VZIhrnJRSiI/edit),
retrieved 2026-09-10. Their original row IDs and values are retained for provenance.
`clinical_donors.csv` contains reported onset/procedure values and conflict notes
checked against official Supplementary Table S5, pages 51–52 of the
[clinical supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41593-022-01095-5/MediaObjects/41593_2022_1095_MOESM1_ESM.pdf).
No PDF or online metadata service is needed for reruns.

- GEO lists GSM6049632–GSM6049642, exactly 11 samples. The extra legacy
  GSM6049643 is excluded because it is absent from this series, not because of QC.
- The 11 sample IDs P1.A through P6.C map to six donor IDs P1–P6. The join checks
  GEO accession and title against curated assay IDs, then tissue and donor keys.
  P2 and P3 are FCD IIb (Table S5); all donors have epilepsy and none is healthy.
- P3.A is a temporal specimen in GEO and Figure 1c, whereas Table S5 describes
  a left frontal procedure. Preserve specimen region and procedure independently.
- P6 is reported as age 4y with seizure onset 4y8m. Both are retained as reported;
  epilepsy duration stays missing. No age correction is inferred.
- Figure 1c specifies occipital cortex/core for P1; the paper's OL/olfactory
  abbreviation is inconsistent. Use the specific GEO/Figure 1c regions.
- Unknown mutation status remains missing. Original curated genetic categories
  are retained in their own columns and do not establish measured mutations.
- The curated study title is a descriptive label; the official publication is
  *Single-cell transcriptomics and surface epitope detection in human brain
  epileptic lesions identifies pro-inflammatory signaling*,
  [DOI 10.1038/s41593-022-01095-5](https://doi.org/10.1038/s41593-022-01095-5).

## Output contract

RNA output: `data/processed/kumar/kumar.h5ad`. All deposited cells and RNA features
are retained; `X` and `layers['counts']` contain identical sparse int32 raw UMI
counts. Source gene IDs, unmodified symbols, barcodes, accessions and filenames
are preserved. Shared `write_study` applies unique gene-symbol indices without
merging features. No QC, normalization, or paper-total-based filtering occurs.

Recoverable Antibody Capture counts are stored separately in
`data/processed/kumar/cite_seq/<GSM>.h5ad`, retaining the 16- versus 18-feature
panels per sample without imputing unmeasured antibodies. Each file is reopened
and checked. No CSV/TSV cell annotation file was supplied, so cell-level author
annotation availability is false; clusters are not inferred.

RNA is reopened and checked against every sample's count fingerprint, all
observation/feature metadata, sample/donor alignment, and the shared validator.
All input SHA256 values are rechecked. `ingestion_audit.json` and
`validation_summary.json` record execution evidence; `sample_mapping.csv` records
the validated join. Existing RNA output causes an explicit refusal to overwrite;
use a separate project root for a fresh ingestion, or run the validator on it.
