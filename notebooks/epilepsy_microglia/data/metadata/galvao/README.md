# Galvao source and ingestion provenance

Source of truth: the Galvao2024 row of `/home/neuro_demo_research/dataset_manifest.csv`,
recorded verbatim as field values in manifest_entry.json. Only Galvao was processed.

- download_sources.json: exact original filenames, URLs and CELLxGENE expected size.
- download_inventory.json: byte sizes, SHA256 and completed integrity checks.
- source_inspection.json: actual matrix dimensions, feature types and barcode checks,
  recorded before ingestion implementation.
- sample_mapping.csv and ingestion_audit.json: explicit file/GEO/CELLxGENE mapping,
  source nuclei and feature counts, shared-count comparisons and exact per-specimen
  RNA fingerprints checked again after reopening the output.
- cellxgene_obs_inspection.csv: unmodified observation metadata extracted for inspection.
- ingestion.log: ingestion run and final validation summary.

Run from the project root with the existing scenv environment:

```bash
python scripts/download/galvao.py
/home/neuro_demo_research/miniconda3/envs/scenv/bin/python scripts/ingest/galvao.py
PYTHONPATH=src /home/neuro_demo_research/miniconda3/envs/scenv/bin/python scripts/validate/validate_h5ad.py data/processed/galvao/galvao.h5ad --study galvao
```

Downloading and ingestion are separate. Downloads resume in metadata staging,
then publish completed originals without overwrite. Ingestion verifies raw SHA256
before and after and refuses an existing output. Re-running the validator is read-only.
The final dataset report and full obs.head() export are in data/processed/galvao/.

Metadata comes directly from GEO characteristics. Normal histology is recorded as
internal histologically normal tissue from an FCD donor, never as healthy control.
FCD subtype is missing for these specimens; mutation status is missing throughout.
CELLxGENE columns carry a cellxgene_ prefix to distinguish author/curator annotations
from standardized GEO metadata. Existing author QC metrics are retained without
recomputing them or using them to select nuclei. Full author var statistics and
embeddings remain in the immutable source H5AD and are not used for ingestion.
