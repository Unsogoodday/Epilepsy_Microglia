# Download stage (pending curated sources)

No download entry point is enabled yet. Populate config/datasets.yaml from the
next curated source list; do not execute URLs archived in legacy/.

Each source will have its exact published filename, URL, required/optional flag
and supplied checksum if available. Use curl with --fail --location --retry and
--continue-at -, downloading into a separate staging directory. Resume only
partial staging files; never append to or overwrite completed files in raw/.
Verify gzip streams, tar members and ZIP CRCs before copying completed originals
with exclusive creation into data/raw/<study>/. Keep original archive names and
bytes. Record sizes, SHA256, source URLs and archive verification results in
metadata. Archives stay compressed in raw/; extraction belongs in processed
scratch space after inspection, with path traversal and link checks.

Inspect matrix dimensions, feature types, barcode keys, metadata tables and
count provenance. Write data/processed/<study>/dataset_report.md only after
inspection. Stop for review before conversion when requested. Never call an
ingestion script automatically from a downloader.
