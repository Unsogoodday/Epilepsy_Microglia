# Preserved reference only — do not run automatically

original/ contains the previous project, moved byte-for-byte including notebooks,
metadata and local caches. original_sha256.json records its pre-move file hashes.
Git preservation commit: 65e132f. Other pre-existing staged edits/deletions in the
parent repository were left intact.

Old flow: environment-dependent external paths → GEO download/extraction →
filename rewriting and CSV/10x conversion → heuristic orientation → concatenation
→ MyGene annotation → QC filtering, normalization/log1p, HVGs → scVI/scANVI and
labeling. The new package deliberately imports none of this code.

Useful leads: 10x MTX triplets may be flat and sample-prefixed; text count tables
need explicit axis checks; gene IDs and symbols are distinct; sample identities
must survive concatenation. The old Kumar table lists GSM6049632–GSM6049643,
replicate labels and clinical fields but lacks per-field evidence. Kumar01A/B
must not be assumed to be different donors or merged without verification.

Unverified old accession leads: Kumar GSE201048, Thrupp GSE153807, Ayhan GSE160189,
Pai GSE140393. The parent datasets.txt also mentions GSE157277 for Pappalardo.
The old Pai notebook flags bulk sequencing and the parent list flags Pappalardo
as FASTQ: verify assay and available count matrices from curated sources.
Galvao and Liu have no study-specific implementation in the inspected project.

Dangerous old behaviors included deleting compressed originals, renaming files
without recompression, guessing orientation from uppercase names, and implicit
concatenation defaults. Preserve these only as historical code, not ingestion.
