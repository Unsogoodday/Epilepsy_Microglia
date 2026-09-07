"""Read-only structural/count checks, not biological QC."""
import numpy as np
from scipy import sparse
from .metadata import OBS_COLUMNS
from .genes import ensembl_mask

STUDIES = {"galvao", "kumar", "ayhan", "pai", "thrupp", "pappalardo", "liu"}


def validate_adata(adata, *, study):
    """Return errors. Missing clinical values are valid; missing columns are not.

    Counts must be declared based on source documentation, not this numeric
    test. Integer-valued expression alone does not establish raw provenance.
    """
    errors = []
    if study not in STUDIES:
        errors.append("Unconfigured study")
    if not adata.n_obs or not adata.n_vars:
        errors.append("Expression matrix is empty")
    for axis, frame in (("obs", adata.obs), ("var", adata.var)):
        if not frame.index.is_unique or frame.index.isna().any():
            errors.append(f"{axis} index must be unique and non-null")
        if not frame.columns.is_unique:
            errors.append(f"{axis} columns must be unique")
    if ensembl_mask(adata.var_names).any():
        errors.append("var_names must use gene symbols, not Ensembl IDs")
    for column in OBS_COLUMNS:
        if column not in adata.obs:
            errors.append(f"Missing obs column: {column}")
    if "study" in adata.obs:
        if adata.obs.study.isna().any() or not adata.obs.study.eq(study).all():
            errors.append("obs.study does not match study")
    for column in ("source_cell_id", "source_file"):
        if column in adata.obs and adata.obs[column].isna().any():
            errors.append(f"obs.{column} must be populated")
    if "source_gene_id" not in adata.var or adata.var["source_gene_id"].isna().any():
        errors.append("var.source_gene_id must preserve original feature identifiers")
    provenance = adata.uns.get("ingestion", {})
    if not isinstance(provenance, dict):
        return errors + ["uns.ingestion must be a mapping"]
    for key in ("source_files", "metadata_sources", "expression_type"):
        if key not in provenance:
            errors.append(f"Missing ingestion provenance: {key}")
    if not len(provenance.get("source_files", [])):
        errors.append("Record at least one source file")
    if not isinstance(provenance.get("expression_type"), str) or not provenance.get("expression_type", "").strip():
        errors.append("expression_type must describe the source scale")
    available = provenance.get("raw_counts_available")
    if not isinstance(available, (bool, np.bool_)):
        errors.append("raw_counts_available must be explicitly true or false")
    elif available:
        if provenance.get("expression_type") != "raw_counts":
            errors.append("Available raw counts must be retained in X")
        if "counts" not in adata.layers:
            errors.append("Raw counts require layers['counts']")
        else:
            counts = adata.layers["counts"]
            values = counts.data if sparse.issparse(counts) else np.asarray(counts)
            if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
                errors.append("Counts must be finite, nonnegative and integer-valued")
            if adata.X is None or (adata.X != counts).sum() != 0:
                errors.append("X and counts layer must match at ingestion")
    else:
        if "counts" in adata.layers or provenance.get("expression_type") == "raw_counts":
            errors.append("Do not label expression as counts when raw counts are unavailable")
    if adata.X is None:
        errors.append("X must contain source expression")
    elif not np.isfinite(adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)).all():
        errors.append("X contains non-finite values")
    return errors
