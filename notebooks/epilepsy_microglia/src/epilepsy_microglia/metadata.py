"""Explicit metadata alignment without clinical inference."""
import pandas as pd

OBS_COLUMNS = (
    "study", "donor_id", "sample_id", "diagnosis", "pathology", "fcd_subtype",
    "control_status", "mutation_status", "brain_region", "assay", "platform",
    "source_accession", "source_sample_id", "source_cell_id", "source_file",
)


def standardize_obs(obs, study):
    """Add missing columns, retaining existing values and original columns.

    Call before changing the observation index. Source cell IDs are copied
    verbatim; sample/donor IDs and clinical values are never guessed.
    """
    result = obs.copy()
    if "study" in result and (result["study"].dropna() != study).any():
        raise ValueError("Existing study values conflict with requested study")
    if "source_cell_id" not in result:
        result["source_cell_id"] = result.index.copy()
    for column in OBS_COLUMNS:
        if column not in result:
            result[column] = pd.Series(pd.Categorical([None] * len(result)), index=result.index)
    result["study"] = study
    return result


def join_metadata(obs, metadata, *, key):
    """Left join on an explicit source key; preserve order and unmatched NA.

    Rename source fields explicitly before joining. Reject null/duplicate
    lookup keys and overlapping columns to avoid accidental overwrites.
    """
    if key not in obs or key not in metadata:
        raise ValueError(f"Missing matching key: {key}")
    if metadata[key].isna().any() or metadata[key].duplicated().any():
        raise ValueError("Metadata matching keys must be non-null and unique")
    overlap = (set(obs) & set(metadata)) - {key}
    if overlap:
        raise ValueError(f"Overlapping metadata columns: {sorted(overlap)}")
    result = obs.join(metadata.set_index(key), on=key, validate="many_to_one")
    if not result.index.equals(obs.index):
        raise ValueError("Metadata join changed observation order")
    return result
