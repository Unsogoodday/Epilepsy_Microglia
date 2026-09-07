"""Validated output only. No source parsing or network access."""
from pathlib import Path
from tempfile import NamedTemporaryFile
import os
from .validation import STUDIES, validate_adata
from .genes import standardize_var_names


def write_study(adata, *, study, project_root):
    """Write the standard output without overwriting an existing artifact.

    Invoke only after source inspection and study-specific ingestion. Raw files
    are read by the study script and never rewritten by this helper.
    """
    if study not in STUDIES:
        raise ValueError("Unknown study")
    standardize_var_names(adata)
    errors = validate_adata(adata, study=study)
    if errors:
        raise ValueError("\n".join(errors))
    root = Path(project_root).resolve()
    destination = root / "data" / "processed" / study / f"{study}.h5ad"
    if destination.resolve() != destination:
        raise ValueError("Output path must not traverse symlinks")
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=destination.parent, suffix=".h5ad", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        adata.write_h5ad(temporary)
        # Hard-link publication is atomic and refuses existing destinations.
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination
