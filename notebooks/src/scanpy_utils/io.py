from pathlib import Path
from collections import defaultdict
import scanpy as sc
import pandas as pd
import shutil, subprocess, gzip

"""
    Codes related to dataset download
    - _download_dataset


"""

from src.files_utils import _download_tar_from_link, _download_from_link
def _download_dataset(
    *,
    link: str,
    download_dir: Path,
    is_tar: bool,
) -> None:
    if is_tar:
        _download_tar_from_link(link, download_dir)
    else:
        _download_from_link(link, download_dir)

"""
    Codes related to dataset download
    - _build_mtx_with_auto_wrap (final wrapper)
    - _wrap_flat_10x_to_subdirs
    - _build_mtx_anndata (intermediate wrapper)
    - _looks_like_10x_dir
    - _standardize_10x_filenames


    
"""

def _build_mtx_with_auto_wrap(
    download_dir: Path,
    raw_dir: Path,
    label: str,
) -> list[Path]:
    """
    Build anndata from mtx regardless of directory structures
    If the directory is flat, create subdirs first
    Then, build anndata from mtx

    """

    # identify subdirectories
    has_subdirs = any(p.is_dir() for p in download_dir.iterdir())

    if not has_subdirs:
        _wrap_flat_10x_to_subdirs(download_dir)

    return _build_mtx_anndata(
        download_dir=download_dir,
        raw_dir=raw_dir,
        label=label,
    )


def _wrap_flat_10x_to_subdirs(
    download_dir: Path,
    sep: str = "_",
) -> Path:
    """
    Convert flat 10x files into per-sample subdirectories.

    Expected flat structure:
        SAMPLEID_matrix.mtx
        SAMPLEID_features.tsv / genes.tsv
        SAMPLEID_barcodes.tsv

    After wrapping:
        download_dir/
            SAMPLEID/
                matrix.mtx
                features.tsv
                barcodes.tsv
    """
    if not download_dir.exists():
        raise FileNotFoundError(f"{download_dir} does not exist")

    files = [p for p in download_dir.iterdir() if p.is_file()]

    buckets: dict[str, list[Path]] = defaultdict(list)

    for f in files:
        if sep not in f.stem:
            continue
        sample_id = f.stem.split(sep, 1)[0]
        buckets[sample_id].append(f)

    for sample_id, paths in buckets.items():
        sample_dir = download_dir / sample_id
        sample_dir.mkdir(exist_ok=True)

        for p in paths:
            target_name = p.name.split(sep, 1)[1]
            shutil.move(p, sample_dir / target_name)
        print(f"Wrapped {p.stem}")

    return download_dir

def _build_mtx_anndata(
    download_dir: Path,
    raw_dir: Path,
    label: str, 
) -> list[Path]:
    """
        Read matrix.mtx / features.tsv / genes.tsv -> .h5ad
        Expects the directory is already classified into subdirs.
        Normalize subdir file names 
        Convert multiple files into AnnData objects
        and save them to raw_dir

        Returns
        ------
        list[Path]
            Paths to written .h5ad files
    """
    if not download_dir.exists():
        raise FileNotFoundError(f"{download_dir} does not exist")

    raw_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []


    for sample_dir in download_dir.iterdir():
        if not sample_dir.is_dir():
            print(f"WARNING : {sample_dir} is not a directory; passing")
            continue

        # set the content of subdir to matrix.mtx.gz, features/genes.mtx.gz, barcodes.mtx.gz
        _standardize_10x_filenames(sample_dir)

        if not _looks_like_10x_dir(sample_dir):
            print(f"WARNING : {sample_dir} is not a 10x directory; passing")
            continue

        sample_id = sample_dir.stem

        ad = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_ids",
        )
        ad.var_names_make_unique()
        out_file = raw_dir / f"{label}_{sample_id}.h5ad"
        ad.write(out_file)
        print(f"Writing {out_file.stem}")

        outputs.append(out_file)

    return outputs   

def _looks_like_10x_dir(p: Path) -> bool:
    if not p.is_dir():
        return False

    for name in ("matrix.mtx", "matrix.mtx.gz"):
        if (p / name).exists():
            return True

    return False

def _standardize_10x_filenames(sample_dir: Path) -> None:
    mapping = {
        "matrix": "matrix.mtx.gz",
        "features": "features.tsv.gz",
        "genes": "genes.tsv.gz",
        "barcodes": "barcodes.tsv.gz",
    }

    for p in sample_dir.iterdir():
        if not p.is_file():
            continue

        for key, target in mapping.items():
            if key in p.name:
                p.rename(sample_dir / target)
                break

"""
    Codes related to CSV data i/o
    - _normalize_text_files
    - _build_csv_anndata (final wrapper)


"""

from src.files_utils import _gunzip_decompress, _tsv_to_csv
def _normalize_text_files(
    download_dir: Path,
    remove_original: bool,
) -> list[Path]:
    """
    Normalize text formats into clean CSV files.
    """
    if not download_dir.exists():
        raise FileNotFoundError(f"{download_dir} does not exist")

    csv_outputs: list[Path] = []

    if download_dir.is_file():
        path = download_dir
        out_dir = download_dir.parent

        path = _gunzip_decompress(path=path, remove_original=remove_original)
        path = _tsv_to_csv(path=path, out_dir=out_dir, remove_original=remove_original)
        csv_outputs.append(path)

    else:
        for path in download_dir.iterdir():
            if not path.is_file():
                continue

            path = _gunzip_decompress(path=path, remove_original=remove_original)
            path = _tsv_to_csv(path=path, out_dir=download_dir, remove_original=remove_original)
            csv_outputs.append(path)

    if not csv_outputs:
        raise RuntimeError(
            f"No CSV-compatible files found in {download_dir}"
        )

    return csv_outputs


def _build_csv_anndata(
    files: list[Path],
    raw_dir: Path,
    label: str,
) -> list[Path]:
    """
        Read .csv -> .h5ad
        Convert multiple CSV files from download_dir into AnnData objects
        and save them into raw_dir.

        Returns
        -------
        list[Path]
            Paths to written .h5ad files
    """
    raw_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []

    for path in files:
        if not path.is_file():
            continue
        if path.suffix.lower() != ".csv":
            continue
        sample_id = path.stem

        ad = sc.read_csv(path)

        out_file = raw_dir / f"{label}_{sample_id}.h5ad"
        ad.write(out_file)
        print(f"Writing {out_file.stem}")

        outputs.append(out_file)

    return outputs


