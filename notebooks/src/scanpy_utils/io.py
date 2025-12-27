from pathlib import Path
import scanpy as sc
import pandas as pd
import shutil, subprocess, gzip

from files_utils import _download_tar_from_link, _download_from_link
def _download_dataset(
    *,
    filename: str,
    link: str,
    download_dir: Path,
    is_tar: bool,
) -> None:
    if is_tar:
        _download_tar_from_link(filename, link, download_dir)
    else:
        _download_from_link(filename, link, download_dir)

from files_utils import _gunzip_decompress, _tsv_to_csv
def _normalize_text_files(
    download_dir: Path,
    remove_original: bool,
) -> list[Path]:
    """
        Normalize text formats into clean CSV file
    """
    if not download_dir.exists():
        raise FileNotFoundError(f"{download_dir} does not exist")
    
    csv_outputs: list[Path] = []

    for path in download_dir.iterdir():
        if not path.is_file():
            continue

        path = _gunzip_decompress(path=path, remove_original=remove_original)
        path = _tsv_to_csv(path=path, out_dir=download_dir, remove_original=remove_original)
        csv_outputs.append[path]

    if not csv_outputs:
        raise RuntimeError(
            f"No CSV-compatible files found in {download_dir}"
        )

    return csv_outputs


def _build_mtx_anndata(
    download_dir: Path,
    raw_dir: Path,
    label: str, 
) -> list[Path]:
    """
        Read matrix.mtx / features.tsv / genes.tsv -> .h5ad
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
            continue

        sample_id = sample_dir.stem

        ad = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_ids",
        )
        ad.var_names_make_unique()
        out_file = raw_dir / f"{label}_{sample_id}.h5ad"
        ad.write(out_file)  

        outputs.append(out_file)

    return outputs   

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

    for path in files():
        if not path.is_file():
            continue
        if path.suffix.lower() != ".csv":
            continue
        sample_id = path.stem

        ad = sc.read_csv(
            path,
            var_names="gene_ids",
        )
        ad.var_names_make_unique()
        out_file = raw_dir / f"{label}_{sample_id}.h5ad"
        ad.write(out_file)

        outputs.append(out_file)

    return outputs


