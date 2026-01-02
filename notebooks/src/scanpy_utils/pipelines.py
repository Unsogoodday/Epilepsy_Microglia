from src.scanpy_utils.io import _download_dataset, _normalize_text_files, _build_mtx_anndata, _build_csv_anndata
from src.scanpy_utils.transform import _fix_orientation
from src.scanpy_utils.annotations import _cleanup_annotation, _map_ensembl_to_symbol, _guardrail_unmapped_hvgs
from src.scanpy_utils.concat import _concat_adata

from pathlib import Path
import scanpy as sc
import numpy as np
import matplotlib as plt

"""
    Pipelines and supporting codes
    Import directly to notebook

    io
    - sc_download_mtx : download mtx type data, fix orientation, write in raw_dir
    - sc_download_tsv : download tsv type data, fix orientation, write in raw_dir

    qc
    - sc_plot_qc_mt_rb_hb : compute mt, rb, hb gene pct and plot
    - sc_build_qc_mask : build qc mask with no. of MAD, mt/rb/hb threshold
    - sc_apply_qc_mask : apply the qc mask
    - sc_plot_qc_distributions : plot qc changes before/after qc mask

    annotation 
    - sc_annotate_mygene : match ENSEMBL -> SYMBOL via mygene

"""

def download_mtx(
    *,
    filename: str,
    link: str,
    download_dir: Path,
    raw_dir: Path,
    label: str,
    is_tar: bool,
) -> None:
    """
    Orchestrator for downloading datasets and building AnnData
    from 10x-style matrix directories.
    """
    _download_dataset(
        filename=filename,
        link=link,
        download_dir=download_dir,
        is_tar=is_tar,
    )
    adata_paths = _build_mtx_anndata(
        download_dir=download_dir,
        raw_dir=raw_dir,
        label=label,
    )
    h5ad_paths = list(h5ad_paths)  
    for f in h5ad_paths:
        print("READING:", f)
        if f.suffix != ".h5ad":
            raise ValueError(f"{f} is not a proper h5ad format")

        name = f.name

        if name.endswith(".raw.h5ad"):
            base = name.removesuffix(".raw.h5ad")
        else:
            base = f.stem  

        out_path = f.with_name(f"{base}.fixed.h5ad")

        adata = sc.read_h5ad(f)
        adata = _fix_orientation(adata)
        adata.write(out_path)
        print("WRITING:", out_path)

def download_tsv_csv(
    *,
    filename: str,
    link: str,
    download_dir: Path,
    raw_dir: Path,
    label: str,
    is_tar: bool,
    remove_original: bool=True,
) -> None:
    """
    Orchestrator for downloading datasets and building AnnData
    from csv/tsv files
    """
    _download_dataset(
        link=link,
        download_dir=download_dir,
        is_tar=is_tar,
    )
    csv_paths = _normalize_text_files(
        download_dir=download_dir,
        remove_original=remove_original,
    )
    h5ad_paths = _build_csv_anndata(
        files=csv_paths,
        raw_dir=raw_dir,
        label=label,
    )
    h5ad_paths = list(h5ad_paths)  
    for f in h5ad_paths:
        print("READING:", f)
        if f.suffix != ".h5ad":
            raise ValueError(f"{f} is not a proper h5ad format")

        name = f.name

        if name.endswith(".raw.h5ad"):
            base = name.removesuffix(".raw.h5ad")
        else:
            base = f.stem  

        out_path = f.with_name(f"{base}.fixed.h5ad")

        adata = sc.read_h5ad(f)
        adata = _fix_orientation(adata)
        adata.write(out_path)
        print("WRITING:", out_path)
    
def load_h5ad(path):
    """
        single file i/o
        input : h5ad path (str)
        return : adata
    """
    path = Path(path)
    adata = sc.read_h5ad(path)
    adata.obs["source_file"] = path.stem
    return adata

from src.files_utils import _extract_tar
def compile_from_tar(tar_path, label="sample", merge="first"):
    """
        directory compilation (tar)
        input : tar.gz path (str)
        return : concatenated adata
    """
    tar_path = Path(tar_path)
    temp_dir = get_temp_dir()
    try:
        _extract_tar(tar_path, temp_dir)
        paths = sorted(temp_dir.glob("*.h5ad"))

        if not paths:
            raise FileNotFoundError("No .h5ad files found in tar")

        adatas = {p.stem: load_h5ad(p) for p in paths}
        adata = _concat_adata(adatas, label=label, merge=merge)
        return adata
    
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def compile_from_dir(dir_path, label="sample", merge="first"):
    """
        directory compilation (h5ad)
        input : dir (.h5ad file directory) path (str)
        return : concatenated adata
    """
    dir_path = Path(dir_path)
    paths = sorted(dir_path.glob("*.h5ad"))

    if not paths:
        raise FileNotFoundError("No .h5ad files found")

    adatas = {p.stem: load_h5ad(p) for p in paths}
    adata = _concat_adata(adatas, label=label, merge=merge)
    return adata

