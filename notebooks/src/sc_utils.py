from pathlib import Path
from env_utils import get_temp_dir, extract_tar
import scanpy as sc
import shutil, subprocess

"""
    Download .h5ad files from link : 
    sc_download_mtx : download mtx-type files
    sc_download_tsv : download tsv-type files
    sc_download_csv : download csv-type files

    Load AnnData object from .h5ad files :
    sc_load_fix_h5ad : single dataset i/o
    sc_compile_from_dir : multiple datasets concatenation
    sc_compile_from_tar : multiple datasets concatenation (packed in tar)

"""

FILE_MAP = {
    "matrix.mtx": "matrix.mtx",
    "barcodes.tsv": "barcodes.tsv",
    "genes.tsv": "features.tsv",
    "features.tsv": "features.tsv",
}

def sc_process_downloaded_10x_mtx(
    download_dir: Path,
    raw_dir: Path,
    label: str, 
    sample_meta: dict,
    cleanup: bool = True
): # import to notebook
    
    groups = group_by_prefix(download_dir)

    for sample_id, files in groups.items():
       sample_dir = normalize_sample_files(files, raw_dir, sample_id)
       out_file = build_and_save_anndata(
        sample_dir,
        raw_dir,
        label,
        sample_id,
        sample_meta
        )
       print(f"Saved {out_file}")
       if cleanup:
        shutil.rmtree(sample_dir)

def download_tar_from_link(
    filename: str,
    link: str,
    download_dir: Path,
) -> None:
    subprocess.run(["curl", "-L", link, "-o", f"{filename}.tar"], check=True)
    subprocess.run(["tar", "-xf", f"{filename}.tar", "-C", download_dir], check=True)
    Path(f"{filename}.tar").unlink()


def normalize_10x_filename(filename: str) -> str | None:
    for key, target in FILE_MAP.items():
        if key in filename:
            if filename.endswith(".gz"):
                return f"{target}.gz"
            return target
    return None


def group_by_prefix(download_dir: Path) -> dict[str, list[Path]]:
    groups: dict[str, list[Path]] = {}
    for f in download_dir.iterdir():
        if not f.is_file():
            continue
        prefix = f.name.split("_")[0]
        groups.setdefault(prefix, []).append(f)
    return groups

def normalize_sample_files(
    files: list[Path],
    raw_dir: Path,
    sample_id: str
):
    sample_dir = raw_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    for f in files:
        normalized = normalize_10x_filename(f.name)
        dest = sample_dir / (normalized or f.name)
        shutil.move(f, dest)
    
    return sample_dir

def build_and_save_anndata(
    sample_dir: Path,
    raw_dir: Path,
    label: str, 
    sample_id: str, 
    sample_meta: dict,
):
    ad = sc.read_10x_mtx(
        sample_dir,
        var_names="gene_ids",
        make_unique=True,
    )
    
    if sample_id in sample_meta:
        for k, v in sample_meta[sample_id].items():
            ad.obs[k] = v if isinstance(v, (list, tuple)) else [v] * ad.n_obs
    else:
        print(f"Warning: no metadata for {sample_id}")
    
    out_file = raw_dir / f"{label}_{sample_id}.h5ad"
    ad.write(out_file)
    return out_file


def sc_download_tsv(filename: str, link, download_dir): # import to notebook


def sc_download_csv(filename: str, link, download_dir): # import to notebook


def tsv_to_csv():


def sc_needs_transpose(adata):
    if adata.X.shape == (adata.n_obs, adata.n_vars):
        return False
    elif adata.X.shape == (adata.n_vars, adata.n_obs):
        return True
    else:
        raise ValueError("AnnData X shape incompatible with obs/var")

def sc_fix_orientation(adata):
    if sc_needs_transpose(adata):
        adata = adata.T.copy()
        adata.var_names_make_unique()
    return adata

def sc_load_fix_h5ad(path): # import to notebook 
    """
        single file i/o
        input : h5ad path (str)
        return : adata
    """
    path = Path(path)
    adata = sc.read_h5ad(path)
    adata = sc_fix_orientation(adata)
    adata.obs["source_file"] = path.stem
    return adata

def sc_concat_adata(adatas, label="sample", merge="first"):
    adata = sc.concat(
        adatas,
        label=label,
        merge=merge
    )
    adata.obs_names_make_unique()
    return adata

def sc_compile_from_tar(tar_path, label="sample", merge="first"): # import to notebook 
    """
        directory compilation (tar)
        input : tar.gz path (str)
        return : concatenated adata
    """
    tar_path = Path(tar_path)
    temp_dir = get_temp_dir()
    try:
        extract_tar(tar_path, temp_dir)
        paths = sorted(temp_dir.glob("*.h5ad"))

        if not paths:
            raise FileNotFoundError("No .h5ad files found in tar")

        adatas = {p.stem: sc_load_fix_h5ad(p) for p in paths}
        adata = sc_concat_adata(adatas, label=label, merge=merge)
        return adata
    
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def sc_compile_from_dir(dir_path, label="sample", merge="first"): # import to notebook 
    """
        directory compilation (h5ad)
        input : dir (.h5ad file directory) path (str)
        return : concatenated adata
    """
    dir_path = Path(dir_path)
    paths = sorted(dir_path.glob("*.h5ad"))

    if not paths:
        raise FileNotFoundError("No .h5ad files found")

    adatas = {p.stem: sc_load_fix_h5ad(p) for p in paths}
    adata = sc_concat_adata(adatas, label=label, merge=merge)
    return adata
