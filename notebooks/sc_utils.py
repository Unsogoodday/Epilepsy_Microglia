from pathlib import Path
from env_utils import get_temp_dir, extract_tar
import scanpy as sc
import shutil

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

def sc_load_fix_h5ad(path):
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

def sc_compile_from_tar(tar_path, label="sample", merge="first"):
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


def sc_compile_from_dir(dir_path, label="sample", merge="first"):
    dir_path = Path(dir_path)
    paths = sorted(dir_path.glob("*.h5ad"))

    if not paths:
        raise FileNotFoundError("No .h5ad files found")

    adatas = {p.stem: sc_load_fix_h5ad(p) for p in paths}
    adata = sc_concat_adata(adatas, label=label, merge=merge)
    return adata
