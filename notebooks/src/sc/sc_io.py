from pathlib import Path
from env_utils import get_temp_dir, extract_tar
import scanpy as sc
import shutil, subprocess

def download_tar_from_link(
    filename: str,
    link: str,
    download_dir: Path,
) -> None:
    download_dir.mkdir(parents=True, exist_ok=True)
    tar_path = download_dir / f"{filename}.tar"
    subprocess.run(["curl", "-L", link, "-o", str(tar_path)], check=True)
    subprocess.run(["tar", "-xf", str(tar_path), "-C", str(download_dir)], check=True)
    tar_path.unlink()


def build_and_save_anndata(
    download_dir: Path,
    raw_dir: Path,
    label: str, 
    sample_meta: dict,
):
    groups: dict[str, list[Path]] = {}
    for sample_dir in download_dir.iterdir():
        if not sample_dir.is_dir():
            continue

        ad = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_ids",
            make_unique=True,
        )

        sample_id = sample_dir.name
        if sample_id in sample_meta:
            for k, v in sample_meta[sample_id].items():
                ad.obs[k] = v if isinstance(v, (list, tuple)) else [v] * ad.n_obs
        else:
            print(f"Warning: no metadata for {sample_id}")
        
        out_file = raw_dir / f"{label}_{sample_id}.h5ad"
        ad.write(out_file)       

def tsv_to_csv():

