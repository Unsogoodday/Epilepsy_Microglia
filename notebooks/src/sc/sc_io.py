from pathlib import Path
from env_utils import get_temp_dir, extract_tar
import scanpy as sc
import pandas as pd
import shutil, subprocess

import subprocess
from pathlib import Path

def download_tar_from_link(
    filename: str,
    link: str,
    download_dir: Path,
) -> None:
    download_dir.mkdir(parents=True, exist_ok=True)
    tar_path = download_dir / filename

    subprocess.run(
        ["curl", "-L", link, "-o", str(tar_path)],
        check=True
    )
    subprocess.run(
        ["tar", "-xf", str(tar_path), "-C", str(download_dir)],
        check=True
    )
    tar_path.unlink()


def download_from_link(
    filename: str,
    link: str,
    download_dir: Path,
) -> None:
    download_dir.mkdir(parents=True, exist_ok=True)
    download_path = download_dir / filename

    subprocess.run(
        ["curl", "-L", link, "-o", str(download_path)],
        check=True
    )

def download_dataset(
    *,
    filename: str,
    link: str,
    download_dir: Path,
    is_tar: bool,
) -> None:
    if is_tar:
        download_tar_from_link(filename, link, download_dir)
    else:
        download_from_link(filename, link, download_dir)


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

def tsv_to_csv(
    filename: str,
    download_dir: Path,
    out_dir=download_dir,
    remove_original=True,
):
    for tsv_path in download_dir.glob("*tsv"):
        basename = os.path.splitext(os.path.basename(tsv_path))[0]
        csv_path = os.path.join(out_dir, f"{basename}.csv")
        try: 
            df = pd.read_csv(tsv_path, sep="\t")
            df.to_csv(csv_path, index=False)
            if(remove_original):
                os.remove(tsv_path)
                print(f"Converted {tsv_path} -> {csv_path}, deleted original")
            else:
                print(f"Converted {tsv_path} -> {csv_path}, kept original")
        except Exception as e:
            print(f"Error processing {tsv_path}: {e}")

