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

from pathlib import Path
import pandas as pd

def tsv_to_csv(
    *,
    download_dir: Path,
    out_dir: Path | None = None,
    remove_original: bool = True,
) -> None:
    """
    Convert all TSV files in download_dir to CSV.

    Parameters
    ----------
    download_dir : Path
        Directory containing .tsv files.
    out_dir : Path | None
        Output directory for .csv files. Defaults to download_dir.
    remove_original : bool
        Whether to delete .tsv files after successful conversion.
    """

    out_dir = out_dir or download_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    for tsv_path in download_dir.glob("*.tsv"):
        csv_path = out_dir / f"{tsv_path.stem}.csv"

        df = pd.read_csv(tsv_path, sep="\t")
        df.to_csv(csv_path, index=False)

        if remove_original:
            tsv_path.unlink()


def check_csv_orientation(
    csv_path: Path,
    n_check: int = 5,
) -> str:
    """
    Infer gene orientation in a CSV file.

    Returns
    -------
    str
        One of {"genes_as_rows", "genes_as_columns", "unknown"}
    """
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)

    df = pd.read_csv(csv_path, index_col=0, nrows=n_check)

    if df.empty:
        raise ValueError(f"{csv_path} is empty or unreadable")

    row_label = str(df.index[0])
    col_label = str(df.columns[0])

    def looks_like_gene(label: str) -> bool:
        return (
            label.startswith("ENSG")
            or label.startswith("MT-")
            or label.replace("-", "").isalnum()
        )

    row_gene_like = looks_like_gene(row_label)
    col_gene_like = looks_like_gene(col_label)

    if row_gene_like and not col_gene_like:
        return "genes_as_rows"
    if col_gene_like and not row_gene_like:
        return "genes_as_columns"

    return "unknown"

from sc_transform import fix_orientation
def csv_to_h5ad(
    *,
    csv_path: Path,
    raw_dir: Path,
    label: str, 
    sample_meta: dict,
):
    for f in csv_path.glob("*.csv"):
        sample_id = f.stem
        ori = check_csv_orientation(f)

        if ori == "genes_as_rows":
            adata = sc.read_csv(f, first_column_names=True)
            adata = fix_orientation(adata)
        elif ori == "genes_as_columns":
            adata = sc.read_csv(f, first_column_names=True)
        else:
            raise ValueError(f"Unknown orientation for {f}")
        
        adata.obs["sample_id"] = sample_id
        
        # I have no idea how to manage the 'saving' part