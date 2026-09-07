from pathlib import Path
from src.env_utils import get_temp_dir
import shutil, subprocess, gzip
import pandas as pd

def is_tar_file(path: Path) -> bool:
    result = subprocess.run(
        ["tar", "-tf", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0

def _extract_tar(src_tar: Path, dst_dir: Path):
    dst_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(src_tar, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue  # skip dirs, symlinks, etc.

            target = dst_dir / Path(member.name).name

            with tar.extractfile(member) as src, open(target, "wb") as dst:
                dst.write(src.read())

def _download_tar_from_link( 
    link: str,
    download_dir: Path,
) -> None:
    download_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        ["curl", "-L", "-O", "-J", link], # or --remote-name --content-disposition
        cwd=download_dir,
        check=True,
    )

    tar_files = [
        p for p in download_dir.iterdir()
        if p.is_file() and is_tar_file(p)
    ]

    if not tar_files:
        raise RuntimeError("No tar archives found")

    for tar_path in tar_files:
        subprocess.run(
            ["tar", "-xf", str(tar_path)],
            cwd=download_dir,
            check=True,
        )


def _download_from_link(
    link: str,
    download_dir: Path,
) -> None:
    download_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        ["curl", "-L", "-O", link],
        cwd=download_dir,
        check=True,
    )

def _gunzip_decompress(
    *,
    path: Path,
    remove_original: bool = True,
) -> Path:
    """
    Decompress a .gz file and return the decompressed path.
    If the file is not gzipped, return it unchanged.
    """
    if path.suffix != ".gz":
        return path

    decompressed = path.with_suffix("")
    if decompressed.exists():
        raise FileExistsError(
            f"Decompressed file already exists: {decompressed}"
        )

    with gzip.open(path, "rb") as f_in, open(decompressed, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    if remove_original:
        path.unlink()

    return decompressed


def _tsv_to_csv(
    *,
    path: Path,
    out_dir: Path | None = None,
    remove_original: bool = True,
) -> Path:
    """
    Convert a single TSV file to CSV and return the CSV path.
    If the file is not a TSV, return it unchanged.
    """
    if path.suffix != ".tsv":
        return path

    out_dir = out_dir or path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / f"{path.stem}.csv"
    if csv_path.exists():
        raise FileExistsError(
            f"CSV already exists: {csv_path}"
        )

    df = pd.read_csv(path, sep="\t")
    df.to_csv(csv_path, index=False)
    print(f"Turning {path.stem}.tsv to csv")

    if remove_original:
        path.unlink()

    return csv_path

