from pathlib import Path
from env_utils import get_temp_dir
import shutil, subprocess, gzip

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


def _download_from_link(
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

    if remove_original:
        path.unlink()

    return csv_path

