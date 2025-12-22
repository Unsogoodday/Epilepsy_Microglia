from pathlib import Path
import os
import tarfile
import glob

def detect_env():
    # Colab
    if (
        "COLAB_RELEASE_TAG" in os.environ
        or os.path.exists("/content")
    ):
        return "colab"

    # code-server/VS Code
    if (
        "VSCODE_IPC_HOOK_CLI" in os.environ
        or os.environ.get("TERM_PROGRAM") == "vscode"
    ):
        return "code-server"
    
    return "local"

def get_paths(PROJECT_NAME):
    env = detect_env()

    if env == "colab":
        base = Path("/content/drive/MyDrive/datas")
    elif env == "code-server":
        base = Path("/home/neuro_demo_research/data_from_drive")
    else:
        base = Path("~/neuro_demo_research").expanduser()

    return {
        "base" : base / PROJECT_NAME,
        "raw" : base / PROJECT_NAME / "raw",
        "processed" : base / PROJECT_NAME / "processed",
        "plots" : base / PROJECT_NAME / "plots",
    }

def get_temp_dir():
    env = detect_env()
    if env == "colab":
        return Path("/content/data")
    elif env == "code-server":
        return Path("/home/neuro_demo_research/temp_data")
    else:
        return Path("~/neuro_demo_research").expanduser()

def extract_tar(src_tar: Path, dst_dir: Path):
    dst_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(src_tar, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue  # skip dirs, symlinks, etc.

            target = dst_dir / Path(member.name).name

            with tar.extractfile(member) as src, open(target, "wb") as dst:
                dst.write(src.read())
