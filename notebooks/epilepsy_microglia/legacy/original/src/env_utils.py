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

def get_paths(DATASET_NAME: str):
    env = detect_env()

    if env == "colab":
        data_root = Path("/content/drive/MyDrive/data/active")
    elif env == "code-server":
        data_root = Path("/home/neuro_demo_research/data/active")
    else:
        data_root = Path.home() / "data" / "active"

    base = data_root / DATASET_NAME

    return {
        "base": base,
        "raw": base / "raw",
        "processed": base / "processed",
        "plots": base / "plots",
    }

def get_temp_dir():
    env = detect_env()
    if env == "colab":
        return Path("/content/data")
    elif env == "code-server":
        return Path("/home/neuro_demo_research/temp_data")
    else:
        return Path("~/neuro_demo_research").expanduser()


