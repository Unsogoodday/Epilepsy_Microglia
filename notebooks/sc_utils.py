from pathlib import Path
from env_utils import detect_env, get_paths
import os
import tarfile
import glob
import scanpy as sc
import anndata as ad

def sc_unpack_tar(RAW_DIR, FCODE):
    # unpack tar.gz to adata 
    env = detect_env()
    if(env == "Colab"):
        TEMP_DIR = "/contents/data"
    elif(env == "code-server"):
        TEMP_DIR = "/home/neuro_demo_research/temp_data"
    else:
        TEMP_DIR = Path("~/neuro_demo_research").expandsuer()
    os.makedirs(TEMP_DIR, exist_ok=True)

    # unpack tar from RAW_DIR
    src = os.path.join(RAW_DIR, f"{FCODE}.tar.gz")
    with tarfile.open(src, "r:gz") as tar:
        for member in tar.getmembers():
            member.name = os.path.basename(member.name)
            tar.extract(memeber, TEMP_DIR)
    
    adatas = {}

    # mark each j5ad file with sample_id
    for file in glob.glob(os.path.join(TEMP_DIR, "*.h5ad")):
        sample_adata = sc.read_h5ad(file)
        sample_adata.var_names_make_unique()
        sample_id = sample_adata.obs["sample_id"].iloc[0]
        adatas[sample_id] = sample_adata

    # concatenate
    if adatas:
        adata = sc.concat(
            adatas.values(),
            label="sample_id",
            keys=adatas.keys(),
            merge="same"
        )
        adata.obs_names_make_unique()
        print(adata.obs["sample_id"].value_counts())
        
        # save adata back to RAW_DIR
        adata.write(os.path.join(RAW_DIR, f"{FCODE}.h5ad"))
    else:
        print("No .h5ad files found in the directory.")

    os.remove(TEMP_DIR)

