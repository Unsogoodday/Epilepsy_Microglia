import scanpy as sc

def concat_adata(adatas, label="sample", merge="first"):
    adata = sc.concat(
        adatas,
        label=label,
        merge=merge
    )
    adata.obs_names_make_unique()
    return adata