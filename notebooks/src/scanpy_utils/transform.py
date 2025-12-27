import scanpy as sc

def _needs_transpose(adata):
    if adata.X.shape == (adata.n_obs, adata.n_vars):
        return False
    elif adata.X.shape == (adata.n_vars, adata.n_obs):
        return True
    else:
        raise ValueError("AnnData X shape incompatible with obs/var")

def _fix_orientation(adata):
    if _needs_transpose(adata):
        adata = adata.T.copy()
        adata.var_names_make_unique()
    return adata