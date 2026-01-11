import scanpy as sc

def _looks_like_gene_ids(names):
    """
        Simple heuristic method to check if a column / row is gene_ids
    """
    return sum(name.isupper() and len(name) < 20 for name in names[:50]) > 20

def _fix_orientation(adata):
    if adata.uns.get("orientation_fixed"):
        return adata

    obs_gene = _looks_like_gene_ids(adata.obs_names)
    var_gene = _looks_like_gene_ids(adata.var_names)

    if obs_gene and not var_gene:
        adata = adata.T.copy()
        adata.uns["orientation_fixed"] = True
        return adata

    adata.uns["orientation_fixed"] = True
    return adata

