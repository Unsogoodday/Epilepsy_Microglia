import scanpy as sc

def _looks_like_gene_ids(names):
    """
        Simple heuristic method to check if a column / row is gene_ids
    """
    return sum(name.isupper() and len(name) < 20 for name in names[:50]) > 20

def _fix_orientation(adata):
    """
        v2 (heuristic approach)
        USE ONLY FOR INTERNALLY CONSISTENT ANNDATA OBJECTS
    """
    if _looks_like_gene_ids(adata.obs_names):
        adata = adata.T.copy()
    return adata