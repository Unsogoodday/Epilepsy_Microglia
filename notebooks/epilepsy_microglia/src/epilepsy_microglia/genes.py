"""Source-grounded gene symbols without dropping or merging features."""
import pandas as pd
from anndata.utils import make_index_unique


def ensembl_mask(names):
    return pd.Index(names).str.match(r'^ENS[A-Z]*G\d+(?:\.\d+)?$')


def standardize_var_names(adata):
    """Replace Ensembl indices using var.gene_symbol; preserve source IDs.

    Missing mappings fail explicitly. Duplicate symbols receive numeric suffixes
    while the unmodified symbols remain in gene_symbol. Expression is untouched.
    """
    mask = ensembl_mask(adata.var_names)
    if not mask.any():
        return
    if 'gene_symbol' not in adata.var:
        raise ValueError('Ensembl var_names require source-grounded var.gene_symbol mappings')
    symbols = adata.var['gene_symbol'].astype('string')
    invalid = symbols.isna() | symbols.str.strip().eq('') | symbols.str.match(r'^ENS[A-Z]*G\d+(?:\.\d+)?$')
    if invalid[mask].any():
        raise ValueError('Missing gene symbols for Ensembl IDs: ' + ', '.join(adata.var_names[mask & invalid.fillna(True).to_numpy(dtype=bool)][:10]))
    names = adata.var_names.to_numpy(dtype=object, copy=True)
    names[mask] = symbols[mask].str.strip().to_numpy()
    if 'source_gene_id' not in adata.var:
        adata.var['source_gene_id'] = adata.var_names.to_numpy(copy=True)
    adata.var_names = make_index_unique(pd.Index(names, name='gene_name'))
