import pandas as pd
import mygene
import scanpy as sc
import AnnData as ad
from scipy.stats import fisher_exact

def cleanup_annotation(
    adata,
    old_id,
    new_id="ensembl_id",
):
    """
    Sanitize Ensembl IDs by removing version suffixes (e.g. ENSG...1 -> ENSG...).

    Input:
        adata.var[old_id] must exist.

    Output:
        Adds adata.var[new_id] with cleaned IDs.
        Row alignment is preserved.
        Returns a copy of adata.
    """
    adata = adata.copy()

    if old_id not in adata.var:
        raise KeyError(f"{old_id} not found in adata.var")

    adata.var[new_id] = (
        adata.var[old_id]
        .astype(str)
        .str.split(".", n=1)
        .str[0]
        .replace("nan", pd.NA)
    )

    return adata


def map_ensembl_to_symbol(
    adata,
    ensembl_key="ensembl_id",
    symbol_key="gene_symbol",
    species="human",
    chunk_size=1000,
):
    """
    Map Ensembl gene IDs to gene symbols using MyGeneInfo.

    Input:
        adata.var[ensembl_key] must exist.

    Output:
        Adds adata.var[symbol_key] (NaN if unmapped).
        Returns a copy of adata.
    """
    adata = adata.copy()

    if ensembl_key not in adata.var:
        raise KeyError(f"{ensembl_key} not found in adata.var")

    mg = mygene.MyGeneInfo()
    genes = adata.var[ensembl_key].astype(str).tolist()

    results = []
    for i in range(0, len(genes), chunk_size):
        chunk = genes[i : i + chunk_size]
        res = mg.querymany(
            chunk,
            scopes="ensembl.gene",
            fields="symbol",
            species=species,
            as_dataframe=False,
        )
        results.extend(res)

    df = (
        pd.DataFrame(results)
        .loc[lambda x: x.get("symbol").notna()]
        .drop_duplicates(subset="query")
        .rename(columns={"query": ensembl_key, "symbol": symbol_key})
        [[ensembl_key, symbol_key]]
    )

    adata.var[symbol_key] = adata.var[ensembl_key].map(
        df.set_index(ensembl_key)[symbol_key]
    )

    return adata

def guardrail_unmapped_hvgs(
    adata,
    hvg_key="highly_variable",
    mapping_key="mapping_status",
    mapped_label="mapped",
    max_unmapped_hvg_fraction=0.05,
    min_unmapped_genes=30,
    or_fail_threshold=0.2,
    pseudocount=0.5,
):
    """
    Guardrail: ensure dropping unmapped genes won't drop too many HVGs
    and unmapped genes aren't strongly enriched for HVGs.

    Returns: (decision, report_dict)
      decision in {"PASS","FAIL"}
    """
    var = adata.var
    if hvg_key not in var or mapping_key not in var:
        raise KeyError("Missing required columns in adata.var")

    is_hvg = var[hvg_key].fillna(False).to_numpy()
    is_mapped = (var[mapping_key] == mapped_label).to_numpy()
    is_unmapped = ~is_mapped

    a = int(np.sum(is_mapped & is_hvg))       # mapped + HVG
    b = int(np.sum(is_mapped & ~is_hvg))      # mapped + not HVG
    c = int(np.sum(is_unmapped & is_hvg))     # unmapped + HVG
    d = int(np.sum(is_unmapped & ~is_hvg))    # unmapped + not HVG

    total_hvgs = a + c
    total_unmapped = c + d

    unmapped_hvg_fraction = (c / total_hvgs) if total_hvgs else 0.0
    odds_ratio = ((a + pseudocount) * (d + pseudocount)) / ((b + pseudocount) * (c + pseudocount))

    fail = False
    reasons = []

    if total_unmapped >= min_unmapped_genes:
        if unmapped_hvg_fraction > max_unmapped_hvg_fraction:
            fail = True
            reasons.append("too_many_hvgs_unmapped")
        if odds_ratio < or_fail_threshold:
            fail = True
            reasons.append("hvg_enriched_in_unmapped")

    decision = "FAIL" if fail else "PASS"
    report = {
        "a_mapped_hvg": a, "b_mapped_not_hvg": b, "c_unmapped_hvg": c, "d_unmapped_not_hvg": d,
        "unmapped_hvg_fraction": unmapped_hvg_fraction,
        "odds_ratio": odds_ratio,
        "decision": decision,
        "reasons": reasons,
    }
    return decision, report
