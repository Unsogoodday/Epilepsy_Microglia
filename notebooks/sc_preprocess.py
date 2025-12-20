import pandas as pd
import mygene
import scanpy as sc
import AnnData as ad
from scipy.stats import fisher_exact

def sc_cleanup_annotation(
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

def sc_guardrail_unmapped_hvgs(
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


def sc_annotate_mygene(
    adata,
    old_id,
): # import to notebook
    """
    High-level orchestration:
    - clean Ensembl IDs
    - map to gene symbols
    - mark mapping status
    - run guardrail
    """
    adata = sc_cleanup_annotation(adata, old_id)
    adata = map_ensembl_to_symbol(adata)

    adata.var["mapping_status"] = np.where(
        adata.var["gene_symbol"].isna(),
        "unmapped",
        "mapped",
    )

    decision, report = sc_guardrail_unmapped_hvgs(adata)

    if decision == "FAIL":
        raise ValueError(f"Gene mapping failed QC: {report}")

    return adata, report


def sc_add_mt_ribo_hb_qc(adata, copy=True):
    """
    Annotate mitochondrial, ribosomal, and hemoglobin genes
    and compute standard Scanpy QC metrics.

    Assumes:
    - adata.X contains raw counts
    - adata.var_names are gene symbols (human-style)
    """
    if copy:
        adata = adata.copy()

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB(?!P)")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        log1p=True
    )

    return adata


def sc_plot_qc_mt_rb_hb(adata, dataset, date):
    sc.pl.violin(
        adata,
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
        ],
        jitter=0.4,
        multi_panel=True,
        save=f"{dataset}_qc_{date}.png",
    )

    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        save=f"{dataset}_qc_{date}.png",
    )


def sc_build_qc_mask(
    adata,
    n_mad,
    mt_threshold, 
    ribo_threshold, 
    hb_threshold, 
):
    def mad_mask(series, n_mad):
        x = np.log1p(series)
        med = np.median(x)
        mad = np.median(np.abs(x - med))
        return (x > med - n_mad * mad) & (x < med + n_mad * mad)

    mask_genes = mad_mask(adata.obs["n_genes_by_counts"], n_mad)
    mask_counts = mad_mask(adata.obs["total_counts"], n_mad)

    mask_mt = adata.obs["pct_counts_mt"] < mt_threshold
    mask_ribo = adata.obs["pct_counts_ribo"] < ribo_threshold
    mask_hb = adata.obs["pct_counts_hb"] < hb_threshold

    qc_mask = mask_genes & mask_counts & mask_mt & mask_ribo & mask_hb

    return qc_mask

def sc_apply_qc_mask(adata, qc_mask, copy=True):
    if copy:
        adata = adata.copy()
    return adata[qc_mask, :]

def sc_plot_qc_distributions(adata, qc_mask=None, bins=100):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].hist(np.log1p(adata.obs["n_genes_by_counts"]), bins=bins)
    axes[0].set_title("log(n_genes_by_counts) - before")

    axes[1].hist(np.log1p(adata.obs["total_counts"]), bins=bins)
    axes[1].set_title("log(total_counts) - before")

    plt.show()

    if qc_mask is not None:
        adata_qc = adata[qc_mask]

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].hist(np.log1p(adata_qc.obs["n_genes_by_counts"]), bins=bins)
        axes[0].set_title("log(n_genes_by_counts) - after")

        axes[1].hist(np.log1p(adata_qc.obs["total_counts"]), bins=bins)
        axes[1].set_title("log(total_counts) - after")

        plt.show()

