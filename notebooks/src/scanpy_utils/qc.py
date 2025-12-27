#미완

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

def sc_annotate_mygene(
    adata,
    old_id,
):
    """
    - clean Ensembl IDs
    - map to gene symbols
    - mark mapping status
    - run guardrail
    """
    adata = cleanup_annotation(adata, old_id)
    adata = map_ensembl_to_symbol(adata)

    adata.var["mapping_status"] = np.where(
        adata.var["gene_symbol"].isna(),
        "unmapped",
        "mapped",
    )

    decision, report = guardrail_unmapped_hvgs(adata)

    if decision == "FAIL":
        raise ValueError(f"Gene mapping failed QC: {report}")

    return adata, report