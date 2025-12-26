from sc_io import download_tar_from_link, build_and_save_anndata
from sc_transform import fix_orientation

from sc_annotations import cleanup_annotation, map_ensembl_to_symbol, guardrail_unmapped_hvgs

import scanpy as sc
import numpy as np
import matplotlib as plt

"""
    Pipelines and supporting codes
    Import directly to notebook

    io
    - sc_download_mtx : download mtx type data
    - sc_download_tsv : download tsv type data
    - sc_download_csv : download csv type data

    qc
    - sc_plot_qc_mt_rb_hb : compute mt, rb, hb gene pct and plot
    - sc_build_qc_mask : build qc mask with no. of MAD, mt/rb/hb threshold
    - sc_apply_qc_mask : apply the qc mask
    - sc_plot_qc_distributions : plot qc changes before/after qc mask

    annotation 
    - sc_annotate_mygene : match ENSEMBL -> SYMBOL via mygene


"""

def sc_download_mtx(
    *,
    filename: str,
    link: str,
    download_dir: Path,
    raw_dir: Path,
    label: str,
    sample_meta: dict,
) -> None:
    """
    Orchestrator for downloading, extracting, and building AnnData
    from 10x-style tar archives.
    """

    download_tar_from_link(
        filename=filename,
        link=link,
        download_dir=download_dir,
    )

    build_and_save_anndata(
        download_dir=download_dir,
        raw_dir=raw_dir,
        label=label,
        sample_meta=sample_meta,
    )

def sc_download_tsv(filename: str, link, download_dir):

def sc_download_csv(filename: str, link, download_dir):

def sc_load_fix_h5ad(path):
    """
        single file i/o
        input : h5ad path (str)
        return : adata
    """
    path = Path(path)
    adata = sc.read_h5ad(path)
    adata = sc_fix_orientation(adata)
    adata.obs["source_file"] = path.stem
    return adata

def sc_compile_from_tar(tar_path, label="sample", merge="first"):
    """
        directory compilation (tar)
        input : tar.gz path (str)
        return : concatenated adata
    """
    tar_path = Path(tar_path)
    temp_dir = get_temp_dir()
    try:
        extract_tar(tar_path, temp_dir)
        paths = sorted(temp_dir.glob("*.h5ad"))

        if not paths:
            raise FileNotFoundError("No .h5ad files found in tar")

        adatas = {p.stem: sc_load_fix_h5ad(p) for p in paths}
        adata = concat_adata(adatas, label=label, merge=merge)
        return adata
    
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def sc_compile_from_dir(dir_path, label="sample", merge="first"):
    """
        directory compilation (h5ad)
        input : dir (.h5ad file directory) path (str)
        return : concatenated adata
    """
    dir_path = Path(dir_path)
    paths = sorted(dir_path.glob("*.h5ad"))

    if not paths:
        raise FileNotFoundError("No .h5ad files found")

    adatas = {p.stem: sc_load_fix_h5ad(p) for p in paths}
    adata = concat_adata(adatas, label=label, merge=merge)
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