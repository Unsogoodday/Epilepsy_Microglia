from scipy.stats import spearmanr
import pandas as pd

def spearman_matrix_with_p(df):
    """
        input : df
        output : spearman rho, pval for each matrix element
    """

    cols = df.columns

    rho = pd.DataFrame(
        np.eye(len(cols)),
        index=cols,
        columns=cols,
    )

    pval = pd.DataFrame(
        np.zeros((len(cols), len(cols))),
        index=cols,
        columns=cols,
    )

    for i, x in enumerate(cols):
        for j, y in enumerate(cols):

            if i == j:
                rho.loc[x, y] = 1
                pval.loc[x, y] = 0
                continue

            valid = df[[x, y]].dropna()

            r, p = spearmanr(
                valid[x],
                valid[y],
            )

            rho.loc[x, y] = r
            pval.loc[x, y] = p

    return rho, pval

from statsmodels.stats.multitest import multipletests

def fdr_symmetric_matrix(pval):
    """
        fdr correction for symmetric matrix (ex. spearman, pearson correlation matrix)
        - input : symmetric matrix of pvals
        - output : qvals
    """

    cols = pval.columns
    qval = pd.DataFrame(
        np.nan,
        index=cols,
        columns=cols,
    )

    pairs = []
    pvalues = []

    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):

            pairs.append((cols[i], cols[j]))
            pvalues.append(pval.iloc[i, j])

    _, qvalues, _, _ = multipletests(
        pvalues,
        method="fdr_bh"
    )

    for (x, y), q in zip(pairs, qvalues):
        qval.loc[x, y] = q
        qval.loc[y, x] = q

    np.fill_diagonal(qval.values, 0)

    return qval

def star(q):
    if q < 0.001:
        return "***"
    elif q < 0.01:
        return "**"
    elif q < 0.05:
        return "*"
    else:
        return ""


def make_annot(rho, qval):
    """
        highlight statistically significant elements with stars
    """

    annot = rho.copy().astype(object)

    for i in rho.index:
        for j in rho.columns:

            if i == j:
                annot.loc[i, j] = "1.00"
            else:
                annot.loc[i, j] = (
                    f"{rho.loc[i, j]:.2f}"
                    f"{star(qval.loc[i, j])}"
                )

    return annot

def pc_loadings_both_sides(adata, pc=1, n=20):
    """
        print PC axis genes
    """
    s = pd.Series(
        adata.varm["PCs"][:, pc - 1],
        index=adata.var_names
    ).sort_values()

    negative = s.head(n).to_frame("loading")
    positive = s.tail(n).sort_values(ascending=False).to_frame("loading")

    return positive, negative

def print_pc_genes(adata, pcs=(1, 2), n=30):
    for pc in pcs:
        s = pd.Series(
            adata.varm["PCs"][:, pc - 1],
            index=adata.var_names
        )

        print(f"\n===== PC{pc} POSITIVE =====")
        print(", ".join(s.nlargest(n).index))

        print(f"\n===== PC{pc} NEGATIVE =====")
        print(", ".join(s.nsmallest(n).index))

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


def aggregate_by_sample(
    adata,
    groupby="sample_id",
    counts_layer="counts",
):
    """
        Generate 'groupby' -level pseudobulk
        - sum raw counts per sample
        - retain all obs keys
        - per sample mean of numerical vals
        - 'mode' for str/categorical vals
        - add obs key ['n_cells_pb']
        input : anndata object
        output : pb
    """
    # ------------------------------------------------------------
    # 1. Use raw counts
    # ------------------------------------------------------------
    X = adata.layers[counts_layer]

    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    # ------------------------------------------------------------
    # 2. Encode sample_id
    # ------------------------------------------------------------
    groups = pd.Categorical(adata.obs[groupby])
    sample_ids = groups.categories
    codes = groups.codes

    n_samples = len(sample_ids)
    n_cells = adata.n_obs

    # group × cell aggregation matrix
    G = sp.csr_matrix(
        (
            np.ones(n_cells, dtype=np.float32),
            (codes, np.arange(n_cells)),
        ),
        shape=(n_samples, n_cells),
    )

    # ------------------------------------------------------------
    # 3. Sum raw counts per sample
    # ------------------------------------------------------------
    X_pb = G @ X

    # ------------------------------------------------------------
    # 4. Aggregate obs while retaining ALL keys
    # ------------------------------------------------------------
    obs = adata.obs.copy()
    obs[groupby] = obs[groupby].astype(str)

    obs_pb = pd.DataFrame(index=sample_ids.astype(str))

    for col in obs.columns:

        # group key itself
        if col == groupby:
            obs_pb[col] = obs_pb.index
            continue

        s = obs[col]

        # Numeric columns -> mean within sample
        if pd.api.types.is_numeric_dtype(s):
            obs_pb[col] = (
                obs.groupby(groupby, observed=True)[col]
                .mean()
                .reindex(obs_pb.index)
            )

        # Categorical / string columns -> mode
        else:
            def get_mode(x):
                x = x.dropna()
                if len(x) == 0:
                    return np.nan
                return x.mode().iloc[0]

            obs_pb[col] = (
                obs.groupby(groupby, observed=True)[col]
                .agg(get_mode)
                .reindex(obs_pb.index)
            )

    # number of cells contributing to each pseudobulk
    obs_pb["n_cells_pb"] = (
        obs.groupby(groupby, observed=True)
        .size()
        .reindex(obs_pb.index)
        .astype(int)
    )

    # ------------------------------------------------------------
    # 5. Construct pseudobulk AnnData
    # ------------------------------------------------------------
    pb = AnnData(
        X=X_pb,
        obs=obs_pb,
        var=adata.var.copy(),
    )

    # Keep raw counts explicitly as well
    pb.layers["counts"] = pb.X.copy()

    return pb