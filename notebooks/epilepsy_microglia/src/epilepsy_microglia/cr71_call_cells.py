#!/usr/bin/env python3

"""
cr71_call_cells.py

Approximate reconstruction of Cell Ranger 7.1-era cell calling from a raw
10x MEX matrix:

    raw_feature_bc_matrix
        -> auto expect-cells estimation
        -> OrdMag initial cell calls
        -> ambient RNA model from barcode ranks 45k-90k
        -> Simple Good-Turing smoothing
        -> EmptyDrops-like multinomial testing
        -> BH FDR
        -> called barcode list

Designed for ONE library / donor at a time.

Input:
    matrix.mtx[.gz]
    barcodes.tsv[.gz]
    features.tsv[.gz]

Output:
    called_barcodes.tsv
    ordmag_barcodes.tsv
    candidate_stats.tsv.gz
    metrics.tsv

Dependencies:
    numpy
    scipy
    pandas
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from scipy.special import gammaln
from scipy import stats


# ---------------------------------------------------------------------
# Defaults intended to mimic Cell Ranger 7.1-era standard 3' GEX
# ---------------------------------------------------------------------

ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100
ORDMAG_QUANTILE = 0.99

MAX_EXPECTED_CELLS = 45_000

AMBIENT_RANK_LOW = 45_000
AMBIENT_RANK_HIGH = 90_000

EMPTYDROPS_MIN_UMI = 500
EMPTYDROPS_MIN_FRAC_MEDIAN = 0.01
EMPTYDROPS_MAX_CANDIDATES = 20_000

FDR_THRESHOLD = 0.01

# 10k is much faster for an initial reproduction.
# Increase to 100k once the behavior looks correct.
NUM_SIMS = 10_000

BOOTSTRAP_SEED = 0
MONTE_CARLO_SEED = 42


# =====================================================================
# I/O
# =====================================================================

def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def open_binary(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_barcodes(path: Path) -> np.ndarray:
    with open_text(path) as f:
        barcodes = np.array(
            [line.rstrip("\n").split("\t")[0] for line in f],
            dtype=object,
        )
    return barcodes


def read_features(path: Path) -> pd.DataFrame:
    rows = []

    with open_text(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            if len(fields) >= 3:
                rows.append((fields[0], fields[1], fields[2]))
            elif len(fields) == 2:
                rows.append((fields[0], fields[1], "Gene Expression"))
            else:
                rows.append((fields[0], fields[0], "Gene Expression"))

    return pd.DataFrame(
        rows,
        columns=["feature_id", "feature_name", "feature_type"],
    )


def read_mex(
    matrix_path: Path,
    barcodes_path: Path,
    features_path: Path,
):
    print(f"[load] matrix:   {matrix_path}")
    print(f"[load] barcodes: {barcodes_path}")
    print(f"[load] features: {features_path}")

    t0 = time.time()

    with open_binary(matrix_path) as f:
        X = mmread(f)

    X = X.tocsc()

    barcodes = read_barcodes(barcodes_path)
    features = read_features(features_path)

    # Standard 10x MEX orientation:
    # rows = features, columns = barcodes.
    if X.shape[1] != len(barcodes):
        if X.shape[0] == len(barcodes):
            print("[load] matrix appears transposed; transposing")
            X = X.T.tocsc()
        else:
            raise ValueError(
                f"Barcode count ({len(barcodes):,}) does not match "
                f"either matrix dimension {X.shape}"
            )

    if X.shape[0] != len(features):
        raise ValueError(
            f"Feature count ({len(features):,}) != matrix rows "
            f"({X.shape[0]:,})"
        )

    # Cell Ranger GEX cell calling uses Gene Expression features.
    if "feature_type" in features.columns:
        gex = features["feature_type"].eq("Gene Expression").to_numpy()

        if not np.all(gex):
            print(
                f"[load] retaining Gene Expression features: "
                f"{gex.sum():,} / {len(gex):,}"
            )
            X = X[gex, :].tocsc()
            features = features.loc[gex].reset_index(drop=True)

    # Raw Cell Ranger matrices should contain integer UMI counts.
    if X.nnz:
        sample = X.data[: min(1_000_000, X.nnz)]
        if not np.allclose(sample, np.round(sample)):
            raise ValueError(
                "Matrix contains non-integer values. "
                "This does not look like a raw UMI matrix."
            )

    X.data = np.asarray(np.round(X.data), dtype=np.int64)

    print(
        f"[load] loaded {X.shape[0]:,} genes × "
        f"{X.shape[1]:,} barcodes; nnz={X.nnz:,}"
    )
    print(f"[load] elapsed: {time.time() - t0:.1f}s")

    return X, barcodes, features


# =====================================================================
# OrdMag
# =====================================================================

def find_within_ordmag(
    counts: np.ndarray,
    baseline_idx,
):
    """
    Reproduce Cell Ranger's "within one order of magnitude" calculation.

    baseline_idx can be scalar or array.
    """
    x = np.sort(counts)

    baseline_idx = np.asarray(baseline_idx, dtype=int)

    baseline = x[-(baseline_idx + 1)]

    cutoff = np.maximum(
        1,
        np.round(0.1 * baseline),
    ).astype(int)

    # Number of elements >= cutoff.
    return len(x) - np.searchsorted(x, cutoff)


def estimate_recovered_cells_ordmag(
    nonzero_counts: np.ndarray,
    max_expected_cells: int = MAX_EXPECTED_CELLS,
):
    """
    Cell Ranger 7.x auto expect-cells approach.

    Search a log2-spaced grid and minimize:

        (OrdMag(x) - x)^2 / x
    """

    grid = np.linspace(
        1,
        np.log2(max_expected_cells),
        2000,
    )

    recovered = np.unique(
        np.round(2 ** grid).astype(int)
    )

    baseline_idx = np.round(
        recovered * (1 - ORDMAG_QUANTILE)
    ).astype(int)

    baseline_idx = np.minimum(
        baseline_idx,
        len(nonzero_counts) - 1,
    )

    called = find_within_ordmag(
        nonzero_counts,
        baseline_idx,
    )

    loss = (called - recovered) ** 2 / recovered

    best = np.argmin(loss)

    return int(recovered[best]), float(loss[best])


def summarize_bootstrapped_top_n(
    top_n_boot: np.ndarray,
    nonzero_counts: np.ndarray,
):
    """
    Approximate Cell Ranger summarize_bootstrapped_top_n().

    It uses the mean bootstrap call count, then attempts to extend the call
    boundary to include equal-count barcodes, unless doing so would add >20%.
    """

    mean_n = float(np.mean(top_n_boot))
    n = int(np.round(mean_n))

    if n <= 0:
        return 0

    sorted_counts = np.sort(nonzero_counts)[::-1]

    n = min(n, len(sorted_counts))

    cutoff = sorted_counts[n - 1]

    i = n - 1

    while (
        i + 1 < len(sorted_counts)
        and sorted_counts[i + 1] == cutoff
    ):
        i += 1

        extra = (i + 1) - n

        if extra > 0.20 * n:
            # Cell Ranger reverts to the original estimate.
            return n

    return i + 1


def call_ordmag(
    totals: np.ndarray,
    bootstrap_samples: int = ORDMAG_NUM_BOOTSTRAP_SAMPLES,
    max_expected_cells: int = MAX_EXPECTED_CELLS,
    seed: int = BOOTSTRAP_SEED,
):
    """
    Auto-estimate expected cells and perform OrdMag calling.
    """

    nonzero = totals[totals > 0].astype(np.int64)

    if len(nonzero) == 0:
        raise ValueError("No non-zero barcodes.")

    rs = np.random.RandomState(seed)

    print(
        f"[ordmag] estimating expect-cells with "
        f"{bootstrap_samples} bootstraps"
    )

    estimates = np.zeros(bootstrap_samples, dtype=float)
    losses = np.zeros(bootstrap_samples, dtype=float)

    for i in range(bootstrap_samples):
        sampled = rs.choice(
            nonzero,
            size=len(nonzero),
            replace=True,
        )

        estimates[i], losses[i] = estimate_recovered_cells_ordmag(
            sampled,
            max_expected_cells=max_expected_cells,
        )

        if (
            i == 0
            or (i + 1) % 10 == 0
            or i + 1 == bootstrap_samples
        ):
            print(
                f"[ordmag] expect bootstrap "
                f"{i + 1}/{bootstrap_samples}: "
                f"{int(estimates[i]):,}"
            )

    recovered_cells = max(
        int(np.round(np.mean(estimates))),
        50,
    )

    print(
        f"[ordmag] auto expect-cells = "
        f"{recovered_cells:,}"
    )
    print(
        f"[ordmag] mean estimation loss = "
        f"{np.mean(losses):.4f}"
    )

    baseline_idx = int(
        np.round(
            recovered_cells *
            (1 - ORDMAG_QUANTILE)
        )
    )

    baseline_idx = min(
        baseline_idx,
        len(nonzero) - 1,
    )

    print(
        f"[ordmag] baseline rank index = "
        f"{baseline_idx}"
    )

    top_n_boot = np.zeros(
        bootstrap_samples,
        dtype=int,
    )

    for i in range(bootstrap_samples):
        sampled = rs.choice(
            nonzero,
            size=len(nonzero),
            replace=True,
        )

        top_n_boot[i] = find_within_ordmag(
            sampled,
            baseline_idx,
        )

        if (
            i == 0
            or (i + 1) % 10 == 0
            or i + 1 == bootstrap_samples
        ):
            print(
                f"[ordmag] calling bootstrap "
                f"{i + 1}/{bootstrap_samples}: "
                f"{top_n_boot[i]:,}"
            )

    top_n = summarize_bootstrapped_top_n(
        top_n_boot,
        nonzero,
    )

    # Stable ordering for reproducibility.
    order = np.argsort(
        totals,
        kind="stable",
    )[::-1]

    initial_idx = np.sort(order[:top_n])

    initial_counts = totals[initial_idx]

    print(f"[ordmag] initial cells = {len(initial_idx):,}")
    print(
        f"[ordmag] initial UMI range = "
        f"{initial_counts.min():,} - "
        f"{initial_counts.max():,}"
    )
    print(
        f"[ordmag] median initial UMI = "
        f"{np.median(initial_counts):.1f}"
    )

    return {
        "initial_idx": initial_idx,
        "recovered_cells": recovered_cells,
        "top_n": top_n,
        "bootstrap_estimates": estimates,
        "bootstrap_top_n": top_n_boot,
    }


# =====================================================================
# Simple Good-Turing
# Ported from the public Cell Ranger implementation
# =====================================================================

class SimpleGoodTuringError(Exception):
    pass


def _averaging_transform(r, nr):
    d = np.concatenate(
        (
            np.ones(1, dtype=int),
            np.diff(r),
        )
    )

    dr = np.concatenate(
        (
            0.5 * (d[1:] + d[:-1]),
            np.array((d[-1],), dtype=float),
        )
    )

    return nr.astype(float) / dr


def _rstest(r, coef):
    return r * np.power(
        1 + 1.0 / r,
        1 + coef,
    )


def simple_good_turing(xr, xnr):
    xr = xr.astype(float)
    xnr = xnr.astype(float)

    xN = np.sum(xr * xnr)

    xnrz = _averaging_transform(
        xr,
        xnr,
    )

    slope, _, _, _, _ = stats.linregress(
        np.log(xr),
        np.log(xnrz),
    )

    if slope >= -1:
        print(
            f"[sgt] slope={slope:.4f} >= -1; "
            "falling back to slope=-1"
        )
        slope = -1.0
    else:
        print(f"[sgt] slope={slope:.4f}")

    xrst = _rstest(
        xr,
        slope,
    )

    xrstrel = xrst / xr

    xrtry = xr == np.concatenate(
        (
            xr[1:] - 1,
            np.zeros(1),
        )
    )

    xrstarel = np.zeros(
        len(xr),
        dtype=float,
    )

    shifted_nr = np.concatenate(
        (
            xnr[1:],
            np.zeros(1),
        )
    )

    xrstarel[xrtry] = (
        (xr[xrtry] + 1)
        / xr[xrtry]
        * shifted_nr[xrtry]
        / xnr[xrtry]
    )

    tursd = np.ones(
        len(xr),
        dtype=float,
    )

    for i in range(len(xr)):
        if xrtry[i]:
            tursd[i] = (
                float(i + 2)
                / xnr[i]
                * np.sqrt(
                    xnr[i + 1]
                    * (
                        1
                        + xnr[i + 1]
                        / xnr[i]
                    )
                )
            )

    xrstcmbrel = np.zeros(
        len(xr),
        dtype=float,
    )

    use_turing = True

    for i in range(len(xr)):
        if not use_turing:
            xrstcmbrel[i] = xrstrel[i]

        elif (
            np.abs(
                xrstrel[i]
                - xrstarel[i]
            )
            * (1 + i)
            / tursd[i]
            > 1.65
        ):
            xrstcmbrel[i] = xrstarel[i]

        else:
            use_turing = False
            xrstcmbrel[i] = xrstrel[i]

    sumpraw = np.sum(
        xrstcmbrel
        * xr
        * xnr
        / xN
    )

    p0 = xnr[0] / xN

    xrstcmbrel = (
        xrstcmbrel
        * (1 - p0)
        / sumpraw
    )

    return xr * xrstcmbrel, p0


def sgt_proportions(frequencies):
    frequencies = np.asarray(
        frequencies,
        dtype=np.int64,
    )

    if len(frequencies) == 0:
        raise ValueError(
            "Input frequency vector is empty."
        )

    if np.any(frequencies <= 0):
        raise ValueError(
            "SGT frequencies must all be > 0."
        )

    freqfreqs = np.bincount(frequencies)

    use_freqs = np.flatnonzero(freqfreqs)

    if len(use_freqs) < 10:
        raise SimpleGoodTuringError(
            f"Too few non-zero frequency classes: "
            f"{len(use_freqs)}"
        )

    rstar, p0 = simple_good_turing(
        use_freqs,
        freqfreqs[use_freqs],
    )

    mapping = dict(
        zip(
            use_freqs,
            rstar,
        )
    )

    rstar_sum = np.sum(
        freqfreqs[use_freqs]
        * rstar
    )

    rstar_i = np.fromiter(
        (
            mapping[f]
            for f in frequencies
        ),
        dtype=float,
        count=len(frequencies),
    )

    pstar = (
        (1 - p0)
        * rstar_i
        / rstar_sum
    )

    if not np.isclose(
        p0 + pstar.sum(),
        1.0,
    ):
        raise RuntimeError(
            "SGT probabilities do not sum to 1."
        )

    return pstar, p0


def estimate_ambient_profile_sgt(
    X: sp.csc_matrix,
    ambient_idx: np.ndarray,
):
    """
    Estimate ambient gene probabilities using SGT.

    Genes with zero ambient counts receive an equal share of p0.
    """

    ambient_counts = np.asarray(
        X[:, ambient_idx].sum(axis=1)
    ).ravel().astype(np.int64)

    zero_features = np.flatnonzero(
        ambient_counts == 0
    )

    nonzero_features = np.flatnonzero(
        ambient_counts > 0
    )

    if len(nonzero_features) == 0:
        raise ValueError(
            "Ambient barcode range contains no counts."
        )

    smoothed, p0 = sgt_proportions(
        ambient_counts[nonzero_features]
    )

    p = np.empty(
        X.shape[0],
        dtype=float,
    )

    if len(zero_features) == 0:
        smoothed /= smoothed.sum()

        p[:] = 0.0
        p[nonzero_features] = smoothed

    else:
        p0_each = p0 / len(zero_features)

        p[:] = p0_each
        p[nonzero_features] = smoothed

    if not np.isclose(
        p.sum(),
        1.0,
    ):
        raise RuntimeError(
            f"Ambient profile sums to {p.sum()}"
        )

    print(
        f"[ambient] genes observed in ambient pool = "
        f"{len(nonzero_features):,}"
    )
    print(
        f"[ambient] genes unobserved in ambient pool = "
        f"{len(zero_features):,}"
    )
    print(
        f"[ambient] SGT unobserved probability mass = "
        f"{p0:.6g}"
    )

    return p


# =====================================================================
# Multinomial likelihood / simulation
# =====================================================================

def eval_multinomial_loglikelihoods(
    X: sp.csc_matrix,
    logp: np.ndarray,
    totals: np.ndarray,
):
    """
    Multinomial log-likelihood for every column in X.
    """

    X = X.tocsc()

    result = np.zeros(
        X.shape[1],
        dtype=float,
    )

    const = gammaln(
        totals + 1
    )

    for j in range(X.shape[1]):
        start = X.indptr[j]
        end = X.indptr[j + 1]

        gene_idx = X.indices[start:end]
        counts = X.data[start:end]

        result[j] = (
            const[j]
            - gammaln(
                counts + 1
            ).sum()
            + (
                counts
                * logp[gene_idx]
            ).sum()
        )

    return result


def incremental_counts_from_sample(
    draws: np.ndarray,
):
    """
    For each categorical draw, return how many times that category
    has appeared up to that point.
    """

    n = len(draws)

    if n == 0:
        return np.empty(
            0,
            dtype=draws.dtype,
        )

    order = np.argsort(
        draws,
        kind="stable",
    )

    sorted_vals = draws[order]

    arange_n = np.arange(n)

    is_new = np.empty(
        n,
        dtype=bool,
    )

    is_new[0] = True

    np.not_equal(
        sorted_vals[1:],
        sorted_vals[:-1],
        out=is_new[1:],
    )

    group_start = np.maximum.accumulate(
        np.where(
            is_new,
            arange_n,
            -1,
        )
    )

    rank_within = (
        arange_n
        - group_start
        + 1
    )

    result = np.empty(
        n,
        dtype=draws.dtype,
    )

    result[order] = rank_within

    return result


def cumulative_multinomial_loglikelihood(
    draws: np.ndarray,
    logp: np.ndarray,
    log_nvals: np.ndarray,
):
    marginal_counts = (
        incremental_counts_from_sample(draws)
    )

    increments = (
        log_nvals[:len(draws)]
        - np.log(marginal_counts)
        + logp[draws]
    )

    return np.cumsum(increments)


def compute_ambient_pvalues(
    ambient_p: np.ndarray,
    candidate_totals: np.ndarray,
    observed_loglk: np.ndarray,
    num_sims: int = NUM_SIMS,
    seed: int = MONTE_CARLO_SEED,
):
    """
    Monte Carlo p-values using the same cumulative-simulation idea used
    by Cell Ranger's public implementation.

    One categorical sequence is drawn per simulation. Its cumulative
    multinomial likelihood provides simulated likelihoods at all
    candidate UMI totals simultaneously.
    """

    candidate_totals = np.asarray(
        candidate_totals,
        dtype=int,
    )

    distinct_n = np.unique(
        candidate_totals
    )

    max_n = int(
        candidate_totals.max()
    )

    print(
        f"[emptydrops] candidate UMI range = "
        f"{candidate_totals.min():,} - "
        f"{candidate_totals.max():,}"
    )
    print(
        f"[emptydrops] distinct UMI totals = "
        f"{len(distinct_n):,}"
    )
    print(
        f"[emptydrops] simulations = "
        f"{num_sims:,}"
    )

    p_cumulative = np.cumsum(
        ambient_p
    )

    # Protect against tiny numerical error.
    p_cumulative[-1] = 1.0

    logp = np.log(
        ambient_p
    )

    n_lookup = np.searchsorted(
        distinct_n,
        candidate_totals,
    )

    lower_count = np.zeros(
        len(candidate_totals),
        dtype=np.int64,
    )

    log_nvals = np.log(
        np.arange(
            1,
            max_n + 1,
        )
    )

    rng = np.random.default_rng(
        seed
    )

    for sim in range(num_sims):

        u = rng.random(
            max_n
        )

        draws = np.searchsorted(
            p_cumulative,
            u,
        )

        cumulative_ll = (
            cumulative_multinomial_loglikelihood(
                draws,
                logp,
                log_nvals,
            )
        )

        sim_ll = cumulative_ll[
            distinct_n - 1
        ]

        lower_count += (
            sim_ll[n_lookup]
            < observed_loglk
        )

        if (
            sim == 0
            or (sim + 1) % 1000 == 0
            or sim + 1 == num_sims
        ):
            print(
                f"[emptydrops] simulation "
                f"{sim + 1:,}/{num_sims:,}"
            )

    pvalues = (
        1 + lower_count
    ) / (
        1 + num_sims
    )

    return pvalues


def bh_adjust(p):
    """
    Benjamini-Hochberg FDR adjustment.
    """

    p = np.asarray(
        p,
        dtype=float,
    )

    n = len(p)

    order = np.argsort(
        p
    )

    ranked = (
        p[order]
        * n
        / np.arange(
            1,
            n + 1,
        )
    )

    adjusted_sorted = np.minimum.accumulate(
        ranked[::-1]
    )[::-1]

    adjusted_sorted = np.minimum(
        adjusted_sorted,
        1.0,
    )

    adjusted = np.empty(
        n,
        dtype=float,
    )

    adjusted[order] = adjusted_sorted

    return adjusted


# =====================================================================
# EmptyDrops / nonambient stage
# =====================================================================

def call_nonambient(
    X: sp.csc_matrix,
    totals: np.ndarray,
    initial_idx: np.ndarray,
    ambient_rank_low: int = AMBIENT_RANK_LOW,
    ambient_rank_high: int = AMBIENT_RANK_HIGH,
    min_umi: int = EMPTYDROPS_MIN_UMI,
    min_frac_median: float = EMPTYDROPS_MIN_FRAC_MEDIAN,
    max_candidates: int = EMPTYDROPS_MAX_CANDIDATES,
    fdr_threshold: float = FDR_THRESHOLD,
    num_sims: int = NUM_SIMS,
    seed: int = MONTE_CARLO_SEED,
):
    order = np.argsort(
        totals,
        kind="stable",
    )[::-1]

    n_bcs = len(totals)

    low = min(
        ambient_rank_low,
        n_bcs,
    )

    high = min(
        ambient_rank_high,
        n_bcs,
    )

    if low >= high:
        raise ValueError(
            f"Cannot construct ambient range "
            f"{low}:{high} from {n_bcs} barcodes."
        )

    # Python slice [45000:90000] corresponds approximately to
    # Cell Ranger's rank 45k-90k ambient range.
    ambient_idx = np.sort(
        order[low:high]
    )

    ambient_umis = totals[
        ambient_idx
    ]

    print(
        f"[ambient] barcode rank range = "
        f"{low:,} - {high:,}"
    )
    print(
        f"[ambient] ambient UMI range = "
        f"{ambient_umis.min():,} - "
        f"{ambient_umis.max():,}"
    )

    ambient_p = estimate_ambient_profile_sgt(
        X,
        ambient_idx,
    )

    # Historical CR7-era candidate threshold:
    #
    # max(
    #     500,
    #     0.01 * median UMI of OrdMag cells
    # )
    #
    median_initial = float(
        np.median(
            totals[initial_idx]
        )
    )

    candidate_min_umi = max(
        int(min_umi),
        int(
            np.round(
                min_frac_median
                * median_initial
            )
        ),
    )

    print(
        f"[emptydrops] median initial UMI = "
        f"{median_initial:.1f}"
    )
    print(
        f"[emptydrops] candidate minimum UMI = "
        f"{candidate_min_umi:,}"
    )

    initial_mask = np.zeros(
        n_bcs,
        dtype=bool,
    )

    initial_mask[
        initial_idx
    ] = True

    # Candidate pool:
    #   not already called by OrdMag
    #   UMI > historical candidate threshold
    #   at most next 20k barcodes
    remaining_order = order[
        ~initial_mask[order]
    ]

    remaining_order = remaining_order[
        totals[remaining_order]
        > candidate_min_umi
    ]

    candidate_idx = remaining_order[
        :max_candidates
    ]

    # Ambient range should normally be far below the candidates.
    # Explicitly remove overlap anyway.
    ambient_set = np.zeros(
        n_bcs,
        dtype=bool,
    )

    ambient_set[
        ambient_idx
    ] = True

    candidate_idx = candidate_idx[
        ~ambient_set[
            candidate_idx
        ]
    ]

    candidate_idx = np.sort(
        candidate_idx
    )

    print(
        f"[emptydrops] candidate barcodes = "
        f"{len(candidate_idx):,}"
    )

    if len(candidate_idx) == 0:
        return {
            "candidate_idx": candidate_idx,
            "pvalue": np.array([]),
            "fdr": np.array([]),
            "nonambient": np.array([], dtype=bool),
            "ambient_idx": ambient_idx,
            "candidate_min_umi": candidate_min_umi,
        }

    Xcand = X[
        :,
        candidate_idx
    ].tocsc()

    candidate_totals = totals[
        candidate_idx
    ].astype(int)

    logp = np.log(
        ambient_p
    )

    print(
        "[emptydrops] computing observed "
        "multinomial likelihoods"
    )

    observed_ll = eval_multinomial_loglikelihoods(
        Xcand,
        logp,
        candidate_totals,
    )

    print(
        "[emptydrops] starting Monte Carlo test"
    )

    pvalues = compute_ambient_pvalues(
        ambient_p=ambient_p,
        candidate_totals=candidate_totals,
        observed_loglk=observed_ll,
        num_sims=num_sims,
        seed=seed,
    )

    fdr = bh_adjust(
        pvalues
    )

    nonambient = (
        fdr <= fdr_threshold
    )

    print(
        f"[emptydrops] rescued cells = "
        f"{nonambient.sum():,}"
    )

    return {
        "candidate_idx": candidate_idx,
        "observed_loglk": observed_ll,
        "pvalue": pvalues,
        "fdr": fdr,
        "nonambient": nonambient,
        "ambient_idx": ambient_idx,
        "candidate_min_umi": candidate_min_umi,
    }


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--matrix",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--barcodes",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--features",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--sample",
        required=True,
    )

    parser.add_argument(
        "--ordmag-only",
        action="store_true",
        help=(
            "Stop after OrdMag. "
            "Useful for first sanity check."
        ),
    )

    parser.add_argument(
        "--bootstrap-samples",
        type=int,
        default=ORDMAG_NUM_BOOTSTRAP_SAMPLES,
    )

    parser.add_argument(
        "--max-expected-cells",
        type=int,
        default=MAX_EXPECTED_CELLS,
    )

    parser.add_argument(
        "--ambient-rank-low",
        type=int,
        default=AMBIENT_RANK_LOW,
    )

    parser.add_argument(
        "--ambient-rank-high",
        type=int,
        default=AMBIENT_RANK_HIGH,
    )

    parser.add_argument(
        "--min-umi",
        type=int,
        default=EMPTYDROPS_MIN_UMI,
    )

    parser.add_argument(
        "--min-frac-median",
        type=float,
        default=EMPTYDROPS_MIN_FRAC_MEDIAN,
    )

    parser.add_argument(
        "--max-candidates",
        type=int,
        default=EMPTYDROPS_MAX_CANDIDATES,
    )

    parser.add_argument(
        "--fdr",
        type=float,
        default=FDR_THRESHOLD,
    )

    parser.add_argument(
        "--num-sims",
        type=int,
        default=NUM_SIMS,
    )

    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=BOOTSTRAP_SEED,
    )

    parser.add_argument(
        "--mc-seed",
        type=int,
        default=MONTE_CARLO_SEED,
    )

    args = parser.parse_args()

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    t0 = time.time()

    X, barcodes, features = read_mex(
        args.matrix,
        args.barcodes,
        args.features,
    )

    print("[counts] calculating UMI totals")

    totals = np.asarray(
        X.sum(axis=0)
    ).ravel().astype(np.int64)

    print(
        f"[counts] non-zero barcodes = "
        f"{np.sum(totals > 0):,} / "
        f"{len(totals):,}"
    )

    print(
        f"[counts] maximum UMI = "
        f"{totals.max():,}"
    )

    # --------------------------------------------------------------
    # OrdMag
    # --------------------------------------------------------------

    ordmag = call_ordmag(
        totals,
        bootstrap_samples=args.bootstrap_samples,
        max_expected_cells=args.max_expected_cells,
        seed=args.bootstrap_seed,
    )

    initial_idx = ordmag[
        "initial_idx"
    ]

    ordmag_barcodes = barcodes[
        initial_idx
    ]

    pd.Series(
        ordmag_barcodes
    ).to_csv(
        args.outdir / "ordmag_barcodes.tsv",
        index=False,
        header=False,
    )

    if args.ordmag_only:
        metrics = pd.DataFrame(
            [{
                "sample": args.sample,
                "raw_barcodes": len(barcodes),
                "nonzero_barcodes": int(
                    np.sum(totals > 0)
                ),
                "auto_expect_cells": ordmag[
                    "recovered_cells"
                ],
                "ordmag_cells": len(initial_idx),
                "median_ordmag_umi": float(
                    np.median(
                        totals[
                            initial_idx
                        ]
                    )
                ),
            }]
        )

        metrics.to_csv(
            args.outdir / "metrics.tsv",
            sep="\t",
            index=False,
        )

        print(
            f"[done] OrdMag only. "
            f"Elapsed: {(time.time() - t0)/60:.1f} min"
        )

        return

    # --------------------------------------------------------------
    # EmptyDrops-style nonambient rescue
    # --------------------------------------------------------------

    ed = call_nonambient(
        X=X,
        totals=totals,
        initial_idx=initial_idx,
        ambient_rank_low=args.ambient_rank_low,
        ambient_rank_high=args.ambient_rank_high,
        min_umi=args.min_umi,
        min_frac_median=args.min_frac_median,
        max_candidates=args.max_candidates,
        fdr_threshold=args.fdr,
        num_sims=args.num_sims,
        seed=args.mc_seed,
    )

    candidate_idx = ed[
        "candidate_idx"
    ]

    rescued_idx = candidate_idx[
        ed["nonambient"]
    ]

    called_idx = np.sort(
        np.unique(
            np.concatenate(
                (
                    initial_idx,
                    rescued_idx,
                )
            )
        )
    )

    # --------------------------------------------------------------
    # Outputs
    # --------------------------------------------------------------

    pd.Series(
        barcodes[
            called_idx
        ]
    ).to_csv(
        args.outdir / "called_barcodes.tsv",
        index=False,
        header=False,
    )

    if len(candidate_idx):
        candidate_df = pd.DataFrame(
            {
                "barcode": barcodes[
                    candidate_idx
                ],
                "total_umi": totals[
                    candidate_idx
                ],
                "log_likelihood": ed[
                    "observed_loglk"
                ],
                "pvalue": ed[
                    "pvalue"
                ],
                "fdr": ed[
                    "fdr"
                ],
                "is_nonambient": ed[
                    "nonambient"
                ],
            }
        )

        candidate_df.to_csv(
            args.outdir
            / "candidate_stats.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    ambient_umis = totals[
        ed["ambient_idx"]
    ]

    metrics = pd.DataFrame(
        [{
            "sample": args.sample,

            "raw_genes": X.shape[0],
            "raw_barcodes": X.shape[1],
            "nonzero_barcodes": int(
                np.sum(totals > 0)
            ),

            "auto_expect_cells": ordmag[
                "recovered_cells"
            ],

            "ordmag_cells": len(
                initial_idx
            ),

            "median_ordmag_umi": float(
                np.median(
                    totals[
                        initial_idx
                    ]
                )
            ),

            "ambient_rank_low": args.ambient_rank_low,
            "ambient_rank_high": args.ambient_rank_high,

            "ambient_max_umi": int(
                ambient_umis.max()
            ),

            "ambient_min_umi": int(
                ambient_umis.min()
            ),

            "candidate_min_umi": ed[
                "candidate_min_umi"
            ],

            "candidate_cells_tested": len(
                candidate_idx
            ),

            "emptydrops_rescued": int(
                ed["nonambient"].sum()
            ),

            "fdr_threshold": args.fdr,

            "num_sims": args.num_sims,

            "final_called_cells": len(
                called_idx
            ),
        }]
    )

    metrics.to_csv(
        args.outdir / "metrics.tsv",
        sep="\t",
        index=False,
    )

    print("")
    print("=" * 60)
    print(f"Sample:             {args.sample}")
    print(
        f"Raw barcodes:       "
        f"{len(barcodes):,}"
    )
    print(
        f"Auto expect-cells:  "
        f"{ordmag['recovered_cells']:,}"
    )
    print(
        f"OrdMag cells:       "
        f"{len(initial_idx):,}"
    )
    print(
        f"Candidates tested:  "
        f"{len(candidate_idx):,}"
    )
    print(
        f"EmptyDrops rescued: "
        f"{len(rescued_idx):,}"
    )
    print(
        f"FINAL CALLED:       "
        f"{len(called_idx):,}"
    )
    print("=" * 60)

    print(
        f"[done] elapsed: "
        f"{(time.time() - t0)/60:.1f} min"
    )


if __name__ == "__main__":
    main()