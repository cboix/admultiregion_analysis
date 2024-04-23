#!/usr/bin/python
"""Utility scripts for correlation class."""

import logging
import numpy as np
from scipy.stats import norm
from scipy import sparse


def triu_mask(A):
    """Get the upper triangular part a matrix."""
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:, None] < r
    return A[mask]


def calculate_correlation(X, center=False):
    """Calculate the correlation of columns of a matrix."""
    logging.info("Calculate correlation of columns of matrix X")
    # TODO: CENTER (+ OTHER WAYS of calc)
    # TODO: Allow numpy corrcoeff
    if center:
        # TODO: See difference with centering X:
        logging.debug("Raw centering not implemented")
        corr = None
    else:
        xTx = X.T.dot(X) / X.shape[0]
        if sparse.issparse(xTx):
            xTx = xTx.toarray()  # TODO: unless x is not sparse
        sd = np.sqrt(np.diag(xTx))
        xTx = xTx / sd[:, np.newaxis]
        corr = xTx / sd[np.newaxis, :]
        del(xTx, sd)
    return(corr)


def calculate_correlation_estimate(U, s, V, power=0, indices=None):
    if type(power) is list:
        corr_est = np.zeros((V.shape[1], V.shape[1]))
        j = 0
        for pw in power:
            j += 1
            corr_est += calculate_single_correlation_estimate(
                U=U, s=s, V=V, power=pw, indices=indices)
        corr_est = corr_est / (j * 1.0)
    else:
        corr_est = calculate_single_correlation_estimate(
            U=U, s=s, V=V, power=power, indices=indices)
    return(corr_est)


def calculate_single_correlation_estimate(U, s, V, power=0, indices=None):
    """Calculate an SVD-derived estimate of the correlation matrix."""
    logging.info("Estimating correlation of columns of matrix X with its SVD")
    # TODO: get proper scaling factor here (d / N ?)
    scale = U.shape[0] * V.shape[1]
    if indices is not None:
        logging.debug("Reducing to subset of indices")
        # TODO: Proper handling of subset.
        # TODO: Return non-scaled
        s = s[indices]
        V = V[indices, :]

    if power != 0:
        logging.debug(f"Using power {power}")
        # Allow any range of s-power, 0.5 is
        # For comparison with Vs^2V.T (approx cov.)
        s_red = np.diag(s**power)
        X_cov = V.T.dot(s_red).dot(V) / scale
    else:
        X_cov = V.T.dot(V) / scale

    X_sd = np.sqrt(np.diag(X_cov))
    cv = X_cov / X_sd[:, np.newaxis]
    corr_est = cv / X_sd[np.newaxis, :]
    del (cv, X_cov, X_sd, scale)
    return(corr_est)


# TODO: This is very expensive, update to compute as we go.
def calculate_correlation_estimate_sd(U, s, V, nperm=25, seed=1):
    """Estimate the sd of correlations by subsampling V."""
    # Estimate the sd of correlations by subsampling V:
    np.random.seed(seed)  # Reproducible sampling
    k = U.shape[1]
    full_ind = np.random.randint(k, size=(k // 2, nperm))
    corr_est = np.zeros((V.shape[1], V.shape[1], nperm))
    # TODO: Perform running computation for reduced memory footprint
    for i in range(nperm):
        logging.debug(i)
        corr_est[:, :, i] = calculate_correlation_estimate(
            U, s, V, indices=full_ind[:, i])

    corr_mean = np.mean(corr_est, axis=-1)
    corr_sd = np.std(corr_est, axis=-1)
    del(corr_est, full_ind)
    return (corr_mean, corr_sd)


def search_FDR_cutoff(t_ij, P, alpha=0.01):
    """Grid search for an FDR cutoff given estimated test statistics."""
    vec = triu_mask(t_ij)  # Upper triangular matrix as vector
    # Valid search space (for which t_ij ~ N(0,1))
    ap = 2 * np.log(np.log(P))
    bp = np.sqrt(4 * np.log(P) - ap)
    res = 0.1
    t = np.array(list(np.arange(2, bp, res)) + [bp])

    def grid_FDR_cutoff(vec, t, alpha):
        est_null = (2 - 2 * norm.cdf(t)) * len(vec)
        sig = np.array([np.sum(vec > x) for x in t])
        efdr = est_null / sig
        t_ind = np.where(efdr < alpha)[0][0]
        return(t_ind, t[t_ind - 1], t[t_ind])

    # Run a hierarchical grid search:
    while res >= 0.01:
        t_ind, t_low, t_high = grid_FDR_cutoff(vec, t, alpha)
        res = res * 0.1
        t = np.arange(t_low, t_high, res)

    return(t_high)
