#!/usr/bin/python
"""Utility scripts for gene graph class."""

import logging
import numpy as np
import pandas as pd
from scipy import sparse
import scipy.stats as st
from scipy.interpolate import SmoothBivariateSpline


def prune_degree(aw, degree_cutoff):
    """Remove nodes with degree of `degree_cutoff` or lower."""
    logging.debug("Removing nodes with degree <=" + str(degree_cutoff))
    aw = aw.tocsr()
    rind = np.array((np.sum(aw > 0, axis=1) > degree_cutoff).T)[0]
    cind = np.array((np.sum(aw > 0, axis=0) > degree_cutoff))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return (aw, ind)


# External / auxiliary functions for graph operations:
def prune_scale(aw, scale=0.90):
    """Prune graph by a relative scaling value."""
    logging.debug("Removing edges below " + str(round(scale*100, 2)) +
                 "\\% of max edge for each node.")
    # TODO: handle 0s for max:
    awm1 = np.max(aw, axis=1).T.toarray()[0]
    awm2 = np.max(aw, axis=0).toarray()[0]
    awm = np.vstack((awm1, awm2))
    awm = np.max(awm, axis=0)
    aw = aw.tocoo()
    scl = np.vstack((awm[aw.row], awm[aw.col]))
    scl = np.max(scl, axis=0)
    pct_data = aw.data / scl
    kind = pct_data > scale
    aw = sparse.coo_matrix(
        (aw.data[kind], (aw.row[kind], aw.col[kind])), shape=aw.shape)
    return aw


def prune_knn(aw, k=50, twodir=True, row_only=False):
    """Prune graph to maximum k for each node."""
    logging.debug("Pruning graph to maximum k-edges for each node")
    aw = aw.tocoo()
    ind = np.argsort(-aw.data)  # Sort correlations
    mk = np.zeros(aw.shape[0], int)  # Kept margin
    keepind = np.zeros(len(ind), int)
    # Add indices in order of strength:
    for i in range(len(ind)):
        r = aw.row[i]
        c = aw.col[i]
        if twodir:
            cond = mk[r] < k or mk[c] < k
        else:
            cond = mk[r] < k and mk[c] < k
        if cond:
            mk[r] += 1
            if twodir or not row_only:
                mk[c] += 1
            keepind[i] = 1
    # Remove edges:
    kind = ind[keepind == 1]
    aw = sparse.coo_matrix((aw.data[kind], (aw.row[kind], aw.col[kind])),
                           shape=aw.shape)
    return aw


def set_zscore(p, n):
    """Set a z-score for p-value p in n tests."""
    # TODO: get proper BY testing correction equation
    z = -st.norm.ppf(p / n)  # If z is none.
    return z


def bin_genes(gmarg, gmarg_cut=0.75):
    if gmarg_cut is not None:
        gmarg[gmarg > gmarg_cut] = gmarg_cut
    lmarg = np.log10(np.array(gmarg))
    lmarg[lmarg < -3] = -3
    # Different ways of determining cutoffs for digitizing:
    # xs = np.linspace(np.min([np.min(lmarg), -3]), 0, 50)
    # xs = np.linspace(np.min([np.min(lmarg), -3]), np.max(lmarg), 50)
    xs = np.linspace(np.min(lmarg), np.max(lmarg), 30)
    dmarg = np.digitize(lmarg, xs)
    bins = np.unique(dmarg)
    return(dmarg, xs, bins)


def make_outliers_zero(x, n_top=25, pct_top=0.01):
    nx = np.prod(x.shape)
    if n_top < (nx * pct_top):
        pct_top = n_top / nx
    qt = np.quantile(x, 1 - pct_top)
    x[x >= qt] = 0
    return(x)


def get_binned_stats(corr, dmarg, bins, zero_outliers=True):
    # Initialize statistics:
    nbin = len(bins)
    smean = np.zeros((nbin, nbin))
    ssd = np.zeros((nbin, nbin))
    scount = np.zeros((nbin, nbin))
    # Populate statistics:
    for i, xi in enumerate(bins):
        rind = (dmarg == xi)
        ci = corr[rind, :]
        for j, xj in enumerate(bins):
            cind = (dmarg == xj)
            cij = ci[:, cind].copy()  # Copy for if we alter by zeroing.
            if zero_outliers:
                # Much faster than remove, but shouldn't affect statistics
                cij = make_outliers_zero(cij, n_top=25)
            smean[i, j] = np.mean(cij)
            ssd[i, j] = np.std(cij)
            # scount[i, j] = len(cij)
            scount[i, j] = cij.shape[0] * cij.shape[1]
    # Turn into sparse matrices with same datasize:
    scount = sparse.coo_matrix(np.log(scount + 2) * (ssd > 0))
    smean = sparse.coo_matrix(smean * (ssd > 0))
    ssd = sparse.coo_matrix(ssd)
    return(smean, ssd, scount)


def smooth_binned_stats(smean, ssd, scount, xs, bins, min_sd=0.01):
    logging.info("Fitting + predicting spline")
    overall_mean = np.mean(smean.data)
    overall_sd = np.mean(ssd.data)
    logging.debug("mean and sd: " + str(overall_mean) + ", " + str(overall_sd))
    ks = 2
    splmean = SmoothBivariateSpline(
        x=bins[smean.row], y=bins[smean.col], z=smean.data - overall_mean,
        kx=ks, ky=ks, w=scount.data,
    )
    splsd = SmoothBivariateSpline(
        x=bins[ssd.row], y=bins[ssd.col], z=ssd.data - overall_sd,
        kx=ks, ky=ks, w=scount.data
    )
    xloc = np.arange(len(xs))
    smean = splmean(xloc, xloc)
    smean = smean + overall_mean
    ssd = splsd(xloc, xloc)
    ssd = ssd + overall_sd
    ssd[ssd < min_sd] = min_sd
    return(smean, ssd, xloc)


def get_binned_cutoffs(corr, gmarg, gmarg_cut=0.75, min_sd=0.01,
                       zero_outliers=True):
    # Bin genes by their sparsity:
    dmarg, xs, bins = bin_genes(gmarg, gmarg_cut=gmarg_cut)
    # Calculate the means / sdevs / counts per combination of bins:
    smean, ssd, scount = get_binned_stats(corr, dmarg, bins,
                                          zero_outliers=zero_outliers)
    # Fit the spline to the bivariate data:
    smean, ssd, xloc = smooth_binned_stats(
        smean, ssd, scount, xs, bins, min_sd=min_sd)
    return(smean, ssd, xloc, dmarg)


def adj_from_bivariate_cutoff(corr, gmarg, z=None, min_sd=0.01,
                              zero_outliers=True):
    """Get thresholded adjacency with a bivariate spline cutoff."""
    # NOTE: Original function, switched to zscore
    # Set the z-score:
    z = set_zscore(p=0.01, n=corr.shape[0]) if z is None else z
    logging.info("Thresholding by sparsity with estimated mean / variance " +
                 "bivariate splines, with z=" + str(np.round(z, 2)))

    # Get the smoothed cutoffs:
    corr = corr - np.diag(np.diag(corr))  # Remove diagonal for adjacency
    smean, ssd, xloc, dmarg = get_binned_cutoffs(
        corr, gmarg, gmarg_cut=0.75, min_sd=min_sd, zero_outliers=zero_outliers)

    # Calculate the z-scored cutoff based on spline predictions:
    zcut = smean + ssd * z  # Correlation cutoff given the zscore
    zcoo = sparse.coo_matrix(zcut)
    zdf = pd.DataFrame({"x": xloc[zcoo.row] + 1,
                        "y": xloc[zcoo.col] + 1,
                        "cutoff": zcoo.data})

    # Merge cutoffs against pre-cut:
    logging.info("Thresholding adjacency matrix")
    min_cut = np.min(zcoo.data)
    logging.debug("Minimum correlation: " + str(round(min_cut, 2)))
    zcorr = corr * (corr >= min_cut)
    aw = sparse.coo_matrix(zcorr)  # TODO: Bypass this creation
    # aind = aw.data >= min_cut
    awcut_df = pd.DataFrame({"x": dmarg[aw.row], "y": dmarg[aw.col],
                             "r": aw.row, "c": aw.col, "dat": aw.data})
    awcut_df = awcut_df.merge(zdf)
    awcut_df = awcut_df[awcut_df.dat >= awcut_df.cutoff]

    # Create thresholded adjacency matrix from rows, cols, data:
    aw = sparse.coo_matrix(
        (awcut_df.dat, (awcut_df.r, awcut_df.c)), shape=aw.shape)
    return (aw, zcut)


def zscore_from_bivariate_cutoff(corr, gmarg, z=None, min_sd=0.01,
                                 zero_outliers=True, keep_all_z=True):
    """Get thresholded adjacency with a bivariate spline cutoff."""
    # Set the z-score:
    z = set_zscore(p=0.01, n=corr.shape[0]) if z is None else z
    logging.info("Thresholding by sparsity with estimated mean / variance" +
                 "bivariate splines, with z=" + str(np.round(z, 2)))

    # Get the smoothed cutoffs:
    corr = corr - np.diag(np.diag(corr))  # Remove diagonal for adjacency
    smean, ssd, xloc, dmarg = get_binned_cutoffs(
        corr, gmarg, gmarg_cut=0.75, min_sd=min_sd, zero_outliers=zero_outliers)

    # Make a df from the binned mean and standard deviation estimates:
    zcut = smean + ssd * z  # Correlation cutoff given the zscore
    min_cut = np.min(zcut)
    logging.info("Min cut: " + str(round(min_cut, 2)))
    smean = sparse.coo_matrix(smean * (ssd > 0))  # Mult to ensure same ind
    ssd = sparse.coo_matrix(ssd)
    statsdf = pd.DataFrame({"x": xloc[smean.row] + 1,
                            "y": xloc[smean.col] + 1,
                            "mean": smean.data,
                            "sd": ssd.data})

    # Turn correlation into i,j,v
    if not keep_all_z:
        # Remove very low correlations for efficiency / space
        cmat = corr * (corr >= (min_cut * .65))
    else:
        cmat = corr
    cmat = sparse.coo_matrix(cmat)
    crdf = pd.DataFrame({"x": dmarg[cmat.row], "y": dmarg[cmat.col],
                         "r": cmat.row, "c": cmat.col, "dat": cmat.data})
    crdf = crdf.merge(statsdf)
    crdf['z'] = (crdf['dat'] - crdf['mean']) / crdf['sd']

    # Create zscored matrix and thresholded adjacency matrix:
    zmat = sparse.coo_matrix((crdf.z, (crdf.r, crdf.c)), shape=cmat.shape)
    crdf = crdf[crdf.z >= z]
    amat = sparse.coo_matrix((crdf.dat, (crdf.r, crdf.c)), shape=cmat.shape)
    del(crdf, cmat)
    return(zmat, amat, zcut)


def adj_from_sd_est(corr, corr_sd, gmarg, z=None, min_sd=0.01):
    """Get thresholded adjacency by estimated standard deviation."""
    # Set a z-threshold:
    z = set_zscore(p=0.01, n=corr.shape[0]) if z is None else z
    logging.info("Thresholding from given sd estimate with z=" +
                 str(np.round(z, 2)))
    # Threshold adjacency matrix
    corr = corr - np.diag(np.diag(corr))
    corr_sd[corr_sd < min_sd] = min_sd
    corr_z = corr / corr_sd
    zcorr = corr * (corr_z >= z)
    aw = sparse.coo_matrix(zcorr)
    return aw
