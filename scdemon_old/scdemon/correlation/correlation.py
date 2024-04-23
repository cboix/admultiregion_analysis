#!/usr/bin/python
"""Class for handling correlations in scdemon."""
# --------------------------------------------
# Class for handling correlations in scdemon:
# Updated: 07/02/21
# --------------------------------------------
from .utils_correlation import (
    calculate_correlation,
    calculate_correlation_estimate,
    calculate_correlation_estimate_sd,
    search_FDR_cutoff
)

import logging
import numpy as np
import pandas as pd
import fbpca
from scipy import sparse

# For adjustments
from scipy.interpolate import UnivariateSpline
from sklearn.linear_model import LinearRegression


# Standalone correlation handler:
class correlation(object):
    """
    Handles correlation computation on raw + decorrelated data.

    Arguments:
        X: Matrix (n x d) for which to calculate the d x d correlation
        U, s, V: SVD of X matrix (not required).
            U is (n x k), s is (k), and V is in (k x d)
        k: Number of SVD components to compute (default = 100)
        calc_raw: Whether to also return correlation on untransformed data

    Functions:
        get_correlation()
            Obtain the overall correlation matrices

        get_correlation_subset(indices)
            Obtain the correlation matrix for a subset of SVD components
    """

    def __init__(self, X, margin, U=None, s=None, V=None,
                 k=100, calc_raw=False, center=False, mean=None):
        """Initialize correlation handling object."""
        # Input:
        self.X = X
        self.margin = margin
        self.U = U
        self.s = s
        self.V = V
        # Arguments:
        self.N = self.X.shape[0]
        self.P = self.X.shape[1]
        self.k = k
        self.calc_raw = calc_raw
        self.center = center
        self.mean = mean
        # Output:
        self.corr_est = None
        self.corr_raw = None
        self.use_v = True
        logging.debug("X: " + str(self.X.shape))
        logging.debug("Margin: " + str(self.margin.shape))

    def setup(self):
        self._calculate_pca()

    def get_correlation(self, keep_first=False, power=0):
        """Get correlation estimates from SVD and raw if calc_raw is True."""
        if self.corr_est is None or self.power != power:
            self._calculate_estimated_correlation(
                keep_first=keep_first, power=power)
        self.corr_raw = self.get_raw_correlation()
        return (self.corr_est, self.corr_raw)

    def get_raw_correlation(self, force=False):
        if force:
            self.calc_raw = True
        if self.corr_raw is None and self.calc_raw:
            self.corr_raw = calculate_correlation(self.X, center=self.center)
        # NOTE: else corr_raw is already None
        return(self.corr_raw)

    def get_correlation_subset(self, indices, power=0):
        """Get a correlation estimate from a subset of SVD columns."""
        logging.info("Calculating correlation from subset of SVD columns")
        corr_subset = calculate_correlation_estimate(
                self.U, self.s, self.V, power=power, indices=indices)
        return corr_subset

    def _calculate_estimated_correlation(self, power=0, keep_first=False):
        self.power = power
        self._calculate_pca()
        if self.use_v:
            # Keep or remove first component (defaults to removing):
            indices = None if keep_first else np.arange(1, self.U.shape[1])
            logging.debug(indices)
            self.corr_est = calculate_correlation_estimate(
                self.U, self.s, self.V, power=power, indices=indices)
        else:
            # TODO: Not correct (transposed) - fix:
            # TODO: Also pass power and keep_first arguments to this:
            self._calculate_data_transformation()
            self.corr_est = calculate_correlation(
                self.Xw, center=self.center)

    def _check_pca(self):
        # Check conditions for re-calculating PCA:
        if self.U is None or self.s is None or \
                self.V is None or self.U.shape[1] < self.k:
            return False
        else:
            return True

    def _calculate_pca(self):
        """Calculate PCA if not computed already."""
        if not self._check_pca():
            logging.info("Calculating PCA of X with fbpca")
            self.U, self.s, self.V = fbpca.pca(
                self.X, k=self.k, raw=not self.center)

    def _calculate_data_transformation(self):
        """Calculate decorrelation transformation of X with SVD."""
        logging.info("Calculating transformation of X with SVD")
        sinv = np.diag(1 / self.s)
        XtU = self.X.T.dot(self.U)
        self.Xw = XtU.dot(sinv.dot(self.U.T)) / self.V.shape[1]

    # TODO: allow both types of adjustments (centralize in graph vs. not
    # in graph types)

    # Perform an adjustment to the correlation estimate
    # empirical t_ij statistics
    def get_adj_correlation(self):
        # Backward compatible for early versions
        self._calculate_var_cov()
        self._adjust_test_statistic()
        self._define_FDR_cutoff()
        return (self.t_ij_adj, self.t_hat)

    # Calculate the variance of the approximate covariance
    # e.g. Var(V_i * V_j^T) for use in calculating T_ij
    # TODO: Working on solution for test statistic, incomplete
    def _calculate_var_cov(self):
        # TODO: Improve, must calc ZCA right now:
        X_z = self.U.dot(self.V)
        X_z2 = X_z ** 2
        s_ij = X_z.T.dot(X_z) / self.N  # TODO: Use VVT estimate
        E1 = X_z2.T.dot(X_z2) / self.N
        E2 = (s_ij) ** 2
        self.theta = E1 - E2
        del (X_z, X_z2, E1, E2)

        # NOTE: Without S, must add this scaling factor (from C^-1)
        self.t_ij = s_ij / np.sqrt(self.theta) * np.sqrt(self.P)
        self.t_ij = self.t_ij - np.diag(np.diag(self.t_ij))

    # TODO: Clean this up:
    def _adjust_test_statistic(self):
        # The T_ij ~ N(0,1) only holds for the non-sparse parts of the data.
        # In the low-sparsity part of the data, we see high variance normal,
        # So we must adjust this statistic.
        # For now we adjust empirically, while we are working on a proper
        # correction.
        # Adjusting statistic:
        zc = sparse.coo_matrix(np.triu(self.t_ij))
        sdf = pd.DataFrame(
            {
                "p": self.margin[zc.row],
                "q": self.margin[zc.col],
                "r": zc.row,
                "c": zc.col,
                "dat": zc.data,
            }
        )
        sdf["pq"] = np.sqrt(sdf.p * sdf.q)

        # Bin and adjust:
        xs = np.linspace(np.min(sdf.pq), np.max(sdf.pq), 50)
        sdf["d_pq"] = np.digitize(sdf.pq, xs)
        u, c = np.unique(sdf["d_pq"], return_counts=True)
        NS = len(u)
        # Calculate the means/sdevs:
        smean = np.zeros(NS)
        ssd = np.zeros(NS)
        for i, ui in enumerate(u):
            # ci = c[i]
            rind = sdf["d_pq"] == ui
            vec = sdf.dat[rind]
            smean[i] = np.mean(vec)
            ssd[i] = np.std(vec)

        # Past halfway, set mean to flat behavior (othw skewed by diff. alt)
        # TODO: Find better way to do this
        ind = np.where(xs > 0.5)[0][0]
        smean[-ind:] = smean[ind]

        # Linear regression for the variance (due to null model):
        regr = LinearRegression()
        mids = (np.array([0] + list(xs[:-1])) + xs) / 2
        regr.fit(mids[:, None], ssd, np.sqrt(c))
        sdf["pred_var"] = regr.predict(sdf.pq[:, None]) ** 2

        # Univariate spline for the mean (inflated values):
        spl = UnivariateSpline(mids[:, None], smean, w=np.sqrt(c), k=2)
        sdf["pred_mean"] = spl(sdf.pq[:, None])

        # Adjust - pretty close to what we want.
        sdf["dat2"] = (sdf["dat"] - sdf["pred_mean"]) / sdf["pred_var"]
        # Turn back into matrix-like
        cu = sparse.coo_matrix(
            (sdf.dat2, (sdf.r, sdf.c)),
            shape=self.t_ij.shape)
        cu = cu.toarray()
        cu = cu + cu.T
        self.t_ij_adj = cu  # Temporary fix

    # TODO move this to graph.py ??
    def _define_FDR_cutoff(self, alpha=0.01):
        self.t_hat = search_FDR_cutoff(self.t_ij_adj, self.P, alpha=alpha)

    # TODO: Estimate the raw correlation SD as well.
    def estimate_corr_sd(self, nperm=25, seed=1):
        corr_mean, corr_sd = calculate_correlation_estimate_sd(
            self.U, self.s, self.V, nperm=nperm, seed=seed)
        return (corr_mean, corr_sd)

    # TODO: Remove later; for comparison (X.T) with non-centered
    def calculate_correlation_centered(self):
        logging.info("Calculating covariance of ZCA transformed X")
        logging.debug("Calculating PCA of X.T with fbpca")
        self.Ucn, self.scn, self.Vcn = fbpca.pca(self.X.T, k=self.k, raw=False)
        s = np.diag(1 / self.scn)
        X = self.X.toarray()
        self.X_zca = self.Vcn.dot(X)
        self.X_zca = self.Vcn.T.dot(s).dot(self.X_zca)
        self.corr_est = np.corrcoef(self.X_zca.T)
        del(X)
