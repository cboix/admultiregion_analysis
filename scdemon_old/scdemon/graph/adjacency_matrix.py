#!/usr/bin/python
"""Class for handling an adjacency matrix."""
# --------------------------------------
# Class for handling an adjacency matrix
# Updated: 06/28/21
# --------------------------------------
from .utils_adjacency import (
    prune_degree,
    prune_scale,
    prune_knn,
    zscore_from_bivariate_cutoff,
    adj_from_sd_est
)

import logging
import numpy as np
from scipy import sparse


class adjacency_matrix(object):
    """
    Turn a correlation matrix into an adjacency matrix.

    Process given correlation matrix into the adjacency matrix.
    """

    def __init__(
        self,
        corr,
        corr_sd=None,
        cutoff=None,
        use_zscore=True,
        labels=None,
        z=None,
        margin=None,
        knn_k=None,  # k-NN
        scale=None,
        zero_outliers=True,
        degree_cutoff=1
    ):
        """Initialize adjacency matrix class."""
        self.corr = corr
        self.corr_sd = corr_sd  # Estimated sd for each pair
        self.labels = labels
        self.indices = np.arange(len(self.labels))
        self.cutoff = cutoff
        self.z = z  # For adjacency cutoff
        self.zero_outliers = zero_outliers  # For getting estimates of cutoffs
        self.use_zscore = use_zscore
        if not self.use_zscore:  # Use cutoff if not using zscore method
            if self.cutoff is None:
                self.cutoff = 0.4

        self.margin = margin  # Margin (fraction non-zero) of dataset
        self.knn_k = knn_k
        self.scale = scale
        self.degree_cutoff = degree_cutoff
        self.adjacency = None
        logging.debug("Margin: " + str(self.margin.shape))
        logging.debug("Corr: " + str(self.corr.shape))

    def get_adjacency(self, keep_all_z=False):
        """Get the adjacency matrix and the final kept labels."""
        if self.adjacency is None:
            self._process_adjacency(keep_all_z=keep_all_z)
        return(self.adjacency, self.labels, self.indices)

    def _process_adjacency(self, keep_all_z=False):
        """Process given correlation matrix into the adjacency matrix."""
        self._create_adjacency_matrix(keep_all_z=keep_all_z)

        # Prune the created matrix:
        if self.scale is not None:
            self._prune_by_scale()
        if self.knn_k is not None:
            self._prune_by_knn()

        # Clean up nodes with low or 0 degree after pruning by other methods:
        self._prune_by_degree()

    def _create_adjacency_matrix(self, keep_all_z=False):
        """Create adjacency matrix from given correlation."""
        logging.info("Creating adjacency matrix from given correlation.")
        if not self.use_zscore:
            self.adjacency = self.corr * (self.corr > self.cutoff)
            self.adjacency = self.adjacency - np.diag(np.diag(self.adjacency))
            self.adjacency = sparse.coo_matrix(self.adjacency)
        else:
            if self.corr_sd is None:
                self.zmat, self.adjacency, self.zcut = \
                    zscore_from_bivariate_cutoff(
                        self.corr, self.margin, z=self.z,
                        zero_outliers=self.zero_outliers,
                        keep_all_z=keep_all_z)
            else:
                self.adjacency = adj_from_sd_est(
                    self.corr, self.corr_sd, self.margin, z=self.z)
        # Save the full and mostly empty adjacency for multiplexing graphs:
        # TODO: Set option to save or not?
        self.full_adjacency = self.adjacency
        self.full_labels = self.labels

    def _prune_by_scale(self):
        """Prune adjacency by a relative cutoff of outgoing edges."""
        logging.info("Pruning by relative cutoff, scale=" + str(self.scale))
        self.adjacency = prune_scale(self.adjacency, scale=self.scale)

    def _prune_by_knn(self):
        """Prune adjacency by a k-NN keeping top k edges."""
        logging.info("Pruning by KNN, K=" + str(self.knn_k))
        self.adjacency = prune_knn(self.adjacency, k=self.knn_k)

    def _prune_by_degree(self):
        """Prune adjacency by the degree of each node."""
        logging.info("Pruning by degree, cutoff=" + str(self.degree_cutoff))
        # Remove nodes with no/very few links:
        self.adjacency, ind = prune_degree(self.adjacency, self.degree_cutoff)
        self.labels = self.labels[ind]
        self.indices = self.indices[ind]
        self.adjacency = sparse.csr_matrix(self.adjacency)
