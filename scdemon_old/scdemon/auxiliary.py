#!/usr/bin/python
"""Auxiliary functions for scdemon."""
# --------------------------------------------
# Auxiliary functions
# Updated: 06/28/21
# --------------------------------------------
import pandas as pd
import numpy as np
import scanpy as sc
import logging


def vcorrcoef(X, y):
    """Vectorized correlation, matrix vs. vector."""
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum(np.array(X - Xm) * (y - ym), axis=1)
    r_den = np.sqrt(np.sum(np.array(X - Xm) ** 2, axis=1)
                    * np.sum((y - ym) ** 2))
    r = r_num / r_den
    return r


def calculate_svd_covar_corr(ut, obsdf, cvlist, cv_mats={}):
    """Calculate covariate correlations with SVD."""
    for covar in cvlist:
        cvcol = obsdf[covar]
        if cvcol.dtype.name == "category":
            covar_dummy = pd.get_dummies(cvcol)
            cv_mats[covar] = np.zeros(
                (ut.shape[0], covar_dummy.shape[1]))
            for i, cv_lvl in enumerate(covar_dummy.columns):
                cv_mats[covar][:, i] = vcorrcoef(
                    ut, covar_dummy[cv_lvl].to_numpy().T)
        else:
            cvcol = cvcol.to_numpy()
            if np.sum(np.isnan(cvcol)) > 0:
                ind = np.where((1 - np.isnan(cvcol)) == 1)[0]
                cv_mats[covar] = vcorrcoef(
                    ut[:, ind], cvcol[ind, np.newaxis].T
                )[:, np.newaxis]
            else:
                cv_mats[covar] = vcorrcoef(
                    ut, cvcol[:, np.newaxis].T)[:, np.newaxis]
    return(cv_mats)


# Simple recipes for example code:

def recipe_preprocess(adata):
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    logging.info("Preprocessed example dataset")


def recipe_annotate(adata):
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    # sc.pl.umap(adata, color="leiden")
    logging.info("Annotated example dataset")


def recipe_full(adata, preprocess, annotate):
    if preprocess:
        recipe_preprocess(adata)
    if annotate:
        recipe_annotate(adata)
