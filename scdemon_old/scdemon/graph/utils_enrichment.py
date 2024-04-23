#!/usr/bin/python
"""Functions for module enrichments on dataFrames."""

import logging
import numpy as np
import pandas as pd

# Plotting:
from matplotlib import pyplot as plt
import seaborn as sns

# For enrichment:
from scipy.stats import hypergeom


def calc_df_enrichment(df, col, genes, module_match):
    """
    Calculate enrichment of modules in dataframe.

    Given a dataframe with genes (df.gene) and some split (df[col]),
    calculate the enrichment of each split in the modules
    """
    # Count genes in each:
    keys = np.flip(np.sort(pd.unique(df[col])))
    cmat = np.zeros((np.max(module_match) + 1, len(keys)), int)
    for i, key in enumerate(keys):
        df_genes = df.gene[df[col] == key]
        ind = np.in1d(genes, df_genes)
        u, c = np.unique(module_match[ind], return_counts=True)
        cmat[u, i] = c

    # Hypergeom for each box:
    rmarg = np.sum(cmat, 1)
    cmarg = np.sum(cmat, 0)
    pmat = np.zeros(cmat.shape)
    rmat = np.zeros(cmat.shape)
    M = np.sum(cmat)  # Total # genes
    stats = []
    for i in range(cmat.shape[0]):
        for j in range(cmat.shape[1]):
            x = cmat[i, j]  # Drawn matching module
            n = rmarg[i]  # Any matching module
            N = cmarg[j]  # Total drawn
            pval = hypergeom.sf(x - 1, M, n, N)
            ratio = (x / N) / (n / M)
            pmat[i, j] = -np.log10(pval)
            rmat[i, j] = ratio
            # Track statistics to return as DataFrame:
            stats.append([i, j, x, n, N, M, ratio, pval])

    # Statistics table:
    statsdf = pd.DataFrame(
        np.array(stats),
        columns=[
            'module',
            'j',
            'x',
            'n',
            'N',
            'M',
            'ratio',
            'p'])
    for column in ['module', 'j', 'x', 'n', 'N', 'M']:
        statsdf[column] = statsdf[column].astype(int)

    statsdf['key'] = keys[statsdf.j.to_numpy()]

    return pmat, rmat, keys, statsdf


def format_enrichment_matrices(rmat, pmat):
    """Format the enrichment matrices for plotting."""
    # Set NaNs and 0 to min (for plotting)
    rmat[np.isnan(rmat)] = 1
    rmat[rmat == 0] = np.min(rmat[rmat != 0])
    # Format the annotation for plotting:
    annot = np.empty(rmat.shape, dtype="<U10")
    annot[:] = ""
    annot[10 ** -pmat < 0.05] = "*"
    annot[10 ** -pmat < 0.01] = "**"
    annot[10 ** -pmat < 0.001] = "***"
    return (rmat, annot)


# TODO: Check works if not using all modules / or genes are missing
def plot_df_enrichment(df, col, genes, module_match, mnames, plotname,
                       title=None):
    """Plot enrichment of modules in dataframe."""
    # Calculate the enrichment:
    if isinstance(df, dict):
        statsdf = None
        pmat, rmat, keys, annot, wr = {}, {}, {}, {}, []
        for dkey in df.keys():
            pmat[dkey], rmat[dkey], keys[dkey], sdf = calc_df_enrichment(
                df[dkey], col, genes, module_match)
            rmat[dkey], annot[dkey] = format_enrichment_matrices(
                rmat[dkey], pmat[dkey])
            wr.append(len(keys[dkey]))
            sdf['dkey'] = dkey
            statsdf = pd.concat([statsdf, sdf])
    else:
        pmat, rmat, keys, statsdf = calc_df_enrichment(
            df, col, genes, module_match)
        rmat, annot = format_enrichment_matrices(rmat, pmat)

    # Plot heatmap of effect size with pvals on top:
    akws = {"ha": "center", "va": "center"}
    if isinstance(df, dict):
        dkey = list(rmat.keys())[0]
        w = 3 + rmat[dkey].shape[1] * len(df)
        h = 1 + 0.2 * (rmat[dkey].shape[0] - 1)
        fig, axs = plt.subplots(1, len(df), figsize=(w, h),
                                gridspec_kw={"width_ratios": wr,
                                             "hspace": 0.025, "wspace": 0.025})
        for i, dkey in enumerate(list(df.keys())):
            yt = False if i > 0 else mnames
            sns.heatmap(
                np.log2(rmat[dkey]), annot=annot[dkey], cmap="RdBu_r",
                yticklabels=yt, xticklabels=keys[dkey],
                ax=axs[i], center=0, cbar=True, fmt="s", annot_kws=akws)
            axs[i].set_aspect(0.5)
            axs[i].set_title(dkey)
    else:
        plt.figure(figsize=(3 + rmat.shape[1],
                            1 + 0.2 * (rmat.shape[0] - 1)))
        ax = plt.gca()
        sns.heatmap(
            np.log2(rmat), annot=annot, cmap="RdBu_r",
            yticklabels=mnames, xticklabels=keys,
            center=0, cbar=True, fmt="s", annot_kws=akws)
        ax = plt.gca()
        ax.set_aspect(0.5)
        if title is not None:
            ax.set_title(title)

    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    logging.info(plotname)
    return(statsdf)
