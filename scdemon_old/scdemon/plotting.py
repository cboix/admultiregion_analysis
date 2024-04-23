#!/usr/bin/python
"""Plotting functions for modules."""
import logging
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

# For graphs
import textwrap


def plot_cell_umap(umat, c, plotname=None, ax=None, title=None,
                   s=1, width=8, axlab=False, cmap='viridis'):
    """Plot single umap with given scores as colors."""
    mx = np.max(umat, axis=0)
    mn = np.min(umat, axis=0)
    mid = (mx + mn) / 2
    if ax is None:
        # Define plotting range:
        rn = mx - mn
        height = width * rn[1] / rn[0]
        plt.figure(figsize=(width, height))
        ax = plt.gca()
    if title is not None:
        tw = textwrap.fill(title, 18)
        # ax.set_title(tw, fontdict={"fontsize": 7})
        ax.text(mid[0], mid[1], tw, fontdict={"fontsize": 7},
                ha='center', va='center')
    # TODO: bring red points to front.
    ax.set_facecolor("white")
    ax.scatter(umat[:, 0], umat[:, 1], s=s, c=c, marker=".",
               edgecolors="none", cmap=plt.get_cmap(cmap))
    # Remove tick labels:
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    if axlab:
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
    if plotname is not None:
        plt.tight_layout()
        plt.savefig(plotname, dpi=350, bbox_inches="tight")
        plt.close()


# TODO: Should width be the full size?
def plot_umap_grid(umat, scores, titles, plotname=None,
                   ind=None, sel=None, width=2, s=0.5):
    """Plot modules scores on UMAP as a grid."""
    nplot = scores.shape[1] if sel is None else len(sel)
    # Determine grid shape: NR/NC
    nr = int(np.round(np.sqrt(nplot) * 0.8))
    nc = int(np.ceil(nplot / (nr * 1.0)))
    # Select specific cells:
    if ind is not None:  # TODO Fix this delimitation with quantiles?
        umat = umat[ind, :]
        scores = scores[ind, :]
    # Define plotting range so each plot preserves aspect:
    mx = np.max(umat, axis=0)
    mn = np.min(umat, axis=0)
    rn = mx - mn
    height = width * rn[1] / rn[0]  # Size of each plot
    # Make subplots:
    fig, axs = plt.subplots(nrows=nr, ncols=nc,
                            gridspec_kw={"hspace": 0.01, "wspace": 0.01},
                            figsize=(width * nc, height * nr))
    for i in range(nr):
        for j in range(nc):
            k = i * nc + j
            ax = axs[i, j]
            if k < nplot:
                k = k if sel is None else sel[k]
                title = None if titles is None else titles[k]
                plot_cell_umap(umat=umat, c=scores[:, k],
                               title=title, s=s, ax=ax)
            else:
                ax.set_facecolor("white")
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
    # Save figure:
    plt.tight_layout()
    plt.savefig(plotname, dpi=350, bbox_inches="tight")
    plt.close()
    logging.info("Plotted grid of modules on UMAP to: " + plotname)


def plot_svd_corr(cvlist, cv_hratio, cv_mats, cv_ticks, svec, plotname,
                  cbar=False):
    """Plot heatmaps of correlation of svd components and covariates."""
    hfull = np.sum(cv_hratio)
    w = 2.5 * len(svec) / 20 + 1.5 + 0.8 * cbar
    h = 1 * hfull / 6 + 0.1 * len(cvlist)
    fig, axs = plt.subplots(len(cv_mats), 1, figsize=(w, h),
                            gridspec_kw={
                                "height_ratios": cv_hratio,
                                "left": 1.5 / w,
                                "right": 0.99 - 0.8 / w * cbar,
                                "top": 1 - 0.1 / h,
                                "bottom": 0.4 / h})
    # Add heatmaps:
    sfact = np.max(svec) / svec
    for i, covar in enumerate(cvlist):
        cmat = cv_mats[covar].T
        cmat = cmat * sfact[np.newaxis, :]
        sns.heatmap(cmat,
                    yticklabels=cv_ticks[covar],
                    xticklabels=(i == len(cvlist) - 1),
                    cmap="RdBu", ax=axs[i], cbar=cbar, center=0)
        axs[i].set_ylabel(covar, fontsize=12, rotation=0,
                          ha="right", va="center")
    # Save figure:
    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    logging.info("Plotted graph to " + plotname)
