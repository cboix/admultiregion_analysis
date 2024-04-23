#!usr/bin/python
"""Rank genes on the precomputed integrated trajectory."""
# --------------------------------------------------------------
# Compute an integrated UMAP across all neuron subtypes:
# Updated: 03/01/22
# --------------------------------------------------------------
import re
import logging
from time import time
from tqdm import tqdm

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing.pool import Pool

from pygam import LinearGAM, GAM, s, f

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)


# Load in the preprocessed data:
# ------------------------------
dbdir = '/home/cboix/data/DEVTRAJ/db/multiRegion/'
h5ad_file = dbdir + 'preprocessed_fordetraj_Exc_log1p.h5ad'
useset = 'testedgenes'
prefstr = '_preprocessed_fordetraj_Exc_log1p_' + useset + "_combat"


adata = sc.read(h5ad_file)

# Also get group2:
llist = adata.obs.leiden.tolist()
adata.obs['group2'] = ['1to2' if x in ['1', '2', '5']
                       else '' for x in llist]

sc.pl.rank_genes_groups_matrixplot(
    adata, n_genes=3, save=prefstr + '_ctalone_leiden_small.pdf',
    colorbar_title='mean z-score', layer='scaled',
    vmin=-2, vmax=2, cmap='RdBu_r')



# Functions to score genes by their GAM fits in either trajectory:
# ----------------------------------------------------------------
# TODO: multiprocessing for this... (have ~ 60-80 nodes)
# 7000 * 2 * 7 / 3600 ~ 27 hrs --> req. parallel or multiprocessing
# TODO: Better trajectory needed
X = adata.obs['C2'].to_numpy()
X = 1 - X / np.max(X)


def predict_gene(adata, gene, group, groupvar='group', pst='C2', ci=.95):
    """Predict GAM fit for a single gene and return CI for PST term."""
    ind = adata.obs[groupvar] == group
    # TODO: speed up getting data.
    X = adata.obs[pst].to_numpy()
    X = 1 - X / np.max(X)
    subX = X[ind]
    y = adata[ind, gene].X.toarray().T[0]
    # print("GOT DATA")
    gam = LinearGAM(s(0))
    # print("BUILT GAM")
    # gam.gridsearch(subX[:, np.newaxis], y, progress=True)
    # print("GRIDSEARCH DONE")
    gam.lam[0][0] = 250
    gam.fit(subX[:, np.newaxis], y)
    # print('LAMBDA', gam.lam[0][0])
    XX = gam.generate_X_grid(term=0)
    pred = gam.predict(XX)
    intervals = gam.partial_dependence(term=0, X=XX, width=.95)[1]
    return(pred, intervals)


def score_se(intervals):
    """Score gene based on CI intervals for the pseudotime term."""
    int_se = (intervals[:, 1] - intervals[:, 0]) / 4
    int_z = np.abs(intervals / int_se[:, np.newaxis])
    maxscore = np.max(int_z.min(1))
    avgscore = np.mean(int_z.min(1))
    return(avgscore, maxscore)


# Run across all genes:
# ---------------------
# dgenes = ['FAIM2', 'PRNP']
dgenes = adata.var_names.tolist()
NG = len(dgenes)
scores = []
pmat = np.zeros((NG, 100, 2))
imat = np.zeros((NG, 100, 2, 2))
groups = ['top', 'bottom']
runs = [[i, j] for i in range(NG) for j in range(2)]


def worker(k, mode='vert', runs=runs):
    """Run a single worker based on runs above, for multiproc."""
    args = runs[k]
    i = args[0]
    j = args[1]
    gene = dgenes[i]
    group = groups[j]
    print(gene, group)
    if mode == 'vert':
        pred, intervals = predict_gene(adata, gene, group,
                                       groupvar='group', pst='C2')
    elif mode == 'horiz':
        pred, intervals = predict_gene(adata, gene, group,
                                       groupvar='group2', pst='C1')
    avgscore, maxscore = score_se(intervals)
    procscores = [gene, group, maxscore, avgscore]
    # print(procscores)
    pmat[i, :, j] = pred
    imat[i, :, j, 0] = intervals[:, 0]
    imat[i, :, j, 1] = intervals[:, 1]
    return(procscores)


ts = time()
scores = [worker(i) for i in tqdm(range(len(runs)))]
prinpstdir = dbdir + 'pseudotime/'
df.to_csv(pstdir + 'gene_scores_integrated_detraj.tsv', sep="\t")

t(time() - ts)

df = pd.DataFrame(scores)
df.columns = ['gene', 'group', 'max', 'avg']
pstdir = dbdir + 'pseudotime/'
df.to_csv(pstdir + 'gene_scores_integrated_detraj.tsv', sep="\t")


# Save IMAT and PMAT (for plotting fits)
# --------------------------------------
pwide = pmat.swapaxes(1, 2).reshape(pmat.shape[0], -1)
np.savetxt(pstdir + 'gene_order.txt', dgenes, fmt="%s")
np.savetxt(pstdir + 'gene_prediction_matrix.tsv', pwide, delimiter="\t")

imat_top = imat[:, :, :, 0]
imat_bot = imat[:, :, :, 1]
itwide = imat_top.swapaxes(1, 2).reshape(imat_top.shape[0], -1)
ibwide = imat_bot.swapaxes(1, 2).reshape(imat_bot.shape[0], -1)
np.savetxt(pstdir + 'gene_ci95top_matrix.tsv', itwide, delimiter="\t")
np.savetxt(pstdir + 'gene_ci95bottom_matrix.tsv', ibwide, delimiter="\t")


# Run across all genes horizontal trajectory:
# TODO: Collapse with above:
# ---------------------
# dgenes = ['FAIM2', 'PRNP']
dgenes = adata.var_names.tolist()
NG = len(dgenes)
scores = []
pmat = np.zeros((NG, 100))
imat = np.zeros((NG, 100, 2))
runs = [i for i in range(NG)]


def worker(k, mode='vert', runs=runs):
    """Run a single worker based on runs above, for multiproc."""
    i = runs[k]
    gene = dgenes[i]
    group = '1to2'
    # print(gene, group)
    if mode == 'vert':
        pred, intervals = predict_gene(adata, gene, group,
                                       groupvar='group', pst='C2')
    elif mode == 'horiz':
        pred, intervals = predict_gene(adata, gene, group,
                                       groupvar='group2', pst='C1')
    avgscore, maxscore = score_se(intervals)
    procscores = [gene, group, maxscore, avgscore]
    # print(procscores)
    pmat[i, :] = pred
    imat[i, :, 0] = intervals[:, 0]
    imat[i, :, 1] = intervals[:, 1]
    return(procscores)


ts = time()
scores = [worker(i, mode='horiz') for i in tqdm(range(len(runs)))]
print(time() - ts)

pstdir = dbdir + 'pseudotime/'
df = pd.DataFrame(scores)
df.columns = ['gene', 'group', 'max', 'avg']
df.to_csv(pstdir + 'gene_scores_integrated_detraj_horiz.tsv', sep="\t")


# Save IMAT and PMAT (for plotting fits)
# --------------------------------------
np.savetxt(pstdir + 'gene_order_horiz.txt', dgenes, fmt="%s")
np.savetxt(pstdir + 'gene_prediction_matrix_horiz.tsv', pmat, delimiter="\t")

imat_top = imat[:, :, 0]
imat_bot = imat[:, :, 1]
np.savetxt(pstdir + 'gene_ci95top_matrix_horiz.tsv', imat_top, delimiter="\t")
np.savetxt(pstdir + 'gene_ci95bottom_matrix_horiz.tsv', imat_bot, delimiter="\t")


# Plot some of these fits:
# ------------------------
pltpref = 'figures/umap' + prefstr + '_ctalone'

def savelocal(fig, suffix, prefix=pltpref, dpi=350):
    """Save plot with current local prefix."""
    plotname = prefix + "_" + suffix + ".png"
    plt.tight_layout()
    plt.savefig(plotname, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(plotname)


X = adata.obs['C2'].to_numpy()
X = 1 - X / np.max(X)

dgenes = ['KCNIP4','NRXN1', 'MT-ND3', 'MT-CO3', 'CALM1', 'CALM3']
suff = 'initial'
# dgenes = ['TARBP1','SLC27A6', 'LPP', 'PKM', 'UBB', 'PCSK1N']
# suff = 'branching'
dgenes = ['PDE5A', 'BEX3', 'ATP1B1', 'NRG3', 'RALYL', 'PDE4D']
suff = 'horizontal'

# NOTE: We are plotting the combat-corrected datapoints.
NG = len(dgenes)
plt.figure()
fig, axs = plt.subplots(2, NG, figsize=(4 * NG * .5, 4))
for i, axvec in enumerate(axs):
    group = ['top', 'bottom'][i]
    gname = ['Left', 'Right'][i]
    ind = adata.obs['group'] == group
    subX = X[ind]
    for j, ax in enumerate(axvec):
        gene = dgenes[j]
        y = adata[ind, gene].X.toarray().T[0]
        gam = LinearGAM(s(0))
        # gam.gridsearch(subX[:, np.newaxis], y, progress=True)
        gam.lam[0][0] = 250
        gam.fit(subX[:, np.newaxis], y)
        # plotting
        XX = gam.generate_X_grid(term=0)
        ax.scatter(subX, y, alpha=0.15, edgecolors='none', marker='.', s=2, c='grey')
        ax.plot(XX[:, 0], gam.predict(XX))
        ax.plot(XX[:, 0], gam.prediction_intervals(XX, width=.95),
                c='r', ls='--')
        if i == 0:
            ax.set_title(gene)
        if j == 0:
            ax.set_ylabel(gname)

savelocal(fig, 'selgenes_linearGAM_' + suff, dpi=450)


# Some multiprocessing options that didn't work well.
# ---------------------------------------------------
# ts = time()
# pool = Pool(10)
# scores = pool.map(worker, range(len(runs)))
# for k in range(len(runs)):
#     worker(runs[k])

# from joblib import Parallel, delayed  # Slow.
# ts = time()
# scores = Parallel(n_jobs=12)(delayed(worker)(k) for k in range(len(runs)))
# print(time() - ts)

# from os import getpid
# def worker(procnum):
#     print('I am number %d in process %d' % (procnum, getpid()))
#     return getpid()

# pool = Pool(processes = 3)
# print(pool.map(worker, range(5)))


# Aggregate scores
# Aggregate predictions

# Save scores and predictions:


# Plot predictions as heatmap, use to rank genes
# Rank genes along temporal time-course
# RTN4, CALM1, AKT1, etc. along bottom pathway
