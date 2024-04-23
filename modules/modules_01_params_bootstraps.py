#!usr/bin/python
"""Compute co-expression modules for multiRegion data using scdemon."""
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# Using data loader + modules + graph handlers
# Created: 06/16/21 // Last updated: 11/27/21
# python $SBINDIR/modules/modules_01_params_bootstraps.py --ct Opc
# -------------------------------------------------------------
import logging
import os
import gc
import pandas as pd
import numpy as np
import scdemon as sm
from adata_handling import adata_loader
import argparse

import socket
domain = socket.getfqdn()
if 'broadinstitute.org' in domain or 'private' in domain:
    topdir = '/home/cboix/data/'
else:
    topdir = '/home/cboix/'

dbdir = topdir + 'DEVTRAJ/db/'
datadir = dbdir + 'multiRegion/'
resultsdir = dbdir + 'multiRegion/modules/'
bootdir = resultsdir + 'bootstrap/'
dedir = datadir + 'dereg/'
imgdir = topdir + 'DEVTRAJ/img/multiRegion/modules/'


# Build argparser so that we can run this from CLI / on batched:
# --------------------------------------------------------------
parser = argparse.ArgumentParser(description='Calculate scdemon modules.')
parser.add_argument('--ct', metavar='ct', type=str,
                    help='celltype')
parser.add_argument('--st', metavar='st', type=str,
                    help='subtype', default=None)
parser.add_argument('--z', metavar='z', type=float,
                    help='zscore', default=4.5)
parser.add_argument('--filt', metavar='filt', type=float,
                    help='filter low expressed genes', default=0.05)
args = parser.parse_args()

# Set logging level:
# NOTE: could also set logging level by argparser
logging.basicConfig(level=logging.INFO)
logging.info(str(args))


# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
scdata = adata_loader(celltype=args.ct, subtype=args.st,
                      dbdir=dbdir, normalize=True,
                      usecombat=False, usehvg=False,
                      remove_mt=False, filter_TH=False)
scdata.load()


# Set up the modules computation object:
# --------------------------------------
max_k = 200
csuff = scdata.csuff + "_z" + str(args.z)
mod = sm.scdemon(scdata.adata, csuff=csuff, imgdir=imgdir, seed=1,
                  svd_k=max_k, filter_expr=args.filt,
                  z=args.z, calc_raw=False)
mod.setup()
kept_genes = mod.genes  # For bootstraps


# Perform a hyperparameter search for k and the power of s:
# NOTE: Save in data frame, to allow pre-computed:
# ---------------------------------------------------------
paramfile = resultsdir + 'paramsearch_' + mod.csuff + '.tsv'
if not os.path.exists(paramfile):
    klist = list(np.arange(30, max_k, 10)) + [max_k]
    ng, nm = mod.get_k_stats(k_list=klist, power=0)
    logging.info(f"ng array - {ng}")
    logging.info(f"nm array - {nm}")
    # Joint score on number of modules and number of genes to pick top k:
    score = np.array(ng) / 100 + np.array(nm)
    paramdf = pd.DataFrame({'k': klist, 'ng': ng, 'nm': nm, 'score': score})
    paramdf.to_csv(paramfile, sep="\t")
else:
    paramdf = pd.read_csv(paramfile, sep="\t")
    logging.info(paramdf.head())
    ng = paramdf.ng.to_numpy()
    nm = paramdf.nm.to_numpy()

ki = np.argmax(paramdf.score.to_numpy())
top_k = paramdf.k.to_numpy()[ki]
logging.info(f"Setting k={top_k}, ng={ng[ki]}, nm={nm[ki]}")


# Set up bootstrap evaluation - individuals and functions:
# --------------------------------------------------------
# NOTE: We will have fewer, but more robust modules here
# An alternative is to run on the full dataset
# --------------------------------------------------------
bs_pct = 0.95
plist = list(scdata.adata.obs['projid'])
rlist = list(scdata.adata.obs['region'])
scdata.adata.obs['pr'] = [str(plist[i]) + "_" + str(rlist[i])
                          for i in range(len(scdata.adata.obs))]
batches = list(pd.unique(scdata.adata.obs['pr']))
bs_size = int(np.floor(bs_pct * len(batches)))


def get_bootstrap_data(scdata, batches, bs_size, kept_genes):
    """Get a bootstrap subset of adata."""
    bs_subset = set(np.random.choice(batches, size=bs_size, replace=False))
    ind = np.array([x in bs_subset for x in scdata.adata.obs.pr])
    logging.info(f"nbatch={len(bs_subset)} and ncell={np.sum(ind)}")
    adata = scdata.adata[ind, kept_genes]
    return(adata)


def get_zscore(adata, svd_k=top_k, power=0):
    """Get the zscored correlation from adata."""
    mod_bs = sm.scdemon(adata, csuff='bootstrap', imgdir=imgdir, seed=1,
                         svd_k=svd_k, filter_expr=0, z=4.5, calc_raw=False)
    mod_bs.setup()
    # Just compute adjacency (no modules, layout):
    mod_bs.make_graph('base', adjacency_only=True, power=power)
    zmat = mod_bs.graphs['base'].adj.zmat
    del(mod_bs, adata)
    gc.collect()
    return(zmat)


def get_bootstrap_zmat(i, outdir, suffix, scdata, batches,
                       bs_size, kept_genes, svd_k, power=0, rerun=False):
    """Get the zscored correlation for one bootstrap."""
    np.random.seed(i)
    logging.info(f"run {i}")
    zmatfile = outdir + f'zmat_{i}_p{power}_' + suffix + '.npz'
    if os.path.exists(zmatfile) and not rerun:
        zmat = np.load(zmatfile, allow_pickle=True)['zmat']
    else:
        adata = get_bootstrap_data(scdata, batches, bs_size, kept_genes)
        zmat = get_zscore(adata, svd_k=top_k, power=power)
        np.savez_compressed(file=zmatfile, zmat=zmat)
    return(zmat)


# Compute and save bootstrapped zscore matrix estimates:
# ------------------------------------------------------
nperm = 10
power_list = list(np.arange(0, 1, .25)) + [1]
for power in power_list:
    ztot = np.zeros((len(kept_genes), len(kept_genes)))
    ztotfile = bootdir + f'ztot_n{nperm}_p{power}' + mod.csuff + '.npz'
    if not os.path.exists(ztotfile):
        for i in range(nperm):
            # Get each of the matrices for the permutations:
            zmat = get_bootstrap_zmat(i, bootdir, mod.csuff,
                                      scdata, batches, bs_size,
                                      kept_genes=kept_genes, svd_k=top_k,
                                      power=power, rerun=False)
            ztot += zmat.toarray()
        # Average and save
        ztot = ztot / nperm
        np.savez_compressed(file=ztotfile, ztot=ztot)
