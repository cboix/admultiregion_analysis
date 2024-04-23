#!usr/bin/python
"""Experiments for co-expression v2."""
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# Using data loader + modules + graph handlers
# Created: 12/01/21
# -------------------------------------------------------------
import logging
import numpy as np
import scdemon as sm
from adata_handling import adata_loader
import argparse

from matplotlib import pyplot as plt
import seaborn as sns

import socket
domain = socket.getfqdn()
if 'broadinstitute.org' in domain or 'private' in domain:
    topdir = '/home/cboix/data/'
else:
    topdir = '/home/cboix/'

dbdir = topdir + 'DEVTRAJ/db/'
datadir = dbdir + 'multiRegion/'
resultsdir = dbdir + 'multiRegion/modules/'
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
parser.add_argument('--res', metavar='res', type=float,
                    help='leiden clustering resolution', default=2)
args = parser.parse_args()

# Set logging level (TODO: set by argparser)
logging.basicConfig(level=logging.INFO)
logging.info(str(args))


args.ct = 'Ast'
args.st = None

# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
scdata = adata_loader(celltype=args.ct, subtype=args.st,
                      dbdir=dbdir,
                      normalize=True,
                      usecombat=False, usehvg=False,
                      remove_mt=False, filter_TH=False)
scdata.load()


# Set up the modules computation object:
# --------------------------------------
max_k = 100
top_k = 100
csuff = scdata.csuff + "_expt_z" + str(args.z)
mod = sm.scdemon(scdata.adata, h5ad_file=scdata.h5ad_file,
                  csuff=csuff, imgdir=imgdir, seed=1,
                  svd_k=max_k, filter_expr=args.filt,
                  z=args.z, calc_raw=False)
mod.setup()
kept_genes = mod.genes  # For bootstraps



# Make graphs from base and from merge:
# -------------------------------------
graph_id = 'base'
k_ind = np.arange(1, top_k)
mod.make_svd_subset_graph(graph_id, k_ind=k_ind, resolution=args.res)

power_list = list(np.arange(0, 1, .25)) + [1]
graph_id = 'merge'
mod.make_merged_graph(graph_id, power_list=power_list,
                      resolution=args.res, keep_all_z=False)
# del(mod.graphs)
# mod.graphs = {}

mod.plot_graph(graph_id, attr="leiden",
               show_labels=True, frac_labels=.5, width=16)

graph_id = 'raw'
mod.make_graph(graph_id, resolution=args.res)


# Calculate the sum of squared errors for clusters:
# -------------------------------------------------
# Not very good scoring for genes:
def score_modules_sse(adata, mlist):
    totsse = 0
    for i in mlist.keys():
        x = mlist[i]
        subadata = adata[:,x]
        # Normalize first to z-scores:
        gsd = subadata.X.std(axis=0)
        gmean = subadata.X.mean(axis=0)
        modmean = ((subadata.X - gmean[None,:]) / gsd[None, :]).mean(axis=1)
        sse = (((subadata.X - gmean[None,:]) / gsd[None,:] -
                modmean[:,None])**2).sum(axis=0)
        sse = np.sum(sse) / subadata.shape[0]
        print(i, sse)
        totsse += sse
    return(totsse)

graph_id = 'base'
mlist = mod.graphs[graph_id].modules['leiden']
kg = mod.graphs[graph_id].graph.vs['name']
out = score_modules_sse(mod.adata, mlist)

graph_id = 'merge'
mlist = mod.graphs[graph_id].modules['leiden']
kg2 = mod.graphs[graph_id].graph.vs['name']
out2 = score_modules_sse(mod.adata, mlist)

graph_id = 'raw'
mlist = mod.graphs[graph_id].modules['leiden']
kg3 = mod.graphs[graph_id].graph.vs['name']
out3 = score_modules_sse(mod.adata, mlist)


# Pretty much equivalent, per gene:
print(out / len(kg))
print(out2 / len(kg2))
print(out3 / len(kg3))

mlist = mod.get_modules(graph_id, print_modules=True)
mod.find_gene(graph_id, 'SLC16A3')
mod.find_gene(graph_id, 'NRXN3')
mod.find_gene(graph_id, 'GAPDH')



# Score the GO enrichments for each gene in each graph:
# -----------------------------------------------------
gpres = mod.get_goterms(graph_id)



