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


args.ct = 'Exc'
args.st = 'ECneurons'

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


# Look at important genes - why are they not captured?
# ----------------------------------------------------
xmean = mod.adata.var['xmean'].to_numpy()
agenes = ['SLC1A2', 'SLC1A3',
          'GRM3', 'LUZP2',
          'MT-ND3', 'MT-ND1',  # Captured
          'GFAP', 'CD44', 'OSMR', 'MAOB',  # Captured 2 / 3
          'VEGFA', 'HILPDA',
          'SYT1', 'KCNIP4']  # Bad, captured
# Visualize V (GFAP
# scdata.adata.mean

mod.adata.var.loc[agenes, :]
aind = np.array([np.where(mod.genes == x)[0][0] for x in agenes])

Vt = mod.cobj.V.T
s = mod.cobj.s
xmean = xmean - np.mean(xmean)
xcol = xmean / np.sqrt(np.sum(xmean**2))


def plot_mean(ax, tck, scmean=1/5, power=0):
    """Plot correlation with mean adjustments. Doesn't work."""
    mat = Vt[aind, :]
    mat = mat + xcol[aind][:, None] * scmean
    s_red = np.diag(s**power)
    X_cov = mat.dot(s_red).dot(mat.T)
    X_sd = np.sqrt(np.diag(X_cov))
    cv = X_cov / X_sd[:, np.newaxis]
    corr_est = cv / X_sd[np.newaxis, :]
    akws = {"ha": "center", "va": "center"}
    sns.heatmap(corr_est, ax=ax,
                xticklabels=tck, yticklabels=tck,
                annot=np.round(10 * corr_est, 1), annot_kws=akws, cbar=False)
    ax.set_title(f"sc={scmean}; power={power}")


fig = plt.figure(figsize=(8, 6),)
ax = plt.gca()
plot_mean(ax, tck=agenes, scmean=1/5, power=0)

plt.tight_layout()
plt.savefig('/home/cboix/test_ast.png')
plt.close()


scm = [0, 1/12, 1/8, 1/5]
powlist = [0, .1, .25, .5]
w = 3 + 4 * 4
h = 1 + 4 * 4
fig, axs = plt.subplots(4, 4, figsize=(w, h),
                        gridspec_kw={"hspace": 0.025, "wspace": 0.025})

for i, scmean in enumerate(scm):
    for j, power in enumerate(powlist):
        print(i, j)
        plot_mean(axs[i, j], tck=False, scmean=scmean, power=power)

# plt.tight_layout()
plt.savefig('/home/cboix/test_ast.png')
plt.close()


# Make graph with this simple adjustment
# -------------------------------------------
xmean = mod.adata.var['xmean'].to_numpy()
xmean = xmean - np.mean(xmean)
xcol = xmean / np.sqrt(np.sum(xmean**2))
V = mod.cobj.V + xcol[None, :] * 0
s = mod.cobj.s
s_red = np.diag(s**power)
X_cov = V.T.dot(V)
X_sd = np.sqrt(np.diag(X_cov))
cv = X_cov / X_sd[:, np.newaxis]
corr_est = cv / X_sd[np.newaxis, :]
mod.cobj.corr_est = corr_est

# Pretty bad, actually!? -- unclear.... need actual benchmarks
graph_id = 'mean'
resolution = 2
mod.make_graph(graph_id, resolution=resolution, z=4)

mod.graphs[graph_id].adj.zcut
mod.find_gene(graph_id, 'SLC1A2')
mod.plot_graph(graph_id, attr="leiden", show_labels=True, width=16)


mod.plot_gene_umap(graph_id, width=16)
if 'X_umap' in mod.adata.obsm:
    mod.plot_umap_grid(graph_id)


# Make graphs from the zscored and from base:
# -------------------------------------------
graph_id = 'base'
resolution = 2
k_ind = np.arange(1, top_k)
mod.make_svd_subset_graph(graph_id, k_ind=k_ind, resolution=resolution)

g1 = mod.graphs[graph_id].graph  # Almost full size, some get pruned
a1 = mod.graphs[graph_id].adj
[x for x in agenes if x in a1.labels]

graph_id = 'p0.5'
resolution = 2
k_ind = np.arange(1, top_k)
mod.make_svd_subset_graph(graph_id, k_ind=k_ind,
                          resolution=resolution, power=0.5)
mod.plot_graph(graph_id, attr="leiden", show_labels=True, width=16)

g1 = mod.graphs[graph_id].graph  # Almost full size, some get pruned
a1 = mod.graphs[graph_id].adj
[x for x in agenes if x in a1.labels]


# Make graphs with multiplexed calling:
# -------------------------------------
graph_id = 'merge'
plist = [0, .5, 1]
plist = list(np.arange(0, 1, .1)) + [1]
mod.z = 4.5
mod.make_merged_graph(graph_id, power_list=plist, resolution=2.5)
# TODO: Make recluster graph possible for merged graph:
# mod.recluster_graph(graph_id,...)

mlist = mod.get_modules(graph_id, print_modules=True)

# 74 clusters, 1873 nodes...? with res = 3 (NOTE: Remove 1-size clusters?)
mod.find_gene(graph_id, 'CD44')
mod.find_gene(graph_id, 'GFAP')
mod.find_gene(graph_id, 'OSMR')
mod.find_gene(graph_id, 'MEIS2')
mod.find_gene(graph_id, 'SLC1A2')
mod.find_gene(graph_id, 'SLC1A3')
mod.find_gene(graph_id, 'VEGFA')
mod.find_gene(graph_id, 'FAU')

# Two views of the modules:
# mod.plot_graph(graph_id, attr="leiden", show_labels=True, width=16)
mod.plot_graph(graph_id, attr="leiden", show_labels=False, width=16)
mod.plot_gene_umap(graph_id, attr="leiden", width=16)

# Modules on cell UMAP:
# mod.plot_umap_grid(graph_id)

# Get the modules/print out:
mlist = mod.get_modules(graph_id, print_modules=False)
# mod.save_modules(graph_id)

# Get functional enrichments for the modules:
gpres = mod.get_goterms(graph_id)
# TODO: Save GO terms?
