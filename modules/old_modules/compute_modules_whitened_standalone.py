#!usr/bin/python
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework:
# Created: 06/16/21
# -------------------------------------------------------------
import glob
import h5py
import re
import gzip
import numpy as np
import pandas as pd
import time
import gc
import os
import sys

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import rcParams
import seaborn as sns

# Single-cell processing:
import scanpy as sc
import anndata

# -----------------------------------------------------
# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# Datafiles:
ct = 'Ast'
st = 'Ast_GRM3'
ct = 'Mic_Immune'
st = ''
ctstr = re.sub("_","/", ct)
ststr = re.sub("_"," ", st)
dbdir = '/home/cboix/data/DEVTRAJ/db/'
datadir = '/home/cboix/data/DEVTRAJ/db/multiRegion/'
pref = 'all_brain_regions_filt_preprocessed_scanpy'
fh5ad = datadir + 'matrices/' + pref + '.majorcelltype.' + ct + '.hdf5'
metafile = datadir + pref + '_norm.final_noMB.cell_labels.tsv.gz'
snapcols = pd.read_csv(dbdir + '/Annotation/snap_colors.tsv', header=None)

# TODO: replace all of this with a loader.

# Metadata:
mdf = pd.read_csv(metafile, sep='\t')
mdf = mdf.loc[mdf.hcelltype == ctstr,:]
if ststr is not '':
    mdf = mdf.loc[mdf.cell_type_high_resolution == ststr,:]

# Only the EC specific subtypes:
if ct == 'Exc' and st == 'EC':
    ecs = ['Exc AGBL1 GPC5', 'Exc COL25A1 SEMA3D', 'Exc DLC1 SNTG2',
           'Exc RELN COL5A2', 'Exc RELN GPC5', 'Exc SOX11 NCKAP5',
           'Exc TOX3 INO80D', 'Exc TOX3 POSTN', 'Exc TOX3 TTC6']
    mdf = mdf.loc[mdf.region == 'EC',:]
    mdf = mdf.loc[mdf.cell_type_high_resolution.isin(ecs),:]
elif ct == 'Exc' and st == 'ECall':
    mdf = mdf.loc[mdf.region == 'EC',:]
    mdf = mdf.loc[mdf['major.celltype'] == 'Exc',:]

kbc = mdf.barcode.tolist()

# Read in data from hdf5 file:
hf = h5py.File(fh5ad, 'r')
mat = hf['matrix'][:]
bcs = hf['barcodes'][:]
genes = hf['genes'][:]
hf.close()

bcs2 = [x.decode() for x in bcs]
genes2 = [x.decode() for x in genes]

# Make adata object:
d = {}
d['obs'] = bcs2
d['var'] = genes2
d['X'] = mat
adata = anndata.AnnData(**d)
adata.obs_names = bcs2
adata.var_names = genes2
adata.var.columns = ['varnames']
adata.obs.columns = ['obsnames']
gc.collect()

# Annotate with metadata (projid in particular):
adata = adata[kbc,:]
adata.obs['barcode'] = adata.obs_names
mdf.index = mdf.barcode
mdf2 = pd.merge(adata.obs, mdf, how='left', left_index=True, right_index=True)
adata.obs['projid'] = mdf2['projid'].astype('category')
adata.obs['celltype'] = mdf2['cell_type_high_resolution'].astype('category')
reg = [re.sub("_.*","", x) for x in adata.obs.barcode.tolist()]
adata.obs['region'] = reg

# -----------------------------------
# Process this cell type with scanpy:
# -----------------------------------
prefstr = '_test_mrad_' + ct + "_" + st

# Filter cells:
sc.pp.filter_genes(adata, min_cells=3)
adata.obs['cell'] = adata.obs_names
sc.pp.filter_cells(adata, min_genes=100)

# Normalize, log1p:
countsafter = None
sc.pp.normalize_per_cell(adata, counts_per_cell_after=countsafter, copy=False)
sc.pp.log1p(adata)

# Can choose to normalize here:
usecombat = False
if usecombat:
    prefstr = prefstr + "_combat"
    sc.pp.combat(adata, key='projid')

# Here, can filter down to highly variable genes if we want:
usehvg = False
n_top_genes = 10000
if usehvg:
    prefstr = prefstr + "_filthvg"
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=True)
    filtgenes = adata.var_names[filter_result.gene_subset]
    adata = adata[:, filtgenes]

# Filter to top expr:
cutoff = 0.1
cmarg = np.mean(adata.X > 0, axis=0)
filtgenes = adata.var_names[cmarg > cutoff]
adata = adata[:, filtgenes]
print(adata.shape)

# TODO: adjust the weights of the network by the % expressed - run a hurdle model on decorr??

# -----------------------------
# Compute PCA, Neighbors, UMAP:
# -----------------------------
# Whiten the X matrix:
sc.tl.pca(adata) # N = 50, default
V = adata.varm['PCs']
s = np.diag(1 / adata.uns['pca']['variance'])
XtV = adata.X.dot(V)
Xw = XtV.dot(s.dot(V.T)) / V.shape[0]

# Properly ??
U = adata.obsm['X_pca']
XtU = adata.X.T.dot(U)
Xtu = XtU.dot(s.dot(U.T)) / V.shape[0]

# Recompute PCA:
# from sklearn.utils.extmath import randomized_svd
import fbpca as fbpca
U, s, Va = fbpca.pca(adata.X, k=200, raw=False)
s = np.diag(1/s) # TODO: check ok.
XtV2 = adata.X.dot(Va.T)
Xw2 = XtV2.dot(s.dot(Va)) / Va.T.shape[0]

XtU2 = adata.X.T.dot(U)
Xtu2 = XtU2.dot(s.dot(U.T)) / Va.shape[0]

# DOUBLE (doesn't work.)
Xuv = Xtu2.T.dot(Va.T)
Xd  = Xuv.dot(s.dot(Va)) / Va.T.shape[0]

# Compute the correlation matrix between genes (5 ways...):
corr_raw = np.corrcoef(adata.X.T)

corr_w = np.corrcoef(Xw.T) # From adata.
corr_w2 = np.corrcoef(Xw2.T) # From fbpca, centered.

corr_u = np.corrcoef(Xtu) # From adata.
corr_u2 = np.corrcoef(Xtu2) # From fbpca

corr_d = np.corrcoef(Xd.T) # Absolute nonsense

g1 = 'APOE'; g2 = 'CLU'
# g1 = 'HILPDA'; g2 = 'IRS2'
# g1 = 'HILPDA'; g2 = 'GAPDH'
i1 = np.where(adata.var_names == g1)[0][0]
i2 = np.where(adata.var_names == g2)[0][0]
corr_w[i1,i2]
corr_w2[i1,i2]
corr_raw[i1,i2]

corr_u[i1,i2]
corr_u2[i1,i2]

# Compare on UMAP:
import umap
uw = umap.UMAP()
er = uw.fit_transform(corr_raw)
# uw = umap.UMAP()
# ew = uw.fit_transform(corr_w)
uw = umap.UMAP()
ew2 = uw.fit_transform(corr_w2)

# Raw:
plt.figure(figsize=(12,12))
plt.scatter(er[:, 0], er[:, 1])
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP of raw', fontsize=24)
plt.savefig('umap' + prefstr + '_raw.png')
plt.close()

# # Plot each:
# plt.figure(figsize=(12,12))
# plt.scatter(ew[:, 0], ew[:, 1])
# plt.gca().set_aspect('equal', 'datalim')
# plt.title('UMAP of whitened', fontsize=24)
# plt.savefig('umap' + prefstr + '_whiten.png')
# plt.close()

# Plot each:
plt.figure(figsize=(12,12))
plt.scatter(ew2[:, 0], ew2[:, 1])
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP of whitened, centered', fontsize=24)
plt.savefig('umap' + prefstr + '_whiten2.png')
plt.close()

# TODO: Plot properties of celltypes and of genes on UMAP;
# Which CT enriched in
# Number of cells

# ---------------------------------------
# Turn all of the correlations to graphs:
# ---------------------------------------
import igraph
from scipy import sparse
import leidenalg as la
from adjustText import adjust_text

visual_style = {}
visual_style["vertex_size"] = 10
visual_style["vertex_color"] = 'slateblue'
visual_style["vertex_frame_color"] = 'slateblue'
visual_style["vertex_label"] = ''
visual_style["edge_color"] = 'lightgrey'
visual_style["bbox"] = (2400, 2400)
visual_style["margin"] = 10
visual_style["edge_width"] = 0.5
# visual_style["vertex_label"] = g.vs["name"]
# visual_style["edge_width"] = [1 + 2 * int(is_formal) for is_formal in g.es["is_formal"]]

def build_adjacency(corr, cutoff, knn=None, rcut=1):
    aw = corr * (corr > cutoff)
    aw = aw - np.diag(np.diag(aw))
    aw = sparse.csr_matrix(aw)
    rind = np.array((np.sum(aw, axis=0) > rcut))[0]
    # Remove nodes with no links:
    aw = aw[rind, :]
    aw = aw[:, rind]
    return(aw, rind)

def prune_scale(aw, cutoff=0.90, rcut=1):
    # TODO: handle 0s for max:
    awm = np.max(aw, axis=1).T.toarray()[0] # Max link
    aw = aw.tocoo()
    scale = awm[aw.row]
    pct_data = aw.data / scale
    kind = pct_data > cutoff
    aw = sparse.coo_matrix((aw.data[kind], (aw.row[kind], aw.col[kind])), shape=aw.shape)
    aw = aw.tocsr()
    # Remove nodes with no/very few links:
    rind = np.array((np.sum(aw, axis=1) > rcut).T)[0]
    cind = np.array((np.sum(aw, axis=0) > rcut))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return(aw, ind)

def prune_knn(aw, k=50, twodir=True, row_only=True):
    # Prune to maximum k for each node:
    aw = aw.tocoo()
    ind = np.where(aw.row < aw.col)[0]
    ind = ind[np.argsort(-aw.data[ind])]
    mk = np.zeros(aw.shape[0], int)
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
    aw = sparse.coo_matrix((aw.data[kind], (aw.row[kind], aw.col[kind])), shape=aw.shape)
    aw = aw.tocsr()
    # Remove nodes with no/very few links:
    rind = np.array((np.sum(aw, axis=1) > 4).T)[0]
    cind = np.array((np.sum(aw, axis=0) > 4))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return(aw, ind)

# TODO: move prefstr out.
def plot_corr_graph(aw, suff, labels=None, k=None,
                    cutoff=None, visual_style=visual_style):
    # Build graph from edge-list from sparse matrix:
    if k is not None:
        print("pruning knn")
        print(np.sum(aw))
        aw, keptind = prune_knn(aw, k=k)
        if labels is not None:
            labels = labels[keptind]
        print(np.sum(aw))
    if cutoff is not None:
        print("pruning by relative cutoff")
        aw, keptind = prune_scale(aw, cutoff=cutoff)
        if labels is not None:
            labels = labels[keptind]
        print(np.sum(aw))
    aw = aw.tocoo()
    ind = aw.row < aw.col
    edg_list = list(tuple(zip(aw.row[ind], aw.col[ind])))
    edg_wt = aw.data[ind]
    gw = igraph.Graph(edg_list, directed=False)
    gw.vs['name'] = labels
    # Delete after labeling:
    todel = [i for i,x in enumerate(gw.degree()) if x < 1]
    gw.delete_vertices(todel)
    layout = gw.layout('fr')
    # Plot:
    plotname = "graph" + prefstr + '_test_' + suff + '.png'
    visual_style["vertex_color"] = 'slateblue'
    visual_style["vertex_frame_color"] = 'slateblue'
    if labels is not None:
        visual_style["vertex_label"] = gw.vs['name']
    else:
        visual_style["vertex_label"] = ''
    igraph.plot(gw, plotname, layout=layout, **visual_style)
    return(gw, layout)


# Plot leiden:
def plot_corr_graph_leiden(gw, layout, suff, pltnames=False,
                           resolution=None, visual_style=visual_style):
    if resolution is None:
        ptn = la.find_partition(gw, la.ModularityVertexPartition)
        plotname = "graph" + prefstr + '_test_' + suff + '_partition.png'
    else:
        plotname = "graph" + prefstr + '_test_' + suff + '_r' + str(resolution) + '_partition.png'
        ptn = la.find_partition(gw, la.RBConfigurationVertexPartition, resolution_parameter=resolution)
    # Colors:
    ptlist = np.array([np.array(x) for x in list(ptn)])
    nc = len(ptlist); j = 1
    print("Found " + str(nc) + " clusters")
    colors = snapcols[j:j+nc].to_numpy().T[0]
    lassign = (np.zeros(len(gw.vs)) - 1).astype(int)
    for i in range(nc):
        lassign[ptlist[i]] = i
    vscols = colors[lassign]
    # Re-write visual_style colors:
    visual_style["vertex_color"] = vscols
    visual_style["vertex_frame_color"] = vscols
    visual_style["vertex_label"] = ''
    # Plot:
    fig = plt.figure(figsize=(24,24))
    ax = plt.gca()
    igraph.plot(gw, layout=layout, target=ax, **visual_style)
    if pltnames:
        llist = list(layout)
        visual_style["vertex_label"] = gw.vs['name']
        texts = [plt.text(llist[i][0], llist[i][1], gw.vs['name'][i],
                          ha='center', va='center', fontsize=12) for i in range(len(llist))]
        adjust_text(texts, lim=25) # Rough adjustment, limit number of iter.
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    return(gw, vscols)

genes = adata.var_names # labels

suff = 'raw'
aw, keptind = build_adjacency(corr_raw, cutoff=0.14)
gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, k=None)
gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)

# suff = 'whiten1'
# aw, keptind = build_adjacency(corr_w, cutoff=0.75)
# gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, k=None)
# gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)

suff = 'whiten2'
aw, keptind = build_adjacency(corr_w2, cutoff=0.4)
gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, k=None)
gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)
# gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, resolution=1.05)

# suff = 'u1'
# aw, keptind = build_adjacency(corr_u, cutoff=0.92)
# gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, k=None)
# gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)

suff = 'u2'
aw, keptind = build_adjacency(corr_u2, cutoff=0.4, rcut=4)
gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, k=None)
gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)

suff = 'u2_knn'
aw, keptind = build_adjacency(corr_u2, cutoff=0.25, rcut=5)
gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, cutoff=.85, k=5)
gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)


suff = 'double'
aw, keptind = build_adjacency(corr_d, cutoff=0.25)
gw, layout = plot_corr_graph(aw, labels=genes[keptind], suff=suff, k=None)
gw, vscols = plot_corr_graph_leiden(gw, layout, suff=suff, pltnames=True)


# Plot whitened UMAP with the whitened leiden colors:
plt.figure(figsize=(12,12))
plt.scatter(ew2[keptind, 0], ew2[keptind, 1], color=vscols)
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP of whitened, centered', fontsize=24)
plt.savefig("umap" + prefstr + "_whiten2_louvain.png")
plt.close()

# Plot raw UMAP with the whitened leiden colors:
plt.figure(figsize=(12,12))
plt.scatter(er[keptind, 0], er[keptind, 1], color=vscols)
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP of whitened, centered', fontsize=24)
plt.savefig("umap" + prefstr + "_raw_louvain_w2.png")
plt.close()

# Plot corr by marg:
cmarg = np.mean(adata.X > 0, axis=0)
mcw = np.median(corr_w2, axis=0)
for i in range(corr_w2.shape[0]):
    c = corr_w2[i,:]
    c = np.sort(c)
    mcw[i] = np.mean(c[-10:-1])

# Slightly higher corr for low-expressed, due to variability issue.
# TODO: Model expression by summary stats of the genes and adjust corr.

# Plot each:
plt.figure(figsize=(12,12))
plt.scatter(cmarg, mcw, s=5)
# plt.gca().set_aspect('equal', 'datalim')
plt.title('Margin of genes vs. median corr', fontsize=24)
plt.savefig("marg_vs_corr" + prefstr + "_whiten2.png")
plt.close()

cmean = np.mean(adata.X, axis=0)
cvar = np.var(adata.X, axis=0)
disp = cvar / cmean

# Plot each:
plt.figure(figsize=(12,12))
plt.scatter(disp, mcw, s=5)
# plt.gca().set_aspect('equal', 'datalim')
plt.title('Dispersion of genes vs. median corr', fontsize=24)
plt.savefig("disp_vs_corr" + prefstr + "_whiten2.png")
plt.close()

dgenes = ['CLDN5','MT-ND3','RPS13','IL6R']
plt.figure();
fig, axs = plt.subplots(1,len(dgenes), figsize=(len(dgenes) * 3,4));
for i, ax in enumerate(axs):
    j = np.where(genes == dgenes[i])[0][0]
    h = ax.hist(corr_w2[j,:])
    t = ax.set_title(dgenes[i]);

plt.tight_layout()
plt.savefig('hist_examples' + prefstr + '.png')
plt.close()

# Alternatively, least sq:
XtV2 = adata.X.dot(Va.T)
Xw2 = XtV2.dot(s.dot(Va)) / Va.T.shape[0]
Va.T.dot(s.dot(Va))


cholsigmainv = np.linalg.cholesky(np.linalg.inv(np.cov(screens.T)))
warped_screens = screens.values @ cholsigmainv
warped_intercept = cholsigmainv.sum(axis=0)

for gene_index in range(len(warped_screens)):
    X = np.stack((warped_intercept, warped_screens[gene_index]), axis=1)
    coef, residues = np.linalg.lstsq(X, ys, rcond=None)[:2]
    df = warped_screens.shape[1] - 2
    GLS_coef[gene_index] = coef[1]
    GLS_se[gene_index] = \
        np.sqrt(np.linalg.pinv(X.T @ X)[1, 1] * residues / df)

# Plot corr by marg:
cmarg = np.mean(adata.X > 0, axis=0)
mcw = np.median(corr_w2, axis=0)
for i in range(corr_w2.shape[0]):
    c = corr_w2[i,:]
    c = np.sort(c)
    mcw[i] = np.mean(c[-10:-1])




