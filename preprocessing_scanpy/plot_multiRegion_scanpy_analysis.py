# !/usr/bin/python
# ----------------------------------------------
# Repeat the multi-region analysis using scanpy:
# ----------------------------------------------
import glob
import h5py
from scipy import sparse
import re
from sklearn.utils import sparsefuncs
# from sklearn.decomposition import IncrementalPCA
# from inc_pca import IncPCA
import gzip
import pickle
import numpy as np
import pandas as pd
import time
import gc
import fbpca
import sys

# For KNN + conn
import umap
from umap.umap_ import nearest_neighbors
from sklearn.utils import check_random_state
from types import MappingProxyType
from umap.umap_ import fuzzy_simplicial_set

# For clustering:
import igraph as ig
import leidenalg
import hdbscan

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

import scanpy.api as sc
import anndata


# Scanpy plotting settings:
sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving

# ----------
# Functions:
# ----------
# Faster than for loop or np.savetxt:
def preformatted_write(mat, f, fmtstring=None, encode=False):
    if fmtstring is None:
        if len(mat.shape) == 1:
            fmtstring = '%g'
        else:
            fmtstring = '\t'.join(['%g']*mat.shape[1])
    fmt = '\n'.join([fmtstring]*mat.shape[0])
    data = fmt % tuple(mat.ravel())
    if encode:
        data = data.encode('utf-8')
    f.write(data)

def gzipped_write(matrix, filename):
    with gzip.open(filename, 'wb') as f:
        preformatted_write(matrix, f, encode=True)

def plot_dimred(mat, filename, redtype=None, title=None,
                values=None, cmap=None, size=.1):
    sns.set(font_scale=1.1)
    # Plot current trace:
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    # Plot % change:
    if title is not None:
        ax.set_title(title)
    ax.set_facecolor('white')
    if values is None:
        plt.scatter(mat[:,0], mat[:,1], s=size)
    else:
        plt.scatter(mat[:,0], mat[:,1], s=size,
                    c=values, cmap=cmap)
    if redtype is not None:
        plt.ylabel(redtype + ' 1')
        plt.xlabel(redtype + ' 2')
    plt.tight_layout()
    fig = plt.gcf()
    fig.savefig(filename, dpi=350, bbox_inches='tight')

def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return(matrix)

# For labels write:
def aggregate_labels_over_matrix(matrix, labels, hdb=False):
    NLAB = np.max(labels) + 1 + 1 * hdb
    avgmat = np.zeros((NLAB, matrix.shape[1]))
    for i in range(NLAB):
        print(i)
        ind = np.where(labels == (i - 1 * hdb))[0]
        avgmat[i,:] = np.array(np.mean(matrix[ind,:], axis=0)[0])[0]
    return(avgmat)

def calc_write_avgmatrix(lblset, hdb=False):
    groupfile = prefix + "." + lblset + '.tsv.gz'
    avgfile = prefix + "." + lblset + '.avg.tsv.gz'
    labels = pd.read_csv(groupfile, header=None, sep="\t").to_numpy().T[0]
    avgmat = aggregate_labels_over_matrix(fullmat, labels, hdb=hdb)
    gzipped_write(avgmat, avgfile)

# ------------------
# Load the datasets:
# ------------------
keptbc = {}
with open('../Annotation/multiRegion_rows.txt','r') as f:
    for line in f:
        line = line.rstrip().split(".")
        if line[0] in keptbc.keys():
            keptbc[line[0]].append(line[1])
        else:
            keptbc[line[0]] = [line[1]]

suffix = "filtered_feature_bc_matrix.h5"
fnames = np.sort(glob.glob("*" + suffix))

# Load in the protein coding genes to subset the data:
anno = pd.read_csv('../Annotation/Gene.v28lift37.annotation.bed',
                   header=None, sep="\t",
                   names=['chr','start','end','strand','ENSG','type','symbol'])
keep_type = ['protein_coding']
pc_symbols = list(anno.symbol[[tp in keep_type for tp in anno['type']]])
pc_ensg = list(anno.ENSG[[tp in keep_type for tp in anno['type']]])
#pc_genes = list(anno.symbol)

# load_data = False
load_data = True

ncells = 0
acell = []
agene = []
# h5file = fnames[0]
fullmat = None
sind = None
allbc = []


# Process anndata for one file:
def read_process_h5(h5file, suffix=suffix):
    region = re.sub("_" + suffix,"", h5file)
    adata = sc.read_10x_h5(h5file)
    print("Finished reading file into adata")
    adata.var_names_make_unique()
    # Subset to only kept individuals:
    bc = list(adata.obs_names)
    keepind = []
    for kb in keptbc[region]:
        keepind = keepind + [i for i,v in enumerate(bc)
                                if re.search('-' + kb + '$',v)]
    keepind = np.sort(np.array(keepind))
    ind = np.array([False] * len(bc))
    ind[keepind] = True
    # Slicing affects the X matrix weirdly, so work on matrices:
    X = adata.X[ind,:]
    vnames = np.array(adata.var_names)
    onames = np.array(adata.obs_names[ind])
    del(adata)
    # Remove cells with greater than 5% ribo and 20% mito:
    ribo_genes = [name for name in list(vnames) if re.match('^RP[0-9]+-', str(name)) or re.match('^RP[SL]', name)]
    mito_genes = [name for name in list(vnames) if re.match('^MT-', name) or re.match('^MTRNR', name)]
    # Calculate mito and ribo pct:
    rind = [np.where(vnames == n)[0][0] for n in ribo_genes]
    mind = [np.where(vnames == n)[0][0] for n in mito_genes]
    gmarg = np.array(np.sum(X, axis=1).T[0])[0]
    rmarg = np.array(np.sum(X[:, rind], axis=1).T[0])[0]
    mmarg = np.array(np.sum(X[:, mind], axis=1).T[0])[0]
    # Percentage:
    rfrac = rmarg / gmarg
    mfrac = mmarg / gmarg
    # Remove cells with greater than 5% ribo and 20% mito:
    mrkeep = ((mfrac < .2) * (rfrac < .05))
    X = X[mrkeep,:]
    onames = onames[mrkeep]
    # Keep only PC genes:
    in_pc = np.array([name in pc_symbols for name in vnames])
    X = X[:, in_pc]
    vnames = vnames[in_pc]
    return(X, vnames, onames)

load_data = True
load_data = False

# ----------------------------------
# Import all matrices one at a time:
# ----------------------------------
if load_data:
    t1 = time.time()
    fullmat = None
    allbc = []
    for filename in fnames:
        X, vnames, onames = read_process_h5(filename)
        region = re.sub("_" + suffix,"", filename)
        allbc = allbc + [region + "_" + b for b in list(onames)]
        if fullmat is None:
            fullmat = X
        else:
            fullmat = sparse.vstack([fullmat, X])
        del(X)
    print(time.time() - t1)
    gc.collect()

# ---------------------------------------
# Preprocess + create the AnnData object:
# ---------------------------------------
prefix = 'all_brain_regions_filt_preprocessed_scanpy'
norm_mat = True
# norm_mat = False
if norm_mat:
    prefix = prefix + "_norm"  # update prefix
    if load_data:
        # Normalize per cell to the median # counts:
        gmarg = np.array(np.sum(fullmat, axis=1).T[0])[0]
        cafter = np.median(gmarg)
        gmarg = gmarg / cafter
        sparsefuncs.inplace_row_scale(fullmat, 1/gmarg)
        newmarg = np.array(np.sum(fullmat, axis=1).T[0])[0]
        print(newmarg[0:5])
        # As seurat, log1p transform the data (inplace)
        np.log1p(fullmat.data, out=fullmat.data)

if load_data:
    d = {}
    d['obs'] = allbc
    d['var'] = vnames
    d['X'] = fullmat
    adata = anndata.AnnData(**d)
    adata.obs_names = allbc
    adata.var_names = vnames
    adata.var.columns = ['varnames']
    adata.obs.columns = ['obsnames']
    gc.collect()


if load_data:
    # Save the full matrix:
    h5ad_file = '/broad/compbio_ce/cboix/multiRegion/matrices/' + \
        prefix + '_fullmatrix.h5ad'
    adata.write(h5ad_file)
    marg_file = '/broad/compbio_ce/cboix/multiRegion/matrices/' + \
        prefix + '_fullmatrix_margin.tsv.gz'
    gmarg = np.array(np.sum(fullmat, axis=1).T[0])[0]
    gzipped_write(gmarg, marg_file)
else:
    # # Read dataset, pre-processed:
    pass
    # h5file = './write/scanpy_preprocessed.h5ad'
    # adata = sc.read_h5ad(h5file)

# For loading the normalized full matrix:
# h5ad_file = './write/' + prefix + '_fullmatrix.h5ad'
# mdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
mdir = './write/'
h5ad_file = mdir + prefix + '_fullmatrix.h5ad'
t1 = time.time()
adata = sc.read_h5ad(h5ad_file)
time.time() - t1

# Read metadata for potential filtering:
metadir = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/metadata/'
metafile = metadir + 'all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv.gz'
metadf = pd.read_csv(metafile, sep="\t")
# Kept barcodes:
kbc = metadf.barcode.to_numpy()

# Reduce to kept barcodes only:
filter_mat = True
if filter_mat:
    prefix = prefix + "_filtered"
    if load_data:
        adata = adata[kbc,:]

if not filter_mat:
    with open(prefix + '_genes.txt', 'w') as f:
        for i in range(len(vnames)):
            out = f.write(vnames[i] + '\n')

# ---------------------
# Pre-process the data:
# ---------------------
# NOTE: ONLY DO IF NOT FILTERED.
from scanpy.preprocessing._deprecated.highly_variable_genes import filter_genes_dispersion
from scanpy.preprocessing import normalize_total

if load_data:
    # only consider genes with more than 10 counts
    # sc.pp.filter_genes(adata, min_counts=10)
    # normalize with total UMI count per cell
    # normalize_total(adata, key_added='n_counts_all')
    n_top_genes = 5000
    filter_result = filter_genes_dispersion(adata.X, flavor='cell_ranger',
                                            n_top_genes=n_top_genes, log=True)
    filtgenes = adata.var_names[filter_result.gene_subset]
    filtgenes
    adata = adata[:, filtgenes]
    gc.collect()
    # NOTE: log1p scaling not done originally with first UMAP run
    sc.pp.log1p(adata)  # log transform: X = log(X + 1) # Already done before!?
    gc.collect()
    # sc.pp.scale(adata) # NOTE CANNOT ACTUALLY SCALE - too big.
    # gc.collect()

# PCA OPTION:
use_rpca = False
# use_rpca = True
k = 50
t1 = time.time()
if use_rpca:
    # Get the fbpca result:
    prefix = prefix + '_rpca'

pcafile = prefix + '.X_pca.tsv.gz'

try:
    X_pca = pd.read_csv(pcafile, header=None, sep="\t").to_numpy()
    adata.obsm['X_pca'] = X_pca
except:
    if use_rpca:
        # Get the fbpca result:
        print("[STATUS] Performing pca with fbpca.pca and k=" + str(k))
        n_iter = 8
        (U, s, Va) = fbpca.pca(adata.X, k=k, raw=False, n_iter=n_iter, l=None)
        sdiag = np.diag(s)
        US = U.dot(sdiag)
        adata.obsm['X_pca'] = US
        gzipped_write(US, pcafile)
    else:
        # Use the scipy / scanpy:
        sc.pp.pca(adata, n_comps=k, zero_center=False)
        gzipped_write(adata.obsm['X_pca'], pcafile)

print('[STATUS] Finished in ' + str(round(time.time() - t1, 2)) + "s")

# run_pp = False
# if run_pp:
#     n_top_genes = 5000
#     filter_result = filter_genes_dispersion(adata.X, flavor='cell_ranger',
#                                             n_top_genes=n_top_genes, log=False)
#     # actually filter the genes, the following is the inplace version of
#     #     adata = adata[:, filter_result.gene_subset]
#     adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
#     normalize_total(adata)  # renormalize after filtering
#     sc.pp.log1p(adata)  # log transform: X = log(X + 1)
#     sc.pp.scale(adata) # TODO can we do this with sparse mat?

# ---------------------------------------------------
# Simple scanpy recipe after our preprocessing + PCA:
# ---------------------------------------------------
h5ad_file = './write/' + prefix + '.h5ad'
if not filter_mat:
    adata = sc.read_h5ad(h5ad_file)
    gc.collect()
    # adata.write(h5ad_file) # Save if first time

barcodes = adata.obs_names
with gzip.open(prefix + '.barcodes.tsv.gz','wb') as f:
    for line in barcodes:
        out = line + '\n'
        a = f.write(out.encode())


t1 = time.time()
sc.pp.neighbors(adata)
print('Neighbors',time.time() - t1)

# sc.tl.louvain(adata)
# niter = 50
niter = 100

#RES = 15
prefstr = '_r' + str(RES) + '_n' + str(niter)
sc.tl.leiden(adata, n_iterations=niter, resolution=RES)
adata.write(h5ad_file)
labels = adata.obs['leiden'].to_numpy().astype(int)
gzipped_write(labels, prefix + '.leiden' + prefstr + '.tsv.gz')
print('Leiden',time.time() - t1)

# sc.tl.paga(adata)
# sc.pl.paga(adata) # Plot
# adata.write(h5ad_file)
# print('PAGA',time.time() - t1)

# If umap exists:
# X_umap = pd.read_csv(prefix + '.umap.tsv.gz', header=None, sep="\t").to_numpy()

#R# adata.obsm['X_umap'] = X_umap

sc.tl.umap(adata)
sc.pl.umap(adata) # Plot umap
# sc.pl.umap(adata, color='leiden') # Plot umap
# sc.pl.umap(adata, color='leiden', legend_loc='on data', frameon=False, save=prefstr + '.png')
# adata.write(h5ad_file)
gzipped_write(adata.obsm['X_umap'], prefix + '.umap.tsv.gz')
print('UMAP',time.time() - t1)

sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.rank_genes_groups(adata, save='.pdf') # Plot
adata.write(h5ad_file)
print('Ranking',time.time() - t1)

# Get the average values for these clusters:
def calc_write_avgmatrix(lblset, hdb=False):
    groupfile = prefix + "." + lblset + '.tsv.gz'
    avgfile = prefix + "." + lblset + '.avg.tsv.gz'
    labels = pd.read_csv(groupfile, header=None, sep="\t").to_numpy().T[0]
    avgmat = aggregate_labels_over_matrix(adata.X, labels, hdb=hdb)
    gzipped_write(avgmat, avgfile)

# NOTE: Need to load the full matrix:
lblset = 'leiden' + prefstr
calc_write_avgmatrix(lblset)



