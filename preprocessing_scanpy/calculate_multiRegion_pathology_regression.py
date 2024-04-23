# !/usr/bin/python
# ----------------------------------------------
# Repeat the multi-region analysis using scanpy:
# ----------------------------------------------
import glob
import h5py
from scipy import sparse
import re
from sklearn.utils import sparsefuncs
from sklearn.decomposition import IncrementalPCA
import gzip
import pickle
from inc_pca import IncPCA
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

# For regression:
from scipy import linalg
from sklearn.linear_model import Lasso, Ridge, ElasticNet
from tqdm import tqdm

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
# Full matrix:
prefix = 'all_brain_regions_filt_preprocessed_scanpy'
prefix = prefix + "_norm"
h5ad_file = './write/' + prefix + '_fullmatrix.h5ad'
adata = sc.read_h5ad(h5ad_file)

# Cell assignments:
path = 'nft'
path = 'plaq_n'
lblset = 'leiden_r5_n50'
lbdf = pd.read_csv(prefix + "." + lblset + ".lbls.tsv",
                   header=None, sep="\t", names=['barcode', 'cluster'])
# Pathology
pathdf = pd.read_csv(prefix + "." + path + ".tsv.gz",
                   header=None, sep="\t", names=['barcode', 'path'])
pathdf = pathdf.merge(lbdf)


# --------------------------------------------------------
# For each specific celltype:
# compute the correlations + regression against pathology:
# --------------------------------------------------------
celltype = 'Microglia'
celltype = 'Oligo'
# celltype = 'Astro'
# celltype = 'Endo'
# celltype = 'OPC'
# celltype = 'Per'

cts = pd.unique(pathdf.cluster)
cts = np.sort([ct for ct in cts if not (re.search('Ex', ct) or re.search('In', ct))])

spathdf = pathdf[pathdf.cluster == celltype]
spathdf = spathdf[(np.isnan(spathdf.path) == False)]
print("Reduced to " + str(spathdf.shape[0]) + " cells for " + celltype)

# Update as ct + path:
if path != 'nft':
    celltype = path + "_" + celltype

# Subset the annotation object and matrix:
sdata = adata[spathdf.barcode,:]
y = spathdf.path.to_numpy()

# 1. Calculate correlation
t1 = time.time()
# Mean center y:
yctr = y - np.mean(y)
ysd = np.std(yctr)
# CSC sparse matrix:
sX = sdata.X.tocsc()
del(sdata)
NR = sX.shape[0]
NC = sX.shape[1]
smean = np.array(sX.sum(0) / NR)[0]

covs = np.zeros(NC)
corrs = np.zeros(NC)
for i in tqdm(range(NC)):
    sdense = sX[:,i].todense() - smean[i]
    covs[i] = (yctr * sdense) / (NR - 1)
    corrs[i] = covs[i] / (ysd * np.std(sdense))

# Top genes:
topind = np.where(corrs > 0.1)
adata.var_names[topind]

# Write out table:
corrdf = pd.DataFrame({'symbol' : adata.var_names,
                       'corrval': corrs})
scdf = corrdf.iloc[np.isnan(corrdf.corrval.to_numpy()) == False]
scdf = scdf.iloc[np.argsort(-scdf.corrval)]
scdf.to_csv(prefix + '.' + lblset + '.corr_' + celltype + '.tsv.gz',
            index=False, sep="\t")
print(time.time() - t1)
gc.collect()


# 2. Perform regression (keep MT and ribo genes - they indicate pathology)
t0 = time.time()
alpha = 1
sparse_lasso = Lasso(alpha=alpha, fit_intercept=False, max_iter=1000)
fit = sparse_lasso.fit(sX, y)
print("Sparse Lasso done in %fs" % (time.time() - t0)) # 27s for 81k in Mic
fit.score(sX, y) # .121 for mic
ypred = fit.predict(sX) # .121 for mic
spathdf['lasso_pred'] = ypred

lcdf = pd.DataFrame({'symbol' : adata.var_names,
                     'coeff': fit.coef_,
                     'fit': 'Lasso'})
lcdf = lcdf.iloc[np.where(lcdf.coeff != 0)]
lcdf = lcdf.iloc[np.argsort(-lcdf.coeff)]


# -----------
# ElasticNet:
# -----------
run_enet = False
if run_enet:
    t0 = time.time()
    alpha = 1
    sparse_enet = ElasticNet(alpha=alpha, l1_ratio=0.5,
                            fit_intercept=False, max_iter=1000)
    fit = sparse_enet.fit(sX, y)
    print("Sparse ElasticNet done in %fs" % (time.time() - t0)) # 121s for 81k in Mic
    fit.score(sX, y) # .168 for mic
    ypred = fit.predict(sX)
    spathdf['enet_pred'] = ypred
    ecdf = pd.DataFrame({'symbol' : adata.var_names,
                        'coeff': fit.coef_,
                        'fit': 'ElasticNet'})
    ecdf = ecdf.iloc[np.where(ecdf.coeff != 0)]
    ecdf = ecdf.iloc[np.argsort(-ecdf.coeff)]
    # All coefficients:
    coeffdf = pd.concat([lcdf, ecdf])
else:
    coeffdf = lcdf

coeffdf.to_csv(prefix + '.' + lblset + '.coeff_' + celltype + '.tsv.gz',
               index=False, sep="\t")

# Final scores:
spathdf.to_csv(prefix + '.' + lblset + '.pathpred_' + celltype + '.tsv.gz',
            index=False, sep="\t")

# TODO: Probably want to plot some of these:
coeff_genes = np.sort(pd.unique(coeffdf.symbol))
geneind = np.array([np.where(adata.var_names == x)[0] for x in coeff_genes]).T[0]

# -------------
# Save as hdf5:
# -------------
array = sX[:, geneind].toarray()
clevel = 4
h5file = prefix + '.' + lblset + '.topgenes_' + celltype + ".hdf5"
hf = h5py.File(h5file, 'w')
hf.create_dataset('matrix', data=array,
                  compression='gzip', compression_opts=clevel)
encgenes = [n.encode() for n in coeff_genes]
hf.create_dataset('genes', data=encgenes,
                  compression='gzip', compression_opts=clevel)
encbc = [n.encode() for n in spathdf.barcode.to_numpy()]
hf.create_dataset('barcodes', data=encbc,
                  compression='gzip', compression_opts=clevel)
hf.close()
del(array)



# ------------------
# Save full as hdf5:
# ------------------
array = sX.toarray()
clevel = 4
h5file = prefix + '.' + lblset + '.full_' + celltype + ".hdf5"
hf = h5py.File(h5file, 'w')
hf.create_dataset('matrix', data=array,
                  compression='gzip', compression_opts=clevel)
encgenes = [n.encode() for n in adata.var_names]
hf.create_dataset('genes', data=encgenes,
                  compression='gzip', compression_opts=clevel)
encbc = [n.encode() for n in spathdf.barcode.to_numpy()]
hf.create_dataset('barcodes', data=encbc,
                  compression='gzip', compression_opts=clevel)
hf.close()
del(array)





# ----------------------------
# Plot the array as a heatmap:
# ----------------------------
array = sX.toarray()
pord = np.argsort(spathdf.path.to_numpy())
cord = np.argsort(corrdf.corrval.to_numpy())
# Reorder the array:
array = array[pord[:, np.newaxis] , cord]


# Plot the full array as a heatmap:
figname = prefix + '.' + lblset + '.full_' + celltype + ".heatmap.png"
sns.set(font_scale=1.1)
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_title(celltype)
ax.set_facecolor('white')
# Important options for heatmap:
ax.axis('off')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_aspect('equal')
cmap = 'viridis'
plt.imshow(array, cmap=cmap, vmin=0, vmax=5)
plt.xlabel('Genes')
plt.ylabel('Cells')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(figname, dpi=350, bbox_inches='tight')


# TODO: Reduce the plotted genes, and split by brain region?

