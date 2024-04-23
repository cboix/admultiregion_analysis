# !/usr/bin/python
# ----------------------------------------------
# Quantify the amount of ambient RNA in each batch
# In particular check for PLP1 cross contam
# Load in celltype annotations to look for incorrect cross-expr
# Use SoupX method (EM)
# ----------------------------------------------
import glob
import h5py
from scipy import sparse
import re
from sklearn.utils import sparsefuncs
import gzip
import pickle
import numpy as np
import pandas as pd
import time
import gc
import fbpca
import sys
from scipy.stats import ttest_ind  # Matrix vs. matrix ttest

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

# ----------------------------------
# Load in the metadata for grouping:
# ----------------------------------
metadir = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/metadata/'
metafile = metadir + 'all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv.gz'
metadf = pd.read_csv(metafile, sep="\t")
# Locations of celltype matrices:
mdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy'

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

# load_data = False
load_data = True
ncells = 0
acell = []
agene = []
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


# ----------------------------------
# Start with analysis of one region:
# ----------------------------------
filename = fnames[1]
print(region, filename)
region = re.sub("_" + suffix,"", filename)
X, vnames, onames = read_process_h5(filename)
regbc = [region + "_" + b for b in list(onames)]

# Make temporary anndata object:
d = {}
d['obs'] = regbc
d['var'] = vnames
d['X'] = X
adata = anndata.AnnData(**d)
adata.obs_names = regbc
adata.var_names = vnames
adata.var.columns = ['varnames']
adata.obs.columns = ['barcode']
adata.obs.index.names = ['index']

# Note adata must be first - othw disorganized.
adata.obs['cellid'] = mdf['label'].astype('category')
adata.obs['genotype'] = mdf['genotype'].astype('str').astype('category')
adata.obs['celltype'] = mdf['celltype'].astype('category')
adata.obs['mouse_id'] = mdf['mouse_id'].astype('category')
adata.obs['timepoint'] = mdf['timepoint'].astype('category')
stg = pd.read_csv(path + 'CKP25/samp_with_transgene.txt', header=None)[0].tolist()
adata.obs['p25dummy'] = [x in stg for x in adata.obs['cell']]

# Subset to the correct barcodes, add metadata:
submeta = metadf.loc[metadf.region == region,:]
submeta.index = submeta['barcode']
submeta.index.names = ['index']
keepbc = submeta.barcode
adata = adata[keepbc,:]
mdf = pd.merge(adata.obs, submeta, how='left',
               left_index=True, right_index=True)
adata.obs['celltype'] = mdf['hcelltype'].astype('category')
adata.obs['cellhigh'] = submeta['cell_type_high_resolution'].astype('category')
adata.obs['rind'] = submeta['rind'].astype('category')


# --------------------------------------------
# Rough pipeline for estimating contamination:
# --------------------------------------------
# 1. Based on hcelltype, find clear marker genes for each cell type:
# --------------------------------------------
NGENES = 50
sc.tl.rank_genes_groups(adata, groupby='celltype', use_raw=True,
                        method='t-test_overestim_var', n_genes=NGENES)

markers = adata.uns['rank_genes_groups']['names']
cts = adata.obs.celltype.unique()
mgenes = np.array([markers[ct].tolist() for ct in cts]).flatten()

# 2. Assign 0-genes to the "wrong" cell types:
# --------------------------------------------
gind = np.nonzero(np.in1d(vnames, mgenes))
gn = vnames[gind]

# Calc quantile expression per ct:
subad = adata[:,mgenes.tolist()]
arr = subad.X.toarray()

qarr = []
for ct in adata.obs['celltype'].unique():
    cind = adata.obs['celltype'] == ct
    print(ct, sum(cind), 'cells')
    qarr.append(np.quantile(arr[cind,:], q=0.9, axis=0).tolist())

# Take locations where 90% of similar cells have no reads:
# NOTE: This approach may not work on MB region!
qmat = np.array(qarr)
del(arr, subad)
gc.collect()

# 3. Use the 0-genes to estimate the level of contam in each cell
# --------------------------------------------
# Calculate the total # reads in each cell:
marg = np.array(np.sum(adata.X, axis=1)).T[0]

# Per-"0-gene" estimates of contam - averaged in the cell
cfrac = np.zeros(adata.shape[0])
for i, ct in enumerate(adata.obs['celltype'].unique()):
    cind = adata.obs['celltype'] == ct
    print(i, ct, sum(cind), 'cells')
    # Keep some genes:
    nonexpr = mgenes[qmat[i,:] == 0]
    subad = adata[cind,nonexpr]
    frac0 = np.array(np.sum(subad.X,axis=1)).T[0]
    cfrac[cind] = frac0 / marg[cind]

100 * np.quantile(cfrac, q=[0,0.05,0.25,0.5,0.75,0.95,1])
# Plot histogram:
sn

# NOTE: CONTAM DEPENDS ON
# (1) GENE's expression in major ct and (2) fraction
# major celltype in the of the original data

# See if certain individuals have higher contam:
imeans = []
for i,pid in enumerate(adata.obs['rind'].unique()):
    cind = adata.obs['rind'] == pid
    imeans.append(np.mean(cfrac[cind]))



# 4. Based on level of contam, estimate the overall baseline for each individ.
# --------------------------------------------


# 5. Repeat this procedure/E-M (?)
# --------------------------------------------
# Non-neg fit of the baseline --> re-est contam --> etc.


# 6. Look at specific genes that are contam
# --------------------------------------------


# 7. Look at the differential genes that are contam in some individuals but not
# others
# --------------------------------------------


# -------------------------------------------------------
# Alternate pipeline using the cell type gene expression:
# -------------------------------------------------------
cellfracs = []
cellsigs = []
ambexpr = np.zeros(adata.shape[1])
for i, ct in enumerate(adata.obs['celltype'].unique()):
    cind = adata.obs['celltype'] == ct
    print(i, ct, sum(cind), 'cells')
    cellfracs.append(sum(cind)/len(cind) * 1.0)
    # Keep some genes:
    subad = adata[cind,:]
    subcounts = np.array(np.sum(subad.X,axis=0))[0]
    cellsigs.append(subcounts / np.sum(subcounts))
    ambexpr = ambexpr + cellsigs[i] * cellfracs[i]

sgenes = vnames[np.argsort(-ambexpr)]


# Fit with NNLS:
from scipy.optimize import nnls
arr = adata.X.toarray()

fit = nnls(arr.T, ambexpr)








# # Quantify amount of each marker gene in "wrong" population:
# micbc = submeta.barcode[submeta.hcluster == 'Mic'].tolist()
# olibc = submeta.barcode[submeta.hcelltype == 'Oli'].tolist()
# astbc = submeta.barcode[submeta.hcelltype == 'Ast'].tolist()
# excbc = submeta.barcode[submeta.hcelltype == 'Exc'].tolist()

# # PLP1 as example:
# i = np.where(vnames == 'PLP1')[0]
# i = np.where(vnames == 'GFAP')[0]
# i = np.where(vnames == 'CSF1R')[0]
# i = np.where(vnames == 'LINGO1')[0]
# mind = np.nonzero(np.in1d(regbc, micbc))[0]
# oind = np.nonzero(np.in1d(regbc, olibc))[0]
# aind = np.nonzero(np.in1d(regbc, astbc))[0]
# eind = np.nonzero(np.in1d(regbc, excbc))[0]

# np.quantile(np.array(X[mind,i]),q=[0.05,0.25,0.5,0.75,0.95,1])
# np.quantile(np.array(X[oind,i]),q=[0.05,0.25,0.5,0.75,0.95,1])
# np.quantile(np.array(X[aind,i]),q=[0.05,0.25,0.5,0.75,0.95,1])
# np.quantile(np.array(X[eind,i]),q=[0.05,0.25,0.5,0.75,0.95,1])

# qmat = np.quantile(X,q=[0.05,0.25,0.5,0.75,0.95,1])
# qmat = np.quantile(X,q=[0.05,0.25,0.5,0.75,0.95,1], axis=1)
# qmat = np.quantile(a=arr,q=[.5], axis=1)


