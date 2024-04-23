# ------------------------------------------------
# Processing + plotting of specific cell datasets:
# Created: 04/04/21
# ------------------------------------------------
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


# -----------------------------------------------------
# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# Datafiles:
ct = 'Mic_Immune'
st = ''
region = 'EC'
prefstr = '_test_deg_' + ct + "_" + st + "_" + region

ctstr = re.sub("_","/", ct)
ststr = re.sub("_","/", st)
datadir = '/home/cboix/data/DEVTRAJ/db/multiRegion/'
pref = 'all_brain_regions_filt_preprocessed_scanpy'
fh5ad = datadir + 'matrices/' + pref + '.majorcelltype.' + ct + '.hdf5'
metafile = datadir + pref + '_norm.final_noMB.cell_labels.tsv.gz'

# Metadata:
mdf = pd.read_csv(metafile, sep='\t')
mdf = mdf.loc[mdf['major.celltype'] == ctstr,:]
if ct == 'Vasc_Epithelia' and (ststr is not ''):
    mdf = mdf.loc[mdf.cell_type_high_resolution == ststr,:]

if region is not None:
    mdf = mdf.loc[mdf.region == region,:]

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
pathlist = ['nft','plaq_d','plaq_n']

# mmetafile = datadir + 'regions_metadata_list.tsv.gz'
rmetafile = datadir + 'multiregion_path_by_region_projid.tsv'
rdf = pd.read_csv(rmetafile, sep='\t')
rdf2 = pd.merge(adata.obs, rdf, how='left', left_index=False, right_index=False)
rdf2.index = rdf2.barcode
adata.obs['nft'] = np.log1p(rdf2['nft_act'].astype('float'))
adata.obs['plaq_d'] = np.log1p(rdf2['plaq_d_act'].astype('float'))
adata.obs['plaq_n'] = np.log1p(rdf2['plaq_n_act'].astype('float'))

adata.obs['bynft'] = adata.obs['nft'] > 2
adata.obs['bynft'] = adata.obs['bynft'].astype('category')
sum(adata.obs['bynft'])

# ------------------
# Save, for testing:
# ------------------
ncell = np.sum(adata.X > 0, 0)
pctcell = ncell / adata.shape[0]
keepgenes = adata.var_names[pctcell > 0.2].tolist()
subadata = adata[:,keepgenes]
# Writing:
adf = pd.DataFrame(subadata.obs)
adf.to_csv('metadata' + prefstr + '.tsv', sep="\t")
adf = pd.DataFrame(subadata.var)
adf.to_csv('genes' + prefstr + '.tsv', sep="\t")
gzipped_write(subadata.X, 'matrix' + prefstr + '.tsv.gz')

# --------------------------------------
# DE genes by various methods in python:
# --------------------------------------
# sc.tl.rank_genes_groups(adata, 'nft', method='t-test_overestim_var', key_added = "t-test_ov")
NG = 50
sc.tl.rank_genes_groups(adata, 'bynft')
sc.pl.rank_genes_groups(adata, n_genes=NG, sharey=False,
                        frameon=False, save=prefstr + '_bynft_ctalone.png')


# Poisson Mixed Model with statsmodels library:
from statsmodels.genmod.bayes_mixed_glm import PoissonBayesMixedGLM as PMM

def one_hot(array):
    unique, inverse = np.unique(array, return_inverse=True)
    onehot = np.eye(unique.shape[0])[inverse]
    return onehot

endog = 1 * np.array(adata.obs['bynft'].tolist())
# Data + intercept:
exog = np.hstack([np.array([1] * adata.shape[0])[:,np.newaxis], adata.X])
# exog_vc
exog_vc = one_hot(np.array(adata.obs['projid'].tolist()))

model = PMM(endog, exog, exog_vc, ident=[0] * exog_vc.shape[1], vcp_p=1, fe_p=2,
            fep_names=None, vcp_names=None, vc_names=None)

result = model.fit_map() # Laplace


random = {"a": '0 + C(Village)',
          "b": '0 + C(Village)*year_cen'}
model = PoissonBayesMixedGLM.from_formula(
                'y ~ year_cen', random, data)
result = model.fit_vb()

