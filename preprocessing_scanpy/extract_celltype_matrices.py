# !/usr/bin/python
# -----------------------------------------
# Extract the celltype read count matrices:
# -----------------------------------------
import glob
import h5py
from scipy import sparse
import re
import gzip
import pickle
import numpy as np
import pandas as pd
import time
import gc
import sys

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

import scanpy as sc
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



# ----------------------------------
# Load in the metadata for grouping:
# ----------------------------------
metadir = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/metadata/'
metafile = metadir + 'all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv.gz'
metadf = pd.read_csv(metafile, sep="\t")

# ---------------------------------------------
# Load the read count (or normalized) datasets:
# ---------------------------------------------
mdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy'
# marg_file = mdir + prefix + '_fullmatrix_margin.tsv.gz'
# marg = pd.read_csv(marg_file, header=None, sep="\t").to_numpy().T[0]

use_norm = True
if use_norm:
    prefix = prefix + "_norm"
    h5ad_file =  './write/' + prefix + '_fullmatrix.h5ad'
else:
    h5ad_file =  mdir + prefix + '_fullmatrix.h5ad'

adata = sc.read_h5ad(h5ad_file)

# Kept barcodes:
kbc = metadf.barcode.to_numpy()

# Separate by major celltype:
ctypes = pd.unique(metadf['major.celltype'])

for celltype in ctypes:
    cellstr = re.sub("/","_",celltype)
    print(celltype)
    submeta = metadf.loc[metadf['major.celltype'] == celltype,:]
    print(submeta.shape[0])
    # Kept barcodes:
    kbc = submeta.barcode.to_numpy()
    # Subset data:
    subad = adata[kbc,:]
    # Make array:
    array = subad.X.toarray()
    clevel = 4
    h5file = mdir + prefix + '.majorcelltype.' + cellstr + ".hdf5"
    hf = h5py.File(h5file, 'w')
    hf.create_dataset('matrix', data=array,
                      compression='gzip', compression_opts=clevel)
    encgenes = [n.encode() for n in adata.var_names]
    hf.create_dataset('genes', data=encgenes,
                      compression='gzip', compression_opts=clevel)
    encbc = [n.encode() for n in subad.obs_names]
    hf.create_dataset('barcodes', data=encbc,
                    compression='gzip', compression_opts=clevel)
    hf.close()
    del(array)
    gc.collect()

# If these matrices already exist, make them the SWMR mode for concurrent reads

# Find all celltype matrices:
mdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy'
mpref = mdir + prefix + ".majorcelltype."
fnames = np.sort(glob.glob(mpref + "*" + '.hdf5'))

celltype = re.sub('.hdf5', '', re.sub(mpref,"",fn))
cellstr = re.sub("/","_",celltype)
clevel = 4

# # Read in the hdf5 file:
# hf = h5py.File(fn, 'r')
# array = hf['matrix'][]
# hf.close()
# hf.create_dataset('matrix', data=array,
#                     compression='gzip', compression_opts=clevel)
# encgenes = [n.encode() for n in adata.var_names]
# hf.create_dataset('genes', data=encgenes,
#                     compression='gzip', compression_opts=clevel)
# encbc = [n.encode() for n in subad.obs_names]
# hf.create_dataset('barcodes', data=encbc,
#                     compression='gzip', compression_opts=clevel)
# hf.close()

for fn in fnames:
    cellstr = re.sub("/","_",celltype)
    print(celltype)
    submeta = metadf.loc[metadf['major.celltype'] == celltype,:]
    print(submeta.shape[0])
    # Kept barcodes:
    kbc = submeta.barcode.to_numpy()
    # Subset data:
    subad = adata[kbc,:]
    # Make array:
    array = subad.X.toarray()
    clevel = 4
    h5file = mdir + prefix + '.majorcelltype.' + cellstr + ".hdf5"
    hf = h5py.File(h5file, 'w')
    hf.create_dataset('matrix', data=array,
                    compression='gzip', compression_opts=clevel)
    encgenes = [n.encode() for n in adata.var_names]
    hf.create_dataset('genes', data=encgenes,
                    compression='gzip', compression_opts=clevel)
    encbc = [n.encode() for n in subad.obs_names]
    hf.create_dataset('barcodes', data=encbc,
                    compression='gzip', compression_opts=clevel)
    hf.close()
    del(array)
    gc.collect()


# -------------------------------------
# Save the strained EC as an HDF5 file:
# -------------------------------------
ecad_file = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/' + \
    'EC_strainedCounts_raw.h5ad'
h5file = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/' + \
    'EC_strainedCounts_raw.hdf5'
adata = sc.read_h5ad(ecad_file)

# Write margin:
marg = np.array(np.sum(adata.X, axis=1).T[0])[0]
marg_file = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/' + \
    'EC_strainedCounts_margin.tsv.gz'
gzipped_write(marg, marg_file)



array = adata.X.toarray()

clevel = 4
hf = h5py.File(h5file, 'w')
hf.create_dataset('matrix', data=array,
                  compression='gzip', compression_opts=clevel)
encgenes = [n.encode() for n in adata.var_names]
hf.create_dataset('genes', data=encgenes,
                  compression='gzip', compression_opts=clevel)
encbc = [n.encode() for n in adata.obs_names]
hf.create_dataset('barcodes', data=encbc,
                  compression='gzip', compression_opts=clevel)
hf.close()
del(array)
gc.collect()


