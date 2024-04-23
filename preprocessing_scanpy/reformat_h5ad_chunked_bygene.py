# !/usr/bin/python
# -----------------------------------------
# Turn the h5ad files into chunked-by-gene hdf5
# TODO: Both raw and normalized
# NOTE: Works at 80G vmem
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

import scanpy.api as sc
import anndata

# Scanpy plotting settings:
sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving

# ---------------------------------------------
# Load the read count (or normalized) datasets:
# ---------------------------------------------
use_norm = True
mdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy'

if use_norm:
    prefix = prefix + "_norm"
    h5ad_file =  './write/' + prefix + '_fullmatrix.h5ad'
else:
    h5ad_file =  mdir + prefix + '_fullmatrix.h5ad'

outfile = mdir + prefix + '_bygene_fullmatrix.swmr.hdf5'
adata = sc.read_h5ad(h5ad_file)


# ---------------------------------------------------------
# Parameters for densifying and writing the data in chunks:
# ---------------------------------------------------------
cs = 10
ncol = adata.shape[1]
chunkshape = (adata.shape[0], cs)
pull = 1000
npull = int(np.floor(ncol / pull) + 1)

# ---------------------------------------------
# Initialize the dataset, write to it in chunks
# ---------------------------------------------
clevel = 4
hf = h5py.File(outfile, 'w', libver='latest')
hd = hf.create_dataset('matrix', shape=adata.shape,
                       chunks=chunkshape, dtype=adata.X.dtype,
                       compression='gzip', compression_opts=clevel)
# Bigger chunks will work better because it is csr format:
# NOTE: The final chunk shape should still stay small
t0 = time.time()
for i in range(npull):
    t1 = time.time()
    print(i, np.round(t1 - t0,2))
    ind = np.array(range((i * pull), np.min([ncol, (i+1) * pull])))
    subad = adata[:,ind]
    chunk = subad.X.toarray()
    hd[:,ind] = chunk

# Add the genes and barcodes:
encgenes = [n.encode() for n in adata.var_names]
hf.create_dataset('genes', data=encgenes,
                  compression='gzip', compression_opts=clevel)
encbc = [n.encode() for n in adata.obs_names]
hf.create_dataset('barcodes', data=encbc,
                  compression='gzip', compression_opts=clevel)

# Change the mode to SWMR for multiple readers
hf.swmr_mode = True # NOTE: Has to come after dataset creation.
hf.close()



