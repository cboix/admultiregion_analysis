#!/usr/bin/python
# --------------------------------------------
# Load multiregion and write each as mtx files
# - write raw chuks, merge after
# Updated 09/08/2023
# --------------------------------------------
import re
import gc
import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
import gzip
import time

from adata_handling import adata_loader


import socket
domain = socket.getfqdn()
if 'broadinstitute.org' in domain or 'private' in domain:
    topdir = '/home/cboix/data/'
else:
    topdir = '/home/cboix/'

dbdir = topdir + 'DEVTRAJ/db/'
sdbdir = dbdir + 'multiRegion/'
ucscdir = sdbdir + 'ucsc_datasets/'
outdir = ucscdir + 'ad-multi-region/'


# Load in the preprocessed dataset (w/ UMAP, etc)
# -----------------------------------------------
scdata = adata_loader(celltype='All',
                      subtype=None,
                      dbdir=dbdir,
                      normalize=True,
                      usecombat=False, usehvg=False,
                      remove_mt=False, filter_TH=False)
scdata.load()

adata_dims = scdata.adata.shape
barcodes = scdata.adata.obs_names.tolist()
features = scdata.adata.var_names.tolist()

# Convert:
matrix = scdata.adata.X.T.tocoo()
datamat = np.c_[matrix.row+1, matrix.col+1, matrix.data]
del(scdata, matrix)
gc.collect()


# Alternatively, perform preformatted_write:
# ------------------------------------------
# Faster than for loop or np.savetxt or pd.csv:
def preformatted_write(mat, f, fmtstring=None, encode=False):
    t0 = time.time()
    if fmtstring is None:
        if len(mat.shape) == 1:
            fmtstring = '%g'
        else:
            fmtstring = '\t'.join(['%g']*mat.shape[1])
    fmt = '\n'.join([fmtstring]*mat.shape[0])
    fmt = fmt + "\n"
    print('Format string - ' + str(round(time.time() - t0,2)))
    data = fmt % tuple(mat.ravel())
    print('Format applied - ' + str(round(time.time() - t0,2)))
    if encode:
        data = data.encode('utf-8')
        print('Format encoded - ' + str(round(time.time() - t0,2)))
    f.write(data)
    print('Data written - ' + str(round(time.time() - t0,2)))

def gzipped_write(matrix, filename, fmtstring=None, gzipped=True):
    if gzipped:
        with gzip.open(filename, 'wb') as f:
            preformatted_write(matrix, f, fmtstring=fmtstring, encode=True)
    else:
        with open(filename, 'w') as f:
            preformatted_write(matrix, f, fmtstring=fmtstring, encode=False)


# Reasonably fast save for features/barcodes:
np.savetxt(outdir + 'features.tsv.gz', features, fmt='%s')
np.savetxt(outdir + 'barcodes.tsv.gz', barcodes, fmt='%s')

def write_gzipped_chunk(chunk, filename, gzipped=True):
    print('[STATUS] Starting write for: ' + filename)
    gzipped_write(chunk, filename, "%d %d %d", gzipped=gzipped)
    print('[STATUS] Finished write for: ' + filename)
    return(filename)

# Chunk the tuples into ~80 chunks to feed to parallel proc.
chunksize = int(5e7)
bks = list(range(0, len(datamat), chunksize))

# Formatting isn't the problem in gzip, writing takes a while.
# When not gzip, writing is immediate, formatting is the whole time
chunksize = int(5e7)
bks = list(range(0, len(datamat), chunksize))
for i, bk in tqdm(enumerate(bks)):
    # Better to write temp files to mtx raw than mtx.gz, then concat + gzip
    filename = outdir + "tmp_matrix_%04d.mtx" % i
    chunk = datamat[bk:(bk + chunksize)]
    result = write_gzipped_chunk(chunk, filename, gzipped=False)

# 1e7: gzipped ~ 50s, non ~ 7s -> 1.34hr
# 5e7: non ~ 39s -> 1.4hr


# Write the header file:
# ----------------------
nnz = len(datamat)
headerfile = outdir + "tmp_matrix_header.mtx"
with open(headerfile, 'w') as f:
    h1 = "%%%%MatrixMarket matrix coordinate %s general\n" % "integer"
    h2 = ("%s %s %s\n" % (adata_dims[1], adata_dims[0], nnz))
    print(h1)
    print(h2)
    f.write(h1)
    f.write(h2)

