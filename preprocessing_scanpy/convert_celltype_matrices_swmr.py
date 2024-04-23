# !/usr/bin/python
# -----------------------------------------
# Extract the celltype read count matrices:
# -----------------------------------------
import glob
import h5py
import re
import gzip
import pickle
import numpy as np
import pandas as pd
import time
import gc
import sys

# If these matrices already exist, make them the SWMR mode for concurrent reads

# Find all celltype matrices:
mdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy'
mpref = mdir + prefix + ".majorcelltype."
fnames = np.sort(glob.glob(mpref + "*" + '.hdf5'))
fnames = list(filter(lambda x:'swmr' not in x, fnames))

clevel = 4

# fn = fnames[0]
for fn in fnames:
    print(fn)
    celltype = re.sub('.hdf5', '', re.sub(mpref,"",fn))
    cellstr = re.sub("/","_",celltype)
    print(celltype)
    # Read in the hdf5 file:
    hf = h5py.File(fn, 'r')
    array = hf['matrix']
    ss = array.shape
    arr = array[:,:]
    print(ss)
    bcs = hf['barcodes'][:]
    gns = hf['genes'][:]
    hf.close()
    # Make hdf5 file in SWMR mode:
    swfile = re.sub('.hdf5', '.swmr.hdf5', fn)
    hf = h5py.File(swfile, 'w', libver='latest')
    hd = hf.create_dataset('matrix', data=arr,
                           compression='gzip', compression_opts=clevel)
    hf.create_dataset('genes', data=gns,
                      compression='gzip', compression_opts=clevel)
    hf.create_dataset('barcodes', data=bcs,
                      compression='gzip', compression_opts=clevel)
    hf.swmr_mode = True # NOTE: Has to come after datasets.
    hf.close()
    del(arr)
    gc.collect()


# # TEST MAKE SWMR:
# f = h5py.File("swmr.h5", 'w', libver='latest')
# arr = np.array([1,2,3,4])
# dset = f.create_dataset("data", chunks=(2,), maxshape=(None,), data=arr)
# f.swmr_mode = True

# TEST READ:
# # Now it is safe for the reader to open the swmr.h5 file
# for i in range(5):
#     new_shape = ((i+1) * len(arr), )
#     dset.resize( new_shape )
#     dset[i*len(arr):] = arr
#     dset.flush()
#     # Notify the reader process that new data has been written

