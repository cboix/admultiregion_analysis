#!/usr/bin/R
# ----------------------
# Turn the cell type matrices into chunked matrices
# Updated: 04/28/2022
# ----------------------
# chunksize=25
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(rhdf5)
library(Matrix)

# Function to estimate end time of a loop:
est.endtime <- function(start, last, done, total){
    curr = proc.time() # Get current time
    # Calculate how much time has passed:
    last.step = (proc.time() - last)[3]
    elapsed = (proc.time() - start)[3]
    cat(paste0(done,'/', total, '\t'))
    cat(paste0(round(last.step,1),'s\t'))
    cat(paste0(round(elapsed,1),'s\t'))
    # Estimate the remaining time:
    # Estimate the end time:
    each.time = elapsed / done
    est.left = (total - done) * each.time
    cat(paste0('Left: ', round(est.left,1),'s\t'))
    fin.time = Sys.time() + est.left
    cat(paste0('Est. ', format(fin.time, "%H:%M:%S"),'\n'))
    return(curr)
}


rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    matdir = paste0(datadir, 'matrices/')
    # matdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')
rdsdir = paste0(matdir, 'rds/')
system(paste('mkdir -p', mtxdir, rdsdir))



# Load in the expression matrix:
# ------------------------------
chunksize = 25

hcts = unique(cellmeta$major.celltype)
for (hct in hcts){
    print(hct)
    celltype = sub("/","_", hct)
    matpref = paste0(matdir, rawpref, '.majorcelltype.', celltype)
    rdspref = paste0(rdsdir, rawpref, '.majorcelltype.', celltype)
    h5file = paste0(matpref, '.hdf5')

    # Extract metadata from hdf5:
    h5f = H5Fopen(h5file)
    genes = h5f$genes
    barcodes = h5f$barcodes
    H5Fclose(h5f)
    ngenes = length(genes)
    
    # Write barcodes and genes:
    saveRDS(barcodes, file=paste0(rdspref, '_barcodes.Rds'))
    saveRDS(genes, file=paste0(rdspref, '_genes.Rds'))

    # Read in chunks of the data:
    nchunk = ceiling(ngenes / chunksize)
    t0 = proc.time()
    t1 = t0
    for (i in 1:nchunk){
        chunkstr = sprintf('_chunk%05d_size%d', i, chunksize)
        rdsfile = paste0(rdspref, chunkstr, '.rds')

        if (!file.exists(rdsfile)){
            ind = ((i - 1) * chunksize + 1):(min(ngenes, i * chunksize))

            h5f = H5Fopen(h5file)
            genes = h5f$genes
            bcs = h5f$barcodes
            h5d = h5f&"matrix"
            mat = h5d[ind,]
            H5Dclose(h5d)
            H5Fclose(h5f)
            rownames(mat) = genes[ind]

            # Save matrix:
            saveRDS(mat, file=rdsfile)
            t1 = est.endtime(t0, t1, i, nchunk)
        }
    }
}


