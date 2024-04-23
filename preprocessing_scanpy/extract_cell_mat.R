#!/usr/bin/R
# ----------------------
# Extract cell matrices:
# Updated: 01/20/2021
# ----------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(rhdf5)
library(Matrix)

# Function to estimate end time of a loop:
est.endtime <- function(start, last, done, total){
    curr = proc.time() # Get current time
    # Calculate how much time has passed:
    last.step = (proc.time() - last)[3]
    elapsed = (proc.time() - start)[3]
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

# --------------------------------
# Load in the barcodes, pathology:
# --------------------------------
hcts = unique(cellmeta$major.celltype)
for (hct in hcts){
    print(hct)
    celltype = sub("/","_", hct)

    rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
    if (dbdir == '~/data/DEVTRAJ/db/') {
        matdir = paste0(datadir, 'matrices/')
        # matdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
    } else {
        matdir = paste0(datadir, 'matrices/')
    }
    mtxdir = paste0(matdir, 'mtx/')
    system(paste('mkdir -p', mtxdir))
    h5file = paste0(matdir, rawpref, '.majorcelltype.', celltype, '.hdf5')

    # Extract metadata from hdf5:
    h5f = H5Fopen(h5file)
    genes = h5f$genes
    barcodes = h5f$barcodes
    H5Fclose(h5f)
    ngenes = length(genes)

    # Make metadata table:
    rownames(cellmeta) = cellmeta$barcode
    pathdf = cellmeta[barcodes,]

    # Split by ct:
    split.var = 'major.celltype'
    # if (celltype %in% c('Exc','Inh', 'Vasc_Epithelia')){ 
    if (celltype %in% c('Vasc_Epithelia')){ 
        split.var = 'cell_type_high_resolution' 
    } else if (celltype == 'Mic_Immune'){
        split.var = 'minor.celltype'
        if (celltype == 'Mic_Immune'){
            pathdf$minor.celltype = 'Mic'
            pathdf$minor.celltype[pathdf$cell_type_high_resolution == 'T cells'] = 'T'
            pathdf$minor.celltype[pathdf$cell_type_high_resolution == 'CAMs'] = 'CAM'
        }
    }
    subtypes = unique(pathdf[[split.var]])

    # Gather + write each region + subtype:
    for (subtype in subtypes){
        for (region in regions[regions != 'MB']){
            ststr = gsub("/","_",gsub(" ","_", subtype))
            matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                             celltype,'.',ststr,'.',region)
            mtxfile = paste0(matpref, '.mtx.gz') 
            rnfile = paste0(matpref, '.rn.tsv.gz') 
            cnfile = paste0(matpref, '.cn.tsv.gz') 
            rdafile = paste0(matpref, '.rda')  # In MM format
            sub.pathdf = pathdf[pathdf[[split.var]] == subtype & pathdf$region == region,] 
            if (!file.exists(rdafile)){
                if (nrow(sub.pathdf) > 1){
                    print(paste("[STATUS] Processing", subtype, 'in',region,'with',nrow(sub.pathdf), 'cells'))
                    bind = match(sub.pathdf$barcode, barcodes)
                    subbcs = barcodes[bind] 
                    t0 = proc.time()
                    t1 = t0
                    # Open handle, extract cells we care about and close:
                    h5f = H5Fopen(h5file)
                    genes = h5f$genes
                    bcs = h5f$barcodes
                    h5d = h5f&"matrix"
                    mat = t(h5d[,bind])
                    H5Dclose(h5d)
                    H5Fclose(h5f)
                    colnames(mat) = genes
                    rownames(mat) = bcs[bind]
                    # Order as model matrix + transpose matrix:
                    mat = mat[sub.pathdf$barcode,]
                    mat = t(Matrix(mat))
                    gcout = gc()
                    t1 = est.endtime(t0, t1, 1, 3)
                    # Save matrices:
                    print("[STATUS] Writing Rda")
                    save(mat, file=rdafile)
                    t1 = est.endtime(t0, t1, 2, 3)
                    print("[STATUS] Writing mtx file")
                    writeMM(mat, file=gzfile(mtxfile)) # TODO: Fix line, not writing
                    # Write barcodes:
                    write.table(rownames(mat), rnfile, quote=F, sep="\t", row.names=F, col.names=F)
                    write.table(colnames(mat), cnfile, quote=F, sep="\t", row.names=F, col.names=F)
                    t1 = est.endtime(t0, t1, 3, 3)
                    rm(mat)
                    gcout=gc()
                }
            }
        }
    }
}

