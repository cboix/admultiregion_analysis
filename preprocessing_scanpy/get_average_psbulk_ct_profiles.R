#!/usr/bin/R
# ----------------------------------------------
# Make average pseudobulk profiles (for website)
# Updated: 04/29/2022
# ----------------------------------------------
# chunksize=25
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(rhdf5)
library(Matrix)



hcts = unique(cellmeta$major.celltype)
psprefix = 'subtype_reg/pseudobulk_data_'

avg.mat = NULL
ncell = NULL
for (subct in hcts){
    ctstr = gsub("/", "_", subct)
    load(paste0(sdbdir, psprefix, ctstr, '.rda'))
    ps.data$meta$pr = with(ps.data$meta, paste0(region, '@', cell_type_high_resolution))
    prs = unique(ps.data$meta$pr)
    tform = make.tform(ps.data$meta$pr, u=prs, norm=FALSE)
    tform = sweep(tform, 1, ps.data$meta$ncell, '*')
    sub.ncell = apply(tform, 2, sum)
    tform = sweep(tform, 2, sub.ncell, '/')
    submat = ps.data$mat %*% tform 
    avg.mat = cbind(avg.mat, submat)
    ncell = c(ncell, sub.ncell)
}

rdsfile = paste0(sdbdir, psprefix, 'allcts_averageprofiles.rds')
saveRDS(avg.mat, file=rdsfile)

ncfile = paste0(sdbdir, psprefix, 'allcts_averageprofiles_ncells.rds')
saveRDS(ncell, file=ncfile)

# NOTE: Missing NRGN for now, need to decide if clean up or relevant CT.

# Mapping:
cmapdf = unique(cellmeta[,c('cell_type_high_resolution','major.celltype')])
cmap = cmapdf$major.celltype
names(cmap) = cmapdf$cell_type_high_resolution
cmfile = paste0(sdbdir, 'celltype_mapping.rds')
saveRDS(cmap, file=cmfile)

