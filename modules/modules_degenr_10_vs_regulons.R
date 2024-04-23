#!/usr/bin/R
# ----------------------------------------------
# Intersect the SCENIC regulons with the modules
# Updated 01/25/2022
# ----------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'regulons_')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)


# Set the run arguments:
# ----------------------
# TODO: Load across all?
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Read in the regulon data:
# ----------------------------
regulonfile = 'Brain_region_regulons_1000.csv'
rwide = read.delim(paste0(sdbdir, regulonfile), header=F, sep=",")

rsets = lapply(1:nrow(rwide), function(i){
    x = as.character(rwide[i, 2:ncol(rwide)])
    x = x[x != ""]
    unique(x)
})

names(rsets) = rwide[,1]
NR = length(rsets)


# Calculate overlaps with each module:
# ------------------------------------
NM = max(coremap) + 1
rtot = rep(0, NR)
mtot = table(coremap)
# jmat = matrix(0, nrow=NR, ncol=NM)
tmat = matrix(0, nrow=NR, ncol=NM)
# TODO: CHANGE to use either...?
use.mapping = coremap

for (i in 1:NR){
    rvec = rsets[[i]]
    rvec = rvec[rvec %in% names(coremap)]
    rtot[i] = length(rvec)
    if (length(rvec) > 0){
        mvec = coremap[rvec]
        mtab = table(mvec)
        mind = as.numeric(names(mtab)) + 1
        tmat[i, mind] = mtab
    }
}

umat = matrix(rep(rtot, NM), byrow=F, nrow=NR) + matrix(rep(mtot, NR), byrow=T, nrow=NR) 
jmat = tmat / (umat - tmat)

for (j in 1:30){
    rind = which.max(jmat[,j])
    cat(j-1, '\t', names(rsets)[rind], '\t', jmat[rind,j], '\t', mmap[j,'mname'], '\n')
    rvec = rsets[[rind]]
    mvec = names(coremap)[coremap == j-1]
    int = rvec[rvec %in% mvec]
    if (length(int) > 10){
        cat(int, '\n\n')
    }
}

