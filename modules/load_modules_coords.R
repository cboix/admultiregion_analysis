#!/usr/bin/R
# ---------------------------------------------
# Load the coordinate system for modules' UMAP:
# Updated 01/25/2023 
# ---------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){ source(paste0(sbindir, 'load_metadata.R')) }
options(width=170)

# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: runset graph_id")
} else {
    runset = args[1]
    graph_id = args[2]
}


if (runset %in% c('Exc', 'Inh', names(exc.sets))){
    submeta = read.delim(ctumap.file, sep="\t")
    submeta = submeta[, c('barcode','C1','C2')]
    names(submeta) = c('barcode','U1','U2')
    rownames(cellmeta) = cellmeta$barcode
    submeta$region = cellmeta[submeta$barcode, 'region']
    submeta$projid = cellmeta[submeta$barcode, 'projid']
    submeta$cell_type_high_resolution = cellmeta[submeta$barcode, 'cell_type_high_resolution']
    rownames(submeta) = submeta$barcode
    submeta = submeta[rownames(scoremat),]
} else if (runset == 'Vasc_Epithelia'){
    submeta = read.delim(ctumap.file)
    rownames(submeta) = submeta$barcode
    submeta = submeta[rownames(scoremat),c('barcode','region','celltype','C1','C2')]
    names(submeta)[3:5] = c('cell_type_high_resolution', 'U1','U2')
} else {
    rownames(cellmeta) = cellmeta$barcode
    submeta = cellmeta[rownames(scoremat),]
}



# Plotting function for each of the module scores:
# ------------------------------------------------
ind = 1:nrow(submeta)
cex = 0.025

# Table of params for each celltype to center + not change aspect ratio:
center.params = list("Ast"=c(0, 8.5, -8.5, 0),
                     "Mic_Immune"=c(1.25, 9.25, -12.5, -4),
                     "Opc"=c(0, 6.5, -8.5, 0),
                     "Oli"=c(-1,11,-14,-2),
                     "HCneurons"=c(14,32,-18,0))

# Set the x and y limits using these:
if (runset %in% names(center.params)){
    cp = center.params[[runset]]
    xlim = c(min(submeta$U1) + cp[1], min(submeta$U1) + cp[2])
    ylim = c(max(submeta$U2) + cp[3], max(submeta$U2) + cp[4])
} else {
    xlim = range(submeta$U1) 
    ylim = range(submeta$U2)
    r = max(c(diff(xlim), diff(ylim))) / 2
    xlim = c(mean(xlim) - r, mean(xlim) + r)
    ylim = c(mean(ylim) - r, mean(ylim) + r)
}

