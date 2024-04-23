#!/usr/bin/R
# -----------------------------------------------
# Functions for processing gprofiler enrichments:
# Updated 12/09/2022
# -----------------------------------------------
library(cbrbase)
library(tidyr)
library(gprofiler2)
library(ComplexHeatmap)
library(circlize)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Functions to process and plot enrichments:
# ------------------------------------------
# TODO: Improve picking top.ids so that we don't get duplicates:
top.ids = function(x, ntop=10, cutoff=0.05){ 
    ord = order(x, decreasing=F)
    ids = head(ord, ntop)
    ids = ids[x[ids] < cutoff] 
    return(ids) }

gpPvalMatrix = function(gpdf, genesets, ntop=8, keep.sets=NULL, thresh=10){
    # Get the enrichment p-values matrix:
    pmat = t(as.matrix(data.frame(gpdf$p_value)))
    rownames(pmat) = NULL
    colnames(pmat) = names(genesets)
    if (!is.null(keep.sets)){
        pmat = pmat[, keep.sets, drop=F]
    }
    cs = colSums(pmat < 0.05)
    pmat = pmat[, cs > 0, drop=F]

    # Subset to the top n terms for each:
    ids = apply(pmat, 2, ntop=ntop, top.ids)
    ids = unique(c(unlist(ids)))
    subpmat = -log10(pmat[ids,,drop=F])
    if (!is.null(thresh)){ subpmat[subpmat > thresh] = thresh }
    rownames(subpmat) = gpdf[ids,'term_name']
    return(subpmat)
}

plotGpPvalMatrix = function(subpmat, pltprefix, ux=1.5){
    # pltmat = t(reord(t(subpmat)))
    pltmat = t(diag.mat2(t(subpmat))[[1]])
    pltmat = pltmat[rev(rownames(pltmat)),]
    rowsplit = colnames(pltmat)[apply(pltmat, 1, which.max)]
    rowsplit = sapply(rowsplit, function(x){ sub("^M","",sub(" .*","",x)) })
    rowsplit = as.numeric(rowsplit)
    gap = 0
    plt = Heatmap(pltmat,
                  cluster_columns=FALSE,
                  cluster_rows=FALSE,
                  cluster_column_slices=FALSE,
                  cluster_row_slices=FALSE,
                  row_split=rowsplit,
                  use_raster=TRUE,
                  border_gp=gpar(color='black'),
                  name='-log10(p)',
                  row_gap=unit(gap, "mm"),
                  column_gap=unit(gap, "mm"),
                  width = ncol(pltmat)*unit(ux, "mm"), 
                  height = nrow(pltmat)*unit(ux, "mm"),
                  col=c('white',colb))

    h = 2.25 + 1 / 15 * nrow(pltmat)
    w = 5 + 1 / 15 * ncol(pltmat)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

