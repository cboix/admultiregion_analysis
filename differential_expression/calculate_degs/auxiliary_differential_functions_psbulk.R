#!/usr/bin/R
# ----------------------------------------
# Auxiliary functions for DE calculations:
# Multiregion (matching AD430)
# Updated: 10/05/22
# ----------------------------------------
library(cbrbase)
library(tidyr)
library(Matrix)


# Functions:
# ----------
# Ensure data matches datasets on both meta/mat:
harmonizePsbulkData = function(ps.data, col='ptype'){
    rownames(ps.data$meta) = ps.data$meta[[col]]
    ps.data$mat = ps.data$mat[,ps.data$meta[[col]]]
    cat('meta:', dim(ps.data$meta),
        '\tmat:', dim(ps.data$mat), '\n')
    return(ps.data)
}


# Subset the counts matrix:
subsetMatrixForDE = function(mat, pathdf, pctcells, pctcut=NULL, selgenes=NULL){
    if (is.null(selgenes)){
        if (!is.null(pctcut)){
            keep.genes = names(pctcells)[pctcells > pctcut]
        }
        # Remove ribosomal genes
        keep.genes = keep.genes[grep("^RP[0-9]*-",keep.genes, invert=TRUE)] 
        keep.genes = keep.genes[grep("^RP[SL]",keep.genes, invert=TRUE)]
    } else {
        keep.genes = selgenes[selgenes %in% rownames(mat)]
    }
    mat = mat[keep.genes, pathdf$barcode]
    print(paste("[STATUS] Subsetting matrix to", 
            paste0(dim(mat), collapse = ' x '),'(g x c)'))
    return(mat)
}


# Subset the ps.data type dataset:
subsetPsbulkMatrixForDE = function(ps.data, pctcut=NULL, 
    selgenes=NULL, path=NULL, col='ptype'){
    if (is.null(selgenes)){
        if (is.null(pctcut)){
            keep.genes = rownames(ps.data$mat)
        } else {
            pctcells = rowSums(ps.data$mat > 0) / ncol(ps.data$mat)
            keep.genes = names(pctcells)[pctcells > pctcut]
        }
        # Remove ribosomal genes
        keep.genes = keep.genes[grep("^RP[0-9]*-",keep.genes, invert=TRUE)] 
        keep.genes = keep.genes[grep("^RP[SL]",keep.genes, invert=TRUE)]
    } else {
        keep.genes = selgenes[selgenes %in% rownames(ps.data$mat)]
    }
    # Subset dataset to samples without NA values in tested variable:
    if (!is.null(path)){
        # TODO: Throw error if path missing:
        ind = which(!is.na(ps.data$meta[[path]]))
        ps.data$meta = ps.data$meta[ind,]
    }
    ps.data$mat = ps.data$mat[keep.genes, ps.data$meta[[col]]]
    print(paste("[STATUS] Subsetting matrix to", 
            paste0(dim(ps.data$mat), collapse = ' x '),'(g x c)'))
    return(ps.data)
}


# Make volcano plot for DE results:
plotVolcano = function(resdf, prefstr, imgpref, ntop=30){
    require(ggplot2)
    require(ggpubr)
    require(ggrepel)
    source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

    labdf = head(resdf[resdf$col != 'NS',], ntop)
    pcols = brewer.pal(12,'Paired')
    devals = c('NS'='grey80', 'Down'=pcols[2], 'Up'=pcols[6])

    gp = ggplot(resdf, aes(log2FoldChange, -log10(pvalue), color=col)) + 
        geom_vline(xintercept=0, lty='dotted', lwd=.5) + 
        geom_point(cex=.25) + 
        geom_text_repel(data=labdf, aes(log2FoldChange, -log10(pvalue), label=gene, color=col), size=2, max.overlaps=20) + 
        scale_color_manual(values=devals, name='DE status')+ 
        scale_y_continuous(expand=c(0,0)) + 
        theme_pubr() + theme(legend.position='none')

    pltprefix = paste0(imgpref, 'volcano_', prefstr)
    saveGGplot(gp, pltprefix, w=5, h=5)
}


