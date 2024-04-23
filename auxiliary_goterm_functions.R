#!/usr/bin/R
# -----------------------------------------------------
# Auxiliary functions for GO term pruning and plotting:
# Updated: 03/22/2022
# -----------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)

library(ComplexHeatmap)
library(circlize)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Functions for enrichments later:
# --------------------------------
sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
sub.sources = c('REAC','WP','KEGG','CORUM')
runSingleGO = function(set, print.res=TRUE){
    gp2.result = gprofiler2::gost(set, organism='hsapiens',
        ordered_query=FALSE, multi_query=FALSE,
        sources = sources, evcodes=TRUE)
    gp2df = gp2.result$result
    gp2df = gp2df[order(gp2df$p_value),]
    if (print.res) {
        print(head(gp2df[gp2df$term_size < 500,], 25))
    }
    return(gp2df)
}


# Function to reduce enrichments by intersections:
jacc.vec = function(x, y){
    int = length(intersect(x,y))
    un = length(x) + length(y) - int
    return(int / un)
}

selGOterms = function(intlist, cutoff=0.75){
    # Separate to lists:
    intlist = lapply(intlist, function(x){ strsplit(x,",")[[1]] })
    NT = length(intlist)
    jmat = matrix(0, nrow=NT, ncol=NT)
    for (i in 1:NT){
        jmat[i,] = sapply(1:NT, function(j){ jacc.vec(intlist[[i]], intlist[[j]]) })
    }
    # Select above cutoff, flag for removal:
    jmat = jmat - diag(diag(jmat))
    jdf = which(jmat >= cutoff, arr.ind=TRUE)
    jdf = jdf[jdf[,1] < jdf[,2],]
    # Greedy approach to remove:
    rmvec = c()
    for (k in 1:nrow(jdf)){
        if (!(jdf[k, 1] %in% rmvec)){
            rmvec = c(rmvec, jdf[k, 2])
        }
    }
    keep.ind = 1:NT
    keep.ind = keep.ind[!(keep.ind %in% rmvec)]
    return(keep.ind)
}


subsetGOterms = function(gp2df, cutoff=0.75, termsize=500, subsetsources=FALSE){
    subdf = gp2df[gp2df$term_size < termsize,]
    if (subsetsources){
        sub.sources = c('REAC','WP','KEGG','CORUM')
        subdf = subdf[subdf$source %in% sub.sources,]
    }
    subdf = subdf[order(subdf$p_value),]
    intlist = subdf$intersection
    keep.ind = selGOterms(intlist, cutoff=cutoff)
    return(subdf[keep.ind,])
}

reduceGeneList = function(x, ntop=8){
    x = strsplit(x, ',')[[1]]
    x = head(x, ntop)
    x = paste(x, collapse=',')
    return(x)
}

plotGObarplot = function(subdf, pltprefix, fill='grey50', ntop=20){
    subdf = unique(subdf[, c('p_value', 'term_name', 'source', 'intersection'), drop=F])
    subdf$lab = paste0(subdf$source, ": ", subdf$term_name)
    subdf = head(subdf, ntop)
    subdf$lab = factor(subdf$lab, levels=rev(subdf$lab))
    # Reduce the intersection term length
    subdf$int.short = sapply(subdf$intersection, reduceGeneList)
    gp = ggplot(subdf, aes(lab, -log10(p_value), label=int.short)) +
        geom_bar(stat='identity', fill=fill) +
        coord_flip() + geom_text(y=0.1, hjust=0) +
        labs(x='Enriched Term') +
        scale_y_continuous(expand=c(0,0)) +
        theme_pubr()
    h = .5 + nrow(subdf) / 4
    w = 10
    saveGGplot(gp, pltprefix, w=w, h=h)
    return(gp)
}



# Functions to process and plot enrichments:
# ------------------------------------------
# TODO: Improve picking top.ids so that we don't get duplicates:
top.ids = function(x, ntop=10, cutoff=0.05){ 
    ord = order(x, decreasing=F)
    ids = head(ord, ntop)
    ids = ids[x[ids] < cutoff] 
    return(ids) }

gpPvalMatrix = function(gpdf, genesets, ntop=8){
    # Get the enrichment p-values matrix:
    pmat = t(as.matrix(data.frame(gpdf$p_value)))
    rownames(pmat) = NULL
    colnames(pmat) = names(genesets)
    cs = colSums(pmat < 0.05)
    pmat = pmat[, cs > 0]

    # Subset to the top n terms for each:
    ids = apply(pmat, 2, ntop=ntop, top.ids)
    ids = unique(c(unlist(ids)))
    subpmat = -log10(pmat[ids,])
    subpmat[subpmat > 10] = 10
    rownames(subpmat) = gpdf[ids,'term_name']
    return(subpmat)
}

pruneWithInt = function(gpdf, intdf, cutoff=0.75){
    if (!is.null(intdf)){
        intdf = intdf[intdf$term_id %in% gpdf$term_id,]
        intdf = intdf[order(intdf$p_value),]
        # Use this to prune the enriched terms:
        sel.ind = selGOterms(intdf$intersection, cutoff=cutoff)
        int.terms = intdf$term_id[sel.ind]
        # Add terms not in intdf:
        all.gpterms = gpdf$term_id
        uq.gpterms = all.gpterms[!(all.gpterms %in% intdf$term_id)]
        kept.terms = c(int.terms, uq.gpterms)
        print(length(kept.terms) / length(all.gpterms))
        return(gpdf[gpdf$term_id %in% kept.terms,])
    } else {
        print("intdf is NULL")
        return(gpdf)
    }
}

diag.mat3 = function(mat, ratio = 0.5, cutoff = 0.25, flip=FALSE){
    broad.ord = colSums(mat > cutoff) > ratio * nrow(mat)
    diag.ord = apply(mat, 2, which.max)
    if (flip){
        max.ord = apply(mat, 2, max)
        ord = order(broad.ord, diag.ord, -max.ord, decreasing = TRUE)
    } else {
        ord = order(broad.ord, diag.ord, max.ord, decreasing = TRUE)
    }
    mat = mat[, ord]
    cto = apply(mat, 2, which.max)
    idx = colSums(mat > cutoff) > ratio * nrow(mat)
    cto[idx] = 0
    return(list(mat, colnames(mat), cto))
}

plotGpPvalMatrix = function(subpmat, pltprefix, cluster_columns=FALSE, use_raster=TRUE){
    # pltmat = t(reord(t(subpmat)))
    pltmat = t(diag.mat3(t(subpmat), flip=TRUE)[[1]])
    pltmat = pltmat[rev(rownames(pltmat)),]
    maxrow = colnames(pltmat)[apply(pltmat, 1, which.max)]
    rowsplit = rep("Down", nrow(pltmat))
    rowsplit[grep("_up$", maxrow)]= 'Up'
    csplit = ifelse(1:ncol(pltmat) %in% grep("_up$", colnames(pltmat)), 'Up','Down')
    ux = 1.5
    plt = Heatmap(pltmat,
                  cluster_columns=cluster_columns,
                  cluster_rows=FALSE,
                  width = ncol(pltmat)*unit(ux, "mm"), 
                  height = nrow(pltmat)*unit(ux, "mm"),
                  row_split=rowsplit,
                  column_split=csplit,
                  use_raster=use_raster,
                  border_gp=gpar(color='black', lwd=.5),
                  name='-log10(p)',
                  col=c('white',colb))

    h = 2.25 + 1 / 15 * nrow(pltmat)
    w = 5 + 1 / 15 * ncol(pltmat)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

orderPvalMatrix = function(pmat, genesets){
    fixord = c(paste0(genesets, '_up'), paste0(genesets, '_down'))
    cn = colnames(pmat)
    cn = fixord[fixord %in% cn]
    pmat = pmat[,cn]
    if ('allregions' %in% genesets){
        colnames(pmat)[cn == 'allregions_up'] = 'All_up'
        colnames(pmat)[cn == 'allregions_down'] = 'All_down'
    }
    return(pmat)
}

plotGpPvalMatrixReduced = function(subpmat, pltprefix, cluster_columns=FALSE, use_raster=TRUE, ux=1.5){
    # Separate up/down functions:
    maxrow = colnames(subpmat)[apply(subpmat, 1, which.max)]
    rowsplit = rep("Down", nrow(subpmat))
    rowsplit[grep("_up$", maxrow)]= 'Up'
    csplit = ifelse(1:ncol(subpmat) %in% grep("_up$", colnames(subpmat)), 'Up','Down')
    upmat = subpmat[rowsplit == 'Up', csplit == 'Up']
    dwmat = subpmat[rowsplit == 'Down', csplit == 'Down']
    colnames(upmat) = sub("_.*", "", colnames(upmat))
    colnames(dwmat) = sub("_.*", "", colnames(dwmat))
    cn = union(colnames(dwmat), colnames(upmat))

    # Make joint matrix:
    pltmat = matrix(0, nrow=nrow(dwmat) + nrow(upmat), ncol=length(cn),
        dimnames=list(c(rownames(upmat), rownames(dwmat)), cn))
    rowsplit = c(rep("Up", nrow(upmat)), rep("Down", nrow(dwmat)))
    names(rowsplit) = rownames(pltmat)
    pltmat[rownames(upmat), colnames(upmat)] = upmat
    pltmat[rownames(dwmat), colnames(dwmat)] = dwmat

    # Reorder:
    pltmat = t(diag.mat3(t(pltmat), flip=TRUE, ratio=.8)[[1]])
    pltmat = pltmat[rev(rownames(pltmat)),]
    rowsplit = rowsplit[rownames(pltmat)]
    plt = Heatmap(pltmat,
                  cluster_columns=cluster_columns,
                  cluster_row_slices=FALSE,
                  cluster_rows=FALSE,
                  width = ncol(pltmat)*unit(ux, "mm"), 
                  height = nrow(pltmat)*unit(ux, "mm"),
                  row_split=rowsplit,
                  use_raster=use_raster,
                  border_gp=gpar(color='black', lwd=.5),
                  name='-log10(p)',
                  col=c('white',colb))

    h = 1 + 1 / 15 * nrow(pltmat)
    w = 4 + 1 / 15 * ncol(pltmat)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

