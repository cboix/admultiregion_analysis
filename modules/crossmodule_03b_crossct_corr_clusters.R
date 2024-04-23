#!/usr/bin/R
# --------------------------------------------------
# Enrichments + plots for the cross-module clusters:
# Updated 03/27/2022
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(ggpubr)
library(gprofiler2)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# plotGpPvalMatrix = function(subpmat, pltprefix, cluster_columns=FALSE){
#     # pltmat = t(reord(t(subpmat)))
#     pltmat = t(diag.mat2(t(subpmat))[[1]])
#     pltmat = pltmat[rev(rownames(pltmat)),]
#     maxrow = colnames(pltmat)[apply(pltmat, 1, which.max)]
#     ux = 1.5
#     plt = Heatmap(pltmat,
#                   cluster_columns=cluster_columns,
#                   cluster_rows=TRUE,
#                   width = ncol(pltmat)*unit(ux, "mm"), 
#                   height = nrow(pltmat)*unit(ux, "mm"),
#                   use_raster=TRUE,
#                   border_gp=gpar(color='black', lwd=.5),
#                   name='-log10(p)',
#                   col=c('white',colb))
#     h = 2.25 + 1 / 15 * nrow(pltmat)
#     w = 5 + 1 / 15 * ncol(pltmat)
#     saveHeatmap(plt, pltprefix, w=w, h=h)
# }


plotGObarplot = function(subdf, pltprefix, fill='grey50', ntop=20, title=NULL){
    subdf = unique(subdf[, c('p_value', 'term_name', 'source', 'intersection'), drop=F])
    subdf$lab = paste0(subdf$source, ": ", subdf$term_name)
    subdf = head(subdf, ntop)
    subdf$lab = factor(subdf$lab, levels=rev(subdf$lab))
    # Reduce the intersection term length
    # subdf$int.short = sapply(subdf$intersection, reduceGeneList)
    # gp = ggplot(subdf, aes(lab, -log10(p_value), label=int.short)) +
    gp = ggplot(subdf, aes(lab, -log10(p_value))) +
        geom_bar(stat='identity', fill=fill) +
        coord_flip() + 
        # geom_text(y=0.1, hjust=0) +
        labs(x='Enriched Term', title=title) +
        scale_y_continuous(expand=c(0,0)) +
        theme_pubr()
    h = .5 + nrow(subdf) / 4
    w = 10
    saveGGplot(gp, pltprefix, w=w, h=h)
    return(gp)
}


# Read in the module clusters:
# ----------------------------
modcls.tsv = paste0(crossdir, 'shared_genes_module_clusters.tsv')
cls.tsv = paste0(crossdir, 'module_cluster_assignments.tsv')
clsdf = read.delim(modcls.tsv, header=T)
clsassign = read.delim(cls.tsv, header=T)
clsdf$cls = paste0('C', clsdf$cls)

acls = unique(clsdf$cls)



# Enrichments for each cluster:
for (i in acls){
    genes = clsdf$gene[clsdf$cls == i]
    topgenes = clsdf$gene[clsdf$count >= 3 & clsdf$cls == i]
    color = clsdf$col[clsdf$cls == i][1]

    cat(i, '\n', topgenes, '\n\n')
}

    gp2.result = gprofiler2::gost(genes, organism='hsapiens',
        ordered_query=FALSE, multi_query=FALSE,
        sources = sources, evcodes=TRUE)
    gpdf = gp2.result$result
    if (!is.null(gpdf)){
        gpdf$nc = nchar(gpdf$term_name)
        gpdf = gpdf[gpdf$term_size < 500,]
        gpdf = gpdf[gpdf$nc < 60,]
        gpdf = gpdf[order(gpdf$p_value),]
        if (nrow(gpdf) > 10){
            gpdf = pruneWithInt(gpdf, gpdf, cutoff=0.5)
        }
        if (nrow(gpdf) > 1){
            pltprefix = paste0(imgpref, 'modules_cluster_enr_barplot_', i)
            plotGObarplot(gpdf, pltprefix, fill=color, ntop=10, title=i)
        }
    }
}

