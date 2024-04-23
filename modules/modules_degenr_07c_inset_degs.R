#!/usr/bin/R
# -----------------------------------------------------------------------
# Plot module insets, coloring/highlighting DEGs for a specific condition
# Updated 12/08/2022
# -----------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(ggrastr)

library(ComplexHeatmap)
library(circlize)
options(width=170)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/pltres/')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}

imgpref = paste0(plotdir, 'modules_indpt_res_', runset, '_')


# Load in and process data and UMAP coordinates:
# ----------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


commandArgs <- function(trailingOnly=TRUE){c(runset)}
source(paste0(sbindir, 'modules/load_modules_umap_coords.R'))


# Load in set of top by p-value functional enrichments for each module:
# ---------------------------------------------------------------------
set = 'coregenes'
toppvals.file = paste0(moddir, 'module_enrichments_toppvals_',
                       set, '_', fullpref, '_small.tsv')
pvalsdf = read.delim(toppvals.file, sep="\t")
pvalsdf$nc = nchar(pvalsdf$term) # For filtering long terms

# Select modules with at least 10 genes:
ctdf = aggregate(gene ~ leiden, nodedf, length)
plt.modules = ctdf$leiden[ctdf$gene >= 10]


# Function to translate edges:
edges.todf = function(edges, lmat){
    e0 = lmat[edges[,1] + 1,]
    e1 = lmat[edges[,2] + 1,]
    data.frame(x0=e0[, 1], y0=e0[, 2],
        x1=e1[, 1], y1=e1[, 2])
}

expand_range = function(rn, scale=.1){
    r = diff(rn) * scale
    rn = c(rn[1] - r, rn[2] + r)
    return(rn)
}

# Plot module network inset:
# --------------------------
# path = 'nft'
for (path in c('nft','plaq_n','plaq_d')){
    print(path)
    sdedf = dedf[(dedf$dkey == path),]
    plot.all = TRUE
    for (i in plt.modules){
        cat(i,':\t')
        if (plot.all){
            # Nodes in module:
            ind = which(nodedf$leiden == i)
            # Connected nodes:
            edges = edgedf[(edgedf$V1 %in% (ind-1)) | (edgedf$V2 %in% (ind-1)),]
            allind = unique(c(unlist(edges))) + 1
        } else {
            # Keep only module nodes and connected nodes with DE results:
            genes = sdedf$gene[sdedf$module == i]
            ind = which(nodedf$gene %in% genes)
            edges = edgedf[(edgedf$V1 %in% (ind-1)) | (edgedf$V2 %in% (ind-1)),]
            allind = unique(c(unlist(edges))) + 1
            allgenes = sdedf$gene[sdedf$gene %in% nodedf[allind,'gene']]
            allind = which(nodedf$gene %in% allgenes)
            edges = edgedf[(edgedf$V1 %in% (allind-1)) & (edgedf$V2 %in% (allind-1)),]
        }
        cat(length(ind), '&', length(allind), '\t')

        # Plot all others as grey
        # Plot only an inset 
        subdf = nodedf[allind,]
        subdf[subdf$leiden != i, 'col'] = 'grey85'
        subdf = merge(subdf, sdedf[,c('gene','logFC_nb','col_nm')], all.x=TRUE)
        subdf$logFC_nb[is.na(subdf$logFC_nb)] = 0
        subdf$col_nm[is.na(subdf$col_nm)] = 0
        subdf$logFC_nb[subdf$leiden != i] = NA
        mx = 0.01
        subdf$logFC_nb[subdf$logFC_nb > mx] = mx
        subdf$logFC_nb[subdf$logFC_nb < -mx] = -mx

        xlim = range(nodedf$L1[ind])
        ylim = range(nodedf$L2[ind])
        r = max(c(diff(xlim), diff(ylim))) / 2
        xlim = c(mean(xlim) - r, mean(xlim) + r)
        ylim = c(mean(ylim) - r, mean(ylim) + r)
        xlim = expand_range(xlim, .2)
        ylim = expand_range(ylim, .2)

        colmap = unique(subdf$col)
        names(colmap) = colmap
        edf = edges.todf(edges, nodedf[,c('L1','L2')])
        limit <- max(abs(subdf$logFC_nb), na.rm=T) * c(-1, 1)
        subdf = subdf[order(abs(subdf$logFC_nb)),]
        subdf = subdf[order(abs(subdf$col_nm)),]
        labdf = subdf[(subdf$col_nm != 0),]
        labdf = labdf[labdf$leiden == i,]
        labdf = labdf[order(abs(labdf$logFC_nb), decreasing=T),]
        labdf = head(labdf, 25)

        gp = ggplot(subdf, aes(L1, L2, col=logFC_nb)) + 
            rasterise(geom_segment(data=edf, aes(x=x0, xend=x1, y=y0, yend=y1), col='grey85', lwd=.25), dpi=450) + 
            geom_point() + 
            geom_text_repel(data=labdf, aes(L1, L2, label=gene), min.segment.length=0.1,
                box.padding=0.01, max.overlaps=30, force=50, force_pull=0.5, segment.color='grey50') +
            scale_color_distiller(palette='RdBu', direction=-1, limit=limit, na.value="grey85") +
            coord_fixed(1, xlim=xlim, ylim=ylim) + 
            theme_pubr() + 
            theme(legend.position=c(0.1,0.1))

        pltprefix = paste0(imgpref, 'graph_degs_', graphpref, '_M', i, '_inset_', path)
        if (plot.all){pltprefix = paste0(pltprefix,'_full')}
        saveGGplot(gp, paste0(pltprefix, '_gp'), h=4, w=4)
    }
}


# Summary number of changes:
# --------------------------
aggdf = aggregate(gene ~ gset + dkey + module, dedf, length)

pcols = brewer.pal(12, 'Paired')
devals = c('Up'=pcols[6], 'Down'=pcols[2], '--' = 'grey85')
gp = ggplot(aggdf, aes(dkey, gene, label=gene, fill=gset)) + 
    facet_wrap(~module) + 
    geom_bar(stat='identity') + 
    geom_text() + 
    scale_fill_manual(values=devals) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() + coord_flip()

pltprefix = paste0(imgpref, 'ndegs_barplot', graphpref)
saveGGplot(gp, paste0(pltprefix, '_gp'), h=8, w=8)

