#!/usr/bin/R
# -----------------------------------------------
# Plot module genes on its top KEGG pathways:
# TODO: pathview > SBGNview (updated w/ reactome)
# Updated 12/21/2022
# -----------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(pathview)
library(tidyr)
library(viridis)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(bindir, 'auxiliary_function_general_repel.R'))

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
for (i in plt.modules){
    cat(i,':\t')
    # Nodes in module:
    # Connected nodes:
    ind = which(nodedf$leiden == i)
    edges = edgedf[(edgedf$V1 %in% (ind-1)) | (edgedf$V2 %in% (ind-1)),]
    allind = unique(c(unlist(edges))) + 1
    cat(length(ind), '&', length(allind), '\t')

    # Plot all others as grey
    # Plot only an inset 
    subdf = nodedf[allind,]
    subdf[subdf$leiden != i, 'col'] = 'grey85'

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

    gp = ggplot(subdf, aes(L1, L2, col=col)) + 
        geom_segment(data=edf, aes(x=x0, xend=x1, y=y0, yend=y1), col='grey85', lwd=.25) + 
        geom_text_repel(data=subdf, aes(L1, L2, label=gene), 
            box.padding=0.01, max.overlaps=10, force=100, force_pull=0.5, segment.color='grey50') +
        geom_point() + 
        scale_color_manual(values=colmap) +
        coord_fixed(1, xlim=xlim, ylim=ylim) + 
        theme_pubr() + 
        theme(legend.position='none')

    pltprefix = paste0(imgpref, 'graph_layout_', graphpref, '_M', i, '_inset')
    saveGGplot(gp, paste0(pltprefix, '_gp'), h=4, w=4)
}

