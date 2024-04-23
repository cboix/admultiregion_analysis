#!/usr/bin/R
# -------------------------------------------------------------
# Replot the modules representations for extended data figures:
# Updated 03/18/2022 to standardize sizes
# -------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


library(tidyr)
library(viridis)
library(igraph)

library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'repr_')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


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


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Plot graph from representation data:
# ------------------------------------
add.edges = function(edges, lmat, col='grey75', lwd=.1){
    e0 = lmat[edges[,1] + 1,]
    e1 = lmat[edges[,2] + 1,]
    segments(x0=e0[, 1], y0=e0[, 2],
             x1=e1[, 1], y1=e1[, 2], 
             col=col, lwd=lwd)
}


w = 1.75
cex = 0.3
sp = 0.1
xlim = range(nodedf$L1) * 1.02
ylim = range(nodedf$L2) * 1.02

mingenes = 4
ngdf = aggregate(gene ~ leiden, nodedf, length)
mod.loc = aggregate(cbind(L1, L2) ~ leiden, nodedf, mean)
mod.loc = merge(mod.loc, ngdf)
mod.loc = mod.loc[mod.loc$gene >= mingenes,]


pltprefix = paste0(imgpref, 'graph_layout_', graphpref)
png(paste0(pltprefix, '_notext.png'), units='in', res=450, width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(0, 4))
plot(nodedf$L1, nodedf$L2, xlim=xlim, ylim=ylim, type='n', axes=F)
add.edges(edgedf, nodedf[,c('L1','L2')], col='grey75', lwd=.1)
points(nodedf$L1, nodedf$L2, col=nodedf$col, pch=19, cex=cex / 2)
dev.off()

pdf(paste0(pltprefix, '_textonly.pdf'), width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(0, 4))
plot(nodedf$L1, nodedf$L2, xlim=xlim, ylim=ylim, type='n', axes=F, ylab='', xlab='')
box() # For aligning:
text(mod.loc$L1, mod.loc$L2, paste0('M', mod.loc$leiden), cex=.35, font=2)
dev.off()


# Plot UMAP from representation data:
# -----------------------------------
w = 1.75
cex = 0.1
sp = 0.1
xlim = range(nodedf$U1) * 1.02
ylim = range(nodedf$U2) * 1.02

# Label modules with at least four genes:
mingenes = 4
ngdf = aggregate(gene ~ leiden, nodedf, length)
mod.loc = aggregate(cbind(U1, U2) ~ leiden, nodedf, mean)
mod.loc = merge(mod.loc, ngdf)
mod.loc = mod.loc[mod.loc$gene >= mingenes,]

pltprefix = paste0(imgpref, 'UMAP_layout_', graphpref)
png(paste0(pltprefix, '_notext.png'), units='in', res=450, width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(0, 4))
plot(nodedf$U1, nodedf$U2, xlim=xlim, ylim=ylim, type='n', axes=F)
points(nodedf$U1, nodedf$U2, col=nodedf$col, pch=19, cex=cex / 2)
dev.off()

pdf(paste0(pltprefix, '_textonly.pdf'), width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(0, 4))
plot(nodedf$U1, nodedf$U2, xlim=xlim, ylim=ylim, type='n', axes=F, ylab='', xlab='')
box() # For aligning:
text(mod.loc$U1, mod.loc$U2, paste0('M', mod.loc$leiden), cex=.35, font=2)
dev.off()




# Print list of genes (extraneous to plotting representations).
# -------------------------------------------------------------
for (i in 0:max(coremap)){
    cat(i, '\t')
    genes = names(coremap)[coremap == i]
    genes = sort(genes)
    cat(genes)
    cat('\n')
}


