#!/usr/bin/R
# ---------------------------------------------------------
# Replot the DEG enrichment table from the modules package:
# Updated 11/27/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
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


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))

# Set a single color scale:
load.colors()
col_fun = colorRamp2(c(-2, 0, 2), c(colrb[90], "white", colrb[10]))


# Turn into a pair of matrices for plotting:
# ------------------------------------------
dfToHeatmapMatrices = function(df, ckey, pkey){
    cmat = pivot.tomatrix(df[,c('module','jointkey',ckey)], 'jointkey', ckey)
    pmat = pivot.tomatrix(df[,c('module','jointkey',pkey)], 'jointkey', pkey)
    cmat[is.na(cmat)] = 1
    pmat[is.na(pmat)] = 1
    cmat[cmat == 0] = 1e-1
    rownames(cmat) = mmap$mname[as.numeric(rownames(cmat)) + 1]
    cmat = cmat[, rev(colnames(cmat))]
    pmat = pmat[, colnames(cmat)]
    csplit = sub(".*\\.", "", colnames(cmat))
    return(list('cmat'=cmat, 'pmat'=pmat, 'column_split'=csplit))
}


ll = dfToHeatmapMatrices(statsdf, 'ratio','p.adj')
plt = Heatmap(log2(ll$cmat),
              col=col_fun,
              use_raster=TRUE,
              column_split=ll$column_split,
              cluster_columns=FALSE,
              cluster_rows=FALSE,
              width = ncol(ll$cmat)*unit(5, "mm"), 
              height = nrow(ll$cmat)*unit(4, "mm"),
              border_gp = gpar(col="black", lty = 1),
              cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                  p = ll$pmat[i,j]
                  ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                  grid.text(ann, x, y)}
)


h = 2.25 + 2.5 / 15 * nrow(ll$cmat)
w = 5 + 2.5 / 15 * ncol(ll$cmat)
pltprefix = paste0(imgpref, 'module_enrheatmap_', fullpref)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()


# Subset to and plot only modules with significant associations:
# --------------------------------------------------------------
topnames = unique(statsdf$mname[(statsdf$p.adj < 0.05) & (statsdf$key != '--')])
subdf = statsdf[statsdf$mname %in% topnames,]
ll = dfToHeatmapMatrices(subdf, 'ratio','p.adj')

plt = Heatmap(log2(ll$cmat),
              col=col_fun,
              use_raster=TRUE,
              column_split=ll$column_split,
              cluster_columns=FALSE,
              # cluster_rows=FALSE,
              width = ncol(ll$cmat)*unit(5, "mm"), 
              height = nrow(ll$cmat)*unit(4, "mm"),
              border_gp = gpar(col="black", lty = 1),
              cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                  p = ll$pmat[i,j]
                  ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                  grid.text(ann, x, y)}
)


h = 2.25 + 2.5 / 15 * nrow(ll$cmat)
w = 5 + 2.5 / 15 * ncol(ll$cmat)
pltprefix = paste0(imgpref, 'module_enrheatmap_', fullpref, '_sigmodules')
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()

# }


