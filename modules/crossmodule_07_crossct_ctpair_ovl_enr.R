#!/usr/bin/R
# --------------------------------------------------
# Calculate the gene overlap and overlap enrichment:
# for pairs of modules across two cell types.
# NOTE: Mainly for OPC / Oli
# Updated 03/26/2022
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(gprofiler2)
library(eulerr)

library(viridis)
library(ggpubr)
library(ggrastr)
library(circlize)
library(ComplexHeatmap)


# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
# TODO: Opt to include neurons or not? Or just add later...
# commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


plotGpPvalMatrix = function(subpmat, pltprefix, cluster_columns=FALSE){
    # pltmat = t(reord(t(subpmat)))
    pltmat = t(diag.mat2(t(subpmat))[[1]])
    pltmat = pltmat[rev(rownames(pltmat)),]
    maxrow = colnames(pltmat)[apply(pltmat, 1, which.max)]
    ux = 1.5
    plt = Heatmap(pltmat,
                  cluster_columns=cluster_columns,
                  cluster_rows=TRUE,
                  width = ncol(pltmat)*unit(ux, "mm"), 
                  height = nrow(pltmat)*unit(ux, "mm"),
                  use_raster=TRUE,
                  border_gp=gpar(color='black', lwd=.5),
                  name='-log10(p)',
                  col=c('white',colb))
    h = 2.25 + 1 / 15 * nrow(pltmat)
    w = 5 + 1 / 15 * ncol(pltmat)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

# Get the gene lists for each runset:
# -----------------------------------
# mingenes = 10
runsets = c('Opc','Oli')
rstr = paste(runsets, collapse='_')
imgpref = paste0(plotdir, 'cross2ct_', rstr)

# Get sets function, default to expanded genes:
getGeneSet <- function(module, uselist=gmlist){
    rs = sub("-.*", "", module)
    num = as.numeric(sub(".*-", "", module))
    coremap = uselist[[rs]]
    x = names(coremap)[coremap == num]
    return(x)
}

# Selected pairs for OPC/Oligodendrocytes:
pltpairs = list(c('Opc-25','Oli-23'), c('Opc-4','Oli-4'), c('Opc-4', 'Oli-1'),
                c('Opc-1','Oli-3'), c('Opc-5','Oli-1'), c('Opc-12','Oli-1'))

for (i in 1:length(pltpairs)){
    # Make the genesets from this pair:
    # ---------------------------------
    pair = pltpairs[[i]]
    pstr = paste(pair, collapse='_')
    genesets = lapply(pair, getGeneSet)
    names(genesets) = pair

    # Plot set sizes:
    pltname = paste0(imgpref, '_eulerovl_', pstr, '.png')
    png(pltname, units='in', res=300, width=2.25, height=1.5)
    par(yaxs='i', xaxs='i', mar=rep(0,4))
    plot(euler(genesets), quantities = TRUE)
    dev.off()
    print(pltname)

    # Add intersection and remove empty sets:
    genesets[['Int']] = intersect(genesets[[1]], genesets[[2]])
    dl = c(lapply(genesets, length))
    keep = names(dl)[dl > 0]
    genesets = genesets[keep]
    core.int = intersect(getGeneSet(pair[1], cmlist),
        getGeneSet(pair[2], cmlist))
    cat(pstr, '\n')
    cat(core.int, '\n')

    gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
        ordered_query=FALSE, multi_query=TRUE,
        sources = sources)
    gpdf = gp2.result$result
    gpdf = gpdf[gpdf$term_size < 500,]

    pltprefix = paste0(imgpref, '_enrplot_', pstr)
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=10)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE)
}

