#!/usr/bin/R
# ------------------------------------
# Plot specific pairs of correlations:
# Updated 02/09/2022
# ------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))

library(tidyr)
library(gprofiler2)

library(viridis)
library(ggpubr)
library(ggrastr)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
resdir = paste0(sdbdir, 'modules/resources/')
plotdir = paste0(imgdir, 'modules/')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
system(cmd)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


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

imgpref = paste0(plotdir, 'module_corrplots_', runset,'_')


# Functions + variables:
# ----------------------
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))
col_log10p = colorRamp2(c(0, 10), c('white','black'))
colg = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(50)
pcols = brewer.pal(12,'Paired')

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



# Load in full data:
# ------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# # Load in the full pseudobulk data for these subtypes (core genes only):
# # NOTE: From modules_degenr_03_plot_pseudobulk_modules.R
# # ------------------------------------------------------
useset = 'coregenes_'
scores.file = paste0(moddir, 'module_pseudobulk_scores_', 
    useset, fullpref, '.tsv.gz')
scdf = read.delim(gzfile(scores.file), sep="\t")
scdf = merge(scdf, mmap, all.x=TRUE)
scdf$mod = paste0('M', scdf$module)
scdf$module = NULL
scdf$mname = NULL

subtype.cols = tcols[sort(unique(scdf$cell_type_high_resolution))]


# Set the comparisons of interest:
# --------------------------------
if (runset == 'Ast'){
    mod.comp = list(c(3, 17), c(0, 12), c(6, 27), c(8,15), c(6,13))
} else if (runset == 'Mic_Immune'){
    mod.comp = list(c(1, 20), c(1, 11), c(1, 15),
        c(1, 2), c(11, 2), c(11, 25))
} else {
    mod.comp = list(c(0, 1))
}

# Get sets function, default to expanded genes:
getGeneSet <- function(module, usemap=coremap){
    x = names(usemap)[usemap == module]
    return(x)
}


for (numpair in mod.comp){
    pair = paste0('M',numpair)
    pstr = paste(pair, collapse='_')
    subdf = scdf[scdf$mod %in% pair,]
    subdf = merge(subdf, unique(metadata[,c('projid','cogdxad')]))
    subdf = spread(subdf, mod, score)
    if (runset == 'Mic_Immune'){
        subdf = subdf[subdf$cell_type_high_resolution != 'T cells',]
    }

    attr = 'cog'
    if (attr == 'cog'){
        use.attr = 'cogdxad'
        use.cols = c('AD'='red', 'CTRL'='grey75')
    } else {
        use.attr = 'cell_type_high_resolution'
        use.cols = subtype.cols
    }

    gp = ggplot(subdf, aes_string(pair[1], pair[2], color=use.attr)) + 
        geom_point(cex=.25) +
        geom_smooth(method='lm', color='black',lwd=.5) + 
        scale_color_manual(values=use.cols) + 
        stat_cor(color='black', cex=3, output.type='text', label.sep='\n', label.y.npc=.95) + 
        theme_pubr() + theme(legend.position='none')
    gp2 = rasterize(gp, layers='Point', dpi=450)

    pltprefix = paste0(imgpref, 'corrplot_', attr, '_', pstr)
    saveGGplot(gp2, pltprefix, w=2, h=2)


    # Make the genesets from this pair:
    # ---------------------------------
    genesets = lapply(numpair, getGeneSet)
    names(genesets) = pair

    gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
        ordered_query=FALSE, multi_query=TRUE,
        sources = sources)
    gpdf = gp2.result$result
    gpdf = gpdf[gpdf$term_size < 500,]

    pltprefix = paste0(imgpref, 'enrplot_', pstr)
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=10)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE)

}


