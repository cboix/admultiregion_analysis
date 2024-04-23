#!/usr/bin/R
# -------------------------------------------------------------
# Replot the modules representations for extended data figures:
# Updated 11/28/2021
# -------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


library(tidyr)
library(viridis)
library(igraph)
library(gprofiler2)

library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'repr_')
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
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Code to look at genes:
# ----------------------
for (i in 0:max(coremap)){
    cat(i, '\t')
    genes = names(coremap)[coremap == i]
    genes = sort(genes)
    cat(genes)
    cat('\n')
}

gene = 'HIF1A'
module = genemap[gene]

module = 11
path = 'cogdxad'

subdf = dedf[(dedf$dkey == path) & (dedf$module == module),]
allgenes = sort(names(coremap)[coremap == module])
upgenes = sort(unique(subdf$gene[subdf$gset == 'Up']))
downgenes = sort(unique(subdf$gene[subdf$gset == 'Down']))

cat('All:\t', allgenes,'\n')
cat('Up:\t', upgenes,'\n')
cat('Down:\t', downgenes, '\n')


sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
gp2.result = gprofiler2::gost(upgenes, organism='hsapiens',
                              ordered_query=FALSE, multi_query=FALSE,
                              sources = sources, evcodes=TRUE)
gpdf = gp2.result$result
if (!is.null(gpdf)){
    gpdf = gpdf[order(gpdf$p_value),]
    # sgpdf = gpdf[gpdf$source %in% c('KEGG','REAC','WP','CORUM'),]
    gpdf$pstr = with(gpdf, paste0(term_name, '\t(', sprintf("%0.1e", p_value), ')\t', intersection))
    out = head(gpdf$pstr, 10)
    cat(paste(out, collapse="\n"))
}





