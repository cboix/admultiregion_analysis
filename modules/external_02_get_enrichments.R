#!/usr/bin/R
# --------------------------------------------------
# Get the GO term enrichments for external datasets:
# Updated 12/01/2023
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(gprofiler2)

library(ComplexHeatmap)
library(circlize)
options(width=170)

source(paste0(sbindir, 'modules/auxiliary_gprofiler_functions.R'))

# Set the run arguments:
# ----------------------
# dataset = 'COVID19_Cell2021'
# graph_id = 'base'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    dataset = args[1]
    graph_id = args[2]
}

# Directories:
extdir = paste0(sdbdir, 'external_datasets/')
dsdir = paste0(extdir, dataset, '/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'external_', dataset, '_')


# Load in the modules data:
# -------------------------
commandArgs <- function(trailingOnly=TRUE){c(dataset, graph_id, FALSE)}
source(paste0(sbindir, 'modules/load_external_modules_data.R'))


# Run GO enrichments on all modules jointly, for plotting:
# --------------------------------------------------------
for (set in c('allgenes', 'coregenes')){
    print(graphpref)
    print(set)
    moduleenr.file = paste0(dsdir, 'module_enrichments_', set, '_', graphpref, '.rda')
    if (!file.exists(moduleenr.file)){
        topmod = mmap$module
        if (set == 'coregenes'){
            genesets = lapply(topmod, function(m){ 
                                  x = names(coremap[(coremap == m)]) 
                                  x[!(is.na(x))]
                                  })
        } else {
            genesets = lapply(topmod, function(m){ 
                                  x = names(genemap[(genemap == m)]) 
                                  x[!(is.na(x))]
                                  })
        }
        names(genesets) = mmap[topmod + 1, 'mname']

        # Remove empty sets:
        dl = c(lapply(genesets, length))
        keep = names(dl)[dl > 0]
        genesets = genesets[keep]

        sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
        NMOD = length(genesets)
        gpdf = c()
        for (i in 1:NMOD){
            print(i)
            gp2.result = try(gprofiler2::gost(genesets[keep[i]], organism='hsapiens',
                    ordered_query=FALSE, multi_query=FALSE, sources=sources))
            if ((class(gp2.result) != 'try-error') & !is.null(gp2.result)){
                g1 = gp2.result$result
                g1$module = i - 1
                gpdf = rbind(gpdf, g1)
            }
        }
        save(gpdf, genesets, file=moduleenr.file)
    }

    # Reduce to the top p-values for figures:
    # ---------------------------------------
    load(moduleenr.file)
    gpdf = gpdf[gpdf$term_size < 500,]
    gpdf = gpdf[gpdf$intersection_size >= 4,]
    gpdf = gpdf[order(gpdf$p_value), ]

    ntop = 10
    ll = lapply(unique(gpdf$module), function(x){
        head(gpdf[gpdf$module == x,], 1) })
    topdf = do.call(rbind, ll)
    pvalsdf = topdf[,c('query','module','p_value','term_size','query_size',
        'intersection_size', 'term_id','source','term_name')]
    names(pvalsdf)[1] = 'mname'

    toppvals.file = paste0(dsdir, 'module_enrichments_toppvals_',set,'_', graphpref, '_small.tsv')
    write.table(pvalsdf, toppvals.file, quote=F, row.names=F, sep="\t")
}

