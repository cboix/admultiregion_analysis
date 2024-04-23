#!/usr/bin/R
# ---------------------------------------------------
# Plot GWAS genes against modules + flag as DE / not:
# Updated 12/21/2022
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
options(width=175)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Directories:
moddir = paste0(sdbdir, 'modules/')
regdir = paste0(sdbdir, 'dereg/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'module_vs_gwas_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])


# Overall statistics for GWAS genes vs. modules:
# ----------------------------------------------
# Number of GWAS genes in core modules vs. total
incore = gwgenes[gwgenes %in% names(coremap)]
infull = gwgenes[gwgenes %in% names(genemap)]

# Mappings:
runsets = c('Ast', 'Mic_Immune', 'Oli', 'Opc', 'Exc', 'Inh', 'Vasc_Epithelia')
rlist = lapply(runsets, function(x){
    coremap = cmlist[[x]]
    incore = gwgenes[gwgenes %in% names(coremap)]
    # Module enrichment for GWAS genes:
    coregw = coremap[incore]
    cdf = data.frame(table(coremap))
    names(cdf) = c('module','nmod')
    cgwdf = data.frame(table(coregw))
    names(cgwdf) = c('module','ngw')

    cdf = merge(cgwdf, cdf, all.y=TRUE)
    cdf$ngw[is.na(cdf$ngw)] = 0
    cdf$ntot = length(coremap)
    cdf$ngwtot = length(incore)
    cdf$p = apply(cdf[,c('ngw','nmod','ngwtot','ntot')], 1, run.hyper)
    cdf = cdf[order(cdf$p),]
    cdf$genes = sapply(cdf$module, function(i){paste(sort(names(coregw)[coregw == i]), collapse=',')})
    cdf$runset = x
    return(cdf)
})
rdf = do.call(rbind, rlist)
rdf = rdf[order(rdf$p),]

# Print top modules / sets:
head(rdf[rdf$nmod >= 10,], 30)



