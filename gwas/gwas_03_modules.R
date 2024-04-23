#!/usr/bin/R
# --------------------------------------------
# Plot the expression patterns for GWAS genes:
# Updated 12/17/2021
# --------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'gwas/')
imgpref = paste0(plotdir, 'gwas_')
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])
topdf = read.delim('Annotation/ADGWAS_topct_multiregion_121721.tsv', sep="\t")


# Load in all gene assignments:
# -----------------------------
modlist.rda = paste0(moddir, 'aggregated_module_gene_mapping.rda')
if (!file.exists(modlist.rda)){
    runlist = c('Ast','Oli','Opc','Mic_Immune','Vasc_Epithelia',
                'Inh','HCneurons', 'ECneurons', 'THneurons', 'CTXneurons')
    graph_id = 'boot'

    genedf = c()
    coredf = c()
    for (runset in runlist){
        commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
        source(paste0(sbindir, 'modules/load_modules_degenr.R'))

        mdf = data.frame(gene=names(genemap), module=genemap, set=runset)
        mdf$mname = mmap[mdf$module + 1,'mname'] 
        genedf = rbind(genedf, mdf)
        
        cdf = data.frame(gene=names(coremap), module=coremap, set=runset)
        cdf$mname = mmap[cdf$module + 1,'mname'] 
        coredf = rbind(coredf, cdf)
    }
    coredf$ct = coredf$set
    coredf$ct[coredf$set %in% names(exc.sets)] = 'Exc'
    genedf$ct = genedf$set
    genedf$ct[genedf$set %in% names(exc.sets)] = 'Exc'
    save(coredf, genedf, file=modlist.rda)
} else {
    load(modlist.rda)
}


# Merge the top-expressed hit with modules:
# -----------------------------------------
corehit = merge(coredf, topdf)

allhit = merge(coredf, topdf[,'gene', drop=F])



# Merge all hits with modules (e.g. ignore expr, only if gene in module core):
# ----------------------------------------------------------------------------


