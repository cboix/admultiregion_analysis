#!/usr/bin/R
# ---------------------------------------------------------------
# Load all of the modules and save pseudobulk data as a resource:
# Updated 01/13/2022
# ---------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


# Directories:
moddir = paste0(sdbdir, 'modules/')
resdir = paste0(moddir, 'resources/')
cmd = paste('mkdir -p', moddir, resdir)
system(cmd)


# Runsets that we can run/have pre-processed results:
# ---------------------------------------------------
graph_id = 'boot'
runlist = c(
    'Ast', 'Mic_Immune', 'Vasc_Epithelia', 
    'Oli', 'Opc', 'Exc', 'Inh', 'All')
# 'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons')
# , 'Glial', 'All')



# Load each of the datasets:
# --------------------------
ps.scores.list = list()
ps.data.list = list()
mod.enr.list = list()
for (runset in runlist){
    cat(runset, '\n')
    commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
    source(paste0(sbindir, 'modules/load_modules_degenr.R'))

    # Load in the full pseudobulk data for these subtypes:
    # NOTE: From modules_degenr_03_plot_pseudobulk_modules.R
    # NOTE: Creating reduced version as a resource for website
    # --------------------------------------------------------
    useset = 'coregenes_'
    runscores.file = paste0(moddir, 'module_pseudobulk_scores_', 
        useset, fullpref, '.tsv.gz')
    scdf = read.delim(gzfile(runscores.file), sep="\t")
    scdf = merge(scdf, mmap)
    scdf$ptype = NULL
    scdf$mname = NULL 
    # Get number of genes (to allow filtering scores with few genes):
    ngenes = table(nodedf$leiden)
    scdf$ng = ngenes[as.character(scdf$module)]
    ps.scores.list[[runset]] = scdf

    # Make as a psbulk dataset:
    # -------------------------
    scdf$set = runset
    scdf = merge(scdf, mmap)
    scdf$ptype = with(scdf, paste(projid, cell_type_high_resolution, region, set, sep="_"))
    # NOTE: Unique may take a while, speed up somehow or precompute as ps.data?
    umeta = unique(scdf[,c('ptype','projid','cell_type_high_resolution','region','set', 'ncell')])
    pmat = pivot.tomatrix(scdf[,c('mname','ptype','score')], 'ptype', 'score')
    pmat = pmat[, umeta$ptype]
    ps.data = list('mat'=pmat, 'meta'=umeta)
    ps.data.list[[runset]] = ps.data

    # Load in set of top by p-value functional enrichments for each module:
    # ---------------------------------------------------------------------
    set = 'coregenes'
    toppvals.file = paste0(moddir, 'module_enrichments_toppvals_',
        set, '_', fullpref, '_small.tsv')
    pvalsdf = read.delim(toppvals.file, sep="\t")
    # pvalsdf$nc = nchar(pvalsdf$term) # For filtering long terms
    pvalsdf = merge(mmap, pvalsdf)
    mod.enr.list[[runset]] = pvalsdf
}


# Save these datasets as resources:
# ---------------------------------
respref = paste0(resdir, 'modules_resource_')
saveRDS(ps.scores.list, file=paste0(respref, 'psbulkScores.Rds'))
saveRDS(ps.data.list, file=paste0(respref, 'psbulkDatasets.Rds'))
saveRDS(mod.enr.list, file=paste0(respref, 'moduleEnrichments.Rds'))
