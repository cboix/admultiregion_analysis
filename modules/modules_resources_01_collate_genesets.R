#!/usr/bin/R
# ---------------------------------------------------------
# Load all of the modules and save gene sets as a resource:
# Updated 01/12/2022
# ---------------------------------------------------------
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
runlist = c('Ast', 'Mic_Immune', 'Vasc_Epithelia', 
    'Oli', 'Opc', 'Exc', 'Inh', 'All')
# , 'Glial', 'All')


# Load each of the datasets:
# --------------------------
nodedf.list = list()
coremap.list = list()
genemap.list = list()
modnames.list = list()
dedf.list = list()
statsdf.list = list()
covar.list = list()
for (runset in runlist){
    cat(runset, '\n')
    commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
    source(paste0(sbindir, 'modules/load_modules_degenr.R'))
    # Add ngenes to modnames dataframe:
    ngenes = table(nodedf$leiden)
    mmap$ngenes = ngenes[as.character(mmap$module)]
    # Enrichments in covariates:
    hgann.file = paste0(moddir, 'module_covariate_hgdf_', 'coregenes_', fullpref, '.tsv')
    covardf = read.delim(hgann.file, sep="\t")
    covardf$set = runset
    # DE tables:
    dein.cols = c('gene','module','dkey','gset', 'logFC_nb','coef_mast', 'log10p_nm')
    deout.cols = c('gene','module', 'advar', 'DE', 'logFC_nb','coef_mast','joint_p')
    if (runset == 'All'){
        dein.cols = c('gene','module','dkey','gset')
        deout.cols = c('gene','module', 'advar', 'DE')
    } 
    dedf.list[[runset]] = dedf[,dein.cols]
    names(dedf.list[[runset]]) = deout.cols
    nodedf = nodedf[,-1]
    nodedf.list[[runset]] = nodedf # Full gene attributes
    coremap.list[[runset]] = coremap # Core genes
    genemap.list[[runset]] = genemap # All assigned genes
    modnames.list[[runset]] = mmap # Module names
    statsdf.list[[runset]] = statsdf # Enrichments
    covar.list[[runset]] = covardf
}


# Save these datasets as resources:
# ---------------------------------
respref = paste0(resdir, 'modules_resource_')
saveRDS(nodedf.list, file=paste0(respref, 'nodedf.Rds'))
saveRDS(coremap.list, file=paste0(respref, 'coremap.Rds'))
saveRDS(genemap.list, file=paste0(respref, 'genemap.Rds'))
saveRDS(modnames.list, file=paste0(respref, 'modnames.Rds'))
saveRDS(statsdf.list, file=paste0(respref, 'statsdf.Rds'))
saveRDS(dedf.list, file=paste0(respref, 'dedf.Rds'))
saveRDS(covar.list, file=paste0(respref, 'covariate_enrichments.Rds'))

# Collate covariates:
hgdf = do.call(rbind, covar.list)
labdf = hgdf[hgdf$cls == 1,]
labdf = merge(labdf, aggregate(p.value ~ covariate + mname + set, labdf, min))
labdf = labdf[labdf$p.value < 1e-2,]

xlfile = 'multiRegion/Supplementary_Table_6_module_covariate_enrichments.xlsx'
tsvfile = 'multiRegion/Supplementary_Table_6_module_covariate_enrichments.tsv'
tsvfile = 'multiRegion/test.tsv'
write.table(labdf, file=tsvfile, quote=F, sep="\t", row.names=F)

require(openxlsx)
write.xlsx(labdf, file=xlfile, rowNames=F)

