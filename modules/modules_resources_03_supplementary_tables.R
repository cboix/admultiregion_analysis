#!/usr/bin/R
# ------------------------------------------------------
# Load all of the modules and save supplementary tables:
# Updated 06/20/2022
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(openxlsx)

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

# Load pre-processed:
respref = paste0(resdir, 'modules_resource_')
coremap.list = readRDS(file=paste0(respref, 'coremap.Rds'))
genemap.list = readRDS(file=paste0(respref, 'genemap.Rds'))

# Save the module lists for supplementary table:
# ----------------------------------------------
xlfile = 'multiRegion/Supplementary_Table_5_gene_expression_modules.xlsx'
rdsfile = 'multiRegion/Supplementary_Table_5_gene_expression_modules.Rds'
dflist = list()
for (runset in runlist){
    print(runset)
    coremap = coremap.list[[runset]]
    genemap = genemap.list[[runset]]
    df = data.frame(gene=names(genemap), module=genemap)
    df$is.core.gene = df$gene %in% names(coremap)
    df$set = runset
    df = df[order(df$gene),]
    df = df[order(!df$is.core.gene),]
    df = df[order(df$module),]
    dflist[[runset]] = df

    if (runset == runlist[[1]]) {
        write.xlsx(df, file=xlfile, sheetName=runset, rowNames=FALSE)
    } else { 
        write.xlsx(df, file=xlfile, sheetName=runset, append=TRUE, rowNames=FALSE)
    }
}

write.xlsx(dflist, file=xlfile)

saveRDS(dflist, file=rdsfile)




