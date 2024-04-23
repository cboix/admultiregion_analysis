#!/usr/bin/R
# ------------------------------------------
# Aggregate the DEGs runs for multiple runs:
# Updated: 01/04/22
# ------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(circlize)
library(ComplexHeatmap)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Get a list of all differential runs:
# -------------------------------------------------------
rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_joint_runlist.tsv'), header=T)
rundf$prefstr = with(rundf, gsub("[()/]", "_", paste(celltype, subtype, region, path, sep="_")))
rundf$setid = with(rundf, gsub("[/]", "_", paste(celltype, subtype, sep="_")))

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
    outfile = paste0(regdir, 'allmethods.', x, '.merged.rda')
    1 * file.exists(outfile) })
table(rundf$merged)
head(rundf[rundf$merged == 0,])

# Reduce to this level (celltype + subtype):
sets = unique(rundf$setid)
length(sets)


# Select runs to use (allregions, main ct):
# -----------------------------------------
for (i in 1:length(sets)){
    set = sets[i]
    cat(set,'\n')
    full.file = paste0(regdir, 'aggregated_fullset.', set, '.rda')
    full.rds = paste0(regdir, 'aggregated_fullset.', set, '.rds')
    full.tsv = paste0(regdir, 'aggregated_fullset.', set, '.tsv.gz')
    if (!file.exists(full.file)){
        # Select runs:
        selrundf = rundf[rundf$setid == set,]

        # Aggregate the DEGs and results across all runs:
        # -----------------------------------------------
        kept.cols = c("gene","pc", "col_nm","log10p_nm","path",
            'region', "logFC_nb","p_nb","padj_nb",
            "coef_mast","p_mast","padj_mast")
        totnsigdf = NULL
        fulldf    = NULL
        for (j in 1:nrow(selrundf)){
            prefstr = selrundf$prefstr[j]
            aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
            if (file.exists(aggrda)){
                load(aggrda)  # Loads aggdf, nsig
                print(prefstr)
                # Concatenate results:
                aggdf$path = nsig[1,'path']
                aggdf$region = nsig[1,'region']
                fulldf = rbind(fulldf, aggdf[,kept.cols])
                # Pad nsig:
                nsigdf = data.frame(nsig)
                if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
                if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
                totnsigdf = rbind(totnsigdf, nsigdf)
            }
        }
        save(fulldf, totnsigdf, file=full.file)

        # Reduced version:
        finaldf = fulldf[,c('gene','path','region','col_nm','log10p_nm','logFC_nb','coef_mast')]
        saveRDS(finaldf, file=full.rds)
        write.table(finaldf, gzfile(full.tsv), quote=F, row.names=F, sep="\t")
    } 
}

setdf = unique(rundf[,c('celltype','subtype','setid')])
setdf$tag = with(setdf, paste0(gsub("_", " ",subtype), " (", gsub("_", "/", celltype), ")"))

setfile = paste0(sdbdir, 'runsets_DEG_analyses.rds')
saveRDS(setdf, file=setfile)




# Repeat for sex*AD DEGs:
# -----------------------
rundf = read.delim(paste0(sdbdir, 'DEG_multiRegion_SI_ACE_runlist.tsv'), header=T)
rundf = rundf[rundf$path == 'msex',]
rundf$prefstr = with(rundf, gsub("[()/]", "_", paste(celltype, subtype, region, path, sep="_")))
rundf$setid = with(rundf, gsub("[/]", "_", paste(celltype, subtype, sep="_")))

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
    outfile = paste0(regdir, 'allmethods.', x, '.merged.rda')
    1 * file.exists(outfile) })
table(rundf$merged)
head(rundf[rundf$merged == 0,])

# Reduce to this level (celltype + subtype):
sets = unique(rundf$setid)
length(sets)


# Select runs to use (allregions, main ct):
# -----------------------------------------
full.file = paste0(regdir, 'aggregated_sexAD.all_celtypes.rda')
full.rds = paste0(regdir, 'aggregated_sexAD.all_celtypes.rds')
full.tsv = paste0(regdir, 'aggregated_sexAD.all_celtypes.tsv.gz')
if (!file.exists(full.file)){
    # Aggregate the DEGs and results across all runs:
    # -----------------------------------------------
    totnsigdf = NULL
    fulldf    = NULL
    for (j in 1:nrow(rundf)){
        prefstr = rundf$prefstr[j]
        aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
        if (file.exists(aggrda)){
            load(aggrda)  # Loads aggdf, nsig
            print(prefstr)
            # Concatenate results:
            aggdf$path = nsig[1,'path']
            aggdf$region = nsig[1,'region']
            fulldf = rbind(fulldf, aggdf)
            # Pad nsig:
            nsigdf = data.frame(nsig)
            if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
            if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
            totnsigdf = rbind(totnsigdf, nsigdf)
        }
    }
    save(fulldf, totnsigdf, file=full.file)

    # Reduced version:
    final.cols = c("path", "region", "gene",
        "logFC_nradAD_nb", "logFC_msex:nradAD_nb",
        "coef_nradAD_mast", "coef_msex:nradAD_mast", 
        "log10p_nradAD_nm", "log10p_msex:nradAD_nm",
        "col_nradAD_nm", "col_msex:nradAD_nm")
    finaldf = fulldf[,final.cols]
    saveRDS(finaldf, file=full.rds)
    write.table(finaldf, gzfile(full.tsv), quote=F, row.names=F, sep="\t")
} 

