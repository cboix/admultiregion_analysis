#!/usr/bin/R
# ---------------------------------------------------
# Calculate enrichments for the final aggregated DEGs
# vs. the final modules (+ save resources for plots).
# Updated 02/18/2021
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
options(width=175)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Directories:
moddir = paste0(sdbdir, 'modules/')
regdir = paste0(sdbdir, 'dereg/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'module_vs_degs_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Load in the DEGs data:
# ----------------------
setfile = paste0(sdbdir, 'runsets_DEG_analyses.rds')
setdf = readRDS(file=setfile)
sets = unique(setdf$setid)
sets = c('Ast'='Ast_Ast','Mic_Immune'='Mic_Immune_Mic',
    'Opc'='Opc_Opc','Oli'='Oli_Oli','Inh'='Inh_Inh',
    'Exc'='Exc_Exc','All'='All_All')


for (use.core in c(TRUE, FALSE)){
    if (use.core){
        uselist = cmlist
        imgpref = paste0(plotdir, 'module_vs_degs_core_')
        dbpref = paste0(moddir, 'modenr_core.on_aggregated_fullset.')
    } else {
        imgpref = paste0(plotdir, 'module_vs_degs_all_')
        dbpref = paste0(moddir, 'modenr_all.on_aggregated_fullset.')
        uselist = gmlist
    }

    # Process each runset:
    for (i in 1:length(sets)){
        set = sets[i]  # DEGs
        runset = names(sets)[i]  # Modules
        cat(set,'\n')
        usemap = uselist[[runset]]

        # Load in data and rank by allregions:
        if (runset == 'All'){
            desuff = 'All_All_allconditions'
            full.rds = paste0(regdir, 'allmethods.', desuff, '.merged.rds')
            fulldf = readRDS(full.rds)

        } else {
            full.rds = paste0(regdir, 'aggregated_fullset.', set, '.rds')
            fulldf = readRDS(full.rds)
        }


        fulldf = fulldf[fulldf$gene %in% names(usemap),]
        fulldf$module = usemap[fulldf$gene]

        # Keep modules with at least 10 genes:
        mingenes = 10
        tab = table(cmlist[[runset]])
        kept.modules = names(tab[tab >= mingenes])
        fulldf = fulldf[fulldf$module %in% kept.modules,]

        # Get summary top DEGs for each condition:
        smdf = agg.rename(gene ~ module + col_nm + path + region,
            fulldf, function(x){
                paste(head(as.character(x), 3), collapse=", ")}, 'genes')

        # Count enrichment occurrences:
        ctdf = agg.rename(gene ~ module + col_nm + path + region, fulldf, length, 'count')
        ctdf = merge(ctdf, agg.rename(gene ~ path + region + col_nm, fulldf, length, 'tot.sign'))
        ctdf = merge(ctdf, agg.rename(gene ~ module + path + region, fulldf, length, 'tot.mod'))
        ctdf = merge(ctdf, agg.rename(gene ~ path + region, fulldf, length, 'total'))

        # Perform hypergeometric test for enrichment:
        df = ctdf[,c('count','tot.sign','tot.mod','total')]
        ctdf$p <- apply(df, 1, run.hyper)
        ctdf$log2FC = with(ctdf, log2((count / tot.sign) / (tot.mod / total)))

        # Prune (remove contam. modules)
        ctdf = ctdf[ctdf$col_nm != 0,]  # Look at only depleted / enriched
        if (runset == 'Ast'){ ctdf = ctdf[ctdf$module != 5,] }

        # Annotate and order:
        ctdf = merge(ctdf, rmap[rmap$runset == runset, c('module','mname')])
        ctdf = merge(ctdf, smdf, all.x=TRUE)
        ctdf$p.adj = p.adjust(ctdf$p, 'BH')
        ctdf = ctdf[order(ctdf$p),]
        show.cols = c('log2FC','p.adj','path','region','col_nm','mname', 'genes')
        head(ctdf[, show.cols], 50)

        # Save this table of enrichments:
        enr.rds = paste0(dbpref, set, '.', runset, '.rds')
        enr.tsv = paste0(dbpref, set, '.', runset, '.tsv.gz')
        saveRDS(ctdf, file=enr.rds)
        write.table(ctdf, file=gzfile(enr.tsv), quote=F, row.names=F, sep="\t")
    }
}

