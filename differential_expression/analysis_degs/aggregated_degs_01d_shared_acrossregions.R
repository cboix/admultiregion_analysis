#!/usr/bin/R
# ----------------------------------------------------------
# Get a list of shared DEGs across regions for each ct/path:
# and make plots for amt shared and consistency
# Updated: 12/20/22
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)

library(ggrepel)
library(ggplot2)
library(ggpubr)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Load regional differential results, process to sets for each cell type:
# -----------------------------------------------------------------------
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
sets = c('Ast_Ast','Exc_Exc','Inh_Inh','Opc_Opc','Mic_Immune_Mic','Oli_Oli')
pctdf = c()
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    shared.degs.rda = paste0(regdir, mstr, '.sharedDElists.rda')
    load(fullaggrda)

    uplist = list()
    dwlist = list()
    nup = c()
    ndw = c()
    for (set in sets){
        print(set)
        setdf = setdflist[[set]]
        setdf = setdf[setdf$col_nm != 0,]
        aggdf = agg.rename(pc ~ gene + col_nm, setdf, length, 'count')
        nup[set] = sum(aggdf$col_nm == 2)
        ndw[set] = sum(aggdf$col_nm == 1)
        # NOTE: Threshold by variability of pc - percent non-zero? 
        uplist[[set]] = aggdf$gene[(aggdf$count >= 4) & (aggdf$col_nm == 2)]
        dwlist[[set]] = aggdf$gene[(aggdf$count >= 4) & (aggdf$col_nm == 1)]
    }
    # Write up and down shared lists:
    save(uplist, dwlist, file=shared.degs.rda)

    # Aggregate statistics on sharing:
    lup = sapply(uplist, length)
    ldw = sapply(dwlist, length) 
    pup = lup / nup
    pdw = ldw / ndw
    sdf = data.frame(
        set=names(pup), path=path,
        lup=lup, ldw=ldw,
        nup=nup, ndw=ndw,
        pup=pup, pdw=pdw)
    pctdf = rbind(pctdf, sdf)
}


# Plot these statistics as a barplot:
gp = ggplot(pctdf, aes(set, pup, fill=set)) + 
    facet_wrap(~path) + 
    geom_bar(stat='identity', color=NA) + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0))+
    labs(y='% of DEGs that are up in 4+ regions', x='Set') +
    theme_pubr() + coord_flip()

pltprefix = paste0(imgpref, 'allmethods_shared4_allmajor')
saveGGplot(gp, pltprefix, h=4, w=6)


# Regression on percent shared by region, pathology, cell type:
# -------------------------------------------------------------
# TODO: (note that on corr / euc. may be better)

