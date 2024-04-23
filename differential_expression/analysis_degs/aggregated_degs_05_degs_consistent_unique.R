#!/usr/bin/R
# -----------------------------------------------------------
# Plot top DEGs that are consistent or unique across regions:
# Updated: 02/17/22
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)

library(ComplexHeatmap)
library(circlize)
print(version)
options(width=175)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Functions:
# ----------
plotSelGenes = function(fulldf, top.genes, zcol, mx, row.split=NULL, raster=FALSE){
    zcol = 'logFC_nb'
    z.top.mat = pivot.tomatrix(fulldf[fulldf$gene %in% top.genes, c('gene','region',zcol)], 'region', zcol)
    col.top.mat = pivot.tomatrix(fulldf[fulldf$gene %in% top.genes, c('gene','region','col_nm')], 'region', 'col_nm')
    z.top.mat[is.na(z.top.mat)] = 0
    col.top.mat[is.na(col.top.mat)] = 0
    reg.cols = c('allregions', 'AG','MT','PFC','EC','HC','TH')
    reg.cols = reg.cols[reg.cols %in% colnames(z.top.mat)]
    z.top.mat = z.top.mat[top.genes, reg.cols]
    col.top.mat = col.top.mat[top.genes, reg.cols]
    col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
    # z.top.mat = z.top.mat * (col.top.mat > 0)
    ux = 1.5
    plt = Heatmap(z.top.mat, name=zcol,
        use_raster=raster,
        col=col_fun,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        row_split = row.split,
        border_gp=gpar(color='black', lwd=.5),
        width = ncol(z.top.mat)*unit(ux, "mm"), 
        height = nrow(z.top.mat)*unit(ux, "mm")
    )
    return(plt)
}


# Get a list of all differential runs:
# -------------------------------------------------------
path = 'cogdxad'

setfile = paste0(sdbdir, 'runsets_DEG_analyses.rds')
setdf = readRDS(file=setfile)
sets = unique(setdf$setid)

# Load shared, for filtering:
sharedup.file = paste0(regdir, 'allmethods.allmajor.', path, '.sharedup.tsv')
shareddw.file = paste0(regdir, 'allmethods.allmajor.', path, '.shareddw.tsv')
up.rep = scan(sharedup.file, 'c', quiet=T)
dw.rep = scan(shareddw.file, 'c', quiet=T)
shared.genes = c(up.rep, dw.rep)

sets = c('Ast_Ast','Mic_Immune_Mic','Opc_Opc','Oli_Oli','Inh_Inh','Exc_Exc')

# Process each set:
for (set in sets){
    cat(set,'\n')
    # Load in data and rank by allregions:
    full.rds = paste0(regdir, 'aggregated_fullset.', set, '.rds')
    fulldf = readRDS(full.rds)
    fulldf = fulldf[fulldf$path == path,]
    mt.genes = unique(fulldf$gene[c(grep("^MTRNR", fulldf$gene), grep("^MT-", fulldf$gene))])
    # fulldf = fulldf[!(fulldf$gene %in% c(shared.genes, mt.genes)),] # Remove shared genes
    # fulldf = fulldf[!(fulldf$gene %in% c(shared.genes, mt.genes)),] # Remove shared genes
    fulldf$zdir = scale(fulldf$logFC_nb) + scale(fulldf$coef_mast)
    fulldf = fulldf[order(-abs(fulldf$zdir)),]
    fulldf = fulldf[order(fulldf$log10p_nm, decreasing=TRUE),]
    ctde = agg.rename(zdir ~ gene, fulldf[fulldf$region != 'allregions' 
        & fulldf$col_nm != 0,], 'length', 'nde')

    # Rank all regions genes aggregate with multiply-DE genes for top consistent:
    subdf = fulldf[fulldf$region == 'allregions' & fulldf$col_nm != 0,]
    subdf$rank = 1:nrow(subdf)
    subdf = merge(subdf, ctde)
    subdf = subdf[order(subdf$rank),]

    if (path %in% c('nft','plaq_n','plaq_d')){
        mx = 0.05 } else { mx = 1 }

    # Top N in 4/6 regions at least
    ntop = 5
    topdf = subdf[subdf$nde >= 4,]
    top.genes = c(head(topdf$gene[topdf$col_nm == 2], ntop),
        head(topdf$gene[topdf$col_nm == 1], ntop))
    # Top N in fewer than 4 regions at least
    ntop = 25
    topdf = subdf[subdf$nde < 4,]
    top.genes2 = c(head(topdf$gene[topdf$col_nm == 2], ntop),
        head(topdf$gene[topdf$col_nm == 1], ntop))
    # Consistent + regional jointly:
    row.split = c(rep('Consistent', length(top.genes)), rep('Regional', length(top.genes2)))
    plt = plotSelGenes(fulldf, c(top.genes, top.genes2), mx=mx, zcol ='logFC_nb', row.split=row.split)
    h = 2.25 + 1 / 15 * length(top.genes2)
    w = 5 + 1 / 15 * (length(reg.nomb) + 1)
    pltprefix = paste0(imgpref, 'allmethods.regional_', path, '.', set, '.top', ntop, '_genes_indv.heatmap')
    saveHeatmap(plt, pltprefix, w, h)

    # Score DEGs, averaging the three neocortex regions:
    # aggdf = fulldf[fulldf$region %in% reg.nomb, c('gene','path','zdir', 'region')]
    # aggdf$mult = ifelse(aggdf$region %in% c('AG','MT','PFC'), 1/3, 1)
    # aggdf$zdir = aggdf$mult * aggdf$zdir
    # aggdf = agg.rename(zdir ~ gene + path, aggdf, sum, 'score')

    # Top N overall:
    ntop = 8
    topdf = subdf
    top.genes = c(head(topdf$gene[topdf$col_nm == 2], ntop),
        head(topdf$gene[topdf$col_nm == 1], ntop))
    # Consistent + regional jointly:
    row.split = c(rep('Up', ntop), rep('Down', ntop))
    plt = plotSelGenes(fulldf, c(top.genes), mx=mx, zcol ='logFC_nb', row.split=row.split)
    h = 2.25 + 1 / 15 * length(top.genes)
    w = 5 + 1 / 15 * (length(reg.nomb) + 1)
    pltprefix = paste0(imgpref, 'allmethods.regional_', path, '.', set, '.top', ntop, '_genes_overall.heatmap')
    saveHeatmap(plt, pltprefix, w, h)
}



