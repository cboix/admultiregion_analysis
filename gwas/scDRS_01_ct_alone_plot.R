#!/usr/bin/R
# ----------------------------------
# Prelim plots scDRS
# Updated 10/26/2022
# ----------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)
library(viridis)
library(ggplot2)
library(ggpubr)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
scddir = paste0(sdbdir, 'scDRS/')
plotdir = paste0(imgdir, 'gwas/')
imgpref = paste0(plotdir, 'scDRS_')
cmd = paste('mkdir -p', plotdir, scddir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Read scDRS summary:
# -------------------
df = read.delim(paste0(scddir, 'group.scores-allGWAS.tsv.gz'), header=T)
for (fdr in c('0.05','0.1','0.2')){
    fdrvar = paste0('pct', fdr)
    df[fdrvar] = df[paste0('n_fdr_', fdr)] / df$n_cell
}
ad.gwas = "PASS_Alzheimers_Jansen2019"

for (ctcol in c('celltype', 'major.celltype')){
    subdf = df[df$ctcol == ctcol,]
    cmat = pivot.tomatrix(subdf[,c('group','gwas','pct0.1')], 'group', 'pct0.1')

    ux = 1.5
    plt = Heatmap(cmat,
        col=colb,
        name='% cells > 0.1',
        width = ncol(cmat)*unit(ux, "mm"), 
        height = nrow(cmat)*unit(ux, "mm"),
        border_gp = gpar(col="black", lty = 1, lwd=.5),
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
    )

    h = 2.25 + 1 / 15 * nrow(cmat)
    w = 5 + 1 / 15 * ncol(cmat)
    pltprefix = paste0(imgpref, 'gwas_heatmap_', ctcol)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}


# Plot only major cell types for AD:
# ----------------------------------
mcts = unique(cellmeta$major.celltype)
mcts = sort(mcts)
subdf = df[df$group %in% mcts,]
subdf = subdf[subdf$gwas == ad.gwas,]
sub.tcols = c(tcols, major.col)[mcts]
subdf = subdf[order(subdf$pct0.1),]
subdf$group = factor(subdf$group, levels=subdf$group)
mcts = rev(levels(subdf$group))

gp = ggplot(subdf, aes(group, pct0.05, fill=group)) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=sub.tcols, 'Celltype') + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(y='% cells with FDR > 0.05', x='Celltype') + 
    theme_pubr() + coord_flip() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'ad_major_overall_barplot')
saveGGplot(gp, pltprefix, w=2.75, h=2)


# Plot microglia subtypes for AD:
# -------------------------------
cts = unique(cellmeta[cellmeta$major.celltype == 'Mic/Immune', 'cell_type_high_resolution'])
cts = c(mcts, sort(cts))
subdf = df[df$group %in% cts,]
subdf = subdf[subdf$gwas == ad.gwas,]
sub.tcols = c(tcols, major.col)[cts]
subdf$group = factor(subdf$group, levels=rev(cts))

gp = ggplot(subdf, aes(group, pct0.05, fill=group)) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=sub.tcols, 'Celltype') + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(y='% cells with FDR > 0.05', x='Celltype') + 
    theme_pubr() + coord_flip() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'ad_overall_barplot')
saveGGplot(gp, pltprefix, w=2.75, h=2.75)


# Load in AD results by region; plot:
# -----------------------------------
addf = read.delim(paste0(scddir, 'group.scores-ADGWAS.tsv.gz'), header=T)

for (fdr in c('0.05','0.1','0.2')){
    fdrvar = paste0('pct', fdr)
    addf[fdrvar] = addf[paste0('n_fdr_', fdr)] / addf$n_cell
    subdf = addf[addf$ctcol %in% c('cthr_region', 'mct_region'),]
    subdf$cthr = sub("_[A-Z]+$","", subdf$group)
    subdf$region = sub('.*_', "", subdf$group)

    mat = pivot.tomatrix(subdf[,c('cthr','region',fdrvar)], 'region', fdrvar)
    cmat = mat[cts,reg.order[-1]]

    ux = 1.5
    plt = Heatmap(cmat,
        # col=colb,
        col=rev(colspec),
        name=paste0('% cells > ', fdr),
        use_raster=FALSE,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        row_split = ifelse(rownames(cmat) %in% mcts, 'Celltype', 'Subtype'),
        width = ncol(cmat)*unit(ux, "mm"), 
        height = nrow(cmat)*unit(ux, "mm"),
        border_gp = gpar(col="black", lty = 1, lwd=.5),
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
    )

    h = .5 + 1 / 15 * nrow(cmat)
    w = 1.5 + 1 / 15 * ncol(cmat)
    pltprefix = paste0(imgpref, 'gwas_adregion_heatmap_', fdr)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}


