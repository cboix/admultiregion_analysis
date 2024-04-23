#!/usr/bin/R
# ------------------------------------------
# Plot pathology measurements by region, ABC
# Updated 11/10/2023
# ------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(qvalue)
library(ComplexHeatmap)
library(circlize)
options(width=170)

# Directories:
plotdir = paste0(imgdir, 'metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Load ABC scores and add metadata:
# ---------------------------------
indmeta_tsv = 'Annotation/metadata_PFC_all_individuals_092520.tsv'
ext.meta = read.delim(indmeta_tsv, header=T)
ext.meta$kept.ind = ifelse(ext.meta$projid %in% kept.individuals, 'Our Cohort\n(48 individ.)', 'ROSMAP\n(Sept. 2020)')
abcdf = read.delim('Annotation/abc_scores_092023.tsv', header=T)
abcdf = merge(abcdf, ext.meta[,c('projid', 'niareagansc', 'cogdx')])

# Make metadata table for our samples:
metadf = unique(metadata[(metadata$projid %in% kept.individuals) & (
        metadata$region %in% reg.nomb), c('projid','rind','region')])
metadf = merge(metadf, abcdf)
metadf = merge(metadf, pqdf, all.x=TRUE)
metadf$region = factor(metadf$region, levels=c('EC','HC','MT','AG','PFC'))
metadf$nia_aa_sc = factor(metadf$nia_aa_sc)


# Plot pathology vs. ABC:
# -----------------------
pathmap = c('nft'='NFT','plaq_d'='Diffuse Plaque', 'plaq_n'='Neuritic Plaque')
for (path in names(pathmap)){
    # gp = ggplot(metadf, aes_string('nia_aa_sc', path, color='region')) + 
    #     scale_color_manual(values=reg.cols) + 
    gp = ggplot(metadf, aes_string('nia_aa_sc', path, fill='region')) + 
        scale_fill_manual(values=reg.cols) + 
        geom_boxplot(outlier.size=.5) + 
        scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) + 
        labs(x='ABC Score', y=pathmap[path]) + 
        theme_pubr()
    pltprefix = paste0(imgpref, 'region_path_abc.', path)
    saveGGplot(gp, pltprefix, w=4, h=3)
}

gp = ggplot(metadf, aes(factor(nia_aa_sc), plaq_n, color=region)) + 
    geom_boxplot(outlier.size=.5) + 
    # scale_y_log10() + 
    scale_color_manual(values=reg.cols) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) + 
    labs(x='ABC Score', y='NFT') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'region_path_abc.nft')
saveGGplot(gp, pltprefix, w=4, h=3)







