#!/usr/bin/R
# --------------------------------------------------------
# Explore AD diagnosis vs. scores from different modules:
# Do different modules explain different parts of AD diag?
# Updated 12/19/2021
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(glasso)
library(cglasso)  # Conditional lasso
library(ggplot2)
library(ggpubr)

# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
# TODO: Opt to include neurons or not? Or just add later...
# commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Turn the score dataframe into a matrix:
# NOTE: At the runset level or at the subtype level:
# --------------------------------------------------
mingenes = 10
modlevel = 'subtype'
modlevel = 'runset'
if (modlevel == 'subtype'){
    mat = pivot.tomatrix(scoredf[scoredf$ng >= mingenes, c('pr','cm','score')], 'pr','score')
    scmat = log(t(mat) + 1e-4)
    ut = 1
} else {
    mat = pivot.tomatrix(runscdf[runscdf$ng >= mingenes, c('pr','rm','score')], 'pr','score')
    scmat = log(t(mat) + 1e-4)
    ut = 3
}


# Make metadata data.frame for the score matrix:
# ----------------------------------------------
mdf = data.frame(pr=rownames(scmat),
                 projid = sapply(rownames(scmat), function(x){sub("-.*","",x)}),
                 region = sapply(rownames(scmat), function(x){sub(".*-","",x)}))
mdf = merge(mdf, metadata[,c('projid','region','rind','msex','age_death', 'cogdxad','nrad')])
rownames(mdf) = mdf$pr
mdf = mdf[rownames(scmat),]
mdf$region = as.character(mdf$region)


# Plot different scores against each other:
# And against AD:
# -----------------------------------------
s1 = 'Mic_Immune-11'   # APOE / C1 etc.
# s1 = 'Mic_Immune-4'   # HIF1A
# s2 = 'Mic_Immune-18'  # Gly
s2 = 'Mic_Immune-25'  # Gly
# s1 = 'Ast-0'   # APOE / C1 etc.
# s2 = 'Ast-6'  # Gly
# s2 = 'Mic_Immune-20'  # IFI44L
mdf$s1 = scmat[,s1]
mdf$s2 = scmat[,s2]

mlong = gather(mdf[,c('nrad','region','s1','s2','projid')], 
               module, val, -region, -projid, -nrad)


# Alone vs. AD:
gp = ggplot(mlong, aes(module, val, fill=nrad)) + 
    geom_boxplot() + 
    labs(x=paste0('s1: ', s1, ', s2: ', s2)) +
    theme_pubr()
pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_scorecomp_boxplot')
ggsave(paste0(pltprefix, '.png'), gp, dpi=400, units='in', width=6, height=5)

# Scatter:
gp = ggplot(mdf, aes(s1, s2, color=nrad)) + 
    geom_point() + 
    labs(x=paste0('s1: ', s1), y=paste0('s2: ', s2)) +
    theme_pubr()
pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_scorecomp_scatter')
ggsave(paste0(pltprefix, '.png'), gp, dpi=400, units='in', width=6, height=5)

# Scatter:
gp = ggplot(mdf, aes(s1 + s2, s2 - s1, color=nrad)) + 
    geom_point() + 
    labs(title=paste0('s1: ', s1, ', s2: ', s2), x='s1 + s2', y='s2 - s1') +
    geom_hline(yintercept=mean(mdf$s2, na.rm=T) - mean(mdf$s1, na.rm=T)) + 
    theme_pubr()
pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_scorecomp_ma')
ggsave(paste0(pltprefix, '.png'), gp, dpi=400, units='in', width=6, height=5)


# Perform lasso to determine best components for predicting AD:
# NOTE: Filtering on DE genes might be better
# -------------------------------------------------------------
library(glmnet)
cind = grep("Mic_Immune", colnames(scmat))
scols = colnames(scmat)[cind]

x = scmat[,scols]
x[is.na(x)] = 0
fit = glmnet(x, mdf$cogdxad, family='binomial')

coefficients(fit)[,0:10]





