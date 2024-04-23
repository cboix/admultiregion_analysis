#!/usr/bin/R
# -----------------------------------------------------------
# Aggregate and plot chunked runs of MAST + RE for DGE
# Plot the comparison between regressions on different pathologies.
# Updated: 02/11/2021
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)
library(MAST)
library(data.table)
library(viridis)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# ------------------
# Load the metadata:
# ------------------
# Data directories:
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    # matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')
ststr = gsub("/","_",gsub(" ","_", subtype))
cellstr = gsub("/","_",gsub(" ","_", celltype))

# --------------------------------------------
# Load from various different regression runs:
# --------------------------------------------
ctlist = c('Ast','Mic_Immune','Oli','Opc','Vasc_Epithelia')
fns = c()
for (ct in ctlist){
    fns = c(fns, list.files(path=regdir, 
                            pattern=paste0(prefix, '.mastlmm_reg.allpath.allreg.major.', ct, '.*.rda')))
}

resdf = c()
for (fn in fns){
    base = sub(".*allreg.major.","", fn)
    ct = sub(".minor.*","", base)
    st = sub(".rda","", sub(".*minor.","", base))
    load(paste0(regdir, fn))
    if (!is.null(alldf)){
        alldf$ct = ct
        alldf$st =st
        cols = c('primerid','coef','padj','path','ct','st')
        resdf = rbind(resdf, alldf[,cols])
    }
}

FCTHRESH=0.02
# FCTHRESH=0.04
pcut = 0.05
resdf$col = 1 * (resdf$padj < pcut) * (1 + 1 * (resdf$coef > 0)) * (abs(resdf$coef) > ifelse(resdf$path == 'nrad', 1,1/50) * FCTHRESH)

nhitdf = aggregate(coef ~ primerid, resdf[resdf$col != 0,], length)
nhitdf = aggregate(coef ~ primerid + path, resdf[resdf$col != 0,], max)
uqgenes = nhitdf$primerid[nhitdf$coef == 1]
uqres = resdf[resdf$primerid %in% uqgenes,]
resdf$abscoef = abs(resdf$coef)
maxdf = aggregate(abscoef ~ primerid + path, resdf[resdf$col != 0,], max)
uqres = merge(resdf, maxdf)

pcols = brewer.pal(12, 'Paired')
nsigdf = agg.rename(coef ~ col + path + ct + st, resdf, length, 'count')
# nsigdf = agg.rename(coef ~ col + path + ct + st, uqres, length, 'count')
nsigdf = nsigdf[nsigdf$path %in% c('nft','plaq_d','plaq_n'),]
nsigdf = spread(nsigdf, col, count, fill=0)
names(nsigdf)[4:6] = c('none','down','up')
# names(nsigdf)[4:5] = c('down','up')

subdf = nsigdf[nsigdf$st %in% c('Opc','Oli','Mic','End','CAMs','Ast'),] 
# TODO: Shared vs. not:

# gplot = ggplot(subdf, aes(st, up, color=path, alpha=path)) + 
subdf = subdf[order(subdf$st),]
subdf$st = factor(subdf$st, levels=rev(unique(subdf$st)))
gplot = ggplot(subdf, aes(st, up)) + 
    facet_wrap(~path, ncol=1) + 
    geom_bar(position='dodge',stat='identity', fill=pcols[6]) + 
    geom_bar(data=subdf, aes(st, -down), position='dodge',stat='identity', fill=pcols[2]) + 
    labs(x='AD variable', y='Number of DEGs') + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, 'mastlmm_comparison_ndeg.allreg.allglial.png'), gplot, dpi=400, units='in', width=2.5,height=5)
ggsave(paste0(imgpref, 'mastlmm_comparison_ndeg.allreg.allglial.pdf'), gplot, width=2.5,height=5)

