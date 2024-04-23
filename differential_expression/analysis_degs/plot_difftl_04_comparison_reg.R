#!/usr/bin/R
# -----------------------------------------------
# Comparison of DE methods - NB / MAST / Wilcoxon
# Updated: 05/25/21
# -----------------------------------------------
# Aggregate number of sign. genes:
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
print(version)

celltype = 'Ast'
subtype = 'Ast'
# subtype = 'Ast_DPP10'
region = 'allregions'

# Data loader:
commandArgs = function(x){ c(celltype, subtype, region)}
source(paste0(bindir, 'multiRegion/load_difftl_data.R'))
gc()

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'compdifftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# -----------------------------
# Load in the regression files:
# -----------------------------
prefstr = paste0(celltype,'_',subtype, '_', region)
fn1list = list.files(path=regdir, pattern=paste0('nebula_ruv.',prefstr, '.*rda'))
fn2list = list.files(path=regdir, pattern=paste0('mast.',prefstr, '.*rda'))
fn3list = list.files(path=regdir, pattern=paste0('wilcoxon',prefstr, '.*rda'))
nsigdf = c(); alldf = c();
for (fn in c(fn1list, fn2list, fn3list)){
    method = sub(paste0(".", prefstr,".*.rda"), "", fn)
    load(paste0(regdir, fn))
    region = nsig[1,'region']
    path = nsig[1,'path']
    ndf = data.frame(nsig)
    ndf$method = method
    if (!('X1' %in% colnames(ndf))){ ndf$X1 = 0 }
    if (!('X2' %in% colnames(ndf))){ ndf$X2 = 0 }
    if (is.null(nsigdf)){
        nsigdf = ndf 
    } else {
        nsigdf = rbind(nsigdf, ndf[,colnames(nsigdf), drop=F])
    }
    fulldf = fulldf[order(fulldf$p),]
    if (method == 'mast'){
        fulldf = fulldf[,c('gene','coef','p','fdr','val','col')]
        names(fulldf) = c('gene','logFC','p','q','pc','col')
    }
    fulldf = fulldf[,c('gene','p','logFC','q','pc','col')]
    fulldf$path = path
    fulldf$region = region
    fulldf$method = method
    fulldf$rank = 1:nrow(fulldf)
    alldf = rbind(alldf, fulldf)
}
names(nsigdf)[names(nsigdf) == 'X1'] = 'ndown'
names(nsigdf)[names(nsigdf) == 'X2'] = 'nup'
nsigdf$nup = as.numeric(nsigdf$nup)
nsigdf$ndown = as.numeric(nsigdf$ndown)

resdf = alldf[alldf$region == 'allregions',]
cat(head(resdf[resdf$col == 1 & resdf$path == 'nrad','gene'], 50), "\n\n")
cat(head(resdf[resdf$col == 2 & resdf$path == 'nrad','gene'], 50), "\n")


path = 'nrad'
resdf = alldf[alldf$region == 'allregions' & alldf$path == path,]
rlong = spread(resdf[,c('gene','logFC','method')], method, logFC)
clong = spread(resdf[,c('gene','col','method')], method, col)
names(clong)[2:3] = paste0(names(clong)[2:3], '_col')
rlong = merge(rlong, clong)

ggplot(rlong, aes(mast, nebula_ruv, col=factor(mast_col + 3 * nebula_ruv_col))) + 
       geom_point() + 
       geom_text_repel(aes(label=gene)) + 
       geom_hline(yintercept=0) + 
       geom_vline(xintercept=0) + 
       theme_pubr()

table(rlong[,c('mast_col','nebula_ruv_col')])




