#!/usr/bin/R
# -------------------------------------------------------
# Extract a specific set of genes across all DEG results:
# Updated: 05/25/21
# -------------------------------------------------------
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

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_extract_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Genes to pull out:
genelist = c('GPC5', 'CDH20','TRPS1','HS6ST3')
mdf = read.delim(paste0(datadir, 'module_df_Ast_leiden_vvt.tsv'), header=T)
k = mdf$leiden[mdf$gene == 'GPC5']
modlist = mdf$gene[mdf$leiden == k]

# -----------------------------
# Load in the regression files:
# -----------------------------
fnlist = list.files(path=regdir, pattern=paste0('nebula_ruv.*rda'))
alldf = c()
for (fn in fnlist){
    load(paste0(regdir, fn))
    base = sub(".rda", "", sub("nebula_ruv.","", fn))
    subdf = fulldf[fulldf$gene %in% c(genelist, modlist),]
    subdf = merge(subdf, nsig[,1:4, drop=F])
    alldf = rbind(alldf, subdf)
}

# -----------------------------------------
# Print some results for GPC5 specifically:
# -----------------------------------------
# GPC5 always down/not DEG in Astrocytes (AG up...)
gene = 'GPC5'
subdf = alldf[alldf$gene == gene,]
table(subdf[,c('celltype','col')])
subdf = alldf[alldf$gene == gene & alldf$celltype == 'Ast',]
table(subdf[,c('region','col')]) # Always down except AG
# Always down except in diffuse plaque:
table(subdf[,c('path','col')]) 
hist(subdf$logFC[subdf$path == 'nrad'], 20)
subdf = subdf[order(subdf$logFC),]
head(subdf[,c('logFC','q','subtype','region','path')], 40)

# -------------------------------------
# Look at the aggregate of the modlist:
# -------------------------------------
subdf = alldf[alldf$gene %in% modlist & alldf$celltype == 'Ast',]
table(subdf[,c('region','col')]) # Always down except AG, some EC (?)
# Always down except in diffuse plaque:
table(subdf[,c('path','col')]) 
hist(subdf$logFC[subdf$path == 'nrad'], 20)
hist(subdf$logFC[subdf$path == 'cogdxad'], 20)
subdf = subdf[order(subdf$logFC),]
head(subdf[,c('logFC','q','subtype','region','path')], 40)




