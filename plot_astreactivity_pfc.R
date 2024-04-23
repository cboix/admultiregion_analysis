#!/usr/bin/R
# ------------------------------------------------
# Explore astrocyte changes in larger data cohort:
# Updated 10/07/2021 
# ------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Matrix)

# library(ggrepel)
# library(ggpmisc)
# library(patchwork)
# library(ComplexHeatmap)
# library(circlize)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/ast_reac/')
imgpref = paste0(plotdir, 'astreac_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# --------------------------
# Load astrocyte snPFC data:
# --------------------------
datadir = 'multiRegion/ast_pfc/'
anndir = 'Annotation/'

# Read in data:
metadata = read.delim(paste0(anndir, 'metadata_PFC_all_individuals_092520.tsv'), header=T)
cellmeta = read.delim(paste0(datadir, 'snPFC.Astro.metadata.txt'), header=T)
cellmeta = cellmeta[,c('projid','batch','barcode','nCount_RNA','nFeature_RNA','percent.mt','decontX_contamination')]

# Matrix ~ 126k x 16k
mat = readRDS(paste0(datadir, 'snPFC.Astro.counts.rds'))

# TODO: Write data out as hdf5 or mtx for reading with python

# Annotate cells with metadata:
# -----------------------------
projids = unique(as.character(cellmeta$projid))
rownames(metadata) = as.character(metadata$projid)
metadata$nrad = 'AD'
metadata$nrad[metadata$niareagansc > 2] = 'CTRL'
metadata$cogdxad = 'AD'
metadata$cogdxad[metadata$cogdx < 4] = 'CTRL'
metadata = metadata[metadata$projid %in% projids,]
metadata$delta = 'N-N'
metadata$delta[metadata$cogdxad == 'AD' & metadata$nrad == 'CTRL'] = 'N-CI'
metadata$delta[metadata$cogdxad == 'CTRL' & metadata$nrad == 'AD'] = 'AD-N'
metadata$delta[metadata$cogdxad == 'AD' & metadata$nrad == 'AD'] = 'AD-CI'
table(metadata$delta)

vars = c('cogdx','nrad','cogdxad','niareagansc', 'braaksc', 'delta')
cellmeta = cbind(cellmeta, metadata[as.character(cellmeta$projid),vars, drop=F])


# Score cells by SLC1A2/1A3 and by reactivity:
# --------------------------------------------
reac.genes = c('GFAP','OSMR','CD44')
glu.genes = c('SLC1A2','SLC1A3')
joint.genes = c('GFAP','OSMR','CD44','SLC1A3')

gfap.score = mat['GFAP',] > 0
gg.score = colSums(mat[c('GFAP', 'SLC1A3'),] > 0) / 2
reac.score = colSums(mat[reac.genes,] > 0) / length(reac.genes)
glu.score = colSums(mat[glu.genes,] > 0) / length(glu.genes)
joint.score = colSums(mat[joint.genes,] > 0) / length(joint.genes)

cellmeta$gfap.score = 1 * (gfap.score[rownames(cellmeta)] == 1)
cellmeta$gg.score = 1 * (gg.score[rownames(cellmeta)] == 1)
cellmeta$reac.score = 1 * (reac.score[rownames(cellmeta)] == 1)
cellmeta$glu.score = 1 * (glu.score[rownames(cellmeta)] == 1)
cellmeta$joint.score = 1 * (joint.score[rownames(cellmeta)] == 1)

aggdf = aggregate(cbind(reac.score, glu.score, joint.score, gfap.score, gg.score) ~ projid + delta, cellmeta, mean)
aggregate(cbind(reac.score, glu.score, joint.score, gfap.score, gg.score) ~ delta, aggdf, mean)

gp = ggplot(aggdf, aes(delta, gfap.score)) + 
    geom_violin(fill='grey90') + 
    geom_boxplot(width=.25, outlier.shape=NA) + 
    geom_jitter(width=.15, cex=.25) + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(title='Only GFAP', x='Pathology + Cognition', y='% of cells with GFAP') + 
    theme_pubr()
ggsave('~/test_gfap.png', dpi=400, gp, units='in', width=5, height=5)

gp = ggplot(aggdf, aes(delta, reac.score)) + 
    geom_violin(fill='grey90') + 
    geom_boxplot(width=.25, outlier.shape=NA) + 
    geom_jitter(width=.15, cex=.25) + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(title='Reactive: GFAP/CD44/OSMR', x='Pathology + Cognition', y='% of cells with all genes') + 
    theme_pubr()
ggsave('~/test_reac.png', dpi=400, gp, units='in', width=5, height=5)

gp = ggplot(aggdf, aes(delta, joint.score)) + 
    geom_violin(fill='grey90') + 
    geom_boxplot(width=.25, outlier.shape=NA) + 
    geom_jitter(width=.15, cex=.25) + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(title='Joint: GFAP/CD44/OSMR + SLC1A3', x='Pathology + Cognition', y='% of cells with all genes') + 
    theme_pubr()
ggsave('~/test_joint.png', dpi=400, gp, units='in', width=5, height=5)

gp = ggplot(aggdf, aes(delta, gg.score)) + 
    geom_violin(fill='grey90') + 
    geom_boxplot(width=.25, outlier.shape=NA) + 
    geom_jitter(width=.15, cex=.25) + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(title='Two markers: GFAP + SLC1A3', x='Pathology + Cognition', y='% of cells with both genes') + 
    theme_pubr()
ggsave('~/test_gfap_glt1.png', dpi=400, gp, units='in', width=5, height=5)



