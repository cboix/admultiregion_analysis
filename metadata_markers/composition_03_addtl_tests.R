#!/usr/bin/R
# ------------------------------------------------------------
# Calculate the fraction differences, by dirichlet regression:
# Updated 04/13/2021 
# Last updated 11/25/2021 to clean up code only
# ------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(qvalue)
library(lme4)
library(emmeans)
library(nnet) # multinom
library(mlogit) # mlogit
library(DirichletReg)


# Directories:
fracdir = paste0(sdbdir, 'fractions/')
plotdir = paste0(imgdir, 'fractions/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, fracdir)
system(cmd)


# Function for plotting odds-ratio matrix:
# ----------------------------------------
plotCoeffPvaluesMatrix = function(cmat, pmat, pltprefix){
    require(ComplexHeatmap)
    plt = Heatmap(cmat,
                  use_raster=TRUE,
                  cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                      p = pmat[i,j]
                      ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                      grid.text(ann, x, y)}
    )
    h = 2.25 + 2.5 / 15 * nrow(cmat)
    w = 4 + 2.5 / 15 * ncol(cmat)
    pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
    print(plt)
    dev.off()
    png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
    print(plt)
    dev.off()
}


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
subsetlist = c('All_minor', 'All', 'OpcOli', unique(cellmeta$major.celltype))
remove.batches = TRUE
suff = '_subset_final_noMB'

metadata$braaksc.ad = 'CTRL'
metadata$braaksc.ad[metadata$braaksc %in% c(5,6)] = 'AD'
metadata$braaksc.ad = factor(metadata$braaksc.ad, levels=c('CTRL','AD'))

covarcols = c('msex', 'age_death', 'pmi')
advars = c('niareagansc','nrad','cogdxad','braaksc.ad', 'ceradsc','braaksc','cogdx')
pathvars = c('nft','plaq_n','plaq_d')


fulldf = c()
subset = 'Ast'

# for (subset in subsetlist){
print(subset)
ststr = gsub("/","_", subset)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))


# Make count matrix for Dirichlet regression:
# -------------------------------------------
cmat = pivot.tomatrix(ctdf[,c('cls','batch','count')], 'cls', 'count')
nmat = cmat / rowSums(cmat)

# Add covariates: 
mdf = data.frame(projid=sub("_.*","", rownames(cmat)), 
                 region=sub(".*_","", rownames(cmat)), 
                 batch=rownames(cmat))
mdf = merge(mdf, metadata[,c('projid','region','rind', 
                             covarcols, advars)], all.x=TRUE)
mdf = merge(mdf, pqdf, all.x=TRUE)
rownames(mdf) = mdf$batch
# Log 1p the pathology measurements:
for (path in pathvars){ mdf[[path]] = log1p(mdf[[path]]) }
mdf = mdf[rownames(nmat),]


# Perform the Dirichlet multinomial regression:
# ---------------------------------------------
counts = as.data.frame(nmat)
counts$counts = DR_data(counts)
data = cbind(counts, mdf)

subsetdf = c()
for (advar in c(advars, pathvars)){
    print(advar)
    drForm = asform(c('counts ~ region +', advar, '+', paste(covarcols, collapse='+')))
    fit = DirichReg(formula=as.formula(drForm), data=data) # , model='alternative') (reference off one)
    # Store results as dataframe:
    cfit = summary(fit)$coef.mat
    vn = rownames(cfit)
    cfit = data.frame(cfit)
    colnames(cfit) = c('est','sd','z','p')
    cfit$variable = vn
    # Add the subtypes:
    nvar = length(u$varnames)
    cfit$cls = rep(u$varnames, rep(nrow(cfit)/nvar, nvar))
    cfit$path = advar
    subsetdf = rbind(subsetdf, cfit)
}

# Subset to just the variables:
kind = ((subsetdf$variable == subsetdf$path) | 
        (subsetdf$variable == paste0(subsetdf$path,'AD')))

pltdf = subsetdf[kind,]
pltdf = pltdf[order(pltdf$p),]
pltdf$p.adj = p.adjust(pltdf$p, 'fdr')
cmat = pivot.tomatrix(pltdf[,c('cls','est','path')], 'path','est')
pmat = pivot.tomatrix(pltdf[,c('cls','p.adj','path')], 'path','p.adj')
pmat = pivot.tomatrix(pltdf[,c('cls','p','path')], 'path','p')

heatpref = paste0(imgpref, 'dirichletreg_', ststr, '_heatmap')
plotCoeffPvaluesMatrix(cmat, pmat, pltprefix=heatpref) 

# NOTE: The p-values and estimates look very odd given effects we see happening othw.


# Fisher's Exact test:
# --------------------
# Run individual-agnostic test:
cf = unique(metadata[,c('projid','region','rind','braaksc','cogdx','niareagansc','msex','age_death','pmi')])
ctdf = merge(ctdf, cf)
ctdf$nrad = 'AD'
ctdf$nrad[ctdf$niareagansc > 2] = 'Control'
ctdf$nrad = factor(ctdf$nrad, c('Control','AD'))

aggdf = aggregate(count ~ cls + major.celltype + region + nrad, ctdf, sum)
c1 = spread(aggdf, nrad, count)
aggdf = aggregate(total ~ cls + major.celltype + region + nrad, ctdf, sum)
t1 = spread(aggdf, nrad, total)
names(t1)[4:5] = paste0(names(t1)[4:5], '.tot')
aggdf = merge(c1, t1)

# Perform test (everything super significant, no accounting for individuals):
aggdf$ft.p = sapply(1:nrow(aggdf), function(i){
                        mat = matrix(unlist(aggdf[i,4:7]),2,2, byrow=T)
                        mat[2,] = mat[2,] - mat[1,]
                        ft = fisher.test(mat)
                        ft$p } )
aggdf$ft.padj = p.adjust(aggdf$ft.p)
aggdf = aggdf[order(aggdf$ft.padj),]


# Mann-Whitney test (non-parametric):
# can account for level of individuals - but doesn't capture multi-region
# -----------------------------------
ctdf$frac = ctdf$count / ctdf$total

# ggplot(ctdf, aes(cls, frac, color=factor(msex))) + 
ggplot(ctdf, aes(cls, frac, color=nrad)) + 
    facet_wrap(~region) + 
    geom_boxplot() + 
    theme_pubr() + 
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust=.5)) + 
    stat_compare_means(method='wilcox.test', label='p.signif')







