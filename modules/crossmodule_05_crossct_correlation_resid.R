#!/usr/bin/R
# ---------------------------------------------------------
# Calculate the cross-ct module correlation from residuals:
# Updated 11/30/2021
# ---------------------------------------------------------
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


# Try regressing signal and then correlation: 
# --------------------------------------------
mdf = data.frame(pr=rownames(scmat),
                 projid = sapply(rownames(scmat), function(x){sub("-.*","",x)}),
                 region = sapply(rownames(scmat), function(x){sub(".*-","",x)}))
mdf = merge(mdf, metadata[,c('projid','region','rind','msex','age_death', 'cogdxad','nrad')])
rownames(mdf) = mdf$pr
mdf = mdf[rownames(scmat),]

mdf$region = as.character(mdf$region)
scres = scmat
for (i in 1:ncol(scmat)){
    mdf$y = scmat[,i]
    ind = !is.na(mdf$y)
    regdf = mdf[ind,]
    if (length(unique(regdf$region)) > 1){
        fit = glm(y ~ region + nrad, regdf, family='gaussian')
    } else {
        fit = glm(y ~ nrad, regdf, family='gaussian')
    } 
    pred = predict(fit, regdf)
    scres[ind,i] = scmat[ind,i] - pred
}

# Calculate and threshold the correlation matrix:
# -----------------------------------------------
modcor = cor(scres, use='pairwise.complete.obs')
modcor[is.na(modcor)] = 0
plt = plotSymMat(modcor, col_fun=col_fun, ut=4)
wdf = edgesFromMat(modcor)

# n = 48 # Effective N may be more between nrow(scres) and 48.
n = nrow(scres)
wdf$p = (1 - wdf$est^2)^(n/2 - 2) / beta(1/2, n/2 - 1)
wdf$p.adj = p.adjust(wdf$p,'fdr')
cutoff = min(wdf$est[wdf$p.adj < 0.01])
plt = plotSymMat(modcor * (modcor >= cutoff), col_fun=col_fun, ut=ut)


# Make edgelist and set list:
# ---------------------------
edgedf = wdf[wdf$est > cutoff,]
ll = prune.edges(edgedf, modlevel)


# Make a network and plot:
# ------------------------
if (modlevel == 'subtype'){cutoff = 0.7} else {cutoff = 0.5}
netlist = make.network(ll$edgedf[,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff)
netlist = set.network.params(netlist)
V(netlist$net)$size = 3
netlist$lb1 = netlist$labels
netlist$labels = sub(".*: ","",netlist$lb1)

pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_correlation_resid_network')
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.15, spacing=.25, adjust=FALSE)
dev.off()


# Cross-subtype only:
ind = substr(ll$edgedf$M1,1,3) != substr(ll$edgedf$M2,1,3) 
netlist = make.network(ll$edgedf[ind,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff)
netlist = set.network.params(netlist, seed=0)
V(netlist$net)$size = 3
netlist$lb1 = netlist$labels
netlist$labels = sub(".*: ","",netlist$lb1)

pltprefix = paste0(imgpref, 'pseudobulk_crossctonly_', 
                   modlevel, '_correlation_resid_network')
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.15, spacing=.25, adjust=FALSE)
dev.off()

