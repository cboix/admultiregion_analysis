#!/usr/bin/R
# ---------------------------------------
# Quantify the cross-module correlations:
# Updated 04/02/2022
# ---------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(ggpubr)
library(gprofiler2)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Turn the score dataframe into a matrix:
# NOTE: At the runset level or at the subtype level:
# --------------------------------------------------
mingenes = 10
modlevel = 'subtype'
modlevel = 'runset'
keep.sets = c('All','Exc','Inh','Ast', 'Mic_Immune','Oli','Opc','Vasc_Epithelia')
scoredf = scoredf[scoredf$runset %in% keep.sets,]
runscdf = runscdf[runscdf$runset %in% keep.sets,]
if (modlevel == 'subtype'){
    mat = pivot.tomatrix(scoredf[scoredf$ng >= mingenes, c('pr','cm','score')], 'pr','score')
    scmat = log(t(mat))
    ut = 1
} else {
    mat = pivot.tomatrix(runscdf[runscdf$ng >= mingenes, c('pr','rm','score')], 'pr','score')
    scmat = log(t(mat))
    ut = 3
}

imgpref = paste0(plotdir, 'cross_allct_')

# Calculate and threshold the correlation matrix:
# -----------------------------------------------
modcor = cor(scmat, use='pairwise.complete.obs')
modcor[is.na(modcor)] = 0
wdf = edgesFromMat(modcor)


# Cutoff for a network:
# ---------------------
r = .4
n = nrow(scmat)
wdf$p = (1 - wdf$est^2)^(n/2 - 2) / beta(1/2, n/2 - 1)
wdf$p.adj = p.adjust(wdf$p,'fdr')
cutoff = min(wdf$est[wdf$p.adj < 0.01])

wdf$S1 = sub("-.*","", wdf$M1)
wdf$S2 = sub("-.*","", wdf$M2)
wdf$sig = wdf$p.adj < 0.01
wdf$sig = wdf$est > .5


# Make a heatmap of percent of significant interactions:
# ------------------------------------------------------
swide = aggregate(sig ~ S1 + S2, wdf, mean)
smat = pivot.tomatrix(swide, 'S2', 'sig')
# smat = smat - diag(diag(smat))
cap = 0.4
smat[smat > cap] = cap

plt = plotSymMat(smat, 
    col_fun=colb, ut=1.5, raster=FALSE)
h = 2 + 1 / 15 * nrow(smat)
w = 3 + 1 / 15 * ncol(smat)
pltprefix = paste0(imgpref, 'percentcorr')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Print the top OPC-Oli pairs and their genes:
# --------------------------------------------
intdf = head(wdf[wdf$S1 == 'Opc' & wdf$S2 == 'Oli',], 20)
intdf$N1 = rmap[intdf$M1, 'mname']
intdf$N2 = rmap[intdf$M2, 'mname']

getGeneSet <- function(module, uselist=gmlist){
    rs = sub("-.*", "", module)
    num = as.numeric(sub(".*-", "", module))
    coremap = uselist[[rs]]
    x = names(coremap)[coremap == num]
    return(x)
}

jaccOvl = function(i, j, uselist=gmlist){
    s1 = getGeneSet(i, uselist=uselist)
    s2 = getGeneSet(j, uselist=uselist)
    jc = length(intersect(s1, s2)) / length(union(s1, s2))
    return(jc)
}

intdf$jacc = sapply(1:nrow(intdf), function(x){
    jaccOvl(intdf$M1[x], intdf$M2[x]) })

# Also get cluster:
cls.tsv = paste0(crossdir, 'module_cluster_assignments.tsv')
clsdf = read.delim(cls.tsv, header=T)
rownames(clsdf) = clsdf$mod

intdf$C1 = clsdf[intdf$M1, 'cls']
intdf$C2 = clsdf[intdf$M2, 'cls']

# Print resulting data.frame:
print(intdf[,c('est','jacc','C1', 'C2','N1','N2')])
