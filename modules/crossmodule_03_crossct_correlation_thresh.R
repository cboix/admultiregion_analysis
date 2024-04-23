#!/usr/bin/R
# ----------------------------------------------------------
# Calculate the cross-ct module correlations and thresholds:
# Updated 11/30/2021
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
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

# ID clusters
require(cba)
nclust = 20
dt <- dist(modcor, 'euclidean')
ht <- hclust(dt, method='ward.D2')
cocl <- order.optimal(dt, ht$merge)$order
reord <- names(cocl)[cocl]
modcor = modcor[reord, reord]
acut <- paste0('C', cutree(ht, nclust)[cocl])

# Add an annotation based on the names and plot:
# ----------------------------------------------
ux = .25
ctann = sub("-.*","", rownames(modcor))
acls = paste0('C', 1:nclust)
acls.cols = snap.cols[1:nclust]
names(acls.cols) = acls

major.col2 = major.col
names(major.col2) = sub("/","_", names(major.col))
major.col2['All'] = 'grey80'

ra = HeatmapAnnotation(
    Set=ctann, Cluster=acut,
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(2, 'mm'),
    gap = unit(0, "mm"),
    col=list(Set=major.col2,
        Cluster=acls.cols), 
    which='row')

gap = 0
plt = Heatmap(modcor, 
    name='Correlation',
    use_raster=TRUE,
    col=col_fun,
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    row_gap=unit(gap, "mm"),
    column_gap=unit(gap, "mm"),
    column_split=acut,
    row_split=acut,
    show_column_names=FALSE,
    show_row_names=FALSE,
    right_annotation=ra,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(modcor)*unit(ux, "mm"), 
    height = nrow(modcor)*unit(ux, "mm")
)

h = 2.25 + 1 / 15 * nrow(modcor)
w = 5 + 1 / 15 * ncol(modcor)
pltprefix = paste0(imgpref, 'corr')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Calculate and plot the same matrix as jaccard overlap:
# ------------------------------------------------------
getGeneSet <- function(module, uselist=gmlist){
    rs = sub("-.*", "", module)
    num = as.numeric(sub(".*-", "", module))
    coremap = uselist[[rs]]
    x = names(coremap)[coremap == num]
    return(x)
}

jmat = modcor * 0
for (i in rownames(jmat)){
    s1 = getGeneSet(i)
    for (j in colnames(jmat)){
        s2 = getGeneSet(j)
        jmat[i,j] = length(intersect(s1, s2)) / length(union(s1, s2))
    }
}

mx = .1
jmat[jmat > mx] = mx

plt = Heatmap(jmat, 
    name='Jaccard distance',
    use_raster=TRUE,
    col=viridis(100),
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    row_gap=unit(gap, "mm"),
    column_gap=unit(gap, "mm"),
    column_split=acut,
    row_split=acut,
    show_column_names=FALSE,
    show_row_names=FALSE,
    right_annotation=ra,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(jmat)*unit(ux, "mm"), 
    height = nrow(jmat)*unit(ux, "mm")
)

h = 2.25 + 1 / 15 * nrow(jmat)
w = 5 + 1 / 15 * ncol(jmat)
pltprefix = paste0(imgpref, 'jaccovl')
saveHeatmap(plt, pltprefix, w=w, h=h)


plt = Heatmap(jmat * modcor, 
    name='Jaccard * Corr.',
    use_raster=TRUE,
    col=rev(colspec),
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    row_gap=unit(gap, "mm"),
    column_gap=unit(gap, "mm"),
    column_split=acut,
    row_split=acut,
    show_column_names=FALSE,
    show_row_names=FALSE,
    right_annotation=ra,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(jmat)*unit(ux, "mm"), 
    height = nrow(jmat)*unit(ux, "mm")
)

h = 2.25 + 1 / 15 * nrow(jmat)
w = 5 + 1 / 15 * ncol(jmat)
pltprefix = paste0(imgpref, 'jacc_corr')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Write out these module clusters:
# --------------------------------
# TODO: print a couple...
# clist = c(9, 5, 7, 18, 19, 16, 10, 8)
clist = c(1:nclust)
clsdf = c()
for (cind in clist){
    submod = rownames(modcor)[acut == paste0('C', cind)]
    subgenes = lapply(submod, uselist=cmlist, getGeneSet)
    subgenes = table(unlist(subgenes))
    subgenes = subgenes[subgenes >= 2]
    subgenes = sort(subgenes, decreasing=T)
    sdf = data.frame(subgenes)
    names(sdf) = c('gene','count')
    sdf$cls = cind
    clsdf = rbind(clsdf, sdf)
    cat(cind, '\tmax:\t', max(subgenes), '\n')
    # cat(head(names(subgenes), 20), '\n\n')
    cat(names(subgenes), '\n\n')
}
clsdf$col = acls.cols[paste0('C', clsdf$cls)]
modcls.tsv = paste0(crossdir, 'shared_genes_module_clusters.tsv')
write.table(clsdf, modcls.tsv, quote=F, row.names=F, sep="\t")

clsassign = data.frame(mod=rownames(modcor), cls=acut, col=acls.cols[acut])
cls.tsv = paste0(crossdir, 'module_cluster_assignments.tsv')
write.table(clsassign, cls.tsv, quote=F, row.names=F, sep="\t")

gene = 'HSP90AA1'
gmlist[['All']][gene]
acut[which(rownames(modcor) == 'All-14')]


# Cutoff for a network:
# ---------------------
r = .4
# n = 48 # Effective N may be more between nrow(scmat) and 48.
n = nrow(scmat)
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
if (modlevel == 'subtype'){cutoff = 0.7} else {cutoff = 0.6}
netlist = make.network(ll$edgedf[,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff)
netlist = set.network.params(netlist)
V(netlist$net)$size = 3

pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_correlation_thresh_network')
w = 4
netlist$lb1 = netlist$labels
netlist$labels = sub(".*: ","",netlist$lb1)
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.15, spacing=.25, adjust=FALSE)
dev.off()


# Cross-subtype only:
ind = substr(ll$edgedf$M1,1,3) != substr(ll$edgedf$M2,1,3) 
netlist = make.network(ll$edgedf[ind,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff)
netlist = set.network.params(netlist, seed=0)
V(netlist$net)$size = 3

pltprefix = paste0(imgpref, 'pseudobulk_crossctonly_', 
                   modlevel, '_correlation_thresh_network')
w = 4
netlist$lb1 = netlist$labels
netlist$labels = sub(".*: ","",netlist$lb1)
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.15, spacing=.25, adjust=FALSE)
dev.off()
