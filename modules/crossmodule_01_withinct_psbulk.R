#!/usr/bin/R
# ------------------------------------------------------------------------
# Plot cross-module scores within the same celltype (at pseudobulk level):
# Updated 02/16/2021
# ------------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(ggrastr)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


imgpref = paste0(plotdir, 'cross1ct_', runset, "_")

# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))

runlist = c('Ast', 'Mic_Immune', 'Vasc_Epithelia', 'Oli','Opc',
    'Inh', 'Exc', 'All')
    # 'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons', 'Glial', 
for (runset in runlist){
imgpref = paste0(plotdir, 'cross1ct_', runset, "_")

# Subset and turn the score dataframe into a matrix:
# --------------------------------------------------
mingenes = 10
subscdf = runscdf[runscdf$ng >= mingenes,]
subscdf = subscdf[subscdf$runset %in% runset,]
scdf = scoredf[scoredf$runset %in% runset & scoredf$ng >= mingenes,]

mmap = unique(scdf[,c('module','mname')])
rownames(mmap) = as.character(mmap$module)
mmap$mname2 = paste0(runset, '-', mmap$module)
mn.mapping = mmap$mname
names(mn.mapping) = mmap$mname2

mat = pivot.tomatrix(subscdf[, c('pr','rm','score')], 'pr','score')
scmat = log(t(mat))
ut = 3


# Calculate and threshold the correlation matrix:
# -----------------------------------------------
modcor = cor(scmat, use='pairwise.complete.obs')
modcor[is.na(modcor)] = 0
modcor2 = modcor
rownames(modcor2) = mn.mapping[rownames(modcor2)]
colnames(modcor2) = mn.mapping[colnames(modcor2)]
plt = plotSymMat(modcor2, col_fun=col_fun, ut=1.5, raster=FALSE)
wdf = edgesFromMat(modcor)

pltprefix = paste0(imgpref, 'corr')
saveHeatmap(plt, pltprefix, w=4, h=4)

}

# Plot in specific order for microglia:
if (runset == 'Mic_Immune'){
mic.mod = c(25, 11, 2, 23, 1, 21, 19, 
    6, 17, 9, 14, 18, 10, 12, 13, 20, 3)
} else {
    mic.mod = sort(unique(scdf$module))
}
mnord = mmap[as.character(mic.mod),'mname']
mnord = paste0(runset, '-', mic.mod)

plt = plotSymMat(modcor[mnord, mnord], col_fun=col_fun, 
    ut=1.5, cluster=FALSE, raster=FALSE)
pltprefix = paste0(imgpref, 'corr_fixed')
saveHeatmap(plt, pltprefix, w=4, h=4)


# # More permissive cutoff:
# NM = ncol(scmat)
# edf = c()
# for (i in 1:NM){
#     for (j in 1:NM){
#         m1 = colnames(scmat)[i]
#         m2 = colnames(scmat)[j]
#         ct = cor.test(scmat[,i], scmat[,j], alternative='greater')
#         edf = rbind(edf, 
#             data.frame(M1=m1, M2=m2, 
#                 est=ct$estimate, p=ct$p.value))
#     }
# }
# edf$p.adj = p.adjust(edf$p, 'fdr')
# min(abs(edf$est[edf$p.adj < 0.01]))


# More stringent cutoff:
r = .4
n = nrow(scmat)
wdf$p = (1 - wdf$est^2)^(n/2 - 2) / beta(1/2, n/2 - 1)
wdf$p.adj = p.adjust(wdf$p,'fdr')
cutoff = min(wdf$est[wdf$p.adj < 0.01])
plt = plotSymMat(modcor * (modcor >= cutoff), col_fun=col_fun, ut=1.5)
pltprefix = paste0(imgpref, 'corrkept')
saveHeatmap(plt, pltprefix, w=4, h=4)


# Make edgelist and set list:
# ---------------------------
modlevel = 'runset'
edgedf = wdf[wdf$est > cutoff,]
nodes = rownames(modcor)

# Remove specific ambient modules:
if (runset == 'Mic_Immune'){
    rm.mod = c(12, 13, 16)
} else if (runset == 'Ast'){
    rm.mod = c(4)
} else {
    rm.mod = c()
}
rm.mod = paste0(runset, '-', rm.mod)
nodes = nodes[!(nodes %in% rm.mod)]

edgedf = edgedf[edgedf$M1 %in% nodes,]
edgedf = edgedf[edgedf$M2 %in% nodes,]
kept.modules = sub(".*-","", nodes)


# Make edges and pie-values:
# --------------------------
ll = prune.edges(edgedf, modlevel, nodes=nodes)

# Make module by subtype score matrix:
scdf$scsum = scdf$score * scdf$ncell
scaggdf = aggregate(cbind(scsum, ncell) ~ module + cell_type_high_resolution, 
                  scdf[scdf$module %in% kept.modules,], sum)
scaggdf$score = with(scaggdf, scsum / ncell)
pvmat = pivot.tomatrix(scaggdf[,c('cell_type_high_resolution','module','score')], 'module','score')
colnames(pvmat) = paste0(runset, '-', colnames(pvmat))

# Make pie chart values:
# pie.values = sweep(pvmat, 2, apply(pvmat, 2, min),'-')
pie.values = pvmat **2
pie.values = sweep(pie.values, 2, colSums(pie.values), '/')
pie.cols = tcols[rownames(pie.values)]
pie.values = lapply(nodes, function(x){pie.values[, x]})
names(pie.values) = nodes


# Make a network and plot:
# ------------------------
w = 1.5
netlist = make.network(ll$edgedf[,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff, label.edge=TRUE)

netlist = set.network.params(netlist, use.lty=FALSE)
V(netlist$net)$size = 15
E(netlist$net)$label.cex = .15


pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_correlation_thresh_network_pie')
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, pie.values=pie.values, pie.cols=pie.cols,
             lbcex=0.15, spacing=.25, adjust=FALSE)
dev.off()
pdf(paste0(pltprefix, '.pdf'), width=w, height=w)
plot.network(netlist, pie.values=pie.values, pie.cols=pie.cols,
             lbcex=0.15, spacing=.25, adjust=FALSE)
dev.off()




# Cross-subtype only:
# -------------------
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



# Look at within a single sub-type:
# ---------------------------------
# Calculate the within-subtype, cross-module correlation:
selfcor.list = list()
modcor = NULL
if (runset == 'Mic_Immune'){
    run.subtypes = subtypes[grep("^Mic",subtypes)]
} else {run.subtypes = subtypes }

for (subtype in run.subtypes){
    subdf = scdf[scdf$cell_type_high_resolution == subtype,]
    smat = pivot.tomatrix(subdf[,c('pr','score','module')], 'module', 'score')
    selfcor.list[[subtype]] = cor(smat)
    if (is.null(modcor)){
        modcor = selfcor.list[[subtype]] 
    } else {
        modcor = modcor + selfcor.list[[subtype]] 
    }
}
modcor = modcor / length(run.subtypes)

rownames(modcor) = mmap$mname[as.numeric(rownames(modcor)) + 1]
plt = Heatmap(modcor, 
              col=col_fun,
              use_raster=TRUE,
              width = ncol(modcor)*unit(5, "mm"), 
              height = nrow(modcor)*unit(5, "mm"),
              border_gp = gpar(col="black", lty = 1))

h = 2.25 + 2.5 / 15 * nrow(modcor)
w = 5 + 2.5 / 15 * ncol(modcor)
pltprefix = paste0(imgpref, 'psbulk_module_correlation_', fullpref)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()


# TODO: How do we pick the top module pairs/how do we decide what to follow up on?
#           --> thresholding correlation
#           --> significance from regressions
#           --> maybe cross-celltype is cleaner
# TODO: Should we plot a within-celltype network of modules (from cell-level?)
# TODO: How do we account for subtype-specific vs. stable modules
# TODO: Should we run this separately for all genes and for up/down DEGs - module signal may clean up, but unclear how to do this systematically
# TODO: Can we figure out the directionality of the correlations 
#           --> anchor some points on the graph (aka least AD, most AD)
#           --> perform some sort of if -> vs. <- (causal inf)
#           --> use three nodes to determine -- use the cell specific signal as the confounding?
#           --> determine ligand-receptor pairs and learn directionality from this


# Look at across subtypes (within same ct):
# -----------------------------------------
# Calculate the cross-subtype, cross-module correlation:
crosscor.list = list()
modcor = NULL
scols = c('pr','score','module')
for (s1 in subtypes){
    s1df = scdf[scdf$cell_type_high_resolution == s1, scols]
    s1mat = pivot.tomatrix(s1df, 'module', 'score')

    for (s2 in subtypes[subtypes != s1]){
        s2df = scdf[scdf$cell_type_high_resolution == s2, scols]
        s2mat = pivot.tomatrix(s2df, 'module', 'score')
        keep.rn = intersect(rownames(s1mat), rownames(s2mat))
        modcor = cor(s1mat[keep.rn,], s2mat[keep.rn,])

        rownames(modcor) = mmap$mname[as.numeric(rownames(modcor)) + 1]
        colnames(modcor) = mmap$mname[as.numeric(colnames(modcor)) + 1]


        # # Same rows:
        # subdf = scdf[scdf$cell_type_high_resolution %in% c(s1, s2),]
        # smat = pivot.tomatrix(subdf[,c('pr','score','module')], 'module', 'score')
        # selfcor.list[[subtype]] = cor(smat)
        # if (is.null(modcor)){
        #     crosscor = selfcor.list[[subtype]] 
        # } else {
        #     crosscor = modcor + selfcor.list[[subtype]] 
        # }
    }
}


modcor = modcor / length(subtypes)

rownames(modcor) = mmap$mname[as.numeric(rownames(modcor)) + 1]
plt = Heatmap(modcor, 
              col=col_fun,
              use_raster=TRUE,
              width = ncol(modcor)*unit(5, "mm"), 
              height = nrow(modcor)*unit(5, "mm"),
              border_gp = gpar(col="black", lty = 1))



# NOTE: UNFINISHED. POSSIBLY BUILD NETWORKS FROM THIS;





