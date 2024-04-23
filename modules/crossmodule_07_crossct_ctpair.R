#!/usr/bin/R
# -------------------------------------------------------------
# Calculate the cross-ct module corr. for a pair of cell types:
# NOTE: Mainly for OPC / Oli
# Updated 03/26/2022
# -------------------------------------------------------------
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


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
# TODO: Opt to include neurons or not? Or just add later...
# commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Subset and turn the score dataframe into a matrix:
# --------------------------------------------------
mingenes = 10
runsets = c('Opc','Oli')
rstr = paste(runsets, collapse='_')
imgpref = paste0(plotdir, 'cross2ct_', rstr)

subscdf = runscdf[runscdf$ng >= mingenes,]
subscdf = subscdf[subscdf$runset %in% runsets,]
scdf = scoredf[scoredf$runset %in% runsets & scoredf$ng >= mingenes,]
scdf$rm = with(scdf, paste0(runset, '-', module))

mat = pivot.tomatrix(subscdf[, c('pr','rm','score')], 'pr','score')
scmat = log(t(mat))
ut = 3


# Calculate and threshold the correlation matrix:
# -----------------------------------------------
modcor = cor(scmat, use='pairwise.complete.obs')
modcor[is.na(modcor)] = 0
plt = plotSymMat(modcor, col_fun=col_fun, ut=1.5)
wdf = edgesFromMat(modcor)

pltprefix = paste0(imgpref, 'corr')
saveHeatmap(plt, pltprefix, w=5, h=5)


# Plot selected pairs in this:
pltpairs = list(c('Opc-25','Oli-23'), c('Opc-4','Oli-4'), c('Opc-4', 'Oli-1'),
                c('Opc-1','Oli-3'), c('Opc-5','Oli-1'), c('Opc-12','Oli-1'), c('Opc-9', 'Oli-7'), c('Opc-7','Oli-7'))

var = 'cogdxad'
for (pair in pltpairs){
    pstr = paste(pair, collapse='_')
    subdf = data.frame(log1p(t(mat)[,pair]))
    pair = names(subdf)
    subdf$ptype = rownames(scmat)
    subdf$projid = sub("-.*","", subdf$ptype)
    subdf$region = sub(".*-","", subdf$ptype)
    if (!(var %in% colnames(subdf))){
        subdf = merge(subdf, unique(metadata[,c('projid','region', var)]))
    }

    gp = ggplot(subdf, aes_string(pair[1], pair[2], color=var)) + 
        geom_point(cex=.25) +
        geom_smooth(method='lm', color='black',lwd=.5) + 
        stat_cor(color='black', cex=3, output.type='text', label.sep='\n', label.y.npc=.95) + 
        theme_pubr() + theme(legend.position='none')
    if (var == 'region'){
        gp = gp + scale_color_manual(values=reg.cols)
    } else {
        gp = gp + scale_color_manual(values=c('CTRL'='grey80','AD'='red'))
    }
    gp2 = rasterize(gp, layers='Point', dpi=450)

    pltprefix = paste0(imgpref, 'corrplot_', var, '_',pstr)
    saveGGplot(gp2, pltprefix, w=2, h=2)
}



# Plot the Opc-Oli network:
# -------------------------
n = nrow(scmat)
wdf$p = (1 - wdf$est^2)^(n/2 - 2) / beta(1/2, n/2 - 1)
wdf$p.adj = p.adjust(wdf$p,'fdr')
cutoff = min(wdf$est[wdf$p.adj < 0.01])
wdf$samect = sub("-.*", "", wdf$M1) == sub("-.*", "", wdf$M2)
plt = plotSymMat(modcor * (modcor >= cutoff), col_fun=col_fun, ut=1.5)
pltprefix = paste0(imgpref, 'corrkept')
saveHeatmap(plt, pltprefix, w=5, h=5)

# Filtered by cell type:
wdf$samect = sub("-.*", "", wdf$M1) == sub("-.*", "", wdf$M2)
wdf$est.filt = wdf$est * (1 - wdf$samect) * (wdf$p.adj < 0.001)
wmat = pivot.tomatrix(wdf[,c('M1','M2','est.filt')], 'M1', 'est.filt')
wmat[is.na(wmat)] = 0
plt = plotSymMat(wmat, col_fun=col_fun, ut=1.5)
pltprefix = paste0(imgpref, 'corrkept_diffct')
saveHeatmap(plt, pltprefix, w=5, h=5)



# Make edgelist and set list:
# ---------------------------
library(plyr)
modlevel = 'runset'
edgedf = wdf[wdf$est.filt > 0,]

# Top corr only:
ntop = 2
edgedf = ldply(unique(edgedf$M1), function(x){
    head(edgedf[edgedf$M1 == x,], ntop)
                })


# Remove specific ambient modules:
nodes = rownames(modcor)
# nodes = unique(c(edgedf$M1, edgedf$M2))
rm.mod = c('Opc-14', 'Oli-5')  # Mic/Immune
nodes = nodes[!(nodes %in% rm.mod)]
nodes = nodes[nodes %in% unique(c(edgedf$M1, edgedf$M2))]
edgedf = edgedf[edgedf$M1 %in% nodes,]
edgedf = edgedf[edgedf$M2 %in% nodes,]


# Make edges and pie-values:
# -------------------------
ll = prune.edges(edgedf, modlevel, nodes=nodes)

# Make module by subtype score matrix:
scdf$scsum = scdf$score * scdf$ncell
scaggdf = aggregate(cbind(scsum, ncell) ~ rm + cell_type_high_resolution, 
    scdf[scdf$rm %in% nodes,], sum)
scaggdf$score = with(scaggdf, scsum / ncell)
pvmat = pivot.tomatrix(scaggdf[,
    c('cell_type_high_resolution','rm','score')], 'rm','score')
pvmat[is.na(pvmat)] = 0

# Make pie chart values:
# pie.values = sweep(pvmat, 2, apply(pvmat, 2, min),'-')
pie.values = pvmat ** 2
pie.values = sweep(pie.values, 2, colSums(pie.values), '/')
pie.cols = tcols[rownames(pie.values)]
pie.values = lapply(nodes, function(x){pie.values[, x]})
names(pie.values) = nodes


# Make a network and plot:
# ------------------------
w = 1.5
netlist = make.network(ll$edgedf[,c('i','j','est')], 'est', 
    ll$mndf, cutoff=cutoff, label.edge=TRUE, symmetric=FALSE)

netlist = set.network.params(netlist, use.lty=FALSE)
V(netlist$net)$size = 15
E(netlist$net)$label = ''

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

