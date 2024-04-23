#!/usr/bin/R
# -----------------------------------------------------------
# Plot the cellular subtype - motif network for visualization:
# Updated: 03/11/2021
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(igraph)
library(tidyr)
library(Matrix)

library(viridis)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/atlas/')
imgpref = paste0(plotdir, 'network_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# --------------------------
# Load in the main datasets:
# --------------------------
# Metadata:
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

# Colors for full:
typelvls = unique(cellmeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
type.cols = c(type.cols, major.col['Inh'], major.col['Exc'])
tsp.type.cols = sapply(type.cols, tsp.col)
stype.cols = type.cols
names(stype.cols) = sub("^Inh ","", names(stype.cols))

# For name mapping:
inh.sts = unique(cellmeta$cell_type_high_resolution[cellmeta$major.celltype=='Inh'])
imap = data.frame(cell_type_high_resolution=inh.sts, C1=gsub("[ ()-]",".",inh.sts), ct.short=sub("^Inh ","", inh.sts))

# Calculate proportions of each celltype:
ctdf = agg.rename(barcode ~ cell_type_high_resolution + region, cellmeta[cellmeta$major.celltype == 'Inh',], length, 'count')
ctwide = spread(ctdf, region, count, fill=0)
ctwide$total = apply(ctwide[,regions[regions !='MB']], 1, sum)
ctfrac = as.matrix(ctwide[,regions[regions!='MB']] / ctwide$total)
rownames(ctfrac) = sub("^Inh ","",ctwide$cell_type_high_resolution)

# -----------------------------
# Read in the adjacency matrix:
# -----------------------------
tab = read.delim('multiRegion/Inh_1000_mean_mscore_SCENIC_projid_average.csv', sep=",")
names(tab)[1] = 'R1'
df = gather(tab, C1, sim, -R1)
df$R1 = sub("_\\(\\+\\)", "", df$R1)
df = merge(df, imap)
all.inh = unique(df$ct.short)
all.mot = unique(df$R1)

# ---------------------------------------------------------
# Plot a graph/extract clusters/modules (community detect):
# ---------------------------------------------------------
ddf = df[,c('R1','ct.short','sim')] # Object for plotting
cutoff = .25
npref = paste0('cutoff', cutoff)
use.nglk = FALSE
if (use.nglk){
    ddf = ddf[abs(ddf$sim) >= cutoff,]
    ddf$COLOR = ifelse(ddf$sim < 0, 'royalblue','indianred')
    npref = paste0(npref, '_negtvlinks')
} else {
    ddf = ddf[ddf$sim >= cutoff,]
    ddf$COLOR = 'grey25'
}

# Remove links if diff sets:
all.nodes=FALSE
if (all.nodes){
    nodes = c(all.inh, all.mot)
    npref = paste0(npref, '_allnodes')
} else {
    nodes = sort(unique(c(ddf$R1, ddf$ct.short)))
}
# Remove edges in opposite direction
dim(ddf)[1]
dim(ddf)[1] / (dim(df)[1])

# Simple network: just the links/points:
sdf = ddf
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
nname = V(net)$name
ncol = ifelse(nname %in% names(stype.cols), stype.cols[nname], 'grey75')
ntype = 1 *(nname %in% all.inh) + 1

# Vertex arguments:
V(net)$color <- ncol # c("steel blue", "orange")[ntype]
V(net)$shape <- c("square", "circle")[ntype]
# V(net)$pch = c(19,15)[ntype]
V(net)$size = c(3,8)[ntype]
V(net)$label = nodes
V(net)$label.cex = .5
V(net)$label.color = 'black'
V(net)$frame.color <- c(NA,'black')[ntype] # vcol
V(net)$label.font = 1
# Edge arguments:
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
E(net)$color = ecol 
elty = rep('dotted', length(sdf$sim))
elty[abs(sdf$sim) >= .5] = 'dashed'
elty[abs(sdf$sim) >= 1] = 'solid'
E(net)$lty = elty
E(net)$width = abs(sdf$sim) * .5 + .5
E(net)$weight = abs(sdf$sim) * 4

# plot(net)
set.seed(2)
l <- layout_with_fr(net, grid='nogrid') # Usually best

w = 5
png(paste0(imgpref, npref, '_inh_motif_regulon_network_basic.png'), res=300, units='in', width=w,height=w)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
dev.off()


# -----------------------
# Pie chart for vertices:
# -----------------------
pie.values = lapply(nodes, function(x){ 
                        if (x %in% rownames(ctfrac)){
                            ctfrac[x,] 
                        } else { 
                            c(1,rep(0,5)) } })

reg.cols = reg.cols[regions[regions !='MB']]
V(net)$shape = c('square',"pie")[ntype]
V(net)$pie = pie.values
V(net)$pie.color = list(reg.cols)
V(net)$pie.border ='black'
V(net)$pie.lty = 1

png(paste0(imgpref, npref, '_inh_motif_regulon_network_pie.png'), res=300, units='in', width=w,height=w)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# legend('bottomright', 
#        legend=c('max -log10p >= 3', 'max -log10p >= 10', 'max -log10p >= 25', 'max -log10p >= 50', 
#                 'Cosine Sim. >= 0.75', 'Cosine Sim. >= 0.85', 'Cosine Sim. >= 0.95'),
#        col='black', lty=c(NA, NA, NA, NA, 'dotted','dashed','solid'),
#        cex=.4, pch=c(19, 19,19,19, NA, NA, NA), pt.cex=c(1, 1.5, 2, 2.5, .75, .75, .75) / 4, 
#        inset=0, box.col=NA)
legend('bottomleft', legend=names(reg.cols), col=reg.cols, ncol=5,
       cex=.4, pch=15, pt.cex=1, inset=0, box.col=NA)
dev.off()




# ------------------------------------------
# Plot the same graph again, but with repel:
# ------------------------------------------
source(paste0('~/ENCODE_DATA/bin/', 'auxiliary_function_general_repel.R'))
V(net)$label = ''
lrange = apply(l, 2, range)
l2 = l
l2 = sweep(l2, 2, lrange[1,], '-')
l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

w = 4
png(paste0(imgpref, npref, '_inh_motif_regulon_network_repel.png'), res=300, units='in', width=w,height=w)
# pdf(paste0(imgpref, npref, '_inh_motif_regulon_network_repel.pdf'), width=w,height=w)
sp = 0.25
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# Repel points:
lbcex=0.5
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         max.overlaps=25,
                         labels=nodes, cex=lbcex, pt.cex=.25)
points(x=rdf$x[ntype == 2], y=rdf$y[ntype==2], col=ncol[ntype==2], pch=15, cex=1)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey50')
# legend('topright', legend=c(paste('Up in', subtype), paste('Down in', subtype)), pch=19, col=c(tsp.col('indianred',.5),tsp.col('royalblue',.5)), bty='n', cex=1)
# text(x=parpos(1,-.025), y=parpos(2,-.98), paste('Expression Corr. in', subtype, '(signif. genes, all regions)'), xpd=TRUE, cex=1, adj=0)
legend('bottomleft', legend=names(reg.cols), col=reg.cols, ncol=5,
       cex=.4, pch=15, pt.cex=1, inset=0, box.col=NA)
dev.off()








# ----------------------------------
# Community detection on this graph:
# ----------------------------------
cls = cluster_louvain(net)
memb = cls$membership
lnet = net
V(lnet)$color = snap.cols[memb + 19]
# Clusters:
cldf = data.frame(gene=nodes, cls=memb, up=(nodes %in% upgene))

png(paste0(imgpref, 'mastlmm_siggenes_graph_repel_louvain_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=7,height=7)
sp = 0.25
par(mar = rep(sp,4))
plot(lnet, layout=l, curved=F)
# Repel points:
lbcex=0.6
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey50')
legend('topright', legend=c(paste('Up in', subtype), paste('Down in', subtype)), pch=19, col=c(tsp.col('indianred',.5),tsp.col('royalblue',.5)), bty='n', cex=1)
text(x=parpos(1,-.025), y=parpos(2,-.98), paste('Expression Corr. in', subtype, '(signif. genes, all regions)'), xpd=TRUE, cex=1, adj=0)
dev.off()


# With background of up/down (not great)
sets = list()
sets[['up']] = upgene[upgene %in% nodes]
sets[['down']] = downgene[downgene %in% nodes]
pal = c(tsp.col('indianred',.75),tsp.col('royalblue',.75))

png(paste0(imgpref, 'mastlmm_siggenes_graph_repel_louvain2_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=7,height=7)
sp = 0.25
par(mar = rep(sp,4))
plot(lnet, layout=l, curved=F,
     mark.groups=sets, mark.border=pal, mark.col=NA, mark.shape=.75, mark.lwd=2)
# Repel points:
lbcex=0.6
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey50')
legend('topright', legend=c(paste('Up in', subtype), paste('Down in', subtype)), pch=19, col=c(tsp.col('indianred',.5),tsp.col('royalblue',.5)), bty='n', cex=1)
text(x=parpos(1,-.025), y=parpos(2,-.98), paste('Expression Corr. in', subtype, '(signif. genes, all regions)'), xpd=TRUE, cex=1, adj=0)
dev.off()


