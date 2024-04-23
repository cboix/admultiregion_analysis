#!/usr/bin/R
# ---------------------------------------
# Differential expression for astrocytes:
# Updated: 05/25/21
# ---------------------------------------
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
# subtype = 'Ast_GRM3'
subtype = 'Ast'
region = 'allregions'

# Data loader:
commandArgs = function(x){ c(celltype, subtype, region)}
source(paste0(bindir, 'multiRegion/load_difftl_data.R'))
gc()

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'ast-difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# -----------------------------
# Load in the regression files:
# -----------------------------
prefstr = paste0(celltype,'_',subtype)
fnlist = list.files(path=regdir, pattern=paste0('nebula_ruv.',prefstr, '.*rda'))
nsigdf = c(); alldf = c();
for (fn in fnlist){
    load(paste0(regdir, fn))
    region = nsig[1,'region']
    path = nsig[1,'path']
    resdf$region = region
    resdf$path = path
    ndf = data.frame(nsig)
    if (!('X1' %in% colnames(ndf))){ ndf$X1 = 0 }
    if (!('X2' %in% colnames(ndf))){ ndf$X2 = 0 }
    if (is.null(nsigdf)){
        nsigdf = ndf 
    } else {
        nsigdf = rbind(nsigdf, ndf[,colnames(nsigdf), drop=F])
    }
    fulldf$path = path
    fulldf$region = region
    # fulldf = fulldf[order(fulldf$col != 0, abs(fulldf$logFC), decreasing=T),]
    fulldf = fulldf[order(fulldf$p),]
    fulldf$rank = 1:nrow(fulldf)
    alldf = rbind(alldf, fulldf)
}
names(nsigdf)[names(nsigdf) == 'X1'] = 'ndown'
names(nsigdf)[names(nsigdf) == 'X2'] = 'nup'
nsigdf$nup = as.numeric(nsigdf$nup)
nsigdf$ndown = as.numeric(nsigdf$ndown)

resdf = alldf[alldf$region == 'allregions',]
cat(head(resdf[resdf$col == 1 & resdf$path == 'nft','gene'], 50))
cat(head(resdf[resdf$col == 2 & resdf$path == 'nrad','gene'], 50))

# ----------------------------------------------------------
# Make the full matrix at aggregate ct + ind + region level:
# ----------------------------------------------------------
mat = mat[,pathdf$barcode]
pathdf$ric = paste0(pathdf$region,'_', pathdf$projid, '_', pathdf$cell_type_high_resolution)
rdf = unique(pathdf[,c('ric','region','projid','cell_type_high_resolution','nrad','cogdxad','nft','plaq_n','plaq_d')])
rownames(rdf) = rdf$ric
rtypes = rdf$ric
tform = make.tform(pathdf$ric, u=rtypes, norm=T)
rctdf = agg.rename(barcode ~ ric, pathdf, length, 'count')
sub.rtypes = rctdf[rctdf$count > 10,'ric']
rind = which(rtypes %in% sub.rtypes)
# rind = 1:nrow(rdf)

# Normalize matrix, aggregate:
fact = (nc / median(nc))[colnames(mat)]
nmat = mat 
nmat@x <- nmat@x / rep.int(fact, diff(nmat@p))
nmat = log1p(nmat)
gc()
# Aggregate:
avg.mat = nmat %*% tform # Takes a while.
# avg.mat = mat %*% tform # Takes a while.
gc()

clsplit = paste0(rdf$cell_type_high_resolution, ' ', rdf$nrad)[rind]
clsplit = paste0(rdf$region, ' ', rdf$nrad)[rind]
ha = HeatmapAnnotation(CT=rdf$cell_type_high_resolution[rind], 
                       Region=rdf$region[rind],
                       nrad=rdf$nrad[rind],
                       nft=rdf$nft[rind],
                       plaq_n=rdf$plaq_n[rind],
                       plaq_d=rdf$plaq_d[rind],
                       cogdxad=rdf$cogdxad[rind],
                       col=list(Region=reg.cols,
                                nrad=colvals[['nrad']],
                                cogdxad=colvals[['cogdxad']],
                                CT=tcols))
# hb = rowAnnotation(Set=pset, col=list(Set=tcols[topec]))

reg = 'allregions'
path = 'nrad'
ntop = 250
for (path in c('nrad','nft','cogdxad','plaq_n')){
    resdf = alldf[alldf$region == reg,]
    resdf = resdf[resdf$path == path,]
    # resdf = resdf[order(resdf$col != 0, abs(resdf$logFC), decreasing=T),]
    labgenes = head(resdf[resdf$col == 2,'gene'], ntop)
    labgenes = labgenes[labgenes %in% rownames(mat)]
    pmat = as.matrix(avg.mat[labgenes,rind])
    pmat = t(scale(t(pmat)))
    # Individuals:
    png(paste0(imgpref, 'indheatmap_',prefstr, '_', path, '_up',ntop,'_',reg,'.png'), res=450, units='in', width=15, height=10)
    draw(Heatmap(pmat,
                 top_annotation=ha,
                 column_split=clsplit,
                 row_names_gp=gpar(fontsize=6),
                 show_column_names=FALSE))
    dev.off()

    # Correlation
    submat = as.matrix(mat[labgenes,pathdf$barcode])
    submat = sweep(submat, 2, pathdf$ncounts / 10000, '/')
    cr = cor(t(submat))
    diag(cr) = 0
    w = 11
    png(paste0(imgpref, 'corheatmap_',prefstr, '_', path, '_up',ntop,'_',reg,'.png'), res=450, units='in', width=w, height=w)
    draw(Heatmap(cr,
            row_names_gp=gpar(fontsize=8),
            column_names_gp=gpar(fontsize=8)))
    dev.off()
}


adind = pathdf$nrad == 'AD'
ctind = pathdf$nrad == 'CTRL'

gene = 'SYNGAP1'
wilcox.test(nmat[gene,ctind],
            nmat[gene,adind])

mean(nmat[gene,ctind])
mean(nmat[gene,adind])

labgenes = head(resdf[resdf$col == 2 & resdf$path == path,'gene'], 10)

x = nmat[labgenes,]
lgenes = c('HILPDA', 'ANGPTL4', 'FOS', 'IRS2', 'APOE',
           'CLU', 'GAPDH', 'PFKP', 'PFKFB3', 'PLOD2',
           'MT1X', 'MT2A', 'AQP1', 'AQP4', 'VEGFA', 'BHLHE40')


x = nmat[lgenes,]
x = data.frame(t(as.matrix(x)))
x$barcode = rownames(x)
x$nrad = pathdf[x$barcode, 'nrad']
x$region = pathdf[x$barcode, 'region']

xlong = gather(x, gene, val, -barcode, -nrad,-region)

ggplot(xlong, aes(nrad, val, fill=nrad)) + 
    facet_grid(gene~region) + 
    geom_violin() + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=colvals[['nrad']]) + 
    theme_pubr()

xlong$expr = xlong$val > 0
a = table(xlong[xlong$region == 'TH',c('gene','expr','nrad')])

# mdf = aggregate(val ~ gene + nrad, xlong[xlong$region == 'TH',], mean)
xlong$projid = pathdf[xlong$barcode,'projid']
xlong$nft = pathdf[xlong$barcode,'nft']
xlong$cogdxad = pathdf[xlong$barcode,'cogdxad']

mdf = aggregate(val ~ gene + nrad + cogdxad + projid, xlong[xlong$region == 'TH',], mean)

ggplot(mdf, aes(nrad, val, fill=nrad)) + 
    facet_wrap(~gene) + 
    geom_violin() + 
    geom_boxplot(width=.25) + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=colvals[['nrad']]) + 
    theme_pubr()

mdf = aggregate(val ~ nrad + cogdxad + projid + region, xlong, mean)

ggplot(mdf, aes(nrad, val, fill=nrad)) + 
    facet_wrap(~region) + 
    geom_violin() + 
    geom_boxplot(width=.25) + 
    # scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=colvals[['nrad']]) + 
    stat_compare_means() + 
    theme_pubr()



mdf = aggregate(val ~ barcode + nrad, xlong[xlong$region == 'TH',], mean)

ggplot(mdf, aes(nrad, val, fill=nrad)) + 
    geom_violin() + 
    geom_boxplot(width=.25) + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=colvals[['nrad']]) + 
    theme_pubr()



# -------------------------------------------

library(uwot)
u = umap(pmat)
udf = data.frame(u)
udf$gene = rownames(pmat)
ggplot(udf, aes(X1, X2, label=gene)) + 
    geom_point() + 
    geom_text_repel() + 
    theme_pubr()

sind = sample(1:nrow(pathdf), 1000, replace=F)
# sind = which(pathdf$ngenes > 2000)
sind = pathdf$barcode[sind]
pmat = as.matrix(nmat[labgenes, sind])
# pmat = scale(pmat)
pmat = t(scale(t(pmat), center=F))

png(paste0(imgpref, 'cellheatmap_',path, '_up',ntop,'_',reg,'.png'), res=450, units='in', width=15, height=7)
draw(Heatmap(pmat,
             # top_annotation=ha,
             col=viridis(100),
             column_split =pathdf[sind,'nrad'],
             row_names_gp=gpar(fontsize=8),
             show_column_names=FALSE))
dev.off()

cs = colSums(nmat[labgenes,])
plot(cs, pathdf$ngenes)


pass0 = (mat[labgenes,] > 0) %*% tform
pmat = as.matrix(pass0)
Heatmap(pmat,
        col=colb,
        top_annotation=ha,
        column_split=rdf$nrad,
        show_column_names=FALSE)


labgenes = head(resdf[resdf$col == 2 & resdf$path == 'plaq_n','gene'], 100)
labgenes = labgenes[labgenes %in% rownames(avg.mat)]

zmat = (mat[labgenes,pathdf$cell_type_high_resolution== 'Ast DCLK1'] > 0)
zz = zmat %*% t(zmat) # Jaccard:
zr = rowSums(zmat)
zr = matrix(rep(zr, length(zr)), length(zr), length(zr))
zu = zr + t(zr)
zj = as.matrix(zz / (zu - zz))
diag(zj) = NA
Heatmap(zj,
        col=colb,
        # top_annotation=ha,
        # column_split=rdf$nrad,
        show_column_names=FALSE)

# TODO: visual is not huge diff; look at coeff, and changes overall?

# Four sets with similar size + see where they are:
s1 = c('MT1X','MT1M','MT1E','FTH1')
s2 = c('ITM2B','ITM2C','CLU','PRNP')
s3 = c('GFAP','VIM','EMP1','CXCL14')
s4 = c('SLC38A2','SLC38A1','GLUL','GLS')
# zmat = (mat[c(s1,s2,s3,s4),] > 0)
zmat = (nmat[c(s1,s2,s3,s4),pathdf$barcode])

zmat.avg = as.matrix(zmat %*% tform)[,rind]
pmat = t(scale(t(zmat.avg)))
Heatmap(pmat,
        column_split=rdf$nrad[rind],
        show_column_names=F,
        top_annotation=ha,
)# ITM2C, ITM2B, PRNP, FTH1, GLS, VIM, CLU are module separate from GFAP, EMP1, SLC38A1 (at ind-level, can't really tell othw...)o

cr = cor(t(pmat))
Heatmap(cr)

pathdf$s1 = colSums(zmat[s1,])
pathdf$s2 = colSums(zmat[s2,])
pathdf$s3 = colSums(zmat[s3,])
pathdf$s4 = colSums(zmat[s4,])

sdf = aggregate(cbind(s1,s2,s3,s4) ~ projid + region + cell_type_high_resolution + nrad + nft + plaq_n, pathdf, mean)

ggplot(sdf, aes(cell_type_high_resolution, log10(s3), fill=nrad)) + 
    geom_violin() + 
    geom_boxplot() + 
    theme_pubr()


# Subset to relevant genes
# Co-expression - pathways
# Comparison to other DE genes (unique genes)

# ----------------------------------------
# Modules from co-expression --> networks:
# ----------------------------------------
reg = 'allregions'
path = 'nrad'
ntop = 250
s1 = c('MT1X','MT1M','MT1E','FTH1')
s2 = c('ITM2B','ITM2C','CLU','PRNP')
s3 = c('GFAP','VIM','EMP1','CXCL14')
s4 = c('SLC38A2','SLC38A1','GLUL','GLS')
slist = c(s1,s2,s3,s4)

resdf = alldf[alldf$region == reg,]
labgenes = head(resdf[resdf$col == 2 & resdf$path == path,'gene'], ntop)
labgenes = labgenes[labgenes %in% rownames(mat)]
labgenes = labgenes[grep("^MT-",labgenes, invert=TRUE)] # Remove MT genes
submat = as.matrix(mat[labgenes,pathdf$barcode])
submat = sweep(submat, 2, pathdf$ncounts / 10000, '/')
# kbcs = pathdf$barcode[pathdf$cell_type_high_resolution == 'Ast GRM3']
cr = cor(t(submat))
diag(cr) = 0
# Heatmap(cr, row_names_gp=gpar(fontsize=8))

# pmat = as.matrix(avg.mat[labgenes,rind])
# pmat = t(scale(t(pmat)))
# # pmat = t(scale(t(pmat), center=FALSE))

# Plot a graph/extract clusters/modules (community detect):
library(igraph)
# cutoff = 0.075
if (subtype == 'Ast_LUZP2'){
    cutoff = 0.12
} else if (subtype == 'Ast_DCLK1'){
    cutoff = 0.1
} else if (subtype == 'Ast_DPP10'){
    cutoff = 0.1
} else { cutoff = 0.08 }

dmat = cr
rn = rownames(dmat)
dmat = data.frame(dmat)
colnames(dmat) = rn
rownames(dmat) = rn
dmat$T1 = rownames(dmat)
ddf  = gather(dmat, T2, sim, -T1)
ddf = ddf[ddf$sim >= cutoff,]
ddf = ddf[ddf$T1 != ddf$T2,]

# Remove links if diff sets:
# ddf$samedir = ((ddf$T1 %in% upgene) & (ddf$T2 %in% upgene)) | ((ddf$T1 %in% downgene) & (ddf$T2 %in% downgene))
# ddf = ddf[ddf$samedir,]
nodes = sort(unique(c(ddf$T1, ddf$T2)))
# Remove edges in opposite direction
ddf$T1 = factor(ddf$T1, levels=nodes)
ddf$T2 = factor(ddf$T2, levels=nodes)
ddf = ddf[as.numeric(ddf$T1) < as.numeric(ddf$T2),]
# Kept pct:
dim(ddf)[1]
dim(ddf)[1] / (dim(dmat)[1]^2)
ddf$COLOR = 'grey25'

# Simple network: just the links/points:
library(igraph)
sdf = ddf
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
vcol = rep(tsp.col('indianred',.5), length(nodes))
# vcol[nodes %in% downgene] = tsp.col('royalblue', .5)
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(net)$size = 2
V(net)$label = nodes
V(net)$label.cex = .5
V(net)$label.color = 'black'
V(net)$color = vcol
V(net)$frame.color <- 'black' # vcol
V(net)$frame.color <- NA
V(net)$pch = 19
E(net)$color = ecol 
elty = rep('dotted', length(sdf$sim))
elty[sdf$sim >= .85] = 'dashed'
elty[sdf$sim >= .95] = 'solid'
E(net)$lty = elty
E(net)$width = sdf$sim * 4
E(net)$weight = sdf$sim  * 4
set.seed(8)
l <- layout_with_fr(net, grid='nogrid') # Usually best

# For repel:
source(paste0('~/ENCODE_DATA/bin/', 'auxiliary_function_general_repel.R'))
V(net)$label = ''
V(net)$size = 2.25
l2 = l
lrange = apply(l, 2, range)
l2 = sweep(l2, 2, lrange[1,], '-')
l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

png(paste0(imgpref, 'crnet_repel',path, '_up',ntop,'_',reg,'.png'), res=300, units='in', width=7,height=7)
sp = 0.25
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# Repel points:
lbcex=0.6
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         max.overlaps=25,
                         labels=nodes, cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey50')
legend('topright', legend=c(paste('Up in', subtype), paste('Down in', subtype)), pch=19, col=c(tsp.col('indianred',.5),tsp.col('royalblue',.5)), bty='n', cex=1)
text(x=parpos(1,-.025), y=parpos(2,-.98), paste('Expression Corr. in', subtype, '(signif. genes, all regions)'), xpd=TRUE, cex=1, adj=0)
dev.off()

# ----------------------------------
# Community detection on this graph:
# ----------------------------------
cls = cluster_louvain(net)
memb = cls$membership
lnet = net
V(lnet)$color = snap.cols[memb + 19]
# Clusters:
cldf = data.frame(gene=nodes, cls=memb)
sapply(unique(memb), function(x){ cat(x, "\n",cldf$gene[cldf$cls == x],"\n")})

png(paste0(imgpref, 'crnet_repel_louvain_',path, '_up',ntop,'_',reg,'.png'), res=300, units='in', width=7,height=7)
sp = 0.25
par(mar = rep(sp,4))
plot(lnet, layout=l, curved=F)
# Repel points:
lbcex=0.6
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         max.overlaps=25,
                         labels=nodes, cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey50')
legend('topright', legend=c(paste('Up in', subtype), paste('Down in', subtype)), pch=19, col=c(tsp.col('indianred',.5),tsp.col('royalblue',.5)), bty='n', cex=1)
text(x=parpos(1,-.025), y=parpos(2,-.98), paste('Expression Corr. in', subtype, '(signif. genes, all regions)'), xpd=TRUE, cex=1, adj=0)
dev.off()


if (subtype == 'Ast_LUZP2'){
    lgenes = c('HILPDA', 'ANGPTL4', 'FOS', 'IRS2', 'APOE',
               'CLU', 'GAPDH', 'PFKP', 'PFKFB3', 'PLOD2',
               'MT1X', 'MT2A', 'AQP1', 'AQP4', 'VEGFA', 'BHLHE40')

    labgenes = lgenes
    labgenes = labgenes[labgenes %in% rownames(mat)]

    labgenes = head(resdf[resdf$col == 2 & resdf$path == path,'gene'], ntop)
    labgenes = labgenes[labgenes %in% rownames(mat)]
    labgenes = labgenes[grep("^MT-",labgenes, invert=TRUE)] # Remove MT genes
    labgenes = unique(c(lgenes, labgenes))

    labgenes = c('CD44','GFAP','HIF1A','F3','MAOB',lgenes)

    tbc = pathdf$barcode[pathdf$region == 'TH']
    adsplit = as.character(pathdf[tbc,'nrad'])
    ha = HeatmapAnnotation(nrad=pathdf[tbc,'nrad'],
                           cogdxad=pathdf[tbc,'cogdxad'],
                           e4=pathdf[tbc,'Apoe_e4'],
                           col=list(nrad=colvals[['nrad']],
                                    e4 = c('no'='white','yes'='grey50'),
                                    cogdxad=colvals[['cogdxad']]))
    pmat = as.matrix(nmat[labgenes,tbc])
    pmat = t(scale(t(pmat)))
    # Cell-level:
    png(paste0(imgpref, 'cellTHheatmap2_',prefstr, '_', path, '_up',ntop,'_',reg,'.png'), res=450, units='in', width=25, height=4.5)
    draw(Heatmap(pmat,
                 top_annotation=ha,
                 column_split=adsplit,
                 # row_names_gp=gpar(fontsize=6),
                 row_split=rownames(pmat) %in% lgenes,
                 show_column_names=FALSE))
    dev.off()


    x = nmat[lgenes,tbc]
    x = data.frame(t(as.matrix(x)))
    x$barcode = rownames(x)
    x$nrad = pathdf[x$barcode, 'nrad']
    x$projid = pathdf[x$barcode, 'projid']
    x$e4 = pathdf[x$barcode, 'Apoe_e4']
    xlong = gather(x, gene, val, -barcode, -nrad,-e4,-projid)

    ggplot(xlong, aes(nrad, val, fill=nrad, color=e4)) + 
        facet_wrap(~gene) + 
        geom_violin() + 
        scale_y_continuous(expand=c(0,0)) + 
        scale_fill_manual(values=colvals[['nrad']]) + 
        theme_pubr()

    mdf = aggregate(val ~ gene + nrad + e4 + projid, xlong, mean)
    mdf = mdf[mdf$gene %in% c('GAPDH','HILPDA','IRS2','PFKFB3','PFKP','PLOD2'),]

    gplot = ggplot(mdf, aes(gene, val, fill=nrad, alpha=e4)) + 
        geom_boxplot(outlier.shape=NA) + 
        geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
        scale_y_continuous(expand=c(0,0)) + 
        scale_fill_manual(values=colvals[['nrad']]) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'metab_genes_', prefstr, '_nrad_boxplots.pdf'), gplot, dpi=400, units='in', width=6,height=4)
}



