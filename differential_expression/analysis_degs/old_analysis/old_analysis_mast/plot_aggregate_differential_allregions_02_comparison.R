#!/usr/bin/R
# -----------------------------------------------------------
# Aggregate and plot chunked runs of MAST + RE for DGE
# Plot the comparison between regressions on different pathologies.
# Updated: 02/11/2021
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)
library(MAST)
library(data.table)
library(viridis)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

celltype = 'Ast'
subtype = 'Ast'
path = 'nrad'

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype subtype region chunksize chunk ascertainment")
} else {        
    celltype = args[1]
    subtype = args[2]
    region = args[3]
    path = args[4]
}

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# ------------------
# Load the metadata:
# ------------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

# Data directories:
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    # matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')
ststr = gsub("/","_",gsub(" ","_", subtype))
cellstr = gsub("/","_",gsub(" ","_", celltype))

# --------------------------------------------
# Load from various different regression runs:
# --------------------------------------------
agg.pref = paste0(regdir, prefix, '.mastlmm_reg.allpath.allreg.major.', cellstr, '.minor.', ststr)
agg.file = paste0(agg.pref, '.tsv.gz') 
agg.rda = paste0(agg.pref, '.rda') 
if (!file.exists(agg.rda)){
    alldf = c()
    sumdf = c()
    for (path in c('nrad','nft','plaq_d','plaq_n')){
        print(paste("[STATUS] Loading regression on", subtype, 'across regions', 'on var', path))
        fpref = paste0(prefix, '.mastlmm_reg.', path, '.allreg.major.', cellstr, '.minor.', ststr)
        # fnlist = list.files(pattern=paste0(fpref,'.*633.*.Rda'), path=regdir)
        fnlist = list.files(pattern=paste0(fpref,'.*.160.*.Rda'), path=regdir)
        if (length(fnlist) > 0){
            print(paste("[STATUS] Loading from",length(fnlist),"files."))
            subdf = c()
            for (fn in fnlist){
                # print(fn)
                load(paste0(regdir,fn))
                regdf$path = path 
                summaryDt$path = path 
                subdf = rbind(subdf, regdf)
                sumdf = rbind(sumdf, summaryDt)
            }
            names(subdf)[2] = 'pvalue'
            subdf = subdf[order(subdf$pvalue),]
            subdf$padj = p.adjust(subdf$pvalue, 'fdr')
            alldf = rbind(alldf, subdf)
        } else {
            print("No files..")
        }
    }
    dim(alldf)
    write.table(alldf, file=gzfile(agg.file), quote=F, sep="\t", col.names=F)
    save(alldf, file=agg.rda)
} else { 
    load(agg.rda)
}

FCTHRESH=0.02 
alldf$col = 1 * (alldf$padj < 0.05) * (1 + 1 * (alldf$coef > 0)) * (abs(alldf$coef) > ifelse(alldf$path == 'nrad', 1,1/50) * FCTHRESH)

nsigdf = data.frame(table(alldf[,c('col','path')]))
nsigdf = spread(nsigdf, col, Freq)
names(nsigdf) = c('path','none','down','up')
tcols = brewer.pal(12, 'Paired')

gplot = ggplot(nsigdf, aes(path, up)) + 
    geom_bar(stat='identity', fill=tcols[6]) + 
    geom_bar(data=nsigdf, aes(path, -down), stat='identity', fill=tcols[2]) + 
    theme_pubr() + coord_flip() + 
    labs(x='AD variable', y='Number of DEGs', title=paste0(celltype, ' (', subtype, ')'))
ggsave(paste0(imgpref, 'mastlmm_comparison_ndeg.allreg.major.', cellstr, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=2.5,height=1.5)


# -----------------------------------------------------------------------------------------------
# Compare which genes are consistently up/down across attributes, which aren't, and which switch:
# -----------------------------------------------------------------------------------------------
rd = alldf[,c('primerid','coef','padj','path')]
alldf$log10p = -log10(alldf$padj)
pwide = spread(alldf[,c('primerid','log10p','path')],path, log10p, fill=0)
cwide = spread(alldf[,c('primerid','coef','path')],path, coef, fill=0)

pmat = as.matrix(pwide[,-1])
rownames(pmat) = pwide[,1]
pcut = 35
Z = 1 * (pmat > pcut)
pmat = pmat[apply(Z,1,sum) > 1,]
pmat[pmat > 100] = 100
# pmat = reord(pmat)

cmat = as.matrix(cwide[,-1])
rownames(cmat) = cwide[,1]
# cmat = cmat[apply(cmat,1,max) > 3,]
cmat = cmat[rownames(pmat),]
cmat = cmat * (1 * (pmat > pcut))
cmat = sweep(cmat, 2, apply(abs(cmat), 2, max), '/')
# cmat[cmat > 100] = 100
# cmat = reord(cmat)

ptstr = paste0(path, '.allreg.major.', cellstr, '.minor.', ststr)
png(paste0(imgpref, 'mastlmm_topgenes_comp_', ptstr, '.png'), res=300, units='in', width=3.4,height=8)
Heatmap(cmat, name='Scaled\ncoef',
        cluster_rows = TRUE, 
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        show_column_names = TRUE)
dev.off()


# ---------------------------------------
# Plot the volcano plot for one of these:
# ---------------------------------------
path = 'nft'
resdf = alldf[alldf$path == path,]
pltdf = resdf

if (path %in% c('nft','plaq_d','plaq_n')){
    scale = 1 / 50
} else { 
    scale = 1
}
if (subtype %in% c('Mic','Ast','Opc','Oli','Exc','Inh')){
    FCTHRESH=0.02 * scale
} else {
    FCTHRESH=0.05 * scale
}
pltdf$col = 1 * (pltdf$padj < 0.05) * (1 + 1 * (pltdf$coef > 0)) * (abs(pltdf$coef) > FCTHRESH)
mx = 0.5
pltdf$coef[pltdf$coef > mx] = mx
pltdf$coef[pltdf$coef < -mx] = -mx
if (subtype == 'Mic'){
    pltdf$pvalue[pltdf$pvalue < 1e-100] = 1e-100
} else if (subtype == 'Ast'){
    pltdf$pvalue[pltdf$pvalue < 1e-150] = 1e-150
}

ntop = sum(pltdf$padj < 0.05 & pltdf$col !=0, na.rm=T)
labcut = 0.05
print(paste("Number of sig. genes:", ntop))
if (ntop > 50){
    labcut = 0.02
}
if (ntop > 100){
    NG = 70
    # Alternatively, plot top each:
    downdf = pltdf[(pltdf$padj < labcut) & pltdf$col == 1,]
    updf = pltdf[(pltdf$padj < labcut) & pltdf$col == 2,]
    labdf = rbind(head(downdf, NG/2), head(updf, NG/2))
} else {
    labdf = pltdf[(pltdf$padj < labcut) & pltdf$col != 0,]
}
labdf = unique(rbind(labdf, pltdf[abs(pltdf$coef) > FCTHRESH * 10 & pltdf$col != 0,]))

gcols = brewer.pal(12, 'Paired')
gplot = ggplot(pltdf, aes(coef, -log10(pvalue), color=factor(col))) + 
    geom_point(cex=.5) + 
    # geom_text_repel(data=labdf, aes(label=primerid), size=2.5, segment.size=.25, max.overlaps=20) + 
    geom_text_repel(data=labdf, aes(label=primerid), size=2, segment.size=.25, max.overlaps=15) + 
    theme_pubr() + 
    theme(legend.position='none') + 
    labs(x='Coefficient (AD / CTRL by NIA-Reagan score)', y='-log10(p-value)', title=paste0(gsub("_"," ",subtype), ' - Differential Genes (All Regions) - ', path)) + 
    # scale_color_manual(values=c('grey75','royalblue','indianred')) + 
    scale_color_manual(values=c('grey75',gcols[1],gcols[5])) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_vline(xintercept=0, lty='dashed') + 
    geom_vline(xintercept=c(-FCTHRESH, FCTHRESH), lty='dashed', col='grey50',lwd=.25) 
# ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.allreg.major.', cellstr, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=6,height=6)
ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.allreg.major.', cellstr, '.minor.', ststr,'_small.png'), gplot, dpi=400, units='in', width=5,height=5)
ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.allreg.major.', cellstr, '.minor.', ststr,'_tiny.png'), gplot, dpi=400, units='in', width=2.5,height=3)
ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.allreg.major.', cellstr, '.minor.', ststr,'_tiny.pdf'), gplot, dpi=400, units='in', width=2.5,height=3)


# -------------------
# Run GO enrichments:
# -------------------
library(gprofiler2)
upgene = pltdf[pltdf$col == 2,'primerid']
downgene = pltdf[pltdf$col == 1,'primerid']
upgene = upgene[!is.na(upgene)]
downgene = downgene[!is.na(downgene)]

udf = gost(upgene)$result
ddf = gost(downgene)$result
adf = gost(c(upgene, downgene))$result
udf = udf[udf$term_size < 1000,]
ddf = ddf[ddf$term_size < 1000,]
adf = adf[adf$term_size < 1000,]

keep.cols = c('p_value','source','term_name')
udf = udf[c(grep("GO",udf$source), grep("REAC",udf$source)), ]
ddf = ddf[c(grep("GO",ddf$source), grep("REAC",ddf$source)),]
adf = adf[c(grep("GO",ddf$source), grep("REAC",adf$source)),]
udf = udf[order(udf$p_value),]
ddf = ddf[order(ddf$p_value),]
adf = adf[order(adf$p_value),]

# Plot up-regulated:
pdf = head(udf, 20)
pdf$term_name = factor(pdf$term_name, levels=rev(unique(pdf$term_name)))
gplot = ggplot(udf, aes(intersection_size / query_size, -log10(p_value), size=term_size, fill=source)) + 
    geom_point(pch=21, alpha=.5) + 
    geom_text_repel(data=pdf, aes(label=term_name), size=2.5, point.padding=.25, segment.size=.25, box.padding=.1) +
    scale_x_continuous(labels=scales::percent) + 
    # scale_y_continuous(expand=c(0,0)) + 
    labs(x='Percent sig. genes in term', y='-log10(p-value)', title=paste0(gsub("_"," ",subtype), ' - Up-regulated Genes (All Regions)')) + 
    theme_pubr() + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'mastlmm_GObubble_up_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=5,height=3)

# Plot down-regulated:
pdf = head(ddf, 20)
pdf$term_name = factor(pdf$term_name, levels=rev(unique(pdf$term_name)))
gplot = ggplot(ddf, aes(intersection_size / query_size, -log10(p_value), size=term_size, fill=source)) + 
    geom_point(pch=21, alpha=.5) + 
    geom_text_repel(data=pdf, aes(label=term_name), size=2.5, point.padding=.25, segment.size=.25, box.padding=.1, max.iter=5000) +
    scale_x_continuous(labels=scales::percent) + 
    # scale_y_continuous(expand=c(0,0)) + 
    labs(x='Percent sig. genes in term', y='-log10(p-value)', title=paste0(gsub("_"," ",subtype), ' - Down-regulated Genes (All Regions)')) + 
    theme_pubr() + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'mastlmm_GObubble_down_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=5,height=3)

cgenes = scan('cholesterol_biosynthesis_genes.txt','c', quiet=T)
resdf[resdf$primerid %in% cgenes,]

# ------------------------------
# Load in data from all regions:
# ------------------------------
agenes = c(upgene, downgene)
keep.reg = regions[regions != 'MB']
amat = c()
barcodes = c()
for (region in keep.reg){
    ststr = gsub("/","_",gsub(" ","_", subtype))
    matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                     celltype,'.',ststr,'.',region)
    rdafile = paste0(matpref, '.rda')  # In Matrix format
    # Load `mat` from rdafile:
    load(rdafile)
    mat = mat[agenes,]
    print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
    barcodes = c(barcodes, colnames(mat))
    genes = rownames(mat)
    ngenes = nrow(mat)
    amat = cbind(amat, mat)
}
rm(mat)
gcout = gc()

# ------------------------------
# Load the appropriate metadata:
# ------------------------------
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

# Normalize:
mat = as.matrix(amat[agenes,])
norm = sweep(mat, 2, marg[colnames(mat),'count'] / 100000,'/')
norm = log(norm + 1)
rm(amat)
gcout = gc()

rownames(cellmeta) = cellmeta$barcode
pathdf = cellmeta[barcodes,]
pathdf = merge(pathdf, metadata[,c('projid','rind','age_death','msex','pmi', 'Apoe_e4', 'cogdx', 'niareagansc')])
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode
pathdf$nrad = 'CTRL'
pathdf$nrad[pathdf$niareagansc %in% c(1,2)] = 'AD'
pathdf$nrad = factor(pathdf$nrad, levels=c('CTRL','AD'))

pathlist = c('nft','plaq_d','plaq_n')
# Get the pathology mapped to each region:
regmap = c('AG','HC','PFC','MT','EC')
names(regmap) = c('ag','hip','mf','mt','ec')
pqdf = NULL
for (path in pathlist){
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    sdf = slong[,c('rind','value','region')]
    names(sdf)[2] = path
    if (is.null(pqdf)){
        pqdf = sdf
    } else {
        pqdf = merge(pqdf, sdf)
    }
}
pathdf = merge(pathdf, pqdf, all.x=TRUE) # Othw. drop out TH region


# --------------------------------------------------------------------------------------
# Plot the expression of the differential genes against each other (together/submodules)
# - Heatmaps (directly
# - Correlation of deg
# - Per-region
# - Per other path/etc. variables
# Ask if the DE gene expression is also significantly affected by other vars:
# --------------------------------------------------------------------------------------
# Order genes:
if (ncol(norm) < 10000){
    rmat = reord(norm)
    # Order cells:
    cellscore = colSums(norm[upgene,]) - colSums(norm[downgene,])
    pathdf$cs = cellscore
    pathdf = pathdf[order(pathdf$cs),]
    pathdf = pathdf[order(pathdf$nrad),]
    # pathdf = pathdf[order(pathdf$region),]
    sbc = pathdf$barcode
    rmat = rmat[,sbc]

    colv = viridis(100)
    png(paste0(imgpref, 'mastlmm_siggenes_heatmap_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=3,height=6)
    par(mar=rep(0.2,4))
    image(rmat, col=colv, useRaster=T, axes=F)
    box(lwd=.25)
    dev.off()
}

# Correlation patterns:
cr = cor(t(norm))
mx = 0.25
rg = rownames(cr)
rgenes = c(rg[rg %in% upgene], rg[rg%in% downgene])
crmat = cr[rgenes,rgenes]
crmat[crmat > mx] = mx
crmat[crmat < -mx] = -mx
rows <- (seq(0,ncol(cr), 50) - .5) / (ncol(cr) -1)

png(paste0(imgpref, 'mastlmm_siggenes_cormap_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=12,height=11.5)
par(mar=c(.25,2,.25,.25))
image(crmat, col=rev(colrb), zlim=c(-mx,mx), useRaster=T, axes=F)
yat=seq(0,1,length.out=nrow(crmat))
text(x=parpos(1,0.001), yat, rownames(crmat), xpd=TRUE, adj=1, cex=.25)
box(lwd=.25)
abline(v=rows, lwd=.2)
abline(h=rows, lwd=.2)
dev.off()

# ---------------------------------------------------------
# Plot a graph/extract clusters/modules (community detect):
# ---------------------------------------------------------
if (subtype == 'CAMs'){
    cutoff = 0.18
} else if (subtype == 'T_cells'){
    cutoff = .24
} else if (subtype == 'End'){ 
    cutoff = 0.18
} else if (subtype == 'Fib'){ 
    cutoff = 0.20
} else if (subtype == 'Ast'){ 
    cutoff = 0.23
} else if (subtype == 'Opc'){ 
    cutoff = 0.19
} else if (subtype == 'Oli'){ 
    cutoff = 0.15
} else { 
    cutoff = 0.15
}

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
ddf$samedir = ((ddf$T1 %in% upgene) & (ddf$T2 %in% upgene)) | ((ddf$T1 %in% downgene) & (ddf$T2 %in% downgene))
ddf = ddf[ddf$samedir,]
all.nodes=FALSE
if (all.nodes){
    nodes = sort(guid)
    npref = paste0(npref, '_allnodes')
} else {
    nodes = sort(unique(c(ddf$T1, ddf$T2)))
}
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
vcol[nodes %in% downgene] = tsp.col('royalblue', .5)
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
E(net)$weight = sdf$sim * 2
set.seed(8)
l <- layout_with_fr(net, grid='nogrid') # Usually best

png(paste0(imgpref, 'mastlmm_siggenes_graph_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=7,height=7)
# pdf(paste0(imgpref, npref, '.pdf'), width=6, height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# plot(net, layout=l, edge.curved=seq(-0.5, 0.5, length = ecount(net)))
dev.off()

# For repel:
source(paste0('~/ENCODE_DATA/bin/', 'auxiliary_function_general_repel.R'))
V(net)$label = ''
V(net)$size = 2.25
l2 = l
lrange = apply(l, 2, range)
l2 = sweep(l2, 2, lrange[1,], '-')
l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

png(paste0(imgpref, 'mastlmm_siggenes_graph_repel_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=7,height=7)
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


# ------------------------------------------------------------------
# Take a look at some top genes by more variables (not mixed model):
# ------------------------------------------------------------------
library(lme4)
library(MASS)
gene = 'HSP90AA1'
# gene = 'VEGFA'
# gene = 'GFAP'
# gene = 'SDCBP'

mx = marg[colnames(norm),'count']
# x = amat[gene,]
x = norm[gene,]

df = data.frame(val=x,barcode=colnames(norm), marg=mx)
df$rind = cellmeta[df$barcode,'rind']
df$projid = cellmeta[df$barcode,'projid']
df$region = cellmeta[df$barcode,'region']
df$ct = cellmeta[df$barcode,'cell_type_high_resolution']
df = merge(df, pqdf, all.x=TRUE)
df = merge(df, unique(metadata[,c('rind','Apoe_e4','msex','niareagansc','cogdx','pmi','age_death')]))
df$nrad = 'CTRL'
df$nrad[df$niareagansc %in% c(1,2)] = 'AD'
df$nrad = factor(df$nrad, levels=c('CTRL','AD'))
# Rescale:
df$nft = df$nft / max(df$nft, na.rm=T)
df$plaq_n = df$plaq_n / max(df$plaq_n, na.rm=T)
df$plaq_d = df$plaq_d / max(df$plaq_d, na.rm=T)
df$age_death = df$age_death / 100
# Normalized:
df$norm = with(df, log(1e3 * (val / marg) + 1))

# PLOT: with / without 0s:
# gplot = ggplot(df, aes(nrad, norm, color=nrad)) + 
#     facet_grid(. ~ ct) + geom_violin() + 
#     stat_compare_means() + theme_pubr()

spread(aggregate(norm ~ ct + nrad, df, mean), nrad, norm)
spread(aggregate(norm ~ region + nrad, df, mean), nrad, norm)

# fit = glm.nb(val ~ offset(log(marg)) + nrad + msex + pmi + age_death + region + ct, df)
path = 'plaq_n'
# fit = glm.nb(asform(c('val ~ offset(log(marg)) + ',path,'+ msex + pmi + age_death + region + ct')), df)
fit = glm.nb(asform(c('val ~ offset(log(marg)) + ',path,'+ msex + pmi + age_death + region')), df)
cfit = coefficients(summary(fit))
cat(round(cfit[path,1],4), -log10(cfit[path,4]), '\n')

# Try to compare with nebula:
library(nebula)
path = 'plaq_n'
genes = c('HSP90AA1','GFAP','RGMA','CRYAB','HMGB1','IFITM3','IFI44L','SLC39A11', 'LINGO1')
genes = c('CTGF','SDCBP','MT2A','HSP90AA1','GPCPD1','RHOJ','VEGFA','LDHA', 'LINGO1')
rownames(df) = df$barcode
df$pr = paste0(df$projid, df$region)
df = df[colnames(norm),]
sdf = df[!is.na(df[[path]]),]
sdf = sdf[order(sdf$pr),]
mat = amat[genes, rownames(sdf)]

# mdx = model.matrix(asform(c('~',path, '+ region + age_death + msex + pmi + ct')), data=sdf)
mdx = model.matrix(asform(c('~',path, '+ region + age_death + msex + pmi')), data=sdf)
# re = nebula(mat, as.character(sdf$projid), pred=mdx, offset=log(sdf$marg)) 
re = nebula(mat, as.character(sdf$pr), pred=mdx, offset=log(sdf$marg)) 
# re = nebula(mat, rep(1, nrow(sdf)), pred=mdx, offset=log(sdf$marg)) 
rdf = re$summary
pv = paste0('p_',path)
ev = paste0('logFC_',path)
rdf[,c('gene',ev, pv)]

# Output files:

# fit = lm(norm ~ nrad + ct + msex + pmi + nft + plaq_n + plaq_d + age_death + Apoe_e4, df)
fit = glm(val ~ offset(log(marg)) + nrad + msex + pmi + age_death + nrad + region + ct, df, family='poisson')
fit = glm(val ~ offset(log(marg)) + nft + msex + pmi + age_death + nrad + region + ct, df, family='poisson')
fit = glm(val ~ offset(log(marg)) + plaq_n * ct + msex + pmi + age_death + nrad + region + ct, df, family='poisson')
fit = glm(val ~ offset(log(marg)) + plaq_n * ct *region + msex + pmi + age_death + nrad + region + ct, df, family='poisson')
# fit = glmer.nb(val ~ offset(log(marg)) + nft + msex + pmi + age_death + nrad + region + ct + (1 | projid), df)

# ------------------------
# Averages by individuals:
# ------------------------
gcout = gc()
rownames(pathdf) = pathdf$barcode
pathdf = pathdf[colnames(norm),]
tform = make.tform(pathdf$rind, norm=T)

# Average normalized matrix
# norm.avg = norm %*% tform
# amat.avg = amat %*% tform
# rmat.avg = reord(norm.avg)
# rmat.avg = t(reord(t(rmat.avg)))

# Scaled for relative expr:
snorm = t(scale(t(norm)))
snorm.avg = snorm %*% tform
# rm(snorm)
# gcout = gc()

smat.avg = scale(t(rmat.avg))
pltmat = snorm.avg
pltmat = t(smat.avg)
# pltmat = t(scale(t(mat.avg)))

cldf$dir = 'Down'
cldf$dir[cldf$up] = 'Up'
cldf = cldf[order(cldf$cls),]
cldf = cldf[order(cldf$up),]
scols = snap.cols[1:max(cldf$cls) + 19]
names(scols) = as.character(1:max(cldf$cls))
pltmat = pltmat[cldf$gene,]
pltmat = pltmat[,order(colnames(pltmat))]
rad = metadata[colnames(pltmat),'niareagansc']
pltmat = pltmat[,order(rad)]
rmeta = metadata
rmeta$nft = NULL
rmeta$plaq_n = NULL
rmeta$plaq_d = NULL
rmeta = merge(rmeta, pqdf, all.x=TRUE)
rownames(rmeta) = rmeta$rind
rmeta = rmeta[colnames(pltmat),]
rmeta$nrad = ifelse(rmeta$niareagansc > 2,"CTRL","AD")
rmeta[rmeta$region == 'TH',pathlist] = 0

# Additional cut
rg = rownames(pltmat)
if (subtype %in% c('Ast','Opc','Mic')){
    pcut = 1e-10
    upgene = resdf[resdf$col == 2  & resdf$fdr < pcut,'primerid']
    downgene= resdf[resdf$col == 1  & resdf$fdr < pcut,'primerid']
    rg = rg[rg %in% c(upgene,downgene)]
    pltmat = pltmat[rg,]
}

# TODO: Diff genes by celltype
if (subtype == 'Mic'){
    labgenes = c('HLA-DRB1','APOC1','C1QC','CD74','IFI44L','MX1','APOE',
                 'UBC','HSPA1A','HIF1A','TLR2','GLDN','SYNDIG1','GLDN',
                 'IL15','IL4R','HDAC4','PTEN','DNAJA4','APPL2','DHFR',
                 'NPL','TMSB4X','RPL19','BIN1','C3','CSF1R','CX3CR1','INPP5D','RASGEF1C',
                 'LINGO1','MT-ND3','MT-CO1','HSP90B1')
} else if (subtype == 'Ast'){
    labgenes = c('MACF1','LINGO1','APPL2','ROBO1','MT1M','GFAP','XAF1','SLC39A11','HMBOX1','AQP4','CD44','FOS','ERBIN', 'CRYAB','HMGB1','HSP90AB1', 'S100B','MT1M','MT1E','MT1X','F3','MT-ND1')
} else if (subtype == 'Opc'){
    labgenes = c('VCAN','FGF14','PREX2','DOCK5','SLC35F1','APOD','GRIA4','FTL','FTH','MT1X','HSPA1A','GPR158','LUZP2','EGFR','HIF3A','SGK1','NRCAM','NCAM2','NRXN','CD81','ARL4C','ASTN2','TNR','MT-ND2','LINGO1','RASGEF1B','CRYAB','FXR','HMGB1','OLIG1','OMG','SLC22A17')
}
labgenes = labgenes[labgenes %in% rownames(pltmat)]
lat = sapply(labgenes, function(x){which(rownames(pltmat) == x)})

# Gene annotation:
rownames(cldf) = cldf$gene
ha = rowAnnotation(Dir = cldf[rg,'dir'], 
                   Cls = cldf[rg,'cls'],
                   Gene = anno_mark(at=lat, labels = labgenes),
                   col = list(Cls=scols,
                              Dir=c('Down'='royalblue','Up'='indianred')))
# Region annotation:
hb = HeatmapAnnotation(Region = rmeta[,'region'], 
                       AD = rmeta[,'nrad'],
                       nft = rmeta[,'nft'],
                       plaq_n = rmeta[,'plaq_n'],
                       plaq_d = rmeta[,'plaq_d'],
                       braak = rmeta[,'braaksc'],
                       cogdx = rmeta[,'cogdx'],
                       col = list(Region=reg.cols,
                                  braak=colvals[['braaksc']],
                                  AD=c('CTRL'='royalblue','AD'='indianred')))


png(paste0(imgpref, 'scaled_heatmap_siggenes_split_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=400, units='in', width=8, height=8)
Heatmap(pltmat, name = "mat", 
        cluster_rows = FALSE, 
        cluster_columns=TRUE,
        use_raster=T,
        right_annotation = ha,
        top_annotation = hb,
        show_row_names=FALSE, 
        row_split = cldf[rg,'dir'], 
        column_split = rmeta$nrad, 
        # column_split = rmeta$braaksc, 
        # column_split = rmeta$region, 
        show_column_names=FALSE)
dev.off()

# Correlation on this matrix:
cr = cor(t(pltmat))
mx = 0.25
rg = rownames(cr)
rgenes = c(rg[rg %in% upgene], rg[rg%in% downgene])
crmat = cr[rgenes,rgenes]
mx = 1
crmat[crmat > mx] = mx
crmat[crmat < -mx] = -mx
rows <- (seq(0,ncol(cr), 50) - .5) / (ncol(cr) -1)

png(paste0(imgpref, 'mastlmm_siggenes_regionlvl_cormap_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=300, units='in', width=12,height=11.5)
par(mar=c(.25,2,.25,.25))
image(crmat, col=rev(colrb), zlim=c(-mx,mx), useRaster=T, axes=F)
yat=seq(0,1,length.out=nrow(crmat))
text(x=parpos(1,0.001), yat, rownames(crmat), xpd=TRUE, adj=1, cex=.25)
box(lwd=.25)
abline(v=rows, lwd=.2)
abline(h=rows, lwd=.2)
dev.off()

# Cell vs. expr heatmap
smat = t(scale(t(norm)))
spathdf = pathdf[,c('barcode','nft','plaq_n')]
spathdf = spathdf[order(spathdf$plaq_n),]
spathdf = spathdf[order(spathdf$nft),]

png(paste0(imgpref, 'annotated_cell_heatmap.png'), res=400, units='in', width=20, height=8)
Heatmap(smat[cldf$gene,spathdf$barcode], name = "mat", 
        cluster_rows = FALSE, 
        cluster_columns=FALSE,
        right_annotation = ha,
        show_column_names=FALSE, show_row_names=FALSE)
dev.off()


ndf = as.data.frame(pltmat)
ndf$gene = rownames(ndf)
ndf = gather(ndf, rind, value, -gene)
ndf = merge(ndf, rmeta)
ndf = merge(ndf, cldf)

aggdf = aggregate(value ~ region + nrad + cls + dir, ndf, mean)
aggdf$nrad = factor(aggdf$nrad, levels=c('CTRL','AD'))

scols = snap.cols[1:max(ndf$cls) + 19]
gplot = ggplot(aggdf , aes(x=region, y=value, fill=factor(cls), alpha=nrad)) + 
# gplot = ggplot(ndf, aes(x=region, y=value, fill=factor(cls), alpha=nrad)) + 
    facet_wrap(dir ~cls, scales='free_y') + 
    geom_bar(stat='identity', position='dodge') + 
    scale_fill_manual(values=scols) + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_alpha_manual(values=c("CTRL" = 0.5, "AD" = 1)) + 
    labs(y='Average Expression', x='Region (light=CTRL / dark=AD)') + 
    theme_pubr() + theme(legend.position='none') 

ggsave(paste0(imgpref, 'mastlmm_graph_subcls_expr_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=7,height=5)


# ----------------------
# Expression by regions:
# ----------------------
# Aggregate averages by regions:
ndf = as.data.frame(norm)
ndf$gene = rownames(ndf)
ndf = gather(ndf, barcode, value, -gene)
ndf = merge(ndf, cldf)
ndf$region = sub("_.*","",ndf$barcode)
ndf$nrad = pathdf[ndf$barcode,'nrad']

aggdf = aggregate(value ~ region + nrad + cls, ndf, mean)
aggdf$nrad = factor(aggdf$nrad, levels=c('CTRL','AD'))

scols = snap.cols[1:max(ndf$cls) + 19]
gplot = ggplot(aggdf, aes(x=region, y=value, fill=factor(cls), alpha=nrad)) + 
    facet_wrap(~cls, scales='free_y') + 
    geom_bar(stat='identity', position='dodge') + 
    scale_fill_manual(values=scols) + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_alpha_manual(values=c("CTRL" = 0.5, "AD" = 1)) + 
    labs(y='Average Expression', x='Region (light=CTRL / dark=AD)') + 
    theme_pubr() + theme(legend.position='none') 
ggsave(paste0(imgpref, 'mastlmm_graph_subcls_expr_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=7,height=5)

# Plot as a heatmap:
ndf$projid = pathdf[ndf$barcode, 'projid']
aggdf = aggregate(value ~ region + nrad + cls + projid + gene, ndf, mean)

awide = spread(aggdf[,c('region','gene','projid','nrad', 'value')], gene, value)
awide = awide[order(awide$region),]
awide = awide[order(awide$nrad),]
amat = as.matrix(awide[,-c(1:3)])

rmat = sweep(amat, 2, colSums(amat) / 100, '/') 
# normalize (unclear..)
rmat = amat
rmat[rmat > 4] = 4
image(rmat, col=colv)

# ----------------------------------
# Plot with the imputed pseudotimes:
# ----------------------------------
impdir = 'multiRegion/rw_top_imputed_scores/'
ibc = scan(paste0(impdir, celltype, '_barcodes.txt'),'c')
ipath = 'nft'
isc = as.numeric(scan(paste0(impdir, celltype, '_imputed_scores_per_region_False_',ipath,'.txt'),'c'))
# Mic_Immune_imputed_scores_per_region_False_plaq_d.txt
# Mic_Immune_imputed_scores_per_region_False_plaq_n.txt
isc = scale(isc)
idf = data.frame(barcode=ibc, score=isc)
idf = idf[idf$barcode %in% colnames(norm),]
gc()

labgenes = c('MACF1','LINGO1','APPL2','ROBO1','MT1M','GFAP','XAF1',
             'SLC39A11','HMBOX1','AQP4','CD44','FOS','ERBIN', 'CRYAB',
             'HMGB1','HSP90AB1', 'S100B','MT1M','MT1E','MT1X','F3','MT-ND1')
# dgenes = c('GFAP', 'DPP10', 'CD44', 'GPC5', 'RGMA',
#     'DGKB', 'NAV2', 'DGKB', 'SLC39A11', 'CRYAB',
#     'HSP90AA1', 'HSPH1', 'FTH1', 'HMGB1', 'SPARCL1', 'CDH20'
labgenes = labgenes[labgenes %in% rownames(norm)]

# lat = sapply(labgenes, function(x){which(rownames(pltmat) == x)})

# Smooth over the pseudotimes:
pdf = c()

for (gene in labgenes){
    print(gene)
    idf$gene = norm[gene,idf$barcode]
    fit = loess(gene ~ score, idf)
    lims = range(idf$score)
    xat = seq(lims[1],lims[2], length.out=50)
    pred = predict(fit, xat)
    pdf = rbind(pdf, data.frame(gene=gene, x=xat, pred=pred))
}
pdf = merge(pdf, data.frame(x=xat, pos=1:length(xat)))
pdf = unique(pdf)

# Plot:
pwide = spread(pdf[,c('pos','pred','gene')], pos, pred)
pmat = as.matrix(pwide[,-1])
rownames(pmat) = pwide$gene
# rownames(cldf) = cldf$gene
pltmat = t(scale(t(pmat)))
pltmat = reord(pltmat, measure='cosine', method='ward.D')
pltmat = pltmat[order(apply(pltmat, 1, which.max), decreasing=T),]

# Gene annotation:
hg = rowAnnotation(Dir = cldf[rownames(pltmat),'dir'], 
                   Cls = cldf[rownames(pltmat),'cls'],
                   # Gene = anno_mark(at=1:nrow(pltmat), labels = rownames(pltmat)),
                   col = list(Cls=scols,
                              Dir=c('Down'='royalblue','Up'='indianred')))


png(paste0(imgpref, 'scaled_heatmap_siggenes_sel_loess_on_', ipath,'_for_', path, '.allreg.major.', celltype, '.minor.', ststr,'.png'), res=400, units='in', width=7, height=8)
Heatmap(pltmat, name = "Loess", 
        # cluster_rows = TRUE, 
        cluster_rows = FALSE, 
        cluster_columns=FALSE,
        # right_annotation = hg,
        show_row_names=TRUE, 
        show_column_names=FALSE)
dev.off()


gene = 'HSPA1A'
gene = 'CD81'
mean(norm[gene,pathdf$nrad == 'AD'], na.rm=T)
mean(norm[gene,pathdf$nrad == 'CTRL'], na.rm=T)

mean(mat[gene,pathdf$plaq_n > 5], na.rm=T)
mean(mat[gene,pathdf$plaq_n < 5], na.rm=T)

# ----------------------------------------
# Load in the pathology pseudotime scores:
# ----------------------------------------
# For EC - change for others:
pathdf = merge(pathdf, unique(metadata[,c('projid','nft_ec','plaq_d_ec','plaq_n_ec')]))
pathdf[pathdf$region == 'TH',c('nft','plaq_d','plaq_n')] = 0 # FOR PLOTTING ONLY
pathdf = merge(pathdf,  unique(metadata[,c('projid','braaksc')]))
rownames(pathdf) = pathdf$barcode
impdir = 'multiRegion/rw_top_imputed_scores/'
ptstr = paste0(path, '.allreg.major.', celltype, '.minor.', ststr)
for (ipath in c('nft','plaq_d','plaq_n')){
    print(ipath)
    ibc = scan(paste0(impdir, celltype, '_barcodes.txt'),'c')
    isc = as.numeric(scan(paste0(impdir, celltype, '_imputed_scores_per_region_False_',ipath,'.txt'),'c'))
    # Mic_Immune_imputed_scores_per_region_False_plaq_d.txt
    # Mic_Immune_imputed_scores_per_region_False_plaq_n.txt
    isc = scale(isc)
    idf = data.frame(barcode=ibc, score=isc)
    idf = idf[idf$barcode %in% colnames(norm),]
    gcout = gc()

    rownames(idf) = idf$barcode
    idf = idf[order(idf$score),]
    rns = idf$barcode
    idf$nrad = pathdf[rns,'nrad']

    snorm = norm[rownames(pltmat),rns]
    snorm = t(scale(t(snorm)))

    ptime = seq(min(idf$score), max(idf$score), length.out=100)
    plasma.col_fun = colorRamp2(ptime, plasma(100))

    hb = HeatmapAnnotation(AD = pathdf[rns,'nrad'],
                           Region = pathdf[rns,'region'],
                           Cogn = pathdf[rns,'cogdx'],
                           NFT = pathdf[rns,'nft'],
                           Braak = pathdf[rns,'braaksc'],
                           PlaqN = pathdf[rns,'plaq_n'],
                           PlaqD = pathdf[rns,'plaq_d'],
                           Pseudotime= anno_simple(idf$score, col = plasma.col_fun),
                           col = list(Region=reg.cols,
                                      Braak=colvals[['braaksc']],
                                      AD=c('CTRL'='royalblue','AD'='indianred')))

    mx = 2
    snorm[snorm > mx] = mx
    snorm[snorm < -mx] = -mx

    png(paste0(imgpref, 'scaled_cell_level_heatmap_siggenes_', ptstr,'_bypstime_',ipath,'.png'), res=400, units='in', width=15, height=8)
    Heatmap(snorm, name = "Scaled\nExpr.", 
            cluster_rows = FALSE, cluster_columns=FALSE,
            right_annotation = ha,
            use_raster=T,
            top_annotation = hb,
            show_row_names=FALSE, 
            row_split = cldf[rg,'dir'], 
            column_split = pathdf[rns,'nrad'], 
            show_column_names=FALSE)
    dev.off()

    # -------------------------------
    # Plot some top genes + GAM fits:
    # -------------------------------
    idf = cbind(idf, t(snorm[labgenes[1:min(c(12,length(labgenes)))],rns]))
    idlong = gather(idf, gene, value, -barcode,-score,-nrad)

    gplot = ggplot(idlong, aes_string('score', 'value', color='nrad')) + 
        facet_wrap(~gene, scales='free_y') + 
        geom_point(alpha=.25, pch='.') + 
        geom_smooth(method='gam') + 
        scale_color_manual(values=c('AD'='indianred','CTRL'='royalblue'), name='Status:') + 
        scale_y_continuous(expand=c(0,0)) + 
        scale_x_continuous(expand=c(0,0)) + 
        labs(x=paste(ipath, 'pseudotime (scaled)'),y='Normalized + Scaled Expression') + 
        theme_pubr() 
    ggsave(paste0(imgpref, 'mastlmm_examplegenes_',ptstr,'_bypstime_',ipath,'.png'), gplot, dpi=400, units='in', width=12,height=10)
    
    gcout = gc()

}








