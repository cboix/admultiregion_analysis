#!/usr/bin/R
# -----------------------------------------------------------
# Aggregate and plot chunked runs of MAST + RE for DGE
# Updated: 01/21/2021
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

celltype = 'Exc'
subtype = 'Exc_TOX3_TTC6'
region = 'EC'
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


# for (region in regions[regions !='MB']){
# sts = c('Inh_ALCAM_TRPM3', 'Inh_CUX2_MSR1', 'Inh_ENOX2_SPHKAP',
#         'Inh_FBN2_EPB41L4A', 'Inh_GPC5_RIT2', 'Inh_L3-5_SST_MAFB',
#         'Inh_L6_SST_NPY', 'Inh_LAMP5_RELN', 'Inh_PAX6_RELN',
#         'Inh_PTPRK_FAM19A1', 'Inh_PVALB_HTR4', 'Inh_RYR3_TSHZ2',
#         'Inh_SGCD_PDE3A', 'Inh_SORCS1_TTN', 'Inh_VIP_ABI3BP',
#         'Inh_VIP_CLSTN2', 'Inh_VIP_THSD7B', 'Inh_VIP_TSHZ2')
# sts = c('End','Fib','Per','SMC') # ,'CPEC')
# sts = c('T_cells','CAMs','Mic')
# for (subtype in sts){

# Load in data:
ststr = gsub("/","_",gsub(" ","_", subtype))
matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                 celltype,'.',ststr,'.',region)
rdafile = paste0(matpref, '.rda')  # In Matrix format
# Load `mat` from rdafile:
load(rdafile)
amat = mat
print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
barcodes = colnames(mat)
genes = rownames(mat)
ngenes = nrow(mat)

# ------------------------------
# Load the appropriate metadata:
# ------------------------------
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

rownames(cellmeta) = cellmeta$barcode
pathdf = cellmeta[barcodes,]
pathdf = merge(pathdf, metadata[,c('projid','rind','age_death','msex','pmi', 'Apoe_e4', 'cogdx', 'niareagansc')])
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode
pathdf$nrad = 'CTRL'
pathdf$nrad[pathdf$niareagansc %in% c(1,2)] = 'AD'
pathdf$nrad = factor(pathdf$nrad, levels=c('CTRL','AD'))

# ----------------------------
# Load the regression results:
# ----------------------------
print(paste("[STATUS] Loading regression on", subtype, 'in',region, 'on var', path, 'with ncell:', ncol(mat)))
fpref = paste0(prefix, '.mastlmm_reg.', path, '.', region, '.major.', celltype, '.minor.', ststr)
fnlist = list.files(pattern=paste0(fpref,'.*.Rda'), path=regdir)
resdf = c()
sumdf= c()
for (fn in fnlist){
    print(fn)
    load(paste0(regdir,fn))
    resdf = rbind(resdf, regdf)
    sumdf = rbind(sumdf, summaryDt)
}
dim(resdf)

Z = 1 * (mat  > 0)
gmarg = rowMeans(Z)
rm(Z); gcout = gc()
# hist(log10(gmarg), 500)
resdf$per = gmarg[resdf$primerid]
names(resdf)[2] = 'pvalue'
resdf = resdf[order(resdf$pvalue),]
resdf$padj = p.adjust(resdf$pvalue, 'fdr')

# resdf[resdf$primerid %in% c('NRCAM','NFASC','SCN8A','SCN9A','KCNQ3','CNTN1','CNTNAP2','ANK3','ADAM22','DLG1','DLG2'),] -> adf
# adf[order(adf$pvalue),c('primerid','coef','pvalue','padj','per')]


FCTHRESH=0.05
resdf$col = 1 * (resdf$padj < 0.05) * (1 + 1 * (resdf$coef > 0)) * (abs(resdf$coef) > FCTHRESH)
pltdf = resdf
mx = 0.5
pltdf$coef[pltdf$coef > mx] = mx
pltdf$coef[pltdf$coef < -mx] = -mx
ntop = sum(pltdf$padj < 0.05 & pltdf$col !=0, na.rm=T)
labcut = 0.05
if (ntop > 50){
    labcut = 0.02
}
labdf = pltdf[(pltdf$padj < labcut) & pltdf$col != 0,]

pcols = brewer.pal(12, 'Paired')
gplot = ggplot(pltdf[pltdf$per > .25,], aes(coef, -log10(pvalue), color=factor(col))) + 
    geom_point(cex=.25) + 
    geom_text_repel(data=labdf, aes(label=primerid), size=2, segment.size=.25, max.overlaps=20) + 
    theme_pubr() + 
    theme(legend.position='none') + 
    labs(x='Coefficient (AD / CTRL by NIA-Reagan score)', y='-log10(p-value)', title=paste0(gsub("_"," ",subtype), ' - Differential Genes')) + 
    scale_color_manual(values=c('grey85',pcols[1],pcols[5])) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_vline(xintercept=c(-FCTHRESH, FCTHRESH), lty='dashed', col='grey50',lwd=.25)  +
    geom_vline(xintercept=0, lty='dashed', lwd=.25)
# ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.', region, '.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=6,height=6)
# ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.', region, '.major.', celltype, '.minor.', ststr,'_small.png'), gplot, dpi=400, units='in', width=5,height=5)
ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.', region, '.major.', celltype, '.minor.', ststr,'_tiny.png'), gplot, dpi=400, units='in', width=2.5,height=3)
ggsave(paste0(imgpref, 'mastlmm_volcano_', path, '.', region, '.major.', celltype, '.minor.', ststr,'_tiny.pdf'), gplot, dpi=400, units='in', width=2.5,height=3)

# }}

# }
# }

# -------------------
# Run GO enrichments:
# -------------------
library(gprofiler2)
upgene = resdf[resdf$col == 2,'primerid']
downgene = resdf[resdf$col == 1,'primerid']
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
pdf$term_name = factor(pdf$term_name, levels=rev(pdf$term_name))
gplot = ggplot(udf, aes(intersection_size / query_size, -log10(p_value), size=term_size, fill=source)) + 
    geom_point(pch=21, alpha=.5) + 
    geom_text_repel(data=pdf, aes(label=term_name), size=2.5, point.padding=.25, segment.size=.25, box.padding=.1) +
    scale_x_continuous(labels=scales::percent) + 
    # scale_y_continuous(expand=c(0,0)) + 
    labs(x='Percent sig. genes in term', y='-log10(p-value)', title=paste0(gsub("_"," ",subtype), ' - Up-regulated Genes (All Regions)')) + 
    theme_pubr() + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'mastlmm_GObubble_up_', path, '.',region, '.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=5,height=3)


# Plot down-regulated:
pdf = head(ddf, 20)
pdf$term_name = factor(pdf$term_name, levels=rev(pdf$term_name))
gplot = ggplot(ddf, aes(intersection_size / query_size, -log10(p_value), size=term_size, fill=source)) + 
    geom_point(pch=21, alpha=.5) + 
    geom_text_repel(data=pdf, aes(label=term_name), size=2.5, point.padding=.25, segment.size=.25, box.padding=.1, max.iter=5000) +
    scale_x_continuous(labels=scales::percent) + 
    # scale_y_continuous(expand=c(0,0)) + 
    labs(x='Percent sig. genes in term', y='-log10(p-value)', title=paste0(gsub("_"," ",subtype), ' - Down-regulated Genes (All Regions)')) + 
    theme_pubr() + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'mastlmm_GObubble_down_', path, '.',region,'.major.', celltype, '.minor.', ststr,'.png'), gplot, dpi=400, units='in', width=5,height=3)




# --------------------------------------------------------------------------------------
# Plot the expression of the differential genes against each other (together/submodules)
# - Heatmaps (directly
# - Correlation of deg # Cond. corr at cell level + plots at cell level
# - Per-region
# - Per other path/etc. variables
# Ask if the DE gene expression is also significantly affected by other vars:
# --------------------------------------------------------------------------------------
agenes = c(upgene, downgene)
mat = as.matrix(amat[agenes,])
norm = sweep(mat, 2, marg[colnames(mat),'count'] / 100000,'/')
norm = log(norm + 1)

# Order genes:
rmat = reord(norm)

# Order cells:
cellscore = colSums(norm[upgene,]) - colSums(norm[downgene,])
pathdf$cs = cellscore
pathdf = pathdf[order(pathdf$cs),]
pathdf = pathdf[order(pathdf$nrad),]
# pathdf = pathdf[order(pathdf$region),]
sbc = pathdf$barcode
rmat = rmat[,sbc]
ptstr = paste0(path, '.', region, '.major.', celltype, '.minor.', ststr)
smat = t(scale(t(rmat)))
clamp.cent = function(mat, mx=2.5){
    mat[mat > mx] = mx
    mat[mat < -mx] = -mx
    return(mat)
}

smat = clamp.cent(smat)

colv = viridis(100)
png(paste0(imgpref, 'mastlmm_siggenes_heatmap_', ptstr ,'.png'), res=300, units='in', width=3,height=6)
par(mar=rep(0.2,4))
image(smat, col=rev(colrb), useRaster=T, axes=F)
box(lwd=.25)
dev.off()

# Correlation patterns:
cr = cor(t(rmat))
mx = 0.25
rg = rownames(cr)
rgenes = c(rg[rg %in% upgene], rg[rg%in% downgene])
crmat = cr[rgenes,rgenes]
crmat[crmat > mx] = mx
crmat[crmat < -mx] = -mx
rows <- (seq(0,ncol(cr), 50) - .5) / (ncol(cr) -1)

png(paste0(imgpref, 'mastlmm_siggenes_cormap_', ptstr, '.png'), res=300, units='in', width=12,height=11.5)
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
cutoff = 0.16
dmat = crmat 
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

# sum(ctodf[as.character(ddf$T1), 'COLOR'] ==
#     ctodf[as.character(ddf$T2), 'COLOR']) /nrow(ddf)

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
E(net)$width = sdf$sim * 3
# E(net)$weight = sdf$sim * .5
E(net)$weight = sdf$sim * 1
# E(net)$weight = sdf$sim * .5
# set.seed(2)
set.seed(8)
l <- layout_with_fr(net, grid='nogrid') # Usually best

png(paste0(imgpref, 'mastlmm_siggenes_graph_', ptstr ,'.png'), res=300, units='in', width=7,height=7)
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

png(paste0(imgpref, 'mastlmm_siggenes_graph_repel_', ptstr ,'.png'), res=300, units='in', width=7,height=7)
sp = 0.25
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# Repel points:
lbcex=0.6
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
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

png(paste0(imgpref, 'mastlmm_siggenes_graph_repel_louvain_', ptstr ,'.png'), res=300, units='in', width=7,height=7)
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



# ------------------------
# Averages by individuals:
# ------------------------
# rm(amat)
# gcout = gc()
rownames(pathdf) = pathdf$barcode
pathdf = pathdf[colnames(norm),]
tform = make.tform(pathdf$rind, norm=T)

# Average normalized matrix
norm.avg = norm %*% tform
mat.avg = mat %*% tform
rmat.avg = reord(norm.avg)
rmat.avg = t(reord(t(rmat.avg)))

# Scaled for relative expr:
snorm = t(scale(t(norm)))
snorm.avg = snorm %*% tform

smat.avg = scale(t(rmat.avg))
pltmat = snorm.avg
pltmat = t(smat.avg)
pltmat = t(scale(t(mat.avg)))

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
rmeta = metadata[colnames(pltmat),]
rmeta$nrad = ifelse(rmeta$niareagansc > 2,"CTRL","AD")
mx = 2.5
pltmat[pltmat > mx] = mx
pltmat[pltmat < -mx] = -mx

# TODO: Diff genes by celltype
labgenes = c('PRNP', 'MAP1A', 'HSP90B1', 'HSP90AB1', 'OLFM1', 'USF2', 
             'HSPA13', 'NPTX1', 'UBB', 'UBA52', 'SV2A', 'BSG', 'PSAP',
             'NCDN', 'APLP1', 'CPNE4', 'ARHGAP22', 'LRRK1', 'MYRIP', 
             'RALGPS1', 'BMPER', 'NRXN3', 'YWHAH')

# labgenes = c('HLA-DRB1','APOC1','C1QC','CD74','IFI44L','MX1','APOE',
#              'UBC','HSPA1A','HIF1A','TLR2','GLDN','SYNDIG1','GLDN',
#              'IL15','IL4R','HDAC4','PTEN','DNAJA4','APPL2','DHFR',
#              'NPL','TMSB4X','RPL19','BIN1','C3','CSF1R','CX3CR1','INPP5D','RASGEF1C',
#              'LINGO1','MT-ND3','MT-CO1','HSP90B1')
labgenes = labgenes[labgenes %in% rownames(pltmat)]
lat = sapply(labgenes, function(x){which(rownames(pltmat) == x)})

# Gene annotation:
ha = rowAnnotation(Dir = cldf[,'dir'], 
                   Cls = cldf[,'cls'],
                   Gene = anno_mark(at=lat, labels = labgenes),
                   col = list(
                              Cls=scols,
                              Dir=c('Down'='royalblue','Up'='indianred')))
# Region annotation:
hb = HeatmapAnnotation(AD = rmeta[,'nrad'],
                       Braak = rmeta[,'braaksc'],
                       Cogn = rmeta[,'cogdx'],
                       # TODO: Diff by region:
                       PlaqN.EC = rmeta[,'plaq_n_ec'],
                       PlaqD.EC = rmeta[,'plaq_d_ec'],
                       NFT.EC = rmeta[,'nft_ec'],
                       col = list(Region=reg.cols,
                                  Braak=colvals[['braaksc']],
                                  AD=c('CTRL'='royalblue','AD'='indianred')))


png(paste0(imgpref, 'scaled_heatmap_siggenes_split_', path, '.', region,'.major.', celltype, '.minor.', ststr,'.png'), res=400, units='in', width=7, height=8)
Heatmap(pltmat, name = "Scaled\nExpr.", 
        cluster_rows = TRUE, cluster_columns=TRUE,
        right_annotation = ha,
        use_raster=T,
        top_annotation = hb,
        show_row_names=FALSE, 
        row_split = cldf$dir, 
        column_split = rmeta$nrad, 
        # column_split = rmeta$braaksc, 
        # column_split = rmeta$region, 
        show_column_names=FALSE)
dev.off()

pdf(paste0(imgpref, 'scaled_heatmap_siggenes_split_braak_', path, '.', region,'.major.', celltype, '.minor.', ststr,'.pdf'), width=7, height=8)
Heatmap(pltmat, name = "Scaled\nExpr.", 
        cluster_rows = TRUE, cluster_columns=TRUE,
        right_annotation = ha,
        use_raster=T,
        top_annotation = hb,
        show_row_names=FALSE, 
        row_split = cldf$dir, 
        # column_split = rmeta$nrad, 
        column_split = rmeta$braaksc, 
        # column_split = rmeta$region, 
        show_column_names=FALSE)
dev.off()

# Also play with computing raw (uncorrected) log2FC
# snorm = norm[rownames(pltmat),]
# snorm = log(mat[rownames(pltmat),] + 1)
snorm = norm[rownames(pltmat),]
snorm = t(scale(t(snorm)))

hb = HeatmapAnnotation(AD = pathdf[,'nrad'],
                       Cogn = pathdf[,'cogdx'],
                       col = list(Region=reg.cols,
                                  Braak=colvals[['braaksc']],
                                  AD=c('CTRL'='royalblue','AD'='indianred')))

png(paste0(imgpref, 'scaled_cell_level_heatmap_siggenes_', path, '.', region,'.major.', celltype, '.minor.', ststr,'.png'), res=400, units='in', width=7, height=8)
Heatmap(snorm, name = "Scaled\nExpr.", 
        cluster_rows = TRUE, cluster_columns=TRUE,
        right_annotation = ha,
        use_raster=T,
        top_annotation = hb,
        show_row_names=FALSE, 
        row_split = cldf$dir, 
        column_split = pathdf$nrad, 
        # column_split = rmeta$braaksc, 
        # column_split = rmeta$region, 
        show_column_names=FALSE)
dev.off()

mean(mat['PRNP',pathdf$nrad == 'AD'], na.rm=T)
mean(mat['PRNP',pathdf$nrad == 'CTRL'], na.rm=T)

table(pathdf[,c('projid','nrad')])

# ----------------------------------------
# Load in the pathology pseudotime scores:
# ----------------------------------------
# For EC - change for others:
pathdf = merge(pathdf, unique(metadata[,c('projid','nft_ec','plaq_d_ec','plaq_n_ec')]))
rownames(pathdf) = pathdf$barcode
impdir = 'multiRegion/rw_top_imputed_scores/'
for (ipath in c('nft','plaq_d','plaq_n')){
    print(ipath)
    ibc = scan(paste0(impdir, celltype, '_barcodes.txt'),'c')
    isc = as.numeric(scan(paste0(impdir, celltype, '_imputed_scores_per_region_False_',ipath,'.txt'),'c'))
    # Mic_Immune_imputed_scores_per_region_False_plaq_d.txt
    # Mic_Immune_imputed_scores_per_region_False_plaq_n.txt
    isc = scale(isc)
    idf = data.frame(barcode=ibc, score=isc)
    idf = idf[idf$barcode %in% colnames(norm),]
    gc()

    rownames(idf) = idf$barcode
    idf = idf[order(idf$score),]
    rns = idf$barcode
    idf$nrad = pathdf[rns,'nrad']

    snorm = norm[rownames(pltmat),rns]
    snorm = t(scale(t(snorm)))

    ptime = seq(min(idf$score), max(idf$score), length.out=100)
    plasma.col_fun = colorRamp2(ptime, plasma(100))

    hb = HeatmapAnnotation(AD = pathdf[rns,'nrad'],
                           Cogn = pathdf[rns,'cogdx'],
                           nft.ec = pathdf[rns,'nft_ec'],
                           plaqn.ec = pathdf[rns,'plaq_n_ec'],
                           plaqd.ec = pathdf[rns,'plaq_d_ec'],
                           Pseudotime= anno_simple(idf$score, col = plasma.col_fun),
                           col = list(Region=reg.cols,
                                      Braak=colvals[['braaksc']],
                                      AD=c('CTRL'='royalblue','AD'='indianred')))

    png(paste0(imgpref, 'scaled_cell_level_heatmap_siggenes_', ptstr,'_bypstime_',ipath,'.png'), res=400, units='in', width=15, height=8)
    Heatmap(snorm, name = "Scaled\nExpr.", 
            cluster_rows = FALSE, cluster_columns=FALSE,
            right_annotation = ha,
            use_raster=T,
            top_annotation = hb,
            show_row_names=FALSE, 
            row_split = cldf$dir, 
            column_split = pathdf[rns,'nrad'], 
            show_column_names=FALSE)
    dev.off()

    # -------------------------------
    # Plot some top genes + GAM fits:
    # -------------------------------
    idf = cbind(idf, t(snorm[labgenes,rns]))
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

}




