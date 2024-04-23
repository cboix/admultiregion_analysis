#!/usr/bin/R
# ----------------------------------------------------------
# Plot the average gene expression for markers for clusters:
# Script used (pre-Sept 2020) for annotating clusters and removing doublets
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)

#' Calculate position relative to par()$usr 
#'
#' @param axis 1 is x; 2 is y;
#' @param shift percent of axis to shift over
#' @return position shifted from start of x or y axis
#' @export
parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}


# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/markers/')
imgpref = plotdir
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# ----------------------------------
# Read the umap coords and barcodes:
# ----------------------------------
datadir = 'multiRegion/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy_norm'
bcfile = paste0(datadir, prefix, '.barcodes.tsv.gz')
barcodes = scan(gzfile(bcfile), 'c')
# bcfile = paste0(datadir,'all_brain_regions_filt_preprocessed_barcodes.txt')
# prefix = 'all_brain_regions_filt_preprocessed_norm_8'
# barcodes = scan(bcfile, 'c')
umapfile = paste0(datadir, prefix, '.umap.tsv.gz')

# Groups to plot:
# lblset = 'groups_hdb'
lblset = 'groups_hdb_lg'
lblset = 'groups_r10'
lblset = 'groups_r5'
# lblset = 'groups_r3'
# lblset = 'groups_r1'
# lblset = 'groups_r2'
# lblset = 'groups_r0.5'
lblset = 'leiden_r5_n50'
db.bc.old = read.delim(paste0(datadir, prefix, '.', lblset, '.dblt.tsv'), header=F)
names(db.bc.old) = 'barcode'
lblset = 'leiden_r15_n100'
lblfile = paste0(datadir, prefix, '.', lblset, '.tsv.gz')
lblavgfile = paste0(datadir, prefix, '.', lblset, '.avg.tsv.gz')

# Read data:
celldf = read.delim(gzfile(umapfile), header=F, sep="\t")
names(celldf) = c('U1','U2')
celldf$barcode = barcodes
celldf$rind = sub("_.*-",".",barcodes)
celldf$region = sub("_.*","",barcodes)
celldf$projid = metadata[celldf$rind,'projid']
celldf$lbl = as.numeric(scan(gzfile(lblfile), 'c')) + 1 # 1 indexing
if (length(grep("hdb", lblset)) > 0){ celldf$lbl = celldf$lbl + 1} # HDB has a -1 category
lblcountdf = aggregate(rind ~ lbl, celldf, length)
lblregcountdf = aggregate(rind~ lbl + region, celldf, length)

# Read gene data:
gfile = paste0(datadir,prefix, '_genes.txt')
genedf = read.delim(gzfile(gfile), header=F)
# names(genedf) = c('symbol','gene')
names(genedf) = c('symbol')
gwide = read.delim(lblavgfile, header=F)
colnames(gwide) = genedf$symbol
gwide = gwide[,which(colSums(gwide) != 0)]

sum(celldf$lbl == 1) # Unassigned.
hlvls = as.character(sort(unique(celldf$lbl)))
lbl.cols = rep(snap.cols,4)[1:length(hlvls)]
names(lbl.cols) = as.character(hlvls)

if (length(grep('hdb', lblset)) > 0){
    ctype = 'HDBSCAN'
} else {
    ctype = 'Leiden'
}

# ------------------------
# Plot the cluster groups:
# ------------------------
NCELL = nrow(celldf)
ind = sample(1:NCELL,NCELL, replace=FALSE)
lbl.loc = aggregate(cbind(U1, U2) ~ lbl, celldf, mean)
cex = 0.025
tsp.lbl.cols = sapply(lbl.cols, tsp.col)

png(paste0(imgpref, 'umap_', lblset, '.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=tsp.lbl.cols[celldf$lbl[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, lbl, cex=.5, xpd=TRUE))
mtext(paste(ctype, 'Clusters'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()

# Read in Hans doublets + compare:
dbrdafile = paste0(datadir, 'predicted_doublets/allregions_predicted_doublets.Rda')
if (!file.exists(dbrdafile)){
    dbdf = c()
    for (region in regions){
        dbfile = paste0(datadir, 'predicted_doublets/', region, '_predicted_doublets.csv')
        rdbdf = read.delim(dbfile, header=T, sep=',', stringsAsFactors=F)
        rdbdf$barcode = paste0(region, '_', rdbdf$index)
        dbdf = rbind(dbdf, rdbdf[,c('barcode','predicted_doublets')])
    }
    flagged = barcodes[!(barcodes %in% dbdf$barcode)]
    dbdf = dbdf[dbdf$predicted_doublets == 'True',]
    save(dbdf, flagged, file=dbrdafile)
} else {
    load(dbrdafile)
}

# Load in Na's predicted doublets
na.db = scan(paste0(datadir, 'na_doublets/na_doublets.tsv'), 'c')

celldf$is.doublet = 1 * (celldf$barcode %in% dbdf$barcode)
celldf$is.doublet[celldf$barcode %in% flagged] = 2
celldf$is.doublet[celldf$barcode %in% na.db] = 3
celldf$is.doublet[celldf$barcode %in% db.bc.old$barcode] = 4

dbcount = aggregate(is.doublet ~ lbl, celldf, function(x){sum(x > 0)})
dbcount = merge(lblcountdf, dbcount)
dbcount$frac = dbcount$is.doublet / dbcount$rind
dbcount = dbcount[order(dbcount$frac, decreasing=T),]

db.cols = rep('grey50',length(hlvls))
names(db.cols) = as.character(hlvls)
db.cols[dbcount$lbl[dbcount$frac > .5]] = 'yellow'

png(paste0(imgpref, 'umap_', lblset, '_w_doublets.png'), units='in', res=450, width=16, height=8)
par(xaxs='i')
par(yaxs='i')
layout(matrix(1:2, nrow=1))
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], 
     col=db.cols[celldf$lbl[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
mtext(paste('Auto labels from', ctype, 'clusters'), side=3,
      cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], 
     col=c('black','yellow','red', 'blue', 'purple')[celldf$is.doublet[ind]+1], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
mtext(paste("HM Doublet calls"), side=3,
      cex=1.5, col='grey25', font=2, line=0.25)
dev.off()



# -------------------------------------------
# Read in basic markers and plot them:
# NOTE: SHOULD COMPUTE AVG WITHIN EACH REGION
# -------------------------------------------
# Overall:
m1df = read.delim('Annotation/broad_ct_markers.tsv', header=T, stringsAsFactors=F)
m2df = read.delim('Annotation/known_markers.tsv', header=T, stringsAsFactors=F)
# Specific:
emdf = read.delim('Annotation/excitatory_markers.tsv', header=F, stringsAsFactors=F)
imdf = read.delim('Annotation/inhibitory_markers.tsv', header=F, stringsAsFactors=F)
# 

# Colors:
cell.map = sapply(names(celltype.col), function(x){
                      if (x %in% c('Inhibitory', 'Excitatory')){ substr(x, 1,2) } else { substr(x, 1,3) }})
cdf = data.frame(col=celltype.col, cell=cell.map[names(celltype.col)])
ct.col = as.character(cdf$col)
names(ct.col) = cdf$cell

# Subset to these genes:
mkdf = data.frame(symbol=c(m1df$gene, m2df$symbol, emdf[,1], imdf[,1], c('AMBP', 'CSPG4', 'PDGFB')),
                  cell=c(m1df$cell, m2df$celltype, rep('Ex', nrow(emdf)), rep('In', nrow(imdf)), rep('Per',3)))
mkdf$symbol = as.character(mkdf$symbol)
mkdf$cell = as.character(mkdf$cell)
mkdf$cell[mkdf$cell == 'Exc'] = 'Ex'
ind = mkdf$cell %in% names(cell.map)
mkdf$cell[ind] = cell.map[mkdf$cell[ind]]
mkdf = unique(mkdf)
mkdf = merge(mkdf, cdf)

# Gene color map:
gcolmap = aggregate(col ~ symbol, mkdf, function(x){ x[1] })
rownames(gcolmap) = gcolmap$symbol
gcolmap$col = as.character(gcolmap$col)

# TODO: Need more markers, OFC (pericyte etc)
markergenes = unique(mkdf$symbol)
mat = as.matrix(gwide[, markergenes])
rownames(mat) = as.character(1:nrow(mat))

reord2d = function(mat, ncol, nrow, measure='euclidean', method='ward.D'){
    require(cba)
    # First, cluster all groups based on their markers:
    dt <- dist(mat, measure)
    ht <- hclust(dt, method = method)
    cocl <- order.optimal(dt, ht$merge)$order
    reord <- names(cocl)[cocl]
    mat = mat[reord, ]
    cols = calc.breaks(ht, nclust=ncol, cocl=cocl)
    # Reorder rows:
    dt <- dist(t(mat), measure)
    ht <- hclust(dt, method = method)
    cocl <- order.optimal(dt, ht$merge)$order
    reord <- names(cocl)[cocl]
    mat = mat[, reord]
    rows = calc.breaks(ht, nclust=nrow, cocl=cocl)
    return(list(mat, rows, cols))
}

ll = reord2d(mat, 20, 20)
mat = ll[[1]]
rows = ll[[2]]
cols = ll[[3]]
# Auto color:

# palette = viridis(100)

w = 10 * nrow(mat) / 130
png(paste0(imgpref, 'avg_nonPEC_markers_groups_', lblset, '.png'), units='in', res=450, width=w, height=6)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
image(mat, col=colb, axes=F)
text(x=seq(0,1, length.out=nrow(mat)), y=parpos(2,.01), 
     rownames(mat), xpd=TRUE, cex=cex, srt=90, adj=1)
text(y=seq(0,1, length.out=ncol(mat)), x=parpos(1,.005), 
     colnames(mat), xpd=TRUE, cex=cex, adj=1, col=gcolmap[colnames(mat),'col'])
abline(h=rows, lwd=.25)
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()

# Add bars on top with # of cells from each.


# ------------------------------------
# Do avg with better markers from PEC:
# ------------------------------------
p1df = read.delim('Annotation/PEC_markers_1.tsv', header=T, stringsAsFactors=F)
p2df = read.delim('Annotation/PEC_markers_2.tsv', header=T, stringsAsFactors=F)
# Allen brain markers:
adf = read.delim(paste0('Annotation/marker_genes_Allen_Brain_Institute.csv'), header=T, sep=",", stringsAsFactors=F)
adf = unique(adf[,c('gene','cluster')])
names(adf) = c('Gene','Cluster')
# Read in joint markers:
mlist = readRDS(paste0(datadir, 'sz_joint_markers.RDS'))$marker.genes
szdf = sapply(names(mlist), function(x){cbind(Gene=mlist[[x]], Cluster=x)})
szdf = do.call(rbind, szdf)
szdf = data.frame(szdf)
# Add to PEC:
p2df = rbind(rbind(p2df, szdf), adf)
p2df = p2df[p2df$Gene %in% colnames(gwide),]
types = unique(p2df$Cluster)
tform = make.tform(p2df$Cluster, norm=T,u=types)

markmat = as.matrix(gwide[, p2df$Gene])
# sort(apply(markmat, 2, max), decreasing=T)
markmat = scale(markmat)
avg.mat = markmat %*% tform
rownames(avg.mat) = as.character(1:nrow(avg.mat))
avg.mat = scale(avg.mat)

ll = reord2d(avg.mat, 30, 25)
avg.mat = ll[[1]]
rows = ll[[2]]
cols = ll[[3]]

fnMatSqrtInverse = function(mA) {
  ei = eigen(mA)
  d = ei$values
      d = (d+abs(d))/2
      d2 = 1/sqrt(d)
      d2[d == 0] = 0
      return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}

w = nrow(avg.mat) / 15
h = ncol(avg.mat) / 15
png(paste0(imgpref, 'avg_all_markers_groups_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,5,0,0))
cex=.4
thresh = 3
plt.mat = avg.mat
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.005), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.0025), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(h=rows, lwd=.25)
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()

# Save matrix for later (plot without doublets)
PEC.avg.mat = avg.mat

prec = fnMatSqrtInverse(cov(avg.mat))
dec.avg.mat = avg.mat %*% prec
colnames(dec.avg.mat) = colnames(avg.mat)
# dec.avg.mat[dec.avg.mat > 3] = 3
tmarg = apply(dec.avg.mat > 1, 1, sum) 

png(paste0(imgpref, 'avg_all_markers_groups_dec_', lblset, '.png'), units='in', res=450, width=w, height=2.5)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
plt.mat = dec.avg.mat
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.01), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1,
     col=ifelse(tmarg > 0, ifelse(tmarg > 1, 'red','black'), 'grey70'))
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()


# ----------------------------------------------------------
# For labeling doublets, collapse inhibitory and excitatory:
# ----------------------------------------------------------
am.ex = apply(avg.mat[,grep('^Ex',colnames(avg.mat))], 1, mean)
am.in = apply(avg.mat[,grep('^In',colnames(avg.mat))], 1, mean)

avg.mat2 = cbind(avg.mat[,c('Oligo', 'OPC', 'Astro', 'Per' ,'Endo', 'Microglia')],
                 'In'=am.in, 'Ex'=am.ex)

tmarg = apply(avg.mat2 > 1, 1, sum) 
# Needs to be relative to how much they're usually joined...

png(paste0(imgpref, 'avg_all_markers_groups_reduced_', lblset, '.png'), units='in', res=450, width=w, height=1)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
plt.mat = avg.mat2
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.01), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1,
     col=ifelse(tmarg > 0, ifelse(tmarg > 1, 'red','black'), 'grey70'))
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()

prec = fnMatSqrtInverse(cov(avg.mat2))
dec.avg.mat = avg.mat2 %*% prec
colnames(dec.avg.mat) = colnames(avg.mat2)
# dec.avg.mat[dec.avg.mat > 3] = 3
tmarg = apply(dec.avg.mat > .5, 1, sum) 

png(paste0(imgpref, 'avg_all_markers_groups_reduced_dec_', lblset, '.png'), units='in', res=450, width=w, height=1)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
plt.mat = dec.avg.mat
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.01), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1,
     col=ifelse(tmarg > 0, ifelse(tmarg > 1, 'red','black'), 'grey70'))
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()



# -----------------------------------------------------------
# Auto assign and label (ignoring the doublet issue for now):
# -----------------------------------------------------------
celldf$celltype = NULL
celldf$col = NULL
celldf$autocls = NULL
autodf = data.frame(lbl = rownames(avg.mat), 
                    autocls = colnames(avg.mat)[apply(avg.mat, 1, which.max)])
autodf$celltype = sapply(autodf$autocls, function(x){ 
                             if ((length(grep('Ex',x)) > 0) ||(length(grep('In',x)) > 0)){
                                 substr(x, 1,2) } else { substr(x, 1,3) }})
t.col = sapply(ct.col, tsp.col)
autodf$col = ct.col[autodf$celltype]
autodf$tspcol = t.col[autodf$celltype]
celldf = merge(celldf, autodf)
celldf$autocls = as.character(celldf$autocls)

id.double = c()
if (lblset == 'groups_hdb_lg'){
    id.double = c('20', '23', '37', '44', '36', '38', '47', '67', '66', 
                  '38', '27', '71', '68', '72', '14', '61', '49', '50', 
                  '60', '9', '16', '42','43', '88','100','32','33', '39','40', 
                  # after check doublets with hans:
                  '131','48','58','70', '137', '144','113')
    celldf$col[celldf$lbl == '1'] = 'grey75'
} else if (lblset == 'leiden_r5_n50'){
    id.double = c('155', '98', '64', '129', '116', '94', '117', '134', '149', '147', 
                  '47', '151', '157', '154', '125', '163', '113', '66', '95', '114', 
                  '118', '131', '120', '141', '140', '126', '158','107', '53')
} else if (lblset == 'leiden_r15_n100'){
    id.double = c(dbcount$lbl[dbcount$frac > .5], 
                  '229', '306', '309', '310', '311', '342', 
                  '326', '331',
                  '328', '127', # Tricky
                  '253', '251', '332', '261', '345', '277',
                  '111', '338','246',
                  '299', '283','278', # Doublets with mixed neurons, endo per
                  '234','242', # In oligos
                  '200', '285','308', '236') # From endo, per
}

# lbl = '248'
# a = lbl.cols[c('291','202', lbl)]
# plot(1:length(a), col=a, cex=10, pch=19)
# autodf[autodf$lbl == lbl,]
# lbl.loc[lbl.loc$lbl == lbl,]
# dbcount[dbcount$lbl == lbl,]
# # TODO: Check markers for 200
# gen = gwide['200',]
# gen = gwide['184',]
# cat(rev(names(tail(unlist(sort(gen)), 50))))

# hm double contains old double too...
hm.double = dbcount$lbl[dbcount$frac > .5]
id = which(celldf$lbl %in% c(hm.double, id.double))
celldf$col[id] = 'yellow'
celldf$tspcol[id] = tsp.col('yellow')
celldf$autocls[id] = 'Doublet'
# Check if any old id doublets left:
notdb = celldf[(celldf$barcode %in% db.bc.old) & (celldf$autocls != 'Doublet'),]
print(dim(notdb))

# Enumeration of doublets:
sum(celldf$autocls == 'Doublet')
sum(celldf$is.doublet > 0)
sum((celldf$autocls == 'Doublet') * (celldf$is.doublet > 0))
sum((celldf$autocls != 'Doublet') * (celldf$is.doublet > 0))
 
# To be conservative, add others to doublets as well:
celldf$autocls[celldf$is.doublet > 0] = 'Doublet'

# Write current doublets:
db.bc = celldf$barcode[celldf$autocls == 'Doublet']
write.table(db.bc, paste0(datadir, prefix, '.', lblset, '.dblt.tsv'), 
            quote=F, row.names=F, sep="\t", col.names=F)

# Write out the labels:
write.table(celldf[,c('barcode','autocls')], paste0(datadir, prefix, '.', lblset, '.lbls.tsv'), 
            quote=F, row.names=F, sep="\t", col.names=F)

# Shorten to label: 
cind = c(grep('^In', celldf$autocls), grep('^Ex', celldf$autocls))
celldf$autocls.short = celldf$autocls
celldf$autocls.short[cind] = sub("_.*$","", sub("^[A-Za-z]*_L","L", celldf$autocls[cind]))

# Write full datasets down:
write.table(celldf, paste0(datadir, prefix, '.', lblset, '.ext.lbl.tsv'), 
            quote=F, row.names=F, sep="\t", col.names=F)
save(celldf, file=paste0(datadir, prefix, '.', lblset, '.ext.lbl.Rda'))


lbl.loc = aggregate(cbind(U1, U2) ~ lbl + autocls, celldf, mean)

NCELL = nrow(celldf)
ind = sample(1:NCELL,NCELL, replace=FALSE)

cex = 0.025
png(paste0(imgpref, 'umap_groups_', lblset, '_autolbl.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], 
     col=celldf$tspcol[ind], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, autocls, cex=.5, xpd=TRUE))
mtext(paste('Auto labels from', ctype, 'clusters'), side=3,
      cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()

db.src.cols = c('black','yellow','red','blue','purple')

png(paste0(imgpref, 'umap_groups_', lblset, '_autolbl_w_doublets.png'), units='in', res=450, width=16, height=8)
par(xaxs='i')
par(yaxs='i')
layout(matrix(1:2, nrow=1))
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], 
     col=celldf$tspcol[ind], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, autocls, cex=.5, xpd=TRUE))
mtext(paste('Auto labels from', ctype, 'clusters'), side=3,
      cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], 
     col=db.src.cols[celldf$is.doublet[ind]+1], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
mtext(paste("HM Doublet calls"), side=3,
      cex=1.5, col='grey25', font=2, line=0.25)
dev.off()


# --------------------------------------------
# Plot the points without the doublet clusters
# --------------------------------------------
ind = which(celldf$autocls != 'Doublet')
ind = sample(ind, length(ind), replace=FALSE)
lbl.loc = aggregate(cbind(U1, U2) ~ lbl, celldf[ind,], mean)
cex = 0.025
tsp.lbl.cols = sapply(lbl.cols, tsp.col)

png(paste0(imgpref, 'umap_', lblset, '_wout_doublets.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=tsp.lbl.cols[celldf$lbl[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, lbl, cex=.5, xpd=TRUE))
mtext(paste(ctype, 'Clusters'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()


lbl.loc = aggregate(cbind(U1, U2) ~ lbl + autocls, celldf[ind,], mean)
# Shorter labels:
cind = c(grep('^In', lbl.loc$autocls), grep('^Ex', lbl.loc$autocls))
lbl.loc$autocls[cind] = sub("_.*$","", sub("^[A-Za-z]*_L","L", lbl.loc$autocls[cind]))
lbl.loc$autocls[-cind] = ''

cex = 0.025
png(paste0(imgpref, 'umap_groups_', lblset, '_autolbl_wout_doublets.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], 
     col=celldf$tspcol[ind], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, autocls, cex=.5, xpd=TRUE))
mtext(paste('Auto labels from', ctype, 'clusters'), side=3,
      cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()


# TODO: Find the bona-fide genes for each group - iterative assign.
# TODO: Diff expr to find top genes for each group...
# - start with means.

# ----------------------------
# Fraction of cells by region:
# ----------------------------
ct.cmap = unique(celldf[celldf$autocls != 'Doublet' & celldf$col != 'grey75', c('autocls', 'col')])
subct.col = ct.cmap$col
names(subct.col) = ct.cmap$autocls

# TODO: only for hdbscan:
countdf = aggregate(lbl ~ region + autocls, celldf[celldf$autocls != 'Doublet',], length)
tot.countdf = aggregate(lbl ~ region, countdf, sum)
names(tot.countdf) = c('region', 'total')
countdf =merge(countdf, tot.countdf)

# Cell fractions per region:
countdf = aggregate(lbl ~ region + autocls, celldf[celldf$autocls != 'Doublet',], length)
cwide = spread(countdf[,c('region','autocls','lbl')], autocls, lbl, fill=0)
cmat = as.matrix(cwide[,-1])
rownames(cmat) = cwide[,1]

au.coldf = unique(celldf[,c('autocls','col')])
au.col = au.coldf$col
names(au.col) = as.character(au.coldf$autocls)

# stackfill.barplot = function(mat, horiz=FALSE){
cmat = sweep(cmat, 1, apply(cmat, 1, sum), '/')
mmat = t(rbind(0, apply(cmat, 1, cumsum)))
NR = nrow(mmat)
NC = ncol(mmat)
xat = seq(0,1,.2)
mmat = mmat[rev(rownames(mmat)),]

w = 5
h = 2
png(paste0(imgpref, 'autocls_fractions_per_region_', lblset, '.png'), units='in', res=450, width=w, height=h)
# Locations:
layout(matrix(1:2, ncol=1), heights=c(.5,5))
sp = 0.1; lsp=2.5
par(xaxs='i')
par(yaxs='i')
par(mar=c(0,lsp, sp, sp))
plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0.5, NR + .5))
red.ct.col = ct.col[c(1,3:9)]
legend('center', legend=names(red.ct.col), col=red.ct.col, 
       pt.cex=1.5, cex=.75, pch=15, bty='n', ncol=8)
par(mar=c(1.5,lsp, sp, sp))
plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0.5, NR + .5))
pad = 0.45
for (i in 1:NR){
    region = rownames(mmat)[i]
    rect(xleft=mmat[i,-NC], xright=mmat[i,-1],
         ybottom=i -pad, ytop=i + pad, col=au.col[colnames(mmat)[-1]],
         border=NA, xpd=TRUE)
    rect(xleft=parpos(1,.1), xright=parpos(1,.01),
         ybottom=i -pad, ytop=i + pad, col=reg.cols[region],
         border=NA, xpd=TRUE)
    text(x=parpos(1,.02), y=i, region, adj=1,  srt=0, xpd=TRUE)
}
axis(1, at=xat, labels=FALSE, lwd=.5)
text(x=xat, y=parpos(2,.12), labels=xat, xpd=TRUE, cex=.8)
dev.off()



# Cell fractions by region, then by individual:
countdf = aggregate(lbl ~ region + rind + autocls, celldf[celldf$autocls != 'Doublet',], length)
cwide = spread(countdf[,c('rind','autocls','lbl')], autocls, lbl, fill=0)
cmat = as.matrix(cwide[,-1])
rownames(cmat) = cwide[,1]

cmat = sweep(cmat, 1, apply(cmat, 1, sum), '/')
mmat = t(rbind(0, apply(cmat, 1, cumsum)))
NR = nrow(mmat)
NC = ncol(mmat)
xat = seq(0,1,.2)
mmat = mmat[rev(rownames(mmat)),]
mlist = lapply(regions, function(x){ mmat[grep(x, rownames(mmat)),]})
heights = c(.5, sapply(mlist, nrow) / 24, .75)

w = 5
h = 3
png(paste0(imgpref, 'autocls_fractions_per_region_individual_', lblset, '.png'), units='in', res=450, width=w, height=h)
# Locations:
layout(matrix(1:9, ncol=1), heights=heights)
sp = 0.1; lsp=1.5
par(xaxs='i')
par(yaxs='i')
par(mar=c(0,lsp, sp, sp))
plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0.5, NR + .5))
red.ct.col = ct.col[c(1,3:9)]
legend('center', legend=names(red.ct.col), col=red.ct.col, 
       pt.cex=1.5, cex=.75, pch=15, bty='n', ncol=8)
for (i in 1:length(mlist)){
    par(mar=c(sp,lsp, sp, sp))
    pad = 0.5
    # Order within?
    sub.mat = mlist[[i]]
    region = regions[i]
    ord = order(sub.mat[,max(grep('Ex',colnames(sub.mat)))])
    sub.mat = sub.mat[ord,]
    plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0.5, nrow(sub.mat) + .5))
    for (j in 1:nrow(sub.mat)){
        rect(xleft=sub.mat[j,-NC], xright=sub.mat[j,-1],
             ybottom=j -pad, ytop=j + pad, col=au.col[colnames(sub.mat)[-1]],
             border=NA, xpd=TRUE)
    }
    rect(xleft=parpos(1,.035), xright=parpos(1,.005),
         ybottom=1 -pad, ytop=nrow(sub.mat) + pad, col=reg.cols[region],
         border=NA, xpd=TRUE)
    text(x=parpos(1,.02), y=(nrow(sub.mat) + 1) / 2, labels=region,
         srt=90, adj=.5, xpd=TRUE, cex=1)
}
par(mar=c(0,lsp, 0, sp))
plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0,1))
segments(x0=0, x1=1, y0=1, y1=1, xpd=TRUE, lwd=.5)
segments(x0=xat, x1=xat, y0=.7, y1=1, xpd=TRUE, lwd=.5)
text(x=xat, y=parpos(2,-.4), labels=xat, xpd=TRUE, cex=.8)
dev.off()



# ------------------------------------------
# Look at the cell fractions by braak stage:
# Pair 1-2, 3-4, 5-6
# ------------------------------------------
covar = 'braaksc'
covstr = sub(" ", "_", sub("\\.", "_", covar))
celldf$braaksc = metadata[celldf$rind, covar]
cmap = colvals[[covar]]

# Cell fractions per region:
ind.countdf = aggregate(lbl ~ region + rind + autocls + braaksc, celldf[celldf$autocls != 'Doublet',], length)
countdf = aggregate(lbl ~ region + autocls + braaksc, ind.countdf, sum)
countdf$mg = paste0(countdf$region, '_', countdf$braaksc)
cwide = spread(countdf[,c('mg','autocls','lbl')], autocls, lbl, fill=0)
cmat = as.matrix(cwide[,-1])
rownames(cmat) = cwide[,1]

# stackfill.barplot = function(mat, horiz=FALSE){
cmat = sweep(cmat, 1, apply(cmat, 1, sum), '/')
mmat = t(rbind(0, apply(cmat, 1, cumsum)))
NR = nrow(mmat)
NC = ncol(mmat)
xat = seq(0,1,.2)
mmat = mmat[rev(rownames(mmat)),]

mlist = lapply(regions, function(x){ mmat[grep(x, rownames(mmat)),]})
heights = c(1, sapply(mlist, nrow) / 6, .75)

w = 5
h = 2.5
png(paste0(imgpref, 'autocls_fractions_per_region_braak_', lblset, '.png'), units='in', res=450, width=w, height=h)
# Locations:
layout(matrix(1:9, ncol=1), heights=heights)
sp = 0.1; lsp=2.25
par(xaxs='i')
par(yaxs='i')
par(mar=c(0,lsp, 0, sp))
plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0,1))
red.ct.col = ct.col[c(1,3:9)]
legend(x=0, y=1.1, legend=names(red.ct.col), col=red.ct.col, 
       pt.cex=1.5, cex=.75, pch=15, bty='n', ncol=8)
legend(x=0, y=0.6, legend=names(cmap), col=cmap, 
       pt.cex=1.5, cex=.75, pch=15, bty='n', ncol=8)
for (i in 1:length(mlist)){
    par(mar=c(sp,lsp, sp, sp))
    pad = 0.5
    # Order within?
    sub.mat = mlist[[i]]
    region = regions[i]
    bks = sub('.*_','', rownames(sub.mat))
    ord = order(bks)
    sub.mat = sub.mat[ord,]
    bks = bks[ord]
    plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0.5, nrow(sub.mat) + .5))
    for (j in 1:nrow(sub.mat)){
        rect(xleft=sub.mat[j,-NC], xright=sub.mat[j,-1],
             ybottom=j -pad, ytop=j + pad, col=au.col[colnames(sub.mat)[-1]],
             border=NA, xpd=TRUE)
        rect(xleft=parpos(1,.025), xright=parpos(1,.005),
             ybottom=j-pad, ytop=j + pad, col=cmap[bks[j]],
             border=NA, xpd=TRUE)
        # text(x=parpos(1,.02), y=j, labels=bks[j], adj=1,  srt=0, xpd=TRUE)
    }
    rect(xleft=parpos(1,.055), xright=parpos(1,.03),
         ybottom=1 -pad, ytop=nrow(sub.mat) + pad, col=reg.cols[region],
         border=NA, xpd=TRUE)
    text(x=parpos(1,.0425), y=(nrow(sub.mat) + 1) / 2, labels=region,
         srt=90, adj=.5, xpd=TRUE, cex=1)
}
par(mar=c(0,lsp, 0, sp))
plot(1,1, type='n', axes=F, xlim=c(0,1), ylim=c(0,1))
segments(x0=0, x1=1, y0=1, y1=1, xpd=TRUE, lwd=.5)
segments(x0=xat, x1=xat, y0=.7, y1=1, xpd=TRUE, lwd=.5)
text(x=xat, y=parpos(2,-.4), labels=xat, xpd=TRUE, cex=.8)
dev.off()

tot.ind.countdf = agg.rename(lbl ~ rind, ind.countdf, sum, 'tot.ind')
ind.countdf  = merge(ind.countdf, tot.ind.countdf)
ind.countdf$ct = substr(ind.countdf$autocls, 1, 2)

gplot = ggplot(ind.countdf[ind.countdf$autocls %in% c('Astro','Oligo','OPC','Endo','Microglia'),], 
               aes(factor(braaksc), lbl/tot.ind, color=factor(braaksc))) + 
    facet_grid(autocls ~ region, scale='free_y') + 
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(width=.25, cex=.05) + 
    scale_color_manual(values=cmap) + 
    scale_y_continuous(labels=scales::percent) + 
    labs(x='Braak Stage', y='Percent of Cells') + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'autocls_fractions_per_region_braak_boxplots_', lblset, '.png'), gplot, dpi=450, units='in', width=6, height=8)


ind.ct.countdf = aggregate(lbl ~ ct + region + rind + tot.ind + braaksc, ind.countdf, sum)

gplot = ggplot(ind.ct.countdf[ind.ct.countdf$ct %in% c('Ex','In'),], 
               aes(factor(braaksc), lbl/tot.ind, color=factor(braaksc))) + 
    facet_grid(ct ~ region, scale='free_y') + 
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(width=.25, cex=.05) + 
    scale_color_manual(values=cmap) + 
    scale_y_continuous(labels=scales::percent) + 
    labs(x='Braak Stage', y='Percent of Cells') + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'autocls_fractions_per_region_braak_boxplots_neu_', lblset, '.png'), gplot, dpi=450, units='in', width=6, height=3)




# --------------------------------------------------------------
# From the identified cell types, which genes mark each cluster:
# --------------------------------------------------------------
# 1) FROM PEC genes
# 2) FROM ALL genes

markmat = as.matrix(gwide[, unique(p2df$Gene)])
rownames(markmat) = as.character(1:nrow(markmat))
markmat = scale(markmat)

ll = reord2d(markmat, 25, 50)
markmat = ll[[1]]
rows = ll[[2]]
cols = ll[[3]]

h = 2.5 * ncol(markmat) / 100
png(paste0(imgpref, 'avg_all_PEC_markers_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
thresh = 3
plt.mat = markmat
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.01), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=.1, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(h=rows, lwd=.25)
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()


# Reduced version:
ind = which(apply(markmat, 2, max) > 4)
plt.mat = markmat[,ind]
h = 2.5 * ncol(plt.mat) / 50
ll = reord2d(plt.mat, 25, 50)
plt.mat = ll[[1]]
rows = ll[[2]]
cols = ll[[3]]

png(paste0(imgpref, 'avg_all_PEC_markers_reduced_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
thresh = 3
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.01), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=.3, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(h=rows, lwd=.25)
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()


# -------------------------------
# Keep only these and replot PEC:
# -------------------------------
p3df = p2df[p2df$Gene %in% colnames(markmat)[ind],]
markmat = as.matrix(gwide[, p3df$Gene])
sort(apply(markmat, 2, max), decreasing=T)
markmat = scale(markmat)
tform = make.tform(p3df$Cluster, norm=T,u=types)
avg.mat = markmat %*% tform
rownames(avg.mat) = as.character(1:nrow(avg.mat))
avg.mat = scale(avg.mat)

ll = reord2d(avg.mat, 25, 5)
avg.mat = ll[[1]]
rows = ll[[2]]
cols = ll[[3]]

png(paste0(imgpref, 'avg_PEC_reduced_markers_', lblset, '.png'), units='in', res=450, width=w, height=2.5)
par(xaxs='i')
par(yaxs='i')
par(mar=c(1,3,0,0))
cex=.4
thresh = 3
plt.mat = avg.mat
plt.mat[plt.mat > thresh] = thresh
image(plt.mat, col=rev(colrb), axes=F, zlim=c(-thresh, thresh))
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.01), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=90, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
abline(h=rows, lwd=.25)
abline(v=cols, lwd=.25)
box(lwd=.5)
dev.off()


# ---------------------------
# Look at one cell type only:
# ---------------------------
# celltype.full = 'Astrocyte'
# celltype = 'Astro'
# celltype.full = 'OPC'
# celltype = 'OPC'
# celltype.full = 'Oligodendrocyte'
# celltype = 'Oligo'
celltype.full = 'Microglia'
celltype = 'Microglia'
sub.lbls = as.character(unique(celldf$lbl[celldf$autocls == celltype]))
not.lbls = hlvls[!(hlvls %in% sub.lbls)]

# 1. Calculate markers for microglia to other cell types:
subexpr = apply(gwide[sub.lbls,],2, mean)
notexpr = apply(gwide[not.lbls,],2, mean)
sdexpr = apply(gwide,2, sd)

# Compare the avg. diff only:
dexpr = subexpr - notexpr
ct.mark = mkdf$symbol[mkdf$cell == cell.map[celltype.full]]
# Markers as significantly differential by SD + expressed at high level
score = (dexpr / sdexpr) * (subexpr > 1) 
plt.genes = names(head(sort(score, decreasing=T), 20))
plt.mat = cbind(notexpr, subexpr)[plt.genes,]
lvlexpr = t(gwide[rev(sub.lbls),])
plt.mat = cbind(lvlexpr, notexpr, subexpr)[plt.genes,]
colnames(plt.mat) = c(colnames(lvlexpr), 'Other cells', paste('All', celltype))

w = nrow(plt.mat) / 8 + .5
h = ncol(plt.mat) / 8 + .5
png(paste0(imgpref, 'overall_markers_', celltype, '_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(2.25,2,.1,.1))
cex=.4
thresh = 3
image(plt.mat, col=colb, axes=F)
axis(1, at=seq(0,1, length.out=nrow(plt.mat)), labels=FALSE, lwd=.25, tck=-0.01)
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.1 * 5 / ncol(plt.mat)), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=45, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
mtext(paste('Learned', celltype, 'Markers'),side=1, cex=.45, line=1.25) 
box(lwd=.5)
dev.off()


# Compare the main cluster to others - Microglia
sublcdf = lblcountdf[lblcountdf$lbl %in% sub.lbls,]
topcls = as.character(sublcdf$lbl[which.max(sublcdf$rind)])

dexpr = lvlexpr[,topcls] - apply(lvlexpr[,sub.lbls[sub.lbls != topcls]],1, mean)
# Markers as significantly differential by SD + expressed at high level
score = abs(dexpr / sdexpr) * (subexpr > 1) 
plt.genes = names(head(sort(score, decreasing=T), 20))
lvlexpr = t(gwide[rev(sub.lbls),])
plt.mat = cbind(lvlexpr, notexpr, subexpr)[plt.genes,]
colnames(plt.mat) = c(colnames(lvlexpr), 'Other cells', paste('All', celltype))

w = nrow(plt.mat) / 8 + .5
h = ncol(plt.mat) / 8 + .5
png(paste0(imgpref, 'ctspecific_markers_', celltype, '_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(2.25,2,.1,.1))
cex=.4
thresh = 3
image(plt.mat, col=colb, axes=F)
axis(1, at=seq(0,1, length.out=nrow(plt.mat)), labels=FALSE, lwd=.25, tck=-0.01)
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.1 * 5 / ncol(plt.mat)), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=45, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
mtext(paste0('Within-', celltype, ' Markers'),side=1, cex=.45, line=1.25) 
box(lwd=.5)
dev.off()

# Comparison of cluster 10 and 80 only for oligodendrocytes:
sublcdf = lblregcountdf[lblregcountdf$lbl %in% sub.lbls & lblregcountdf$region == 'MB',]
sublcdf = sublcdf[sublcdf$lbl != topcls,]
mbcls = as.character(sublcdf$lbl[which.max(sublcdf$rind)])

mainexpr = lvlexpr[,topcls] 
mbexpr = lvlexpr[,mbcls]
dexpr = mainexpr - mbexpr
# Markers as significantly differential by SD + expressed at high level
score = abs(dexpr / sdexpr) * ((mbexpr + mainexpr) / 2> 1) 
plt.genes = names(head(sort(score, decreasing=T), 20))
lvlexpr = t(gwide[rev(sub.lbls),])
plt.mat = cbind(mbexpr, mainexpr)[plt.genes,]
colnames(plt.mat) = c(paste0('MB ', cell.map[celltype.full], ' (',  mbcls, ')'), 
                      paste0('Main ', cell.map[celltype.full], ' (', topcls, ')'))

w = nrow(plt.mat) / 8 + .5
h = ncol(plt.mat) / 8 + .5
png(paste0(imgpref, 'ctspecific_markers_mbcomp_', celltype, '_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(2.25,2.5,.1,.1))
cex=.4
thresh = 3
image(plt.mat, col=colb, axes=F)
axis(1, at=seq(0,1, length.out=nrow(plt.mat)), labels=FALSE, lwd=.25, tck=-0.01)
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.1 * 5 / ncol(plt.mat)), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=45, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
mtext(paste0('Main vs. MB ', celltype, ' Markers'),side=1, cex=.45, line=1.25) 
box(lwd=.5)
dev.off()



# ------------------------
# Plot the cluster groups:
# ------------------------
sub.celldf = celldf[celldf$lbl %in% sub.lbls,]
ind = which(celldf$autocls != 'Doublet')
ind = sample(ind, length(ind), replace=FALSE)
lbl.loc = aggregate(cbind(U1, U2) ~ lbl, celldf, mean)
cex = 0.015
xlim = range(sub.celldf$U1)
ylim = range(sub.celldf$U2)
if(diff(ylim) > diff(xlim)){
    radius = diff(ylim) / 2
    xlim = c(mean(xlim) - radius, mean(xlim) + radius)
} else {
    radius = diff(xlim) / 2
    ylim = c(mean(ylim) - radius, mean(ylim) + radius)
}
auto.loc = aggregate(cbind(U1, U2) ~ lbl + autocls, celldf, mean)
auto.loc = auto.loc[auto.loc$lbl %in% sub.lbls,]

scale = 3
w = 3 * scale
h = 1 * scale * 1.1
png(paste0(imgpref, 'subct_panels_', celltype, '_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
layout(matrix(1:3, nrow=1))
sp = 0.2
par(mar=c(sp,sp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=lbl.cols[celldf$lbl[ind]], 
     pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.1225 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(auto.loc, text(U1, U2, lbl, cex=1))
mtext(paste(ctype, 'Clusters'), side=3, cex=1, col='grey25', font=2, line=0.25)
par(mar=c(sp,sp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=reg.cols[celldf$region[ind]], 
     pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.1225 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
mtext(paste('Brain Regions'), side=3, cex=1, col='grey25', font=2, line=0.25)
legend('right', legend=names(reg.cols), col=reg.cols, 
       pt.cex=1.5, cex=.75, pch=19, bty='n', ncol=2)
par(mar=c(sp,sp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=celldf$col[ind], 
     pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.1225 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
# with(auto.loc, text(U1, U2, autocls, cex=1))
mtext(paste('Celltype Labels'), side=3, cex=1, col='grey25', font=2, line=0.25)
dev.off()



# Comparison of cluster 20 and 29 only for oligodendrocytes:
sublcdf = lblregcountdf[lblregcountdf$lbl %in% sub.lbls & lblregcountdf$region == 'MB',]
sublcdf = sublcdf[sublcdf$lbl != topcls,]
mbcls = as.character(sublcdf$lbl[which.max(sublcdf$rind)])

# sdexpr = apply(gwide,2, sd)
if (celltype == 'Astro') {
    ct1 = '73'; ct2 = '78' 
} else if (celltype == 'Oligo'){
    ct1 = '20'; ct2 = '29'
} else if (celltype == 'OPC'){
    ct1 = '85'; ct2 = '100'
} else if (celltype == 'Microglia'){
    ct1 = '75'; ct2 = '83'
}
mainexpr = lvlexpr[,ct1] 
mbexpr = lvlexpr[,ct2]
dexpr = mainexpr - mbexpr
# Markers as significantly differential by SD + expressed at high level
# score = abs(dexpr / sdexpr) * (subexpr > 1) 
score = abs(dexpr / sdexpr) * ((mbexpr + mainexpr) / 2 > 1)
plt.genes = names(head(sort(score, decreasing=T), 20))
plt.mat = cbind(mbexpr, mainexpr)[plt.genes,]
colnames(plt.mat) = c(paste0('MB-2 ', cell.map[celltype.full], ' (',  ct2, ')'), 
                      paste0('MB-1 ', cell.map[celltype.full], ' (', ct1, ')'))

w = nrow(plt.mat) / 8 + .5
h = ncol(plt.mat) / 8 + .5
png(paste0(imgpref, 'ctspecific_markers_mbcomp_internal_', celltype, '_', lblset, '.png'), units='in', res=450, width=w, height=h)
par(xaxs='i')
par(yaxs='i')
par(mar=c(2.25,2.5,.1,.1))
cex=.4
thresh = 3
image(plt.mat, col=colb, axes=F)
axis(1, at=seq(0,1, length.out=nrow(plt.mat)), labels=FALSE, lwd=.25, tck=-0.01)
text(x=seq(0,1, length.out=nrow(plt.mat)), y=parpos(2,.1 * 5 / ncol(plt.mat)), 
     rownames(plt.mat), xpd=TRUE, cex=cex, srt=45, adj=1)
text(y=seq(0,1, length.out=ncol(plt.mat)), x=parpos(1,.005), 
     colnames(plt.mat), xpd=TRUE, cex=cex, adj=1) # , col=gcolmap[colnames(mat),'col'])
mtext(paste0('MB-1 vs. MB-2 ', celltype, ' Markers'),side=1, cex=.45, line=1.25) 
box(lwd=.5)
dev.off()

