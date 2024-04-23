#!/usr/bin/R
# -----------------------------------------------------
# Plot umap with the batches and with indivual metadata
# TODO: Which are the HC and PFC individuals
# -----------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = paste0(plotdir, 'umap_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

metadata$rind = rownames(metadata)

# ----------------------------------
# Read the umap coords and barcodes:
# ----------------------------------
datadir = 'multiRegion/'
# prefix = 'all_brain_regions_filt_preprocessed_norm'
# bcfile = paste0(datadir,'all_brain_regions_filt_preprocessed_barcodes.txt')
# barcodes = scan(bcfile, 'c')
prefix = 'all_brain_regions_filt_preprocessed_scanpy_norm'
bcfile = paste0(datadir, prefix, '.barcodes.tsv.gz')
barcodes = scan(gzfile(bcfile), 'c')
umapfile = paste0(datadir, prefix, '.umap.tsv.gz')

# Read data:
celldf = read.delim(gzfile(umapfile), header=F, sep="\t")
names(celldf) = c('U1','U2')
celldf$rind = sub("_.*-",".",barcodes)
celldf$region = sub("_.*","",barcodes)
celldf$projid = metadata[celldf$rind,'projid']
celldf$barcode = barcodes

no.db = TRUE
if (no.db){
    lblset = 'leiden_r5_n50'
    db.bc = scan(paste0(datadir, prefix, '.', lblset, '.dblt.tsv'),'c')
    celldf = celldf[!(celldf$barcode %in% db.bc),]
    imgpref = paste0(imgpref, 'dbrm_')
}


# ------------------------
# Plot the region batches:
# ------------------------
cex = 0.025

png(paste0(imgpref, 'basic.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1, celldf$U2, col='grey85', 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col=tsp.col('grey85', alpha=0.5), border=NA, lwd=.5, xpd=TRUE)
mtext('UMAP', side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()


NCELL = nrow(celldf)
ind = sample(1:NCELL,NCELL, replace=FALSE)
treg.cols = sapply(reg.cols, tsp.col)

png(paste0(imgpref, 'batch.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=treg.cols[celldf$region[ind]], 
     pch=19, cex=cex, axes=F)
# legend('topleft', legend=names(reg.cols), col=reg.cols, 
# pt.cex=2.5, cex=1.25, pch=19, bty='n', ncol=2)
legend('topleft', legend=reg.long[names(reg.cols)][1:4], col=reg.cols[1:4], 
       pt.cex=2.5, cex=1.2, pch=19, bty='n', ncol=1)
legend('topright', legend=reg.long[names(reg.cols)][5:7], col=reg.cols[5:7], 
       pt.cex=2.5, cex=1.2, pch=19, bty='n', ncol=1)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
mtext('UMAP by Brain Region', side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()


# -----------------------------------
# Plot the region batches separately:
# -----------------------------------
xlim = range(celldf$U1)
ylim = range(celldf$U2)
png(paste0(imgpref, 'batch_panels.png'), units='in', res=450, width=10, height=5.5)
par(xaxs='i')
par(yaxs='i')
layout(matrix(1:8, nrow=2, byrow=TRUE))
for (i in 1:length(reg.cols)){
    ind = which(celldf$region == regions[i])
    sp = 0.1
    bsp = 1.5
    par(mar=c(bsp,bsp,2,sp))
    plot(celldf$U1[ind], celldf$U2[ind], 
         col=treg.cols[celldf$region[ind]], 
         pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
    rect(xleft=par()$usr[1], xright=par()$usr[2],
         ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
         ytop=par()$usr[4] + 0.1225 * diff(par()$usr[3:4]), 
         col=reg.cols[regions[i]], border=NA, lwd=.5, xpd=TRUE)
    mtext(reg.long[regions[i]], side=3, cex=1.25, col='grey15', font=2, line=0.3)
    if (i == 1){
        mtext('UMAP 1', side=1, line=0.25, cex=.8)
        mtext('UMAP 2', side=2, line=0, cex=.8)
    }
}
dev.off()


# -------------------------------------
# In each region, plot the individuals:
# -------------------------------------
indvs = sort(unique(celldf$projid))
ind.cols = snap.cols[1:length(indvs)]
names(ind.cols) = as.character(indvs)

png(paste0(imgpref, 'indbatch_panels.png'), units='in', res=450, width=10, height=5.5)
par(xaxs='i')
par(yaxs='i')
layout(matrix(1:8, nrow=2, byrow=TRUE))
cex = 0.01
for (i in 1:length(reg.cols)){
    region = regions[i]
    ind = which(celldf$region == region)
    ind = sample(ind,length(ind), replace=FALSE) # Shuffle for plotting:
    sp = 0.1
    bsp = 1.5
    par(mar=c(bsp,bsp,2,sp))
    plot(celldf$U1[ind], celldf$U2[ind], 
         col=ind.cols[as.character(celldf$projid[ind])], 
         pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
    rect(xleft=par()$usr[1], xright=par()$usr[2],
         ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
         ytop=par()$usr[4] + 0.1225 * diff(par()$usr[3:4]), 
         col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
    # mtext(region, side=3, cex=1.5, col='grey15', font=2, line=0.25)
    mtext(reg.long[regions[i]], side=3, cex=1.25, col='grey15', font=2, line=0.3)
    if (i == 1){
        mtext('UMAP 1', side=1, line=0.25, cex=.8)
        mtext('UMAP 2', side=2, line=0, cex=.8)
    }
}
par(mar=rep(0.1,4))
plot(1,1,type='n', axes=F, ylab='', xlab='')
mtext('Individuals:', side=3, line=0)
legend('top', legend=names(ind.cols), col=ind.cols, 
       pt.cex=1.5, cex=.8, pch=19, bty='n', ncol=4)
dev.off()

