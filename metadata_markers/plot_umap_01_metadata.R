#!/usr/bin/R
# --------------------------------------------
# Plot umap with different metadata variables:
# TODO: Which are the HC and PFC individuals
# --------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)

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

# Remove doublet calls:
no.db = TRUE
if (no.db){
    lblset = 'leiden_r5_n50'
    db.bc = scan(paste0(datadir, prefix, '.', lblset, '.dblt.tsv'),'c')
    celldf = celldf[!(celldf$barcode %in% db.bc),]
    imgpref = paste0(imgpref, 'dbrm_')
}

# Number of cells + limits:
NCELL = nrow(celldf)
full.ind = sample(1:NCELL,NCELL, replace=FALSE)
xlim = range(celldf$U1)
ylim = range(celldf$U2)

# --------------------------------------
# Plot first the categorical covariates:
# --------------------------------------
# TODO: ALSO PLOT LEVELS like NFT per region later.
covariates = c('niareagansc','braaksc','amyloid','tangles','cogdx', 'sex','apoe_genotype')
covar = 'niareagansc'

for (covar in covariates){
    covstr = sub(" ", "_", sub("\\.", "_", covar))

    vars = metadata[celldf$rind, covar]
    uqv = unique(vars)
    if (length(uqv) < 10){ 
        vartype = 'categorical' 
        if (covar %in% names(colvals)){
            cmap = colvals[[covar]]
            cols = cmap[as.character(vars)]
        }
    } else { 
        vartype = 'numeric' 
    }

    cex = 0.02
    png(paste0(imgpref, 'covariate_', covar, '.png'), units='in', res=450, width=8, height=8)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    bsp = 1.5
    par(mar=c(bsp,bsp,2,sp))
    plot(celldf$U1[full.ind], celldf$U2[full.ind], col=cols[full.ind], 
         pch=19, cex=cex, axes=F)
    legend('topleft', legend=names(cmap), col=cmap, 
           pt.cex=2.5, cex=1.25, pch=19, bty='n', ncol=2)
    rect(xleft=par()$usr[1], xright=par()$usr[2],
         ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
         ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
         col='grey85', border=NA, lwd=.5, xpd=TRUE)
    mtext(covar, side=3, cex=1.5, col='grey25', font=2, line=0.25)
    mtext('UMAP 1', side=1, line=0.25, cex=1.25)
    mtext('UMAP 2', side=2, line=0, cex=1.25)
    dev.off()

    # Panel by panel:
    png(paste0(imgpref, 'covariate_', covar, '_panels.png'), units='in', res=450, width=10, height=5.5)
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
        plot(celldf$U1[ind], celldf$U2[ind], col=cols[ind], 
             pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
        rect(xleft=par()$usr[1], xright=par()$usr[2],
             ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
             ytop=par()$usr[4] + 0.1225 * diff(par()$usr[3:4]), 
             col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
        mtext(reg.long[region], side=3, cex=1.25, col='grey15', font=2, line=0.3)
        if (i == 1){
            mtext('UMAP 1', side=1, line=0.25, cex=.8)
            mtext('UMAP 2', side=2, line=0, cex=.8)
        }
    }
    par(mar=rep(0.1,4))
    plot(1,1,type='n', axes=F, ylab='', xlab='')
    mtext(covar, side=3, line=0)
    legend('top', legend=names(cmap), col=cmap, 
           pt.cex=2, cex=1.25, pch=19, bty='n', ncol=2)
    dev.off()
}


# ------------------------------------------
# Plot pathology scores on the regions umap:
# ------------------------------------------
regmap = c('AG','HC','PFC','MT','EC')
names(regmap) = c('ag','hip','mf','mt','ec')
for (path in c('nft','plaq_d','plaq_n')){
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    celldf$path = slong[celldf$rind,'value']

    path.zlim = range(log10(celldf$path + 1), na.rm=T)
    palette = sapply(viridis(100), tsp.col)
    # palette = sapply(rev(colrb), tsp.col)
    palette = sapply(colr, tsp.col)
    col_fun = function(x, pal=palette){
            bin <- cut(x, seq(path.zlim[1], path.zlim[2], length.out=length(palette)), include.lowest=T) 
            palette[bin] 
    }

    # Panel by panel:
    mvdf = aggregate(value ~ region, slong, function(x){mean(x) + median(x)})
    regord = mvdf[order(mvdf$value, decreasing=T), 'region']
    regord = c('EC','HC','MT','PFC','AG')

    h = 2.5
    w = 11
    png(paste0(imgpref, 'pathology_', path, '_panels.png'), units='in', res=450, width=w, height=h)
    par(xaxs='i')
    par(yaxs='i')
    layout(matrix(1:5, nrow=1, byrow=TRUE), heights=rep(1.1, 1), widths=rep(1, 5), TRUE)
    cex = 0.01
    for (i in 1:length(regord)){
        region = regord[i]
        ind = which(celldf$region == region)
        ind = sample(ind,length(ind), replace=FALSE) # Shuffle for plotting:
        sp = 0.1
        bsp = 1.5
        par(mar=c(sp,sp,2,sp))
        plot(celldf$U1[ind], celldf$U2[ind], col=col_fun(log10(celldf$path[ind] + 1)),
             pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
        rect(xleft=par()$usr[1], xright=par()$usr[2],
             ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
             ytop=par()$usr[4] + 0.225 * diff(par()$usr[3:4]), 
             col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
        # mtext(region, side=3, cex=1.5, col='grey15', font=2, line=0.25)
        mtext(reg.long[region], side=3, cex=1.25, col='grey15', font=2, line=0.3)
    }
    dev.off()

    pathfile = paste0(datadir, prefix, '.', path, '.tsv.gz')
    write.table(celldf[,c('barcode','path')], gzfile(pathfile), col.names=F, row.names=F, sep="\t", quote=F)
}


