#!/usr/bin/R
# ----------------------------------------------------------
# Plot simple predictions and top genes vs. pathology
# Last updated mid-2020 (exploratory analysis)
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(rhdf5)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/markers/')
imgpref = paste0(plotdir, 'pred_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# ----------------------------------
# Read the umap coords and barcodes:
# ----------------------------------
datadir = 'multiRegion/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy_norm'
bcfile = paste0(datadir, prefix, '.barcodes.tsv.gz')
barcodes = scan(gzfile(bcfile), 'c')
umapfile = paste0(datadir, prefix, '.umap.tsv.gz')

# Groups to plot:
# lblset = 'groups_hdb'
lblset = 'leiden_r5_n50'
lblfile = paste0(datadir, prefix, '.', lblset, '.tsv.gz')
lblavgfile = paste0(datadir, prefix, '.', lblset, '.avg.tsv.gz')
cellfile = paste0(datadir, prefix, '.', lblset, '.lbls.tsv')

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
lbldf = read.delim(cellfile, header=F, stringsAsFactors=F)
names(lbldf) = c('barcode','celltype')
celldf = merge(celldf, lbldf)
write.table(celldf, gzfile(paste0(datadir, prefix, '.', lblset, '.meta.tsv.gz')),
            quote=F, row.names=F, sep="\t", col.names=T)


sum(celldf$lbl == 1) # Unassigned.
hlvls = as.character(sort(unique(celldf$lbl)))
lbl.cols = rep(snap.cols,3)[1:length(hlvls)]
names(lbl.cols) = as.character(hlvls)

if (length(grep('hdb', lblset)) > 0){
    ctype = 'HDBSCAN'
} else {
    ctype = 'Leiden'
}

# -------------------------------
# Load in and plot for each cell:
# -------------------------------
clist = c('Astro','Oligo','OPC', 'Microglia', 'Endo', 'Per')
celltype = 'Microglia'
celltype = 'Oligo'
path = 'nft' # TODO: Generalize to other pathology
# path = 'plaq_n'
go.list = list()
gsea.list = list()
for (celltype in clist){
    print(celltype)
    # Read in pathology and coefficients:
    if (path != 'nft'){ celltype = paste0(path, '_', celltype) }
    coeffdf = read.delim(paste0(datadir, prefix, '.', lblset, '.coeff_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
    corrdf = read.delim(paste0(datadir, prefix, '.', lblset, '.corr_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
    pathdf = read.delim(paste0(datadir, prefix, '.', lblset, '.pathpred_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
    h5file = paste0(datadir, prefix, '.', lblset, '.topgenes_', celltype, '.hdf5')
    pathdf$region = sub('_.*','',pathdf$barcode)
    pathdf = merge(pathdf, celldf, all.x=TRUE)

    lcdf = coeffdf[coeffdf$fit == 'Lasso',]
    lcdf = lcdf[abs(lcdf$coeff) > 0.1,]

    only.goe4=TRUE
    if (!only.goe4){
        gplot = ggplot(pathdf, aes(path, lasso_pred)) + 
            facet_wrap(~region, nrow=1) + 
            geom_point(alpha=.5) + 
            geom_smooth(method='lm') + 
            labs(x='Individual + Region pathology (NFTs)', y='Predicted Pathology') + 
            theme_pubr()
        ggsave(paste0(imgpref, celltype, '_path_vs_lassopath_', lblset, '.png'), gplot, dpi=450, units='in', width=8, height=3)

        # --------------------------
        # Plot the top correlations:
        # --------------------------
        NTOP = 25
        mat = as.matrix(corrdf$corrval[rev(1:NTOP)])
        rownames(mat) = corrdf$symbol[rev(1:NTOP)]

        png(paste0(imgpref, celltype, '_top', NTOP, '_corr_', lblset, '.png'), units='in', res=450, width=1.5, height=6)
        sp = 0.1
        par(mar=c(sp, 5,sp,sp))
        image(t(mat), col=colrb, zlim=c(-max(mat), max(mat)), axes=F)
        text(x=parpos(1, .05),
             y=seq(0,1,length.out=nrow(mat)),
             labels=rownames(mat), xpd=TRUE, adj=1)
        text(x=parpos(1, -.5),
             y=seq(0,1,length.out=nrow(mat)),
             labels=round(mat,3), xpd=TRUE, cex=.8)
        box(lwd=.5)
        dev.off()


        # ----------------------------
        # Plot the lasso coefficients: 
        # ----------------------------
        NR = nrow(lcdf)
        mat = as.matrix(lcdf$coeff[rev(1:NR)])
        rownames(mat) = lcdf$symbol[rev(1:NR)]
        thresh = 1
        zmat = mat
        mat[mat > thresh] = thresh
        mat[mat < -thresh] = -thresh

        png(paste0(imgpref, celltype, '_top', NTOP, '_lasso_coeff_', lblset, '.png'), units='in', res=450, width=1.5, height=6 * NR / 25)
        sp = 0.1
        par(mar=c(sp, 5,sp,sp))
        image(t(mat), col=colrb, zlim=c(-max(abs(mat)), max(abs(mat))), axes=F)
        text(x=parpos(1, .05),
             y=seq(0,1,length.out=nrow(mat)),
             labels=rownames(mat), xpd=TRUE, adj=1)
        text(x=parpos(1, -.5),
             y=seq(0,1,length.out=nrow(mat)),
             labels=round(zmat,3), xpd=TRUE, cex=.8)
        box(lwd=.5)
        dev.off()

        # -------------------------------------------
        # Plot the pathology scores for this one set:
        # -------------------------------------------
        xlim = range(pathdf$U1)
        ylim = range(pathdf$U2)
        if(diff(ylim) > diff(xlim)){
            radius = diff(ylim) / 2
            xlim = c(mean(xlim) - radius, mean(xlim) + radius)
        } else {
            radius = diff(xlim) / 2
            ylim = c(mean(ylim) - radius, mean(ylim) + radius)
        }

        path.zlim = range(log10(pathdf$path + 1), na.rm=T)
        palette = sapply(colr, tsp.col)
        col_fun = function(x, pal=palette){
                bin <- cut(x, seq(path.zlim[1], path.zlim[2], length.out=length(palette)), include.lowest=T) 
                palette[bin] 
        }
        regord = c('EC','HC','MT','PFC','AG')

        h = 2.5
        w = 11
        png(paste0(imgpref, celltype, '_pathology_', path, '_panels.png'), units='in', res=450, width=w, height=h)
        par(xaxs='i')
        par(yaxs='i')
        layout(matrix(1:5, nrow=1, byrow=TRUE), heights=rep(1.1, 1), widths=rep(1, 5), TRUE)
        cex = 0.01
        for (i in 1:length(regord)){
            region = regord[i]
            ind = which(pathdf$region == region)
            ind = sample(ind,length(ind), replace=FALSE) # Shuffle for plotting:
            sp = 0.1
            bsp = 1.5
            par(mar=c(sp,sp,2,sp))
            plot(pathdf$U1[ind], pathdf$U2[ind], col=col_fun(log10(pathdf$path[ind] + 1)),
                 pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
            rect(xleft=par()$usr[1], xright=par()$usr[2],
                 ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
                 ytop=par()$usr[4] + 0.225 * diff(par()$usr[3:4]), 
                 col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
            mtext(reg.long[region], side=3, cex=1.25, col='grey15', font=2, line=0.3)
        }
        dev.off()


        # -------------------------------
        # Plot the predicted pathologies:
        # -------------------------------
        pathdf$lasso_pred = pathdf$lasso_pred - min(pathdf$lasso_pred)
        path.zlim = range(pathdf$lasso_pred, na.rm=T)
        palette = sapply(colr, tsp.col)
        col_fun = function(x, pal=palette){
                bin <- cut(x, seq(path.zlim[1], path.zlim[2], length.out=length(palette)), include.lowest=T) 
                palette[bin] 
        }
        regord = c('EC','HC','MT','PFC','AG')

        h = 2.5
        w = 11
        png(paste0(imgpref, celltype, '_lasso_pathology_', path, '_panels.png'), units='in', res=450, width=w, height=h)
        par(xaxs='i')
        par(yaxs='i')
        layout(matrix(1:5, nrow=1, byrow=TRUE), heights=rep(1.1, 1), widths=rep(1, 5), TRUE)
        cex = 0.01
        for (i in 1:length(regord)){
            region = regord[i]
            ind = which(pathdf$region == region)
            ind = sample(ind,length(ind), replace=FALSE) # Shuffle for plotting:
            sp = 0.1
            bsp = 1.5
            par(mar=c(sp,sp,2,sp))
            plot(pathdf$U1[ind], pathdf$U2[ind], col=col_fun(pathdf$lasso_pred[ind]),
                 pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
            rect(xleft=par()$usr[1], xright=par()$usr[2],
                 ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
                 ytop=par()$usr[4] + 0.225 * diff(par()$usr[3:4]), 
                 col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
            mtext(reg.long[region], side=3, cex=1.25, col='grey15', font=2, line=0.3)
        }
        dev.off()


        # -----------------------------------
        # Plot the values for specific genes:
        # -----------------------------------
        kept.genes = head(lcdf$symbol,5)
        # Slice the matrix for these genes:
        h5f = H5Fopen(h5file)
        genes = h5f$genes
        bcs = h5f$barcodes
        kept.ind = which(genes %in% kept.genes)
        # Open handle, extract genes we care about and close:
        h5d = h5f&"matrix"
        mat = t(h5d[kept.ind,])
        H5Dclose(h5d)
        H5Fclose(h5f)
        colnames(mat) = genes[kept.ind]
        rownames(mat) = bcs

        # Order as pathdf:
        mat = mat[pathdf$barcode,]

        # ---------------------
        # Plot one gene (APOE):
        # ---------------------
        if (length(grep('Microglia', celltype)) > 0){
            gene = 'APOE'
            x = mat[,gene]
            # x[x > 3] = 3
            gene.zlim = range(x, na.rm=T)
            palette = sapply(colr, tsp.col)
            # palette = sapply(viridis(100), tsp.col)
            col_fun = function(x, pal=palette){
                    bin <- cut(x, seq(gene.zlim[1], gene.zlim[2], length.out=length(palette)), include.lowest=T) 
                    palette[bin] 
            }

            h = 2.5
            w = 11
            png(paste0(imgpref, celltype, '_singlegene_', gene, '_panels.png'), units='in', res=450, width=w, height=h)
            par(xaxs='i')
            par(yaxs='i')
            layout(matrix(1:5, nrow=1, byrow=TRUE), heights=rep(1.1, 1), widths=rep(1, 5), TRUE)
            cex = 0.01
            for (i in 1:length(regord)){
                region = regord[i]
                ind = which(pathdf$region == region)
                ind = sample(ind,length(ind), replace=FALSE) # Shuffle for plotting:
                # ind = order(x, decreasing=F)
                sp = 0.1
                bsp = 1.5
                par(mar=c(sp,sp,2,sp))
                plot(pathdf$U1[ind], pathdf$U2[ind], col=col_fun(x[ind]),
                    pch=19, cex=cex, axes=F, xlim=xlim, ylim=ylim)
                rect(xleft=par()$usr[1], xright=par()$usr[2],
                    ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
                    ytop=par()$usr[4] + 0.225 * diff(par()$usr[3:4]), 
                    col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
                mtext(reg.long[region], side=3, cex=1.25, col='grey15', font=2, line=0.3)
            }
            dev.off()
        }


        # --------------------
        # Plot multiple genes:
        # --------------------
        coeff.genes = head(lcdf$symbol,5)
        NG = length(coeff.genes)
        # For mic:
        if (celltype == 'Microglia'){
            cut.ylim = c(mean(ylim), ylim[2])
        } else { cut.ylim = ylim }

        hscale = (.5 + .5 * (celltype != 'Microglia'))
        h = 2.4 * (1 + 1 / 1.1 * (NG - 1)) * hscale
        w = 11
        png(paste0(imgpref, celltype, '_coeffgenes_top',NG, '_panels.png'), units='in', res=450, width=w, height=h)
        par(xaxs='i')
        par(yaxs='i')
        layout(matrix(1:(6 * NG), nrow=NG, byrow=TRUE), heights=c(1.1 * hscale, rep(1 * hscale, NG-1)), widths=c(.3,rep(1, 5)), TRUE)
        cex = 0.01
        for (j in 1:NG){
            gene = coeff.genes[j]
            x = mat[,gene]
            gene.zlim = range(x, na.rm=T)
            palette = sapply(colr, tsp.col)
            col_fun = function(x, pal=palette){
                bin <- cut(x, seq(gene.zlim[1], gene.zlim[2], length.out=length(palette)), include.lowest=T) 
                palette[bin] 
            }
            # Gene only header:
            par(mar=c(sp,sp,sp + 2 * (j==1),sp))
            plot(1,1, type='n', axes=F)
            rect(xleft=parpos(1, -.6), xright=par()$usr[2],
                 ybottom=par()$usr[4], ytop=par()$usr[3], 
                 col='grey85', border=NA, lwd=.5, xpd=TRUE)
            text(x=parpos(1, -.8), y=parpos(2,-.5), gene, cex=1.5, col='grey15', font=2, srt=90, adj=.5)
            for (i in 1:length(regord)){
                region = regord[i]
                ind = which(pathdf$region == region)
                ind = sample(ind,length(ind), replace=FALSE) # Shuffle for plotting:
                # ind = order(x, decreasing=F)
                sp = 0.1
                bsp = 1.5
                par(mar=c(sp,sp,sp + 2 * (j==1),sp))
                plot(pathdf$U1[ind], pathdf$U2[ind], col=col_fun(x[ind]),
                     pch=19, cex=cex, axes=F, xlim=xlim, ylim=cut.ylim)
                if (j==1){
                    rect(xleft=par()$usr[1], xright=par()$usr[2],
                         ybottom=par()$usr[4] + 0.01 * diff(par()$usr[3:4]),
                         ytop=par()$usr[4] + 0.325 * diff(par()$usr[3:4]), 
                         col=reg.cols[region], border=NA, lwd=.5, xpd=TRUE)
                    mtext(reg.long[region], side=3, cex=1.25, col='grey15', font=2, line=0.3)
                }
            }
        }
        dev.off()
    }

    # ------------------------
    # Simple view of pathways:
    # ------------------------
    # 1. Pathways in NFT vs in PLAQ for glial cells
    library(gprofiler2)
    glist = corrdf$symbol[corrdf$corrval > .1]
    # cat(paste(glist, concatenate=' '))
    gtest= gost(glist)
    gtab = gtest$result
    gtab = gtab[order(gtab$p_value),]
    # Save into list:
    go.list[[paste0(celltype, '_', path)]] = gtab

    # Ordered:
    glist = corrdf$symbol
    gtest= gost(glist[1:1000], ordered_query=TRUE)
    gsea.tab = gtest$result
    gsea.tab = gsea.tab[order(gsea.tab$p_value),]
    # gsea.tab[gsea.tab$source == 'REAC',]
    gsea.list[[paste0(celltype, '_', path)]] = gsea.tab


    # 2. Differential use of genes in APOEe4 vs not (in particular APOE?)
    # kept.genes = head(lcdf$symbol,5)
    kept.genes = head(lcdf$symbol,100)
    # Slice the matrix for these genes:
    h5f = H5Fopen(h5file)
    genes = h5f$genes
    bcs = h5f$barcodes
    kept.ind = which(genes %in% kept.genes)
    # Open handle, extract genes we care about and close:
    h5d = h5f&"matrix"
    mat = t(h5d[kept.ind,])
    H5Dclose(h5d)
    H5Fclose(h5f)
    colnames(mat) = genes[kept.ind]
    rownames(mat) = bcs
    mat = mat[pathdf$barcode,]

    # By e4 (only - could also do by region)
    pathdf$has.e4 = metadata[pathdf$rind,'Apoe_e4'] == 'yes'
    i3 = which(!pathdf$has.e4)
    i4 = which(pathdf$has.e4)
    p3 = pathdf$path[i3]
    p4 = pathdf$path[i4]
    gdf = c()
    for (gene in kept.genes){
        gdf = rbind(gdf, data.frame(symbol=gene, c3=cor(p3, mat[i3,gene]), c4=cor(p4, mat[i4,gene])))
    }
    gdf$diff = gdf$c4 - gdf$c3

    # Overall: 
    gmat = as.matrix(gdf[,-1])
    rownames(gmat) = gdf[,1]
    gmat = gmat[rev(1:nrow(gmat)),]
    w = 0.5 + ncol(gmat) / 1.5
    h = 0.5 + nrow(gmat) / 3

    png(paste0(imgpref, celltype, '_', path, '_e4_compare_corr_', lblset, '.png'), units='in', res=450, width=w, height=h)
    sp = 0.1
    par(mar=c(sp, 5,3,sp))
    image(t(gmat), col=rev(colrb), zlim=c(-max(abs(gmat)), max(abs(gmat))), axes=F)
    text(x=parpos(1, .03),
         y=seq(0,1,length.out=nrow(gmat)),
         labels=rownames(gmat), xpd=TRUE, adj=1, cex=.9)
    text(y=parpos(2, -1.012), x=seq(0,1,length.out=3),
             labels=c('No e4', 'e4', 'Diff.'), xpd=TRUE, adj=.5)
    ind = which(gmat !=0,arr.ind=T)
    vals = round(gmat[ind], 2)
    xat = seq(0,1,length.out=3)
    yat = seq(0,1,length.out=nrow(gmat))
    text(x=xat[ind[,2]], y=yat[ind[,1]],
         labels=vals, xpd=TRUE, cex=.8)
    mtext(celltype, side=3, line=2)
    box(lwd=.5)
    dev.off()
}

