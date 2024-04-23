#!/usr/bin/R
# ------------------------------------------------
# Plot UMAPs for each of the celltypes separately:
# Updated 11/04/2021 
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)

# Directories:
plotdir = paste0(imgdir, 'markers/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir)
system(cmd)


# List of chosen runs corresponding to each cell type:
# ----------------------------------------------------
runlist = list()
runlist[['Exc']] = 'log1p_combat_filthvg1000_nomt'
runlist[['Inh']] = 'log1p_combat_filthvg1000_nomt'

rownames(cellmeta) = cellmeta$barcode


# Plot the UMAPs for each celltype
# --------------------------------
for (celltype in names(runlist)){
    print(celltype)
    cellsuff = paste0(celltype, '_', runlist[[celltype]])
    cellpref = paste0(imgpref, 'umap_highres_cols_subset_', cellsuff)
    fname = paste0(sdbdir, 'metadata/multiregion_majorctUMAP_', cellsuff, '.tsv')
    df = read.delim(fname, header=T, sep="\t")
    df = df[, c('barcode','C1','C2')]

    # Merge the metadata in: 
    df$region = cellmeta[df$barcode, 'region']
    df$projid = cellmeta[df$barcode, 'projid']
    df$cthr = cellmeta[df$barcode, 'cell_type_high_resolution']


    # ------------------------------------
    # Plot with own colors for full types:
    # ------------------------------------
    ind = sample(1:nrow(df), nrow(df), replace=F)
    celltype.loc = aggregate(cbind(C1, C2) ~ cthr, df, mean)
    cex = 0.025
    ctdf = agg.rename(barcode ~ cthr, df, length, 'count')
    ctdf$frac = ctdf$count / sum(ctdf$count)
    ctdf$pct = round(ctdf$frac * 100, 2)
    celltype.loc = merge(celltype.loc, ctdf)

    w = 4
   # No text resource:
    png(paste0(cellpref, '_notext.png'), units='in', res=450, width=w, height=w)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    par(mar=rep(sp, 4))
    plot(df$C1[ind], df$C2[ind], col=tsp.tcols[df$cthr[ind]], pch=19, cex=cex, axes=F)
    dev.off()

    # Text-only resource:
    pdf(paste0(cellpref, '_textonly.pdf'), width=w, height=w)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    par(mar=rep(sp, 4))
    plot(df$C1[ind], df$C2[ind], type='n', col=tsp.tcols[df$cthr[ind]], pch=19, cex=cex, axes=F)
    box()
    with(celltype.loc, text(C1, C2, paste0(cthr, ' ', count, ' (', pct, '%)'), cex=.75, xpd=TRUE, font=2))
    dev.off()


    # Plot with region colors:
    # ------------------------
    w = 4
   # No text resource:
    png(paste0(cellpref, '_notext_region.png'), units='in', res=450, width=w, height=w)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    par(mar=rep(sp, 4))
    plot(df$C1[ind], df$C2[ind], col=reg.cols[df$region[ind]], pch=19, cex=cex, axes=F)
    dev.off()

    # Plot 2x3 for regions:
    # ---------------------
    w = 4
    sp = 0.1
   # No text resource:
    png(paste0(cellpref, '_notext_region_2x3.png'), units='in', res=450, width=w * 3, height=w * 2)
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    layout(matrix(1:6, nrow=2, byrow=TRUE))
    for (region in reg.order[-1]){
        subind = ind[df$region[ind] == region]
        plot(df$C1[subind], df$C2[subind], col=reg.cols[df$region[subind]], pch=19, cex=cex, axes=F)
    }
    dev.off()
}

