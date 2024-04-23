#!/usr/bin/R
# ------------------------------------------------------
# Second pass of refining metadata:
# Final clusters from Hans, names for exc, removal of MB
# Updated 10/11/2020 with names / cluster 172 changes
# ------------------------------------------------------
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


final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
if (!file.exists(final.rdafile)){
    # Load old meta (cellmeta):
    load(paste0(datadir, prefix, '.final.cell_labels.Rda'))

    # -------------------------------------
    # Load in second pass Hans assignments:
    # -------------------------------------
    metadir = 'multiRegion/metadata/'
    fnames = list.files(path=metadir, pattern='Metadata_.*.csv')
    newmeta = c()
    for (fnam in fnames){
        print(fnam)
        df = read.delim(paste0(metadir, fnam), sep=",")
        print(dim(df))
        df = df[, c(names(cellmeta), 'cell_type_high_resolution')]
        newmeta = rbind(newmeta, df)
    }

    # Size of newmeta:
    print(dim(newmeta))
    # Size of old metadata:
    print(dim(cellmeta))
    print(dim(cellmeta[cellmeta$region != 'MB',]))
    # Update the neuronal subclusters based on these labels?
    newmeta$neuronal.layer[newmeta$neuronal.layer == 'nan'] = NA

    # -----------------
    # Save annotations:
    # -----------------
    celldf = newmeta
    ind = which(!is.na(celldf$full.exttype))
    celldf$col = major.col[celldf$hcelltype]
    celldf$tsp.col = tsp.major.col[celldf$hcelltype]
    keep.cols = c("lbl", "U1", "U2", "barcode", "rind", "region", "projid", "is.doublet", 
                  "col", "tspcol", "hcluster", "hcelltype", "hsubclass", 
                  "major.celltype", "minor.celltype", "neuronal.layer", "inh.subtype", 
                  "neuronal.exttype", "full.exttype",'cell_type_high_resolution')  
    write.table(celldf[ind,keep.cols], gzfile(paste0(datadir, prefix, '.final_noMB.cell_labels.tsv.gz')), 
                quote=F, row.names=F, sep="\t", col.names=T)
    cellmeta = celldf[ind,keep.cols]
    save(cellmeta, file=final.rdafile)
} else { 
    load(final.rdafile)
}


# ------------------------------------
# Plot with own colors for full types:
# ------------------------------------
ind = 1:nrow(cellmeta)
typelvls = unique(cellmeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
tsp.type.cols = sapply(type.cols, tsp.col)
celltype.loc = aggregate(cbind(U1, U2) ~ cell_type_high_resolution, cellmeta, mean)
cex = 0.025

png(paste0(imgpref, 'umap_highres_cols_', lblset, '_wout_doublets.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(cellmeta$U1[ind], cellmeta$U2[ind], col=tsp.type.cols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
# with(lbl.loc, text(U1, U2, lbl, cex=.5, xpd=TRUE, font=1))
with(celltype.loc, text(U1, U2, cell_type_high_resolution, cex=.75, xpd=TRUE, font=2))
mtext(paste('Cell Types (High-res)'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()

# -------------------------------
# Plot with new UMAP coordinates:
# -------------------------------
udf = read.delim(gzfile(paste0('multiRegion/',prefix, '_filtered.umap.tsv.gz')), header=F, sep="\t")
bcs = scan(gzfile(paste0('multiRegion/',prefix, '_filtered.barcodes.tsv.gz')), 'c', quiet=T)
names(udf) = c('UMAP1', 'UMAP2')
rownames(udf) = bcs

cellmeta$UMAP1 = udf[cellmeta$barcode,'UMAP1']
cellmeta$UMAP2 = udf[cellmeta$barcode,'UMAP2']

ind = 1:nrow(cellmeta)
typelvls = unique(cellmeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
tsp.type.cols = sapply(type.cols, tsp.col, alpha=.25)
celltype.loc = aggregate(cbind(UMAP1, UMAP2) ~ cell_type_high_resolution, cellmeta, mean)
cex = 0.02

png(paste0(imgpref, 'umap_final_highres_cols_', lblset, '_wout_doublets.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.type.cols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
# with(lbl.loc, text(UMAP1, UMAP2, lbl, cex=.5, xpd=TRUE, font=1))
with(celltype.loc, text(UMAP1, UMAP2, cell_type_high_resolution, cex=.75, xpd=TRUE, font=2))
mtext(paste('Cell Types (High-res)'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()

# No text resource:
png(paste0(imgpref, 'umap_final_highres_cols_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.type.cols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
dev.off()

# Text-only resource:
pdf(paste0(imgpref, 'umap_final_highres_cols_', lblset, '_wout_doublets_textonly.pdf'), width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], type='n', col=tsp.type.cols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
box()
with(celltype.loc, text(UMAP1, UMAP2, cell_type_high_resolution, cex=.75, xpd=TRUE, font=2))
dev.off()


# Brain regions resource:
regnmb = regions[regions != 'MB']
tsp.reg.cols = sapply(reg.cols, tsp.col, alpha=.5)
lab.reg.cols = reg.cols[regnmb]
names(lab.reg.cols) = c('Ang. Gyrus', 'Ent. Cortex', 'Hippocampus',
                       'Mid-Tmp Ctx', 'Pre-Ftl Ctx', 'Thalamus')

png(paste0(imgpref, 'umap_final_highres_cols_', lblset, '_wout_doublets_notext_regions.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.reg.cols[cellmeta$region[ind]], 
     pch=19, cex=cex, axes=F)
legend('topleft', legend=names(lab.reg.cols), col=lab.reg.cols, 
           pt.cex=2.5, cex=1, pch=15, bty='n', ncol=2)
dev.off()

w = 2
png(paste0(imgpref, 'umap_final_highres_cols_', lblset, '_wout_doublets_notext_regions_small.png'), units='in', res=450, width=w, height=w)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
# plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.reg.cols[cellmeta$region[ind]], pch='.', cex=.25, axes=F)
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=reg.cols[cellmeta$region[ind]], pch='.', cex=.25, axes=F)
dev.off()


# --------------------------------------
# Improve the colors for the cell types:
# --------------------------------------
# 1. Same color scale for each type
# 2. Diverse colors for each type
# 3. Brighter exc, pastel inh
tcols = type.cols
tnam = names(tcols)
# tcols[grep("^Ast", tnam)] = brewer.pal(5, 'Reds')[2:5] # OK for now
tcols[grep("^Oli", tnam)] = brewer.pal(5, 'Oranges')[c(2,4)]
tcols[c(grep("^Mic", tnam), which(tnam %in% c('CAMs','T cells')))] = c(brewer.pal(6, 'Purples')[3:6], 'goldenrod1','firebrick')
tcols[c(grep("^OPC", tnam))] = colorRampPalette(c('white','saddlebrown'))(5)[3:5]

tcols['CA1 pyramidal cells'] = 'grey80'
tcols['DG granule cells'] = '#6699CC'
tcols['Inh PTPRK FAM19A1'] = 'grey85'
tcols['Inh VIP THSD7B'] = '#44AA99'
a = tcols['Inh PVALB SULF1']
tcols['Inh PVALB SULF1'] = tcols['Inh ENOX2 SPHKAP']
tcols['Inh ENOX2 SPHKAP'] = a

a = tcols['Exc L5-6 RORB LINC02196']
tcols['Exc L5-6 RORB LINC02196'] = tcols['Exc NRGN']
tcols['Exc NRGN'] = a
tcols['Exc SOX11 NCKAP5'] = '#DDCC77'
tcols['Exc L2-3 CBLN2 LINC02306'] = brewer.pal(9, 'Blues')[5]
tcols['Exc L3-5 RORB PLCH1'] = '#DA6D7D' 

tsp.tcols = sapply(tcols, tsp.col)

png(paste0('~/test_legend.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
par(mar=rep(0,4))
plot(1,1,type='n', axes=F)
legend('center', legend=names(tcols), col=tcols, 
           pt.cex=2.5, cex=1, pch=15, bty='n', ncol=2)
dev.off()


png(paste0(imgpref, 'umap_final_highres_tcols_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.tcols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
dev.off()

save(tcols,file='Annotation/multiregion_celltypes_colors.Rda')
coldf = data.frame(color=tcols)
write.table(coldf,file='Annotation/multiregion_celltypes_colors.tsv', quote=F, sep="\t")


# ------------------------------
# Merged cols for microglia, OPC
# ------------------------------
mcols = tcols
mcols[grep("^Mic", tnam)] = major.col['Mic/Immune']
mcols[grep("^OPC", tnam)] = major.col['Opc']
tsp.mcols = sapply(mcols, tsp.col)

png(paste0(imgpref, 'umap_final_highres_merged_cols_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.mcols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
dev.off()

# Highlighting Exc subsets:
set = 'Exc'
cts = unique(cellmeta$cell_type_high_resolution[cellmeta$major.celltype == set])
gcols = tcols
gcols[!(names(gcols) %in% cts)] = 'grey95'
tsp.gcols = sapply(gcols, tsp.col)
png(paste0(imgpref, 'umap_final_highres_highlight_exc_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.gcols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
dev.off()

# Highlighting Inh subsets:
set = 'Inh'
cts = unique(cellmeta$cell_type_high_resolution[cellmeta$major.celltype == set])
gcols = tcols
gcols[!(names(gcols) %in% cts)] = 'grey95'
tsp.gcols = sapply(gcols, tsp.col)
png(paste0(imgpref, 'umap_final_highres_highlight_inh_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=tsp.gcols[cellmeta$cell_type_high_resolution[ind]], 
     pch=19, cex=cex, axes=F)
dev.off()



# --------------------
# Fully merged colors:
# --------------------
png(paste0(imgpref, 'umap_final_highres_major_cols_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=major.col[cellmeta$hcelltype[ind]], pch=19, cex=cex, axes=F)
dev.off()

w = 2
png(paste0(imgpref, 'umap_final_highres_major_cols_', lblset, '_wout_doublets_notext_small.png'), units='in', res=450, width=w, height=w)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(cellmeta$UMAP1[ind], cellmeta$UMAP2[ind], col=major.col[cellmeta$hcelltype[ind]], pch='.', cex=.25, axes=F)
dev.off()


# ------------------------
# Make the fraction plots:
# ------------------------
margdf = agg.rename(barcode ~ region, cellmeta, length, 'tot')
ctdf = aggregate(barcode ~ region + hcelltype, cellmeta, length)
ctdf = merge(ctdf, margdf)
totdf = aggregate(barcode ~ hcelltype, cellmeta, length)
totdf$tot = sum(totdf$barcode)
totdf$region = 'All'
ctdf = rbind(ctdf, totdf)
ctdf$fracloc = ctdf$barcode / ctdf$tot
ctdf$frac = paste0(round(ctdf$barcode / ctdf$tot * 100,1),'%')

ctdf$region = factor(ctdf$region, levels= c('All','EC','HC','AG','MT','PFC','TH'))
gplot = ggplot(ctdf, aes(region, barcode, fill=hcelltype)) + 
    geom_bar(position='fill', stat='identity') + 
    # geom_text(data=ctdf, aes(region, fracloc, label=frac)) + 
    scale_fill_manual(values=major.col) + 
    theme_pubr() + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    labs(x='Region', y='Percentage')

ggsave(paste0(imgpref, 'fractions_final_hct_major_cols_', lblset, '.pdf'), gplot, dpi=450, units='in', width=5, height=6)

# ----------------------------------------------
# Make the fraction plot for specific celltypes:
# ----------------------------------------------
set = 'Inh'
margdf = agg.rename(barcode ~ region, cellmeta[cellmeta$major.celltype == set,], length, 'tot')
ctdf = aggregate(barcode ~ region + cell_type_high_resolution, cellmeta[cellmeta$major.celltype == set,], length)
ctdf = merge(ctdf, margdf)
totdf = aggregate(barcode ~ cell_type_high_resolution, cellmeta[cellmeta$major.celltype == set,], length)
totdf$tot = sum(totdf$barcode)
totdf$region = 'All'
ctdf = rbind(ctdf, totdf)
ctdf$fracloc = ctdf$barcode / ctdf$tot
ctdf$frac = paste0(round(ctdf$barcode / ctdf$tot * 100,1),'%')

ctdf$region = factor(ctdf$region, levels= c('All','EC','HC','AG','MT','PFC','TH'))
gplot = ggplot(ctdf, aes(region, barcode, fill=cell_type_high_resolution)) + 
    geom_bar(position='fill', stat='identity') + 
    # geom_text(data=ctdf, aes(region, fracloc, label=frac)) + 
    scale_fill_manual(values=tcols) + 
    theme_pubr() + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    labs(x='Region', y='Percentage')

ctdf = ctdf[ctdf$region != 'All',]
ctdf$region = factor(ctdf$region, levels= c('EC','HC','AG','MT','PFC','TH'))
gplot = ggplot(ctdf, aes(cell_type_high_resolution, barcode, fill=region)) + 
    geom_bar(position='fill', stat='identity') + 
    # geom_text(data=ctdf, aes(region, fracloc, label=frac)) + 
    coord_flip() + 
    scale_fill_manual(values=reg.cols) + 
    theme_pubr() + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    labs(x='Region', y='Percentage')
ggsave(paste0(imgpref, 'fractions_final_',sub("/","_", set),'_regions_', lblset, '.pdf'), gplot, dpi=450, units='in', width=5, height=6)


# --------------------------------------------------
# Update the metadata with the new umap coordinates:
# --------------------------------------------------
if (!is.null(cellmeta$UMAP1)){
    cellmeta$U1 = cellmeta$UMAP1
    cellmeta$U2 = cellmeta$UMAP2
    cellmeta$UMAP1 = NULL
    cellmeta$UMAP2 = NULL
}

ind = which(!is.na(cellmeta$full.exttype))
cellmeta = cellmeta[ind,]
write.table(cellmeta, gzfile(paste0(datadir, prefix, '.final_noMB.cell_labels.tsv.gz')), 
            quote=F, row.names=F, sep="\t", col.names=T)
save(cellmeta, file=final.rdafile)

