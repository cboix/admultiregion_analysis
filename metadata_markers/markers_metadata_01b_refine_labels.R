#!/usr/bin/R
# ----------------------------------------------------------
# Plot the average gene expression for markers for clusters:
# Updated 09/17/2020 with final pass of doublet idenfitication.
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


# --------------------------
# Load in Hans' assignments:
# --------------------------
edf = read.delim('multiRegion/hans_excitatory_markers.tsv', header=T)
idf = read.delim('multiRegion/hans_inhibitory_markers.tsv', header=T)
ndf = read.delim('multiRegion/hans_nonneuronal_markers.tsv', header=T)
edf$lbl = sub(".*_","", edf$cls_lbl)
idf$lbl = sub(".*_","", idf$cls_lbl)
ndf$lbl = sub(".*_","", ndf$cls_lbl)
edf[edf$cls_lbl == 'Ex_224_Stressed_neurons',] = c('Ex_224', 'Exc_L2-3_Stressed', '224')
edf[edf$cls_lbl == 'Ex_113_L6_CT',] = c('Ex_113', 'Exc_L6_CT', '113')
idf[idf$lbl %in% c('231','288'),'hcluster'] = 'Inh_PVALB'
# Additional filtering:
edf[edf$lbl %in% c('239','14','334','344','298','297','291'),'hcluster'] = 'Doublet'
edf$hcelltype = 'Exc'
idf$hcelltype = 'Inh'
ndf$hcelltype = ndf$hsubclass
idf$hsubclass[idf$hsubclass == 'Doublet'] = NA

# Additional fixes (may want to recluster - Hans?)
# Non 257 (mixture of multiple types of cells)


hcdf = rbind(edf[,c('lbl','hcluster','hcelltype')],
             idf[,c('lbl','hcluster','hcelltype')],
             ndf[,c('lbl','hcluster','hcelltype')])

celldf = merge(celldf, hcdf, all.x=TRUE)
celldf = merge(celldf, idf[,c('hsubclass','lbl')], all.x=TRUE)

# --------------------
# Finalize the labels:
# --------------------
celldf$hcelltype[celldf$hcluster == 'Doublet'] = NA
celldf$hcluster[celldf$hcluster == 'Doublet'] = NA
celldf$hcelltype[celldf$hcelltype == 'single_projid'] = NA
celldf$hcluster[celldf$hcluster == 'single_projid'] = NA
celldf$major.celltype = celldf$hcelltype
celldf$minor.celltype = celldf$hcelltype

# Relabel the vascular / epithelia:
vind = which(celldf$hcelltype == 'Vasc/Epithelia') 
celldf$minor.celltype[vind] = celldf$hcluster[vind]
celldf$minor.celltype[celldf$minor.celltype == 'choroid_plexus_epithelial_cells'] = 'CPEC'
celldf$minor.celltype[celldf$minor.celltype == 'endothelial_cells'] = 'End'
celldf$minor.celltype[celldf$minor.celltype == 'ependymal_cells'] = 'Epd'
celldf$minor.celltype[celldf$minor.celltype == 'fibroblasts'] = 'Fib'
celldf$minor.celltype[celldf$minor.celltype == 'pericytes'] = 'Per'
celldf$minor.celltype[celldf$minor.celltype == 'smooth_muscle_cells'] = 'SMC'

# Relabel the vascular / epithelia:
vind = which(celldf$hcelltype == 'Mic/Immune') 
celldf$minor.celltype[vind] = celldf$hcluster[vind]
celldf$minor.celltype[celldf$minor.celltype == 'CAMs'] = 'CAM'
celldf$minor.celltype[celldf$minor.celltype == 'Mic'] = 'Mic'
celldf$minor.celltype[celldf$minor.celltype == 'T_cells'] = 'T'
print(table(celldf[, 'minor.celltype']))
print(sum(table(celldf[, 'minor.celltype'])))

# -------------------------------------------
# Add the neuronal markers on the next level:
# -------------------------------------------
# Layer markers where available:
celldf$neuronal.layer = ""
vind = which(celldf$minor.celltype %in% c('Exc','Inh')) 
celldf$neuronal.layer[vind] = sub("_.*","", sub(".*_L([0-9])","L\\1", celldf$hcluster[vind]))
lind = grep("^L[0-9]", celldf$neuronal.layer, invert=T)
celldf$neuronal.layer[lind] = ""
print(table(celldf[, 'neuronal.layer']))

# Neuronal subtype for inhibitory neurons:
celldf$inh.subtype = ""
hind = !is.na(celldf$hsubclass)
celldf$inh.subtype[hind] = celldf$hsubclass[hind]
print(table(celldf[, 'inh.subtype']))

# Neuronal subtype for other neurons:
celldf$neuronal.exttype = celldf$inh.subtype
eind = which(celldf$minor.celltype == 'Exc') 
celldf$neuronal.exttype[eind] = celldf$hcluster[eind]
eind = which(celldf$minor.celltype == 'Exc' & celldf$neuronal.layer != "") 
celldf$neuronal.exttype[eind] = celldf$neuronal.layer[eind]
celldf$neuronal.exttype[celldf$hcluster == 'DG_granule_cells'] = 'Granule'
celldf$neuronal.exttype[celldf$hcluster == 'Inh_L1-4_LAMP5_DUSP4_Rosehip'] = 'Rosehip'
celldf$neuronal.exttype[celldf$hcluster == 'CA1_pyramidal_cells'] = 'CA1_Pyr'
celldf$neuronal.exttype[celldf$hcluster == 'CA2_CA3_pyramidal_cells'] = 'CA23_Pyr'
celldf$neuronal.exttype[celldf$hcluster == 'Subcubiculum'] = 'Sbc'
celldf$neuronal.exttype[celldf$hcluster == 'Ex-NRGN'] = 'NRGN'
table(celldf$neuronal.exttype)

# Aggregate the minor cell type with the neuronal exttype:
vind = which(celldf$minor.celltype %in% c('Exc','Inh')) 
celldf$full.exttype = celldf$minor.celltype 
celldf$full.exttype[vind] = celldf$neuronal.exttype[vind]

table(celldf$full.exttype)
sum(table(celldf$full.exttype))

# Finally, update doublet calculations:
celldf$is.doublet = is.na(celldf$full.exttype)

# -----------------
# Save annotations:
# -----------------
ind = which(!is.na(celldf$full.exttype))

celldf$col = major.col[celldf$hcelltype]
celldf$tsp.col = tsp.major.col[celldf$hcelltype]

keep.cols = c("lbl", "U1", "U2", "barcode", "rind", "region", "projid", "is.doublet", 
              "col", "tspcol", "hcluster", "hcelltype", "hsubclass", 
              "major.celltype", "minor.celltype", "neuronal.layer", "inh.subtype", 
              "neuronal.exttype", "full.exttype")  

write.table(celldf[ind,keep.cols], gzfile(paste0(datadir, prefix, '.final.cell_labels.tsv.gz')), 
            quote=F, row.names=F, sep="\t", col.names=T)
write.table(celldf[,keep.cols], gzfile(paste0(datadir, prefix, '.final_with_db.cell_labels.tsv.gz')),
            quote=F, row.names=F, sep="\t", col.names=T)

cellmeta = celldf[ind,keep.cols]
cellmeta.full = celldf[,keep.cols]
save(cellmeta, file=paste0(datadir, prefix, '.final.cell_labels.Rda'))
save(cellmeta.full, file=paste0(datadir, prefix, '.final_with_db.cell_labels.Rda'))

# --------------------------------------------
# Plot all of the updated annotations:
# Plot the points without the doublet clusters
# --------------------------------------------
ind = which(!is.na(celldf$full.exttype))
ind = sample(ind, length(ind), replace=FALSE)
celltype.loc = aggregate(cbind(U1, U2) ~ full.exttype, celldf[ind,], mean)
lbl.loc = aggregate(cbind(U1, U2) ~ lbl, celldf[ind,], mean)
cex = 0.025
tsp.lbl.cols = sapply(lbl.cols, tsp.col)

png(paste0(imgpref, 'umap_fulltypes_', lblset, '_wout_doublets.png'), units='in', res=450, width=8, height=8)
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
with(lbl.loc, text(U1, U2, lbl, cex=.5, xpd=TRUE, font=1))
with(celltype.loc, text(U1, U2, full.exttype, cex=.75, xpd=TRUE, font=2))
mtext(paste(ctype, 'Clusters'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()


# ------------------------------------
# Plot with own colors for full types:
# ------------------------------------
typelvls = unique(celldf$full.exttype)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
tsp.type.cols = sapply(type.cols, tsp.col)

png(paste0(imgpref, 'umap_fulltypes_cols_', lblset, '_wout_doublets.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=tsp.type.cols[celldf$full.exttype[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, lbl, cex=.5, xpd=TRUE, font=1))
with(celltype.loc, text(U1, U2, full.exttype, cex=.75, xpd=TRUE, font=2))
mtext(paste(ctype, 'Clusters'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()

# ------------------------------------------
# Plot with own colors for major cell types:
# ------------------------------------------
png(paste0(imgpref, 'umap_fulltypes_major_cols_', lblset, '_wout_doublets.png'), units='in', res=450, width=8, height=8)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(bsp,bsp,2,sp))
plot(celldf$U1[ind], celldf$U2[ind], col=tsp.major.col[celldf$hcelltype[ind]], 
     pch=19, cex=cex, axes=F)
rect(xleft=par()$usr[1], xright=par()$usr[2],
     ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
     ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
     col='grey85', border=NA, lwd=.5, xpd=TRUE)
with(lbl.loc, text(U1, U2, lbl, cex=.5, xpd=TRUE, font=1))
with(celltype.loc, text(U1, U2, full.exttype, cex=.75, xpd=TRUE, font=2))
mtext(paste(ctype, 'Clusters'), side=3, cex=1.5, col='grey25', font=2, line=0.25)
mtext('UMAP 1', side=1, line=0.25, cex=1.25)
mtext('UMAP 2', side=2, line=0, cex=1.25)
dev.off()


