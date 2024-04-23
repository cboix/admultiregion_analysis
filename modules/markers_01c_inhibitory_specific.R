#!/usr/bin/R
# ----------------------------------------------------------
# Plot basic marker differences between inhibitory subtypes:
# Updated 12/03/2021
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(viridis)
# library(ggplot2)
# library(ggrepel)
# library(ggpubr)
# library(ggpmisc)
# library(patchwork)

# library(ComplexHeatmap)
# library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'markers/')
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
remove.batches = TRUE
suff = '_subset_final_noMB'

subset = 'Inh'
imgpref = paste0(plotdir, subset, '_')
ststr = gsub("/","_", subset)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
subtypes = unique(ctdf$cls)
celltypes = unique(ctdf$major.celltype)


# Load in the full ast data for these subtypes:
# ---------------------------------------------
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
if (!file.exists(psdata.rda)){
    ps.data = load_pseudobulk_dataset(subset, subtypes, reg.nomb)
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }


# Load in the specific UMAPs for Inhibitory:
# ------------------------------------------
celltype = 'Inh'
cellsuff = paste0(celltype, '_log1p_combat_filthvg1000_nomt')
ctumap.file = paste0(sdbdir, 'metadata/multiregion_majorctUMAP_', cellsuff, '.tsv')

submeta = read.delim(ctumap.file, sep="\t")
submeta = submeta[, c('barcode','C1','C2')]
names(submeta) = c('barcode','U1','U2')
rownames(cellmeta) = cellmeta$barcode
submeta$region = cellmeta[submeta$barcode, 'region']
submeta$projid = cellmeta[submeta$barcode, 'projid']
submeta$cell_type_high_resolution = cellmeta[submeta$barcode, 'cell_type_high_resolution']
rownames(submeta) = submeta$barcode


# Get the non-pseudobulk data for plotting on a cell-level:
# ---------------------------------------------------------
nmat = load_full_dataset(celltypes, subtypes, reg.nomb)
submeta = submeta[colnames(nmat),]


# Plot some of these genes on the Inhibitory-only UMAP:
# -----------------------------------------------------
ind = 1:nrow(submeta)
cex = 0.01
xlim = range(submeta$U1)
ylim = range(submeta$U2)

palette = viridis(100)
col_fun = function(x, pal=palette){
                bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
                palette[bin]  }

pgenes = c('GRM3','MEIS2','FOXP2','CHRM2','LAMP5','SST','PVALB','PAX6','VIP')
for (gene in pgenes){
    print(gene)
    x = log(nmat[gene,submeta$barcode] + 1)
    ind = order(x)
    png(paste0(imgpref, 'umap_gene_', gene,'_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=3, height=3)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    par(mar=rep(sp, 4))
    plot(submeta$U1[ind], submeta$U2[ind], col=col_fun(x[ind]), ylim=ylim, xlim=xlim, pch=19, cex=cex, axes=F)
    mtext(gene, side=1, line=-1, font=2, cex=2)
    dev.off()
}

# TODO: replot without text

