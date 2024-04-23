#!/usr/bin/R
# ------------------------------------------------------
# Compare our microglia definitions to Tuddenham markers
# Updated 11/17/2023
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(viridis)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)

library(ComplexHeatmap)
library(circlize)
print(version)
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'metadata/')
imgpref = paste0(plotdir, 'microglia_review_')
cmd = paste('mkdir -p', plotdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))



# Load all microglia data at cell-level:
# --------------------------------------
celltype = 'Mic'
commandArgs = function(x){ c('Mic_Immune', 'Mic', 'allregions')}
source(paste0(bindir, 'multiRegion/load_difftl_data.R'))


# Load the Tuddenham et al markers:
# ---------------------------------
markdf = read.delim(paste0(sdbdir, 'tuddenham_microglia_upmarkers.tsv'))
markdf = markdf[markdf$gene %in% rownames(mat),] # Subset to expr data
types = unique(markdf$up_type)
ll = lapply(types, function(x){
    df = markdf[markdf$up_type == x,]
    df$rank = 1:nrow(df)
    df
})
markdf = do.call(rbind, ll)

# Make unique by ranking, break ties by sum_logFC
markdf = merge(markdf, aggregate(rank ~ gene, markdf, min))
markdf = merge(markdf, aggregate(sum_logFC ~ gene, markdf, max))
markdf = markdf[order(markdf$rank),]

# Top N genes for these substates:
NTOP = 4
ll = lapply(types, function(x){
    head(markdf[markdf$up_type == x,], NTOP)
})
markdf = do.call(rbind, ll)
table(markdf$rank)
kept.genes = markdf$gene

gtmap = markdf$up_type
names(gtmap) = markdf$gene


# Average signal for these genes:
# -------------------------------
marg = colSums(mat)
submat = mat[kept.genes,]
submat = sweep(submat, 2, median(marg) / marg, '*')
submat = log1p(submat)

# Average signals:
submeta = cellmeta[colnames(submat),]
tform = make.tform(submeta$cell_type_high_resolution, norm=TRUE)
avg.mat = t(submat %*% tform) # Average signal by subtype:
pct.mat = t((submat > 0) %*% tform) # Percent expr per subtype
avg.mat = as.matrix(avg.mat)
pct.mat = as.matrix(pct.mat)

# Column split:
colsplit = gtmap[colnames(avg.mat)]
uqcol = length(unique(colsplit))

# Heatmap and bubble plot:
# mx = .5
# col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))

ux = 1.5
pltmat = sweep(avg.mat, 2, apply(avg.mat, 2, max), '/')
ht1 = Heatmap(pltmat, 
    use_raster=FALSE, 
    name='Row-max\nAvg. log1p(expr)',
    col=col1,
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    cluster_column_slices=FALSE,
    column_split=colsplit,
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    width=(3/4 * (uqcol-1) + ncol(avg.mat)) * unit(ux, 'mm'),
    height=nrow(avg.mat) * unit(ux, 'mm'),
)

ht2 = Heatmap(avg.mat, 
    use_raster=FALSE, 
    name='Avg. log1p(expr)',
    col=viridis(100),
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    cluster_column_slices=FALSE,
    column_split=colsplit,
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    width=(3/4 * (uqcol-1) + ncol(avg.mat)) * unit(ux, 'mm'),
    height=nrow(avg.mat) * unit(ux, 'mm'),
)

ht = ht1 %v% ht2

pltprefix = paste0(imgpref, 'dejager_markers_avgexpr_heatmap')
h = 1.5 + 1 / 15 * (nrow(avg.mat) * 2)
w = 1.5 + 1 / 15 * ncol(avg.mat)
saveHeatmap(ht, pltprefix, w=w, h=h)

