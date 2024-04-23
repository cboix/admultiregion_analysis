#!/usr/bin/R
# ---------------------------------------------------------
# Calculate and plot the module enrichments for covariates:
# Updated 03/18/2022 (standardized figure size)
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


library(tidyr)
library(viridis)
library(ggpubr)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Load in the full pseudobulk data for these subtypes:
# NOTE: From modules_degenr_03_plot_pseudobulk_modules.R
# ------------------------------------------------------
useset = 'coregenes_'
scores.file = paste0(moddir, 'module_pseudobulk_scores_', useset, fullpref, '.tsv.gz')
scdf = read.delim(gzfile(scores.file), sep="\t")

# For All, aggregate at major celltype level:
if (runset == 'All'){
    ctmap = unique(cellmeta[,c('cell_type_high_resolution','major.celltype')])
    scdf = merge(scdf, ctmap)
    scdf$totscore = scdf$score * scdf$ncell
    runscdf = aggregate(cbind(totscore, ncell) ~ major.celltype + projid + region + mname, scdf, sum)
    runscdf$score = runscdf$totscore / runscdf$ncell
    runscdf$ptype = with(runscdf, paste0(projid, '_', major.celltype, '_', region))
    scdf = runscdf
    cls = 'major.celltype'
    ct.cols = major.col
} else {
    cls = 'cell_type_high_resolution'
    ct.cols = tcols[plt.subtypes]
}


# Turn back into separate matrix and metadata:
# --------------------------------------------
mod.mat = pivot.tomatrix(scdf[,c('mname','ptype','score')], 'ptype','score')
umeta = unique(scdf[,c('ptype','projid','region',cls,'ncell')])
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(mod.mat),]


# Select only modules with a minimum number of genes:
# ---------------------------------------------------
ngdf = aggregate(gene ~ leiden, nodedf, length)
names(ngdf) = c('module','ngene')
ngdf = merge(ngdf, mmap)

mingenes = 10
kept.modules = ngdf$mname[ngdf$ngene >= mingenes]
kept.modules = rownames(mod.mat)[rownames(mod.mat) %in% kept.modules]
mod.mat = mod.mat[kept.modules,]


# Make annotation from the metadata table:
# ----------------------------------------
ux = 1.5
fixwidth = 125

full.projids = sort(unique(cellmeta$projid))
projid.cols = snap.cols[1:48]
names(projid.cols) = full.projids
nc.col_fun = colorRamp2(range(umeta$ncell), c("white", "indianred")) 

# Make annotation from pseudobulk data:
clsplit = umeta$nrad
ha = HeatmapAnnotation(CT=umeta[[cls]], 
                       Region=umeta$region,
                       Braak=umeta$braaksc,
                       nrad=umeta$nrad,
                       cogdxad=umeta$cogdxad,
                       ncell=umeta$ncell,
                       e4=umeta$Apoe_e4,
                       projid=as.character(umeta$projid),
                       annotation_name_gp = gpar(fontsize=5),
                       simple_anno_size = unit(ux / 1.25, 'mm'),
                       gap = unit(0, "mm"),
                       col=list(CT=ct.cols,
                                Region=reg.cols,
                                ncell=nc.col_fun,
                                Braak=colvals[['braaksc']],
                                nrad=colvals[['nrad']],
                                cogdxad=colvals[['cogdxad']],
                                projid=projid.cols,
                                e4=c('no'='grey90','yes'='slateblue')
                                ))


# Load the annotation table (enrichments in covariates):
# NOTE: from modules_degenr_04_metadata_enrichments.R
# ------------------------------------------------------
useset = 'coregenes_'  # Score the module by the core genes only.
hgann.file = paste0(moddir, 'module_covariate_hgenr_', useset, fullpref, '.tsv')
labdf = read.delim(hgann.file, sep="\t")

# Reorder as matrix:
rownames(labdf) = labdf$mname
labdf = labdf[rownames(mod.mat),]

# Make annotation from covariate enrichments:
hb = rowAnnotation(CT=labdf[[cls]], 
    region=labdf$region,
    braak=labdf$braaksc,
    cogdx=labdf$cogdx,
    nrad=labdf$nrad,
    cogdxad=labdf$cogdxad,
    e4=labdf$Apoe_e4,
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(ux / 1.25, 'mm'),
    gap = unit(0, "mm"),
    col=list(CT=ct.cols,
        region=reg.cols,
        braak=colvals[['braaksc']],
        nrad=colvals[['nrad']],
        cogdxad=colvals[['cogdxad']],
        cogdx=colvals[['cogdx']],
        e4=c('no'='grey90','yes'='slateblue')
        ), 
    na_col='white')



# Score all modules for (a) all genes and (b) tested DE genes:
# ------------------------------------------------------------
pltmat = log1p(mod.mat)[, umeta$ptype]
pltmat = t(scale(t(pltmat)))
plt = Heatmap(pltmat, 
              name='Module\nscore\n(scaled)', 
              use_raster=TRUE,
              top_annotation=ha, 
              column_split=clsplit, 
              width=unit(fixwidth, 'mm'),  # Fixed width
              height=nrow(pltmat) * unit(ux, 'mm'),
              border_gp=gpar(lwd=.5, color='black'),
              row_dend_width = unit(ux * 2, "mm"),
              column_dend_height = unit(ux * 2, "mm"),
              show_column_names=FALSE,
              right_annotation=hb)


h = 5 + nrow(pltmat) / 15
w = 5 + (fixwidth / ux) / 15
pltprefix = paste0(imgpref, 'module_pseudobulk_scaled_annotated_', fullpref)
saveHeatmap(plt, pltprefix, w=w, h=h)

