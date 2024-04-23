#!/usr/bin/R
# ---------------------------------------------------
# Select modules from multiregion set for Sumaiya
# 1. Remove/aggregate the redundant clusters
# 2. Score based on average gene expression per cell type
# - Remove ambient RNA / contam
# - Keep cell type specific
# - Cut down to top N genes
# Updated 10/13/2023
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
options(width=175)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Directories:
moddir = paste0(sdbdir, 'modules/')
regdir = paste0(sdbdir, 'dereg/')
srdir = paste0(sdbdir, 'subtype_reg/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'module_selection_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the module 
# -----------------------------------------------
project = 'snPFC'
if (project == 'multiRegion'){
    source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))
} else {
    graph_id = 'boot'
    cts = c('Ast','Exc','Inh','Mic','Opc','Oli','Vas')
    cmlist = list()
    set_proj('AD430', 'snPFC')
    for (runset in cts){
        # Load in modules annotation:
        commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
        source(paste0(sbindir, 'modules/load_modules_degenr.R'))
        cmlist[[runset]] = coremap
    }
}



# Average expression of each modules in each set:
# -----------------------------------------------
dfll = lapply(names(cmlist), function(x){
    df = data.frame(
        gene=names(cmlist[[x]]),
        module=paste0(x, '-', cmlist[[x]]))
    rownames(df) = NULL
    return(df) })
df = do.call(rbind, dfll)
df$ind = 1

mingenes = 10 # Keep modules with at least 10 genes:
tform = pivot.tomatrix(df, 'module','ind')
tform[is.na(tform)] = 0
marg = apply(tform, 2, sum)
tform = sweep(tform, 2, marg, '/')
tform = tform[,marg >= mingenes]


# Load in the average cell type expression:
avgmat = readRDS(paste0(srdir, 'pseudobulk_data_allcts_averageprofiles.rds'))
ind = grep("@SMC", colnames(avgmat), invert=TRUE)
kept_genes = intersect(rownames(tform), rownames(avgmat))
avgmat = avgmat[kept_genes,ind]
modmat = t(tform[kept_genes,]) %*% avgmat
df = df[df$gene %in% kept_genes,]


# Get modules x expression matrix
# -------------------------------
ctmap = unique(cellmeta[,c('major.celltype', 'cell_type_high_resolution')])
rownames(ctmap) = ctmap$cell_type_high_resolution

module.ct = sub("-.*", "", rownames(modmat))
module.ct[module.ct == 'Mic'] = 'Mic_Immune'
module.ct[module.ct == 'Vas'] = 'Vasc_Epithelia'
ind = which(module.ct != 'HCneurons')
module.ct = module.ct[ind]
modmat = modmat[ind,]

expr.subct = sub(".*@","", colnames(modmat))
expr.ct = sub("/","_", ctmap[expr.subct, 'major.celltype'])

# Plot modules expression:
ux = 1.5
pltmat = sweep(modmat, 1, apply(modmat, 1, max), '/')
ht = Heatmap(
    as.matrix(pltmat),
    column_split=expr.ct,
    row_split=module.ct,
    width = ncol(pltmat)*unit(ux, "mm"), 
    height = nrow(pltmat)*unit(ux, "mm"),
    border_gp = gpar(col="black", lwd = .5)
)

h = 1 + nrow(pltmat) / 15
w = 2 + ncol(pltmat) / 15
pltprefix = paste0(imgpref, 'modexpr_heatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)


# Merge at cell type level and score modules:
# -------------------------------------------
tform = make.tform(expr.ct, norm=T)
ctmodmat = modmat %*% tform
# ctmodmat = ctmodmat[,1:6]
wm = colnames(ctmodmat)[apply(ctmodmat, 1, which.max)]
# good correspondence:

ha = HeatmapAnnotation(match = ifelse(wm == module.ct, 'Yes', 'No'),
    which='row', col=list(match=c('No'='white', 'Yes'='slateblue')))

# Plot modules expression:
ux = 1.5
pltmat = sweep(ctmodmat, 1, apply(ctmodmat, 1, max), '/')
ht = Heatmap(
    as.matrix(pltmat),
    # column_split=expr.ct,
    left_annotation=ha,
    row_split=module.ct,
    width = ncol(pltmat)*unit(ux, "mm"), 
    height = nrow(pltmat)*unit(ux, "mm"),
    border_gp = gpar(col="black", lwd = .5)
)

h = 1 + nrow(pltmat) / 15
w = 2 + ncol(pltmat) / 15
pltprefix = paste0(imgpref, 'modexpr_ct_heatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)


# Cut down modules based on size, etc.
# ------------------------------------
kept.modules = rownames(modmat)[wm == module.ct]
MAXGENES=50 

# Size of these modules:
ngdf = agg.rename(gene ~ module, df, length, 'ngene')
ngdf = ngdf[ngdf$module %in% kept.modules,]
ngdf$ct = sub("-.*", "", ngdf$module)
ngdf$ct[ngdf$ct == 'Mic'] = 'Mic_Immune'
ngdf$ct[ngdf$ct == 'Vas'] = 'Vasc_Epithelia'
sum(ngdf$ngene <= MAXGENES)

# For modules with too many genes:
tform = make.tform(expr.ct, norm=T)
ctavgmat = avgmat %*% tform
seldf = c()
for (module in kept.modules){
    genes = df$gene[df$module == module]
    ct = ngdf$ct[ngdf$module == module]
    ng = ngdf$ngene[ngdf$module == module]

    # Sort by expr
    sub.expr = sort(ctavgmat[genes, ct], decreasing=T)
    top.genes = head(names(sub.expr), MAXGENES)
    cat(module, '\tng:', ng, 'to', length(top.genes), 
        'genes \n', top.genes, '\n\n')
    seldf = rbind(seldf, 
        data.frame(gene=top.genes, module=module, celltype=ct))
}

write.table(seldf, paste0(sdbdir, 'selected_modules_', project, 'AD_101323.tsv'), quote=F, sep="\t", row.names=F)



# 2. Score based on average gene expression per cell type
# - Remove ambient RNA / contam
# - Keep cell type specific
# - Cut down to top N genes


# Read in the module clusters (for multiregion deg)
# ----------------------------
modcls.tsv = paste0(crossdir, 'shared_genes_module_clusters.tsv')
cls.tsv = paste0(crossdir, 'module_cluster_assignments.tsv')
clsdf = read.delim(modcls.tsv, header=T)
clsassign = read.delim(cls.tsv, header=T)
clsdf$cls = paste0('C', clsdf$cls)
acls = unique(clsdf$cls)


# Flag modules by clusters:
# -------------------------
# TODO: Keep core genes from these clusters as one module:
# C10 is cholesterol, merge?
flag.clust = c('C5', 'C9', 'C10', 'C18', 'C19', 'C20')
rm.clust = 'C16' # Ambient - TODO: CHECK

# 1. Remove/aggregate the redundant clusters
# 2. Score based on average gene expression per cell type
# - Remove ambient RNA / contam
# - Keep cell type specific
# - Cut down to top N genes








