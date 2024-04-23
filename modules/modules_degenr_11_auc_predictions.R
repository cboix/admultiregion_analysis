#!/usr/bin/R
# --------------------------------------------------
# Ask whether modules predict AD and other variables 
# at sample & cell type level
# Updated 01/26/2022
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))

library(tidyr)
library(viridis)
library(PRROC)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
resdir = paste0(sdbdir, 'modules/resources/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'modules_predAD_')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Set the run arguments:
# ----------------------
# TODO: Load across all?
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Heatmap plotting functions:
# ---------------------------
plotEffSizeHeatmap = function(cmat, pmat, ux, 
                              col.split=NULL, row.split=NULL, 
                              cluster=TRUE, topann=NULL){
    plt = Heatmap(cmat,
                  col=col_fun,
                  use_raster=TRUE,
                  column_split=col.split,
                  row_split=row.split,
                  cluster_columns=cluster,
                  cluster_rows=cluster,
                  cluster_row_slices=cluster,
                  top_annotation=topann,
                  width = ncol(cmat)*unit(ux, "mm"), 
                  height = nrow(cmat)*unit(ux, "mm"),
                  border_gp = gpar(col="black", lty = 1, lwd=.5),
                  cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                      p = pmat[i,j]
                      if (p < 0.05){ grid.text('*', x, y, gp=gpar(fontsize=gridtxt.fs))} }
    )
    return(plt)
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# TODO: Load pseudobulk data, and compute these scores just on DE genes?
# Load in the full pseudobulk data for these subtypes:
# ----------------------------------------------------
psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
if (!file.exists(psdata.rda)){
    ps.data = load_pseudobulk_dataset(celltype, subtypes, region.set)
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }

umeta = ps.data$meta
umeta = unique(umeta)
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(ps.data$mat),]
# Remove very low abundance batches + will use for weight
umeta = umeta[umeta$ncell > 10,] 


# Filter modules by min genes:
# ----------------------------
ngdf = aggregate(gene ~ leiden, nodedf, length)
names(ngdf) = c('module','ng')
mmap = merge(mmap, ngdf, all.x=TRUE)
mingenes = 10
kept.mnames = mmap$mname[mmap$ng >= mingenes]


# Score all modules for (a) core genes and (b) tested DE genes:
# ------------------------------------------------------------
subset.de = FALSE
modules = sort(unique(coremap))
kept.genes = rownames(ps.data$mat)
kept.genes = kept.genes[kept.genes %in% names(coremap)]
if (subset.de){
    advar = 'cogdxad'
    subdedf = dedf[dedf$dkey == advar,]
    rownames(subdedf) = subdedf$gene
    de.genes = subdedf$gene[subdedf$gset != '--'] 
    # de.genes = dedf$gene[(dedf$dkey == advar) & (dedf$gset == 'Up')]
    kept.genes = kept.genes[kept.genes %in% de.genes]
    kept.eff = subdedf[kept.genes, 'logFC_nb']
}

# Conversion matrix:
tform = make.tform(coremap[kept.genes], u=modules, norm=TRUE)
if (subset.de){
    mod.mat = t(tform) %*% sweep(ps.data$mat[kept.genes,], 1, kept.eff, '*')
} else {
    mod.mat = t(tform) %*% ps.data$mat[kept.genes,]
}
mod.mat = as.matrix(mod.mat)
rownames(mod.mat) = mmap$mname[1:nrow(mod.mat)]

# Remove non-scored (if using DE genes):
kept.mod = which(!is.na(rowSums(mod.mat)))
mod.mat = mod.mat[kept.mod,]


# Aggregate at the individual + region level as well:
# ---------------------------------------------------
mod.data = list('mat'=mod.mat[,umeta$ptype], 'meta'=umeta)
if (runset == 'Mic_Immune'){
    # Subset to microglia and aggregate:
    mod.data.mic = mod.data
    mind = grep("^Mic", mod.data.mic$meta$cell_type_high_resolution)
    mod.data.mic$meta = mod.data.mic$meta[mind,]
    mod.data.mic$mat = mod.data.mic$mat[, mod.data.mic$meta$ptype]
    ind.data = aggregatePsbulkIndRegion(mod.data.mic)
} else {
    ind.data = aggregatePsbulkIndRegion(mod.data)
}


# Perform AUC calculations for each module alone:
# ----------------------------------------------
# TODO: For multiple (binary) classes?
advar = 'cogdxad'
oracle = as.numeric(ind.data$meta[[advar]]) - 1

NM = nrow(ind.data$mat)
AUCstat = rep(0, NM)
AUPRCstat = rep(0, NM)
for (i in 1:NM){
    x = as.numeric(ind.data$mat[i,])
    fg <- x[oracle == 1]
    bg <- x[oracle == 0]
    # Flip if reverse:
    if (mean(fg) < mean(bg)){
        fg <- x[oracle == 0]
        bg <- x[oracle == 1]
    }
    # ROC Curve    
    roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    AUPRCstat[i] = pr$auc.integral
    AUCstat[i] = roc$auc
}
audf = data.frame(mname=rownames(ind.data$mat),
    AUC=AUCstat, AUPRC=AUPRCstat)
audf = audf[order(audf$AUPRC, decreasing=T),]


# Look at R.sq increase from combining scores:
# TODO: Also look at AIC / BIC and deviance.
# --------------------------------------------
# First run each module separately:
# scdf = ind.data$meta
scdf = mod.data$meta
scdf = merge(scdf, unique(metadata[,c('projid','region','rind')]))
scdf = merge(scdf, pqdf, all.x=TRUE)
scdf$cogdxad = factor(scdf$cogdxad, levels=c('CTRL','AD'))
scdf$is.ad = (scdf$cogdxad == 'AD')
scdf$nft = log1p(scdf$nft)
scdf$plaq_n = log1p(scdf$plaq_n)
scdf$plaq_d = log1p(scdf$plaq_d)

resdf = c()
for (path in c('is.ad','nft','plaq_n','plaq_d')){
    for (i in 1:NM){
        # mn = rownames(ind.data$mat)[i]
        # scdf$x = ind.data$mat[i,scdf$pr]
        mn = rownames(mod.data$mat)[i]
        scdf$x = mod.data$mat[i,scdf$ptype]
        # form = asform(c(path, '~ x * region')) # TODO add mixed effects
        # form = asform(c(path, '~ x * cell_type_high_resolution')) # TODO add mixed effects
        # form = asform(c(path, '~ x + cell_type_high_resolution * region')) # TODO add mixed effects
        form = asform(c('x ~', path, '+ cell_type_high_resolution * region')) # TODO add mixed effects
        fit = lm(form, scdf, weights=log(scdf$ncell + 1))
        cfit = data.frame(coefficients(summary(fit)))
        colnames(cfit) = c('Est','SE','t','p')
        resdf = rbind(resdf, data.frame(
                mname = mn, path = path,
                est = cfit[2, 'Est'], p = cfit[2, 'p'],
                rsq = summary(fit)$r.squared))
    }
}
resdf = resdf[order(resdf$rsq, decreasing=T),]
resdf$rsq = abs(resdf$rsq) * sign(resdf$est)

# Save these results for collating later:
rsq.file = paste0(moddir, 'rsq_vs_pathology_measures_', runset, '.rda')
rsq.tsv = paste0(moddir, 'rsq_vs_pathology_measures_', runset, '.tsv.gz')
saveRDS(resdf, file=rsq.file)
write.table(resdf, file=gzfile(rsq.tsv), quote=F, row.names=F, sep="\t")


# Plot the r-squared effect sizes as a heatmap:
# ---------------------------------------------
resdf = resdf[resdf$mname %in% kept.mnames,]
resdf$p.adj = p.adjust(resdf$p, 'BH')
cmat = pivot.tomatrix(resdf[,c('mname','path','est')], 'path','est')
pmat = pivot.tomatrix(resdf[,c('mname','path','p.adj')], 'path','p.adj')

cn = rownames(reord(t(abs(cmat))))
zmat = 1 * (abs(cmat[,cn]) * (pmat[, cn] < 0.05)) * 200
ll = diag.mat2(t(zmat[, cn]))
rn = rev(ll[[2]])

mx = 0.1
col_fun = colorRamp2(c(-mx, 0, mx), c('blue', "white", 'red'))
ht = plotEffSizeHeatmap(cmat[rn, cn], pmat[rn, cn], ux=1.5, cluster=FALSE)
h = 3 + 1 / 15 * nrow(cmat)
w = 5 + 1 / 15 * ncol(cmat)
pltprefix = paste0(imgpref, 'rsq_vs_pathology_measures_', runset, '_heatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)

# Signif only:
ind = which(apply(pmat < 0.05, 1, sum) > 0)

cn = rownames(reord(t(abs(cmat))))
zmat = 1 * (abs(cmat[,cn]) * (pmat[, cn] < 0.05)) * 200
ll = diag.mat2(t(zmat[ind, cn]), ratio=1.1)
rn = rev(ll[[2]])


ht = plotEffSizeHeatmap(cmat[ind,], pmat[ind,], ux=1.5)
h = 3 + 1 / 15 * length(ind)
pltprefix = paste0(imgpref, 'rsq_vs_pathology_measures_', runset, '_sig_heatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)



# Joint predictions:
# ------------------

kept.mod = c(1, 11, 20)
kept.mod = 1:27 - 1
kept.names = mmap[mmap$module %in% kept.mod,'mname']
scdf = data.frame(t(ind.data$mat[kept.names,]))
colnames(scdf) = paste0('M', kept.mod)
scdf = cbind(scdf, ind.data$meta)
covars = '+msex + Apoe_e4 + region'
form = asform(c('is.ad~', paste(paste0('M', kept.mod), collapse='+'), covars)) # TODO add mixed effect
scdf$is.ad = scdf$cogdxad == 'AD'

fit = glm(form, scdf, family='gaussian')
cfit = data.frame(coefficients(summary(fit)))
colnames(cfit) = c('Est','SE','t','p')
cfit$module = sub("M", "", rownames(cfit))
cfit = merge(cfit, mmap, all.x=TRUE)
cfit = cfit[order(cfit$p),]

with(summary(fit), 1 - deviance/null.deviance)



fit.all = lm(is.ad ~ M11 + region, scdf)
summary(fit.all)$r.squared

fit.all = lm(is.ad ~ M1 + M11 + M20 + region, scdf)
summary(fit.all)$r.squared



# Joint score:
score = rep(0, ncol(ind.data$mat))
for (name in kept.names){
    x = as.numeric(ind.data$mat[i,])
    fg <- x[oracle == 1]
    bg <- x[oracle == 0]
    # Flip if reverse:
    if (mean(fg) < mean(bg)){
        score = score + scale(-x)
    } else {
        score = score + scale(x)
    }
}

fg <- score[oracle == 1]
bg <- score[oracle == 0]

# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

AUPRCstat[i] = pr$auc.integral
AUCstat[i] = roc$auc
audf = data.frame(mname=rownames(ind.data$mat),
    AUC=AUCstat, AUPRC=AUPRCstat)

# For microglia, show specified modules together


# For astrocytes, also show specified modules together


