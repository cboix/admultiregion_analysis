#!/usr/bin/R
# -----------------------------------------------------
# Plot top marker overlaps between multiple cell types 
# when comparing thalamus <> neocortex subtypes
# Updated 06/29/2023
# -----------------------------------------------------
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
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'markers/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Compute marker differences (pseudobulk) between 
# neocortex and thalamus within each major cell type:
# ---------------------------------------------------
sets = unique(cellmeta$major.celltype)
fulldf = c()
for (subset in sets){
    print(subset)
    ststr = gsub("/","_", subset)
    source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
    if (!file.exists(psdata.rda)){
        ps.data = load_pseudobulk_dataset(subset, subtypes, reg.nomb)
        save(ps.data, file=psdata.rda)
    } else { load(psdata.rda) }

    # Merge to region-level comparison between thalamus and neocortex:
    ps.data$meta$region[ps.data$meta$region %in% c('AG', 'MT', 'PFC')] = 'Neocortex'
    ps.data$meta = ps.data$meta[ps.data$meta$region %in% c('Neocortex', 'TH'),]
    ps.data$mat = ps.data$mat[, ps.data$meta$ptype]
    ps.data$meta$pr = with(ps.data$meta, paste0(projid, '-', region))
    tform = make.tform(ps.data$meta$pr, norm=FALSE)
    tform = sweep(tform, 1, ps.data$meta$ncell, '*')
    tform = sweep(tform, 2, apply(tform, 2, sum), '/')

    comp.mat = ps.data$mat %*% tform 
    comp.df = aggregate(ncell ~ pr + projid + region, ps.data$meta, sum)
    comp.mat = comp.mat[, comp.df$pr]

    # Simple test to rank genes:
    expr.cut = 1
    mean.expr = apply(comp.mat, 1, mean)
    comp.mat = comp.mat[mean.expr > expr.cut,]
    ind.neo = comp.df$region == 'Neocortex'
    ind.th = comp.df$region == 'TH'
    pvals = rep(1, nrow(comp.mat))
    for (i in 1:nrow(comp.mat)){
        gene = rownames(comp.mat)[i]
        x.neo = comp.mat[gene, ind.neo]
        x.th = comp.mat[gene, ind.th]
        wt = wilcox.test(x.neo, x.th)
        pvals[i] = wt$p.value
    }

    rankdf = data.frame(gene=rownames(comp.mat), p=pvals, set=subset)
    rankdf$mn.neo = apply(comp.mat[,ind.neo], 1, mean)
    rankdf$mn.th = apply(comp.mat[,ind.th], 1, mean)
    rankdf$padj = p.adjust(rankdf$p, 'BH')
    rankdf$sig = with(rankdf, ifelse(padj < 0.05, ifelse(mn.neo > mn.th, 'Neocortex', 'Thalamus'), 'NS'))
    rankdf$gset = paste0(rankdf$set, '-', rankdf$sig)
    rankdf = rankdf[order(rankdf$p),]
    fulldf = rbind(fulldf, rankdf)
}


# Compute overlaps for the (signficant) marker gene sets:
# -------------------------------------------------------
fulldf$logFC = log2(fulldf$mn.neo / fulldf$mn.th)
fulldf$abs.logFC = abs(fulldf$logFC)
fulldf$log10p = -log10(fulldf$p)

NS = length(sets)
ovl.mat = matrix(0, nrow=NS, ncol=NS, dimnames=list(sets, sets))
jacc.mat = ovl.mat 
for (set.neo in sets){
    gset.neo = fulldf$gene[(fulldf$set == set.neo) & (fulldf$sig == 'Neocortex')]
    for (set.th in sets){
        gset.th = fulldf$gene[(fulldf$set == set.th) & (fulldf$sig == 'Thalamus')]
        ovl.set = intersect(gset.neo, gset.th)
        ovl.mat[set.neo, set.th] = length(ovl.set)
        jacc.mat[set.neo, set.th] = length(ovl.set) / length(union(gset.neo, gset.th))
    } 
}



# Plot the number of overlapping genes:
# -------------------------------------
pltmat = ovl.mat - diag(NA *diag(ovl.mat))
mx = max(pltmat, na.rm=T)
col_fun = colorRamp2(c(0, mx), c("white", "red"))
ux = 1.5
plt = Heatmap(pltmat, name='# of overlapping\nmarker genes',
    use_raster=TRUE,
    col=col_fun,
    # column_split=csplit,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(pltmat)*unit(ux * 2, "mm"), 
    height = nrow(pltmat)*unit(ux, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
        ann = abs(pltmat[i,j])
        if (!is.na(ann)){
        grid.text(ann, x, y, gp=gpar(fontsize=5)) 
        }
    })

h = 1 + 1 / 15 * nrow(pltmat)
w = 1.5 + 1 / 15 * ncol(pltmat) * 2
pltprefix = paste0(imgpref, 'overlap_th_neocortex_separation_heatmap')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Plot jaccard too:
# -----------------
# mx = max(jacc.mat)
# col_fun = colorRamp2(c(0, mx), c("white", "red"))
pltmat = jacc.mat - diag(NA *diag(jacc.mat))
ux = 1.5
plt = Heatmap(pltmat, name='Jaccard (marker\ngene overlap)',
    use_raster=FALSE,
    # col=col_fun,
    # column_split=csplit,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(pltmat)*unit(ux, "mm"), 
    height = nrow(pltmat)*unit(ux, "mm")
    )

h = 1 + 1 / 15 * nrow(pltmat)
w = 1.5 + 1 / 15 * ncol(pltmat)
pltprefix = paste0(imgpref, 'jaccard_th_neocortex_separation_heatmap')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Report the top overlapping genes
# for top pairs of overlaps:
# --------------------------------
jdf = data.frame(jacc.mat, check.names=FALSE)
jdf$set.neo = rownames(jdf)
jdf = gather(jdf, set.th, jaccard, -set.neo)
jdf = jdf[order(jdf$jaccard, decreasing=TRUE),]
jdf$genes = ''
for (i in 1:nrow(jdf)){
    subdf = fulldf[(fulldf$set == jdf$set.neo[i]) & (fulldf$sig == 'Neocortex'),]

    gset.neo = fulldf$gene[(fulldf$set == jdf$set.neo[i]) & (fulldf$sig == 'Neocortex')]
    gset.th = fulldf$gene[(fulldf$set == jdf$set.th[i]) & (fulldf$sig == 'Thalamus')]
    ovl.set = intersect(gset.neo, gset.th)
    jdf$genes[i] = paste(head(ovl.set), collapse=',')

    # TODO: Make sure genes in order?
}



