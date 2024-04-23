#!/usr/bin/R
# ----------------------------------
# Put together the full GWAS figure:
# Updated 02/21/2022
# ----------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'gwas/')
imgpref = paste0(plotdir, 'gwasfull_')
cmd = paste('mkdir -p', plotdir, srdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Heatmap plotting functions:
# ---------------------------
plotDEgenesHeatmap = function(cmat, pmat, ux, col.split=NULL, row.split=NULL, cluster=TRUE, topann=NULL){
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
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lty = 1, lwd=.5),
        cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
            p = pmat[i,j]
            if (p < 0.05){ grid.text('*', x, y, gp=gpar(fontsize=gridtxt.fs))} }
    )
    return(plt)
}


plotExprHeatmap = function(emat, ux, row.split=NULL, cluster=FALSE, raster=TRUE, topann=NULL){
    plt = Heatmap(emat,
        col=viridis(50),
        name='expr',
        use_raster=raster,
        column_split=sub("-.*","", colnames(emat)),
        row_split=row.split,
        cluster_columns=cluster,
        cluster_column_slices=cluster,
        cluster_row_slices=cluster,
        cluster_rows=cluster,
        top_annotation=topann,
        width = ncol(emat)*unit(ux / 2, "mm"), 
        height = nrow(emat)*unit(ux, "mm"),
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lty = 1, lwd=.5),
    )
    return(plt)
}


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])


# Load pseudobulk data for all runsets:
# -------------------------------------
gw.ps.rda = paste0(srdir, 'pseudobulk_data_gwasgenes_regionaverages.rda')
load(gw.ps.rda)  # norm.exprmat full.exprmat 
topdf = read.delim('Annotation/ADGWAS_topct_multiregion_121721.tsv', sep="\t")


# Load DEG data for all runsets:
# ------------------------------
# region = 'allregions'
# gw.reg.file = paste0(regdir, 'aggregated_allres.gwgenes.', region, '.rda')


# Get the compiled results for NIA-Reagan regional DEGs:
# ------------------------------------------------------
path = 'nrad'
mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)


# Turn DEGs into a single data frame:
# -----------------------------------
dedf = c()
sets = names(setdflist)
for (set in sets){
    setdf = setdflist[[set]]
    if (!is.null(setdf)){
        setdf$set = set
        dedf = rbind(dedf, setdf)
    }
}
gwdedf = dedf[(dedf$gene %in% gwgenes),]


# Match the DEGs with the expression data:
# ----------------------------------------
path = 'nrad'
col_fun = colorRamp2(c(-.25, 0, .25), c('blue', "white", 'red'))

rnlist = list()
reglist = c('allregions', reg.order[-1])
for (region in reglist){
    subdf = gwdedf[gwdedf$region == region,]
    subdf$lp = ifelse(subdf$col_nm > 0, subdf$log10p_nm, 0)
    cmat = pivot.tomatrix(subdf[,c('gene','set','logFC_nb')], 'set', 'logFC_nb')
    pmat = pivot.tomatrix(subdf[,c('gene','set','lp')], 'set', 'lp')
    pmat = 10**(-pmat)
    cmat[is.na(cmat)] = 0
    pmat[is.na(pmat)] = 1

    full.cmat = matrix(0, nrow=nrow(norm.exprmat), ncol=ncol(cmat),
        dimnames=list(rownames(norm.exprmat), colnames(cmat)))
    full.pmat = matrix(1, nrow=nrow(norm.exprmat), ncol=ncol(pmat),
        dimnames=list(rownames(norm.exprmat), colnames(pmat)))
    full.cmat[rownames(cmat),] = cmat
    full.pmat[rownames(pmat),] = pmat
    prsuff = paste0(region, "_", path)

    # Plot GWAS genes DE tested as heatmap:
    # -------------------------------------
    ux = 2
    cn = colnames(full.cmat)
    col.split = ifelse(1:length(cn) %in% grep("Vasc",cn),'Vasculature','All')
    htde = plotDEgenesHeatmap(full.cmat, full.pmat,
        col.split=col.split, ux=ux, cluster=FALSE)
    h = 2.25 + 1 / 15 * nrow(full.cmat)
    w = 5 + 1 / 15 * ncol(full.cmat)
    pltprefix = paste0(imgpref, 'gwas_deresource_heatmap_', prsuff)
    saveHeatmap(htde, pltprefix, w=w, h=h)

    # Join the GWAS heatmap and the DE heatmap:
    # -----------------------------------------
    htexp = plotExprHeatmap(norm.exprmat, ux=ux, cluster=FALSE)
    ht = htde + htexp
    h = 3 + 1 / 15 * nrow(full.cmat)
    w = 5 + 1 / 15 * (ncol(full.cmat) + ncol(norm.exprmat) / 2)
    pltprefix = paste0(imgpref, 'jointdegwas_allgwgenes_heatmap_', prsuff)
    saveHeatmap(ht, pltprefix, w=w, h=h)

    # Plot significant genes only in the joint heatmaps:
    # --------------------------------------------------
    NTOP = 15
    nonve = grep("Vasc",colnames(full.cmat), invert=TRUE)
    issig = which(apply(full.pmat[, nonve], 1, min) < 0.05)
    icmat = full.cmat[issig, nonve]
    ipmat = full.pmat[issig, nonve]
    rmat = t(reord(t(icmat)))
    cn = colnames(rmat)

    zmat = 1 * (abs(icmat[,cn]) * (ipmat[, cn] < 0.05))
    zmarg = apply(zmat, 1, sum)
    thresh = tail(head(sort(zmarg, decreasing=T), NTOP),1)
    rn = names(zmarg[zmarg >= thresh])
    ll = diag.mat2(t(zmat[rn, cn]))
    rn = rev(ll[[2]])
    col.split = ifelse(1:length(cn) %in% grep("Vasc",cn),'Vasculature','All')
    rnlist[[region]] = rn  # Important for full plot later

    ux = 1.5
    htde = plotDEgenesHeatmap(icmat[rn, cn], ipmat[rn, cn],
        col.split=col.split, ux=ux, cluster=FALSE)

    htexp = plotExprHeatmap(norm.exprmat[rn,], ux=ux, cluster=FALSE)
    ht = htde + htexp
    h = 3 + 1 / 15 * length(rn)
    w = 5 + 1 / 15 * (ncol(full.cmat) + ncol(norm.exprmat) / 2)
    pltprefix = paste0(imgpref, 'jointdegwas_siggenes_heatmap_', prsuff)
    saveHeatmap(ht, pltprefix, w=w, h=h)
}


# Plot a merged version now:
# --------------------------
# Make matrices:
# subdf = full.gwdedf[full.gwdedf$region == region,]
subdf = gwdedf[gwdedf$path == path,]
subdf = subdf[subdf$region %in% reglist,]
subdf$set2 = paste0(subdf$set, '@', subdf$region)

subdf$lp = ifelse(subdf$col_nm > 0, subdf$log10p_nm, 0)
cmat = pivot.tomatrix(subdf[,c('gene','set2','logFC_nb')], 'set2', 'logFC_nb')
pmat = pivot.tomatrix(subdf[,c('gene','set2','lp')], 'set2', 'lp')
pmat = 10**(-pmat)
cmat[is.na(cmat)] = 0
pmat[is.na(pmat)] = 1

sort(apply(pmat < 0.05, 1, sum))
sum(apply(pmat < 0.05, 1, sum) > 0)
topdf = aggregate(lp ~ gene + set, subdf, max)
topdf = topdf[topdf$lp > -log10(0.05),]
topdf = topdf[grep("Vasc", topdf$set, invert=T),]
topdf = topdf[order(topdf$gene),]
tp = sort(table(topdf$gene))
sort(table(topdf$set))


full.cmat = matrix(0, nrow=nrow(norm.exprmat), ncol=ncol(cmat),
    dimnames=list(rownames(norm.exprmat), colnames(cmat)))
full.pmat = matrix(1, nrow=nrow(norm.exprmat), ncol=ncol(pmat),
    dimnames=list(rownames(norm.exprmat), colnames(pmat)))
full.cmat[rownames(cmat),] = cmat
full.pmat[rownames(pmat),] = pmat

# Margins and reordering:
colmarg = apply(full.pmat < 0.05, 2, sum)
dwmarg = apply((full.pmat < 0.05) * (full.cmat < 0), 2, sum)
upmarg = apply((full.pmat < 0.05) * (full.cmat > 0), 2, sum)
mmarg = cbind(dwmarg, upmarg)

sets = unique(sub("@.*", "", colnames(full.cmat)))
sdf = expand.grid(region=reglist, set=sets)
cn = paste0(sdf$set, "@", sdf$region)
cn = cn[cn %in% colnames(full.cmat)]
full.cmat = full.cmat[,cn]
full.pmat = full.pmat[,cn]

# Margin for top expr:
wind = apply(norm.exprmat, 1, which.max)
emarg = data.frame(table(colnames(norm.exprmat)[wind]))
emarg = merge(emarg, data.frame(Var1=colnames(norm.exprmat)), all.y=TRUE)
emarg$Freq[is.na(emarg$Freq)] = 0
rownames(emarg) = emarg$Var1
emarg = emarg[colnames(norm.exprmat),]


# Subset to top genes and reorder matrices:
# -----------------------------------------
top.genes = unique(do.call(c, rnlist))
full.cmat = full.cmat[top.genes,]
full.pmat = full.pmat[top.genes,]

nonve = grep("Vasc",colnames(full.cmat), invert=TRUE)
icmat = full.cmat[, nonve]
ipmat = full.pmat[, nonve]
# rmat = t(reord(t(icmat)))
rmat = icmat
cn = colnames(rmat)

zmat = 1 * (abs(icmat[,cn]) * (ipmat[, cn] < 0.05))
ll = diag.mat2(t(zmat[, cn]))
row.split = sub("@.*", "", cn[ll[[3]]])
rn = rev(ll[[2]])
col.split = gsub("_", "\n", sub("@.*","", cn))

htbr = HeatmapAnnotation(NDE=anno_barplot(mmarg[cn,], bar_width=1, 
        height=unit(.5, 'cm'), border=FALSE, gp=gpar(fill=c(colrb[90], colrb[10]), lty=0)))

httopexpr = HeatmapAnnotation(NTOP=anno_barplot(emarg[colnames(norm.exprmat),'Freq'], bar_width=1, 
        height=unit(.5, 'cm'), border=FALSE, gp=gpar(fill='grey70', lty=0)))


# Plot heatmaps:
# --------------
ux = 1.5
htde = plotDEgenesHeatmap(icmat[rn, cn], ipmat[rn, cn],
    col.split=col.split, ux=ux, cluster=FALSE, topann=htbr)

htexp = plotExprHeatmap(norm.exprmat[rn,], 
    ux=ux, cluster=FALSE, raster=TRUE, topann=httopexpr)
ht = htde + htexp
h = 3 + 1 / 15 * length(rn)
w = 5 + 1 / 15 * (ncol(full.cmat) + ncol(norm.exprmat) / 2)
pltprefix = paste0(imgpref, 'jointdegwas_siggenes_heatmap_region')
saveHeatmap(ht, pltprefix, w=w, h=h)


# Make the full GWAS heatmap for extended data:
# ---------------------------------------------
full.cmat = matrix(0, nrow=nrow(norm.exprmat), ncol=ncol(cmat),
    dimnames=list(rownames(norm.exprmat), colnames(cmat)))
full.pmat = matrix(1, nrow=nrow(norm.exprmat), ncol=ncol(pmat),
    dimnames=list(rownames(norm.exprmat), colnames(pmat)))
full.cmat[rownames(cmat),] = cmat
full.pmat[rownames(pmat),] = pmat

sets = unique(sub("@.*", "", colnames(full.cmat)))
sdf = expand.grid(region=reglist, set=sets)
cn = paste0(sdf$set, "@", sdf$region)
cn = cn[cn %in% colnames(full.cmat)]

icmat = full.cmat[,cn]
ipmat = full.pmat[,cn]
cn = colnames(icmat)
zmat = abs(icmat[,cn])
ll = diag.mat2(t(zmat[, cn]))
row.split = sub("@.*", "", cn[ll[[3]]])
rn = rev(ll[[2]])
rn = rownames(norm.exprmat)
col.split = gsub("_", "\n", sub("@.*","", cn))

htbr = HeatmapAnnotation(NDE=anno_barplot(mmarg[cn,], bar_width=1, 
        height=unit(.5, 'cm'), border=FALSE, gp=gpar(fill=c(colrb[90], colrb[10]), lty=0)))
httopexpr = HeatmapAnnotation(NTOP=anno_barplot(emarg$Freq, bar_width=1, 
        height=unit(.5, 'cm'), border=FALSE, gp=gpar(fill='grey70', lty=0)))

ux = 1.5
htde = plotDEgenesHeatmap(icmat[rn,cn], ipmat[rn, cn],
    col.split=col.split, ux=ux, cluster=FALSE, topann=htbr)

htexp = plotExprHeatmap(norm.exprmat[rn,], 
    ux=ux, cluster=FALSE, raster=TRUE, topann=httopexpr)
ht = htde + htexp
h = 3 + 1 / 15 * length(rn)
w = 5 + 1 / 15 * (ncol(full.cmat) + ncol(norm.exprmat) / 2)
pltprefix = paste0(imgpref, 'jointdegwas_allgenes_heatmap_region')
saveHeatmap(ht, pltprefix, w=w, h=h)

