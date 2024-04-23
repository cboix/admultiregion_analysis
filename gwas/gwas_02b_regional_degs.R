#!/usr/bin/R
# --------------------------------------------------
# Plot GWAS gene intersections with aggregated DEGs:
# Updated: 12/16/21
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

library(ComplexHeatmap)
library(circlize)
options(width=170)

print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'degwas_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


plotDEgenesHeatmap = function(cmat, pmat, col.split, ux, cluster=TRUE){
    plt = Heatmap(cmat,
        col=col_fun,
        use_raster=TRUE,
        column_split=col.split,
        cluster_columns=cluster,
        cluster_rows=cluster,
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
}


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])


# Get the compiled results for NIA-Reagan regional DEGs:
# ------------------------------------------------------
path = 'nrad'
mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)


# Turn DEGs into a single data frame:
# -----------------------------------
dedf = c()
for (set in names(setdflist)){
    setdf = setdflist[[set]]
    if (!is.null(setdf)){
        setdf$set = set
        dedf = rbind(dedf, setdf)
    }
}

gwdedf = dedf[(dedf$gene %in% gwgenes),]


# Make matrices to plot DE genes as heatmap:
# ------------------------------------------
reglist = c('allregions', reg.order[-1])
for (region in reglist){
    subdf = gwdedf[gwdedf$region == region,]
    subdf$lp = ifelse(subdf$col_nm > 0, subdf$log10p_nm, 0)
    cmat = pivot.tomatrix(subdf[,c('gene','set','logFC_nb')], 'set', 'logFC_nb')
    pmat = pivot.tomatrix(subdf[,c('gene','set','lp')], 'set', 'lp')
    if (path %in% c('nft','plaq_n','plaq_d')){
        col_fun = colorRamp2(c(-.025, 0, .025), c('blue', "white", 'red'))
    } else {
        col_fun = colorRamp2(c(-.25, 0, .25), c('blue', "white", 'red'))
    }

    pmat = 10**(-pmat)
    cmat[is.na(cmat)] = 0
    pmat[is.na(pmat)] = 1


    # Plot GWAS genes DE tested as heatmap:
    # -------------------------------------
    ux = 2
    col.split = ifelse(1:ncol(cmat) %in% grep("Vasc",colnames(cmat)),'Vasculature','All')
    plt = plotDEgenesHeatmap(cmat, pmat, col.split, ux=ux)

    h = 1 + 1 / 15 * nrow(cmat)
    w = 2 + 1 / 15 * ncol(cmat)
    pltprefix = paste0(imgpref, 'gwas_desimple_heatmap_', region, '_', path)
    saveHeatmap(plt, pltprefix, w=w, h=h)


    # Only plot signif:
    # ----------------------------------------------------
    issig = which(apply(pmat, 1, min) < 0.05)
    icmat = cmat[issig,]
    ipmat = pmat[issig,]
    rmat = t(reord(t(icmat)))
    cn = colnames(rmat)
    zmat = 1 * (abs(icmat[,cn]) * (ipmat[, cn] < 0.05))
    ll = diag.mat2(t(zmat))
    rn = rev(ll[[2]])

    plt = plotDEgenesHeatmap(icmat[rn, cn], ipmat[rn, cn],
        col.split=NULL, ux=ux, cluster=FALSE)

    h = 1 + 1 / 15 * nrow(icmat[rn, cn])
    w = 2 + 1 / 15 * ncol(icmat[rn, cn])
    pltprefix = paste0(imgpref, 'gwas_desimple_heatmap_signif_', region, '_', path)
    saveHeatmap(plt, pltprefix, w=w, h=h)


    # Only plot the non-vasculature cells and only signif:
    # ----------------------------------------------------
    nonve = grep("Vasc",colnames(cmat), invert=TRUE)
    issig = which(apply(pmat[, nonve], 1, min) < 0.05)
    icmat = cmat[issig, nonve]
    ipmat = pmat[issig, nonve]
    rmat = t(reord(t(icmat)))
    cn = colnames(rmat)
    zmat = 1 * (abs(icmat[,cn]) * (ipmat[, cn] < 0.05))
    zmarg = apply(zmat, 1, sum)
    # zmarg = apply(abs(icmat), 1, sum)
    thresh = sort(zmarg, decreasing=T)[20]
    rn = names(zmarg[zmarg >= thresh])
    ll = diag.mat2(t(zmat[rn, cn]))
    rn = rev(ll[[2]])

    plt = plotDEgenesHeatmap(icmat[rn, cn], ipmat[rn, cn], col.split=NULL, ux=ux, cluster=FALSE)

    h = 1 + 1 / 15 * nrow(icmat)
    w = 2 + 1 / 15 * ncol(icmat)
    pltprefix = paste0(imgpref, 'gwas_desimple_heatmap_nonve_', region, '_', path)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}



# Get the compiled results for all regions DE runs (03):
# ------------------------------------------------------
for (region in c('allregions', reg.nomb)){
    full.file = paste0(regdir, 'aggregated_allres.', region, '.rda')
    load(full.file)

    # Turn DEGs into a single data frame:
    # -----------------------------------
    dedf = c()
    for (set in names(setdflist)){
        setdf = setdflist[[set]]
        if (!is.null(setdf)){
            setdf$set = set
            dedf = rbind(dedf, setdf)
        }
    }

    gwdedf = dedf[(dedf$gene %in% gwgenes),]
    gwdedf$abseff = abs(gwdedf$logFC_nb)
    gwdedf$region = region
    mdf = aggregate(log10p_nm ~ gene + set, gwdedf, max)
    topdedf = merge(gwdedf, mdf)
    # Break top-level p-value ties with effect (higher always for nrad/cogdxad)
    mdf = aggregate(abseff ~ gene + set, topdedf, max)
    topdedf = merge(topdedf, mdf)

    gw.reg.file = paste0(regdir, 'aggregated_allres.gwgenes.', region, '.rda')
    save(gwdedf, topdedf, file=gw.reg.file)
}


