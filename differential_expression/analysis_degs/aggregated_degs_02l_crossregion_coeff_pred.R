#!/usr/bin/R
# --------------------------------------------------
# Starting from all DEG results,
# - Predict each region from others + from pathology
# Updated: 11/02/23
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(patchwork)
print(version)
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


# Functions:
# ----------
load_psbulk_matrix = function(set){
    ststr = sub('_[A-Za-z]+$', '', set)
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
    load(psdata.rda)
    # If Mic, remove T, CAMs from psbulk data:
    if (set == 'Mic_Immune_Mic'){
        ind = !(ps.data$meta$cell_type_high_resolution %in% c('T cells', 'CAMs'))
        ps.data$meta = ps.data$meta[ind,]
        ps.data$mat = ps.data$mat[,rownames(ps.data$meta)]
    }
    # Merge to sample-level (individual x region):
    ps.data = aggregate_psbulk_samplelevel(ps.data)
    return(ps.data)
}

pivot.tomatrix.v2 = function(df, index, key, value){
    require(dplyr)
    df = df[,c(index, key, value)]
    wide = spread(df, key, value)
    mat = as.matrix(wide[, -1])
    rownames(mat) = wide[, 1]
    return(mat)
}

# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
degcols = c('NS'='grey90','Down'=col.paired[2],'Up'=col.paired[6])
denumcols = c('0'='grey90','1'=col.paired[2],'2'=col.paired[6])


# Load in DEGs across all regions, path vars, sets
# ------------------------------------------------
alldf = c()
kept.cols = c('gene','col_nm','path','region', 'logFC_nb', 'p_nb', 'coef_mast', 'p_mast')
for (path in pathlist){
    cat(path, '\n')
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)

    for (set in keep.sets){
        cat('-', set, '\n')
        setdf = setdflist[[set]][, kept.cols]
        setdf$set = set
        alldf = rbind(alldf, setdf)
    }
}


# Run analysis for each set separately:
# -------------------------------------
# lvar = 'logFC_nb'
# lvar = 'coef_mast'
lvar = 'avg_mn'
cvar = 'col_nm'
rescale_mast = sd(alldf$logFC_nb) / sd(alldf$coef_mast)
alldf$avg_mn = (alldf$logFC_nb + alldf$coef_mast * rescale_mast) / 2

path = 'nft'
coeffdf = c()
cs.fulldf = c()
for (set in keep.sets){
    suff = paste0(set, '.', path, '.', lvar)
    subdf = alldf[(alldf$set == set) & (alldf$path == path),]
    df = spread(subdf[,c('region','gene', lvar)], 'region', lvar)
    cdf = spread(subdf[,c('region','gene', cvar)], 'region', cvar)
    names(cdf)[2:ncol(cdf)] = paste0('col_', names(cdf)[2:ncol(cdf)])
    df = merge(df, cdf, all=TRUE)

    regset = reg.nomb[reg.nomb %in% subdf$region]
    regorder = c('EC', 'HC', 'TH', 'AG', 'MT', 'PFC')
    regorder = regorder[regorder %in% regset]
    resdf = c()
    seldf = c()
    updf = c()
    ntop = 20 # GGplot scatters
    nsig = 10 # Heatmaps
    for (region in regset){
        cat(region, '\n')
        nonregion = regset[regset != region]
        formvec = c(region, '~', paste(nonregion, collapse=' + '))
        form = asform(formvec)
        print(form)

        fit = glm(form, data=df, family='gaussian')
        cfit = as.data.frame(coefficients(summary(fit)))
        names(cfit) = c('Est', 'SE', 't', 'p')
        cfit$var = rownames(cfit)
        cfit$region = region
        cfit$path = path
        resdf = rbind(resdf, cfit)
        # TODO: Get summary of fit quality

        # Get consistent / inconsistent genes:
        # ------------------------------------
        # TODO: Prediction with missingness?
        # TODO: get SE from orig analysis?
        predvar = paste0('pred_', region)
        crvar = paste0('col_', region)
        df[[crvar]] = factor(df[[crvar]])
        df[[predvar]] = predict(fit, df)
        rdf = df[!is.na(df[[predvar]]),]
        rdf = rdf[!is.na(rdf[[region]]),]
        rdf$resid = rdf[[region]] - rdf[[predvar]]

        sig.rdf = rdf[rdf[[crvar]] != 0,]
        sig.rdf = sig.rdf[order(sig.rdf[[predvar]]),]
        top.pred = c(head(sig.rdf$gene, ntop),
            tail(sig.rdf$gene, ntop))

        sig.rdf = sig.rdf[order(sig.rdf$resid),]
        top.resid = c(head(sig.rdf$gene, ntop),
            tail(sig.rdf$gene, ntop))
        top.genes = unique(c(top.pred, top.resid))
        labdf = rdf[rdf$gene %in% top.genes,]

        # Select specific genes for later:
        sig.rdf = sig.rdf[order(abs(sig.rdf$resid), decreasing=T),]
        ind = ifelse(sig.rdf[[crvar]] == 1, sig.rdf$resid < 0, sig.rdf$resid > 0)
        selgenes = head(sig.rdf[ind,'gene'], nsig)
        selrdf = data.frame(gene=selgenes, region = region)
        seldf = rbind(seldf, selrdf)
        # Up only:
        upgenes = head(sig.rdf[ind & (sig.rdf[[region]] > 0),'gene'], nsig)
        uprdf = data.frame(gene=upgenes, region = region)
        updf = rbind(updf, uprdf)

        # TODO: Consistent vs. different
        # TODO: Make sure no flips
        g1 = ggplot(rdf, aes_string(region, predvar, col=crvar)) + 
            geom_point(cex=.25) + 
            # geom_smooth(method='lm', color='black') +
            # stat_cor(color='black') + 
            geom_hline(yintercept=0, lty='dashed') + 
            geom_vline(xintercept=0, lty='dashed') + 
            scale_color_manual(values=denumcols) + 
            geom_text_repel(data=labdf, aes_string(region, predvar, label='gene', col=crvar), cex=2,
                box.padding=0.01, max.overlaps=20, force=20, force_pull=0.5, min.segment.length=0.02) + 
            geom_abline(intercept=0, slope=1, lty='dotted') + 
            labs(x=paste0(lvar,' in ', region), y='Prediction from other regions') +
            theme_pubr() + theme(legend.position='none')

        # TODO: Consistent vs. different
        # TODO: Make sure no flips
        g2 = ggplot(rdf, aes_string(region, 'resid', col=crvar)) + 
            geom_point(cex=.25) + 
            # geom_smooth(method='lm', color='black') +
            # stat_cor(color='black') + 
            geom_hline(yintercept=0, lty='dashed') + 
            geom_vline(xintercept=0, lty='dashed') + 
            scale_color_manual(values=denumcols) + 
            geom_text_repel(data=labdf, aes_string(region, 'resid', label='gene', col=crvar), cex=2,
                box.padding=0.01, max.overlaps=20, force=20, force_pull=0.5, min.segment.length=0.02) + 
            geom_abline(intercept=0, slope=1, lty='dotted') + 
            labs(x=paste0(lvar,' in ', region), y='Residual from prediction') +
            theme_pubr() + theme(legend.position='none')

        garr = g1 + g2
        pltprefix = paste0(imgpref, 'scatter_pred.', region, '.', suff)
        saveGGplot(garr, pltprefix, w=8, h=4)

        # TODO: TABLE of consistent vs. inconsistent
    }
    coeffdf = rbind(coeffdf, resdf)


    # Plot coefficients:
    # ------------------
    cmat = pivot.tomatrix.v2(resdf, 'var', 'region', 'Est')
    pmat = pivot.tomatrix.v2(resdf, 'var', 'region', 'p')

    mx = 1
    col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
    ux = 1.5
    ht = Heatmap(cmat, 
        use_raster=FALSE, 
        name=suff,
        col=col_fun,
        cluster_columns=FALSE, 
        cluster_rows=FALSE,
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lwd=.5),
        width=ncol(cmat) * unit(ux, 'mm'),
        height=nrow(cmat) * unit(ux, 'mm'),
        cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
            p = pmat[i,j]
            if (!is.na(p) & (p < 0.01)){
                grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.1))
            }
        })

    pltprefix = paste0(imgpref, 'crossregion_pred_est_heatmap.', suff)
    h = 1 + 1 / 15 * nrow(cmat)
    w = 2 + 1 / 15 * ncol(cmat)
    saveHeatmap(ht, pltprefix, w=w, h=h)


    # Plot consistent genes:
    # ----------------------
    predcols = paste0("pred_", regset)
    sigcols = paste0("col_", regset)
    resid = df[, regset] - df[, predcols]
    # NOTE: ALT
    # df$abs_mean_lvar = apply(abs(df[,regset]), 1, mean)
    # df$mean_resid = apply(abs(resid), 1, mean)
    df$mean_lvar = apply(df[,regset], 1, mean)
    df$abs_mean_lvar = abs(df$mean_lvar)
    df$mean_resid = apply(abs(resid), 1, mean)
    df$nsig = apply(df[,sigcols] != 0, 1, sum)
    sdf = df[!is.na(df$mean_resid),]

    # Order:
    fit = glm(mean_resid ~ abs_mean_lvar, data=sdf, family='gaussian')
    sdf$pred_mean_resid = predict(fit, sdf)
    sdf$diff_resid_pred = sdf$mean_resid - sdf$pred_mean_resid
    sdf = sdf[order(sdf$diff_resid_pred),]
    ntop = 20
    top.genes = c(head(sdf$gene[sdf$nsig > 0], ntop),
        tail(sdf$gene[sdf$nsig > 0], ntop))
    labdf = sdf[sdf$gene %in% top.genes,]
    consistent.genes = head(sdf$gene[sdf$nsig > 0], nsig)
    up.consistent.genes = head(sdf$gene[(sdf$nsig > 0) & (sdf$mean_lvar > 0)], nsig)

    gp = ggplot(sdf, aes(abs_mean_lvar, mean_resid, col=nsig)) + 
        geom_point(cex=.25) + 
        geom_smooth(method='lm') + 
        geom_hline(yintercept=0, lty='dashed') + 
        geom_vline(xintercept=0, lty='dashed') + 
        geom_text_repel(data=labdf, aes(abs_mean_lvar, mean_resid, label=gene), cex=2,
            box.padding=0.01, max.overlaps=20, force=20, force_pull=0.5, min.segment.length=0.02) + 
        labs(x=paste0('|Mean ',lvar,'|'), y='Mean |residual|') +
        theme_pubr() + theme(legend.position='none')
    pltprefix = paste0(imgpref, 'scatter_meansd_pred.', suff)
    saveGGplot(gp, pltprefix, w=4, h=4)

    gp = ggplot(sdf, aes(abs_mean_lvar, mean_resid - pred_mean_resid, col=nsig)) + 
        geom_point(cex=.25) + 
        geom_smooth(method='lm') + 
        geom_hline(yintercept=0, lty='dashed') + 
        geom_vline(xintercept=0, lty='dashed') + 
        geom_text_repel(data=labdf, aes(abs_mean_lvar, mean_resid - pred_mean_resid, label=gene), cex=2,
            box.padding=0.01, max.overlaps=20, force=20, force_pull=0.5, min.segment.length=0.02) + 
        labs(x=paste0('|Mean ',lvar,'|'), y='Mean |residual| - predicted') +
        theme_pubr() + theme(legend.position='none')
    pltprefix = paste0(imgpref, 'scatter_meansd_pred_diff.', suff)
    saveGGplot(gp, pltprefix, w=4, h=4)


    # Plot consistent and specific genes:
    # -----------------------------------
    for (plot.uponly in c(TRUE, FALSE)){
        if (plot.uponly){
            con.seldf = rbind(updf, data.frame(gene=up.consistent.genes, region='Consistent'))
            con.seldf$set = set
            pltprefix = paste0(imgpref, 'consistentspecific_topgenes_heatmap.', suff, '.uponly')
        } else {
            con.seldf = rbind(seldf, data.frame(gene=consistent.genes, region='Consistent'))
            con.seldf$set = set
            pltprefix = paste0(imgpref, 'consistentspecific_topgenes_heatmap.', suff)
        }

        # Make matrices for plotting:
        sdf = merge(df, con.seldf)
        rownames(sdf) = paste0(sdf$gene, ':', sdf$region)
        cs.fulldf = rbind(cs.fulldf, con.seldf)
        cmat = as.matrix(sdf[, regorder])
        pmat = as.matrix(sdf[, paste0('col_', regorder)])

        # Order in same way:
        rowreg = sub(".*:","", rownames(cmat))
        rowreg = factor(rowreg, levels=c('Consistent', regorder))
        roword = order(rowreg)
        cmat = cmat[roword,]
        pmat = pmat[roword,]
        rmap = paste0(1:length(levels(rowreg)), '-', levels(rowreg))
        names(rmap) = levels(rowreg)
        rowsplit = rmap[sub(".*:","", rownames(cmat))]

        mx = .5
        if (path %in% c('nft','plaq_n', 'plaq_d')){ mx = mx / 20 }
        col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
        ux = 1.5
        ht = Heatmap(cmat, 
            use_raster=FALSE, 
            name=lvar,
            col=col_fun,
            column_title=suff,
            cluster_columns=FALSE, 
            row_split=rowsplit, 
            cluster_rows=TRUE,
            cluster_row_slices=FALSE,
            row_dend_width = unit(.25, "cm"),
            column_dend_height = unit(.25, "cm"),
            row_dend_gp = gpar(lwd=.5),
            column_dend_gp = gpar(lwd=.5),
            border_gp = gpar(col="black", lwd=.5),
            width=ncol(cmat) * unit(ux, 'mm'),
            height=nrow(cmat) * unit(ux, 'mm'),
            cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                p = pmat[i,j]
                if (p != 0){
                    grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.25))
                }
            })

        h = .5 + 1 / 15 * nrow(cmat)
        w = 1.5 + 1 / 15 * ncol(cmat)
        saveHeatmap(ht, pltprefix, w=w, h=h)
    }

    # TODO: Save consistent / selected genes:

}


