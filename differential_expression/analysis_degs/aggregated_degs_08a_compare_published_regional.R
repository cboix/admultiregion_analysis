#!/usr/bin/R
# --------------------------------------
# Comparison of published with own DEGs:
# - regional DEGs in this analysis
# Updated: 07/27/23
# --------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'comppubAD_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Load published and own DEGs:
# ----------------------------
# Load published DEGs:
source(paste0(sbindir, 'differential_expression/analysis_degs/load_published_AD_DEGs.R'))

# Count of up/down/ns on these:
ndf = merge(agg.rename(gene ~ sig + run + region, all.degs.df, length, 'ngene'),
    agg.rename(gene ~ run + region + only.sig, all.degs.df, length, 'ntot'))

# Regional: NIA-Reagan
path = 'nrad'
mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)
sets = names(setdflist)
reglist = unique(totnsigdf$region)

rmsets = c('Mic_Immune_T', 'Mic_Immune_CAM',
    'Vasc_Epithelia_CPEC', 'Vasc_Epithelia_Epd', 'Vasc_Epithelia_Per',
    'Vasc_Epithelia_SMC', 'Vasc_Epithelia_Fib')
sets = sets[!(sets %in% rmsets)]


# Plot overlap to explore consistency:
# ------------------------------------
all.cor.df = c()
for (set in sets){
    for (region in reglist){
        cat(set, '\t', region, '\n')
        setdf = setdflist[[set]]
        setdf = setdf[setdf$region == region, ]
        if (nrow(setdf) > 25){
            setdf$set = set
            lfc.map = setdf$logFC_nb
            names(lfc.map) = setdf$gene
            subdf = all.degs.df[all.degs.df$gene %in% setdf$gene,]
            subdf$logFC_nb = lfc.map[subdf$gene]
            runs = unique(subdf$run)
            out = sapply(runs, 
                function(x){
                    rundf = subdf[subdf$run == x,]
                    if (nrow(rundf) > 25){
                        ct = cor.test(rundf$lfc, rundf$logFC_nb, method='spearman')
                        return(data.frame(est=ct$estimate, p=ct$p.value, n=nrow(rundf)))
                    }
                }
            )
            cordf = data.frame(do.call(rbind, out))
            cordf$run = rownames(cordf)
            cordf = cordf[order(cordf$p),]
            cordf$set = set
            cordf$region = region
            all.cor.df = rbind(all.cor.df, cordf)
        }
    }
}

all.cor.df = merge(all.cor.df, unique(all.degs.df[,c('run', 'celltype', 'study', 'path')]))
all.cor.df = merge(all.cor.df, pub.cellmap, all.x=TRUE)
all.cor.df$sp = with(all.cor.df, paste0(set, '@', region))
all.cor.df = all.cor.df[all.cor.df$set %in% sets,]
# all.cor.df = all.cor.df[all.cor.df$study != 'Leng_2021',]
all.cor.df$padj = p.adjust(all.cor.df$p, method='BH')


# Plot results from correlation analysis:
# ---------------------------------------
all.cor.df$run2 = with(all.cor.df, paste0(study, '-', celltype, ':', path, '@', major.celltype))
cmat = pivot.tomatrix(all.cor.df[,c('run2', 'sp', 'est')], 'sp','est')
pmat = pivot.tomatrix(all.cor.df[,c('run2', 'sp', 'padj')], 'sp','padj')

dfind = cbind(apply(pmat, 2, which.min), 1:ncol(pmat))
pmat2 = pmat * 0 + 1
pmat2[dfind] = 0

roword = order(dfind[,1])
mx = 1
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))

ux = 1.5
pltmat = cmat
pltmat = reord(pltmat)
pltmat = diag.mat2(pltmat)[[1]]
csplit = gsub("_", "\n", sub("@.*","", colnames(pltmat)))
rsplit = sub(".*@","", rownames(pltmat))
cmat2 = cmat[rownames(pltmat), colnames(pltmat)]
pmat2 = pmat2[rownames(pltmat), colnames(pltmat)]
pmat3 = pmat[rownames(pltmat), colnames(pltmat)]
# pltmat = sweep(pltmat, 2, apply(pltmat, 2, max, na.rm=T), '/')
pltmat = sweep(pltmat, 1, apply(pltmat, 1, max, na.rm=T), '/')
ht = Heatmap(pltmat, 
    use_raster=TRUE, 
    col=col_fun,
    column_split=csplit,
    cluster_column_slices=FALSE,
    cluster_row_slices=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    row_split=rsplit,
    border_gp=gpar(color='black', lwd=.5),
    width=ncol(cmat) * unit(ux, 'mm'),
    height=nrow(cmat) * unit(ux, 'mm'), 
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        p = pmat3[i,j]
        cc = cmat2[i,j]
        # if (!is.na(p) & (p < 0.001) & (cc > .3)){
        if (!is.na(p) & (p < 0.001) & (cc > .2)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs * 1.5))
        }
    })

pltprefix = paste0(imgpref, 'byregion_', path, 'cortest_heatmap')
h = 1 + 1 / 15 * nrow(cmat)
w = 2 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)


pltmat = -log10(pmat)
pltmat[pltmat > 50] = 50
pltmat = reord(pltmat)
pltmat = diag.mat2(pltmat)[[1]]
csplit = sub("@.*","", colnames(pltmat))
rsplit = sub(".*@","", rownames(pltmat))
cmat2 = cmat[rownames(pltmat), colnames(pltmat)]
pmat2 = pmat2[rownames(pltmat), colnames(pltmat)]
pmat3 = pmat[rownames(pltmat), colnames(pltmat)]
# pltmat = sweep(pltmat, 2, apply(pltmat, 2, max, na.rm=T), '/')
pltmat = sweep(pltmat, 1, apply(pltmat, 1, max, na.rm=T), '/')
ht = Heatmap(pltmat, 
    use_raster=FALSE, 
    col=viridis(100),
    column_split=gsub("_", "\n", csplit),
    cluster_column_slices=FALSE,
    cluster_row_slices=FALSE,
    row_split=rsplit,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    border_gp=gpar(color='black', lwd=.5),
    width=ncol(pltmat) * unit(ux, 'mm'),
    height=nrow(pltmat) * unit(ux, 'mm'), 
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        plab = pltmat[i,j]
        p = pmat3[i,j]
        cc = cmat2[i,j]
        if (!is.na(p) & (p < 0.001) & (cc > .2)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs * 1.5, col=ifelse(plab > .6, 'black', 'white')))
        }
    })

pltprefix = paste0(imgpref, 'byregion_', path, 'cortest_pval_heatmap')
h = 1.5 + 1 / 15 * nrow(pltmat)
w = 1.5 + 1 / 15 * ncol(pltmat)
saveHeatmap(ht, pltprefix, w=w, h=h)






# NO GUARANTEE OF OVERLAP!
# Assumes same background - wrong.

        # uplist = setdf$gene[setdf$col_nm == 2]
        # dwlist = setdf$gene[setdf$col_nm == 1]
        # all.degs.df$is.up = 1*(all.degs.df$gene %in% uplist)
        # all.degs.df$is.dw = 1*(all.degs.df$gene %in% dwlist)
        # aggdf = aggregate(cbind(is.up, is.dw) ~ run + study + celltype + region + sig + only.sig, all.degs.df, sum)
        # aggdf = merge(aggdf, ndf)  # Add aggregate stats
        # aggdf$n.up = length(uplist)
        # aggdf$n.dw = length(dwlist)
        # aggdf$ntest = nrow(setdf)
        # aggdf$stat.up = with(aggdf, log((is.up/n.up)/(ngene/ntot)))
        # aggdf$stat.dw = with(aggdf, log((is.dw/n.dw)/(ngene/ntot)))
        # head(aggdf[order(aggdf$stat.up, decreasing=T),], 20)


dw.hgdf = aggdf[,c('is.dw','ngene','ntot', 'n.dw')] # So wrong! FIX
up.hgdf = aggdf[,c('is.up','ngene','ntot', 'n.up')]
dw.hgdf = aggdf[,c('is.dw','ngene','ntot', 'n.dw')]
aggdf$p.up = apply(up.hgdf, 1, run.hyper)
aggdf$p.dw = apply(dw.hgdf, 1, run.hyper)
# up.hgdf[order(up.hgdf$p),]
head(aggdf[order(aggdf$p.up),], 20)

# degdf = rbind(degdf, setdf)
degs[[paste0(set, "_up")]] = setdf$gene[setdf$col_nm == 2]
degs[[paste0(set, "_down")]] = setdf$gene[setdf$col_nm == 1]


# TODO: Also evaluate discrete (enr) analysis:

# TODO: Also evaluate 
