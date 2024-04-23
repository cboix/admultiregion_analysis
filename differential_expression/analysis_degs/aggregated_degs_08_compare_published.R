#!/usr/bin/R
# -----------------------------------------
# Comparison of published with own DEGs:
# - major cell types, across pathology vars
# Updated: 07/27/23
# -----------------------------------------
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
ndf = merge(agg.rename(gene ~ sig + run + path, all.degs.df, length, 'ngene'),
    agg.rename(gene ~ run + path + only.sig, all.degs.df, length, 'ntot'))

# Load own, overall DEGs:
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
load(fullaggrda)
sets = names(setdflist)
pathlist = unique(totnsigdf$path)

rmsets = c('Mic_Immune_T', 'Mic_Immune_CAM',
    'Vasc_Epithelia_CPEC', 'Vasc_Epithelia_Epd', 'Vasc_Epithelia_Per',
    'Vasc_Epithelia_SMC', 'Vasc_Epithelia_Fib')
sets = sets[!(sets %in% rmsets)]


# Plot overlap to explore consistency:
# ------------------------------------
all.cor.df = c()
all.enr.df = c()
for (set in sets){
    for (path in pathlist){
        cat(set, '\t', path, '\n')
        setdf = setdflist[[set]]
        setdf = setdf[setdf$path == path, ]
        setdf$set = set

        lfc.map = setdf$logFC_nb
        sig.map = setdf$col_nm
        names(lfc.map) = setdf$gene
        names(sig.map) = setdf$gene
        subdf = all.degs.df[all.degs.df$gene %in% setdf$gene,]
        subdf$logFC_nb = lfc.map[subdf$gene]
        subdf$sig.nm = sig.map[subdf$gene]

        # Perform correlation test for directionality of shared reported genes
        runs = unique(subdf$run)
        out = sapply(runs, 
            function(x){
                rundf = subdf[subdf$run == x,]
                if (nrow(rundf) > 25){
                    ct = cor.test(rundf$lfc, rundf$logFC_nb, method='pearson')
                    return(data.frame(est=ct$estimate, p=ct$p.value, n=nrow(rundf)))
                }
            }
        )
        cordf = data.frame(do.call(rbind, out))
        cordf$run = rownames(cordf)
        cordf = cordf[order(cordf$p),]
        cordf$set = set
        cordf$path = path
        all.cor.df = rbind(all.cor.df, cordf)

        # Also perform an enrichment test - overlap of down/down and up/up within overlapped DEGs
        sig.subdf = subdf[subdf$sig != 'NS',] # Remove NS so all comparable
        aggdf = merge(agg.rename(gene ~ run + sig + sig.nm, sig.subdf, length, 'nint'),
            agg.rename(gene ~ run + sig.nm, sig.subdf, length, 'nown'))
        aggdf = merge(aggdf, agg.rename(gene ~ run + sig, subdf, length, 'npub'))
        aggdf = merge(aggdf, agg.rename(gene ~ run, sig.subdf, length, 'ntot'))
        aggdf = rbind(aggdf[(aggdf$sig == 'Down') & (aggdf$sig.nm == 1),],
            aggdf[(aggdf$sig == 'Up') & (aggdf$sig.nm == 2),])
        if (nrow(aggdf) > 0){
            hgdf = aggdf[,c('nint','nown','npub', 'ntot')]
            aggdf$p = apply(hgdf, 1, run.hyper)
            aggdf$padj = p.adjust(aggdf$p, method='BH')
            aggdf = aggdf[order(aggdf$p),]
            aggdf$lfc = with(aggdf, log((nint/nown)/(npub/ntot)))
            aggdf$set = set
            aggdf$path = path
            all.enr.df = rbind(all.enr.df, aggdf)
        }
    }
}

all.enr.df = merge(all.enr.df, unique(all.degs.df[,c('run', 'celltype', 'study')]))
all.enr.df = merge(all.enr.df, pub.cellmap, all.x=TRUE)
all.enr.df$sp = with(all.enr.df, paste0(set, '@', path))
all.enr.df = all.enr.df[all.enr.df$study != 'Leng_2021',]
all.enr.df$padj = p.adjust(all.enr.df$p, method='BH')

all.cor.df = merge(all.cor.df, unique(all.degs.df[,c('run', 'celltype', 'study')]))
all.cor.df = merge(all.cor.df, pub.cellmap, all.x=TRUE)
all.cor.df$sp = with(all.cor.df, paste0(set, '@', path))
all.cor.df = all.cor.df[!(all.cor.df$study %in% c('Leng_2021')),]
all.cor.df$padj = p.adjust(all.cor.df$p, method='BH')


# Plot results from correlation analysis:
# ---------------------------------------
all.cor.df$run2 = with(all.cor.df, paste0(study, '-', celltype, '@', major.celltype))
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
        if (!is.na(p) & (p < 0.001) & (cc > .3)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs * 1.5))
        }
    })

pltprefix = paste0(imgpref, 'allregions_cortest_heatmap')
h = 1 + 1 / 15 * nrow(cmat)
w = 2 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)


pltmat = -log10(pmat)
pltmat[pltmat > 200] = 200
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
    # cluster_rows=FALSE,
    # show_row_names=FALSE)
    border_gp=gpar(color='black', lwd=.5),
    width=ncol(pltmat) * unit(ux, 'mm'),
    height=nrow(pltmat) * unit(ux, 'mm'), 
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    # column_title=paste0('Modules (', setstr, ')'),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        plab = pltmat[i,j]
        p = pmat3[i,j]
        cc = cmat2[i,j]
        # p = pmat[i,j]
        if (!is.na(p) & (p < 0.001) & (cc > .3)){
            # ann = ifelse(p < 0.1, ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***', '**'), '*'), '.'),'')
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs * 1.5, col=ifelse(plab > .6, 'black', 'white')))
        }
    })

pltprefix = paste0(imgpref, 'allregions_cortest_pval_heatmap')
h = 1.5 + 1 / 15 * nrow(pltmat)
w = 1.5 + 1 / 15 * ncol(pltmat)
saveHeatmap(ht, pltprefix, w=w, h=h)






# Plot results from enrichment analysis:
# ---------------------------------------
# Aggregate up + down pvalues
all.enr.df$run2 = with(all.enr.df, paste0(study, '-', celltype, '@', major.celltype))
l.enr.df = aggregate(lfc ~ sp + run2, all.enr.df, mean)
# NOTE: not completely ok, because missingness of some combinations 
# TODO: Add all down/up with all sp + run2
p.enr.df = aggregate(padj ~ sp + run2, all.enr.df, function(x){ 
    if (length(x) > 1){
    pstat = -2 * sum(log(x))
    p = pchisq(pstat, df=length(x), lower.tail=FALSE)
    } else { p = x}
    return(p)
    })
cmat = pivot.tomatrix(l.enr.df, 'sp','lfc')
pmat = pivot.tomatrix(p.enr.df, 'sp','padj')
pmat[is.na(pmat)] = 1
cmat[is.na(cmat)] = 0

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
        if (!is.na(p) & (p < 0.001) & (cc > .3)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs * 1.5))
        }
    })

pltprefix = paste0(imgpref, 'allregions_enrtest_heatmap')
h = 1 + 1 / 15 * nrow(cmat)
w = 2 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)


