#!/usr/bin/R
# -------------------------------------------------
# Are region scores different
# Updated: 06/23/23
# -------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(ggrastr)
print(version)
options(width=170)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Load modules:
# -------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Load regional genes for each pathology measure:
# -----------------------------------------------
pathlist = c('plaq_d', 'plaq_n', 'nft')
keep.sets = c('Mic_Immune_Mic', 'Ast_Ast', 'Oli_Oli', 'Opc_Opc')
kept.cols = c('gene','col_nm','path','region')
alldf = c()
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)
    for (set in keep.sets){
        print(set)
        setdf = setdflist[[set]][, kept.cols]
        setdf = setdf[setdf$col_nm != 0,]
        setdf$set = set
        setdf$path = path
        alldf = rbind(alldf, setdf)
    }
}


# Annotate with relevant modules:
# -------------------------------
modmap = list('Oli' = c(7), 'Opc' = c(7,9,24), 'Ast'=c(6, 27), 'Mic_Immune'=c(11, 25))

resdf = c()
for (set in keep.sets){
    use.core = TRUE
    setstr = sub("_[A-Za-z]+$","", set)
    cat(set, setstr, '\n')
    coremap = cmlist[[setstr]]
    genemap = gmlist[[setstr]]
    if (use.core){ usemap = coremap } else { usemap = genemap }

    for (path in pathlist){
        subdf = alldf[(alldf$set == set) & (alldf$path == path),]
        subdf$module = usemap[subdf$gene]
        subdf = subdf[!is.na(subdf$module),]
        kept.genes = unique(subdf$gene)

        # Count # genes per module:
        ndf = data.frame(table(usemap))
        names(ndf) = c('module', 'nmod')
        mod.cutoff = 10 # Modules with at least 10 genes
        kept.modules = ndf$module[ndf$nmod >= mod.cutoff]

        # Enrichment:
        hgdf = agg.rename(gene ~ module + col_nm + region, subdf, length, 'nint')
        hgdf = merge(hgdf, agg.rename(gene ~ col_nm + region, subdf, length, 'nsig'))
        hgdf = merge(hgdf, ndf)
        hgdf = hgdf[hgdf$nmod >= mod.cutoff,]
        hgdf$ntot = sum(usemap %in% kept.modules)

        # Run hypergeometric tests:
        # -------------------------
        pvdf = hgdf[,c('nint','nmod', 'nsig','ntot')]
        pout <- apply(pvdf, 1, run.hyper)
        hgdf$p = pout
        hgdf$lp = -log10(pout)
        hgdf = hgdf[order(hgdf$p),]
        hgdf$padj = p.adjust(hgdf$p,'BH') # Correct by BH
        hgdf$log2FC = with(hgdf, log2((nint / nsig) / (nmod / ntot)))
        hgdf$set = setstr
        hgdf$path = path
        resdf = rbind(resdf, hgdf[hgdf$module %in% modmap[[setstr]],])
    }
}
resdf = merge(resdf, aggregate(p ~ region + module + set + path, resdf, min))
resdf$region = sub('allregions', 'All', resdf$region)
resdf$pr = with(resdf, paste0(region, ':', path))
resdf$tag = with(resdf, paste0(set, '-', module))
resdf$signed.log2FC = resdf$log2FC * (resdf$col_nm - 1.5) * 2


# Plot these modules:
# -------------------
cmat = pivot.tomatrix(resdf[,c('tag','pr','signed.log2FC')], 'pr', 'signed.log2FC')
pmat = pivot.tomatrix(resdf[,c('tag','pr','padj')], 'pr', 'padj')

cmat[is.na(cmat)] = 0
pmat[is.na(pmat)] = 1
cnord = expand.grid(region= reg.order[-length(reg.order)], path=pathlist)
cn = with(cnord, paste0(region, ":", path))
cmat = cmat[,cn]
pmat = pmat[,cn]

mx = 3
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))

ux = 1.5
plt = Heatmap(cmat,
    use_raster=FALSE,
    name='signed.log2FC',
    col=col_fun,
    cluster_columns=FALSE,
    cluster_column_slices=FALSE,
    cluster_rows=TRUE,
    column_split=sub(".*:","", colnames(cmat)),
    width = ncol(cmat)*unit(ux, "mm"), 
    height = nrow(cmat)*unit(ux, "mm"),
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        p = pmat[i,j]
        if (!is.na(p)){
            ann = ifelse(p < 0.1, ifelse(p < 0.05, '*', '.'),'')
            grid.text(ann, x, y,gp=gpar(fontsize=gridtxt.fs))
        }
    })

pltprefix = paste0(imgpref, 'opcoli_enr_heatmap')
h = 1 + 1 / 15 * nrow(cmat)
w = 1.5 + 1 / 15 * ncol(cmat) 
saveHeatmap(plt, pltprefix, w=w, h=h)

