#!/usr/bin/R
# --------------------------------------------------------
# Plot the aggregated DE results across different methods:
# Number of DEGs for msex * AD interaction model
# Updated: 11/15/23
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)

library(ggplot2)
library(ggrepel)
library(ggrastr)
library(ggpubr)
options(width=170) 
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Get a list of all differential runs:
# -------------------------------------------------------
rundf = read.delim(paste0(sdbdir, 'DEG_multiRegion_SI_ACE_runlist.tsv'), header=T)
rundf$prefstr = with(rundf, gsub("[()/]", "_", paste(celltype, subtype, region, path, sep="_")))
rundf = rundf[rundf$path == 'msex',]


# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
    outfile = paste0(regdir, 'allmethods.', x, '.merged.rda')
    1 * file.exists(outfile) })
table(rundf$merged)
head(rundf[rundf$merged == 0,])


# Select runs to use (allregions, main ct):
# -----------------------------------------
path = 'msex'
print(path)
selrundf = rundf[rundf$path == path,]
selrundf$setid = with(selrundf, paste0(celltype, '_', subtype))
sets = unique(selrundf$setid)
sets = c(sets[grep("^Exc_",sets, invert=TRUE)], 'Exc_Exc')
selrundf = selrundf[selrundf$setid %in% sets,]
print(sets)


# Aggregate the DEGs and results across all regional runs for NFT:
# ----------------------------------------------------------------
mstr = paste0('allmethods.', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
if (!file.exists(fullaggrda)){
    kept.cols = c("gene","pc", "col_nm","log10p_nm",
        "path","region",
        "logFC_nb","p_nb","padj_nb","col_nb",
        "coef_mast","p_mast","padj_mast","col_mast")
    setdflist = lapply(sets, function(x){})
    names(setdflist) = sets
    totnsigdf = NULL
    for (i in 1:nrow(selrundf)){
        prefstr = selrundf$prefstr[i]
        setid = selrundf$setid[i]
        aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
        if (file.exists(aggrda)){
            load(aggrda)  # Loads aggdf, nsig
            print(nsig)
            # cat(nsig, '\n')
            # Concatenate results:
            aggdf$path = nsig[1,'path']
            aggdf$region = nsig[1,'region']
            setdflist[[setid]] = rbind(setdflist[[setid]], aggdf)
            # Pad nsig:
            nsigdf = data.frame(nsig)
            if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
            if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
            totnsigdf = rbind(totnsigdf, nsigdf)
        }
    }
    save(setdflist, totnsigdf, file=fullaggrda)
} else {
    load(fullaggrda)
}


# Plot total number of significant genes as a heatmap:
# ----------------------------------------------------
pathstr = unique(totnsigdf$eff)
totnsigdf$X1 = as.numeric(totnsigdf$X1)
totnsigdf$X2 = as.numeric(totnsigdf$X2)
for (pstr in pathstr){
    umat = pivot.tomatrix(totnsigdf[totnsigdf$eff == pstr,c('region','subtype','X2')], 'region','X2')
    dmat = pivot.tomatrix(totnsigdf[totnsigdf$eff == pstr,c('region','subtype','X1')], 'region','X1')
    umat[is.na(umat)] = 0
    dmat[is.na(dmat)] = 0
    regset = c('allregions','EC','HC','TH','AG','MT','PFC')
    umat = umat[,regset]
    dmat = dmat[,regset]
    colnames(umat) = paste0(colnames(umat), '_up')
    colnames(dmat) = paste0(colnames(dmat), '_down')

    bmat = cbind(umat, -dmat)
    mx = max(abs(bmat))

    # Column, row splits:
    csplit = ifelse(1:ncol(bmat) %in% grep("_up", colnames(bmat)), 'Up','Down')
    keep.rows = c('Ast','Mic','Oli','Opc','Inh','Exc')
    rn = keep.rows[keep.rows %in% rownames(bmat)]
    bmat = bmat[rn,]

    col_fun = colorRamp2(c(-300, 0, 300), c("blue", "white", "red"))
    plt = Heatmap(bmat, name='Number\nof DEGs',
        use_raster=TRUE,
        col=col_fun,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        column_split=csplit,
        border_gp=gpar(color='black', lwd=.5),
        width = ncol(bmat)*unit(4.5, "mm"), 
        height = nrow(bmat)*unit(2, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
            ann = abs(bmat[i,j])
            grid.text(ann, x, y, gp=gpar(fontsize=5)) })

    h = 2.25 + 1 / 15 * nrow(bmat)
    w = 5 + 1 / 15 * ncol(bmat)
    pltprefix = paste0(imgpref, 'allmethods_ndeg.regional_', path, '.', pstr, '.heatmap')
    saveHeatmap(plt, pltprefix, w=w, h=h)
}



# Plot side-by-side volcano plots for each of the main cell types:
# ----------------------------------------------------------------
pcols = brewer.pal(12,'Paired')
decols = c('0'='grey80','1'=pcols[1],'2'=pcols[5])
totnsigdf[totnsigdf$region == 'allregions',]

for (set in sets){
    setdf = setdflist[[set]]
    subdf = setdf[setdf$region == 'allregions',]
    peff = colnames(subdf)[grep('^p_.*_nb$', colnames(subdf))]
    pathstr = sub("_nb$", "", sub("^p_","", peff))

    for (pstr in pathstr){
        subdf$pvar = subdf[[paste0('p_', pstr, '_nb')]]
        subdf$lvar = subdf[[paste0('logFC_', pstr, '_nb')]]
        subdf$cvar = subdf[[paste0('col_', pstr, '_nm')]]
        subdf$lp = -log10(subdf$pvar)
        subdf = subdf[order(subdf$cvar, decreasing=F),]
        labdf = subdf[subdf$cvar != 0,]

        # Plot as a volcano plot:
        gp = ggplot(subdf, aes(lvar, lp, color=factor(cvar))) + 
            rasterise(geom_point(cex=.25), dpi=450) + 
            geom_text_repel(data=labdf, aes(lvar, lp, label=gene, color=factor(cvar)), 
                box.padding=0.075, cex=3, max.overlaps=10, min.segment.length=0.1, force=5, force_pull=0.1) + 
            scale_color_manual(values=decols) + 
            labs(x='log Fold Change', y='log10 p-value', title=paste0(set, ' - ', path, '(', pstr, ')')) + 
            geom_vline(xintercept=0, lty='dotted') + 
            scale_y_continuous(expand=c(0,0)) + 
            theme_pubr() + theme(legend.position='none')
        pltprefix = paste0(imgpref, 'volcano_', path , '.', pstr, '.', set)
        saveGGplot(gp, pltprefix, w=4, h=4)
    }
}


