#!/usr/bin/R
# --------------------------------------------------------
# Plot the aggregated DE results across different methods:
# Raw counts of degs
# Overlap between conditions
# Updated: 03/22/22
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
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
rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_runlist.tsv'), header=T)
rundf = rbind(rundf, read.delim(paste0(sdbdir, 'nebula_wRUV_excitatory_subsets_runlist.tsv'), header=T))
rundf$prefstr = with(rundf, paste(celltype, subtype, region, path, sep="_"))

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
                            length(list.files(path=regdir, pattern=paste0('allmethods.', x, '.merged.rda'))) })
table(rundf$merged)
head(rundf[rundf$merged == 0,])


# Select runs to use (allregions, main ct):
# -----------------------------------------
selrundf = rundf[rundf$region == 'allregions',]
selrundf$setid = with(selrundf, paste0(celltype, '_', subtype))
sets = unique(selrundf$setid)
print(sets)


# Aggregate the DEGs and results across all runs:
# -----------------------------------------------
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
if (!file.exists(fullaggrda)){
    kept.cols = c("gene","pc", "col_nm","log10p_nm","path",
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
            cat(nsig, '\n')
            # Concatenate results:
            aggdf$path = nsig[1,'path']
            setdflist[[setid]] = rbind(setdflist[[setid]], aggdf[,kept.cols])
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
totnsigdf$X1 = as.numeric(totnsigdf$X1)
totnsigdf$X2 = as.numeric(totnsigdf$X2)
umat = pivot.tomatrix(totnsigdf[,c('path','subtype','X2')], 'path','X2')
dmat = pivot.tomatrix(totnsigdf[,c('path','subtype','X1')], 'path','X1')
umat[is.na(umat)] = 0
dmat[is.na(dmat)] = 0
colnames(umat) = paste0(colnames(umat), '_up')
colnames(dmat) = paste0(colnames(dmat), '_down')

bmat = cbind(umat, -dmat)
mx = max(abs(bmat))
csplit = ifelse(1:ncol(bmat) %in% grep("_up", colnames(bmat)), 'Up','Down')

col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
plt = Heatmap(bmat, name='Number\nof DEGs',
              use_raster=TRUE,
              col=col_fun,
              column_split=csplit,
              border_gp=gpar(color='black', lwd=.5),
              width = ncol(bmat)*unit(4.5, "mm"), 
              height = nrow(bmat)*unit(2, "mm"),
              cell_fun = function(j, i, x, y, w, h, fill) {
                  ann = abs(bmat[i,j])
                  grid.text(ann, x, y, gp=gpar(fontsize=5)) })

h = 2.25 + 1 / 15 * nrow(bmat)
w = 5 + 1 / 15 * ncol(bmat)
pltprefix = paste0(imgpref, 'allmethods_ndeg.allregions.heatmap')
saveHeatmap(plt, pltprefix, w=w, h=h)

