#!/usr/bin/R
# --------------------------------------------------------
# Plot the aggregated DE results across different methods:
# - Adding DESeq2 results, + make aggregate datasets
# Updated: 11/02/23
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
print(version)
options(width=170)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Get a list of all differential runs:
# ------------------------------------
rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_joint_runlist.tsv'), header=T)
rundf$prefstr = with(rundf, gsub("[()/]", "_", 
        paste(celltype, subtype, region, path, sep="_")))

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
    outfile = paste0(regdir, 'allmethods_ds2.', x, '.merged.rda')
    1 * file.exists(outfile) })
table(rundf$merged)
head(rundf[rundf$merged == 0,], 20)


# Select runs to use (allregions, main ct):
# -----------------------------------------
path = 'nrad'
for (path in unique(rundf$path)){
    print(path)
    selrundf = rundf[rundf$path == path,]
    selrundf$setid = with(selrundf, paste0(celltype, '_', subtype))
    sets = unique(selrundf$setid)
    sets = c(sets[grep("^Exc_",sets, invert=TRUE)], 'Exc_Exc')
    selrundf = selrundf[selrundf$setid %in% sets,]
    print(sets)


    # Aggregate the DEGs and results across all regional runs for NFT:
    # ----------------------------------------------------------------
    mstr = paste0('allmethods_ds2.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    if (!file.exists(fullaggrda)){
        kept.cols = c("gene","pc", "col_nm","log10p_nm",
            "path","region",
            "logFC_nb","p_nb","padj_nb","col_nb",
            "coef_mast","p_mast","padj_mast","col_mast",
            "log2FC_ds","p_ds","padj_ds","col_ds")

        setdflist = lapply(sets, function(x){})
        names(setdflist) = sets
        totnsigdf = NULL
        for (i in 1:nrow(selrundf)){
            prefstr = selrundf$prefstr[i]
            setid = selrundf$setid[i]
            aggrda = paste0(regdir, 'allmethods_ds2.', prefstr, '.merged.rda')
            if (file.exists(aggrda)){
                load(aggrda)  # Loads aggdf, nsig
                cat(nsig, '\n')
                # Concatenate results:
                aggdf$path = nsig[1,'path']
                aggdf$region = nsig[1,'region']
                setdflist[[setid]] = rbind(setdflist[[setid]], aggdf[,kept.cols])
                # Pad nsig:
                nsigdf = data.frame(nsig)
                if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
                if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
                nsigdf$run = 'Nebula_MAST'

                # Get number of significant DEseq results:
                # NOTE: Separate from the cell-level models:
                densig = spread(as.data.frame(table(aggdf$col_ds)), Var1, Freq)
                names(densig) = paste0('X', names(densig))
                if (!("X1" %in% colnames(densig))){ densig$X1 = 0}
                if (!("X2" %in% colnames(densig))){ densig$X2 = 0}
                densig$run = 'DESeq2'
                cn = colnames(nsigdf)[!(colnames(nsigdf) %in% colnames(densig))]
                densig = merge(nsigdf[,cn], densig)
                cat('DESeq2:\t', as.character(densig[1,]), '\n')
                # Add to output:
                nsigdf = rbind(nsigdf, densig[,colnames(nsigdf)])
                totnsigdf = rbind(totnsigdf, nsigdf)
            }
        }
        save(setdflist, totnsigdf, file=fullaggrda)
    } else {
        load(fullaggrda)
    }


    # Plot total number of significant genes as a heatmap:
    # NOTE: Plotting only DESeq2 here, as others are done in original script
    # ----------------------------------------------------------------------
    totnsigdf = totnsigdf[totnsigdf$run == 'DESeq2',]
    totnsigdf$X1 = as.numeric(totnsigdf$X1)
    totnsigdf$X2 = as.numeric(totnsigdf$X2)
    umat = pivot.tomatrix(totnsigdf[,c('region','subtype','X2')], 'region','X2')
    dmat = pivot.tomatrix(totnsigdf[,c('region','subtype','X1')], 'region','X1')
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
    pltprefix = paste0(imgpref, 'allmethods_ndeg.regional_', path, '.heatmap')
    saveHeatmap(plt, pltprefix, w=w, h=h)

}

