#!/usr/bin/R
# --------------------------------------------------------
# Plot the aggregated DE results across different methods:
# for specific examples (e.g. P2RY12 for reviewer)
# Updated: 11/30/23
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)

library(ggrepel)
library(ggplot2)
library(ggpubr)
print(version)
options(width=170)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Load all of the runs:
# ---------------------
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
sets = c('Ast_Ast','Exc_Exc','Inh_Inh','Opc_Opc','Mic_Immune_Mic','Oli_Oli')
kept.cols = c('gene','col_nm','path','region', 'logFC_nb', 'set')

fulldf = c()
for (path in pathlist){
    print(path)
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)
    for (set in sets){
        setdf = setdflist[[set]]
        setdf$set = set
        fulldf = rbind(fulldf, setdf[, kept.cols])
    }
}


# Subset to a specific gene, celltype: 
# -------------------------------------
set = 'Mic_Immune_Mic'
gene = 'P2RY12'

subdf = fulldf[(fulldf$set == set) & (fulldf$gene == gene),]
subdf = subdf[subdf$path %in% c('nft','plaq_n','plaq_d'),]
cmat = pivot.tomatrix(subdf[,c('path','region','logFC_nb')], 'region','logFC_nb')
pmat = pivot.tomatrix(subdf[,c('path','region','col_nm')], 'region','col_nm')

ux = 1.5
ht = Heatmap(cmat, 
    use_raster=FALSE, 
    name='logFC_nb',
    # col=col_fun,
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
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
        if (!is.na(p) & (p != 0)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.25))
        }
    })

pltprefix = paste0(imgpref, 'regionalDEG_example.', set, '.', gene)
h = 1.5 + 1 / 15 * nrow(cmat)
w = 1.5 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)

