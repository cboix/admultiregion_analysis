#!/usr/bin/R
# ---------------------------------------------------------
# Calculate and plot the module enrichments for covariates:
# Updated 11/28/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


library(tidyr)
library(viridis)
library(ggpubr)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)


# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Load in the full pseudobulk data for these subtypes:
# NOTE: From modules_degenr_03_plot_pseudobulk_modules.R
# ------------------------------------------------------
useset = 'coregenes_'
scores.file = paste0(moddir, 'module_pseudobulk_scores_', 
    useset, fullpref, '.tsv.gz')
scdf = read.delim(gzfile(scores.file), sep="\t")

# Merge rind to use as indexing:
scdf = merge(scdf, metadata[,c('projid','region','rind')])

# Merge major cell type for 'All'
if (runset == 'All'){
    ctmap = unique(cellmeta[,c('major.celltype','cell_type_high_resolution')])
    scdf = merge(scdf, ctmap)
}


# Run hypergeometric tests to score each attribute vs. each module:
# -----------------------------------------------------------------
mnames = sort(unique(scdf$mname))
covarcols = c('msex','Apoe_e4', 'age_bin', 'apoe_genotype')
advars = c('niareagansc','nrad','cogdxad','braaksc56', 'ceradsc','braaksc','cogdx')
# pathvars = c('nft','plaq_n','plaq_d')
# TODO: Add pathvars and cognition variables

metadata$age_bin = metadata$age_death > 85
metadata$braaksc56 = 'CTRL'
metadata$braaksc56[metadata$braaksc %in% c(5,6)] = 'AD'
metadata$braaksc56 = factor(metadata$braaksc56, levels=c('CTRL','AD'))

# Merge onto counts table:
scdf = merge(scdf, unique(metadata[,c('rind', covarcols, advars)]), all.x=TRUE)
# ctdf = merge(ctdf, pqdf, all.x=TRUE)
# ctdf$nft = log1p(ctdf$nft)
# ctdf$plaq_n = log1p(ctdf$plaq_n)
# ctdf$plaq_d = log1p(ctdf$plaq_d)

covarlist = c('cell_type_high_resolution', advars, 'region',
              'msex','Apoe_e4','age_bin','apoe_genotype')
if (runset == 'All'){ covarlist = c('major.celltype', covarlist) }
hgdf = NULL
for (mn in mnames){
    cat(mn, '\n')
    subdf = scdf[scdf$mname == mn,]
    subdf$z = scale(subdf$score)
    subdf$dummy = 1 * (subdf$z > 1)
    covar = 'cell_type_high_resolution'
    for (covar in covarlist){
        lvls = unique(subdf[[covar]])
        # For each level, test:
        for (lv in lvls){
            draws = table(subdf[, 'dummy'])
            subhits = table(subdf[subdf[[covar]] == lv, 'dummy'])
            hits = draws * 0
            hits[names(subhits)] = subhits
            df = cbind(q=hits, draw=draws,
                       m=nrow(subdf[subdf[[covar]] == lv,]),
                       N=nrow(subdf))
            pout <- apply(df, 1, run.hyper)
            df = data.frame(df, p.value=pout, cls=rownames(df), covariate=covar, level=lv, mname=mn)
            rownames(df) = NULL
            hgdf = rbind(hgdf, df)
        }
    }
}
hgdf = hgdf[order(hgdf$p),]
head(hgdf, 20)


# Plot the enrichments:
# ---------------------
hgdf$log10p = -log10(hgdf$p.value)
hgdf$log10p[hgdf$p.value < 10^-12] <- 12  # Cap inf
hgdf$cl = paste0(hgdf$covariate, '\n', hgdf$level)
hgwide = spread(hgdf[hgdf$cls == 1, c('mname','covariate','cl','log10p')], mname, log10p)
hgwide$cl = factor(hgwide$cl, unique(hgdf$cl))
hgwide = hgwide[order(hgwide$cl),]
hgmat = as.matrix(hgwide[,-c(1,2)])
CUTOFF=5
hgmat[hgmat > CUTOFF] <- CUTOFF
hgmat[hgmat < 2] <- 0
rownames(hgmat) = hgwide$cl
mcols = colnames(hgmat)
mcols = mmap$mname[mmap$mname %in% mcols]
hgmat = hgmat[,mcols]

plt = Heatmap(t(hgmat),
              col=colb, 
              use_raster=T,
              border_gp=gpar('black'),
              cluster_rows=FALSE,
              column_split=hgwide$covariate)


# Turn into an annotation-style object by keeping top label for each:
# -------------------------------------------------------------------
labdf = hgdf[hgdf$cls == 1,]
labdf = merge(labdf, aggregate(p.value ~ covariate + mname, labdf, min))
labdf = labdf[labdf$p.value < 1e-2,]
labdf = aggregate(level ~ mname + covariate, labdf, function(x){x[1]})
labdf = spread(labdf, covariate, level)
labdf[is.na(labdf)] = 'NA'

# Save the annotation table:
hgann.file = paste0(moddir, 'module_covariate_hgenr_', useset, fullpref, '.tsv')
write.table(labdf, hgann.file, quote=F, row.names=F, sep="\t")

hgdf.file = paste0(moddir, 'module_covariate_hgdf_', useset, fullpref, '.tsv')
hgdf$cl = NULL
write.table(hgdf, hgdf.file, quote=F, row.names=F, sep="\t")

