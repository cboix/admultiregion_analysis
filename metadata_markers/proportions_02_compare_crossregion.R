#!/usr/bin/R
# ------------------------------------------------------------------
# Compare the proportions across regions for the main paper figures:
# Updated 12/01/2021 
# ------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(lme4)
library(emmeans)
library(ComplexHeatmap)
library(circlize)

# Directories:
fracdir = paste0(sdbdir, 'fractions/')
plotdir = paste0(imgdir, 'metadata/')
imgpref = paste0(plotdir, 'metadata_prop_')
cmd = paste('mkdir -p', plotdir, fracdir)
system(cmd)


# Function for plotting odds-ratio matrix:
# ----------------------------------------
# Set a single color scale:
load.colors()
col_fun = colorRamp2(c(-1, 0, 1), c(colrb[90], "white", colrb[10]))
# col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))

plotCoeffPvaluesMatrix = function(cmat, pmat, pltprefix, raster=TRUE){
    require(ComplexHeatmap)
    plt = Heatmap(cmat,
                  col=col_fun,
                  use_raster=raster,
                  cluster_columns=FALSE,
                  width = ncol(cmat)*unit(2, "mm"), 
                  height = nrow(cmat)*unit(5, "mm"),
                  cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                      p = pmat[i,j]
                      ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                      grid.text(ann, x, y)}
    )
    h = 2.25 + 2.5 / 15 * nrow(cmat)
    w = 4 + 2.5 / 15 * ncol(cmat)
    pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
    print(plt)
    dev.off()
    png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
    print(plt)
    dev.off()
}



# Test for each subtype whether there are more cells in any specific region:
# --------------------------------------------------------------------------
subsetlist = c('All_minor', 'All', 'OpcOli', unique(cellmeta$major.celltype))
remove.batches = FALSE # Need to keep SV2C, for example
suff = '_subset_final_noMB'

fulldf = c()
for (subset in subsetlist){
    cat(subset)
    ststr = gsub("/","_", subset)
    qbreg.file = paste0(fracdir, 'quasibinom_compareregions_', ststr, '.tsv')
    if (!file.exists(qbreg.file)){

        # Load in and process data (saves to matrices):
        commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
        source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))

        # Run quasi-binomial regression and evaluate contrasts with emmeans:
        # ------------------------------------------------------------------
        p.cut = 0.05
        subsetdf = c()
        for (region in reg.nomb){
            cat("\t", region)
            # Prepare dummy variable and run regression:
            ctdf$is.region = ctdf$region == region
            gform = asform(c('cbind(count, other) ~ is.region * cls + cls * msex * age_death'))
            fit = glm(gform, ctdf, family='quasibinomial')

            # Get contrasts:
            emform = asform(c('revpairwise ~ is.region|cls'))
            emm1 <- emmeans(fit, specs=emform)
            cdf = as.data.frame(rbind(summary(emm1$contrasts, infer = TRUE, type = 'response')))
            cdf = cdf[order(cdf$odds.ratio),]
            cdf$p.adj = p.adjust(cdf$p.value, 'fdr')
            cdf$region = region
            cdf$col = 1 * (cdf$p.adj < p.cut) + 2 * (cdf$odds.ratio > 1)
            subsetdf = rbind(subsetdf, cdf)
        }

        # Save table:
        write.table(subsetdf, qbreg.file, quote=F, row.names=F, sep="\t")
    } else {
        subsetdf = read.delim(qbreg.file, sep="\t")
    }
    cat("\tDone.\n")
    # Plot table, but adjust p-value for all tests:
    subsetdf$p.adj = p.adjust(subsetdf$p.value, 'fdr')
    cmat = pivot.tomatrix(subsetdf[,c('cls','odds.ratio','region')], 'region','odds.ratio')
    pmat = pivot.tomatrix(subsetdf[,c('cls','p.adj','region')], 'region','p.adj')

    heatpref = paste0(imgpref, 'oddratios_compareregions_', ststr, '_heatmap')
    plotCoeffPvaluesMatrix(log2(cmat), pmat, pltprefix=heatpref, raster=FALSE) 
}


