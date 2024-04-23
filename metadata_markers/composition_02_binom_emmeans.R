#!/usr/bin/R
# --------------------------------------------------------
# Calculate composition differences using binomial emmeans
# Run contrasts on 2-level categorical vars, plot heatmap
# of effect sizes (log2 OR) and adjusted p-values
# As in:
# https://github.com/vals/Blog/tree/master/201127-cell-count-glm
# Updated 11/25/2021 
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(emmeans)
library(qvalue)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

# Directories:
fracdir = paste0(sdbdir, 'fractions/')
plotdir = paste0(imgdir, 'fractions/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, fracdir)
system(cmd)


# Function for plotting odds-ratio matrix:
# ----------------------------------------
# Set a single color scale:
load.colors()
# col_fun = colorRamp2(c(-1, 0, 1), c(colrb[90], "white", colrb[10]))
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))

plotCoeffPvaluesMatrix = function(cmat, pmat, pltprefix){
    require(ComplexHeatmap)
    plt = Heatmap(cmat,
                  col=col_fun,
                  use_raster=TRUE,
                  cluster_columns=FALSE,
                  width = ncol(cmat)*unit(5, "mm"), 
                  height = nrow(cmat)*unit(3, "mm"),
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


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
subsetlist = c('All_minor', 'All', 'OpcOli', unique(cellmeta$major.celltype))
remove.batches = TRUE
suff = '_subset_final_noMB'

fulldf = c()
for (subset in subsetlist){
    print(subset)
    ststr = gsub("/","_", subset)

    # Load in and process data (saves to matrices):
    commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
    source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))

    qbreg.file = paste0(fracdir, 'quasibinom_contrasts_', ststr, '.tsv')

    if (!file.exists(qbreg.file)){
        # Run quasi-binomial regression and evaluate contrasts with emmeans:
        # ------------------------------------------------------------------
        p.cut = 0.05
        pathlist = c('nrad','cogdxad','braaksc.ad', 
                     'any.nft','any.plaq_n','any.plaq_d')
        subsetdf = c()
        for (pathval in pathlist){
            print(pathval)
            lmdf = ctdf[!is.na(ctdf[[pathval]]),]
            gform = asform(c('cbind(count, other) ~ ', pathval, '* cls + cls * region * msex'))
            fit = glm(gform, lmdf, family='quasibinomial')

            emform = asform(c('revpairwise ~', pathval, '|cls'))
            emm1 <- emmeans(fit, specs=emform)
            cdf = as.data.frame(rbind(summary(emm1$contrasts, infer = TRUE, type = 'response')))
            cdf = cdf[order(cdf$odds.ratio),]
            cdf$p.adj = p.adjust(cdf$p.value, 'fdr')
            cdf$path = pathval
            cdf$col = 1 * (cdf$p.adj < p.cut) + 2 * (cdf$odds.ratio > 1 )
            subsetdf = rbind(subsetdf, cdf)

            sigcols = c('0' ='grey75','1'='royalblue','2'='grey75', '3'='indianred')
            labdf = cdf[cdf$p.adj < p.cut,]
            gplot = ggplot(cdf, aes(x=odds.ratio, y=cls, col=factor(col))) +
                geom_point() +
                geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend =cls)) +
                geom_text(data=cdf, aes(x=max(labdf$odds.ratio) * 1.1, y=cls, 
                                        label=sprintf("OR=%0.2f p=%0.1e", odds.ratio, p.adj)), col='black') +
                scale_color_manual(values=sigcols) + 
                geom_vline(xintercept=1, lty='dashed') + 
                theme_pubr() + theme(legend.position = 'none') +
                labs(x = cdf$contrast[1], y='Subtype')
            h = nrow(cdf) / 30 * 8 + 1
            ggsave(paste0(imgpref, 'oddratios_', ststr, '_by_', pathval,'.pdf'), gplot, dpi=450, units='in', width=4, height=h)
            ggsave(paste0(imgpref, 'oddratios_', ststr, '_by_', pathval,'.png'), gplot, dpi=450, units='in', width=4, height=h)
        }

        # Save table:
        write.table(subsetdf, qbreg.file, quote=F, row.names=F, sep="\t")

        # TODO: plot with splits for Exc

        # Plot table, but adjust p-value for all tests:
        subsetdf$p.adj = p.adjust(subsetdf$p.value, 'fdr')
        cmat = pivot.tomatrix(subsetdf[,c('cls','odds.ratio','path')], 'path','odds.ratio')
        pmat = pivot.tomatrix(subsetdf[,c('cls','p.adj','path')], 'path','p.adj')

        heatpref = paste0(imgpref, 'oddratios_', ststr, '_heatmap')
        plotCoeffPvaluesMatrix(log2(cmat), pmat, pltprefix=heatpref) 
    } else {
        subsetdf = read.delim(qbreg.file, header=T)

        # Plot table, but adjust p-value for all tests:
        subsetdf$p.adj = p.adjust(subsetdf$p.value, 'fdr')
        cmat = pivot.tomatrix(subsetdf[,c('cls','odds.ratio','path')], 'path','odds.ratio')
        pmat = pivot.tomatrix(subsetdf[,c('cls','p.adj','path')], 'path','p.adj')

        heatpref = paste0(imgpref, 'oddratios_', ststr, '_heatmap')
        plotCoeffPvaluesMatrix(log2(cmat), pmat, pltprefix=heatpref) 
    }
    subsetdf$subset = subset
    fulldf = rbind(fulldf, subsetdf)
}


# Plot all of these together:
fulldf$cls2 = paste0(fulldf$subset, ":", fulldf$cls)
fulldf$p.adj = p.adjust(fulldf$p.value, 'fdr')
cmat = pivot.tomatrix(fulldf[,c('cls2','odds.ratio','path')], 'path','odds.ratio')
pmat = pivot.tomatrix(fulldf[,c('cls2','p.adj','path')], 'path','p.adj')

heatpref = paste0(imgpref, 'oddratios_majorct_heatmap')
plotCoeffPvaluesMatrix(log2(cmat), pmat, pltprefix=heatpref) 
