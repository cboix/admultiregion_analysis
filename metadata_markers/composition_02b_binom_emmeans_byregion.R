#!/usr/bin/R
# --------------------------------------------------------
# Calculate composition differences using binomial emmeans
# Run contrasts on 2-level categorical vars, plot heatmap
# of effect sizes (log2 OR) and adjusted p-values
#
# Here, we compare region differences in major cell types
# Updated 12/20/2021 
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

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Function for plotting odds-ratio matrix:
# ----------------------------------------
# Set a single color scale:
load.colors()
# col_fun = colorRamp2(c(-1, 0, 1), c(colrb[90], "white", colrb[10]))
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))

plotCoeffPvaluesMatrix = function(cmat, pmat, pltprefix){
    require(ComplexHeatmap)
    ux = 1.5
    plt = Heatmap(cmat,
                  col=col_fun,
                  use_raster=TRUE,
                  cluster_columns=FALSE,
                  cluster_rows=FALSE,
                  width = ncol(cmat)*unit(ux * 1.5, "mm"), 
                  height = nrow(cmat)*unit(ux, "mm"),
                  border_gp = gpar(col='black', lwd=.5),
                  cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                      p = pmat[i,j]
                      ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                      grid.text(ann, x, y, gp=gpar(fontsize=5))}
    )
    h = 2.25 + 1 / 15 * nrow(cmat)
    w = 4 + 1 / 15 * ncol(cmat)
    saveHeatmap(plt, pltprefix, w=w, h=h)
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
    qbreg.file = paste0(fracdir, 'quasibinom_contrasts_', ststr, '_byregion.tsv')

    # Add the early-mid + mid-late braak stages:
    ctdf$braaksc.stages = ifelse(ctdf$braaksc <= 4, ifelse(ctdf$braaksc <= 2, 'Early', 'Mid'), 'Late')
    table(ctdf$braaksc, ctdf$braaksc.stages) # Check ok
    # ctdf$braaksc.stages = ifelse(ctdf$braaksc <= 4, ifelse(ctdf$braaksc <= 2, 'Early', 'Mid'), 'Late')

    # Set NA to not eval late or early stages:
    ctdf$braaksc.em = ctdf$braaksc.stages
    ctdf$braaksc.ml = ctdf$braaksc.stages
    ctdf$braaksc.em[ctdf$braaksc.em == 'Late'] = NA
    ctdf$braaksc.ml[ctdf$braaksc.ml == 'Early'] = NA

    countfile = paste0(fracdir, 'counts_', ststr, '.mat.tsv')
    metafile = paste0(fracdir, 'counts_', ststr, '.meta.tsv')
    if (!file.exists(countfile)){
        # Make counts and metadata, separately:
        ctdf$pr = with(ctdf, paste0(projid, '_', region))
        cmat = pivot.tomatrix(ctdf[,c('pr', 'cls', 'count')], 'cls', 'count')
        write.table(cmat, countfile, quote=F, sep="\t")
        cmeta = ctdf
        cmeta[, c('cls', 'major.celltype', 'count', 'other', 'total')] = NULL
        cmeta = unique(cmeta)
        rownames(cmeta) = cmeta$pr
        cmeta = cmeta[rownames(cmat),]
        write.table(cmeta, metafile, quote=F, row.names=F, sep="\t")
    }

    pathlist = c('nrad','cogdxad','braaksc.ad', 'braaksc.em', 'braaksc.ml')
    # 'any.nft','any.plaq_n','any.plaq_d')
    if (!file.exists(qbreg.file)){
        fulldf = c()
        for (region in reg.nomb){
            print(region)
            # Run quasi-binomial regression and evaluate contrasts with emmeans:
            # ------------------------------------------------------------------
            p.cut = 0.05
            subsetdf = c()
            for (pathval in pathlist){
                print(pathval)
                # Remove NAs: Early for ml, Late for em (braak staging)
                lmdf = ctdf[!is.na(ctdf[[pathval]]),]  
                lmdf = lmdf[lmdf$region == region,]
                gform = asform(c('cbind(count, other) ~ ', pathval, '* cls + cls * msex'))
                fit = glm(gform, lmdf, family='quasibinomial')

                emform = asform(c('revpairwise ~', pathval, '|cls'))
                emm1 <- emmeans(fit, specs=emform)
                cdf = as.data.frame(rbind(summary(emm1$contrasts, infer = TRUE, type = 'response')))
                cdf = cdf[order(cdf$odds.ratio),]
                cdf$p.adj = p.adjust(cdf$p.value, 'fdr')
                cdf$path = pathval
                cdf$col = 1 * (cdf$p.adj < p.cut) + 2 * (cdf$odds.ratio > 1 )
                cdf$region = region
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
            fulldf = rbind(fulldf, subsetdf)
        }
        # Save table:
        write.table(fulldf, qbreg.file, quote=F, row.names=F, sep="\t")
    } else {
        fulldf = read.delim(qbreg.file, header=T)
    }


    # Plot table:
    fulldf$p.adj = p.adjust(fulldf$p.value, 'fdr')
    fulldf$cls2 = paste0(fulldf$cls, ":", fulldf$region)

    # Set low counts to NA:
    nregdf = aggregate(count ~ cls + region, ctdf, sum)
    ct.cutoff = 500
    fulldf = merge(fulldf, nregdf, all.x=TRUE)

    for (path in pathlist){
        subdf = fulldf[fulldf$path == path,]
        subdf$p.adj = p.adjust(subdf$p.value, 'fdr')
        subdf$odds.ratio[subdf$count < ct.cutoff] = NA
        subdf = subdf[!is.na(subdf$odds.ratio),]
        subdf = subdf[order(subdf$p.value),]
        cmat = pivot.tomatrix(subdf[,c('cls','odds.ratio','region')], 'region','odds.ratio')
        pmat = pivot.tomatrix(subdf[,c('cls','p.adj','region')], 'region','p.adj')
        pmat[is.na(pmat)] = 1

        heatpref = paste0(imgpref, 'oddratios_', ststr, '_heatmap_byregion_', path)
        plotCoeffPvaluesMatrix(log2(cmat), pmat, pltprefix=heatpref) 
    }

}

