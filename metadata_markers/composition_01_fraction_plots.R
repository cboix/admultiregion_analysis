#!/usr/bin/R
# ---------------------------------------------------
# Plot composition differences - boxplots + barplots:
# Updated 11/25/2021 
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)

library(ggplot2)
library(ggrepel)
library(ggpubr)

# Directories:
plotdir = paste0(imgdir, 'fractions/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir)
system(cmd)


subsetlist = c('All_minor', 'All', unique(cellmeta$major.celltype))
remove.batches = TRUE
suff = '_subset_final_noMB'

for (subset in subsetlist){
    print(subset)
    ststr = gsub("/","_", subset)

    use.cols = c(tcols, major.col)

    # Load in and process data (saves to matrices):
    commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
    source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))


    # Plot agnostic to AD across regions, etc:
    aggdf = aggregate(count ~ cls + region, ctdf, sum)
    totdf = aggregate(count ~ cls, ctdf, sum)
    totdf$region = 'All'
    aggdf = rbind(aggdf, totdf[,colnames(aggdf)])
    aggdf$region = factor(aggdf$region, levels=reg.order)

    gplot = ggplot(aggdf, aes(region, count, fill=cls)) + 
        geom_bar(position='fill', stat='identity') + 
        scale_fill_manual(values=use.cols) + 
        theme_pubr() + 
        scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
        labs(x='Region', y='Percentage') + theme(legend.position='none')
    ggsave(paste0(imgpref, 'fractions_', ststr, suff, '.png'), gplot, dpi=450, units='in', width=5, height=6)
    ggsave(paste0(imgpref, 'fractions_', ststr, suff, '.pdf'), gplot, dpi=450, units='in', width=5, height=6)


    # Plot w.r.t. to AD (NIA-Reagan 1-2 vs. 3-4):
    aggdf = aggregate(count ~ cls + region + nrad, ctdf, sum)
    totdf = aggregate(count ~ cls + nrad, ctdf, sum)
    totdf$region = 'All'
    aggdf = rbind(aggdf, totdf[,colnames(aggdf)])
    aggdf$region = factor(aggdf$region, levels=reg.order)

    gplot = ggplot(aggdf, aes(nrad, count, fill=cls)) + 
        facet_wrap(~region, nrow=1) + 
        geom_bar(position='fill', stat='identity') + 
        scale_fill_manual(values=use.cols) + 
        theme_pubr() + 
        scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
        labs(x='Region', y='Percentage') + theme(legend.position='none')
    ggsave(paste0(imgpref, 'fractions_', ststr, suff, '_ad.pdf'), gplot, dpi=450, units='in', width=7, height=6)
    ggsave(paste0(imgpref, 'fractions_', ststr, suff, '_ad.png'), gplot, dpi=450, units='in', width=7, height=6)


    # Plot by individual / region:
    gplot = ggplot(ctdf, aes(region, count / total, alpha=nrad, fill=cls)) + 
        facet_wrap(~cls, nrow=1) + 
        geom_boxplot(outlier.shape=NA) + 
        scale_fill_manual(values=use.cols) + 
        geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
        theme_pubr() + 
        scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
        stat_compare_means(hide.ns=TRUE, label='p.format') + 
        labs(x='Region', y='Percentage') + theme(legend.position='none')
    ggsave(paste0(imgpref, 'boxplots_', ststr, suff, '.png'), gplot, dpi=450, units='in', width=8, height=4)
    ggsave(paste0(imgpref, 'boxplots_', ststr, suff, '.pdf'), gplot, dpi=450, units='in', width=8, height=4)

}

