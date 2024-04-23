#!/usr/bin/R
# --------------------------------------------------------
# Score the pseudobulks by modules, plot against pathology
# Updated 01/18/2023
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
options(width=175)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Directories:
moddir = paste0(sdbdir, 'modules/')
regdir = paste0(sdbdir, 'dereg/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'modorder_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Score modules vs. stages (by region)
# ------------------------------------
runset = 'Ast'
region = 'PFC'
subdf = runscdf[runscdf$runset == runset,]
mods = c(28, 0, 24, 7, 27, 15, 6, 1, 11, 12)
subdf = subdf[subdf$module %in% mods,]
subdf$projid = sub("-.*","", subdf$pr)
subdf$region = sub(".*-","", subdf$pr)
if (region == 'HC'){
    ptail = 'hip'
} else if (region == 'PFC'){
    ptail = 'mf'
} else {
    ptail = tolower(region)
}
regpath = paste0(c('nft_','plaq_n_', 'plaq_d_'), ptail)
advars = c('cogdxad','braaksc','nrad','niareagansc', regpath)
subdf = merge(subdf,
    unique(metadata[metadata$region == region, c('projid',advars)]))
aggdf = merge(agg.rename(score ~ rm + region, subdf, mean, 'mean'),
    agg.rename(score ~ rm + region, subdf, sd, 'sd'))
subdf = merge(subdf, aggdf)
regdf = subdf[subdf$region == 'EC',]


# Just PFC:
for (path in regpath){
    print(path)
    regdf$path = log1p(regdf[[path]])
    gp = ggplot(regdf, aes(path, score, color=rm)) +
        facet_wrap(~rm, scales='free') + 
        geom_point(cex=.5) + 
        labs(x=paste0('log1p(', path, ')'), y='Module score') + 
        geom_smooth() + 
        theme_pubr()
    pltprefix = paste0(imgpref, runset, '_scatter_score_', region, '_', path)
    saveGGplot(gp, pltprefix, h=5, w=6)

    gp = ggplot(regdf, aes(path, score / mean, color=rm)) +
        geom_point(cex=.5) + 
        labs(x=paste0('log1p(', path, ')'), y='Module score') + 
        geom_smooth(se=FALSE) + 
        theme_pubr()
    pltprefix = paste0(imgpref, runset, '_scatter_score_stack_', region, '_', path)
    saveGGplot(gp, pltprefix, h=5, w=6)
}



# Aggregate + summarize:
# ----------------------
subdf$wscore = subdf$score * subdf$ncell / subdf$mean
aggdf = aggregate(cbind(ncell, wscore) ~ module + rm + projid, subdf, sum)
aggdf$score = aggdf$wscore / aggdf$ncell
aggdf = merge(aggdf, aggregate(cbind(tangles, amyloid) ~ projid, metadata, mean))

gp = ggplot(aggdf, aes(log1p(tangles), score, color=rm)) +
    geom_point(cex=.5) + 
    labs(y='Module score') + 
    geom_smooth(se=FALSE) + 
    theme_pubr()
pltprefix = paste0(imgpref, runset, '_scatter_score_stack_tangles')
saveGGplot(gp, pltprefix, h=5, w=6)

gp = ggplot(aggdf, aes(log1p(amyloid), score, color=rm)) +
    geom_point(cex=.5) + 
    labs(y='Module score') + 
    geom_smooth(se=FALSE) + 
    theme_pubr()
pltprefix = paste0(imgpref, runset, '_scatter_score_stack_amyloid')
saveGGplot(gp, pltprefix, h=5, w=6)




# Can we summarize these curves in one statistic (e.g. peak?):
# Then plot stats for each region
# ------------------------------------------------------------





# Global, across regions:
gp = ggplot(subdf, aes(gpath, score, color=region)) +
    facet_wrap(~rm, scales='free') + 
    geom_point() + 
    geom_smooth() + 
    theme_pubr()
pltprefix = paste0(imgpref, runset, '_scatter_score_gpath')
saveGGplot(gp, pltprefix, h=5, w=6)


# Mean locs:
# ----------
# All samp:
subdf$wscore = subdf$score * subdf$ncell
aggdf = aggregate(cbind(ncell, wscore) ~ module + rm + projid, subdf, sum)
aggdf$score = aggdf$wscore / aggdf$ncell

# By region:


# Plot heatmap:


# Plot overlapping densities



# Boxplots:


