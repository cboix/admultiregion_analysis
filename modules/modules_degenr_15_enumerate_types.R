#!/usr/bin/R
# ------------------------------------
# Enumerate types of modules captured:
# Updated 03/23/2022
# ------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(PRROC)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
options(width=150)

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
resdir = paste0(sdbdir, 'modules/resources/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'module_panels_')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
system(cmd)


# List of arguments to load:
# --------------------------
graph_id = 'boot'
runlist = c('Ast','Oli','Opc','Mic_Immune',
    'Vasc_Epithelia', 'Inh', 'Exc', 'All')


full.hgdf = c()
for (runset in runlist){
    # Load in metadata for modules:
    # -----------------------------
    print(runset)
    commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, FALSE, FALSE)}
    source(paste0(sbindir, 'modules/load_modules_degenr.R'))

    # Get hypergeometric enrichment-based covariate scores:
    # -----------------------------------------------------
    useset = 'coregenes_'
    hgdf.file = paste0(moddir, 'module_covariate_hgdf_', useset, fullpref, '.tsv')
    print(hgdf.file)
    hgdf = read.delim(hgdf.file, header=T)
    hgdf$log10p = -log10(hgdf$p.value)

    unique(hgdf$covariate)
    sub.hgdf = merge(hgdf, aggregate(p.value ~ covariate + level + mname, 
            hgdf, min))
    sub.hgdf$log10p = (2 * sub.hgdf$cls - 1) * sub.hgdf$log10p
    sub.hgdf$runset = runset
    full.hgdf = rbind(full.hgdf, sub.hgdf)
}

mingenes = 10
full.hgdf$ng = with(full.hgdf, sub(".*\\(","", sub(" gene.*","", mname)))
full.hgdf$ng = as.numeric(full.hgdf$ng)

# Count number of modules:
uqdf = unique(full.hgdf[,c('runset','mname','ng')])
totdf = agg.rename(mname ~ runset, uqdf[uqdf$ng >= mingenes,], length, 'total')

# Plot number of modules per runset:
major.col2 = major.col
names(major.col2) = sub("/","_", names(major.col))
major.col2['All'] = 'grey80'
gp = ggplot(totdf, aes(runset, total, fill=runset)) + 
    geom_bar(stat='identity') + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=major.col2) + 
    theme_pubr() + coord_flip() + theme(legend.position='right')
pltprefix = paste0(imgpref, 'celltype_total_barplot')
saveGGplot(gp, pltprefix, w=8, h=3)


# Get significant + break ties for enrichment double counting:
# ------------------------------------------------------------
sig.hgdf = full.hgdf[abs(full.hgdf$log10p) > 3,]
sig.hgdf = sig.hgdf[sig.hgdf$cls == 1,]
sig.hgdf = merge(sig.hgdf, aggregate(p.value ~ covariate + mname + runset + ng, sig.hgdf, min))
# sig.hgdf = merge(sig.hgdf, aggregate(cls ~ covariate + mname + runset, sig.hgdf, max))
sig.hgdf$log2FC = with(sig.hgdf, log2((q / draw) / (m /  N)))
subsig.hgdf = sig.hgdf[sig.hgdf$ng >= mingenes,]

# Barplot of enrichments (by level and by runset), merge cell types:
hg.covars = c('msex','Apoe_e4','region','cell_type_high_resolution',
    'nrad','braaksc56','cogdxad')
aggdf = agg.rename(p.value ~ level + covariate, subsig.hgdf, length, 'count')
aggdf$level[aggdf$covariate == 'cell_type_high_resolution'] = 'Merge'
aggdf = aggregate(count ~ level + covariate, aggdf, sum)
aggdf = aggdf[aggdf$covariate %in% hg.covars,]

ctdf = agg.rename(p.value ~ runset + covariate, subsig.hgdf, length, 'count')
ctdf = ctdf[ctdf$covariate %in% hg.covars,]

cv.cols = c(colvals[['cogdxad']], colvals[['Apoe_e4']], reg.cols,'Merge'='grey80','0'='pink', '1'='lightblue', 'no'='grey80', 'yes'='slateblue')
gp = ggplot(aggdf, aes(covariate, count, fill=level)) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=cv.cols) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() + coord_flip() + theme(legend.position='right')
pltprefix = paste0(imgpref, 'covar_barplot')
saveGGplot(gp, pltprefix, w=8, h=3)

ctdf = merge(ctdf, totdf)
major.col2 = major.col
names(major.col2) = sub("/","_", names(major.col))
major.col2['All'] = 'grey50'
gp = ggplot(ctdf, aes(covariate, count, fill=runset)) + 
    geom_bar(stat='identity', position='fill') + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=major.col2) + 
    theme_pubr() + coord_flip() + theme(legend.position='right')
pltprefix = paste0(imgpref, 'celltype_barplot')
saveGGplot(gp, pltprefix, w=8, h=3)


# Assign to top enrichment:
subdf = subsig.hgdf[subsig.hgdf$covariate %in% c('major.celltype', hg.covars),]
topdf = merge(subdf, aggregate(p.value ~ mname + runset, subdf, min))
topdf = topdf[order(topdf$p.value),]
topdf = merge(topdf, aggregate(covariate ~ mname + runset, topdf, function(x){head(x,1)}))
topdf = topdf[order(topdf$p.value),]

aggdf = agg.rename(p.value ~ covariate + runset, topdf, length, 'count')
enrdf = agg.rename(count ~ runset, aggdf, sum, 'enr')
enrdf = merge(enrdf, totdf)
enrdf$count = enrdf$total - enrdf$enr
enrdf$covariate = 'None'
aggdf = rbind(aggdf, enrdf[, colnames(aggdf)])
# pmat = pivot.tomatrix(aggdf, 'covariate', 'count')

# Plot number of modules with max val per runset:
# major.col2 = major.col
# names(major.col2) = sub("/","_", names(major.col))
# major.col2['All'] = 'grey80'
gp = ggplot(aggdf, aes(runset, count, fill=covariate)) + 
    geom_bar(stat='identity') + 
    scale_y_continuous(expand=c(0,0)) + 
    # scale_fill_manual(values=major.col2) + 
    theme_pubr() + coord_flip() + theme(legend.position='right')
pltprefix = paste0(imgpref, 'celltype_total_barplot')
saveGGplot(gp, pltprefix, w=8, h=3)


# Give the most associated one for each (by sig. log2FC):
# -------------------------------------------------------
# min.hgdf = merge(subsig.hgdf, aggregate(p.value ~ covariate, subsig.hgdf, min))
max.hgdf = merge(subsig.hgdf, aggregate(log2FC ~ covariate, subsig.hgdf, max))

