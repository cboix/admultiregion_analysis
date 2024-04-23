#!/usr/bin/R
# ---------------------------------------------------
# Plot composition differences in a heatmap for Exc:
# Updated 03/16/2022
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)

library(ggplot2)
library(ggrepel)
library(ggpubr)

library(ComplexHeatmap)
library(circlize)

# Directories:
plotdir = paste0(imgdir, 'fractions/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


subsetlist = c('All_minor', 'All', unique(cellmeta$major.celltype))
remove.batches = TRUE
suff = '_subset_final_noMB'

# for (subset in subsetlist){
subset = 'Exc'
print(subset)
ststr = gsub("/","_", subset)
use.cols = c(tcols, major.col)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))


# Calculate average fractions in early AD:
# ----------------------------------------
early.sets = c(0,1,2)
aggdf = aggregate(cbind(count, total) ~ cls + projid + region + braaksc + nrad, ctdf, sum)
aggdf$frac = aggdf$count / aggdf$total
nondf = aggregate(frac ~ cls + region, aggdf[aggdf$braaksc %in% early.sets,], mean)

# Select sets for plotting and subset aggregated:
kdf = merge(data.frame(cls=exc.sets[['CTXneurons']]), data.frame(region=c('AG','MT','PFC')))
kdf = rbind(kdf, data.frame(cls=exc.sets[['ECneurons']], region=c('EC')))
kdf = rbind(kdf, data.frame(cls=exc.sets[['HCneurons']], region=c('HC')))
kdf = rbind(kdf, data.frame(cls=exc.sets[['THneurons']], region=c('TH')))

aggdf = merge(kdf, aggdf)
nondf = merge(kdf, nondf)
aggdf$cr = with(aggdf, paste0(cls,'@',region))
nondf$cr = with(nondf, paste0(cls,'@',region))

# nondf$frac = nondf$count / nondf$total
rownames(nondf) = nondf$cr

aggdf$frac = aggdf$count / aggdf$total
amat = pivot.tomatrix(aggdf[,c('cr','projid','frac')], 'projid', 'frac')


# Metadata for individuals:
# -------------------------
umeta = unique(metadata[metadata$region == 'PFC',
    c('projid','nrad','cogdxad','cogdx', 'niareagansc', 'age_death',
        'msex', 'braaksc', 'Apoe_e4','gpath','tangles','amyloid', 'cogn_global_lv')])
rownames(umeta) = umeta$projid
umeta = umeta[as.character(colnames(amat)),]

umeta = umeta[order(umeta$braaksc),]
umeta = umeta[order(-umeta$niareagansc),]
amat = amat[,as.character(umeta$projid)]

age.col_fun = colorRamp2(range(umeta$age_death), c("white", "slateblue")) 
pmi.col_fun = colorRamp2(c(2, 15), c("white", "indianred")) 
gpath.col_fun = colorRamp2(c(0, max(umeta$gpath)), c("white", "indianred")) 
gcog.col_fun = colorRamp2(c(min(umeta$cogn_global_lv), max(umeta$cogn_global_lv)),
    c("red", 'white')) 
# mat.col_fun = colorRamp2(c(0, max(pmat, na.rm=T)), c("white", "blue")) 

# Make metadata annotation:
ux = 1.5
ha = HeatmapAnnotation(
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(ux, 'mm'),
    Sex=ifelse(umeta$msex == 0, 'female','male'), 
    # PMI=umeta$pmi,
    Age=umeta$age_death,
    Apoe_e4=umeta$Apoe_e4,
    GPath=umeta$gpath,
    Braak=umeta$braaksc,
    AD=umeta$niareagansc,
    Cognition=umeta$cogdx,
    GlobalCog=umeta$cogn_global_lv,
    col=list(AD=colvals[['niareagansc']],
        Apoe_e4=c('no'='grey95','yes'='grey70'),
        Age=age.col_fun,
        GPath=gpath.col_fun,
        GlobalCog=gcog.col_fun,
        # PMI=pmi.col_fun,
        Braak=colvals[['braaksc']],
        Cognition=colvals[['cogdx']],
        Sex=colvals[['sex']]), which='column')


# Plot heatmap:
# -------------
mx = 4
col_fun = colorRamp2(c(-mx, 0, mx), c('blue', "white", 'red'))
pltmat = amat
marg = nondf[rownames(pltmat), 'frac']
pltmat = log2(sweep(pltmat, 1, marg, '/'))
pltmat[is.infinite(pltmat)] = -mx

ux = 1.5
ht = Heatmap(pltmat, 
    col=col_fun, 
    width=ncol(pltmat) * unit(ux, 'mm'),
    height=nrow(pltmat) * unit(ux, 'mm'),
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    border_gp=gpar(lwd=.5, color='black'),
    bottom_annotation=ha,
    row_split = sub('.*@', '', rownames(pltmat))
)

w = 3.5 + ncol(pltmat) / 15
h = 3.5 + nrow(pltmat) / 15
pltprefix = paste0(imgpref, 'fractions_', ststr, suff, '_heatmap.png')
saveHeatmap(ht, pltprefix, w=w, h=h)


# Plot cortical neuron %s as boxplots, by late AD:
# ------------------------------------------------
ctxdf = ctdf[ctdf$region %in% c('AG','MT','PFC'),]
ctxdf = ctxdf[ctxdf$cls %in% exc.sets[['CTXneurons']],]
gplot = ggplot(ctxdf, aes(cls, count / total, fill=braaksc.ad)) + 
    facet_wrap(~region, scales='free', nrow=1) + 
    geom_boxplot(outlier.shape=NA) + 
    scale_fill_manual(values=colvals[['cogdxad']]) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.5) +
    theme_pubr() + coord_flip() + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    stat_compare_means(hide.ns=TRUE, label='p.format', label.y.npc=0.8) + 
    labs(x='Subtype', y='Percentage') + theme(legend.position='none')
pltprefix = paste0(imgpref, 'boxplots_', ststr, suff, '_ctx_braakscad.png')
saveGGplot(gplot, pltprefix, w=15, h=4)

