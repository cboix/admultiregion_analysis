#!/usr/bin/R
# ------------------------------------------------------------------------
# Plot the base region proportions figures for the atlas part of the paper
# Updated 11/04/2021 
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(qvalue)
library(lme4)

# Directories:
plotdir = paste0(imgdir, 'metadata/')
imgpref = paste0(plotdir, 'metadata_prop_')
cmd = paste('mkdir -p', plotdir)
system(cmd)


# Load in the abundance data:
# ---------------------------
subset = 'Ast'
remove.batches = FALSE

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))


# Aggregate and order the raw + normalized counts:
# ------------------------------------------------
cdf = aggregate(count ~ cls + region, ctdf, sum)
cwide = spread(cdf, region, count)
cmat = as.matrix(cwide [-1])
rownames(cmat) = cwide[,1] 
norm = cmat / rowSums(cmat)

cls.lvls = rev(rownames(reord(norm)))
cdf$cls = factor(cdf$cls, levels=cls.lvls)

tot.cdf = aggregate(count ~ cls, ctdf, sum)
tot.cdf$frac = tot.cdf$count / sum(tot.cdf$count)
tot.cdf$lbl = with(tot.cdf, paste0(count, ' (', round(frac * 100,2),'%)'))
tot.cdf$cls = factor(tot.cdf$cls, levels=cls.lvls)

h = nrow(cmat) / 30 * 8 + 1

# Plot simple total percentage barplot:
# -------------------------------------
g1 = ggplot(cdf, aes(cls, count, fill=region)) + 
    geom_bar(stat='identity', color=NA, position='fill') +
    scale_fill_manual(values=reg.cols) + 
    theme_pubr() + theme(legend.position = 'none') + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    coord_flip()
ggsave(paste0(imgpref, 'proportions_overall_', clsstr, '.png'), g1, dpi=450, units='in', width=5, height=h)
ggsave(paste0(imgpref, 'proportions_overall_', clsstr, '.pdf'), g1, dpi=450, units='in', width=5, height=h)


# Plot simple total counts barplot:
# ---------------------------------
g2 = ggplot(tot.cdf, aes(cls, count, label=lbl, fill=cls)) + 
    geom_bar(stat='identity', color=NA) +
    geom_text() +
    scale_fill_manual(values=tcols) + 
    theme_pubr() + theme(legend.position = 'none') + 
    scale_y_continuous(expand=c(0,0), labels=scales::comma) + 
    coord_flip()
ggsave(paste0(imgpref, 'counts_overall_', clsstr, '.png'), g2, dpi=450, units='in', width=5, height=h)
ggsave(paste0(imgpref, 'counts_overall_', clsstr, '.pdf'), g2, dpi=450, units='in', width=5, height=h)


# Plot joint resources:
# ---------------------
g3 = ggplot(tot.cdf, aes(cls, 1, label=lbl, fill=cls)) + 
    geom_bar(stat='identity', color=NA, position='fill') +
    geom_text() +
    scale_fill_manual(values=tcols) + 
    theme_pubr() + theme(legend.position = 'none') + 
    scale_y_continuous(expand=c(0,0)) + 
    coord_flip()

# Resource:
garr = ggarrange(g3,g1,g2, ncol=3)
ggsave(paste0(imgpref, 'lbls_overall_', clsstr, '.png'), garr, dpi=450, units='in', width=12, height=h)
ggsave(paste0(imgpref, 'lbls_overall_', clsstr, '.pdf'), garr, dpi=450, units='in', width=12, height=h)



# Numbers for overall metadata/paper:
# -----------------------------------
table(cellmeta$major.celltype)
cthrdf = unique(cellmeta[,c('major.celltype','cell_type_high_resolution')])
table(cthrdf$major.celltype)

# Counts per region:
mat = table(cellmeta[,c('major.celltype','region')])

# Percent neuronal:
ctreg = c('AG','MT','PFC')
sum(mat['Exc',ctreg] + mat['Inh',ctreg]) / sum(mat[,ctreg])
(mat['Exc',] + mat['Inh',]) / colSums(mat)

# Percent cycling microglia:
submeta = cellmeta[cellmeta$major.celltype == 'Mic/Immune',]
micmat = table(submeta$cell_type_high_resolution)
micmat / sum(micmat) * 100




