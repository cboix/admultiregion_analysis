#!/usr/bin/R
# ------------------------------------------------------
# Calculate metadata statistics for the reviews response
# - stats on amyloid / nft burden 
# - stats on staging
# - stats on pathol diagnosis
# - comparison to Montine ABC score
# (all of these on ours vs. full ROSMAP cohort)
# Updated 09/19/2023
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Directories:
plotdir = paste0(imgdir, 'metadata/')
imgpref = paste0(plotdir, 'metadata_comparison_')
cmd = paste('mkdir -p', plotdir)
system(cmd)
options(width=170)


# Load extended metadata:
# -----------------------
indmeta_tsv = 'Annotation/metadata_PFC_all_individuals_092520.tsv'
ext.meta = read.delim(indmeta_tsv, header=T)
ext.meta$kept.ind = ifelse(ext.meta$projid %in% kept.individuals, 'Our Cohort\n(48 individ.)', 'ROSMAP\n(Sept. 2020)')

coh.meta = ext.meta[ext.meta$projid %in% kept.individuals,]

# Braak, niareagan, cerad, etc.
# -----------------------------
int.vars = c('braaksc', 'niareagansc', 'ceradsc', 'cogdx')

df = ext.meta[,c('projid','study', int.vars, 'kept.ind')]
df = gather(df, var, value, -projid, -study, -kept.ind)
df = df[!is.na(df$value),]

gp = ggplot(df, aes(value, fill=kept.ind)) + 
    facet_grid(kept.ind~var, scales='free_y') + 
    geom_bar() + 
    labs(x='Staging or Score', y='Number of individuals') +
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'staging_barplot')
saveGGplot(gp, pltprefix, h=2.75, w=6)


# Quantitative variables across full cohort:
# ------------------------------------------
quant.vars = c('gpath', 'nft', 'amyloid', 'tangles', 'plaq_n', 'plaq_d') # TODO: by region?

df = ext.meta[,c('projid','study', quant.vars, 'kept.ind')]
df = gather(df, var, value, -projid, -study, -kept.ind)
df = df[!is.na(df$value),]

gp = ggplot(df, aes(value, fill=kept.ind)) + 
    facet_grid(kept.ind~var, scales='free') + 
    geom_density() + 
    labs(x='value (density or burden)', y='Number of individuals') +
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() + theme(legend.position='none')

pltprefix = paste0(imgpref, 'overallquant_density')
saveGGplot(gp, pltprefix, h=3, w=8)


# Plaque burden in each region
# ----------------------------


# ABC scores versus NIA-Reagan, cogdx:
# ------------------------------------
abcdf = read.delim('Annotation/abc_scores_092023.tsv', header=T)
abcdf = merge(abcdf, ext.meta[,c('projid', 'niareagansc', 'cogdx')])
abcmat = t(as.matrix(abcdf[,c('a_score', 'b_score', 'c_score', 'nia_aa_sc')]))
scmat = t(as.matrix(abcdf[,c('niareagansc', 'cogdx')]))

pcols = brewer.pal(12, 'Paired')
ux = 1.5
ht = Heatmap(abcmat,
    col=c(pcols[c(2,1,5,6)]),
    cluster_column_slices=FALSE,
    column_split=ifelse(scmat['niareagansc',] < 3, 'NIA-Reagan 1-2', 'NIA-Reagan 3-4'),
    cluster_rows=FALSE,
    border_gp=gpar(color='black', lwd=.5),
    width=ncol(abcmat) * unit(ux, 'mm'),
    height=nrow(abcmat) * unit(ux, 'mm'), 
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
)

pltprefix = paste0(imgpref, 'abc_scores_heatmap')
h = 1 + 1 / 15 * nrow(abcmat)
w = 2 + 1 / 15 * ncol(abcmat)
saveHeatmap(ht, pltprefix, w=w, h=h)



# Plot our samples in context of full cohort


