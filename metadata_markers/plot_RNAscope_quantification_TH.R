#!/usr/bin/R
# ---------------------------------------------------
# Plot the RNAscope quantification for the thalamus: 
# Updated 12/07/2021 
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(viridis)
library(pzfx)  # For reading prism data
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)
library(patchwork)


# Directories:
fracdir = paste0(sdbdir, 'fractions/')
plotdir = paste0(imgdir, 'fractions/')
imgpref = paste0(plotdir, 'RNAscope_')
cmd = paste('mkdir -p', plotdir, fracdir)
system(cmd)


# List and read in tables from prism file:
# ----------------------------------------
pzfile = paste0(fracdir, 'FOXP2_MEIS2_multiregion.pzfx')
tables = pzfx_tables(pzfile)
pzfx_tables(pzfile)

df1 = read_pzfx(pzfile, table=1)
df2 = read_pzfx(pzfile, table=2)
names(df1) = c('TH','PFC')
names(df2) = c('TH','PFC')

df1$table = tables[1]
df2$table = tables[2]

df = rbind(df1, df2)


# Reshape tables, remove NAs (ragged table), and plot:
# ----------------------------------------------------
df = gather(df, region, transcripts, -table)
df = df[!(is.na(df$transcripts)),]


gp = ggplot(df, aes(region, transcripts, fill=region)) + 
    facet_wrap(~table) + 
    scale_fill_manual(values=c(reg.cols[c('TH','PFC')])) + 
    geom_boxplot(outlier.shape=NA) + 
    # geom_violin(scale='width') + 
    geom_jitter(width=.15, cex=1, height=0) + 
    stat_compare_means() + 
    scale_y_continuous(expand=c(0,0)) + 
    labs(x='Region', y='# Transcripts per GAD2+ cell') + 
    theme_pubr()

ggsave(paste0(imgpref, 'TH_quantification_boxplots.png'), gp, dpi=450, units='in', width=2.25, height=3.5)
ggsave(paste0(imgpref, 'TH_quantification_boxplots.pdf'), gp, dpi=450, units='in', width=2.25, height=3.5)





