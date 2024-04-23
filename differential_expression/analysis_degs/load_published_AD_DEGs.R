#!/usr/bin/R
# ---------------------------------------------
# Load published AD DEGs from previous studies:
# Updated: 07/27/23
# ---------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
# source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
options(width=175)


# List of all studies:
# --------------------
pddir = paste0(sdbdir, 'Published_AD_studies_DEGs/')
deg.meta = read.delim(paste0(pddir, 'metadata.tsv'), header=T)
dirs = list.files(path=pddir)
dirs = dirs[dirs != 'metadata.tsv']

all.degs.df = c()
for (study in dirs){
    stdir = paste0(pddir, study, '/')
    fns = list.files(stdir)
    deg.line = deg.meta[deg.meta$study == study,]
    cols = unlist(deg.line[c('gene','lfc','padj')])

    for (fn in fns){
        df = read.delim(paste0(stdir, fn), sep=",", header=T, check.names=F)
        cat(study, '\t', fn, '\t', nrow(df), '\n')  # TODO: ngenes
        subdf = df[,cols]
        names(subdf) = c('gene', 'lfc', 'padj')
        if (is.na(deg.line$celltype)){
            ct = sub(".csv", "", sub('_.*', '', fn))
            subdf$celltype = ct
        } else {
            subdf$celltype = df[[deg.line$celltype]]
        }
        if (is.na(deg.line$path)){
            subdf$path = 'AD'
        } else {
            subdf$path = df[[deg.line$path]]
        }
        subdf$study = study
        all.degs.df = rbind(all.degs.df, subdf)
    }
}


# Assign significance, check if study only has significant genes reported:
# ------------------------------------------------------------------------
all.degs.df$sig = with(all.degs.df, ifelse(padj < 0.05, ifelse(lfc < 0, 'Down', 'Up'), 'NS'))
mdf = as.data.frame(table(all.degs.df[,c('study', 'sig')]))
mdf = mdf[mdf$sig == 'NS',]
mdf$only.sig = (mdf$Freq == 0)
all.degs.df = merge(all.degs.df, mdf[,c('study', 'only.sig')])

all.degs.df$run = with(all.degs.df, paste0(study, '-', celltype, '@', path))
# table(all.degs.df[,c('run', 'sig')])


# Map all cell types to main classes:
# -----------------------------------
ctdf = unique(all.degs.df[,c('run','celltype')])
ctdf$major.celltype = sapply(ctdf$celltype, function(x){
    substr(tolower(x), 1,3)})
ctdf$major.celltype[ctdf$major.celltype %in% c('ex','ec:','ex0','ex1','exc', 'glu', 'neu')] = 'Excitatory'
ctdf$major.celltype[ctdf$major.celltype %in% c('in', 'inh', 'gab')] = 'Inhibitory'
ctdf$major.celltype[ctdf$major.celltype %in% c('cho','epe','end','per')] = 'Vasc/Epithelia'
ctdf$major.celltype[ctdf$major.celltype %in% c('ast','asc')] = 'Astrocyte'
ctdf$major.celltype[ctdf$major.celltype %in% c('oli','odc')] = 'Oligodendrocyte'
ctdf$major.celltype[ctdf$major.celltype %in% c('opc')] = 'OPC'
ctdf$major.celltype[ctdf$major.celltype %in% c('mg', 'mic', 'cam', 'mg1','mg0','mg4')] = 'Microglia/Immune'

pub.cellmap = ctdf[,c('run','major.celltype')]

