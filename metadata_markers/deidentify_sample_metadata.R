#!/usr/bin/R
# ------------------------------------------------------
# Deidentify the supplementary metadata for multiregion:
# Updated: 12/13/23
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

# Libraries:
library(tidyr)

options(width=170)


# Read in clinical mapping:
# -------------------------
ad.dbdir = sub("DEVTRAJ", "AD430",dbdir)
mapdf = read.delim(paste0(ad.dbdir, 'AD427_subject_projid_mapping.tsv'), header=T)
submap = mapdf$subject
names(submap) = mapdf$projid


# Deidentify metadata:
meta = metadata[metadata$rind %in% cellmeta$rind,]
meta$subject = submap[as.character(meta$projid)]
meta$projid = NULL
meta$pathAD = ifelse(meta$nrad == 'AD', 'AD', 'non-AD')
meta$clinicalDiag = ifelse(meta$cogdxad == 'AD', 'AD dementia', 'no dementia')

keep.cols = c('library_id', 'subject', 'region', 
    'msex', 'age_death', 'pmi', 'pathAD','clinicalDiag')
meta = meta[,keep.cols]

cuts = c(70, 75, 80, 85, 90)
age.binned = as.character(cut(meta$age_death, breaks=cuts))
age.binned [meta$age_death >= 90] = '90+' 
meta$age_death = age.binned
meta$pmi = round(meta$pmi, 0)

write.table(meta, 'multiregion_sample_metadata_deidentified.tsv', quote=F, sep="\t", row.names=F)

