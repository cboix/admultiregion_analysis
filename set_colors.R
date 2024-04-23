#!/usr/bin/R
# ------------------------------------------
# Load metadata for the multiregion project:
# ------------------------------------------
library(cbrbase)
library(RColorBrewer)
library(tidyr)
set_proj('DEVTRAJ')

# Arguments:
args = commandArgs()
if (length(args) > 0){
    project = args[1]
} else { 
    project = 'RNA_Regions'
}

# ---------------
# Metadata files:
# ---------------
mfiles = list.files(path='Annotation', pattern='^metadata_[A-Z]*.tsv')
regions = sub(".tsv", "", sub("metadata_", "",mfiles))
reg.nomb = regions[regions != 'MB']
NREGIONS = length(regions)
reg.order = c('All','AG','MT','PFC','EC','HC','TH')

kept.individuals = scan('Annotation/multiRegion_individuals.txt', 'c')

# --------------------
# Colors for plotting:
# --------------------
# t-SNE colors (preliminary)
snap.cols = scan('Annotation/snap_colors.tsv', 'c')
j = 38
cls.cols = snap.cols[j:(j+15)]
celltype.col = c(cls.cols[c(6, 4, 4, 3, 9, 11, 12, 13)], '#ffa340', cls.cols[2], '#FFE33D', 'grey75')
celltype.col = c(cls.cols[c(6, 4, 4, 3, 9, 11, 12,16)], '#ffa340', cls.cols[2], 'grey80', 'grey75')
names(celltype.col) = c('Astrocyte','Neuronal', 'Excitatory','Inhibitory','Microglia', 'OPC','Pericyte','VLM','Oligodendrocyte','Endothelial', 'Undefined/Dying', 'Other')

major.col = celltype.col[c(1:6,9, 8)]
names(major.col) = c('Ast','Neu','Exc','Inh','Mic/Immune','Opc','Oli', 'Vasc/Epithelia')
tsp.major.col = sapply(major.col, tsp.col)


# Covariate colors:
colvals = list()
colvals[['sex']] = c("female" = "pink", "male" = "lightblue")
colvals[['amyloid.level']] =  c("low" = "royalblue", "high" = "indianred")
colvals[['tangles.level']] =  c("low" = "royalblue", "high" = "indianred")
colvals[['nrad']] =  c("CTRL" = "royalblue", "AD" = "indianred")
colvals[['cogdxad']] =  c("CTRL" = "royalblue", "AD" = "indianred")
colvals[['braaksc.ad']] =  c("CTRL" = "royalblue", "AD" = "indianred")
colvals[['braaksc.early']] =  c("CTRL" = "royalblue", "AD" = "indianred")
# colvals[['apoe_genotype']] =  c("23" = "lightgrey", "33" = "lightgrey", "34" = "darkgrey", "44" = "darkgrey", "0" = "black")
colvals[['apoe_genotype']] =  c("23" = "lightblue", '24' = 'slateblue', "33" = "grey80", "34" = "grey60", "44" = "grey40", "0" = "grey50")
colvals[['sex.and.amyloid']] =  c("female, high amyloid" = "indianred", "female, low amyloid" = "pink", "male, high amyloid" = "royalblue", "male, low amyloid" = "lightblue")
col.paired =  brewer.pal(n=12,name="Paired")
colvals[['niareagansc']] = c('1' = col.paired[6], '2' = col.paired[5], 
                             '3' = col.paired[1], '4' = col.paired[2])
colvals[['cogdx']] = c('1' = col.paired[2], '2' = col.paired[1], 
                       '3' = col.paired[3], '4' = col.paired[5], 
                       '5' = col.paired[6], '6' = 'grey')
colvals[['braaksc']] = c(col.paired[1:6], 'grey80')
names(colvals[['braaksc']]) = as.character(c(1:6,0))

# Colors for regions:
reg.cols = brewer.pal(7, 'Set3')
names(reg.cols) = regions
reg.cols[2] = brewer.pal(12, 'Set3')[12]

# Long names:
reg.long = c('Angular Gyrus', 'Entorhinal Cortex','Hippocampus', 
             'Mammillary Body', 'Mid-Temporal Cortex', 'Prefrontal Cortex','Thalamus')
names(reg.long) = regions

# Colors for individuals:
indvs = sort(kept.individuals)
ind.cols = snap.cols[1:length(indvs)]
names(ind.cols) = as.character(indvs)

# Other useful colors:
colb <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100)
colr <- colorRampPalette(brewer.pal(n=7,name="Reds"))(100)
colrb <- colorRampPalette(brewer.pal(n=7,name="RdBu"))(100)
col.spec <- colorRampPalette(brewer.pal(n=11,name="Spectral"))(100)
load.colors()  # From cbrbase


colr <- colorRampPalette(c('blue','white','red'))(100)


# Cell type colors:
load('Annotation/multiregion_celltypes_colors.Rda')
tsp.tcols = sapply(tcols, tsp.col)
