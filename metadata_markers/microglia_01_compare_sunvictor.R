#!/usr/bin/R
# ------------------------------------------------------
# Compare our microglia definitions to Sun/Victor et al.
# Updated 11/17/2023
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(viridis)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)

library(ComplexHeatmap)
library(circlize)
print(version)
options(width=170)

# Directories:
plotdir = paste0(imgdir, 'metadata/')
imgpref = paste0(plotdir, 'microglia_review_')
cmd = paste('mkdir -p', plotdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Load the Sun/Victor metadata:
# -----------------------------
ad.dbdir = sub("DEVTRAJ","AD430", dbdir)
micdir = paste0(ad.dbdir, 'snPFC/ucsc_datasets/')
micdf = readRDS(paste0(micdir, "ROSMAP.Microglia.6regions.metadata.rds"))

coldf = read.delim(paste0(micdir, 'microglia-states/microglia_states_colors.tsv'), header=F)
mic.cols = paste0("#", coldf$V2)
names(mic.cols) = coldf$V1

# Subset to microglia in multi-region dataset:
micdf = micdf[grep("^SM_", micdf$batch, invert=T),]
submeta = cellmeta[grep("^Mic", cellmeta$cell_type_high_resolution),]

# Map region names:
brmap = c('AngularGyrus'='AG', 'EntorhinalCortex'='EC',
    'Hippocampus'='HC', 'MidtemporalCortex'='MT',
    'PFC'='PFC', 'Thalamus'='TH')
micdf$region = brmap[micdf$brainRegion]

# Map the dataset definitions to each other:
micdf$bc = sub("-[0-9]+$", "", micdf$barcode)
micdf$dataset = paste0(micdf$region, '-', sub("^[A-Z]+-", "", micdf$barcode))
submeta$bc = sub("-[0-9]+$", "", sub("^[A-Z]+_", "", submeta$barcode))

# Rind vs. dataset:
micdf$rind = ''
for (region in reg.nomb){
    print(region)
    reg.rind = unique(submeta$rind[submeta$region == region])
    reg.dataset = unique(micdf$dataset[micdf$region == region])

    jmat = matrix(0, nrow=length(reg.rind), ncol=length(reg.dataset), dimnames=list(reg.rind, reg.dataset))
    for (rind in reg.rind){
        mr.bcs = submeta$bc[submeta$rind == rind]
        for (dset in reg.dataset){
            mic.bcs = micdf$bc[micdf$dataset == dset]
            int = length(intersect(mr.bcs, mic.bcs))
            union= length(unique(c(mr.bcs, mic.bcs)))
            jmat[rind,dset] = int / union
        }
    }
    dset.map = apply(jmat, 2, which.max)
    dset.map[] = rownames(jmat)[dset.map]
    reg.ind = micdf$region == region
    micdf$rind[reg.ind] = dset.map[micdf$dataset[reg.ind]]
}


# Update microglia barcodes:
# --------------------------
micdf$barcode = paste0(micdf$region, '_', micdf$bc, '-', sub("^[A-Z]+\\.", "", micdf$rind))
mean(micdf$barcode %in% submeta$barcode) # 89% recovery

# Merge datasets:
df = merge(micdf[,c('barcode','rind','cellstate','region')],
    submeta[,c('barcode','rind','cell_type_high_resolution','region', 'U1','U2')])


# Alignment between subtypes:
# ---------------------------
mapdf = as.data.frame(table(df[,c('cellstate','cell_type_high_resolution')]))

# Barplot (heatmap isn't necessary)
gp = ggplot(mapdf, aes(cell_type_high_resolution, Freq, fill=cellstate)) + 
    geom_bar(stat='identity', position='fill') +
    scale_fill_manual(values=mic.cols) + 
    scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
    labs(y='% of cells', x='Microglia subtype') +
    coord_flip() + theme_pubr() 
pltprefix = paste0(imgpref, "pctbreakdown_micstates_barh")
saveGGplot(gp, pltprefix, w=6, h=3.5)


