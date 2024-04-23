#!/usr/bin/R
# ---------------------------------------------------------
# Script to load proportion data for plotting + composition
# Updated 11/04/2021 
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(bindir, 'multiRegion/load_metadata.R'))
}

library(tidyr)
library(reshape2)


# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: subset remove.batches")
} else {
    subset = args[1]
    remove.batches= args[2]
}


# Select the appropriate runs for abundance estimation:
# -----------------------------------------------------
# Pick the appropriate cell type resolution column:
if (subset == 'All'){
    clsval = 'major.celltype'
} else if (subset == 'All_minor'){
    clsval = 'minor.celltype'
} else {
    clsval = 'cell_type_high_resolution'
}


# For Excitatory neuron subtypes, pre-subset:
if (subset %in% names(exc.sets)){
    celltype = 'Exc'
    region = exc.set.regions[[subset]]
    subtypes = exc.sets[[subset]]
    region.set = region
} else { 
    region.set = reg.nomb 
}


# Get the unique combinations from metadata:
# ------------------------------------------
cmap = unique(cellmeta[,c('major.celltype',clsval)])
names(cmap)[2] = 'cls'
combdf = expand.grid(cls=unique(cellmeta[[clsval]]),
                     region=region.set, 
                     projid=unique(cellmeta$projid))


# Count up by individual and region in addition to clsval
ctdf = agg.rename(asform(c('barcode ~ projid + region +', clsval)), cellmeta, length, 'count')
names(ctdf)[3] = 'cls'
ctdf = merge(ctdf, cmap)
combdf = merge(combdf, cmap)
ctdf = merge(ctdf, combdf, all.y=TRUE)
ctdf$count[is.na(ctdf$count)] = 0

# Remove problematic batches for effect estimation:
# -------------------------------------------------
# TODO: Only remove if estimating changes:
if (remove.batches){
    # Remove region and individual hugely inflated for fibroblasts:
    ctdf = ctdf[!(ctdf$projid == '50106280' & ctdf$region == 'HC'),]
    if (clsval == 'cell_type_high_resolution'){
        # Remove problematic batch
        ctdf = ctdf[ctdf$cls != 'Exc SV2C LINC02137',] 
    }
}


# Aggregate to total and get number of other cells:
totdf = agg.rename(count ~ projid + region, ctdf, sum, 'total')
ctdf = merge(ctdf, totdf) # Removes a couple missing projid x region comb.
ctdf$other = ctdf$total - ctdf$count


# For subsetting relative to all cells:
if (subset == 'OpcOli'){
    sts = unique(cellmeta[cellmeta$major.celltype %in% c('Opc','Oli'), clsval])
    ctdf = ctdf[ctdf$cls %in% sts,]
    clsstr = paste0(clsval,'_', sub("/","_",subset))
} else if (subset %in% names(exc.sets)){
    ctdf = ctdf[ctdf$cls %in% subtypes,]
    clsstr = paste0(clsval,'_', sub("/","_",subset))
} else if (substr(subset, 1,3) != 'All'){
    sts = unique(cellmeta[cellmeta$major.celltype == subset, clsval])
    ctdf = ctdf[ctdf$cls %in% sts,]
    clsstr = paste0(clsval,'_', sub("/","_",subset))
} else { 
    clsstr = clsval 
}

# Add other attributes:
attr.cols = c('projid','region','rind','braaksc', 'cogdx',
              'niareagansc','msex','age_death','pmi', 'nrad','cogdxad')
cf = unique(metadata[,attr.cols])
ctdf = merge(ctdf, cf)
ctdf = merge(ctdf, pqdf, all.x=TRUE)
ctdf$projid = factor(ctdf$projid)


# For regressions:
ctdf$braaksc.ad = 'CTRL'
ctdf$braaksc.ad[ctdf$braaksc > 4] = 'AD'
ctdf$braaksc.ad = factor(ctdf$braaksc.ad, levels=c('CTRL','AD'))
ctdf$any.nft = ctdf$nft > 0
ctdf$any.plaq_n = ctdf$plaq_n > 0
ctdf$any.plaq_d = ctdf$plaq_d > 0
