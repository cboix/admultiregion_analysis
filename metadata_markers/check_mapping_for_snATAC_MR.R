#!/usr/bin/R
# -------------------------------------------------------------
# Check the barcode / sample mapping for the snATAC-MR project:
# Updated 05/12/2023
# -------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
options(width=170)

# Directories:
plotdir = paste0(imgdir, 'metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir)
system(cmd)


# Read in the integration:
# ------------------------
groupdir = '/net/bmc-lab5/data/kellis/group/'
zpdir = 'Zunpeng_Sharing/MR_RNA/Integration/'
metafile = paste0(groupdir, zpdir, 'MR_snRNA_integrated_1_7M_cells_metadata.tsv.gz')
mdf = read.delim(gzfile(metafile), header=T)
mdf = mdf[,-1]

# Process mdf data:
regmap = c('AG', 'TH', 'EC', 'HC', 'MFC', 'MT', 'PFC')
names(regmap) = c('Angular_gyrus', 'Anterior_thalamus', 'Entorhinal_cortex', 
    'Hippocampus', 'Medialfrontal_cortex', 'Midtemporal_cortex', 'Prefrontal_cortex')
mdf$region = regmap[mdf$batch]

# Merge cellmeta with library ID, check mdf:
cellmeta = merge(cellmeta, metadata[,c('library_id', 'rind')])
mdf$in.meta = mdf$Sample %in% unique(cellmeta$library_id)
table(mdf[,c('region', 'in.meta')])

# Extract cell barcodes in each:
mdf$CB = sub("-.*","", sub(".*#", "", mdf$Sample_barcode))
cellmeta$CB = sub("-.*","", sub(".*_", "", cellmeta$barcode))

# Save for sharing:
filepref = paste0(groupdir, 'cboix/AD_multiRegion/', 'cell_metadata_with_library_id_051223')
saveRDS(cellmeta, file=paste0(filepref, '.Rds'))
write.table(cellmeta, file=gzfile(paste0(filepref, '.tsv.gz')), quote=F, sep="\t", row.names=F)


# Look at quality of mapping alignment:
# -------------------------------------
lids = unique(cellmeta$library_id)
err.lids = c()
statdf = c()
tot.accounted = 0
for (lid in lids){
    reg = cellmeta$region[cellmeta$library_id == lid][1]
    ccb = cellmeta$CB[cellmeta$library_id == lid]
    mcb = mdf$CB[mdf$Sample == lid]
    int = length(intersect(mcb, ccb))
    statdf = rbind(statdf, 
        data.frame(library_id=lid, region=reg, 
            cb.pct=int / length(ccb), int.pct=int / length(mcb)))
    cat(lid, reg, int / length(ccb), int / length(mcb), '\n')
    if (int / length(ccb) < .5){ err.lids = c(err.lids, lid) }
    tot.accounted = tot.accounted + int
}
print(err.lids) # 4 in my data not in integration
print(tot.accounted) # 1.322M accounted for, of 1.353M

statdf = statdf[order(statdf$cb.pct),]


# For error library IDs, check barcodes against all of integration barcodes:
# --------------------------------------------------------------------------
# Match barcodes to cellmeta from MR project:
errdf = c()
for (lid in err.lids){
    reg = cellmeta$region[cellmeta$library_id == lid][1]
    ccb = cellmeta$CB[cellmeta$library_id == lid]
    mm = mdf$CB %in% ccb
    tab = sort(table(mdf$Sample[mm]), decreasing=T)
    errdf = rbind(errdf, 
        data.frame(library_id=lid, orig=length(ccb), sum=sum(mm), max=tab[1], 
            match=names(tab[1]), match.pct=tab[1] / sum(mm)))
}
# Conclusion: these four IDs have been excluded from the MR integration
          # library_id orig sum max     match  match.pct
# D19-12350   D17-9566 1724  81   3 D19-12350 0.03703704
# D19-12338   D17-9567 1703  59   4 D19-12338 0.06779661
# D19-8375    D17-9572 1413  33   2  D19-8375 0.06060606
# D19-8367   D19-12336 7665 229   5  D19-8367 0.02183406




