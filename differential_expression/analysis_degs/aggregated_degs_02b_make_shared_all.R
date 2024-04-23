#!/usr/bin/R
# -----------------------------------------------------------
# Make the set of shared DEGs "All_All" for every combination
# Updated: 03/27/22
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)


# Read in the full RDS results across all cell types:
# ---------------------------------------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc",
    "Oli_Oli", 'Inh_Inh','Exc_Exc')
# pathlist = unique(totnsigdf$path)

# Requires running:
# aggregated_degs_04_aggregate_regions.R
fulldf = c()
for (set in keep.sets){
    print(set)
    full.rds = paste0(regdir, 'aggregated_fullset.', set, '.rds')
    setdf = readRDS(full.rds)
    setdf$set = set
    fulldf = rbind(fulldf, setdf)
}


# Reduce to a merged data.frame (for "All_All")
# ---------------------------------------------
cdf = agg.rename(coef_mast ~ gene + col_nm + path + region, 
    fulldf, length, 'count')
cdf$count = cdf$count + .1 * cdf$col_nm  # Break ties
cdf = merge(cdf, aggregate(count ~ gene + path + region, cdf, max))
cdf = cdf[cdf$count >= 3, ]
cdf$count = round(cdf$count)
head(cdf[order(cdf$count, decreasing=T),], 50)


# Save this merged table:
# -----------------------
desuff = 'All_All_allconditions'
fname = paste0(regdir, 'allmethods.', desuff, '.merged.tsv.gz')
fname.rds = paste0(regdir, 'allmethods.', desuff, '.merged.rds')

write.table(cdf, gzfile(fname), quote=F, sep="\t", row.names=F)
saveRDS(cdf, file=fname.rds)

