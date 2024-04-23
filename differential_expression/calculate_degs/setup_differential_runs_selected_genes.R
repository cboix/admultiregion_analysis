#!/usr/bin/R
# -------------------------------------
# Set up a matrix for the differential 
# expression runs using Nebula + RUV
# Only for the selected gene sets
# Updated: 02/23/22
# -------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))


# ----------------------------------------------------
# List all of the cell type contexts that we will use:
# ----------------------------------------------------
df = data.frame(celltype='Mic_Immune', subtype=c('Mic'))

cts = c('Ast','Opc')  # ,'Oli','Inh','Exc')
df = rbind(df, data.frame(celltype=cts, subtype=cts))


# -----------------------------------------------
# Merge the regions + AD ascertainment variables:
# -----------------------------------------------
pathlist = c('nrad','nft','plaq_n','plaq_d','cogdxad', 'braaksc.early','braaksc.ad')
reglist = c('allregions', reg.nomb)
pdf = expand.grid(region=reglist, path=pathlist)

# Remove TH from nft/plaq_n/plaq_d:
pdf = pdf[!(pdf$path %in% c('nft','plaq_n','plaq_d') &
            pdf$region == 'TH'),]

# Merge all:
testdf = merge(df, pdf)

testdf$geneset = 'metabolic_genes'

# Write out:
write.table(testdf, paste0(datadir, 'nebula_wRUV_selgenes_runlist.tsv'), 
            quote=F, row.names=F, sep="\t")

