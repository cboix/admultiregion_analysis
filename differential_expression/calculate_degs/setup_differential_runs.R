#!/usr/bin/R
# -------------------------------------
# Set up a matrix for the differential 
# expression runs using Nebula + RUV
# Updated: 11/24/21
# -------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))


# ----------------------------------------------------
# List all of the cell type contexts that we will use:
# ----------------------------------------------------
df = data.frame(celltype='Mic_Immune',subtype=c('Mic','CAM','T'))

vasc.sts = unique(cellmeta$cell_type_high_resolution[cellmeta$major.celltype == 'Vasc/Epithelia'])
df = rbind(df, data.frame(celltype='Vasc_Epithelia', subtype=vasc.sts))

cts = c('Ast','Opc','Oli','Inh','Exc')
df = rbind(df, data.frame(celltype=cts, subtype=cts))

# -----------------------------------------------
# Merge the regions + AD ascertainment variables:
# -----------------------------------------------
pathlist = c('nrad','nft','plaq_n','plaq_d','cogdxad')
reglist = c('allregions', reg.nomb)
pdf = expand.grid(region=reglist, path=pathlist)

# Remove TH from nft/plaq_n/plaq_d:
pdf = pdf[!(pdf$path %in% c('nft','plaq_n','plaq_d') &
            pdf$region == 'TH'),]

# Merge all:
testdf = merge(df, pdf)

# Write out:
write.table(testdf, paste0(datadir, 'nebula_wRUV_runlist.tsv'), 
            quote=F, row.names=F, sep="\t")


# ---------------------------------------------------------
# Setup the differential runs for specific neuron subtypes:
# ---------------------------------------------------------
# Select neuronal subtypes in EC, HC, and TH:
rdf = agg.rename(barcode ~ region + cell_type_high_resolution, cellmeta[cellmeta$major.celltype == 'Exc',], length, 'count')
rdf = spread(rdf, region, count, fill=0)
rmat = as.matrix(rdf[,-1])
rownames(rmat) = rdf[,1]

df = c()
kept = c()
for (reg in c('EC','HC','TH')){
    top = apply(rmat,1, max) == rmat[,reg]
    top = names(top)[top]
    top = top[top!= 'Exc SV2C LINC02137']
    kept = c(kept, top)
    top = gsub(" ", "_", top)
    df = rbind(df, data.frame(celltype='Exc', subtype=top, region=reg))
}

# Evaluate the others in the neocortex only:
other = rownames(rmat)[!(rownames(rmat) %in% kept)]
other = gsub(" ", "_", other)
df = rbind(df, data.frame(celltype='Exc', subtype=other, region='neocortex'))

# Merge all:
pathlist = unique(c(pathlist, 'braaksc.early', 'braaksc.ad'))
testdf = merge(df, data.frame(path=pathlist))
testdf = testdf[order(testdf$region),]

# Remove TH from nft/plaq_n/plaq_d:
testdf = testdf[!(testdf$path %in% c('nft','plaq_n','plaq_d') &
                  testdf$region == 'TH'),]

# Write out:
write.table(testdf, 
            paste0(datadir, 'nebula_wRUV_excitatory_subsets_runlist.tsv'), 
            quote=F, row.names=F, sep="\t")

