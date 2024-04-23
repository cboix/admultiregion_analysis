#!/usr/bin/R
# -----------------------------------------------
# Create files for the MR UCSC site, minus matrix
# Updated: 09/08/23
# -----------------------------------------------
library(cbrbase)
library(Matrix)
library(gplots)
library(RColorBrewer) # colorRampPalette
library(gplots) # col2hex

set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

options(width=170)

# Directories:
udir = paste0(sdbdir, 'ucsc_datasets/')
outdir = paste0(udir, 'ad-multi-region/')
cmd = paste('mkdir -p', udir, outdir)
system(cmd)


# For deidentification:
# ---------------------
ad430.dbdir = sub('DEVTRAJ', 'AD430', dbdir)
mapdf = read.delim(paste0(ad430.dbdir, 'AD427_subject_projid_mapping.tsv'), header=T)
indmap = as.character(mapdf$subject)
names(indmap) = as.character(mapdf$projid)


# Make the metadata and UMAP tables:
# ----------------------------------
# Subset to matrix barcodes:
bcs = scan(gzfile(paste0(outdir, 'barcodes.tsv.gz')), 'c', quiet=T)
cellmeta = cellmeta[cellmeta$barcode %in% bcs,]

# Remove UMAP coords from metadata, make these separate:
umapdf = cellmeta[,c('barcode', 'U1', 'U2')]
cellmeta$U1 = NULL
cellmeta$U2 = NULL
names(umapdf)[1] = c('cellName')

# Reduce cell metadata to simple annotations:
cellmeta$subject = indmap[as.character(cellmeta$projid)]
keep.cols = c('barcode', 'region', 'subject', 'major.celltype', 'cell_type_high_resolution')
cellmeta = cellmeta[,keep.cols]
names(cellmeta) = c('cellName', 'Region' ,'Individual', 'Major_Cell_Type', 'Sub_Cell_Type')

# Write out files:
write.table(cellmeta, gzfile(paste0(outdir, 'metadata.tsv.gz')), quote=F, row.names=F, sep="\t")
write.table(umapdf, gzfile(paste0(outdir, 'umap.coords.tsv.gz')), quote=F, row.names=F, sep="\t")



# Save the cell colors:
# ---------------------
mct = unique(cellmeta$Major_Cell_Type)
cct = unique(cellmeta$Sub_Cell_Type)
coldf = data.frame(tag=c(mct, cct), col=c(major.col[mct], tcols[cct]))
coldf = rbind(coldf, data.frame(tag=reg.nomb, col=reg.cols[reg.nomb]))
rownames(coldf) = NULL

ind = grep("^#", coldf$col, invert=T)
coldf$col[ind] = col2hex(coldf$col[ind])
coldf$col = sub("#","",coldf$col)
coldf$col = toupper(coldf$col)

write.table(coldf, paste0(outdir, 'multiregion_colors.tsv'), quote=F, col.names=F, row.names=F, sep="\t")


# Make marker genes and quickGenes:
# ---------------------------------
mark.rds = paste0(sdbdir, 'subtype_reg/aggregated_markers_dataframe.Rds')
df = readRDS(mark.rds)

# Subset to protein-coding genes only
pcgenes = anno$symbol[anno$type == 'protein_coding']
df = df[df$gene %in% pcgenes,]
df = df[order(df$est, decreasing=T),]
df = df[order(df$p),]

ntop = 5
ll = lapply(unique(df$celltype), function(x){
    sdf = df[df$celltype == x,]
    return(head(sdf, ntop)) })
markdf = do.call(rbind, ll)
names(markdf)[1] = 'cluster'

markdf = markdf[,c('cluster', 'gene', 'est')]
write.table(markdf, paste0(outdir, 'markers_in.tsv'), quote=F, row.names=F, sep="\t")

# Quick genes:
qg = sapply(unique(markdf$cluster), function(x){
    subdf = markdf[markdf$cluster == x,]
    subdf = subdf[order(subdf$est, decreasing=T),]
    head(subdf$gene, 2)
    })
qg = data.frame(symbol=c(qg))
write.table(qg, paste0(outdir, 'quickGenes.tsv'), quote=F, row.names=F, sep="\t")

