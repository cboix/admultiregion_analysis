#!/usr/bin/R
# -------------------------------------------------
# Collate lists of markers for pseudobulk plotting:
# Updated 04/29/22
# -------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')  # TODO: Change dir for these analyses
cmd = paste('mkdir -p', srdir)
system(cmd)


# Load markers for all cells except Exc/Inh (using Hans' markers):
# ----------------------------------------------------------------
spattern = paste0('^difftl_pseudobulk_markers_.*_.*.rda')
fns = list.files(path=srdir, pattern=spattern)
fns = fns[grep("_Exc_", fns, invert=T)]
fns = fns[grep("_Inh_", fns, invert=T)]

markerdf = NULL
for (fn in fns){
    a = load(paste0(srdir, fn))
    df = est.regdf[est.regdf$var == 'is.stTRUE', c('st','symbol','p','Est')]
    names(df) = c('celltype','gene','p','est')
    markerdf = rbind(markerdf, df)
}

excdf = read.delim(paste0(srdir, 'Table_S1_Exc_markers.tsv'))
inhdf = read.delim(paste0(srdir, 'Table_S1_Inh_markers.tsv'))
excdf = excdf[,c('cluster','gene','p_val','avg_log2FC')]
inhdf = inhdf[,c('cluster','gene','p_val','avg_log2FC')]
names(excdf) = c('celltype','gene','p','est')
names(inhdf) = c('celltype','gene','p','est')

markerdf = rbind(markerdf, excdf)
markerdf = rbind(markerdf, inhdf)

markerdf = markerdf[order(markerdf$est, decreasing=T), ]
markerdf = markerdf[order(markerdf$p), ]
markerdf = markerdf[markerdf$est > 0,]


# Save full markerdf:
# -------------------
fullfile = paste0(srdir, 'aggregated_markers_dataframe.Rds')
saveRDS(markerdf, file=fullfile)


# Save top markers (for faster selection/plotting):
# -------------------------------------------------
# Subset to protein coding genes, remove MT genes:
genes.file = paste0(sdbdir,
    'all_brain_regions_filt_preprocessed_scanpy_norm',
    '_genes.txt')
genes = scan(genes.file, 'c', quiet=T)
markerdf = markerdf[markerdf$gene %in% genes,]
markerdf = markerdf[grep("^MT-", markerdf$gene, invert=TRUE),]

ctlist = lapply(unique(markerdf$celltype), function(x){
    df = markerdf[markerdf$celltype == x,]
    head(df, 10)
    })

topdf = do.call(rbind, ctlist)
topfile = paste0(srdir, 'top_aggregated_markers_dataframe.Rds')
saveRDS(topdf, file=topfile)

