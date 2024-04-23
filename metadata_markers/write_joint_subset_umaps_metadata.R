#!/usr/bin/R
# ----------------------------------------
# Write the merged + updated subtype UMAPs
# Updated 01/07/2022
# ----------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


# Create the C1/C2 (subtype UMAP) columns:
# ----------------------------------------
cellmeta$C1 = cellmeta$U1
cellmeta$C2 = cellmeta$U2
rownames(cellmeta) = cellmeta$barcode

for (runset in c('Exc','Inh', 'Vasc_Epithelia')){
    ct = sub("_", "/", runset)
    if (runset %in% c('Exc', 'Inh')){
        cellsuff = paste0(runset, '_log1p_combat_filthvg1000_nomt')
        ctumap.file = paste0(sdbdir, 'metadata/multiregion_majorctUMAP_', cellsuff, '.tsv')
    } else {
        # Mainly for vasculature, which othw can't be subsetted
        ctumap.file = paste0(sdbdir, 'multiregion_majorctUMAP_', runset, '_combat.tsv')
    } 
    # Match metadata to UMAP coordinates:
    # -----------------------------------
    if (runset %in% c('Exc', 'Inh', names(exc.sets))){
        submeta = read.delim(ctumap.file, sep="\t")
        submeta = submeta[, c('barcode','C1','C2')]
        rownames(cellmeta) = cellmeta$barcode
    } else if (runset == 'Vasc_Epithelia'){
        submeta = read.delim(ctumap.file)
        submeta = submeta[,c('barcode','C1','C2')]
    }
    cellmeta[submeta$barcode, 'C1'] = submeta$C1
    cellmeta[submeta$barcode, 'C2'] = submeta$C2
}

# Save updated:
cellmeta.rds = 'Annotation/multiregion_cell_metadata.Rds'
if (!file.exists(cellmeta.rds)){
    sel.cols = c('barcode','rind','region','projid','U1','U2','C1','C2',
        # TODO: Do we need all of these?
                 'major.celltype','minor.celltype','neuronal.layer','inh.subtype',
                 'neuronal.exttype','full.exttype','cell_type_high_resolution')
    saveRDS(cellmeta[,sel.cols], file=cellmeta.rds)
}


