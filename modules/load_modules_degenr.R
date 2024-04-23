#!/usr/bin/R
# ---------------------------------------------------------
# Replot the DEG enrichment table from the modules package:
# Updated 11/27/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(bindir, 'multiRegion/load_metadata.R'))
}

library(tidyr)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: runset graph_id load.repr load.scores")
} else {
    runset = args[1]
    graph_id = args[2]
    if (length(args) > 2){
        load.repr = args[3]  # Load representation or not
    } else { load.repr = FALSE }
    if (length(args) > 3){
        load.scores = args[4]
    } else { load.scores = FALSE }
}


# Centralized options for specific run sets:
# ------------------------------------------
if (runset %in% names(exc.sets)){
    celltype = 'Exc'
    region = exc.set.regions[[runset]]
    subtypes = exc.sets[[runset]]
    plt.subtypes = subtypes
    if (runset == 'ECneurons'){
        subtype = 'Exc TOX3 TTC6'
    } else {
        subtype = subtypes[1]
    }
    modsuff = paste0('_', runset, '_', paste0(region, collapse="_"))
    region.set = region
} else {
    region = 'allregions'
    region.set = reg.nomb
    modsuff = ''
    celltype = runset
    if (celltype == 'Mic_Immune'){
        celltype = 'Mic/Immune'
        subtype = 'Mic'
    } else if (celltype == 'Vasc_Epithelia'){
        celltype = 'Vasc/Epithelia'
        subtype = 'End'
    } else {
        subtype = celltype
    }
    subtypes = unique(cellmeta[cellmeta$major.celltype == celltype, 'cell_type_high_resolution'])
    plt.subtypes = subtypes
    if (celltype == 'All'){
        plt.subtypes = unique(cellmeta$cell_type_high_resolution)
    }
    if (celltype %in% c('Opc','Oli')) { 
        subtypes = celltype
    } else if (celltype == 'Mic/Immune'){
        subtypes = c('Mic','CAM','T')
    }
}

# Directories:
moddir = paste0(sdbdir, 'modules/')


# Load the modules x DEG enrichment table for a cell type:
# --------------------------------------------------------
ctstr = gsub("/","_", celltype)
ststr = gsub("/","_", gsub(" ","_", subtype))

# TODO: Set this more robustly (e.g. Mic_Immune and Vasc_Epithelia)
if (runset == 'CTXneurons'){
    depref = paste0(ctstr, '_', ststr, '_neocortex_')
} else {
    depref = paste0(ctstr, '_', ststr, '_', region, '_')
}
modpref = paste0(ctstr, modsuff, '_log1p_z4.5')
fullpref = paste0(depref, modpref, '_', graph_id)
graphpref = paste0(modpref, "_", graph_id)


# Load the mapping between genes and modules:
# -------------------------------------------
assigndf = read.delim(paste0(moddir, 'allgene_assignments_', graphpref, '.tsv'), header=T)
genemap = assigndf$module
names(genemap) = assigndf$gene


# Load the statistics for the graph/modules/DEGs combination:
# -----------------------------------------------------------
stats.file = paste0(moddir, 'enrstats_allpath_allmethods.', fullpref, '.tsv')
if (file.exists(stats.file)){
    statsdf = read.delim(stats.file, header=T)
    # Adjust p-value + add joint key:
    statsdf$p.adj = p.adjust(statsdf$p, 'fdr')
    statsdf$jointkey = with(statsdf, paste0(key, '.', dkey))
    # Mapping from module to short name:
    mmap = unique(statsdf[,c('module','mname')])
    mmap = mmap[order(mmap$module), ]
}

# Load in the aggregated differential genes used for this enrichment:
de.file = paste0(moddir, 'allpath_allmethods.', depref, 'usedDEGs.tsv')
if (file.exists(de.file)){
    dedf = read.delim(de.file, header=T)
    names(dedf)[1:3] = c('index1', 'dkey', 'index2')
    dedf$module = genemap[dedf$gene]
}


# For celltypes where we use a different UMAP:
# --------------------------------------------
if (runset %in% c('Exc', 'Inh', names(exc.sets))){
    cellsuff = paste0(celltype, '_log1p_combat_filthvg1000_nomt')
    ctumap.file = paste0(sdbdir, 'metadata/multiregion_majorctUMAP_', cellsuff, '.tsv')
} else {
    # Mainly for vasculature, which othw can't be subsetted
    ctumap.file = paste0(sdbdir, 'multiregion_majorctUMAP_', runset, '_combat.tsv')
} 

# mod.df.file = paste0(moddir, 'module_df_', modpref, '_leiden_', graph_id, '.tsv')
# moddf = read.delim(mod.df.file, sep="\t")


# Load representation data if needed:
# -----------------------------------
if (load.repr){
    nodefile = paste0(moddir, 'graph_nodes_info_', graphpref, '.tsv')
    nodedf = read.delim(nodefile, sep="\t")
    coremap = nodedf$leiden
    names(coremap) = nodedf$gene

    edgefile = paste0(moddir, 'graph_edgelist_', graphpref, '.tsv')
    edgedf = read.delim(edgefile, sep=" ", header=F)
} 


# Load the scores as well (can be quite large):
# ---------------------------------------------
if (load.scores){
    scorefile = paste0(moddir, 'module_cell_scores_', graphpref, '.tsv.gz')
    scoredf = read.delim(gzfile(scorefile), sep="\t", header=T)[,-1]
    scoremat = as.matrix(scoredf[,-1])
    rownames(scoremat) = scoredf$bc
    rm(scoredf)
}
