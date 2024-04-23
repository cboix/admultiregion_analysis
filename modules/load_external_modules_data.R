#!/usr/bin/R
# --------------------------------------------------------
# Load data for the modules computed on external datasets:
# Updated 12/01/2023
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')

library(tidyr)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: dataset graph_id load.repr load.scores")
} else {
    dataset = args[1]
    graph_id = args[2]
    if (length(args) > 2){
        load.scores = args[3]
    } else { load.scores = FALSE }
}

# Directories:
extdir = paste0(sdbdir, 'external_datasets/')
dsdir = paste0(extdir, dataset, '/')

# Prefixes:
modpref = paste0(dataset, '_z4.5')
graphpref = paste0(modpref, '_', graph_id)


# Load the mapping between genes and modules:
# -------------------------------------------
assigndf = read.delim(paste0(dsdir, 'allgene_assignments_', graphpref, '.tsv'), header=T)
genemap = assigndf$module
names(genemap) = assigndf$gene


# Load module names:
# ------------------
namefile = paste0(dsdir, 'module_names_', graphpref, '.tsv')
mmap = read.delim(namefile, header=F)
names(mmap) = 'mname'
mmap$module = (1:nrow(mmap) - 1)


# Load representation data and core genes:
# ----------------------------------------
nodefile = paste0(dsdir, 'graph_nodes_info_', graphpref, '.tsv')
nodedf = read.delim(nodefile, sep="\t")
coremap = nodedf$leiden
names(coremap) = nodedf$gene

edgefile = paste0(dsdir, 'graph_edgelist_', graphpref, '.tsv')
edgedf = read.delim(edgefile, sep=" ", header=F)


# Load the scores as well (can be quite large):
# ---------------------------------------------
if (load.scores){
    scorefile = paste0(dsdir, 'module_cell_scores_', graphpref, '.tsv.gz')
    scoredf = read.delim(gzfile(scorefile), sep="\t", header=T)[,-1]
    scoremat = as.matrix(scoredf[,-1])
    rownames(scoremat) = scoredf$bc
    rm(scoredf)
}
