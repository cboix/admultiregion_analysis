#!/usr/bin/R
# ------------------------------------------------------------------------
# Load the cross-module scores across all celltypes (at pseudobulk level):
# Updated 11/30/2021
# ------------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){ source(paste0(sbindir, 'load_metadata.R')) }


crossdir = paste0(sdbdir, 'crossmodule/')
cmd = paste('mkdir -p', crossdir)
system(cmd)


# Set the run arguments:
# ----------------------
# Default arguments:
graph_id = 'boot'
runlist = c('Ast','Oli','Opc','Mic_Immune','Vasc_Epithelia',
            'Inh','Exc', 'All') 
# 'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons')
fullscores.rda = paste0(crossdir, 'aggregated_crossct_psbulk_scores.rda')
fullscores.file = paste0(crossdir, 'aggregated_crossct_psbulk_scores.tsv.gz')

if (!file.exists(fullscores.rda)){
    gmlist = list()
    cmlist = list()
    scoredf = NULL
    for (runset in runlist){
        print(paste("Loading", runset))
        # Load in runset data and prefixes:
        # ---------------------------------
        commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
        source(paste0(sbindir, 'modules/load_modules_degenr.R'))


        # Load in the full pseudobulk data for these subtypes:
        # NOTE: From modules_degenr_03_plot_pseudobulk_modules.R
        # ------------------------------------------------------
        useset = 'coregenes_'
        runscores.file = paste0(moddir, 'module_pseudobulk_scores_', 
                                useset, fullpref, '.tsv.gz')
        scdf = read.delim(gzfile(runscores.file), sep="\t")
        # scdf = merge(scdf, mmap, all.x=TRUE)
        scdf = merge(scdf, mmap)
        scdf$pr = with(scdf, paste0(projid, '-', region))
        scdf$cm = with(scdf, paste0(cell_type_high_resolution, '-', module))
        scdf$runset = runset

        # Get number of genes:
        ngenes = table(nodedf$leiden)
        scdf$ng = ngenes[as.character(scdf$module)]

        # Concatenate all:
        scoredf = rbind(scoredf, scdf)

        # Save the genelists:
        gmlist[[runset]] = genemap
        cmlist[[runset]] = coremap
    }

    save(scoredf, gmlist, cmlist, file=fullscores.rda)
    write.table(scoredf, gzfile(fullscores.file), quote=F, row.names=F, sep="\t")
} else { load(fullscores.rda) }


# Remove flagged modules for ambient RNA:
# ---------------------------------------
# scoredf = scoredf[!((scoredf$module %in% c(4)) & (scoredf$runset == 'Ast')),]
# scoredf = scoredf[!((scoredf$module %in% c(2,12)) & (scoredf$runset == 'Mic_Immune')),]
# scoredf = scoredf[!((scoredf$module %in% c(1,11)) & (scoredf$runset == 'Opc')),]
# scoredf = scoredf[!((scoredf$module %in% c(4)) & (scoredf$runset == 'Oli')),]

col_fun = colorRamp2(c(-1,0,1), c("blue", "white", "red"))

# Aggregate score matrix for runset level:
# ----------------------------------------
scoredf$totscore = scoredf$score * scoredf$ncell

runscdf = aggregate(cbind(totscore, ncell) ~ runset + module + pr + ng, scoredf, sum)
runscdf$score = runscdf$totscore / runscdf$ncell
runscdf$rm = with(runscdf, paste0(runset, '-', module))


# Module mappings:
# ----------------
rmap = unique(scoredf[,c('module','mname','runset')])
rmap$rname = with(rmap, paste0(runset, '-', module))
rownames(rmap) = rmap$rname

smap = unique(scoredf[,c('module','mname','cell_type_high_resolution')])
smap$rname = with(smap, paste0(cell_type_high_resolution, '-', module))
rownames(smap) = smap$rname


# For making the graph (with smap, rmap):
# ---------------------------------------
prune.edges = function(edgedf, modlevel, nodes=NULL){
    if (modlevel == 'subtype'){
        edgedf = edgedf[grep("MKI67",edgedf$M1, invert=TRUE),]
        edgedf = edgedf[grep("MKI67",edgedf$M2, invert=TRUE),]
    }
    if (is.null(nodes)){
        nodes = sort(unique(c(edgedf$M1, edgedf$M2)))
    }
    nodemap = 1:length(nodes)
    names(nodemap) = nodes
    edgedf$i = nodemap[edgedf$M1]
    edgedf$j = nodemap[edgedf$M2]
    mndf = data.frame(ind=1:length(nodes), node=nodes)

    if (modlevel == 'runset'){
        mndf$mname = rmap[nodes,'mname']
        mndf$runset = rmap[nodes,'runset']
        mndf$name = paste0(mndf$runset, "_", sapply(mndf$mname, function(x){sub(" [(].*[)]", "", x)}))
        mcol = major.col[sapply(mndf$runset, function(x){gsub("_", "/", x)})]
        mndf$col = ifelse(mndf$runset %in% names(exc.sets), major.col['Exc'], mcol)
    } else {
        mndf$mname = smap[nodes,'mname']
        mndf$celltype = smap[nodes,'cell_type_high_resolution']
        mndf$name = paste0(mndf$celltype, "_", sapply(mndf$mname, function(x){sub(" [(].*[)]", "", x)}))
        mndf$col = tcols[mndf$celltype]
    }
    return(list(edgedf=edgedf, mndf=mndf))
}


