#!/usr/bin/R
# ---------------------------------------------------------
# Replot the DEG enrichment table from the modules package:
# Updated 11/27/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(gprofiler2)

library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Functions to process and plot enrichments:
# ------------------------------------------
# TODO: Improve picking top.ids so that we don't get duplicates:
top.ids = function(x, ntop=10, cutoff=0.05){ 
    ord = order(x, decreasing=F)
    ids = head(ord, ntop)
    ids = ids[x[ids] < cutoff] 
    return(ids) }

gpPvalMatrix = function(gpdf, genesets, ntop=8){
    # Get the enrichment p-values matrix:
    pmat = t(as.matrix(data.frame(gpdf$p_value)))
    rownames(pmat) = NULL
    colnames(pmat) = names(genesets)
    cs = colSums(pmat < 0.05)
    pmat = pmat[, cs > 0]

    # Subset to the top n terms for each:
    ids = apply(pmat, 2, ntop=ntop, top.ids)
    ids = unique(c(unlist(ids)))
    subpmat = -log10(pmat[ids,])
    subpmat[subpmat > 10] = 10
    rownames(subpmat) = gpdf[ids,'term_name']
    return(subpmat)
}

plotGpPvalMatrix = function(subpmat, pltprefix){
    # pltmat = t(reord(t(subpmat)))
    pltmat = t(diag.mat2(t(subpmat))[[1]])
    pltmat = pltmat[rev(rownames(pltmat)),]
    rowsplit = colnames(pltmat)[apply(pltmat, 1, which.max)]
    rowsplit = sapply(rowsplit, function(x){ sub("^M","",sub(" .*","",x)) })
    rowsplit = as.numeric(rowsplit)
    plt = Heatmap(pltmat,
                  cluster_columns=FALSE,
                  cluster_rows=FALSE,
                  row_split=rowsplit,
                  use_raster=TRUE,
                  border_gp=gpar(color='black'),
                  name='-log10(p)',
                  col=c('white',colb))

    h = 2.25 + 2.5 / 15 * nrow(pltmat)
    w = 5 + 2.5 / 15 * ncol(pltmat)
    pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
    print(plt)
    dev.off()
    png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
    print(plt)
    dev.off()
    return(plt)
}


# Identify/print the top associations:
# ------------------------------------
subdf = statsdf[statsdf$p.adj < 0.05, c('mname','key','dkey','p.adj','ratio', 'module')]
subdf = subdf[order(subdf$p.adj),]
print(tibble(subdf)[1:10,])  # Just to see top associations:

# Print the top DE genes driving each enrichment:
for (i in 1:nrow(subdf)){
    m = subdf$module[i]
    key = subdf$key[i]
    path = subdf$dkey[i]

    cat('\n')
    cat(as.character(subdf[i,]), '\n')
    sdedf = dedf[(dedf$module == m) & (dedf$gset == key) & (dedf$dkey == path),]
    cat(head(sdedf$gene, 25),'\n')
    sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
    gp2.result = gprofiler2::gost(sdedf$gene, organism='hsapiens',
                                  ordered_query=FALSE, multi_query=FALSE,
                                  sources = sources)
    gpdf = gp2.result$result
    if (!is.null(gpdf)){
        gpdf = gpdf[order(gpdf$p_value),]
        sgpdf = gpdf[gpdf$source %in% c('KEGG','REAC','WP','CORUM'),]
        cat(head(sgpdf$term_name,10), '\n')
        cat(head(gpdf[gpdf$term_size < 1500,'term_name'],10),'\n')
        # print(head(sgpdf$term_name,10))
        # print(head(gpdf[gpdf$term_size < 1500,'term_name'],10))
    }
}


# Run GO enrichments on all modules jointly, for plotting:
# run both on DE-tested genes vs. on allgenes:
# --------------------------------------------------------
for (set in c('deonly', 'allgenes', 'coregenes')){
    print(set)
    moduleenr.file = paste0(moddir, 'module_enrichments_',set,'_', fullpref, '.rda')
    # if (!file.exists(moduleenr.file)){
        topmod = mmap$module
        if (set == 'deonly'){
            genesets = lapply(topmod, function(m){
                                  x = unique(dedf[(dedf$module == m),'gene'])
                                  x[!(is.na(x))]
                                  })
        } else if (set == 'coregenes'){
            genesets = lapply(topmod, function(m){ 
                                  x = names(coremap[(coremap == m)]) 
                                  x[!(is.na(x))]
                                  })
        } else {
            genesets = lapply(topmod, function(m){ 
                                  x = names(genemap[(genemap == m)]) 
                                  x[!(is.na(x))]
                                  })
        }
        names(genesets) = mmap[topmod + 1, 'mname']

        # Remove empty sets:
        dl = c(lapply(genesets, length))
        keep = names(dl)[dl > 0]
        genesets = genesets[keep]

        sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
        gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
                                      ordered_query=FALSE, multi_query=TRUE,
                                      sources = sources)
        save(gp2.result, genesets, file=moduleenr.file)

        gpdf = gp2.result$result
        gpdf = gpdf[gpdf$term_size < 1000,]
        pltprefix = paste0(imgpref, 'module_enrichments_', fullpref, '_', set)
        subpmat = gpPvalMatrix(gpdf, genesets, ntop=8)
        plt = plotGpPvalMatrix(subpmat, pltprefix)

        gpdf = gp2.result$result
        gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
        pltprefix = paste0(imgpref, 'module_enrichments_', fullpref,'_', set, '_src')
        subpmat = gpPvalMatrix(gpdf, genesets, ntop=10)
        plt = plotGpPvalMatrix(subpmat, pltprefix)

        gpdf = gp2.result$result
        gpdf = gpdf[gpdf$term_size < 500,]
        pltprefix = paste0(imgpref, 'module_enrichments_', fullpref, '_', set, '_small')
        subpmat = gpPvalMatrix(gpdf, genesets, ntop=5)
        plt = plotGpPvalMatrix(subpmat, pltprefix)

        # Save the top p-values, for extended data figures:
        # -------------------------------------------------
        load(moduleenr.file)
        gpdf = gp2.result$result
        gpdf = gpdf[gpdf$term_size < 500,]

        # Get the enrichment p-values matrix:
        pmat = t(as.matrix(data.frame(gpdf$p_value)))
        rownames(pmat) = gpdf$term_name
        colnames(pmat) = names(genesets)

        # Subset to the top n terms for each:
        ntop = 10
        ids = apply(pmat, 2, ntop=ntop, top.ids)

        # Remove empty sets:
        dl = c(lapply(ids, length))
        keep = names(dl)[dl > 0]
        ids = ids[keep]

        # Get the top-pvalues for each:
        pvalsdf = NULL
        for (i in 1:length(ids)){
            mn = names(ids)[i]
            mnids = ids[[mn]]
            pn = rownames(pmat)[mnids]
            pvals = pmat[mnids, mn]
            # TODO: add order column
            pdf = data.frame(p=pvals, term=pn, mname=mn)
            pvalsdf = rbind(pvalsdf, pdf)
        }
        pvalsdf = pvalsdf[order(pvalsdf$p),]

        # Not very nice, lets report differently
        # library(ggpubr)
        # ggplot(pvalsdf, aes(term, -log10(p))) + 
        #     facet_wrap(~mname, scales='free') + 
        #     geom_bar(stat='identity') + 
        #     theme_pubr() + coord_flip()

        toppvals.file = paste0(moddir, 'module_enrichments_toppvals_',set,'_', fullpref, '_small.tsv')
        write.table(pvalsdf, toppvals.file, quote=F, row.names=F, sep="\t")
    # }
}
