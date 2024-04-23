#!/usr/bin/R
# ----------------------------------------------------------
# Plot enrichments for the collated + aggregated DE results:
# Shared and specific pathways up and down.
# Updated: 03/22/22
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


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

plotGpPvalMatrix = function(subpmat, pltprefix, cluster_columns=FALSE){
    # pltmat = t(reord(t(subpmat)))
    pltmat = t(diag.mat2(t(subpmat))[[1]])
    pltmat = pltmat[rev(rownames(pltmat)),]
    maxrow = colnames(pltmat)[apply(pltmat, 1, which.max)]
    rowsplit = rep("Down", nrow(pltmat))
    rowsplit[grep("_up$", maxrow)]= 'Up'
    csplit = ifelse(1:ncol(pltmat) %in% grep("_up$", colnames(pltmat)), 'Up','Down')
    ux = 1.5
    plt = Heatmap(pltmat,
                  cluster_columns=cluster_columns,
                  cluster_rows=FALSE,
                  width = ncol(pltmat)*unit(ux, "mm"), 
                  height = nrow(pltmat)*unit(ux, "mm"),
                  row_split=rowsplit,
                  column_split=csplit,
                  use_raster=TRUE,
                  border_gp=gpar(color='black', lwd=.5),
                  name='-log10(p)',
                  col=c('white',colb))

    h = 2.25 + 1 / 15 * nrow(pltmat)
    w = 5 + 1 / 15 * ncol(pltmat)
    saveHeatmap(plt, pltprefix, w=w, h=h)
}


pruneWithInt = function(gpdf, intdf, cutoff=0.75){
    intdf = intdf[intdf$term_id %in% gpdf$term_id,]
    intdf = intdf[order(intdf$p_value),]
    # Use this to prune the enriched terms:
    sel.ind = selGOterms(intdf$intersection, cutoff=cutoff)
    int.terms = intdf$term_id[sel.ind]
    # Add terms not in intdf:
    all.gpterms = gpdf$term_id
    uq.gpterms = all.gpterms[!(all.gpterms %in% intdf$term_id)]
    kept.terms = c(int.terms, uq.gpterms)
    print(length(kept.terms) / length(all.gpterms))
    return(gpdf[gpdf$term_id %in% kept.terms,])
}


# Aggregate the DEGs and results across all runs:
# -----------------------------------------------
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
load(fullaggrda)
sets = names(setdflist)


# Reduce to a specific AD variable + keep major glial only:
# ---------------------------------------------------------
path = 'cogdxad'
reduce.sets = TRUE
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
pathlist = unique(totnsigdf$path)

for (path in pathlist){
    degdf = NULL
    degs = list()
    for (set in keep.sets){
        setdf = setdflist[[set]]
        setdf = setdf[setdf$path == path, ]
        setdf$set = set
        degdf = rbind(degdf, setdf)
        degs[[paste0(set, "_up")]] = setdf$gene[setdf$col_nm == 2]
        degs[[paste0(set, "_down")]] = setdf$gene[setdf$col_nm == 1]
    }


    # Reduce to a merged data.frame (for "All")
    # -----------------------------------------
    cdf = aggregate(path ~ gene + col_nm, degdf, length)
    # Break ties:
    cdf$path = cdf$path + .1 * cdf$col_nm
    cdf = merge(cdf, aggregate(path ~ gene, cdf, max))
    cdf = cdf[cdf$path >= 3, ]
    cdf$path = path
    desuff = 'All_All_allregions'
    fname = paste0(regdir, 'allmethods.', desuff, '_', path, '.merged.tsv.gz')
    write.table(cdf, gzfile(fname), quote=F, sep="\t", row.names=F)


    # Separate DEGs into shared / specific:
    # -------------------------------------
    sigdf = degdf[degdf$col_nm != 0,]
    repdf = agg.rename(pc ~ col_nm + gene, sigdf, length, 'ntimes')

    pdf(paste0(imgpref, 'barplot_ndeg_shared_', path ,'.pdf'), width=6, height=6)
    barplot(table(repdf$ntimes), border=F, ylab='Number of DEG', xlab='Number of cell types')
    dev.off()


    # Plot on ndeg the number of genes that are shared:
    # -------------------------------------------------
    # Define shared as up or down in at least 3 cell types::
    nset = 3
    up_rep = repdf[repdf$ntimes >= nset & repdf$col_nm == 2,'gene']
    down_rep = repdf[repdf$ntimes >= nset & repdf$col_nm == 1,'gene']
    sigdf$intop = sigdf$gene %in% c(up_rep, down_rep)
    print(length(c(up_rep, down_rep)))
    mtop = aggregate(intop ~ set, sigdf, mean)

    sharedup.file = paste0(regdir, 'allmethods.allmajor.', path, '.sharedup.tsv')
    shareddw.file = paste0(regdir, 'allmethods.allmajor.', path, '.shareddw.tsv')
    write.table(up_rep, sharedup.file, quote=F, row.names=F, sep="\t", col.names=F)
    write.table(down_rep, shareddw.file, quote=F, row.names=F, sep="\t", col.names=F)

    # Number of significant (full)
    nsigdf = spread(agg.rename(pc ~ col_nm + set, sigdf, length, 'ndeg'), col_nm, ndeg)
    names(nsigdf)[2:3] = c('down','up')
    subdf = nsigdf
    subdf = subdf[order(subdf$set),]
    subdf$set = factor(subdf$set, levels=rev(unique(subdf$set)))

    # Number of significant (specific)
    nsigdf2 = spread(agg.rename(pc ~ col_nm + intop + set, sigdf, length, 'ndeg'), col_nm, ndeg)
    names(nsigdf2)[3:4] = c('down','up')
    subdf2 = nsigdf2[nsigdf2$intop == FALSE,]
    pcols = brewer.pal(12, 'Paired')

    # Plot 
    gplot = ggplot(subdf, aes(set, up)) + 
        geom_bar(position='dodge',stat='identity', fill=pcols[5]) + 
        geom_bar(data=subdf, aes(set, -down), position='dodge',stat='identity', fill=pcols[1]) + 
        geom_bar(data=subdf2, aes(set, up), position='dodge',stat='identity', fill=pcols[6]) + 
        geom_bar(data=subdf2, aes(set, -down), position='dodge',stat='identity', fill=pcols[2]) + 
        geom_text(data=subdf, aes(set, -down, label=down), color='black') + 
        geom_text(data=subdf, aes(set, up, label=up), color='black') + 
        labs(x='Cell type', y='Number of DEGs') + 
        theme_pubr() + coord_flip()
    ggsave(paste0(imgpref, 'ndeg_barplot_bytop_', path, '.png'), gplot, dpi=400, units='in', width=5,height=3)
    ggsave(paste0(imgpref, 'ndeg_barplot_bytop_', path, '.pdf'), gplot, width=5,height=3)


    # Get enrichments for shared and specific DEGs:
    # ---------------------------------------------
    use.shared = TRUE
    suff = paste0(path, ifelse(use.shared, '_wshared', '_all'))
    # gliaenr.file = paste0(regdir, 'aggregatedDEGs_glia_enrichments_', suff, '.rda')
    fullenr.file = paste0(regdir, 'aggregatedDEGs_full_enrichments_', suff, '.rda')
    if (!file.exists(fullenr.file)){
        if (use.shared){
            genesets = list('shared_up'=up_rep, 
                'shared_down'=down_rep)
        } else { genesets = list() }

        # Add specific:
        for (set in unique(sigdf$set)){
            if (use.shared){
                subsetdf = sigdf[(sigdf$set == set) & (sigdf$intop == FALSE),]
            } else {
                subsetdf = sigdf[(sigdf$set == set),]
            }

            genesets[[paste0(set, '_up')]] = subsetdf$gene[subsetdf$col_nm == 2]
            genesets[[paste0(set, '_down')]] = subsetdf$gene[subsetdf$col_nm == 1]
        }

        # Remove empty sets:
        dl = c(lapply(genesets, length))
        keep = names(dl)[dl > 0]
        genesets = genesets[keep]

        sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
        gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
            ordered_query=FALSE, multi_query=TRUE,
            sources = sources)


        # Get all of the genes in each of the term intersections now:
        # -----------------------------------------------------------
        # Get the genes in the term intersections (slow):
        allgenes = unique(unlist(genesets))
        gp2.res.int = gprofiler2::gost(allgenes, organism='hsapiens',
            ordered_query=FALSE, multi_query=FALSE, user_threshold=1,
            sources = sources, evcodes=TRUE)
        intdf = gp2.res.int$result

        # Get the top N genes per term:
        # degs = dedf$gene
        # intdf$topgenes = sapply(intdf$intersection, function(x, ntop=3){
        #     x = strsplit(x, ',')[[1]]
        #     x = head(degs[degs %in% x], ntop)
        #     return(paste0(x, collapse=','))
        #     })

        # Mapping term ids to genes:
        # mapgenes = intdf$topgenes
        # names(mapgenes) = intdf$term_id
        # gp2df$topgenes = mapgenes[gp2df$term_id]

        save(gp2.result, intdf, file=fullenr.file)
    } else { load(fullenr.file) }
    
    gpdf = gp2.result$result
    gpdf = gpdf[gpdf$term_size < 1000,]
    gpdf = pruneWithInt(gpdf, intdf)
    pltprefix = paste0(imgpref, 'aggregatedDEGs_glia_enrichments_', suff)
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=8)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE)

    gpdf = gp2.result$result
    gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
    gpdf = pruneWithInt(gpdf, intdf)
    pltprefix = paste0(imgpref, 'aggregatedDEGs_glia_enrichments_', suff, '_src')
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=10)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE)

    gpdf = gp2.result$result
    gpdf = gpdf[gpdf$term_size < 500,]
    gpdf = pruneWithInt(gpdf, intdf, cutoff=0.5)
    pltprefix = paste0(imgpref, 'aggregatedDEGs_glia_enrichments_', suff, '_small')
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=5)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE)

}
