#!/usr/bin/R
# ----------------------------------------------------------
# Plot enrichments for the collated + aggregated DE results:
# - Compare enrichments in same celltype, across regions.
# Updated: 01/30/23
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
options(width=170)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
remove.shared = TRUE
run.intersections = TRUE

path = 'nrad'
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)

    # Calculate the set shared in each region:
    # ----------------------------------------
    kept.cols = c('gene','col_nm','path','region')
    alldf = c()
    for (set in keep.sets){
        print(set)
        setdf = setdflist[[set]][, kept.cols]
        setdf = setdf[setdf$col_nm != 0,]
        setdf$set = set
        alldf = rbind(alldf, setdf)
    }
    cdf = agg.rename(path ~ gene + col_nm + region, alldf, length, 'count')
    cdf = cdf[cdf$count >= 3, ]
    # head(cdf[order(cdf$count, decreasing=T),], 50)

    # TODO: Ordering of regions (not critical)
    for (set in keep.sets){
        print(set)
        suff = paste0(path, '.', set)
        fulldf = setdflist[[set]][, kept.cols]
        fulldf = fulldf[fulldf$col_nm != 0,]
        fulldf$de = ifelse(fulldf$col == 1, 'down', ifelse(fulldf$col == 2, 'up', ''))
        fulldf$pr = paste0(fulldf$region, '_', fulldf$de)
        if (remove.shared){ 
            suff = paste0('nonshared_', suff) 
            fulldf = merge(fulldf, cdf, all.x=TRUE)
            fulldf = fulldf[is.na(fulldf$count),]
        }

        rn = unique(fulldf$region)
        groupnames = unique(fulldf$pr)
        genesets = lapply(groupnames, function(x){
            ind = fulldf$pr == x
            fulldf$gene[ind] })
        names(genesets) = groupnames

        # Get enrichments for shared and specific DEGs:
        # ---------------------------------------------
        # TODO: Save no int
        if (run.intersections){ suff = paste0('withint_', suff) }
        fullenr.file = paste0(regdir, 'regionalDEGs_full_enrichments_', suff, '.rda')
        if (!file.exists(fullenr.file)){
            # Remove empty sets:
            dl = c(lapply(genesets, length))
            keep = names(dl)[dl > 0]
            genesets = genesets[keep]

            print('--GO enrichments')
            sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
            gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
                ordered_query=FALSE, multi_query=TRUE,
                sources = sources)

            # Get all of the genes in each of the term intersections now:
            # -----------------------------------------------------------
            # Get the genes in the term intersections (slow):
            if (run.intersections){
                print('--Term intersections')
                allgenes = unique(unlist(genesets))
                gp2.res.int = gprofiler2::gost(allgenes, organism='hsapiens',
                    ordered_query=FALSE, multi_query=FALSE, user_threshold=1,
                    sources = sources, evcodes=TRUE)
                intdf = gp2.res.int$result
            } else {
                intdf = NULL
            }
            save(gp2.result, intdf, file=fullenr.file)
        } else { load(fullenr.file) }

        # Region order should be fixed:
        NTOP = 3
        dereg.order = c('allregions', 'EC','HC','TH','AG','MT','PFC')
        gp2.result$result$nc = nchar(gp2.result$result$term_name)

        for (tag in c('', '_small', '_src')){
            gpdf = gp2.result$result
            gpdf = gpdf[gpdf$nc < 50,]
            if (tag == '_small'){
                gpdf = gpdf[gpdf$term_size < 500,]
            } else if (tag == '_src'){
                gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
            } else {
                gpdf = gpdf[gpdf$term_size < 1000,]
            }
            gpdf = pruneWithInt(gpdf, intdf)
            pltprefix = paste0(imgpref, 'regionalDEGs_enrichments_', NTOP, "_", suff, tag)
            subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
            subpmat = orderPvalMatrix(subpmat, dereg.order)
            plt = plotGpPvalMatrixReduced(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)
        }
    }

}

