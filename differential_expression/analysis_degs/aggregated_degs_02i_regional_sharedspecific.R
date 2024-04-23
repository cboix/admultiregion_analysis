#!/usr/bin/R
# -------------------------------------------------
# Plot stats on DEGs with shared/specific profiles:
# Updated: 06/23/23
# -------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(ggrastr)
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
run.intersections = FALSE
degcols = c('NS'='grey90','Down'=col.paired[2],'Up'=col.paired[6])
denrcols = c('NS'='grey90',
    'Down (1-2)'=col.paired[2],
    'Down (3+)'=col.paired[1],
    'Down (multi-CT)'='slateblue4',
    'Up (1-2)'=col.paired[6],
    'Up (3+)'=col.paired[5],
    'Up (multi-CT)'='brown4')


# Load in the annotated DEG results:
# ----------------------------------
path = 'nrad'
mstr = paste0('allmethods.regional_', path)
outpref = paste0(regdir, mstr, '.merged.sharedspecific')
allresdf = readRDS(paste0(outpref, '.rds'))


# Aggregate to plot number of genes:
ngdf = agg.rename(gene ~ set + denr + de + region, allresdf[allresdf$de != 'NS',], length, 'ngene')
# ngmat = pivot.tomatrix(ngdf[,c('ngene','set','denr')], 'set', 'ngene')

# Stats:
ngdf = merge(ngdf, agg.rename(ngene ~ set + region, ngdf, sum, 'ntot'))
ngdf$frac = ngdf$ngene / ngdf$ntot
pctdf = aggregate(frac ~ denr, ngdf, mean)


gp = ggplot(ngdf, aes(region, ngene, fill=denr)) + 
    facet_grid(set ~ de) + 
    geom_bar(stat='identity', position='stack') + 
    scale_y_continuous(expand=c(0,0), labels=scales::comma) + 
    scale_fill_manual(values=denrcols) + 
    theme_pubr()
pltprefix = paste0(imgpref, 'regional_shared_stats')
saveGGplot(gp, pltprefix, w=6, h=6)


allresdf$type = sub(".* \\(", "\\(", allresdf$denr)


# Get enrichments for only region-specific DEGs:
# ----------------------------------------------
remove.shared = TRUE
run.intersections = FALSE
for (set in keep.sets){
    for (tag in c('1-2', '3+', 'multi-CT')){
        cat(set, tag,'\n')
        suff = paste0(tag, 'only_', path, '.', set)
        type = paste0('(', tag, ')')
        fulldf = allresdf[allresdf$set == set,]
        fulldf = fulldf[(fulldf$de != 'NS') & (fulldf$type == type),]
        fulldf$de = tolower(fulldf$de)
        fulldf$pr = paste0(fulldf$region, '_', fulldf$de)

        rn = unique(fulldf$region)
        groupnames = unique(fulldf$pr)
        genesets = lapply(groupnames, function(x){
            ind = fulldf$pr == x
            fulldf$gene[ind] })
        names(genesets) = groupnames

        # Get enrichments:
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
            intdf = NULL
            save(gp2.result, intdf, file=fullenr.file)
        } else { load(fullenr.file) }

        # Region order should be fixed:
        NTOP = 2
        dereg.order = c('allregions', 'EC','HC','TH','AG','MT','PFC')
        gp2.result$result$nc = nchar(gp2.result$result$term_name)

        # for (tag in c('', '_small', '_src')){
        for (tag in c('_small')){
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


# Enrichments of fully shared genes:
# ----------------------------------
suff = paste0('multi-CT_', path, '.allCT')
fulldf = allresdf[allresdf$type == '(multi-CT)',]
fulldf = fulldf[(fulldf$de != 'NS'),]
fulldf = unique(fulldf[,c('gene','de','region')])
fulldf$de = tolower(fulldf$de)
fulldf$pr = paste0(fulldf$region, '_', fulldf$de)

groupnames = unique(fulldf$pr)
genesets = lapply(groupnames, function(x){
    ind = fulldf$pr == x
    fulldf$gene[ind] })
names(genesets) = groupnames

# Get enrichments:
run.intersections = TRUE
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
# for (tag in c('_small')){
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


