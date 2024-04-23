#!/usr/bin/R
# ----------------------------------------------------------
# Plot DEGs with highlighted enrichments 
# for the collated + aggregated DE results:
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

    allresdf = NULL
    for (set in keep.sets){
        print(set)
        suff = paste0(path, '.', set)
        fulldf = setdflist[[set]]
        fulldf$de = ifelse(fulldf$col_nm == 1, 'Down', ifelse(fulldf$col_nm == 2, 'Up', 'NS'))
        fulldf$pr = paste0(fulldf$region, '_', fulldf$de)
        fulldf = merge(fulldf, cdf, all.x=TRUE) # Shared within regions across cell types

        # Count up shared between regions in this cell type:
        rdf = agg.rename(path ~ gene + col_nm, fulldf, length, 'nreg')
        fulldf = merge(fulldf, rdf, all.x=TRUE)
        fulldf$nregion = ifelse(is.na(fulldf$count), ifelse(fulldf$nreg >= 3, '3+', '1-2'), 'multi-CT')
        fulldf$denr = paste0(fulldf$de, ' (', fulldf$nregion, ')')
        fulldf$denr[fulldf$de == 'NS'] = 'NS'
        fulldf$set = set
        # table(fulldf$denr)
        allresdf = rbind(allresdf, fulldf)

        # Volcano plots, without pathways labeled:
        # Labelled which points are region-unique vs. shared in CT
        # Also added genes that are shared across cell types
        # --------------------------------------------------------
        for (region in unique(fulldf$region)){
            subdf = fulldf[fulldf$region == region,]
            subdf = subdf[order(subdf$p_nb),]
            NTOP = 20
            labdf = rbind(head(subdf[subdf$col_nm == 2,], NTOP),
                head(subdf[subdf$col_nm == 1,], NTOP))

            # NOTE: Plotting with Nebula parameters for consistency
            subdf = subdf[order(subdf$col_nm, decreasing=F),]
            gp = ggplot(subdf, aes(logFC_nb, -log10(p_nb), color=denr)) + 
                scale_color_manual(values=denrcols, name='DE status:') +
                rasterize(geom_point(cex=.25), dpi=450) +
                scale_y_continuous(expand=c(0,0)) + 
                geom_vline(xintercept=0, lty='dotted') + 
                geom_text_repel(data=labdf, aes(logFC_nb, -log10(p_nb), label=gene, color=denr), 
                    box.padding=0.075, cex=3, max.overlaps=30, min.segment.length=0.1) + 
                labs(x='log Fold Change', y='log10 p-value', title=region) + 
                theme_pubr() + theme(legend.position='none')
            pltprefix = paste0(imgpref, 'regional_volcano.', suff, '.', region)
            saveGGplot(gp, pltprefix, w=4, h=4)
        }
    }

    outpref = paste0(regdir, mstr, '.merged.sharedspecific')
    saveRDS(allresdf, file=paste0(outpref, '.rds'))
    write.table(allresdf, file=gzfile(paste0(outpref, '.tsv.gz')), quote=F, sep="\t", row.names=F)
}


#         # TODO: PLOT HIGHLIGHTING PATHWAYS?
#         # Get enrichments for shared and specific DEGs:
#         # ---------------------------------------------
#         # TODO: Save no int
#         if (run.intersections){ suff = paste0('withint_', suff) }
#         fullenr.file = paste0(regdir, 'regionalDEGs_full_enrichments_', suff, '.rda')
#         # save(gp2.result, intdf, file=fullenr.file)
#         load(fullenr.file) 
#     }
# }

