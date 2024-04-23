#!/usr/bin/R
# ----------------------------------------------------------
# Plot the per-sample scores of regional DEGs (from nrad)
# - Compare across regions, staging.
# Updated: 02/27/23
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
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


load_psbulk_matrix = function(set){
    ststr = sub('_[A-Za-z]+$', '', set)
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
    load(psdata.rda)
    # If Mic, remove T, CAMs from psbulk data:
    if (set == 'Mic_Immune_Mic'){
        ind = !(ps.data$meta$cell_type_high_resolution %in% c('T cells', 'CAMs'))
        ps.data$meta = ps.data$meta[ind,]
        ps.data$mat = ps.data$mat[,rownames(ps.data$meta)]
    }
    # Merge to sample-level (individual x region):
    ps.data = aggregate_psbulk_samplelevel(ps.data)
    return(ps.data)
}


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
remove.shared = TRUE
run.intersections = FALSE
degcols = c('NS'='grey90','Down'=col.paired[2],'Up'=col.paired[6])
denrcols = c('NS'='grey90',
    'Down (1-2)'=col.paired[2], 'Down (3+)'=col.paired[1], 'Down (multi-CT)'='slateblue4',
    'Up (1-2)'=col.paired[6], 'Up (3+)'=col.paired[5], 'Up (multi-CT)'='brown4')
region.cols = c('allregions'='grey85', reg.cols[reg.nomb])

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

        # Load pseudobulk data:
        ps.data = load_psbulk_matrix(set)

        # By all genes in region
        scdf = c()
        for (region in unique(fulldf$region)){
            subdf = fulldf[fulldf$region == region,]
            uplist = subdf[subdf$col_nm == 2,'gene']
            dwlist = subdf[subdf$col_nm == 1,'gene']
            ngenes = length(uplist) + length(dwlist)
            sc = (colSums(ps.data$mat[uplist,, drop=F]) - colSums(ps.data$mat[dwlist,, drop=F])) / ngenes
            df = data.frame(score=sc, pr=names(sc), deregion=region, path=path)
            scdf = rbind(scdf, df)
        }
        scdf = merge(scdf, ps.data$meta)
        scdf = merge(scdf, unique(metadata[metadata$region == 'PFC',c('projid', 'nrad', 'cogn_global_lv')]))
        scdf = merge(scdf, unique(metadata[,c('projid', 'region', 'rind')]))

        # Plot these scores as violins:
        gp = ggplot(scdf, aes(region, score, fill=nrad)) + 
            facet_wrap(~deregion, scales='free_y',nrow=2) + 
            scale_fill_manual(values=c('AD'='red','CTRL'='grey85')) +
            geom_violin() + theme_pubr()
        pltprefix = paste0(imgpref, 'regional_crossDEscores_violin.', set)
        saveGGplot(gp, pltprefix, w=9, h=6)

        # Plot these scores on plaque / trajectory:
        scdf2 = merge(scdf, pqdf) # No measurements for TH
        gp = ggplot(scdf2, aes(log1p(plaq_n), score, color=region)) + 
            facet_wrap(~deregion, scales='free_y') + 
            scale_color_manual(values=c(region.cols)) +
            geom_point() + geom_smooth(method='loess') + 
            theme_pubr()
        pltprefix = paste0(imgpref, 'regional_crossDEscores_scatter.plaq_n.', set)
        saveGGplot(gp, pltprefix, w=8, h=8)

        gp = ggplot(scdf2, aes(log1p(nft), score, color=region)) + 
            facet_wrap(~deregion, scales='free_y') + 
            scale_color_manual(values=c(region.cols)) +
            geom_point() + geom_smooth(method='loess') + 
            theme_pubr()
        pltprefix = paste0(imgpref, 'regional_crossDEscores_scatter.nft.', set)
        saveGGplot(gp, pltprefix, w=8, h=8)
    }
}


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Collect all of the DEGs at a major cell type level:
# ---------------------------------------------------
alldf = c()
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)

    kept.cols = c('gene','col_nm','path','region')
    for (set in keep.sets){
        setdf = setdflist[[set]][, kept.cols]
        setdf = setdf[setdf$col_nm != 0,]
        setdf$set = set
        alldf = rbind(alldf, setdf)
    }
}


# Subset DEGs to one major cell type, annotate with modules:
# ----------------------------------------------------------
runset = 'Mic_Immune'
set = 'Mic_Immune_Mic'
setdf = alldf[alldf$set == set,]
ps.data = load_psbulk_matrix(set)
keep.mod = as.numeric(names(nmod)[nmod >= 10])


# Subset to conditions:
# ---------------------
use.any.region = FALSE
use.core = TRUE 
if (!use.any.region){ setdf = setdf[setdf$region == 'allregions',] }
if (use.core){
    genemap = cmlist[[runset]]
} else {
    genemap = gmlist[[runset]]
}
setdf$module = genemap[setdf$gene]
setdf = setdf[!is.na(setdf$module),]


# Get either all DE or only DE in allregions, and score modules:
# --------------------------------------------------------------
scdf = c()
for (mn in keep.mod){
    # TODO: Also allow allgenes vs. core genes
    subdf = setdf[setdf$module == mn,]
    uplist = unique(subdf$gene[subdf$col_nm == 2]) # NOTE: May overlap
    dwlist = unique(subdf$gene[subdf$col_nm == 1])

    ngenes = length(uplist) + length(dwlist)
    sc = (colSums(ps.data$mat[uplist,, drop=F]) - colSums(ps.data$mat[dwlist,, drop=F])) / ngenes
    df = data.frame(score=sc, pr=names(sc), module=mn)
    scdf = rbind(scdf, df)
}
scdf = merge(scdf, ps.data$meta)
scdf = merge(scdf, unique(metadata[metadata$region == 'PFC',c('projid', 'nrad', 'cogn_global_lv')]))
scdf = merge(scdf, unique(metadata[,c('projid', 'region', 'rind')]))

# Plot these scores as violins:
gp = ggplot(scdf, aes(region, score, fill=nrad)) + 
    facet_wrap(~module, scales='free_y') + 
    scale_fill_manual(values=c('AD'='red','CTRL'='grey85')) +
    geom_violin() + theme_pubr()
pltprefix = paste0(imgpref, 'regional_moduleDEscores_violin.', set)
saveGGplot(gp, pltprefix, w=12, h=12)

# Plot these scores on plaque / trajectory:
scdf2 = merge(scdf, pqdf) # No measurements for TH
gp = ggplot(scdf2, aes(log1p(plaq_n), score, color=region)) + 
    facet_wrap(~module, scales='free_y') + 
    scale_color_manual(values=c(region.cols)) +
    geom_point() + geom_smooth(method='loess') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'regional_moduleDEscores_scatter.plaq_n.', set)
saveGGplot(gp, pltprefix, w=8, h=8)

gp = ggplot(scdf2, aes(log1p(nft), score, color=region)) + 
    facet_wrap(~module, scales='free_y') + 
    scale_color_manual(values=c(region.cols)) +
    geom_point() + geom_smooth(method='loess') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'regional_moduleDEscores_scatter.nft.', set)
saveGGplot(gp, pltprefix, w=8, h=8)





# Repeat these plots, but by subset of genes (GO or modules)
# ----------------------------------------------------------
# If modules - any genes that are DE for any condition (nrad, plaq, etc.) (in any region???)
# First score each region by each other


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

