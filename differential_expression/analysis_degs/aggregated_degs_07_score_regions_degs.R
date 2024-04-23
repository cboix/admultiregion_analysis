#!/usr/bin/R
# ----------------------
# Score regions by DEGs:
# Updated: 03/29/22
# ----------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(circlize)
library(ComplexHeatmap)

library(progress)
print(version)
options(width=175)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_regionpred_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Arguments for runs:
# -------------------
# Parameters:
path = 'plaq_d'
ntop = 250
use.weighted = TRUE
optstr = paste0('_', path, '_n', ntop)
if (use.weighted){ optstr = paste0(optstr, '_weighted') }
# TODO: imgpref for these parameters:

sets = c('Ast'='Ast_Ast','Mic_Immune'='Mic_Immune_Mic',
    'Opc'='Opc_Opc', 'Oli'='Oli_Oli','Inh'='Inh_Inh',
    'Exc'='Exc_Exc', 'All'='All_All')
# TODO: Fix "All" to have the merged coefficient (for weighted)

score.rds = paste0(srdir, 'region_descores_', optstr, '.rds')
if (!file.exists(score.rds)){
    scoredf = c()
    for (i in 1:length(sets)){
        subset = sets[i]
        runset = names(sets)[i]
        print(runset)

        # Load pseudobulk data and aggregate at the sample level:
        # -------------------------------------------------------
        source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
        psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
        load(psdata.rda)
        ps.data = aggregate_psbulk_samplelevel(ps.data)


        # Write this data out for python models:
        # --------------------------------------
        mat.tsv = paste0(srdir, 'pseudobulk_data_samplelvl_', runset, '_matrix.tsv.gz')
        meta.tsv = paste0(srdir, 'pseudobulk_data_samplelvl_', runset, '_metadata.tsv.gz')
        if (!file.exists(meta.tsv)){
            print(paste0("Writing matrices for: ", runset))
            meta.cols = c('projid', 'region', 'rind', 'cogdxad', 'nrad',
                'age_death','msex','Apoe_e4','pmi')
            umeta = merge(ps.data$meta, unique(metadata[, meta.cols]))
            umeta = merge(umeta, pqdf, all.x=TRUE)
            pmat = as.matrix(ps.data$mat[, umeta$pr])
            write.table(pmat, gzfile(mat.tsv), quote=F, sep="\t")
            write.table(umeta, gzfile(meta.tsv), quote=F, sep="\t", row.names=F)
        }


        # Load the DEGs:
        # --------------
        if (runset == 'All'){
            desuff = 'All_All_allconditions'
            full.rds = paste0(regdir, 'allmethods.', desuff, '.merged.rds')
        } else {
            full.rds = paste0(regdir, 'aggregated_fullset.', subset, '.rds')
        }

        dedf = readRDS(full.rds)
        dedf = dedf[dedf$path == path & dedf$region == 'allregions',]
        if (runset == 'All'){ dedf = dedf[order(dedf$count, decreasing=T),] }
        tdf = head(dedf[dedf$col_nm == 2,], ntop)
        bdf = head(dedf[dedf$col_nm == 1,], ntop)
        genes = c(tdf$gene, bdf$gene)
        if (use.weighted & runset != 'All') {
            weights = c(tdf$logFC_nb, bdf$logFC_nb)
        } else {
            weights = c(rep(1, nrow(bdf)), rep(-1, nrow(bdf))) / 100
        }

        # Score matrix:
        pmat = ps.data$mat[genes,]
        scores = colSums(sweep(pmat, 1, weights, '*'))
        df = cbind(ps.data$meta, score=scores, ct=runset)
        df$zscore = (df$score - mean(df$score)) / sd(df$score)
        scoredf = rbind(scoredf, df)
    }
    saveRDS(scoredf, file=score.rds)
} else {
    scoredf = readRDS(score.rds)
}


# Plot these scores as heatmaps:
# ------------------------------
scoredf$cr = with(scoredf, paste0(ct, '_', region))
cmat = pivot.tomatrix(scoredf[,c('cr','projid','zscore')], 'cr','zscore')
csplit = sub('_.*', '', colnames(cmat))

# Regional metadata matrix:
plong = merge(pqdf, unique(metadata[,c('rind', 'region','projid')]))
plong$rind = NULL
plong = gather(plong, path, val, -projid, -region)
plong$pr = paste0(plong$region, "_", plong$path)
pwide = spread(plong[,c('pr','projid','val')], pr, val)
pmat = as.matrix(pwide[,-1])
rownames(pmat) = pwide$projid
pmat = pmat[as.character(rownames(cmat)),]

# Add the metadata for the individual-level:
umeta = unique(metadata[metadata$region == 'PFC',
    c('projid','nrad','cogdxad','cogdx', 'niareagansc', 'age_death',
        'msex', 'braaksc', 'Apoe_e4','gpath','tangles','amyloid', 'cogn_global_lv')])
rownames(umeta) = umeta$projid
umeta = umeta[as.character(rownames(cmat)),]

mat.col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")) 
age.col_fun = colorRamp2(range(umeta$age_death), c("white", "slateblue")) 
pmi.col_fun = colorRamp2(c(2, 15), c("white", "indianred")) 
gpath.col_fun = colorRamp2(c(0, max(umeta$gpath)), c("white", "indianred")) 
gcog.col_fun = colorRamp2(c(min(umeta$cogn_global_lv), max(umeta$cogn_global_lv)),
    c("red", 'white')) 

# Make metadata annotation:
ux = 1.5
ha = HeatmapAnnotation(
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(ux, 'mm'),
    Sex=ifelse(umeta$msex == 0, 'female','male'), 
    gap=unit(0,'mm'),
    Age=umeta$age_death,
    Apoe_e4=umeta$Apoe_e4,
    GPath=umeta$gpath,
    Braak=umeta$braaksc,
    AD=umeta$niareagansc,
    Cognition=umeta$cogdx,
    GlobalCog=umeta$cogn_global_lv,
    col=list(AD=colvals[['niareagansc']],
        Apoe_e4=c('no'='grey95','yes'='grey70'),
        Age=age.col_fun,
        GPath=gpath.col_fun,
        GlobalCog=gcog.col_fun,
        Braak=colvals[['braaksc']],
        Cognition=colvals[['cogdx']],
        Sex=colvals[['sex']]), which='row')

cht = Heatmap(cmat, 
    use_raster=TRUE, 
    column_split=csplit,
    # cluster_columns=FALSE,
    border_gp=gpar(color='black', lwd=.5),
    col=mat.col_fun,
    width=ncol(cmat) * unit(ux, 'mm'),
    height=nrow(cmat) * unit(ux, 'mm'), 
    show_row_names=FALSE)

pqht = Heatmap(pmat, 
    use_raster=TRUE, 
    cluster_columns=FALSE,
    border_gp=gpar(color='black', lwd=.5),
    column_split=sub("[A-Z]+_","", colnames(pmat)),
    col=viridis(50),
    right_annotation=ha,
    width=ncol(pmat) * unit(ux, 'mm'),
    height=nrow(pmat) * unit(ux, 'mm'), 
    show_row_names=FALSE)

plt = cht + pqht
# plt = pqht + cht

pltprefix = paste0(imgpref, 'region_descores_heatmap', optstr)
w = 2 + (ncol(cmat) + ncol(pmat)) / 15
h = 2 + nrow(cmat) / 15 
saveHeatmap(plt, pltprefix, w=w, h=h)


