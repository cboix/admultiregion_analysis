#!/usr/bin/R
# ------------------------------------------------------------
# Make a figure to show our cohort along multimodal reference.
# and add ABC scores to it as well
# Updated 11/29/2023
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Directories:
plotdir = paste0(imgdir, '/metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Assemble plotting dataset for MR cohort only:
# ---------------------------------------------
indmeta_tsv = 'Annotation/metadata_PFC_all_individuals_092520.tsv'
ext.meta = read.delim(indmeta_tsv, header=T)
ext.meta = ext.meta[ext.meta$projid %in% kept.individuals,]
abcdf = read.delim('Annotation/abc_scores_092023.tsv', header=T)
ext.meta = merge(ext.meta, abcdf) 

# Assemble metadata:
umeta = unique(metadata[,c('projid', 'niareagansc','braaksc','cogdx',
                           'Apoe_e4','pmi','msex', 'age_death')])
umeta = umeta[umeta$projid %in% kept.individuals,]
umeta = merge(umeta, ext.meta[,c('projid','tangles','amyloid','gpath','cogn_global_lv', 'nft','plaq_n','plaq_d', 'a_score','b_score','c_score','nia_aa_sc')])
rownames(umeta) = umeta$projid

# Get the pathology mapped to each region:
pqdf = merge(pqdf, unique(metadata[,c('rind','projid')]))
sub.pqdf = pqdf[pqdf$projid %in% kept.individuals,]
sub.pqdf$rind = NULL
plong = gather(sub.pqdf, path, val, -projid, -region)
plong$pr = paste0(plong$region, "_", plong$path)
pwide = spread(plong[,c('pr','projid','val')], pr, val)
pmat = as.matrix(pwide[,-1])
rownames(pmat) = pwide$projid

# Make a preliminary ordering of the individuals:
umeta = umeta[order(umeta$gpath, decreasing=F),]
umat = t(pmat[as.character(umeta$projid),])
ureg = sub("_.*","", rownames(umat))
upath = sub("[A-Z]*_","", rownames(umat))
freg = factor(ureg, levels=c('EC','HC','AG','MT','PFC'))
ind = order(freg)
ureg = ureg[ind]
upath = upath[ind]
umat = umat[ind,]

age.col_fun = colorRamp2(range(umeta$age_death), c("white", "slateblue")) 
# pmi.col_fun = colorRamp2(c(2, 15), c("white", "indianred")) 
gpath.col_fun = colorRamp2(c(0, max(umeta$gpath)), c("white", "indianred")) 
mat.col_fun = colorRamp2(c(0, max(umat, na.rm=T)), c("white", "royalblue")) 
pcols = brewer.pal(12, 'Paired')
abc_score.cols = c(pcols[c(2,1,5,6)])
names(abc_score.cols) = as.character(0:3)

# Make matrix and annotations:
clsplit = umeta$region
ux = 1.5
ha = HeatmapAnnotation(Sex=ifelse(umeta$msex == 0, 'female','male'), 
                       # PMI=umeta$pmi,
                       Age=umeta$age_death,
                       Apoe_e4=umeta$Apoe_e4,
                       #
                       Global_Pathology=umeta$gpath,
                       Braak_Stage=umeta$braaksc,
                       NIA_Reagan=umeta$niareagansc,
                       Cognitive_Diagnosis=umeta$cogdx,
                       #
                       ABC_Score=umeta$nia_aa_sc,
                       a_score=umeta$a_score,
                       b_score=umeta$b_score,
                       c_score=umeta$c_score,
                       # 
                       annotation_name_gp = gpar(fontsize=5),
                       simple_anno_size = unit(ux, 'mm'),
                       gap=unit(0,'mm'),
                       col=list(NIA_Reagan=colvals[['niareagansc']],
                                Apoe_e4=c('no'='grey95','yes'='grey70'),
                                Age=age.col_fun,
                                Global_Pathology=gpath.col_fun,
                                # PMI=pmi.col_fun,
                                Braak_Stage=colvals[['braaksc']],
                                Cognitive_Diagnosis=colvals[['cogdx']],
                                ABC_Score=abc_score.cols, 
                                a_score=abc_score.cols,
                                b_score=abc_score.cols,
                                c_score=abc_score.cols,
                                Sex=colvals[['sex']]))

udsplit = upath
hb = rowAnnotation(Region=ureg, col=list(Region=reg.cols),
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(ux, 'mm'),
    gap=unit(0,'mm'))

# png(paste0(imgpref, 'individual_metadata_heatmap.png'), res=400, units='in', width=11, height=5)
scale = 5/6
ht = Heatmap(umat, 
    name='Path.\nDensity', 
    col=mat.col_fun,
    use_raster=TRUE,
    top_annotation=ha, 
    right_annotation=hb,
    row_split=udsplit,
    column_split=clsplit, 
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    width=ncol(umat) * unit(ux, 'mm') * scale,
    height=nrow(umat) * unit(ux, 'mm'),
)


pltprefix = paste0(imgpref, 'individual_metadata_heatmap_abcscore')
h = 1 + 1 / 15 * nrow(umat)
w = 5 + 1 / 15 * ncol(umat) * scale
saveHeatmap(ht, pltprefix, w=w, h=h)

