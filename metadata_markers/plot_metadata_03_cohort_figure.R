#!/usr/bin/R
# ------------------------------------------------------------
# Make a figure to show our cohort along multimodal reference.
# Also plot the full ROSMAP cohort within that figure (alt)
# Updated 04/26/2021
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(qvalue)
library(ComplexHeatmap)
library(circlize)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# ---------------------------------------------
# Assemble plotting dataset for MR cohort only:
# ---------------------------------------------
emeta = read.delim('Annotation/metadata_PFC_all_individuals_092520.tsv', header=T)
load(file=paste0(datadir, prefix, '.final_noMB.cell_labels.Rda'))
ddf = unique(cellmeta[,c('projid','region')]) # scDatasets

# Assemble metadata:
umeta = unique(metadata[,c('projid', 'niareagansc','braaksc','cogdx',
                           'Apoe_e4','pmi','msex', 'age_death')])
umeta = umeta[umeta$projid %in% ddf$projid,]
umeta = merge(umeta, emeta[,c('projid','tangles','amyloid','gpath','cogn_global_lv', 'nft','plaq_n','plaq_d')])
rownames(umeta) = umeta$projid

# Get the pathology mapped to each region:
pqdf = NULL
regmap = c('AG','HC','PFC','MT','EC')
names(regmap) = c('ag','hip','mf','mt','ec')
rmeta = merge(emeta, data.frame(region=regmap))
pathlist = c('nft','plaq_d','plaq_n')
for (path in pathlist){
    vars = colnames(rmeta)[grep(path, colnames(rmeta))]
    vars = vars[vars != path]
    submeta = unique(rmeta[,c('projid','region', vars)])
    slong = gather(submeta, path, value, -projid, -region)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    if (is.null(pqdf)){
        pqdf = slong[,c('projid','value','region')]
        names(pqdf)[2] = path
    } else { 
        sub.pqdf = slong[,c('projid','value','region')]
        names(sub.pqdf)[2] = path
        pqdf = merge(pqdf, sub.pqdf)
    }
}
plong = gather(pqdf, path, val, -projid, -region)
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

# Make matrix and annotations:
clsplit = umeta$region
ha = HeatmapAnnotation(Sex=ifelse(umeta$msex == 0, 'female','male'), 
                       # PMI=umeta$pmi,
                       Age=umeta$age_death,
                       Apoe_e4=umeta$Apoe_e4,
                       #
                       GPath=umeta$gpath,
                       Braak=umeta$braaksc,
                       AD=umeta$niareagansc,
                       Cognition=umeta$cogdx,
                       col=list(AD=colvals[['niareagansc']],
                                Apoe_e4=c('no'='grey95','yes'='grey70'),
                                Age=age.col_fun,
                                GPath=gpath.col_fun,
                                # PMI=pmi.col_fun,
                                Braak=colvals[['braaksc']],
                                Cognition=colvals[['cogdx']],
                                Sex=colvals[['sex']]))

udsplit = upath
hb = rowAnnotation(Region=ureg, col=list(Region=reg.cols))

# png(paste0(imgpref, 'individual_metadata_heatmap.png'), res=400, units='in', width=11, height=5)
pdf(paste0(imgpref, 'individual_metadata_heatmap.pdf'), width=11, height=5)
Heatmap(umat, name='Path.\nDensity', 
        # col=colb,
        col=mat.col_fun,
        use_raster=TRUE,
        top_annotation=ha, 
        column_split=clsplit, 
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_column_names=FALSE,
        row_split=udsplit,
        right_annotation=hb
)
dev.off()

# table(ddf$region)
