#!/usr/bin/R
# ------------------------------------------------------
# Compare OPC fraction to OPC genes at pseudobulk level:
# Updated 02/17/2021
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)
library(patchwork)
library(emmeans)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'markers/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Run composition analysis for each of these subsets:
# ---------------------------------------------------
remove.batches = TRUE
suff = '_subset_final_noMB'

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c('All', remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
subtypes = unique(ctdf$cls)
celltypes = unique(ctdf$major.celltype)


# Load in the full opc data for these subtypes:
# ---------------------------------------------
subset = 'Opc'
ststr = gsub("/","_", subset)

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
if (!file.exists(psdata.rda)){
    ps.data = load_pseudobulk_dataset(subset, subtypes, reg.nomb)
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }


# Aggregate matrix at the pseudobulk sample level:
# ------------------------------------------------
pmat = ps.data$mat
umeta = ps.data$meta

umeta$pr = paste0(umeta$projid, '_', umeta$region)
ptypes = unique(umeta$pr)
tform = make.tform(umeta$pr, u=ptypes, norm=F)
tform = sweep(tform, 1, umeta$ncell, '*')
tform = sweep(tform, 2, apply(tform, 2, sum), '/')
pmat = pmat %*% tform
umeta = aggregate(ncell ~ region + projid + pr, umeta, sum)
pmat = pmat[, umeta$pr]

# Further annotate the pseudo-bulk metadata:
# ------------------------------------------
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
# Remove very low abundance batches + will use for weight
umeta = umeta[umeta$ncell > 10,] 
rownames(umeta) = umeta$pr
pmat = pmat[, umeta$pr]


# Perform regression of OPC expression vs. OPC count:
# ---------------------------------------------------
subct = ctdf[ctdf$cls == 'Opc',]
subct$pr = paste0(subct$projid, '_', subct$region)
pmat = pmat[, subct$pr]

gform = asform(c('lfrac ~ x + region * msex'))
# subct$lfrac = log10(subct$count / subct$other)
subct$lfrac = log10(subct$count / subct$total)

resdf = c()
for (gene in rownames(pmat)){
    subct$x = pmat[gene,]
    fit = lm(gform, subct)
    cfit = data.frame(coefficients(summary(fit)))
    colnames(cfit) = c('Est','SE','t','p')
    cfit = cfit['x',, drop=F]
    cfit$gene = gene
    resdf = rbind(resdf, cfit)
}
resdf = resdf[order(resdf$p), ]

gene = 'MARCKS'
gene = 'NAV2'
subct$x = pmat[gene,]
gp = ggplot(subct, aes(x, count / total)) + 
    geom_point(cex=.5) + 
    geom_smooth(method='lm') + 
    stat_cor(output.type='text', label.sep='\n') + 
    labs(x=gene) +
    theme_pubr()
pltprefix = paste0(imgpref, 'opcfrac_',gene, '_scatter')
saveGGplot(gp, pltprefix, 5, 4)


genes = c('EGFR','MARCKS','OLIG1','OLIG2','VCAN','BCAN','SAMHD1','APOD','NFKBIA','SLAIN1','MT3','ZBTB20')
genes = c('XYLT1', 'TRIO','NAV2','SGCD','GPR158','SORCS1','CAMK2D')
resdf[resdf$gene %in% genes,]

