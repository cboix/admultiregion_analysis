#!/usr/bin/R
# ---------------------------------------------------------
# Plot basic marker differences between astrocyte subtypes:
# Updated 11/26/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)
library(patchwork)

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

subset = 'Ast'
ststr = gsub("/","_", subset)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
subtypes = unique(ctdf$cls)
celltypes = unique(ctdf$major.celltype)


# Load in the full ast data for these subtypes:
# ---------------------------------------------
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
if (!file.exists(psdata.rda)){
    ps.data = load_pseudobulk_dataset(subset, subtypes, reg.nomb)
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }


# Further annotate the pseudo-bulk metadata:
# ------------------------------------------
pmat = ps.data$mat
umeta = ps.data$meta
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(pmat),]
# Remove very low abundance batches + will use for weight
umeta = umeta[umeta$ncell > 100,] 


# Compare the GRM3 astrocyte markers to the MEIS2 Inhibitory markers:
# -------------------------------------------------------------------
srpref = paste0(srdir, 'difftl_pseudobulk_markers_')
processRegdf = function(df, txt){
    df = df[df$var == 'is.stTRUE',]
    df = df[order(df$p),]
    df = df[order(-sign(df$Est)),]
    df$p.adj = p.adjust(df$p)
    df$rank = 1:nrow(df)
    varlist = c('Est','t', 'p','p.adj','rank')
    df = df[,c(varlist, 'symbol')]
    names(df) = c(paste0(txt, '_', varlist), 'symbol')
    return(df)
}

# Load each:
st = 'Ast GRM3'
sub.ststr = gsub("/","_",gsub(" ","_", st))
subtype.diff.rda = paste0(srpref, ststr, '_', sub.ststr, '.rda')
load(subtype.diff.rda)
ast.regdf = processRegdf(est.regdf, 'ast')

cutoff = 100
st = 'Thalamus'
for (st in c('Thalamus','SST','VIP','PVALB','LAMP5','PAX6')){
    sub.ststr = gsub("/","_",gsub(" ","_", st))
    subtype.diff.rda = paste0(srpref, 'Inh_', sub.ststr, '.rda')
    load(subtype.diff.rda)
    inh.regdf = processRegdf(est.regdf, 'inh')

    regdf = merge(inh.regdf, ast.regdf)
    cross.mat = table((regdf$inh_rank < cutoff),  (regdf$ast_rank < cutoff))
    ft = fisher.test(cross.mat)
    print(cross.mat)
    cat(st, '\tpval: ', ft$p.value, '\n')
    regdf$col = (regdf$inh_rank < cutoff) + (regdf$ast_rank < cutoff) 
    labdf = regdf[(regdf$inh_rank < cutoff) | (regdf$ast_rank < cutoff),]
    regdf = regdf[order(-regdf$col), ]
}



gp = ggplot(regdf, aes(inh_t, ast_t, col=factor(col))) + 
    geom_point() + 
    scale_color_manual(values=c('0'='grey90','1'='grey50','2'='red')) + 
    geom_text_repel(data=labdf, aes(inh_t, ast_t, label=symbol)) + 
    labs(x=paste0('Inhibitory (', st, ') t-statistic (dark = top 100 markers)'), y='Astrocyte t-statistic (dark = top 100 markers)') + 
    theme_pubr() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'ast_inh_', st, '_top', cutoff, 'markers_comparison')
saveGGplot(gp, pltprefix, w=8, h=6)



ggplot(regdf, aes(-log10(inh_p), -log10(ast_p))) + 
    geom_point() + 
    geom_text_repel(data=labdf, aes(-log10(inh_p), -log10(ast_p), label=symbol)) + 
    theme_pubr()



