#!/usr/bin/R
# --------------------------------------------------------------------------
# Calculate the fraction differences, by logistic regression + mixed effects
# Updated 10/23/2020 
# --------------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(lme4)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = plotdir
outpref = 'multiRegion/fracreg/'
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

asform = function(x){ as.formula(paste0(x, collapse='')) }

chunksize = 4 # 570 chunks
chunk = as.integer(Sys.getenv('SGE_TASK_ID'))

agg.rename = function(formula, data, FUN, name){
    agg.df = aggregate(formula, data, FUN)
    names(agg.df)[ncol(agg.df)] = name
    return(agg.df)
}

# --------------------------------------
# Load in the final metadata (cellmeta):
# --------------------------------------
load(file=paste0(datadir, prefix, '.final_noMB.cell_labels.Rda'))
# Colors:
load('Annotation/multiregion_celltypes_colors.Rda')
type.cols = tcols
type.cols = c(type.cols, major.col['Inh'], major.col['Exc'])
tsp.type.cols = sapply(type.cols, tsp.col)

# Order for full.exttype:
odf = unique(cellmeta[,c('major.celltype','full.exttype')])
odf = odf[order(odf$full.exttype),]
odf$major.celltype = factor(odf$major.celltype, levels=c('Ast','Mic/Immune','Oli','Opc','Vasc/Epithelia', 'Exc','Inh'))
odf = odf[order(odf$major.celltype),]
odf$full.exttype = factor(odf$full.exttype, levels=odf$full.exttype)

# Order for cell_type_high_resolution:
codf = unique(cellmeta[,c('major.celltype','cell_type_high_resolution')])
codf = codf[order(codf$cell_type_high_resolution),]
codf$major.celltype = factor(codf$major.celltype, levels=c('Ast','Mic/Immune','Oli','Opc','Vasc/Epithelia', 'Exc','Inh'))
codf = codf[order(codf$major.celltype),]
codf$cell_type_high_resolution = factor(codf$cell_type_high_resolution, levels=codf$cell_type_high_resolution)

# Get the per-region pathology scores:
pqdf = read.delim(paste0(datadir, 'region_pathology_scores.tsv'), header=T) 
colnames(pqdf)[3:5] = paste0(colnames(pqdf)[3:5],"_regional")
# Make nft diff vars:
pqdf$nft_log1p = log(pqdf$nft_regional + 1)
sdf = merge(agg.rename(nft_regional ~ region, pqdf, sd, 'nft_sd'),
            agg.rename(nft_regional ~ region, pqdf, mean, 'nft_mean'))
pqdf = merge(pqdf, sdf, all.x=TRUE)
pqdf$nft_zscore = (pqdf$nft_regional - pqdf$nft_mean) / pqdf$nft_sd

metadata = merge(metadata, pqdf, all.x=TRUE)
metadata$nrad = metadata$niareagansc <= 2
metadata$cogad = metadata$cogdx > 3

# Create all test cases:
test.vars = c('nrad','cogad','nft','plaq_n','plaq_d', 'nft_regional','plaq_n_regional','plaq_d_regional', 'nft_zscore','nft_log1p')
ct.var = 'cell_type_high_resolution'
reg.vars = unique(cellmeta$region)
cell.vars = unique(cellmeta[[ct.var]])
optdf = expand.grid(cell=cell.vars, region=reg.vars, var=test.vars)

# Run on specific cases:
ind = (chunksize * (chunk - 1) + 1):min(c(nrow(optdf), chunk * chunksize))
for (i in ind){
    cell = as.character(optdf$cell[i])
    test.var = as.character(optdf$var[i])
    region = as.character(optdf$region[i])
    cellstr = gsub("/","_", gsub(" ", "_", cell))
    outfile = paste0(outpref, 'glmm_logreg_', region, '_', test.var,'_', cellstr, '.tsv')
    if (!file.exists(outfile)){
        cat("[STATUS] Running:", '\t', region, '\t', test.var, '\t', cell, '\t')
        # Set up testing data:
        subdf = cellmeta[cellmeta$region == region,]
        subdf = merge(subdf, metadata[,c('rind',test.var, 'msex','age_death','pmi', 'projid')])
        subdf$is.ct = subdf[[ct.var]] == cell
        # Run logistic regression with mixed effects:
        fit = try(glmer(asform(c('is.ct ~',test.var,'+ msex + age_death + pmi + (1|rind)')), subdf, family=binomial))
        if (class(fit) != 'try-error'){
            cdf = data.frame(coefficients(summary(fit)))
            names(cdf) = c('est','se','z','p')
            cdf$reg = region
            cdf$coef = rownames(cdf)
            cdf$ct = cell
            cdf$var = test.var
            if (test.var %in% c('nrad','cogad')){
                line = cdf[cdf$coef == paste0(test.var,'TRUE'),]
            } else {
                line = cdf[cdf$coef == test.var,]
            }
            cat(round(line$est,2), '\t', round(line$p,2), '\n')
            write.table(cdf, file=outfile, row.names=F, col.names=T, quote=F, sep="\t")
        } else {
            cat("NA\tNA\n")
        }
    } else { 
        cat("[STATUS] Skipping ", i, ", finished already\n")
    }
} 

