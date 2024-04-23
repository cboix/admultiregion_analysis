#!/usr/bin/R
# --------------------------------------------------------
# Calculate whether modules at psbulk level are assoc. 
# with compositional changes at the cell / subtype level.
# NOTE: Particularly interested in Mic subtypes, OPC/Oli
# Nothing major - probably not the right way to go about asking this question.
# Updated 02/16/2021
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(emmeans)
options(width=175)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'module_vs_comp_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
# TODO: Opt to include neurons or not? Or just add later...
# commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Run the regressions on each compositional set:
# ----------------------------------------------
reg.tsv = paste0(moddir, 'modules_vs_composition_regresults.tsv')
reg.rds = paste0(moddir, 'modules_vs_composition_regresults.Rds')

if (!file.exists(reg.rds)){

    mingenes = 10
    # subset = 'Mic/Immune'
    runsets = c('Ast','Opc','Oli','Mic/Immune','Vasc/Epithelia','Inh', 'Exc')
    fulldf = c()
    for (subset in runsets){
        sstr = sub('/', '_', subset)
        print(sstr)

        # Load compositional data:
        commandArgs <- function(trailingOnly=TRUE){c(subset, TRUE)}
        source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))

        # Question: which modules predict which subtypes
        # also look only cross cell types?
        rmat = pivot.tomatrix(runscdf[runscdf$ng >= mingenes, c('pr','rm','score')], 'rm','score')
        ctdf$pr = with(ctdf, paste0(projid, '-', region))

        # Regression 
        lvls = sort(unique(ctdf$cls))
        ctdf$cls = factor(ctdf$cls, levels=lvls)
        subsetdf = c()
        gform = asform(c('cbind(count, other) ~ x * cls + cls * region * msex'))
        for (mod in colnames(rmat)){
            ctdf$x = rmat[ctdf$pr, mod]
            fit = try(glm(gform, ctdf, family='quasibinomial'))
            if (class(fit) != 'try-error'){
                cfit = coefficients(summary(fit))
                cfit = data.frame(cfit[grep("^x",rownames(cfit)),])
                names(cfit) = c('Est','SE','t','p')
                cat(mod, '\t', min(cfit$p), '\n')
                cfit$var = rownames(cfit)
                cfit$rm = mod
                subsetdf = rbind(subsetdf, cfit)
            }
        }
        subsetdf = subsetdf[order(subsetdf$p),]
        subsetdf$subset = sstr
        rownames(subsetdf) = NULL
        subsetdf$var[subsetdf$var == 'x'] = as.character(lvls[1])
        subsetdf$var = sub("^x:cls", "", subsetdf$var)
        fulldf = rbind(fulldf, subsetdf)
    }
    fulldf$mod.subset = sub("-.*", "", fulldf$rm)
    fulldf$in.subset = (fulldf$subset == fulldf$mod.subset)

    # Save results:
    write.table(fulldf, reg.tsv, quote=F, row.names=F, sep="\t")
    saveRDS(fulldf, file=reg.rds)
} else {
    fulldf = readRDS(reg.rds)
}


# Subset to not in same subtype, sort, adjust p-values:
# -----------------------------------------------------
subdf = fulldf[!fulldf$in.subset, ]
subdf = subdf[order(subdf$p),]
subdf = subdf[abs(subdf$Est) < 1e6,]
subdf$p.adj = p.adjust(subdf$p, 'BH')

head(subdf, 20)
head(subdf[!(subdf$subset %in% c('Exc', 'Inh', 'Vasc_Epithelia')),], 40)
# subdf = subdf[!(subdf$subset %in% c('Exc', 'Inh', 'Vasc_Epithelia')),]
# subdf = subdf[!(subdf$rm %in% c('Ast-12', 'Ast-24', 'Ast-19','Ast-7')),]


# 1. Astrocyte FOS + metallostasis are both linked to increased DUSP1 %
head(subdf[subdf$subset == 'Mic_Immune',])
# 2. 
head(subdf[subdf$subset == 'Exc' & subdf$var != 'Exc TOX3 INO80D',c('Est','p.adj','var','rm')], 40)
# -8.4795276 0.0000000000140122032426240            Exc TOX3 TTC6            Opc-27
# 0.2454900 0.0000000000735207797374313       Exc L3-4 RORB CUX2             Opc-5
# 0.3093390 0.0000000001897510276606802         Exc NXPH1 RNF220             Opc-5

head(subdf[subdf$subset == 'Oli',c('Est','p.adj','var','rm')], 40)
# 3. Oligodendrocytes:
# -0.85708629 0.000000396516729 Oli OPALIN             Ast-0
# -1.27053791 0.000014683570923 Oli OPALIN            Ast-10
# -1.04900543 0.000002971034266 Oli OPALIN Vasc_Epithelia-14 (may be contam SMC CD81)
# 1408 -0.44606411 0.000206172666681 Oli OPALIN     Mic_Immune-21 # TPT1 vs. OPALIN
# 1409 -0.45013931 0.000241623882211 Oli OPALIN      Mic_Immune-8

head(subdf[subdf$subset == 'Opc',c('Est','p.adj','var','rm')], 20)
head(subdf[subdf$subset == 'Ast',c('Est','p.adj','var','rm')], 20)



# Plot some of the top results as scatter plots:
# NOTE: Color by region to show if it's driven by region
# ------------------------------------------------------
# Load compositional data:
subset = 'Ast'
sstr = sub('/', '_', subset)
commandArgs <- function(trailingOnly=TRUE){c(subset, TRUE)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
ctdf$pr = with(ctdf, paste0(projid, '-', region))

# Question: which modules predict which subtypes
# also look only cross cell types?
mod = 'Oli-25'
ctdf$x = rmat[ctdf$pr, mod]
gp = ggplot(ctdf, aes(x, count / total, color=region)) + 
    facet_wrap(~cls, scales='free') + 
    geom_point(cex=.5) + 
    scale_color_manual(values=reg.cols) + 
    # geom_smooth(method='lm', color='black') + 
    labs(x=mod) + 
    geom_smooth(method='lm') + 
    stat_cor(label.sep="\n", output.type='text', color='black') + 
    theme_pubr() + theme(legend.position='none')
pltprefix = paste0(imgpref, sstr, '_', mod)
saveGGplot(gp, pltprefix, w=7.5, h=4)



