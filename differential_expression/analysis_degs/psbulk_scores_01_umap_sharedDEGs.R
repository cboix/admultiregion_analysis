#!/usr/bin/R
# ----------------------------------------------------------
# Score all samples (psbulk) based on shared degs (01e)
# Updated: 01/09/23
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(uwot)

library(ggrepel)
library(ggplot2)
library(ggpubr)
print(version)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'sampscore_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


# Load regional differential results, process to sets for each cell type:
# -----------------------------------------------------------------------
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
sets = c('Ast_Ast','Exc_Exc','Inh_Inh','Opc_Opc','Mic_Immune_Mic','Oli_Oli')
pctdf = c()
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    shared.degs.rda = paste0(regdir, mstr, '.sharedDElists.rda')
    load(shared.degs.rda) # uplist, dwlist

    # Load pseudobulk data:
    # set = 'Ast_Ast'
    ststr = sub('_[A-Za-z]+$', '', set)
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
    load(psdata.rda)

    # Merge to sample-level (individual x region):
    ps.data = aggregate_psbulk_samplelevel(ps.data)

    # Subset to shared DEGs 
    shared.degs = c(uplist[[set]], dwlist[[set]])
    demat = ps.data$mat[shared.degs,]

    # UMAP from demat:
    dt = dist(t(demat))  # Euclidean distance

    nn = 10; mdist = 0.05; set.seed(0);
    u = umap(dt, n_neighbors=nn, min_dist=mdist, verbose=F,
        repulsion_strength=.25)
    mdf = t(sapply(rownames(u), function(x){strsplit(x, "_")[[1]]}))
    udf <- data.frame(X1=u[,1],X2=u[,2], region=mdf[,2], projid=mdf[,1], pr=rownames(u))
    udf = merge(udf, aggregate(gpath ~ projid + niareagansc, metadata, mean))
    pqdf = merge(pqdf, unique(metadata[,c('rind','projid')]))
    udf = merge(udf, pqdf, all.x=TRUE)

    gp = ggplot(udf, aes(X1, X2, color=region, fill=log1p(nft))) + 
        geom_point(pch=21) + 
        # scale_color_manual(values=colvals$niareagansc) +
        scale_color_manual(values=reg.cols) +
        scale_fill_viridis(na.value='white') +
        theme_pubr()

    pltprefix = paste0(imgpref, 'euclidean_umap_sharedgenes_', set, '_', path)
    saveGGplot(gp, pltprefix, h=5, w=6)


    # Summary scores from this:
    upmat = ps.data$mat[uplist[[set]],]
    dwmat = ps.data$mat[dwlist[[set]],]
    sampscore = (apply(upmat, 2, sum) - apply(dwmat, 2, sum)) / length(shared.degs)

    df = data.frame(score=sampscore, pr=names(sampscore))
    df = merge(df, udf, all.x=TRUE)
    df = df[order(df$score),]
    df$pr = factor(df$pr, levels=df$pr)
    df$path = log1p(df[[path]])

    # gp = ggplot(df, aes(pr, score, color=region, fill=log1p(nft))) + 
    gp = ggplot(df, aes(pr, score, fill=path)) + 
        geom_bar(stat='identity') + 
        scale_y_continuous(expand=c(0,0)) + 
        # scale_color_manual(values=reg.cols) +
        scale_fill_viridis(na.value='white', name=paste0('log1p(', path,')')) +
        theme_pubr() + theme(legend.position=c(.1,.6)) +
        theme(axis.text.x=element_blank())

    pltprefix = paste0(imgpref, 'barplot_score_sharedgenes_', set, '_', path)
    saveGGplot(gp, pltprefix, h=3, w=10)


    gp = ggplot(df, aes(score, path, color=region)) + 
        geom_smooth(method='gam') + 
        geom_point() + 
        scale_color_manual(values=reg.cols) +
        labs(x='Score', y=path) + 
        theme_pubr()
    pltprefix = paste0(imgpref, 'scatter_score_sharedgenes_', set, '_', path)
    saveGGplot(gp, pltprefix, h=5, w=6)



}


