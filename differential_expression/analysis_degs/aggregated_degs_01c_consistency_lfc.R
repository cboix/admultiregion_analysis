#!/usr/bin/R
# --------------------------------------------------------
# Plot the aggregated DE results across different methods:
# Consistency of DEG directionality by different variables
# Updated: 12/19/22
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
options(width=170)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Function for getting log fold-change vector:
getLFC = function(df, path){
    df = df[df$path == path,]
    lfc = df$logFC_nb
    names(lfc) = df$gr
    return(lfc)
}

# Function to process glm output:
process.fit = function(fit, set){
    cfit = data.frame(coefficients(summary(fit)))
    names(cfit) = c('Est','SE','t','p')
    cfit$var = rownames(cfit)
    cfit$set = set
    rownames(cfit) = NULL
    return(cfit)
}

# LFC comparison across runs:
# ---------------------------
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
sets = c('Ast_Ast','Exc_Exc','Inh_Inh','Opc_Opc','Mic_Immune_Mic','Oli_Oli')
allmat = c()
fitdf = c()
for (set in sets){
    print(set)
    fulldf = c()
    for (path in pathlist){
        mstr = paste0('allmethods.regional_', path)
        fullaggrda = paste0(regdir, mstr, '.merged.rda')
        load(fullaggrda)
        setdf = setdflist[[set]]
        kept.cols = c('gene','col_nm','path','region', 'logFC_nb')
        fulldf = rbind(fulldf, setdf[, kept.cols])
    }
    fulldf$gr = with(fulldf, paste0(gene, '_', region))
    fulldf$gs = paste0(fulldf$gene, '_', set)
    fulldf$pr = with(fulldf, paste0(path, ':', region))

    # Path-level similarity:
    NP = length(pathlist)
    rmat = matrix(0, nrow=NP, ncol=NP, dimnames=list(pathlist, pathlist))
    for (p1 in pathlist){
        l1 = getLFC(fulldf, p1)
        for (p2 in pathlist){
            l2 = getLFC(fulldf, p2)
            mn = names(l1)[names(l1) %in% names(l2)]
            rmat[p1, p2] = cor(l1[mn], l2[mn])
        }
    }

    # UMAP:
    fulldf$lfc = fulldf$logFC_nb
    ind = fulldf$path %in% c('nft','plaq_n','plaq_d')
    fulldf$lfc[ind] = fulldf$lfc[ind] * 25
    mat = pivot.tomatrix(fulldf[,c('gs','pr','lfc')], 'pr', 'lfc')
    dt = dist(t(mat), 'euclidean')

    nn = 5; mdist = 0.05;
    u = umap(dt, n_neighbors=nn, min_dist=mdist, verbose=F,
        repulsion_strength=.25)
    mdf = t(sapply(rownames(u), function(x){strsplit(x, ":")[[1]]}))
    udf <- data.frame(X1=u[,1],X2=u[,2], region=mdf[,2], path=mdf[,1], pr=rownames(u))

    gp = ggplot(udf, aes(X1, X2, label=path, color=region, pch=path)) + 
        geom_point() + 
        geom_text_repel() + 
        scale_color_manual(values=reg.cols) +
        theme_pubr()

    pltprefix = paste0(imgpref, 'allmethods_euclidean_lfc_umap_', set)
    saveGGplot(gp, pltprefix, h=5, w=6)

    # Perform regression on correlation:
    cmat = cor(mat, use='pairwise.complete.obs')
    ddf = as.data.frame(cmat)
    ddf$M1 = rownames(ddf)
    ddf = gather(ddf, M2, cor, -M1)
    ddf$M1 = factor(ddf$M1, levels=unique(ddf$M1))
    ddf$M2 = factor(ddf$M2, levels=unique(ddf$M1))
    ddf = ddf[as.numeric(ddf$M1) < as.numeric(ddf$M2),]
    # Separate region/path
    ddf$R1 = sub(".*:", "", ddf$M1)
    ddf$R2 = sub(".*:", "", ddf$M2)
    ddf$P1 = sub(":.*", "", ddf$M1)
    ddf$P2 = sub(":.*", "", ddf$M2)
    ddf$region = ddf$R1 == ddf$R2
    ddf$path   = ddf$P1 == ddf$P2

    # Not best family for cor but ok
    fit1 = glm(cor ~ region + path, ddf, family='gaussian')
    fit2 = glm(cor ~ region*R1 + path*P1, ddf, family='gaussian')
    cfit1 = process.fit(fit1, set)
    cfit2 = process.fit(fit2, set)
    cfit1$run = 'simple'
    cfit2$run = 'interact'
    fitdf = rbind(fitdf, rbind(cfit1, cfit2))

    if (is.null(allmat)){
        allmat = mat
    } else {
        allmat = rbind(allmat, mat[, colnames(allmat)])
    }
}


# Plot the coefficients for regression on logFC_nb correlations:
# --------------------------------------------------------------
sdf = fitdf[fitdf$run == 'simple',]

# Plot these statistics as a barplot:
gp = ggplot(sdf, aes(var, Est, fill=set)) + 
    geom_bar(stat='identity', position='dodge', color=NA) + 
    scale_y_continuous(expand=c(0,0))+
    labs(y='Estimated coeff. for pred. correlation', x='Variable') +
    theme_pubr()

pltprefix = paste0(imgpref, 'allmethods_regr_corr_allmajor')
saveGGplot(gp, pltprefix, h=4, w=6)


# Plot interaction model statistics as a barplot:
sdf = fitdf[fitdf$run == 'interact',]
gp = ggplot(sdf, aes(var, Est, fill=set)) + 
    geom_bar(stat='identity', position='dodge', color=NA) + 
    scale_y_continuous(expand=c(0,0))+
    labs(y='Estimated coeff. for pred. correlation', x='Variable') +
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

pltprefix = paste0(imgpref, 'allmethods_regr_corr_interaction_allmajor')
saveGGplot(gp, pltprefix, h=5, w=10)



# Plot UMAP for all cell types jointly:
# -------------------------------------
dt = dist(t(allmat), 'euclidean')
nn = 5; mdist = 0.05;
u = umap(dt, n_neighbors=nn, min_dist=mdist, verbose=F,
    repulsion_strength=.25)
mdf = t(sapply(rownames(u), function(x){strsplit(x, ":")[[1]]}))
udf <- data.frame(X1=u[,1],X2=u[,2], region=mdf[,2], path=mdf[,1], pr=rownames(u))

gp = ggplot(udf, aes(X1, X2, label=path, color=region, pch=path)) + 
    geom_point() + 
    geom_text_repel() + 
    scale_color_manual(values=reg.cols) +
    theme_pubr()

pltprefix = paste0(imgpref, 'allmethods_euclidean_lfc_umap_allmajor')
saveGGplot(gp, pltprefix, h=5, w=6)


