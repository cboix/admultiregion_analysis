#!/usr/bin/R
# ----------------------------------
# Starting from all DEG results,
# - Build region-specificity scores
# - Regression of lFC across results
# NOTE: mostly preliminary, not finished analysis
# Updated: 03/02/23
# ----------------------------------
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
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))
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
    'Down (1-2)'=col.paired[2],
    'Down (3+)'=col.paired[1],
    'Down (multi-CT)'='slateblue4',
    'Up (1-2)'=col.paired[6],
    'Up (3+)'=col.paired[5],
    'Up (multi-CT)'='brown4')

path = 'nrad'
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)


    # Calculate the set shared in each region:
    # ----------------------------------------
    kept.cols = c('gene','col_nm','path','region', 'logFC_nb', 'p_nb')
    alldf = c()
    for (set in keep.sets){
        print(set)
        setdf = setdflist[[set]][, kept.cols]
        # setdf = setdf[setdf$col_nm != 0,]
        setdf$set = set
        alldf = rbind(alldf, setdf)
    }

    set = 'Mic_Immune_Mic'
    ps.data = load_psbulk_matrix(set)
    fulldf = alldf[alldf$set == set,]
    degs = unique(fulldf$gene[fulldf$col_nm != 0])
    covars = c('msex','Apoe_e4','pmi','age_death')
    ps.data$meta = merge(ps.data$meta, unique(metadata[metadata$region == 'PFC',c('projid', 'nrad', covars)]))
    ps.data$mat = ps.data$mat[,ps.data$meta$pr]
    mat = as.matrix(ps.data$mat[degs,])


    # TODO: Redo this scoring, with DESeq2?

    # Test genes:
    resdf = c()
    ndf = c()
    pb = txtProgressBar(min = 0, max = length(degs), initial = 0, style=3) 
    for (i in 1:length(degs)){
        setTxtProgressBar(pb,i)
        gene = degs[i]
        ps.data$meta$val = mat[gene,]
        fit = glm(val ~ nrad * region + msex + pmi + age_death, ps.data$meta, 
            family='gaussian', weights=log(ps.data$meta$ncell))
        cfit = as.data.frame(coefficients(summary(fit)))
        names(cfit) = c('Est','SE','t','p')
        cfit$var = rownames(cfit)
        cfit$gene = gene
        int = cfit['nradAD', 'Est']
        ind = grep('nradAD:', cfit$var)
        df = data.frame(var=c('nrad:regionAG', cfit[ind,'var']), Est=c(int, int + cfit[ind,'Est']), gene=gene)
        resdf = rbind(resdf, cfit)
        ndf = rbind(ndf, df)
    }
    close(pb)


    # Average expression of these genes:
    vdf = data.frame(avgexpr=rowMeans(mat))
    vdf$gene = rownames(vdf)

    # Statistic on regression estimates:
    sdf = agg.rename(Est ~ gene, ndf, sd, 'sd.int.Est')
    sdf = merge(sdf, agg.rename(Est ~ gene, ndf, function(x){sd(x)/mean(x)}, 'cov.int.Est')) # coefficient of variation:
    sdf = merge(sdf, agg.rename(Est ~ gene, ndf, function(x){max(abs(x))}, 'max.int.Est'))
    sdf = merge(sdf, vdf)
    sdf$stat = abs(sdf$sd.int.Est / sdf$avgexpr)
    # sdf$stat = abs(sdf$cov.int.Est / sdf$avgexpr)
    sdf = sdf[order(sdf$stat, decreasing=T),]

    head(sdf, 20) # Is pulling out most significant genes, should control for this?? also PCDH9/PLP1
    tail(sdf)
    # Relative to their intercept?

    sdf[sdf$gene == 'APOC1',]
    ndf[ndf$gene == 'APP',] # TH?









    sdf = merge(sdf, agg.rename(Est ~ gene, resdf[grep('nradAD:', resdf$var),], max, 'max.int.Est'))
    # Compare the nradAD coef to each of the regions 
    # + to variance of coeffs for each region 
    sdf = agg.rename(Est ~ gene, resdf[grep('nradAD', resdf$var),], sd, 'sd.int.Est')
    sdf = merge(sdf, agg.rename(Est ~ gene, resdf[grep('nradAD:', resdf$var),], max, 'max.int.Est'))
    sdf = merge(sdf, resdf[resdf$var == 'nradAD', c('Est','SE','p','gene')])
    sdf$stat = log2(abs(sdf$sd.int.Est / sdf$Est))
    
    sdf = sdf[order(sdf$stat, decreasing=T),]
    topdf = sdf[abs(sdf$Est + sdf$max.int.Est) > 0.2,]


    topdf = topdf[order(topdf$sd.int.Est, decreasing=T),]

    gene = 'SKAP2'
    sdf[sdf$gene == gene,]
    resdf[resdf$gene == gene,]



    subdf = fulldf[fulldf$gene == gene,]
    # fit = glm(logFC_nb ~ 0 + region * set, subdf, family='gaussian')
    f1 = glm(logFC_nb ~ 1, subdf, family='gaussian')
    f2 = glm(logFC_nb ~ 1 + region, subdf, family='gaussian')
    anova(f1, f2)

    # df = expand.grid(region=unique(subdf$region), set=unique(subdf$set))
    # df$pred = predict(fit, df)


    # TODO: Load in the original naming splits
    cdf = agg.rename(path ~ gene + col_nm + region, alldf, length, 'count')
    cdf = cdf[cdf$count >= 3, ]
    # head(cdf[order(cdf$count, decreasing=T),], 50)

}
