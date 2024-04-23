#!/usr/bin/R
# -----------------------------------------------
# Plot some example genes with regional patterns:
# Updated 06/21/2023
# -----------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


library(tidyr)
library(viridis)
library(ggplot2)
library(ggpubr)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'gwas/')
imgpref = paste0(plotdir, 'gwas_')
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])
topdf = read.delim('Annotation/ADGWAS_topct_multiregion_121721.tsv', sep="\t")


# Load pseudobulk data for all runsets:
# -------------------------------------
reg.ps.rda = paste0(srdir, 'pseudobulk_data_all_regionaverages.rda')
load(reg.ps.rda)  # avg.mat


# Score genes by region specificity in CAMs / Mic:
# ------------------------------------------------
runset = 'Mic_Immune'
graph_id = 'boot'
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))

psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
load(psdata.rda)

cts = unique(ps.data$meta$cell_type_high_resolution)
kept.genes = intersect(gwgenes, rownames(ps.data$mat))

for (ct in cts){
    ind = ps.data$meta$cell_type_high_resolution == ct
    mdf = ps.data$meta[ind,]
    mat = ps.data$mat[kept.genes,ind]

    avdf = NULL
    resdf = NULL
    for (gene in kept.genes){
        mdf$val = mat[gene,]
        fit = lm(val ~ region, mdf)
        cfit = coefficients(summary(fit))
        cfit = data.frame(cfit)
        names(cfit) = c('Est','SE','t','p')
        cfit$var = rownames(cfit)
        cfit$gene = gene
        resdf = rbind(resdf, cfit)

        av = anova(fit)
        adf = as.data.frame(av)[1,,drop=FALSE]
        names(adf) = c('df','ss','ms','F','p')
        adf$gene = gene
        avdf = rbind(avdf, adf)
    }
    avdf = merge(avdf, topdf)
    avdf = avdf[order(avdf$p),]
    head(avdf, 20)
}


# Merge expression across microglia + macrophages:
# ------------------------------------------------
ind = ps.data$meta$cell_type_high_resolution != 'T cells'
ind = ps.data$meta$cell_type_high_resolution != 'T cells'
ps.data$meta$pr = with(ps.data$meta, paste0(projid, '_', region))
tform = make.tform(ps.data$meta$pr[ind])
tform = sweep(tform, 1, ps.data$meta$ncell[ind], '*')
tform = sweep(tform, 2, apply(tform, 2, sum), '/')

# Sample matrix + metadata:
samp.mat = ps.data$mat[,ind] %*% tform
samp.mdf = aggregate(ncell ~ pr + region + projid, ps.data$meta[ind,], sum)
samp.mat = samp.mat[, samp.mdf$pr]


# Score genes across this merged dataset:
# ---------------------------------------
mat = samp.mat[kept.genes,]
avdf = NULL
resdf = NULL
for (gene in kept.genes){
    samp.mdf$val = samp.mat[gene,]
    fit = lm(val ~ region, samp.mdf)
    cfit = coefficients(summary(fit))
    cfit = data.frame(cfit)
    names(cfit) = c('Est','SE','t','p')
    cfit$var = rownames(cfit)
    cfit$gene = gene
    resdf = rbind(resdf, cfit)

    av = anova(fit)
    adf = as.data.frame(av)[1,,drop=FALSE]
    names(adf) = c('df','ss','ms','F','p')
    adf$gene = gene
    avdf = rbind(avdf, adf)
}
avdf = merge(avdf, topdf)
avdf = avdf[order(avdf$p),]

mic.avdf = avdf[avdf$ct == 'Mic_Immune',]
mic.avdf$p.adj = p.adjust(mic.avdf$p, 'BH')
head(mic.avdf, 20)
hist(mic.avdf$p) # Looks ok.


# Plot these top prioritized genes:
# ---------------------------------
subdf = mic.avdf[mic.avdf$p.adj < 0.05,]
topgenes = subdf$gene
subdf$gene = factor(subdf$gene, levels=rev(topgenes))

# Which is the top region:
subresdf = resdf[resdf$gene %in% topgenes,]
subresdf = subresdf[subresdf$var != '(Intercept)',]
subresdf = merge(subresdf, aggregate(Est ~ gene, subresdf, max))
subresdf$region = ifelse(subresdf$Est < 0, 'AG', sub("region","", subresdf$var))
subdf = merge(subdf, subresdf[,c('region','gene')])
subdf$plab = with(subdf, ifelse(p.adj < 0.05, ifelse(p.adj < 0.01, ifelse(p.adj < 0.001, '***', '**'), '*'), '.'))
subdf$lab = with(subdf, paste(region, plab))

gp = ggplot(subdf, aes(gene, F, fill=region, label=lab)) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=reg.cols) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_text(aes(y=4), adj=0) + 
    labs(y='Anova F-statistic', x='GWAS gene') +
    theme_pubr() + coord_flip() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'ad_regional_genes_microglia')
saveGGplot(gp, pltprefix, w=2.25, h=2)


# Load the NIA-Reagan regional DEGs:
# ----------------------------------
path = 'nrad'
mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)

kept.cols = c('gene','col_nm','path','region')
set = 'Mic_Immune_Mic'
setdf = setdflist[[set]][, kept.cols]

regdf = setdf[setdf$gene %in% topgenes,]
regmat = pivot.tomatrix(regdf[,c('gene','col_nm','region')], 'region', 'col_nm')
regmat[is.na(regmat)] = 0
regmat = regmat[, reg.order[-1]]


# Plot the underlying expression patterns for these top genes:
# ------------------------------------------------------------
pltmat = avg.mat[topgenes, paste0("Mic_Immune-", reg.order[-1])]
pltmat = sweep(pltmat, 1, apply(pltmat, 1, max), '/')
pltmat = as.matrix(pltmat)
colnames(pltmat) = sub(".*-","",colnames(pltmat))

pmat = pltmat * 0
pmat[rownames(regmat), colnames(regmat)] = regmat

ux = 1.5
plt = Heatmap(pltmat,
              col=viridis(50),
              name='Expr.',
              use_raster=FALSE,
              cluster_columns=FALSE,
              cluster_rows=FALSE,
              width = ncol(pltmat)*unit(ux, "mm"), 
              height = nrow(pltmat)*unit(ux, "mm"),
              border_gp = gpar(col="black", lty = 1, lwd=.5),
              cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                  if (pmat[i,j] != 0){
                      grid.text('*', x, y,gp=gpar(fontsize=gridtxt.fs))
                  }
              })

h = 1 + 1 / 15 * nrow(pltmat)
w = 1.5 + 1 / 15 * ncol(pltmat)
pltprefix = paste0(imgpref, 'ad_regional_genes_microglia_heatmap')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Plot these genes as boxplots:
# -----------------------------
samp.mdf = merge(samp.mdf, unique(metadata[,c('projid','nrad','cogdxad','braaksc')]))
samp.mat = samp.mat[, samp.mdf$pr]

pltdf = NULL
pltgenes = c('APOE', 'SORL1', 'MS4A4A')
for (gene in pltgenes){
    samp.mdf$val = samp.mat[gene,]
    subdf = samp.mdf[samp.mdf$ncell >= 25,]
    subdf$gene = gene
    qt.th = quantile(subdf$val, .995)
    subdf$val[subdf$val > qt.th] = qt.th
    pltdf = rbind(pltdf, subdf)
}

# NOTE: Truncated outliers
pltdf$region = factor(pltdf$region, levels=reg.order[-1])

# TODO: ADD THE NIA-Reagan DEG scores:
nr.cols = c('CTRL'='grey80', 'AD'='indianred')
gp = ggplot(pltdf, aes(region, val, fill=nrad)) + 
    facet_wrap(~ gene, ncol=1, scales='free_y') + 
    scale_fill_manual(values=nr.cols) + 
    geom_boxplot(outlier.size=.5) + 
    labs(x='Region', y='Expr') +
    theme_pubr() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'ad_regional_genes_microglia_boxplots')
saveGGplot(gp, pltprefix, w=2.5, h=4)


