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
library(ggplot2)
library(ggpubr)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = plotdir
outpref = 'multiRegion/fracreg/'
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

asform = function(x){ as.formula(paste0(x, collapse='')) }

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
test.vars = c('nrad','cogad','nft','plaq_n','plaq_d', 'nft_regional','plaq_n_regional','plaq_d_regional','nft_zscore','nft_log1p')
ct.var = 'cell_type_high_resolution'
reg.vars = unique(cellmeta$region)
cell.vars = unique(cellmeta[[ct.var]])
optdf = expand.grid(cell=cell.vars, region=reg.vars, var=test.vars)

# Run on specific cases:
cdf = c()
for (i in 1:nrow(optdf)){
    cell = as.character(optdf$cell[i])
    test.var = as.character(optdf$var[i])
    region = as.character(optdf$region[i])
    cellstr = gsub("/","_", gsub(" ", "_", cell))
    outfile = paste0(outpref, 'glmm_logreg_', region, '_', test.var,'_', cellstr, '.tsv')
    if (file.exists(outfile)){
        df = read.delim(outfile, header=T, sep="\t", stringsAsFactors=F)
        cdf = rbind(cdf, df)
    }
}

# Calculate the overall fractions to decide which tests are kept out
# TODO: Plot left out ones as grey/NA
totdf = agg.rename(barcode ~ region, cellmeta, length, 'total')
ctdf = aggregate(barcode ~ cell_type_high_resolution + region, cellmeta, length)
ctdf = merge(ctdf, totdf)
ctdf$frac = ctdf$barcode / ctdf$total
mxdf = agg.rename(frac ~ cell_type_high_resolution, ctdf, max, 'max')
ctdf = merge(ctdf, mxdf)
ctdf$ratio = ctdf$max / ctdf$frac
ctdf = ctdf[ctdf$ratio < 10,]
ctdf = ctdf[,c('region','cell_type_high_resolution')]
names(ctdf) = c('reg','ct')


# Plot an example:
test.var = 'nft_regional'
for (tvar in test.vars){
    if (tvar %in% c('nrad','cogad')){
        test.var = paste0(tvar, 'TRUE')
    } else { test.var = tvar }
    ndf = cdf[cdf$coef == test.var,]
    ndf$log10p = -log10(ndf$p)
    ndf = ndf[ndf$est > -5,]
    ndf = ndf[order(ndf$p),]

    subdf = ndf[ndf$p < 0.05,]
    lims = c(-2.25,2.25)  /20
    mx = max(abs(subdf$est))
    lims = c(-mx, mx)
    gp = ggplot(subdf, aes(reg, ct, color=est)) + 
        geom_point(pch=15, size=10) + 
        scale_color_distiller(palette = 'RdBu', limits=lims) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'frac_sign_gg_', test.var, '.png'), gp, dpi=450, units='in', width=5, height=10)

    # Plot the figures in a boxplot manner:
    cts = as.character(unique(codf$major.celltype))
    countdf = aggregate(cell_type_high_resolution ~ major.celltype, codf, length)
    hs = c(1.5, countdf$cell_type_high_resolution)
    ws = 5
    nr = length(hs)
    nc = length(ws)

    if (sum(ndf$reg == 'TH') > 1){
        reglevels = c('EC','HC','TH','AG','MT','PFC')
    } else {
        reglevels = c('EC','HC','AG','MT','PFC')
    }
    ndf$reg = factor(ndf$reg, levels=reglevels)
    ndf = merge(ndf, ctdf)
    ndf$q = p.adjust(ndf$p, 'fdr')

    png(paste0(imgpref, 'frac_sign_heat_', test.var, '.png'), res=450, units='in', width=2, height=5)
    layout(matrix(1:(nr * nc), nrow=nr, ncol=nc), widths=hs, heights=hs)
    sp = 0.1; rsp=7;
    par(xaxs='i')
    par(yaxs='i')
    # Title here:
    par(mar=c(0, rsp, 0,0))
    image(as.matrix(1:length(reglevels)), col=reg.cols[reglevels], axes=F, useRaster=TRUE)
    text(y=parpos(2, -.5), x=parpos(1,.015),
         labels=tvar, xpd=TRUE, cex=1, font=2, adj=1)
    text(y=parpos(2, -.5), x=seq(0,1, length.out=length(reglevels)),
         labels=reglevels, xpd=TRUE, cex=.75, font=2, adj=.5)
    for (ct in cts) {
        subct = as.character(codf[codf$major.celltype == ct, 'cell_type_high_resolution'])
        subct = subct[subct %in% ndf$ct]
        subdf = ndf[ndf$ct %in% subct,]
        subdf$ct = factor(as.character(subdf$ct), levels=subct)
        subdf$log10p = -log10(subdf$p)
        # pwide = spread(subdf[,c('reg','ct','log10p')], reg, log10p, fill=0)
        rwide = spread(subdf[,c('reg','ct','est')], reg, est, fill=0)
        rmat = as.matrix(rwide[,-1])
        mx = max(abs(subdf$est))
        lims = c(-mx, mx)
        sigdf = subdf[subdf$p < 0.05,]
        # Image: 
        par(mar=c(sp, rsp, sp,sp))
        image(t(rmat), axes=F, useRaster=TRUE, col=rev(colrb), zlim=lims)
        box(lwd=.25)
        # Text:
        text(x=parpos(1, 0.015), y=seq(0,1, length.out=length(subct)),
             labels=subct, xpd=TRUE, cex=.5, adj=1)
        xat=seq(0,1, length.out=length(reglevels))
        yat=seq(0,1, length.out=length(subct))
        # Effect sizes 
        if (nrow(sigdf) > 0){
            within(sigdf, text(xat[as.numeric(reg)], yat[as.numeric(ct)], 
                               labels=ifelse(log10p > 3, '***', ifelse(log10p > 2,'**', '*')), xpd=TRUE))
        }
    }
    dev.off()

}


head(ndf[,c('est','reg','ct','log10p')], 20)


# TODO: Plot these coefficients, also remove ones with fractions significantly below other regions (region-specific ones)

# ---------------------------------------------
# Plot PCA of compositional data 
# both hcelltype and cell_type_high_resolution:
# ---------------------------------------------
library(compositions)
library(factoextra)

totdf = agg.rename(barcode ~ rind, cellmeta, length, 'total')
ctdf = aggregate(barcode ~ rind + cell_type_high_resolution, cellmeta, length)
hcdf = aggregate(barcode ~ rind + hcelltype, cellmeta, length)
ctdf = merge(ctdf, totdf)
hcdf = merge(hcdf, totdf)
ctdf$frac = ctdf$barcode / ctdf$total
hcdf$frac = hcdf$barcode / hcdf$total

# Matrices for PCA:
hwide = spread(hcdf[,c('hcelltype','rind','frac')], hcelltype, frac)
hmat = as.matrix(hwide[,-1])
rownames(hmat) = hwide[,1]
hcmat = as.matrix(clr(hmat))
groups = factor(sub("\\..*","", rownames(hcmat)), levels=names(reg.cols))

res.pca <- prcomp(hcmat, scale = TRUE)
png(paste0(imgpref, 'frac_pca_prop.png'), res=450, units='in', width=6, height=5.5)
fviz_pca_biplot(res.pca, repel = TRUE, pch=19,
                geom="point",
                col.ind = groups,
                palette = reg.cols,
                col.var = "#2E9FDF", # Variables color
                addEllipses = TRUE, # Concentration ellipses
                ellipse.type = "confidence"
)
dev.off()



# Matrices for PCA:
cwide = spread(ctdf[,c('cell_type_high_resolution','rind','frac')], cell_type_high_resolution, frac)
cmat = as.matrix(cwide[,-1])
rownames(cmat) = cwide[,1]
ccmat = as.matrix(clr(cmat))
groups = factor(sub("\\..*","", rownames(ccmat)), levels=names(reg.cols))

res.pca <- prcomp(ccmat, scale = TRUE)
png(paste0(imgpref, 'frac_pca_prop_highres.png'), res=450, units='in', width=6, height=5.5)
fviz_pca_ind(res.pca, repel = TRUE, pch=19,
             geom="point",
             col.ind = groups,
             palette = reg.cols,
             # col.var = "#2E9FDF", # Variables color
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence"
)
dev.off()

rownames(metadata) = metadata$rind 
nrc = metadata[rownames(ccmat), 'niareagansc']
nrad = rep('CTRL',length(nrc))
nrad[nrc <=2] = 'AD'

png(paste0(imgpref, 'frac_pca_prop_highres_ad.png'), res=450, units='in', width=6, height=5.5)
fviz_pca_ind(res.pca, repel = TRUE, pch=19,
             geom="point",
             col.ind = nrad,
             palette = c('indianred', 'royalblue'),
             # col.var = "#2E9FDF", # Variables color
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence"
)
dev.off()


rr = sweep(res.pca$rotation,2, res.pca$sdev, '*')
rdf = data.frame(region=groups, PC1 = res.pca$x[,'PC1'], PC2 =res.pca$x[,'PC2'])
adf = data.frame(ct=rownames(rr), PC1 = rr[,'PC1'], PC2 = rr[,'PC2'])

ggplot(rdf, aes(PC1, PC2, color=region)) + 
    geom_point() + 
    geom_segment(data=adf, aes(x=0, y=0, xend=PC1, yend=PC2), color='darkgrey') + 
    geom_text(data=adf, aes(PC1, PC2, label=ct), color='darkgrey') + 
    scale_color_manual(values=reg.cols) + 
    theme_pubr() + 
    theme(legend.position='right')








