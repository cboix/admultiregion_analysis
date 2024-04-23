#!/usr/bin/R
# -----------------------------------------------
# Plot number of differential genes as a heatmap:
# Updated: 09/28/21
# -----------------------------------------------
# Aggregate number of sign. genes:
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
print(version)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# -------------------------------------------------------
# Load in the # significant tsv files across ALL regions:
# -------------------------------------------------------
# prefstr = paste0(celltype,'_',subtype)
method = 'nebula_ruv'
fnlist = list.files(path=regdir, pattern=paste0(method, '.*_allregions_.*.nsig.tsv'))

nsigdf = c(); alldf = c();
for (fn in fnlist){
    ndf = read.delim(paste0(regdir, fn), header=T)
    if (!('X1' %in% colnames(ndf))){ ndf$X1 = 0 }
    if (!('X2' %in% colnames(ndf))){ ndf$X2 = 0 }
    if (is.null(nsigdf)){
        nsigdf = ndf 
    } else {
        nsigdf = rbind(nsigdf, ndf[,colnames(nsigdf), drop=F])
    }
}
nsigdf$subtype[nsigdf$subtype == TRUE] = 'T cells'
names(nsigdf)[names(nsigdf) == 'X1'] = 'ndown'
names(nsigdf)[names(nsigdf) == 'X2'] = 'nup'
nsigdf$nup = as.numeric(nsigdf$nup)
nsigdf$ndown = as.numeric(nsigdf$ndown)
nsigdf = nsigdf[grep('Ast_', nsigdf$subtype, invert=T),]


# Plot number of significant genes:
# ---------------------------------
pcols = brewer.pal(12, 'Paired')
gplot = ggplot(nsigdf, aes(path, nup)) + 
    facet_wrap(subtype~region) + 
    geom_bar(position='dodge',stat='identity', fill=pcols[6]) + 
    geom_bar(data=nsigdf, aes(path, -ndown), position='dodge',stat='identity', fill=pcols[2]) + 
    labs(x='AD variable', y='Number of DEGs') + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, method, '_ndeg.allregions.png'), gplot, dpi=400, units='in', width=8,height=5)
ggsave(paste0(imgpref, method, '_ndeg.allregions.pdf'), gplot, dpi=400, units='in', width=8,height=5)


# Plot as a heatmap:
# ------------------
uwide = spread(nsigdf[,c('path','subtype','nup')], path, nup, fill=0)
dwide = spread(nsigdf[,c('path','subtype','ndown')], path, ndown, fill=0)

umat = as.matrix(uwide[,-1])
dmat = as.matrix(dwide[,-1])
rownames(umat) = uwide[,1]
rownames(dmat) = dwide[,1]
colnames(umat) = paste0(colnames(umat), '_up')
colnames(dmat) = paste0(colnames(dmat), '_down')


bmat = cbind(umat, -dmat)
mx = max(abs(bmat))
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
pdf(paste0(imgpref, method, '_ndeg.allregions.heatmap.pdf'), width=7, height=6)
Heatmap(bmat, name='Number\nof DEGs',
        use_raster=TRUE,
        col=col_fun,
        column_split=ifelse(1:ncol(bmat) %in% grep("_up", colnames(bmat)), 'Up','Down'),
        cell_fun = function(j, i, x, y, w, h, fill) {
            ann = abs(bmat[i,j])
            grid.text(ann, x, y)
        })
dev.off()



# Repeat, loading # DEGs for each of the regions:
# -----------------------------------------------
# prefstr = paste0(celltype,'_',subtype)
method = 'nebula_ruv'
fnlist = list.files(path=regdir, pattern=paste0(method, '.*.nsig.tsv'))

nsigdf = c(); alldf = c();
for (fn in fnlist){
    ndf = read.delim(paste0(regdir, fn), header=T)
    if (!('X1' %in% colnames(ndf))){ ndf$X1 = 0 }
    if (!('X2' %in% colnames(ndf))){ ndf$X2 = 0 }
    if (is.null(nsigdf)){
        nsigdf = ndf 
    } else {
        nsigdf = rbind(nsigdf, ndf[,colnames(nsigdf), drop=F])
    }
}
nsigdf$subtype[nsigdf$subtype == TRUE] = 'T cells'
names(nsigdf)[names(nsigdf) == 'X1'] = 'ndown'
names(nsigdf)[names(nsigdf) == 'X2'] = 'nup'
nsigdf$nup = as.numeric(nsigdf$nup)
nsigdf$ndown = as.numeric(nsigdf$ndown)
nsigdf = nsigdf[grep('Ast_', nsigdf$subtype, invert=T),]
nsigdf = nsigdf[grep('Exc_', nsigdf$subtype, invert=T),]
nsigdf = nsigdf[grep('^CA', nsigdf$subtype, invert=T),]
nsigdf = nsigdf[grep('^DG', nsigdf$subtype, invert=T),]
nsigdf = nsigdf[nsigdf$region != 'allregions',]
nsigdf = nsigdf[nsigdf$region != 'neocortex',]

# Plot as a heatmap:
# ------------------
uwide = spread(nsigdf[,c('path','region','subtype','nup')], path, nup, fill=0)
dwide = spread(nsigdf[,c('path','region','subtype','ndown')], path, ndown, fill=0)

umat = as.matrix(uwide[,-c(1:2)])
dmat = as.matrix(dwide[,-c(1:2)])
rownames(umat) = paste0(uwide[,1], '-', uwide[,2])
rownames(dmat) = paste0(dwide[,1], '-', dwide[,2])
colnames(umat) = paste0(colnames(umat), '_up')
colnames(dmat) = paste0(colnames(dmat), '_down')


bmat = cbind(umat, -dmat)
reg = sub("-.*", "", rownames(bmat))
rownames(bmat) = sub(".*-","", rownames(bmat))
mx = max(abs(bmat))
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
# pdf(paste0(imgpref, method, '_ndeg.allregions.heatmap.pdf'), width=7, height=6)
Heatmap(bmat, name='Number\nof DEGs',
        use_raster=TRUE,
        col=col_fun,
        row_split=reg,
        column_split=ifelse(1:ncol(bmat) %in% grep("_up", colnames(bmat)), 'Up','Down'),
        cell_fun = function(j, i, x, y, w, h, fill) {
            ann = abs(bmat[i,j])
            grid.text(ann, x, y)
        })
# dev.off()







# Annotate sides with number of cells

# Correlation + p-value from regression accounting for covariates
pdf(paste0(imgpref, 'correlation_regrpvals_',analysis,'heatmap.pdf'), width=7, height=6)
Heatmap(cr1mat,
        use_raster=TRUE,
        cluster_columns=FALSE, 
        cluster_rows=FALSE, 
        column_split=ifelse(colnames(cmat) %in% tophc, 'Hippocampus','Entorhinal Cortex'),
        row_split=ifelse(colnames(cmat) %in% tophc, 'Hippocampus','Entorhinal Cortex'),
        rect_gp = gpar(type = "none"), column_dend_side = "bottom",
        cell_fun = function(j, i, x, y, w, h, fill) {
            cond1 = (as.numeric(x) > 1 - as.numeric(y) + 1e-6)
            cond2 = (i <= length(topec) & j > length(topec))
            cond3 = (i > length(topec) & j <= length(topec))
            if((cond1 | cond2) & (!cond3)) {
                grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                if (i != j){
                    p = p2mat[i,j] # Adjusted p-values from regression.
                    ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                    grid.text(ann, x, y)
                }
            }
        })
dev.off()



# ---------------------------------
# Plot number of significant genes:
# ---------------------------------
pcols = brewer.pal(12, 'Paired')
gplot = ggplot(nsigdf, aes(path, nup)) + 
    facet_wrap(~region, ncol=1) + 
    geom_bar(position='dodge',stat='identity', fill=pcols[6]) + 
    geom_bar(data=nsigdf, aes(path, -ndown), position='dodge',stat='identity', fill=pcols[2]) + 
    labs(x='AD variable', y='Number of DEGs') + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, 'nb_ruv_ndeg.', prefstr,'.png'), gplot, dpi=400, units='in', width=2.5,height=5)
ggsave(paste0(imgpref, 'nb_ruv_ndeg.', prefstr,'.pdf'), gplot, dpi=400, units='in', width=2.5,height=5)


# --------------------------
# Plot a small volcano plot:
# --------------------------
for (path in unique(alldf$path)){
    pltdf = alldf[alldf$path == path,]
    pltdf = pltdf[order(pltdf$q),]
    pcols = brewer.pal(12,'Paired')
    ntop = sum(pltdf$padj < 0.05 & pltdf$col !=0, na.rm=T)
    labcut = 0.05
    print(paste("Number of sig. genes:", ntop))
    if (ntop > 50){ labcut = 0.02 }
    if (ntop > 100){
        NG = 70
        # Alternatively, plot top each:
        downdf = pltdf[(pltdf$padj < labcut) & pltdf$col == 1,]
        updf = pltdf[(pltdf$padj < labcut) & pltdf$col == 2,]
        labdf = rbind(head(downdf, NG/2), head(updf, NG/2))
    } else {
        labdf = pltdf[(pltdf$padj < labcut) & pltdf$col != 0,]
    }
    FCTHRESH = min(abs(pltdf$logFC[pltdf$col > 0]))

    gplot = ggplot(pltdf, aes(logFC , -log10(q), color=factor(col))) + 
        geom_point(cex=.1) + 
        geom_text_repel(data=labdf, aes(logFC, -log10(q), label=gene, color=factor(col)), size=2.5, segment.size=.25, max.overlaps=20) + 
        scale_color_manual(values=c('grey80',pcols[1],pcols[5])) + 
        scale_y_continuous(expand=c(0,0)) + 
        geom_vline(xintercept=0, lty='dashed') + 
        geom_vline(xintercept=c(-FCTHRESH, FCTHRESH), lty='dashed', col='grey50',lwd=.25) +
        theme_pubr() + theme(legend.position='none')
    ggsave(paste0(imgpref, 'volcano_',prefstr,'_', path, '_nebulaRUV.png'),gplot, units='in', dpi=450, width=4, height=5)
    ggsave(paste0(imgpref, 'volcano_',prefstr,'_', path, '_nebulaRUV.pdf'),gplot, units='in', dpi=450, width=4, height=5)
}

# -------------------------------------------------------------
# Get and plot overlaps between variables for the same regions:
# -------------------------------------------------------------
region = 'HC'
regs = unique(alldf$region)
for (region in regs){
    print(region)
    resdf = alldf[alldf$region == region,]
    # Sig. genes:
    sgenes = unique(resdf$gene[resdf$col != 0])
    resdf = resdf[resdf$gene %in% sgenes,]
    labdf = resdf[resdf$col != 0,]
    lwide = spread(labdf[,c('path','gene','logFC')], path, logFC, fill=0)
    lmat = as.matrix(lwide[,-1])
    rownames(lmat) = lwide[,1]
    smat = sweep(lmat,2,apply(abs(lmat), 2, max),'/')
    # Heatmap(smat)
    sind = which(apply(smat != 0, 1, sum) > 1)
    png(paste0(imgpref,'nb_ruv_heatmap_gt1.',prefstr,'.reg_', region ,'.png'), res=400, units='in',width=3, height=2 + 9 * length(sind) / 75)
    draw(Heatmap(smat[sind,], col=col_fun))
    dev.off()
}

# --------------------------------------------------------
# Get and plot overlaps between regions for the same vars:
# --------------------------------------------------------
vars = unique(alldf$path)
for (path in vars){
    print(path)
    resdf = alldf[alldf$path == path,]
    # Sig. genes:
    sgenes = unique(resdf$gene[resdf$col != 0])
    resdf = resdf[resdf$gene %in% sgenes,]
    labdf = resdf[resdf$col != 0,]
    lwide = spread(labdf[,c('region','gene','logFC')], region, logFC, fill=0)
    lmat = as.matrix(lwide[,-1])
    rownames(lmat) = lwide[,1]
    smat = sweep(lmat,2,apply(abs(lmat), 2, max),'/')
    # Heatmap(smat)
    sind = which(apply(smat != 0, 1, sum) > 1)
    png(paste0(imgpref,'nb_ruv_heatmap_gt1.',prefstr,'.path_', path,'.png'), res=400, units='in',width=3.5, height=2 + 9 * length(sind) / 75)
    draw(Heatmap(smat[sind,], col=col_fun))
    dev.off()
}


