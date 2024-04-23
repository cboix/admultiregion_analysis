#!/usr/bin/R
# --------------------------------------------------------------
# Nebula as in minimal example, but using the precomp. matrices:
# Updated: 05/25/21
# --------------------------------------------------------------
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

celltype = 'Oli'
subtype = ''

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

# -----------------------------
# Load in the regression files:
# -----------------------------
prefstr = paste0(celltype,'_',subtype)
fnlist = list.files(path=regdir, pattern=paste0('nebula_ruv.',prefstr, '.*_allregions_.*rda'))
nsigdf = c(); alldf = c();
for (fn in fnlist){
    load(paste0(regdir, fn))
    region = nsig[1,'region']
    path = nsig[1,'path']
    resdf$region = region
    resdf$path = path
    ndf = data.frame(nsig)
    if (!('X1' %in% colnames(ndf))){ ndf$X1 = 0 }
    if (!('X2' %in% colnames(ndf))){ ndf$X2 = 0 }
    if (is.null(nsigdf)){
        nsigdf = ndf 
    } else {
        nsigdf = rbind(nsigdf, ndf[,colnames(nsigdf), drop=F])
    }
    fulldf$path = path
    fulldf$region = region
    fulldf = fulldf[order(fulldf$col != 0, abs(fulldf$logFC), decreasing=T),]
    fulldf$rank = 1:nrow(fulldf)
    alldf = rbind(alldf, fulldf)
}
names(nsigdf)[names(nsigdf) == 'X1'] = 'ndown'
names(nsigdf)[names(nsigdf) == 'X2'] = 'nup'
nsigdf$nup = as.numeric(nsigdf$nup)
nsigdf$ndown = as.numeric(nsigdf$ndown)

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


