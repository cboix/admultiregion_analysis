#!/usr/bin/R
# -----------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Last updated mid-2020 (exploratory analysis)
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(rhdf5)
library(glmnet)
library(lme4)
library(cba)

# For model plots:
library(sjPlot)
library(sjmisc)
library(patchwork)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/markers/')
imgpref = paste0(plotdir, 'pred_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# ------------------------
# Load pathology measures:
# ------------------------
pathfile = paste0(datadir, prefix, '.nft.tsv.gz')
nfdf = read.delim(pathfile, header=F)
names(nfdf) = c('barcode','nft')
pathfile = paste0(datadir, prefix, '.plaq_n.tsv.gz')
pqdf = read.delim(pathfile, header=F)
names(pqdf) = c('barcode','plaq_n')
pathfile = paste0(datadir, prefix, '.plaq_d.tsv.gz')
pq2df = read.delim(pathfile, header=F)
names(pq2df) = c('barcode','plaq_d')
pqdf = merge(pqdf, pq2df)
pqdf = merge(pqdf, nfdf)
print(dim(pqdf))


# -------------------------------
# Load in and plot for each cell:
# -------------------------------
clist = c('Astro','Oligo','OPC', 'Microglia', 'Endo', 'Per')
celltype = 'Oligo'
celltype = 'Microglia'
path = 'nft' # TODO: Generalize to other pathology
# path = 'plaq_n'
for (celltype in clist){
    print(celltype)
    ctpref = paste0(imgpref, celltype, '_', path, '_')
    # Read in pathology and coefficients:
    if (path != 'nft'){ celltype = paste0(path, '_', celltype) }
    lblset = 'leiden_r5_n50'
    coeffdf = read.delim(paste0(datadir, prefix, '.', lblset, '.coeff_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
    corrdf = read.delim(paste0(datadir, prefix, '.', lblset, '.corr_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
    pathdf = read.delim(paste0(datadir, prefix, '.', lblset, '.pathpred_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
    h5file = paste0(datadir, prefix, '.', lblset, '.topgenes_', celltype, '.hdf5')
    pathdf$region = sub('_.*','',pathdf$barcode)
    pathdf = merge(pathdf, celldf, all.x=TRUE)
    lcdf = coeffdf[coeffdf$fit == 'Lasso',]
    lcdf = lcdf[abs(lcdf$coeff) > 0.1,]
    pathdf = merge(pathdf, pqdf)

    # Get the expression matrix for top genes (testing lmm models):
    # kept.genes = head(lcdf$symbol,5)
    # Slice the matrix for these genes:
    h5attr = h5ls(h5file)
    ngenes = as.numeric(h5attr[5][h5attr[2] == 'genes'])
    chunk = 500
    if (ngenes > 100){
    h5f = H5Fopen(h5file)
    genes = h5f$genes
    bcs = h5f$barcodes
    # Open handle, extract genes we care about and close:
    h5d = h5f&"matrix"
    mat = t(h5d[])
    H5Dclose(h5d)
    H5Fclose(h5f)
    colnames(mat) = genes
    rownames(mat) = bcs

    # Order as pathdf:
    mat = mat[pathdf$barcode,]

    # Reorder the genes:
    dt = dist(t(mat))
    ht <- hclust(dt, method = 'ward.D')
    cocl <- order.optimal(dt, ht$merge)$order
    reord <- names(cocl)[cocl]

    scmat = scale(mat)

    # Visualize these top genes against pathology variables:
    # Add region 
    # Add APOE status
    ord = order(pathdf$path) 
    png(paste0(ctpref, 'topgenes_heatmap_pathord.png'), units='in', res=450, width=8, height=4.5)
    sp = 0.1
    par(mar=c(sp,2,sp,sp))
    image(scmat[ord,reord], axes=F, col=viridis(100))
    text(y=seq(0,1,length.out=ncol(mat)), x=parpos(1, 0.002),
         labels=reord, xpd=TRUE, srt=0, cex=.3, adj=1)
    dev.off()

    # Plot genes:
    df = data.frame(mat)
    df$barcode = rownames(df)
    df = gather(df, symbol, value, - barcode)
    df = merge(df, pathdf[,c('barcode','region','projid', 'plaq_n', 'nft','plaq_d')])
    df$projid = factor(df$projid)

    asform = function(x){ as.formula(paste0(x, collapse='')) }

    # gene = 'HLA.DRB1'
    # gene = 'APOE'
    gene = 'PRKCB'
    subdf = df[df$symbol == gene,]
    m0 = glm(asform(c(path, '~ region')), data=subdf)
    m1 = glm(asform(c(path, '~ value + region')), data=subdf)
    m2 = glm(asform(c(path, '~ value + region + value*region')), data=subdf)
    av = anova(m0, m1, m2)
    pchisq(av$Deviance, df=av$Df, lower.tail=FALSE)

    plot.gene=FALSE
    gene = 'PRKCB'
    df1 = c()
    df2 = c()
    dfp = c()
    regpdf = c()
    use.re=FALSE
    if (use.re) { rpref = paste0(ctpref, 'withrandeff_') } else { rpref = ctpref } 
    pathlist = c('nft','plaq_d','plaq_n')
    for (gene in genes){
        print(gene)
        gene = sub("-","\\.", gene)
        subdf = df[df$symbol == gene,]
        cl1 = list()
        cl2 = list()
        pvals = matrix(0, nrow=2, ncol=3, dimnames=list(NULL,pathlist))
        for (path in pathlist){
            if (use.re){
                l0 = lmer(asform(c(path, '~ region + (1|projid)')), data=subdf)
                l1 = lmer(asform(c(path, '~ value + region + (1|projid)')), data=subdf)
                l2 = lmer(asform(c(path, '~ region + value:region + (1|projid)')), data=subdf)
            } else {
                l0 = lm(asform(c(path, '~ region')), data=subdf)
                l1 = lm(asform(c(path, '~ value + region')), data=subdf)
                l2 = lm(asform(c(path, '~ region + value:region')), data=subdf)
            }
            av = anova(l0, l1, l2)
            # Plot interactions terms visually:
            if (plot.gene){
                g1 = plot_model(l2, title=paste('Modeling', path, 'by', gene))
                if (use.re){
                    g2 = plot_model(l2, type = "int", title=paste('Prediction of', path, 'from', gene))
                } else {
                    g2 = plot_model(l2, type = "pred", title=paste('Prediction of', path, 'from', gene), terms=c('value','region'))
                }
                g3 = ggplot(subdf, aes_string('value', path)) + facet_wrap(~region) + 
                    geom_point(alpha=0.5) + geom_smooth(method='lm')
                subdf$vbins = cut(subdf$value, c(0,1,2,3,4,6), include.lowest=T)
                g3 = ggplot(subdf, aes_string('vbins', path)) + 
                    facet_wrap(~region, scales='free_y') + geom_boxplot()
                gplot = g1 + g2
                ggsave(paste0(rpref, 'gene_', gene, '_', path, '_coeff_int_plot.png'), gplot,
                       dpi=450, units='in', width=7, height=4)
            }
            # Save coeff:
            cl1[[path]] = coefficients(l1)
            cl2[[path]] = coefficients(l2)
            # Extended coefficients for interaction term:
            cvals = summary(l2)$coefficients
            cvals = as.data.frame(cvals[grep('value', rownames(cvals)),])
            cvals$gene = gene
            cvals$path = path
            cvals$region = sub(':value', '', sub('region', '', rownames(cvals)))
            regpdf = rbind(regpdf, cvals)
            # P-values in general:
            df1 = rbind(df1, cbind(gene=gene, path=path, cl1[[path]]))
            df2 = rbind(df2, cbind(gene=gene, path=path, cl2[[path]]))
            dfp = rbind(dfp, data.frame(gene=gene, path=path, p1=av$Pr[2], p2=av$Pr[3]))
            pvals[,path] = av$Pr[-1]
        }
    }
    regpdf$log10p = -log10(regpdf[,'Pr(>|t|)'])


    dfp[order(dfp[,3]),]
    dfp[order(dfp[,4]),]

    df1 = dfp[,c('gene','path','p1')]
    df2 = dfp[,c('gene','path','p2')]
    names(df2)[3] = 'p1'
    df2$path = paste0(df2$path, '\n(int)')
    dfboth = rbind(df1, df2)
    pwide = spread(dfboth, path, p1)
    pmat = as.matrix(pwide[,-1])
    rownames(pmat) = pwide$gene
    pmat = reord(pmat)
    pmat = pmat[order(pmat[, 'nft'], decreasing=T),]
    pmat = -log10(pmat)
    pmat = t(pmat)

    scale = 0.6
    h = 0.5 + ncol(pmat) * scale * .3
    w = 0.5 + nrow(pmat) * scale
    png(paste0(ctpref, 'overall_lm_pvals.png'), units='in', res=450, width=w, height=h)
    sp = 0.1
    par(mar=c(sp, 6,3,sp))
    mxr = 50
    plt.mat = pmat
    plt.mat[plt.mat > mxr] = mxr
    image(plt.mat, col=colb, zlim=c(0, mxr), axes=F)
    text(x=parpos(1, .02),
         y=seq(0,1,length.out=ncol(plt.mat)),
         labels=colnames(plt.mat), xpd=TRUE, adj=1, cex=.9)
    text(y=parpos(2, -1.012), x=seq(0,1,length.out=nrow(plt.mat)),
             labels=rownames(plt.mat), xpd=TRUE, adj=c(.5,0))
    ind = which(plt.mat != 0,arr.ind=T)
    vals = round(pmat[ind])
    xat = seq(0,1,length.out=nrow(plt.mat))
    yat = seq(0,1,length.out=ncol(plt.mat))
    text(x=xat[ind[,1]], y=yat[ind[,2]],
         labels=vals, xpd=TRUE, cex=.8, col=ifelse(vals > 35,'grey90','black'))
    box(lwd=.5)
    dev.off()


    # ---------------------------------------------------------------------
    # Plot the region-specific coefficients and significance for each gene:
    # Plot as a heatmap, cluster by which genes show different effects by region.
    # ---------------------------------------------------------------------
    path = 'nft'
    pmlist = list()
    emlist = list()
    for (path in pathlist){
        srdf = regpdf[regpdf$path == path, c('gene','log10p', 'Estimate', 'region')]
        ewide = spread(srdf[,c('gene','Estimate','region')], region, Estimate)
        pwide = spread(srdf[,c('gene','log10p','region')], region, log10p)
        emat = as.matrix(ewide[,-1])
        pmat = as.matrix(pwide[,-1])
        rownames(emat) = ewide$gene
        rownames(pmat) = pwide$gene
        regst = c('EC','HC','MT','AG','PFC')
        emat = reord(emat)[,regst]
        pmat = pmat[rownames(emat),regst]
        emat = t(emat)
        pmat = t(pmat)

        scale = 0.6
        h = 0.5 + ncol(emat) * scale * .3
        w = 0.5 + nrow(emat) * scale
        png(paste0(ctpref, 'overall_', path, '_perregion_estimates.png'), units='in', res=450, width=w, height=h)
        sp = 0.1
        par(mar=c(sp, 6,3,sp))
        plt.mat = emat
        mxr = max(abs(plt.mat))
        mxr = 5
        plt.mat[plt.mat > mxr] = mxr
        plt.mat[plt.mat < -mxr] = -mxr
        image(plt.mat, col=rev(colrb), zlim = c(-mxr, mxr), axes=F)
        text(x=parpos(1, .02),
             y=seq(0,1,length.out=ncol(plt.mat)),
             labels=colnames(plt.mat), xpd=TRUE, adj=1, cex=.9)
        text(y=parpos(2, -1.012), x=seq(0,1,length.out=nrow(plt.mat)),
                 labels=rownames(plt.mat), xpd=TRUE, adj=c(.5,0))
        ind = which(plt.mat != 0,arr.ind=T)
        vals = ifelse(pmat[ind] > 2, ifelse(pmat[ind] > 3, '***', '**'),'')
        xat = seq(0,1,length.out=nrow(plt.mat))
        yat = seq(0,1,length.out=ncol(plt.mat))
        text(x=xat[ind[,1]], y=yat[ind[,2]],
             labels=vals, xpd=TRUE, cex=.8, col=ifelse(vals > 35,'grey90','black'))
        box(lwd=.5)
        dev.off()
         
        pmlist[[path]] = pmat
        emlist[[path]] = emat
    }


    scale = 0.6
    h = 0.5 + ncol(emat) * scale * .2
    w = 0.5 + nrow(emat) * scale * 1.5
    png(paste0(ctpref, 'overall_allpath_perregion_estimates.png'), units='in', res=450, width=w, height=h)
    sp = 0.1
    layout(matrix(1:4, ncol=4), widths=c(.5, 1,1,1))
    par(mar=c(sp, sp,3,sp))
    geneord = colnames(pmlist[['nft']])
    for (path in pathlist){
        plt.mat = emlist[[path]][,geneord]
        pmat = pmlist[[path]][,geneord]
        mxr = 5
        plt.mat[plt.mat > mxr] = mxr
        plt.mat[plt.mat < -mxr] = -mxr
        if (path == pathlist[[1]]){
            image(plt.mat, col='white', axes=F)
            text(x=parpos(1, -1), y=seq(0,1,length.out=ncol(plt.mat)),
                 labels=colnames(plt.mat), xpd=TRUE, adj=1, cex=.9)
        }
        image(plt.mat, col=rev(colrb), zlim = c(-mxr, mxr), axes=F)
        text(y=parpos(2, -1.002), x=seq(0,1,length.out=nrow(plt.mat)),
             labels=rownames(plt.mat), xpd=TRUE, adj=c(.5,0), font=2)
        ind = which(plt.mat != 0,arr.ind=T)
        vals = ifelse(pmat[ind] > 2, ifelse(pmat[ind] > 3, '***', '**'),'')
        xat = seq(0,1,length.out=nrow(plt.mat))
        yat = seq(0,1,length.out=ncol(plt.mat))
        text(x=xat[ind[,1]], y=yat[ind[,2]],
             labels=vals, xpd=TRUE, cex=.8, col=ifelse(vals > 35,'grey90','black'))
        mtext(toupper(path), side=3, line=1.5)
        box(lwd=.5)
    }
    dev.off()


    # Look at which dataframes are 

    # Plot out the coefficients:
    mat1 = as.matrix(do.call(rbind, cl1))
    mat2 = as.matrix(do.call(rbind, cl2))

    scale = 0.75
    h = 0.5 + ncol(mat2) * scale
    w = 0.5 + nrow(mat2) * scale
    png(paste0(ctpref, 'gene_', gene, '_lmm_coeff.png'), units='in', res=450, width=w, height=h)
    sp = 0.1
    par(mar=c(sp, 5,3,sp))
    mxr = max(abs(mat2))
    image(mat2, col=rev(colrb), zlim=c(-mxr, mxr), axes=F)
    text(x=parpos(1, .03),
         y=seq(0,1,length.out=ncol(mat2)),
         labels=sub(":", " by \n", colnames(mat2)), xpd=TRUE, adj=1, cex=.9)
    text(y=parpos(2, -1.012), x=seq(0,1,length.out=nrow(mat2)),
             labels=rownames(mat2), xpd=TRUE, adj=.5)
    ind = which(mat2 != 0,arr.ind=T)
    vals = round(mat2[ind], 2)
    xat = seq(0,1,length.out=nrow(mat2))
    yat = seq(0,1,length.out=ncol(mat2))
    text(x=xat[ind[,1]], y=yat[ind[,2]],
         labels=vals, xpd=TRUE, cex=.8)
    mtext(gene, side=3, line=1, font=2)
    box(lwd=.5)
    dev.off()


    mat2 = as.matrix(df2[,-c(1:2)])
    rownames(mat2) = paste0(df2[,1], '_', df2[,2])
    mat2 = reord(mat2)
    mat2 = t(mat2)
    # image(t(mat2), col=colrb)

    scale = 0.75
    h = 0.5 + ncol(mat2) * scale / 10
    w = 0.5 + nrow(mat2) * scale
    png(paste0(ctpref, 'gene_', gene, '_dfp_coeff.png'), units='in', res=450, width=w, height=h)
    sp = 0.1
    par(mar=c(sp, 5,3,sp))
    mxr = max(abs(mat2))
    image(mat2, col=rev(colrb), zlim=c(-mxr, mxr), axes=F)
    text(x=parpos(1, .01),
         y=seq(0,1,length.out=ncol(mat2)),
         labels=colnames(mat2), xpd=TRUE, adj=1, cex=.5)
    text(y=parpos(2, -1.01), x=seq(0,1,length.out=nrow(mat2)),
             labels=sub(":", " by \n", rownames(mat2)), xpd=TRUE, adj=.5, cex=0.7)
    dev.off()


    # ------------------------------
    # Pathology vs. gene expression:
    # ------------------------------
    print(gene)
    gene = 'QDPR'
    gene = 'HSPA1A'
    for (gene in genes){
        gene = sub("-","\\.", gene)
        subdf = df[df$symbol == gene,]
        l2 = lm(asform(c('value ~ (', paste(pathlist, collapse='+'), ') * region')), data=subdf)

        g3 = plot_model(l2, type = "pred", title=paste('Prediction of', path, 'from', gene), terms=c('region',pathlist))


        dat <- ggeffects::ggpredict(model = l2, terms = c('region',pathlist), ci.lvl = 0.95, type = 'fe')
        dat = as.data.frame(dat)
        colnames(dat)[1] = 'region' 
        colnames(dat)[6:8] = pathlist



        dat <- ggeffects::ggpredict(model = l2, terms = c(pathlist, 'region'), ci.lvl = 0.95, type = 'fe')
        dat = as.data.frame(dat)
        colnames(dat)[1] = 'nft' 
        colnames(dat)[6:8] = c(pathlist[2:3], 'region')

        ggplot(dat, aes(nft,predicted, pch=factor(plaq_n), lty=factor(plaq_n), color=factor(plaq_d))) + 
            facet_wrap(~region) + 
            # geom_point() + 
            geom_line() + 
            theme_pubr()
            # geom_segment(aes(




        ggsave(paste0(rpref, 'gene_', gene, '_multipath_coeff_int_plot.png'), dpi=450, units='in', width=7, height=8)
    }






    # Try with multiple genes:
    subdf = df[df$symbol %in% c('APOE','RASGEF1B', 'MT.ND3','DPYD','SPP1'),]
    l0 = lmer(asform(c(path, '~ region + (1|projid)')), data=subdf)
    l1 = lmer(asform(c(path, '~ value*symbol + region + (1|projid)')), data=subdf)
    l2 = lmer(asform(c(path, '~ value*symbol + region + value*symbol*region + (1|projid)')), data=subdf)
    av = anova(l0, l1, l2)

    plot_model(l2, show.values=T)
    plot_model(l2, type = "int")

    # ----------------------
    # Experiments with STAN:
    # ----------------------
    require("rstanarm", quietly = TRUE)
    library(ggplot2)
    theme_set(theme_sjplot())

    ppred = rstantools::posterior_predict(m)

    gene = 'APOE'
    subdf = df[df$symbol == gene,]
    m <- stan_glmer(asform(c(path, '~ value + region + value * region + (1|projid)')), 
                    data = subdf, chains = 1)

    # g1 = plot_model(m, title=paste('Modeling', path, 'by', gene))
    g1 = plot_model(m, 
                    bpe = "mean", bpe.style = "dot", prob.inner = .4, prob.outer = .95, 
                    bpe.color = 'black',
                    # dot.size = 3,
                    line.size = 1.5,
                    title=paste('Stan model,', path, 'by',gene))
    g2 = plot_model(m, type = "int", title=paste('Prediction of', path, 'from', gene))
    gplot = g1 + g2
    ggsave(paste0(ctpref, 'gene_', gene, '_', path, '_stan_coeff_int_plot.png'), gplot,
           dpi=450, units='in', width=7, height=4)

    # default model
    plot_model(m)
    # same model, with mean point estimate, dot-style for point estimate
    # and different inner/outer probabilities of the HDI

    slong = gather(subdf[,c('projid','region', pathlist)], pathvar, pathval, -projid, -region)
    slong = merge(slong, subdf)

        
}


