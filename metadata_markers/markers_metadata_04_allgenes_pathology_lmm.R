#!/usr/bin/R
# -----------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Last updated mid-2020 (exploratory analysis)
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(rhdf5)
library(lme4)

celltype = 'Microglia'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need celltype")
} else {        
    celltype = args[1]
}

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/markers/')
imgpref = paste0(plotdir, 'pred_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
pathlist = c('nft','plaq_d','plaq_n')

# ------------------------
# Load pathology measures:
# ------------------------
all.path.rda = paste0(datadir, prefix, '.allpath.Rda')
if (!file.exists(all.path.rda)){
    pqdf = NULL
    for (i in 1:length(pathlist)){
        path = pathlist[i]
        pathfile = paste0(datadir, prefix, '.', path, '.tsv.gz')
        pdf = read.delim(pathfile, header=F)
        names(pdf) = c('barcode',path)
        if (is.null(pqdf)){ pqdf = pdf } else { pqdf = merge(pqdf, pdf) }
    }
    print(dim(pqdf))
    save(pqdf, file=all.path.rda)
} else {
    load(all.path.rda) 
}


# -------------------------------
# Load in and plot for each cell:
# -------------------------------
path = 'nft'
ctpref = paste0(imgpref, celltype, '_', path, '_')
lblset = 'leiden_r5_n50'
h5file = paste0(datadir, prefix, '.', lblset, '.full_plaq_n_', celltype, '.hdf5')
regfile = paste0(datadir, prefix, '.', lblset, '.full_pathregression.', celltype, '.Rda')
regtsv = paste0(datadir, prefix, '.', lblset, '.full_pathregression.', celltype, '.tsv.gz')
# Read in pathology and coefficients:
if (path != 'nft'){ celltype = paste0(path, '_', celltype) }
pathdf = read.delim(paste0(datadir, prefix, '.', lblset, '.pathpred_', celltype, '.tsv.gz'), header=T, stringsAsFactors=F)
pathdf$region = sub('_.*','',pathdf$barcode)
pathdf = merge(pathdf, celldf, all.x=TRUE)
pathdf = merge(pathdf, pqdf)

# Get the expression matrix for top genes (testing lmm models):
# Slice the matrix for these genes:
h5attr = h5ls(h5file)
ngenes = as.numeric(h5attr[5][h5attr[2] == 'genes'])

if (!file.exists(regfile)){
    chunksize = 1000
    nchunk = floor(ngenes / chunksize) + 1
    regpdf = c()
    for (i in 1:nchunk){
        print(i)
        ind = (1 + (i-1) * chunksize):min(c(i * chunksize, ngenes))
        h5f = H5Fopen(h5file)
        genes = h5f$genes
        bcs = h5f$barcodes
        # Open handle, extract genes we care about and close:
        h5d = h5f&"matrix"
        mat = t(h5d[ind,])
        H5Dclose(h5d)
        H5Fclose(h5f)
        colnames(mat) = genes[ind]
        rownames(mat) = bcs

        # Order as pathdf, merge with pathology, metadata:
        print("[STATUS] Merging data:")
        mat = mat[pathdf$barcode,]
        df = data.frame(mat)
        df$barcode = rownames(df)
        df = gather(df, symbol, value, - barcode)
        rm(mat)

        print("[STATUS] Data merged with expression values")
        t0 = proc.time()
        for (j in 1:length(ind)){
            t1 = proc.time()
            gene = genes[ind[j]]
            cat(ind[j], gene, '\t')
            gene = sub("-","\\.", gene)
            # subdf = df[df$symbol == gene,]
            baseind = 1:nrow(pathdf)
            geneind = nrow(pathdf) * (j-1) + baseind
            subdf = df[geneind,]
            subdf = cbind(subdf, pathdf[,c('barcode','region', pathlist, 'projid')])
            if (var(subdf$value) > 0){
                cl1 = list()
                cl2 = list()
                # pvals = matrix(0, nrow=2, ncol=3, dimnames=list(NULL,pathlist))
                for (path in pathlist){
                    cat(path,'\t')
                    # l2 = lm(asform(c(path, '~ region + value:region')), data=subdf)
                    l2 = lmer(asform(c('value ~ region + ', path, ':region + (1|projid)')), data=subdf)
                    cvals = summary(l2)$coefficients
                    cind = grep(path, rownames(cvals))
                    if (length(cind) > 0){
                        cvals = as.data.frame(cvals[cind,,drop=F])
                        cvals$gene = gene
                        cvals$path = path
                        cvals$region = sub(paste0(':', path), '', sub('region', '', rownames(cvals)))
                        regpdf = rbind(regpdf, cvals)
                    }
                }
            }
            t2 = (proc.time() - t1)[3]
            ttot = (proc.time() - t0)[3]
            names(t2) = NULL
            names(ttot) = NULL
            cat(paste0(round(t2,1),'s\t'))
            cat(paste0(round(ttot,1),'s\n'))
        }
        # rdf = regpdf[order(regpdf$t, decreasing=T),]
    }

    # TODO: Get pval from the t-stats
    regpdf$log10p = -log10(regpdf[,'Pr(>|t|)'])

    save(regpdf, file=regfile)
    write.table(regpdf, file=gzfile(regtsv), quote=F, row.names=F, col.names=T, sep="\t")
} else {
    load(regfile)
}




# ---------------------------------------------------------------------
# Plot the region-specific coefficients and significance for each gene:
# Plot as a heatmap, cluster by which genes show different effects by region.
# ---------------------------------------------------------------------
# Only genes with at least one significant hit: 
regpdf$padj = p.adjust(regpdf[['Pr(>|t|)']])
regpdf$log10adj = -log10(regpdf$padj)
gdf = aggregate(log10adj ~ gene, regpdf, max)
kgene = gdf$gene[gdf$log10adj > 5]
sdf = regpdf[regpdf$gene %in% kgene,]

path = 'nft'
pmlist = list()
emlist = list()
for (path in pathlist){
    srdf = sdf[sdf$path == path, c('gene','log10adj', 'Estimate', 'region')]
    ewide = spread(srdf[,c('gene','Estimate','region')], region, Estimate)
    pwide = spread(srdf[,c('gene','log10adj','region')], region, log10adj)
    emat = as.matrix(ewide[,-1])
    pmat = as.matrix(pwide[,-1])
    emat[is.na(emat)] = 0
    pmat[is.na(pmat)] = 0
    rownames(emat) = ewide$gene
    rownames(pmat) = pwide$gene
    regst = c('EC','HC','MT','AG','PFC')
    emat = emat[,regst]
    pmat = pmat[,regst]
    colnames(pmat) = paste0(path, '_', colnames(pmat))
    colnames(emat) = paste0(path, '_', colnames(emat))
    emat = t(emat)
    pmat = t(pmat)
    emlist[[path]] = emat
    pmlist[[path]] = pmat
}


amat = do.call(rbind, emlist)
pmat = do.call(rbind, pmlist)
# rownames(amat) = paste0(sapply(pathlist, times=5, rep), "_", regst)
# rownames(pmat) = paste0(sapply(pathlist, times=5, rep), "_", regst)
plim = 5
amat[pmat < plim] = 0
amat = reord(t(amat))
pmat = t(pmat[colnames(amat),])

mxr = 5
amat[amat < -mxr] = -mxr
amat[amat > mxr] = mxr
mmarg = apply(amat, 1, max)
ind = mmarg >= 4
# Present across:
nhit = apply(abs(amat) > 0, 1, sum)
NTOP = 50
topgenes = names(head(sort(nhit, decreasing=T), NTOP))
ind = rownames(amat)[rownames(amat) %in% topgenes]
# ind = 1:nrow(amat) 
plt.mat = t(amat[ind,])
sub.pmat = t(pmat[ind,])


w = (10 + nrow(plt.mat)) / 8
h = (10 +  ncol(plt.mat)) / 8
png(paste0(ctpref, 'top',NTOP,'_consistent_genes.png'), units='in', res=450, width=w, height=h)
sp = 0.1
par(mar=c(5, 5,sp,sp))
image(plt.mat, col=rev(colrb), zlim=c(-mxr, mxr), axes=F, useRaster=T)
text(x=parpos(1, 0.01), y=seq(0,1,length.out=ncol(plt.mat)),
     labels=colnames(plt.mat), xpd=TRUE, adj=1, cex=.9)
text(y=parpos(2, 0.01), x=seq(0,1,length.out=nrow(plt.mat)),
     labels=rownames(plt.mat), xpd=TRUE, srt=90,adj=1, font=1)
abline(v=seq(par()$usr[1], par()$usr[2], length.out=4))
box()
dev.off()


# Per path hit:
nhitdf = data.frame(nft=apply(abs(amat[,1:5]) > 0, 1, sum),
                    plaq_d=apply(abs(amat[,6:10]) > 0, 1, sum),
                    plaq_n=apply(abs(amat[,11:15]) > 0, 1, sum),
                    total=nhit,
                    gene=rownames(amat))
nhitdf = nhitdf[order(nhitdf$total, decreasing=T),]

# Plot mainly PLAQ, not NFT:
plaqdf = nhitdf[nhitdf$nft == 0,]
topgenes = head(plaqdf$gene,NTOP)
ind = rownames(amat)[rownames(amat) %in% topgenes]
plt.mat = t(amat[ind,])
sub.pmat = t(pmat[ind,])

w = (10 + nrow(plt.mat)) / 8
h = (10 +  ncol(plt.mat)) / 8
png(paste0(ctpref, 'top',NTOP,'_consistent_genes_plaq.png'), units='in', res=450, width=w, height=h)
sp = 0.1
par(mar=c(5, 5,sp,sp))
image(plt.mat, col=rev(colrb), zlim=c(-mxr, mxr), axes=F, useRaster=T)
text(x=parpos(1, 0.01), y=seq(0,1,length.out=ncol(plt.mat)),
     labels=colnames(plt.mat), xpd=TRUE, adj=1, cex=.9)
text(y=parpos(2, 0.01), x=seq(0,1,length.out=nrow(plt.mat)),
     labels=rownames(plt.mat), xpd=TRUE, srt=90,adj=1, font=1)
abline(v=seq(par()$usr[1], par()$usr[2], length.out=4))
box()
dev.off()


# Plot mainly PLAQ, not NFT:
nftdf = nhitdf[nhitdf$nft >= nhitdf$total - 2,]
topgenes = head(nftdf$gene,NTOP)
ind = rownames(amat)[rownames(amat) %in% topgenes]
plt.mat = t(amat[ind,])
sub.pmat = t(pmat[ind,])

w = (10 + nrow(plt.mat)) / 8
h = (10 +  ncol(plt.mat)) / 8
png(paste0(ctpref, 'top',NTOP,'_consistent_genes_nft.png'), units='in', res=450, width=w, height=h)
sp = 0.1
par(mar=c(5, 5,sp,sp))
image(plt.mat, col=rev(colrb), zlim=c(-mxr, mxr), axes=F, useRaster=T)
text(x=parpos(1, 0.01), y=seq(0,1,length.out=ncol(plt.mat)),
     labels=colnames(plt.mat), xpd=TRUE, adj=1, cex=.9)
text(y=parpos(2, 0.01), x=seq(0,1,length.out=nrow(plt.mat)),
     labels=rownames(plt.mat), xpd=TRUE, srt=90,adj=1, font=1)
abline(v=seq(par()$usr[1], par()$usr[2], length.out=4))
box()
dev.off()



# -------------------------------
# Go back to plot specific genes:
# -------------------------------
# For model plots:
library(sjPlot)
library(sjmisc)
library(patchwork)
library(ggplot2)
library(ggpubr)

# keep.genes = c('CDKN2B','ATP8B4','FTL')
geneset = 'SLC'
geneset = 'TNC'
geneset = 'RYR2'
if (geneset == 'APOE'){
    keep.genes = c('APOE','MT-ND3','HLA-DRB1')
} else if (geneset == 'HSP'){ 
    keep.genes = c('HSP90AA1','HSPA1A','HSPA1B')
} else if (geneset == 'SLC'){ 
    keep.genes = c('SLC39A11', 'SLC26A3','SLC27A6')
} else if (geneset == 'TNC'){ 
    keep.genes = c('TNC', 'DSCAM','CEBPD')
} else if (geneset == 'RYR2'){ 
    keep.genes = c('RYR2', 'LINGO2','ANLN')
}
h5f = H5Fopen(h5file)
genes = h5f$genes
bcs = h5f$barcodes
ind = which(genes %in% keep.genes)
# Open handle, extract genes we care about and close:
h5d = h5f&"matrix"
kmat = t(h5d[ind,])
H5Dclose(h5d)
H5Fclose(h5f)
colnames(kmat) = genes[ind]
rownames(kmat) = bcs

# Order as pathdf, merge with pathology, metadata:
kmat = kmat[pathdf$barcode,]
kdf = data.frame(kmat)
kdf$barcode = rownames(kdf)
kdf = gather(kdf, symbol, value, - barcode)
rm(kmat)

print("[STATUS] Data merged with expression values")
kpdf = c()
t0 = proc.time()
glist = list()
gplot = NULL
for (gene in keep.genes){
    t1 = proc.time()
    cat(gene, '\t')
    gene = sub("-","\\.", gene)
    subdf = kdf[kdf$symbol == gene,]
    subdf = cbind(subdf, pathdf[,c('barcode','region', pathlist)])
    for (path in pathlist){
        cat(path,'\t')
        l2 = lm(asform(c(path, '~ region + value:region')), data=subdf)
        cvals = summary(l2)$coefficients
        cind = grep('value', rownames(cvals))
        if (length(cind) > 0){
            cvals = as.data.frame(cvals[cind,,drop=F])
            cvals$gene = gene
            cvals$path = path
            cvals$region = sub(':value', '', sub('region', '', rownames(cvals)))
            kpdf = rbind(kpdf, cvals)
        }
        glist[[paste0(gene, '_', path)]] = plot_model(l2, type = "pred", title=paste('Prediction of', path, 'from', gene), terms=c('value','region')) + theme_pubr()
    }
}


garr = ggarrange(plotlist=glist, nrow=length(keep.genes), ncol=3, common.legend=TRUE)
scale = 3
ggsave(paste0(ctpref, 'gene_examples_', geneset, '_intpred_plot.png'), garr,
       dpi=450, units='in', width=4 * scale, height=3 * scale)



# -----------------------------------------------------------
# Plot the actual values, against a different scoring metric:
# -----------------------------------------------------------
keep.genes = c('HSP90AA1','HSPA1A','HSPA1B', 'SLC39A11', 'SLC26A3','SLC27A6')
h5f = H5Fopen(h5file)
genes = h5f$genes
bcs = h5f$barcodes
ind = which(genes %in% keep.genes)
# Open handle, extract genes we care about and close:
h5d = h5f&"matrix"
kmat = t(h5d[ind,])
H5Dclose(h5d)
H5Fclose(h5f)
colnames(kmat) = genes[ind]
rownames(kmat) = bcs

# Order as pathdf, merge with pathology, metadata:
kmat = kmat[pathdf$barcode,]
kdf = data.frame(kmat)
kdf$barcode = rownames(kdf)
kdf = gather(kdf, symbol, value, - barcode)
kdf = merge(kdf, pathdf[,c('barcode','region', pathlist)])
kdf$on.val = kdf$value > 1

# Genes are not correlated, except HSP genes:
crmat = cor(kmat)
zmat = 1 * (kmat > 0)
aind = grep("^AG", rownames(zmat))

g1 = 'HSPA1A'
g2 = 'SLC27A6'
tmat = table(zmat[aind,g1], zmat[aind,g2])
# Expected:
tot = sum(tmat)
e = (sum(tmat[,2]) / tot) * (sum(tmat[2,]) / tot) * tot
print(tmat[2,2] - e)
print((tmat[2,2] - e)^2 / e)

fisher.test(tmat)



reg.cols2 = col.paired[seq(2,10,2)]
names(reg.cols2) = c('EC','HC','AG','PFC','MT')

ggplot(kdf, aes(region, plaq_n, fill=region, alpha=factor(on.val))) + 
    facet_wrap(~symbol) + 
    geom_boxplot() + 
    scale_fill_manual(values=reg.cols2) + 
    theme_pubr()

# NOTE: Need to check if results are robust to removing individuals:
kdf2 = kdf[kdf$plaq_d < 20,]

ggplot(kdf2[kdf2$symbol == 'HSPA1A',], aes(value, plaq_d, color=region)) + 
    facet_wrap(~region) + 
    geom_point(alpha=.25) + 
    geom_smooth(method = 'lm') + 
    scale_color_manual(values=reg.cols2) + 
    theme_pubr()

# -------------------------------------------------------------------
# Jack-knife validation of the directions (leave out 1/8 of samples):
# -------------------------------------------------------------------
# Random dropout of the projids - will help validate significance of the coefficients/directionality:
kpdf = c()
for (gene in keep.genes){
    t1 = proc.time()
    cat(gene, '\t')
    gene = sub("-","\\.", gene)
    subdf = kdf[kdf$symbol == gene,]
    subdf = cbind(subdf, pathdf[,c('barcode','region', 'projid', pathlist)])
    projids = sort(unique(subdf$projid))
    nproj = length(projids)
    for (path in pathlist){
        cat(path,'\t')
        NPERM = 25
        for (k in 1:NPERM){
            cat(k)
            set.seed(k)
            # Jack-knife sampling:
            rand.pids = sample(projids, nproj - 6, replace=FALSE)
            pdf = subdf[subdf$projid %in% rand.pids,]
            l2 = lm(asform(c(path, '~ region + value:region')), data=pdf)
            cvals = summary(l2)$coefficients
            cind = grep('value', rownames(cvals))
            if (length(cind) > 0){
                cvals = as.data.frame(cvals[cind,,drop=F])
                cvals$gene = gene
                cvals$path = path
                cvals$perm = k
                cvals$region = sub(':value', '', sub('region', '', rownames(cvals)))
                kpdf = rbind(kpdf, cvals)
            }
        }
    }
}

# aggregate(Estimate ~ gene + path + region, kpdf, range)
rndf = aggregate(Estimate ~ gene + path + region, kpdf, function(x){quantile(x, c(0.005, .075, .5, .975, 0.995))})
# Keep only those with consistent stats:
rind = which((rndf$Estimate[,1] * rndf$Estimate[,5]) > 0)
rndf$est = 0
rndf$est[rind] = rndf$Estimate[rind,3]

# ----------------------------------------------------------------------------
# Plot specific genes against each other in each region, colored by pathology:
# ----------------------------------------------------------------------------
mdf = data.frame(kmat)
mdf$barcode = rownames(mdf)
mdf = merge(mdf, pathdf[,c('barcode','region', 'U1','U2',pathlist)])

library(viridis)

# ggplot(mdf, aes(SLC26A3, HSPA1A, color=plaq_d)) + 
ggplot(mdf, aes(SLC26A3, SLC27A6, alpha=plaq_n)) + 
    facet_wrap(~region) + 
    geom_jitter() + 
    geom_smooth(method='lm') + 
    scale_color_viridis() + 
    # scale_fill_manual(values=reg.cols2) + 
    theme_pubr()

# mdf[mdf$

# -------------------------------------------------------------------------
# Plot on the umap, demonstrate that it is not a cell-type specific effect:
# -------------------------------------------------------------------------
mdf$col = tsp.col('grey90')
mdf$col[mdf$SLC26A3 > 0] = tsp.col(col.paired[2])
mdf$col[mdf$SLC27A6 > 0] = tsp.col(col.paired[4])
mdf$col[mdf$SLC27A6 > 0 & mdf$SLC26A3 > 0] = tsp.col(col.paired[8])

xlim = range(mdf$U1)
ylim = range(mdf$U2)
png('~/test.png', units='in', res=450, width=4 * 5, height=4)
par(xaxs='i')
par(yaxs='i')
layout(matrix(1:5, ncol=5))
sp = 0.1
bsp = 1.5
cex = 0.025
for (region in names(reg.cols2)){
    ind = which(mdf$region == region)
    ind = sample(ind, length(ind), replace=FALSE)
    par(mar=c(bsp,bsp,2,sp))
    plot(mdf$U1[ind], mdf$U2[ind], col=mdf$col, 
         xlim = xlim, ylim = ylim,
         pch=19, cex=cex, axes=F)
    rect(xleft=par()$usr[1], xright=par()$usr[2],
         ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
         ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
         col='grey85', border=NA, lwd=.5, xpd=TRUE)
    # with(celltype.loc, text(U1, U2, full.exttype, cex=.5, xpd=TRUE))
    mtext(region, side=3, cex=1.5, col='grey25', font=2, line=0.25)
    mtext('UMAP 1', side=1, line=0.25, cex=1.25)
    mtext('UMAP 2', side=2, line=0, cex=1.25)
}
dev.off()



# ---------------------------------------
# Regress, controlling by each pathology:
# ---------------------------------------
kpdf = c()
for (gene in keep.genes){
    t1 = proc.time()
    cat(gene, '\t')
    gene = sub("-","\\.", gene)
    subdf = kdf[kdf$symbol == gene,]
    subdf = cbind(subdf, pathdf[,c('barcode','region', 'projid', pathlist)])
    projids = sort(unique(subdf$projid))
    nproj = length(projids)
    for (path in pathlist){
        cat(path,'\t')
        NPERM = 25
        npath = pathlist[pathlist != path]
        l2 = lm(asform(c(path, '~', npath[1], '+', npath[2], '+ region + value:region')), data=subdf)
        l2 = lm(asform(c(path, '~', npath[1], ':region +', npath[2], ':region + region + value:region')), data=subdf)
        cvals = summary(l2)$coefficients
        cind = grep('value', rownames(cvals))
        if (length(cind) > 0){
            cvals = as.data.frame(cvals[cind,,drop=F])
            cvals$gene = gene
            cvals$path = path
            cvals$region = sub(':value', '', sub('region', '', rownames(cvals)))
            kpdf = rbind(kpdf, cvals)
        }
        # plot_model(l2, type = "pred", title=paste('Prediction of', path, 'from', gene), terms=c('value','region')) + theme_pubr()
    }
}
kpdf$log10p = -log10(kpdf[,'Pr(>|t|)'])
kpdf[kpdf$log10p > 4,c('gene','path','region','Estimate', 'log10p')]
kpdf[kpdf$log10p > 4,]














