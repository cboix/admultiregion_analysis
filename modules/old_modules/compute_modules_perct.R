#!/usr/bin/R
# --------------------------
# Compute cell type modules:
# Updated: 05/25/21
# --------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)

# For nebula: 
library(nebula)
library(DESeq2)
library(RUVSeq)
library(qvalue)

# For plotting
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(ggrepel)
print(version)

celltype = 'Ast'
subtype = 'Ast'
region = 'allregions'
path = 'nft'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    print("No arguments supplied: celltype subtype region chunksize chunk ascertainment")
    # stop("No arguments supplied: celltype subtype region chunksize chunk ascertainment")
} else {        
    celltype = args[1]
    subtype = args[2]
    region = args[3]
    path = args[4]
}

# Data loader:
commandArgs = function(x){ c(celltype, subtype, region)}
source(paste0(bindir, 'multiRegion/load_difftl_data.R'))
gcout = gc()

prefstr = paste(celltype, subtype, region, path, sep="_")

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'modules_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# -----------------------------------------------
# 1. Compute unbiased modules across the full matrix
# 2. Regress out region, subtype, and other signal and recompute
# -----------------------------------------------
# outtsv = paste0(regdir, 'modules.', prefstr,'.tsv.gz')
outrda = paste0(regdir, 'modules.', prefstr,'.rda')
# nsigfile = paste0(regdir, 'modules.', prefstr,'.nsig.tsv')
if (!file.exists(outrda)){
    # --------------
    # Subset matrix:
    # --------------
    # pctcut = 0.1
    pctcut = 0.2
    keep.genes = names(pctcells)[pctcells > pctcut]
    keep.genes = keep.genes[grep("^RP[0-9]*-",keep.genes, invert=TRUE)] # Remove ribosomal genes
    keep.genes = keep.genes[grep("^RP[SL]",keep.genes, invert=TRUE)]
    mat = mat[keep.genes, pathdf$barcode]
    print(paste("[STATUS] Subsetting matrix to", 
                paste0(dim(mat), collapse = ' x '),'(g x c)'))
    gcout = gc()

    # Testing regression approach:
    g1 = 'CLU'; g2 = 'APOE'
    g1 = 'CLU'; g2 = 'HIF3A'
    g1 = 'HILPDA'; g2 = 'IRS2'
    g1 = 'HILPDA'; g2 = 'GAPDH'
    g1 = 'HILPDA'; g2 = 'SLC6A11'
    x1 = mat[g1,]
    x2 = mat[g2,]

    pfact = 10000 / pathdf$ncounts
    submat = mat[c(g1,g2),pathdf$barcode]
    submat = sweep(submat, 2, pfact, '*')

    # Un-normalized
    cor(x1, x2) 

    # Normalized (much better):
    cor(submat[1,], submat[2,])

    # Norm + log1p
    cor(log1p(submat[1,]), log1p(submat[2,]))
    cor(log1p(submat[1,]), log1p(submat[2,]), method='spearman')

    # OPT A: regress beforehand:
    subdf = pathdf[,c('region','cell_type_high_resolution')]
    subdf$x1 = log1p(submat[1,])
    subdf$x2 = log1p(submat[2,])
    f1 = lm(x1 ~ region * cell_type_high_resolution, subdf)
    f2 = lm(x2 ~ region * cell_type_high_resolution, subdf)
    subdf$p1 = predict(f1, subdf)
    subdf$p2 = predict(f2, subdf)

    cor(subdf$p1, subdf$x1)
    cor(subdf$p2, subdf$x2)

    cor(subdf$p1, subdf$p2) # Removed comp
    cor(subdf$x1 - subdf$p1, subdf$x2 - subdf$p2) # Correlation on residuals
    subdf$r1 = subdf$x1 - subdf$p1
    subdf$r2 = subdf$x2 - subdf$p2
    cor(subdf$r1,subdf$r2) # Correlation on residuals
    cor(subdf$r1,subdf$r2, method='spearman') # Spearman is less stable on residuals

    # Very inflated p-values due to covariance - could decorr. here:
    f12 = lm(x1 ~ x2 + region * cell_type_high_resolution, subdf)
    coefficients(summary(f12))

    # Testing regression approach:
    g1 = 'CLU'; g2 = 'APOE'
    g1 = 'CLU'; g2 = 'HIF3A'
    g1 = 'HILPDA'; g2 = 'IRS2'
    g1 = 'HILPDA'; g2 = 'GAPDH'
    g1 = 'HILPDA'; g2 = 'SLC6A11'
    x1 = mat[g1,]
    x2 = mat[g2,]

    # For corr
    submat = as.matrix(mat[,pathdf$barcode])
    submat = sweep(submat, 2, pathdf$ncounts / 10000, '/')

    # Regress out subtype + 
    library(Rfast)
    adjmat = 0 * submat
    adjfill = rowSums(adjmat)
    mdx = model.matrix(~ region * cell_type_high_resolution, subdf)
    for (i in 1:nrow(mat)){
        g2 = rownames(mat)[i]
        if (i %% 10 == 0){print(i)}
        f2 = lmfit(y=log1p(submat[g2,]), x=mdx)
        adjmat[g2,] = log1p(submat[g2,]) - mdx %*% f2$be
    }

    # Alternatively, compute the Mahalanobis statistic 



    g1 = 'APOE'
    subdf = pathdf[,c('region','cell_type_high_resolution')]
    subdf$x1 = log1p(submat[g1,])
    f1 = lm(x1 ~ region * cell_type_high_resolution, subdf)
    adjmat[g1,] = subdf$x1 - predict(f1, subdf)
    adjfill[g1] = 1

    for (g2 in rownames(mat)){
        corlist = cor(log1p(submat[g1,]), log1p(submat[g2,]))
        if (adjfill[g2] == 0) {
            subdf$x2 = log1p(submat[g2,])
            f2 = lm(x2 ~ region * cell_type_high_resolution, subdf)
            adjmat[g2,] = subdf$x2 - predict(f2, subdf)
            adjfill[g2] = 1
        }
    }

    # ----------------------------
    # Basic correlation + jaccard:
    # ----------------------------

    # TODO: Diff corr?
    # g500 = rownames(mat)[1:500]
    g500 = rownames(mat) # Only for cluster
    submat = as.matrix(mat[g500,pathdf$barcode])
    submat = sweep(submat, 2, pathdf$ncounts / 10000, '/')
    cr = cor(t(submat))
    # diag(cr) = 0

    # Plot this:
    w = 12
    png(paste0(imgpref, 'corheatmap_',prefstr, '_full.png'), res=450, units='in', width=w, height=w)
    draw(Heatmap(cr,
            row_names_gp=gpar(fontsize=4),
            column_names_gp=gpar(fontsize=4)))
    dev.off()

    # Jaccard - faster but not normalized:
    zmat = mat > 0
    # Intersections:
    imat = zmat %*% t(zmat)
    imat = as.matrix(imat)

    NG = nrow(zmat)
    mz = rowSums(zmat)
    mz = matrix(rep(mz, NG), nrow=NG, ncol=NG)
    jmat = imat / (mz + t(mz) - imat)

    # Plot this:
    w = 12
    png(paste0(imgpref, 'jaccheatmap_',prefstr, '_full.png'), res=450, units='in', width=w, height=w)
    draw(Heatmap(jmat,
            row_names_gp=gpar(fontsize=4),
            column_names_gp=gpar(fontsize=4)))
    dev.off()

    save(cr, jmat, file=outrda)

    library(uwot)
    udf = data.frame(umap(cr))
    udf$gene = rownames(cr)

    gp = ggplot(udf, aes(X1, X2, label=gene)) + 
        geom_point(cex=.5) + 
        geom_text_repel(size=2.5) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'umap_',prefstr,'_from_cr.png'), gp, units='in', dpi=450, width=8, height=8)


    # Plot a graph/extract clusters/modules (community detect):
    library(igraph)
    # TODO: develop more robust method for building network
    # cutoff = 0.08
    cutoff = 0.16

    dmat = cr
    rn = rownames(dmat)
    dmat = data.frame(dmat)
    colnames(dmat) = rn
    rownames(dmat) = rn
    dmat$T1 = rownames(dmat)
    ddf  = gather(dmat, T2, sim, -T1)
    ddf = ddf[ddf$sim >= cutoff,]
    ddf = ddf[ddf$T1 != ddf$T2,]

    # Remove links if diff sets:
    # ddf$samedir = ((ddf$T1 %in% upgene) & (ddf$T2 %in% upgene)) | ((ddf$T1 %in% downgene) & (ddf$T2 %in% downgene))
    # ddf = ddf[ddf$samedir,]
    nodes = sort(unique(c(ddf$T1, ddf$T2)))
    # Remove edges in opposite direction
    ddf$T1 = factor(ddf$T1, levels=nodes)
    ddf$T2 = factor(ddf$T2, levels=nodes)
    ddf = ddf[as.numeric(ddf$T1) < as.numeric(ddf$T2),]
    # Kept pct:
    dim(ddf)[1]
    dim(ddf)[1] / (dim(dmat)[1]^2)
    ddf$COLOR = 'grey25'

    # Simple network: just the links/points:
    sdf = ddf
    net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
    vcol = rep(tsp.col('indianred',.5), length(nodes))
    # vcol[nodes %in% downgene] = tsp.col('royalblue', .5)
    ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
    V(net)$size = 2
    V(net)$label = nodes
    V(net)$label.cex = .5
    V(net)$label.color = 'black'
    V(net)$color = vcol
    V(net)$frame.color <- 'black' # vcol
    V(net)$frame.color <- NA
    V(net)$pch = 19
    E(net)$color = ecol 
    elty = rep('dotted', length(sdf$sim))
    elty[sdf$sim >= .85] = 'dashed'
    elty[sdf$sim >= .95] = 'solid'
    E(net)$lty = elty
    E(net)$width = sdf$sim * 4
    E(net)$weight = sdf$sim  * 4
    set.seed(8)
    l <- layout_with_fr(net, grid='nogrid') # Usually best

    cls = cluster_louvain(net)
    memb = cls$membership
    lnet = net
    V(net)$color = snap.cols[memb + 19]
    # Clusters:
    cldf = data.frame(gene=nodes, cls=memb)
    out = sapply(unique(memb), function(x){ cat(x, "\n",cldf$gene[cldf$cls == x],"\n")})

    png(paste0(imgpref, 'crnet_basic_',prefstr, '_fromcorr.png'), res=300, units='in', width=7,height=7)
    sp = 0.25
    par(mar = rep(sp,4))
    lnet = net
    # V(lnet)$label = ''
    plot(lnet, layout=l, curved=F, label.cex=.25)
    dev.off()

    # For repel:
    source(paste0('~/data/ENCODE_DATA/bin/', 'auxiliary_function_general_repel.R'))
    V(net)$label = ''
    V(net)$size = 2.25
    l2 = l
    lrange = apply(l, 2, range)
    l2 = sweep(l2, 2, lrange[1,], '-')
    l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

    png(paste0(imgpref, 'crnet_repel_',prefstr, '_fromcorr.png'), res=300, units='in', width=7,height=7)
    sp = 0.25
    par(mar = rep(sp,4))
    plot(net, layout=l, curved=F)
    # Repel points:
    lbcex=0.45
    rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                             xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                             hjust=.5, vjust=.5, seed=1, max.iter=50000,
                             max.overlaps=25,
                             labels=nodes, cex=lbcex, pt.cex=.25)
    text(x=rdf$x, y=rdf$y, labels=rdf$lab,
         srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
    segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey50')
    # legend('topright', legend=c(paste('Up in', subtype), paste('Down in', subtype)), pch=19, col=c(tsp.col('indianred',.5),tsp.col('royalblue',.5)), bty='n', cex=1)
    # text(x=parpos(1,-.025), y=parpos(2,-.98), paste('Expression Corr. in', subtype, '(signif. genes, all regions)'), xpd=TRUE, cex=1, adj=0)
    dev.off()




    # --------------------
    # Regression approach:
    # --------------------
    g1 = 'CLU'
    g2 = 'APOE'
    x1 = mat[g1,]
    x2 = mat[g2,]

    submat = mat[c(g1,g2),pathdf$barcode]
    submat = sweep(submat, 2, pathdf$ncounts / 10000, '/')



    # Remove TH if running allregions + nft/plaq_n/plaq_d:
    if (path %in% c('nft','plaq_n','plaq_n')){
        pathdf = pathdf[pathdf$region != 'TH',]
        mat = mat[,pathdf$barcode]
    }

    # --------------------------------------------
    # Run RUV on the data at the individual level:
    # --------------------------------------------
    # Make the individual-aggregate matrix:
    NRUV = 10
    print(paste0("[STATUS] Running RUV for N=", NRUV))
    if (length(unique(pathdf$region)) > 1){
        pathdf$rp = paste0(pathdf$region, "_", pathdf$projid)
        rpids = as.character(unique(pathdf$rp))
        tform = make.tform(pathdf$rp, u=rpids)
        indvar = 'rp'
    } else {
        pids = as.character(unique(pathdf$projid))
        tform = make.tform(pathdf$projid, u=pids)
        indvar = 'projid'
    }
    data_ind = mat %*% tform 

    # Make the aggregate design matrix:
    if (indvar == 'rp'){
        uqcols = c(indvar,'region',path)
        dform = asform(c('~',path, '+ region'))
    } else {
        uqcols = c(indvar,path)
        dform = asform(c('~',path))
    }
    uqobs = unique(pathdf[,uqcols])
    rownames(uqobs) = uqobs[[indvar]]
    uqobs = uqobs[colnames(data_ind),]
    design = model.matrix(dform, data=uqobs)

    # DESeq2 object
    d_e = DGEList(data_ind, genes=rownames(data_ind))
    keep = rowSums(cpm(d_e)>1) >= 3
    d_e = d_e[keep, , keep.lib.sizes=FALSE]
    d_e = calcNormFactors(d_e, method="TMM")
    d_e = estimateGLMCommonDisp(d_e, design)
    d_e = estimateGLMTagwiseDisp(d_e, design)
    fit1 = glmFit(d_e, design)
    res1 = residuals(fit1, type="deviance")
    ruv_cov = RUVr(round(d_e$counts), 
                    as.character(rownames(d_e$counts)), 
                    k=NRUV, res1)

    # Merge the learned factors back into the data.frame:
    uqobs = cbind(uqobs, ruv_cov$W)
    pathdf = merge(pathdf, uqobs, all.x=TRUE)
    # Re-order to original
    rownames(pathdf) = pathdf$barcode
    pathdf = pathdf[colnames(mat),]

    # -----------------------------------------------
    # Run NEBULA with the RUV results and other vars:
    # -----------------------------------------------
    ruvw = paste0("W_", 1:NRUV)
    flist = c('~', path)
    flist = c(flist, ' + age_death + msex + pmi')
    flist = c(flist, ' + cpg + ngenes')
    if (length(unique(pathdf$region)) > 1){
        flist = c(flist, '+ region')
    }
    # flist = c(flist, ' + pctMT')
    if (length(unique(pathdf$cell_type_high_resolution)) > 1){
        flist = c(flist, '+ cell_type_high_resolution')
    }
    flist = c(flist, " + ", paste(ruvw, collapse=" + "))
    nb.form = asform(flist)
    print(nb.form)
    mdx = model.matrix(nb.form, data=pathdf)

    if (path %in% c('cogdxad','nrad')){
        pathstr = paste0(path, 'AD')
    } else { 
        pathstr = path 
    }
    leff = paste0('logFC_', pathstr)
    peff = paste0('p_', pathstr)

    chunksize=500
    nchunk = ceiling(nrow(mat) / chunksize)
    fulldf = c()
    offset = log10(pathdf$ncounts)
    t0 = proc.time()
    t1 = t0
    for (chunk in 1:nchunk){
        print(chunk)
        ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
        submat = mat[ind,]
        re = nebula(submat, as.character(pathdf[[indvar]]), 
                    pred=mdx, offset=offset, model='PMM') 
        rdf = re$summary
        resdf = rdf[order(rdf[[peff]]),c('gene',leff, peff)]
        names(resdf) = c('gene','logFC','p')
        fulldf = rbind(fulldf, resdf)
        t1 = est.endtime(t0, t1, chunk, nchunk)
    }
    fulldf = fulldf[order(fulldf$p),]
    fulldf$padj = p.adjust(fulldf$p, 'fdr')
    fulldf$q = qvalue(fulldf$p)$q

    pcut = 0.05
    fulldf$col = 1 * (fulldf$q < pcut) * (2 - 1 * (fulldf$logFC < 0))
    fulldf$pc = pctcells[fulldf$gene]
    labdf = fulldf[fulldf$col != 0,]

    pcols = brewer.pal(12,'Paired')
    gplot = ggplot(fulldf, aes(logFC, -log10(p), color=factor(col))) + 
        geom_point(cex=.25) + 
        geom_text_repel(data=labdf, aes(logFC, -log10(p), label=gene, color=factor(col)), size=2, max.overlaps=20) + 
        scale_color_manual(values=c('grey80',pcols[1],pcols[5])) + 
        scale_y_continuous(expand=c(0,0)) + 
        theme_pubr() + theme(legend.position='none')
    ggsave(paste0(imgpref, 'volcano_',prefstr,'_nebula_frommat.png'),gplot, units='in', dpi=450, width=6, height=6)
    ggsave(paste0(imgpref, 'volcano_',prefstr,'_nebula_frommat.pdf'),gplot, units='in', dpi=450, width=6, height=6)

    # Write out nsig: 
    ctvec = c(celltype=celltype, subtype=subtype, 
              region=region, path=path, pctcut=pctcut, pcut=pcut)
    nsig = t(c(ctvec, table(fulldf$col)))
    write.table(nsig, nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(outtsv), quote=F, row.names=F, sep="\t")
    save(resdf, fulldf, nsig, file=outrda)
} else {
    print("[STATUS] Regression output files already exist")
}

