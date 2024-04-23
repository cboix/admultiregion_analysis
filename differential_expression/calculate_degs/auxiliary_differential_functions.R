#!/usr/bin/R
# ----------------------------------------
# Auxiliary functions for DE calculations:
# Updated: 02/23/22
# ----------------------------------------
library(cbrbase)

library(tidyr)
library(Matrix)

# For nebula: 
library(nebula)
library(DESeq2)
library(RUVSeq)
library(qvalue)


# Same as plotting settings
saveGGplot = function(gp, pltprefix, w, h){
    require(ggplot2)
    ggsave(paste0(pltprefix, '.pdf'), gp, units='in', dpi=450, width=w, height=h)
    ggsave(paste0(pltprefix, '.png'), gp, units='in', dpi=450, width=w, height=h)
    print(pltprefix)
}


# Functions:
# ----------
# Subset the counts matrix:
subsetMatrixForDE = function(mat, pathdf, pctcells, pctcut=NULL, selgenes=NULL){
    if (is.null(selgenes)){
        if (!is.null(pctcut)){
            keep.genes = names(pctcells)[pctcells > pctcut]
        }
        keep.genes = keep.genes[grep("^RP[0-9]*-",keep.genes, invert=TRUE)] # Remove ribosomal genes
        keep.genes = keep.genes[grep("^RP[SL]",keep.genes, invert=TRUE)]
    } else {
        keep.genes = selgenes[selgenes %in% rownames(mat)]
    }
    mat = mat[keep.genes, pathdf$barcode]
    print(paste("[STATUS] Subsetting matrix to", 
                paste0(dim(mat), collapse = ' x '),'(g x c)'))
    return(mat)
}


# Run RUV on the sample-level pseudobulk:
runRUVpsbulk = function(pathdf, mat, path){
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

    # Return metadata with the RUV factors:
    return(pathdf)
}


# Make the regression formula and design matrix:
makeRegFormula = function(pathdf, path, nruv, int.var=NULL){
    flist = c('~', path)
    if (!is.null(int.var)){ flist = c(flist, '*', int.var) }
    flist = c(flist, ' + age_death + msex + pmi')
    flist = c(flist, ' + cpg + ngenes')
    if (length(unique(pathdf$region)) > 1){
        flist = c(flist, '+ region')
    }
    # flist = c(flist, ' + pctMT')  # Not used
    if (length(unique(pathdf$cell_type_high_resolution)) > 1){
        flist = c(flist, '+ cell_type_high_resolution')
    }
    if (nruv > 0){
        ruvw = paste0("W_", 1:nruv)
        flist = c(flist, " + ", paste(ruvw, collapse=" + "))
    }

    nb.form = asform(flist)
    print(nb.form)
    mdx = model.matrix(nb.form, data=pathdf)

    # Set the ad variable string:
    if (path %in% c('cogdxad','nrad', 'braaksc.early', 'braaksc.ad')){
        pathstr = paste0(path, 'AD')
    } else { 
        pathstr = path 
    }
    if (!is.null(int.var)){
        if (int.var %in% c('cogdxad','nrad', 'braaksc.early', 'braaksc.ad')){
            intstr = paste0(int.var, 'AD')
        } else { 
            intstr = int.var
        }
        pathstr = c(intstr, paste0(pathstr,':',intstr))
    }
    leff = paste0('logFC_', pathstr)
    peff = paste0('p_', pathstr)
    return(list(nb.form=nb.form, mdx=mdx,
            leff=leff, peff=peff, pathstr=pathstr))
}


# Functions for MAST, specifically:
# ---------------------------------
# For normalizing to log(norm + 1)
makeNormMat = function(mat, csm){
    norm = as.matrix(mat)
    norm = sweep(norm, 2, csm / 10000,'/')
    norm = log(norm + 1)
    return(norm)
}


# Function for MAST:
runSingleChunkMAST = function(mat, covmat, fdata, 
                              ind, csm, nb.form, pathstr){
    # Make matrix and chunked sca object:
    norm = makeNormMat(mat[ind,], csm)
    sca = FromMatrix(exprsArray=norm, cData=covmat, 
                     fData=fdata[ind,,drop=F], check_sanity=TRUE)
    # Calculate regression, run LRT, refit:
    zlm.obj = zlm(nb.form, sca)
    # NOTE: if multiple pathstr, this will be multiple x slower
    summaryCond <- summary(zlm.obj, doLRT=pathstr) 
    # MAST results:
    summaryDt <- summaryCond$datatable
    resdf = NULL
    for (pstr in pathstr){
        fcHurdle <- merge(
            # Hurdle P values
            summaryDt[contrast==pstr & component=='H',
                .(primerid, `Pr(>Chisq)`)],
            # logFC coefficients
            summaryDt[contrast==pstr & component=='logFC',
                .(primerid, coef, ci.hi, ci.lo)], by='primerid')
        chunkdf = data.frame(fcHurdle)
        if (length(pathstr) == 1){
            colnames(chunkdf) = c('gene','p','coef','ci.hi','ci.lo')
        } else {
            colnames(chunkdf) = c('gene',
                paste0(c('p','coef','ci.hi','ci.lo'), '_', pstr))
        }
        if (is.null(resdf)){
            resdf = chunkdf
        } else {
            resdf = merge(resdf, chunkdf, all=TRUE)
        }
    }
    return(resdf)
}


