#!/usr/bin/R
# -------------------------------------------------------------------------
# Run MAST + wilcoxon (Nebula in other script) using the precomp. matrices:
# Updated: 10/31/23 to add interaction
# -------------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(Matrix)

# For nebula: 
library(MAST)
library(nebula)
library(DESeq2)
library(RUVSeq)
library(qvalue)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
print(version)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    print("No arguments supplied: runfile task")
} else {        
    runfile = args[1]
    task = as.numeric(args[2])
}


# Functions for DE:
cbindir = paste0(sbindir, 'differential_expression/calculate_degs/')
source(paste0(cbindir, 'auxiliary_differential_functions.R'))


# Read run table and set run options:
# -----------------------------------
rundf = read.delim(runfile, header=T, stringsAsFactors=FALSE)
celltype = as.character(rundf[task,1])
subtype  = as.character(rundf[task,2])
region   = as.character(rundf[task,3])
path     = as.character(rundf[task,4])
int = (path == 'msex') # Only interaction model for now


# Set the directories:
# --------------------
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)


# Set the run outputs: 
# ---------------------
prefstr = gsub("[()/]","_", paste(celltype, subtype, region, path, sep="_"))
masttsv = paste0(regdir, 'mast.', prefstr,'.tsv.gz')
mastrda = paste0(regdir, 'mast.', prefstr,'.rda')
mast.nsigfile = paste0(regdir, 'mast.', prefstr,'.nsig.tsv')
wxtsv = paste0(regdir, 'wilcoxon.', prefstr,'.tsv.gz')
wxrda = paste0(regdir, 'wilcoxon.', prefstr,'.rda')
wx.nsigfile = paste0(regdir, 'wilcoxon.', prefstr,'.nsig.tsv')
pcut = 0.05
print(prefstr)


# Check if already computed. If not, load data:
if (!file.exists(mastrda) || !file.exists(wxrda)){
    # Load dataset:
    commandArgs = function(x){ c(celltype, subtype, region)}
    source(paste0(bindir, 'multiRegion/load_difftl_data.R'))

    # Subset matrix to genes in 20%+ of cells:
    pctcut = 0.2
    mat = subsetMatrixForDE(mat, pathdf, pctcells=pctcells, pctcut=pctcut)

    # Remove TH if running allregions + nft/plaq_n/plaq_d:
    if (path %in% c('nft','plaq_n','plaq_d')){
        pathdf = pathdf[pathdf$region != 'TH',]
        mat = mat[,pathdf$barcode]
    }

    # For normalize to log(TPM + 1):
    csm = colSums(mat)
}


# Compute the MAST results:
if (!file.exists(mastrda)){
    # Run RUV on the data at the individual level:
    # --------------------------------------------
    pathdf = runRUVpsbulk(pathdf, mat, path)

    # Put together the regression formula:
    # ------------------------------------
    if (int){ int.var = 'nrad' } else { int.var = NULL }
    ll = makeRegFormula(pathdf, path=path, nruv=10, int.var=int.var)

    print("[STATUS] Making MAST covariate matrix")
    fdata = data.frame(primerid=rownames(mat))
    covmat = data.frame(ll$mdx)
    covmat$wellKey = pathdf$barcode
    covmat$indiv = pathdf$projid
    covmat$ncounts = pathdf$ncounts
    # Due to way in which MAST takes arguments:
    covmat[path] = pathdf[,path]
    if (!is.null(int.var)){ covmat[int.var] = pathdf[,int.var] }
    if (length(unique(pathdf$region)) > 1){
        covmat$region = pathdf$region
    }
    if (length(unique(pathdf$cell_type_high_resolution)) > 1){
        covmat$cell_type_high_resolution = pathdf$cell_type_high_resolution
    }


    # Run MAST on normalized, chunked data:
    # -------------------------------------
    print("[STATUS] Running MAST")
    chunksize=500
    nchunk = ceiling(nrow(mat) / chunksize)
    fulldf = c()
    t0 = proc.time()
    t1 = t0
    for (chunk in 1:nchunk){
        ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
        print(paste(chunk, length(ind)))
        chunkrda = paste0(regdir, 'mast.', prefstr,'.', 
                          chunksize, '_', chunk, '.rda')
        if (!file.exists(chunkrda)){
            chunkdf = runSingleChunkMAST(mat=mat, covmat=covmat, fdata=fdata, 
                                         ind=ind, csm=csm, nb.form=ll$nb.form, pathstr=ll$pathstr)
            save(chunkdf, file=chunkrda)
        } else {
            load(chunkrda)
        }
        fulldf = rbind(fulldf, chunkdf)
        t1 = est.endtime(t0, t1, chunk, nchunk)
    }


    # Process the full results:
    # -------------------------
    fulldf$val = pctcells[fulldf$gene]
    if (length(ll$pathstr) > 1){
        for (pstr in ll$pathstr){
            lvar = paste0('coef_', pstr)
            pvar = paste0('p_', pstr)
            qvar = paste0('fdr_', pstr)
            cvar = paste0('col_', pstr)
            fulldf[[qvar]] = p.adjust(fulldf[[pvar]], 'fdr')
            fulldf = fulldf[order(fulldf[[pvar]]),]
            fulldf[[cvar]] = ((fulldf[[qvar]] < pcut) *
                (2 - 1 * (fulldf[[lvar]] < 0)))
        }
        fulldf = fulldf[fulldf$val > 0.2,]
        fulldf = fulldf[fulldf$gene != 'MALAT1',]
    } else {
        fulldf$fdr = p.adjust(fulldf$p, 'fdr')
        fulldf = fulldf[order(fulldf$p),]
        fulldf = fulldf[!is.na(fulldf$coef),]
        fulldf$col = (fulldf$fdr < pcut) * (2 - 1 * (fulldf$coef < 0))
        fulldf = fulldf[fulldf$val > 0.2,]
        fulldf = fulldf[fulldf$gene != 'MALAT1',]

        # Plot as a volcano plot:
        # -----------------------
        ntop = 50
        labdf = rbind(head(fulldf[fulldf$col == 1,], ntop),
            head(fulldf[fulldf$col == 2,], ntop))
        pcols = brewer.pal(12, 'Paired')
        gplot = ggplot(fulldf, aes(coef, -log10(fdr), color=factor(col))) + 
            geom_point(cex=.1) + theme_pubr() + 
            scale_color_manual(values=c('grey80',pcols[1], pcols[5])) + 
            geom_text_repel(data=labdf, aes(coef, -log10(fdr), color=factor(col), label=gene), size=2, max.overlaps=20) + 
            labs(x=paste('logFC ', ll$pathstr), y='log10 adj. p-value') + 
            geom_vline(xintercept=0, lty='dashed') + 
            geom_hline(yintercept=2, lty='dotted') +
            scale_y_continuous(expand=c(0,0)) + 
            theme(legend.position='none')
        pltprefix = paste0(imgpref, 'volcano_', prefstr, '_mast_frommat')
        saveGGplot(gplot, pltprefix, w=6, h=6)
    }

    # Write out the results:
    # ----------------------
    # Write out nsig: 
    if (length(ll$pathstr) > 1){
        nsig = c()
        for (pstr in ll$pathstr){
            cvar = paste0("col_", pstr)
            ndf = table(fulldf[,cvar])
            ndf = as.data.frame(ndf)
            ndf = spread(ndf, Var1, Freq)
            ndf$eff = pstr
            if (!('1' %in% colnames(ndf))){ ndf[['1']] = 0 }
            if (!('2' %in% colnames(ndf))){ ndf[['2']] = 0 }
            cn = c('0','1','2','eff')
            nsig = rbind(nsig[,cn, drop=F], ndf[,cn, drop=F])
        }
    } else {
        nsig = table(fulldf[,'col'])
        nsig = as.data.frame(nsig)
        nsig = spread(nsig, Var1, Freq)
    }
    ctvec = data.frame(celltype=celltype, subtype=subtype, 
        region=region, path=path, 
        pctcut=pctcut, pcut=pcut)
    nsig = merge(ctvec, nsig)
    write.table(nsig, mast.nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)


    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(masttsv), quote=F, row.names=F, sep="\t")
    save(fulldf, nsig, file=mastrda)
} else {
    print("[STATUS] MAST regression output files already exist")
}



# Compute differential genes by wilcoxon test as well:
# NOTE: Do not run if interaction
if ((!file.exists(wxrda)) && (!int)){
    print("[STATUS] Running Wilcoxon test")
    norm = as.matrix(mat)
    # NOTE: May not be able to fit into memory as matrix if huge.
    norm = sweep(norm, 2, colSums(norm) / 10000,'/')
    norm = log(norm + 1)

    # Running both at individual level and cell level, aggregate:
    if (region == 'allregions'){
        pathdf$pr = paste0(pathdf$projid, "_", pathdf$region)
        udf = unique(pathdf[,c('pr', path)])
        rownames(udf) = udf$pr
        pids = as.character(unique(udf$pr))
        tform = make.tform(pathdf$pr, u=pids, norm=TRUE)
    } else {
        udf = unique(pathdf[,c('projid',path)])
        rownames(udf) = udf$projid
        pids = as.character(unique(pathdf$projid))
        tform = make.tform(pathdf$projid, u=pids, norm=TRUE)
    }
    avg.nmat = norm[,pathdf$barcode] %*% tform

    # adid = which(pathdf$Apoe_e4 == 'yes')
    # ctid = which(pathdf$Apoe_e4 == 'no')
    # adid.ind = which(udf$Apoe_e4 == 'yes')
    # ctid.ind = which(udf$Apoe_e4 == 'no')

    # Splits for different variables:
    if (path %in% c('cogdxad','nrad', 'braaksc.early', 'braaksc.ad')){
        adid = which(pathdf[path] == 'AD')
        ctid = which(pathdf[path] == 'CONTROL')
        adid.ind = which(udf[path] == 'AD')
        ctid.ind = which(udf[path] == 'CONTROL')
    } else {
        adid = which(pathdf[path] > 0)
        ctid = which(pathdf[path] == 0)
        adid.ind = which(udf[path] > 0)
        ctid.ind = which(udf[path] == 0)
    }

    # Wilcoxon at cell level:
    fulldf = c()
    for (gene in keep.genes){
        # Cell level:
        xi = norm[gene,adid]
        xo = norm[gene,ctid]
        cp = wilcox.test(xi, xo)$p.value
        # Individual level:
        xi.ind = avg.nmat[gene,adid.ind]
        xo.ind = avg.nmat[gene,ctid.ind]
        ip = wilcox.test(xi.ind, xo.ind)$p.value
        # Collate:
        fulldf = rbind(fulldf, data.frame(gene=gene,
                                          AD = mean(xi), CTRL=mean(xo), 
                                          # e4=mean(xi), e3=mean(xo), 
                                          p_cell=cp, p_ind=ip))
    }

    fulldf$path = path
    fulldf$logFC = log2(fulldf$AD / fulldf$CTRL)
    # fulldf$logFC = log2(fulldf$e4 / fulldf$e3)
    fulldf$fdr = p.adjust(fulldf$p_cell, 'fdr')
    fulldf$fdr_ind = p.adjust(fulldf$p_ind, 'fdr')
    fulldf = fulldf[order(abs(fulldf$logFC), decreasing=T),]
    fulldf = fulldf[order(fulldf$p_cell),]
    fulldf$col = 1 * (fulldf$fdr < 0.05) * (2 - 1 * (fulldf$logFC < 0))
    fulldf$col_ind = 1 * (fulldf$fdr_ind < 0.05) * (2 - 1 * (fulldf$logFC < 0))
    ntop = 50
    labdf = rbind(head(fulldf[fulldf$col == 1,], ntop),
                  head(fulldf[fulldf$col == 2,], ntop))

    pcols = brewer.pal(12, 'Paired')
    gplot = ggplot(fulldf, aes(logFC, -log10(fdr), color=factor(col))) + 
        geom_point(cex=.1) + theme_pubr() + 
        scale_color_manual(values=c('grey80',pcols[1], pcols[5])) + 
        geom_text_repel(data=labdf, aes(logFC, -log10(fdr), color=factor(col), label=gene), size=2, max.overlaps=20) + 
        labs(x=paste('logFC ', ll$pathstr), y='log10 adj. p-value') + 
        geom_vline(xintercept=0, lty='dashed') + 
        geom_hline(yintercept=2, lty='dotted') +
        scale_y_continuous(expand=c(0,0)) + 
        theme(legend.position='none')
    pltprefix = paste0(imgpref, 'volcano_', prefstr, '_wilcoxon_frommat')
    saveGGplot(gplot, pltprefix, w=6, h=6)

    # Write out the results:
    # ----------------------
    # Write out number of significant genes:
    ctvec = c(celltype=celltype, subtype=subtype, 
              region=region, path=path, pctcut=pctcut, pcut=pcut)
    nsig = t(c(ctvec, table(fulldf$col)))
    write.table(nsig, wx.nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(wxtsv), quote=F, row.names=F, sep="\t")
    save(fulldf, nsig, file=wxrda)
} else {
    print("[STATUS] Wilcoxon regression output files already exist")
}

print("[STATUS] Finished script")
