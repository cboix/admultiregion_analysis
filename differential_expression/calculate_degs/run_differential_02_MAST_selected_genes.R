#!/usr/bin/R
# --------------------------------
# Run MAST on selected genes only:
# Updated: 02/23/22
# --------------------------------
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
geneset  = as.character(rundf[task,5])


# Set the directories:
# --------------------
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
gsdir = paste0(datadir , 'genesets/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, gsdir, regdir)
system(cmd)


# Set the run outputs: 
# ---------------------
prefstr = gsub("[()/]","_", paste(celltype, subtype, region, path, geneset, sep="_"))
masttsv = paste0(regdir, 'mast.', prefstr,'.tsv.gz')
mastrda = paste0(regdir, 'mast.', prefstr,'.rda')
mast.nsigfile = paste0(regdir, 'mast.', prefstr,'.nsig.tsv')
pcut = 0.05
print(prefstr)


# Read the relevant geneset:
# --------------------------
selgenes = scan(paste0(gsdir, geneset, '.txt'), 'c', quiet=TRUE)
print(paste("[STATUS] Read set:", geneset, 'with', length(selgenes), 'genes'))
print(head(selgenes))


# Check if already computed. If not, load data:
if (!file.exists(mastrda)){
    # Load dataset:
    commandArgs = function(x){ c(celltype, subtype, region)}
    source(paste0(bindir, 'multiRegion/load_difftl_data.R'))

    # Subset matrix to the selected genes:
    mat = subsetMatrixForDE(mat, pathdf, pctcells=pctcells, selgenes=selgenes)

    # Remove TH if running allregions + nft/plaq_n/plaq_d:
    if (path %in% c('nft','plaq_n','plaq_d')){
        pathdf = pathdf[pathdf$region != 'TH',]
        mat = mat[,pathdf$barcode]
    }

    # For normalize to log(TPM + 1):
    csm = colSums(mat)

    # Run RUV on the data at the individual level:
    # --------------------------------------------
    pathdf = runRUVpsbulk(pathdf, mat, path)

    # Put together the regression formula:
    # ------------------------------------
    ll = makeRegFormula(pathdf, path=path, nruv=10)

    print("[STATUS] Making MAST covariate matrix")
    fdata = data.frame(primerid=rownames(mat))
    covmat = data.frame(ll$mdx)
    covmat$wellKey = pathdf$barcode
    covmat$indiv = pathdf$projid
    covmat$ncounts = pathdf$ncounts
    # Due to way in which MAST takes arguments:
    covmat[path] = pathdf[,path]
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
        chunkrda = paste0(regdir, 'mast.', prefstr, '.', 
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
    fulldf$fdr = p.adjust(fulldf$p, 'fdr')
    fulldf = fulldf[order(fulldf$p),]
    fulldf = fulldf[!is.na(fulldf$coef),]
    fulldf$val = pctcells[fulldf$gene]
    fulldf$col = (fulldf$fdr < pcut) * (2 - 1 * (fulldf$coef < 0))
    fulldf = fulldf[fulldf$val > 0.2,]
    fulldf = fulldf[fulldf$gene != 'MALAT1',]
    ntop = 50
    labdf = rbind(head(fulldf[fulldf$col == 1,], ntop),
                  head(fulldf[fulldf$col == 2,], ntop))

    # Plot as a volcano plot:
    # -----------------------
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

    # Write out the results:
    # ----------------------
    # Write out number of significant genes:
    ctvec = c(celltype=celltype, subtype=subtype, 
              region=region, path=path, geneset=geneset, pcut=pcut)
    nsig = t(c(ctvec, table(fulldf$col)))
    write.table(nsig, mast.nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(masttsv), quote=F, row.names=F, sep="\t")
    save(fulldf, nsig, file=mastrda)
} else {
    print("[STATUS] MAST regression output files already exist")
}


print("[STATUS] Finished script")
