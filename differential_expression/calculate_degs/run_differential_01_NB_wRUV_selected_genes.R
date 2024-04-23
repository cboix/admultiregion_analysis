#!/usr/bin/R
# ----------------------------------
# Run Nebula on selected genes only:
# Updated: 02/23/22
# ----------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(Matrix)

# For nebula: 
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
regdir = paste0(datadir,'dereg/')
gsdir = paste0(datadir , 'genesets/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir, gsdir)
system(cmd)


# Set the run outputs: 
# ---------------------
prefstr = gsub("[()/]","_", paste(celltype, subtype, region, path, geneset, sep="_"))
outtsv = paste0(regdir, 'nebula_ruv.', prefstr,'.tsv.gz')
outrda = paste0(regdir, 'nebula_ruv.', prefstr,'.rda')
nsigfile = paste0(regdir, 'nebula_ruv.', prefstr,'.nsig.tsv')
pcut = 0.05
print(prefstr)

# Read the relevant geneset:
# --------------------------
if (geneset == 'FA'){
    selgenes = c('ADCY8', 'HILPDA', 'PFKP', 'GFAP', 'AQP4')
} else {
    selgenes = scan(paste0(gsdir, geneset, '.txt'), 'c', quiet=TRUE)
}
print(paste("[STATUS] Read set:", geneset, 'with', length(selgenes), 'genes'))
print(head(selgenes))


# Check if already computed:
if (!file.exists(outrda)){
    # Run differential expression using Nebula + RUV:
    # -----------------------------------------------
    # Load data and subset matrix:
    commandArgs = function(x){ c(celltype, subtype, region)}
    source(paste0(bindir, 'multiRegion/load_difftl_data.R'))

    # Subset matrix to selected genes:
    mat = subsetMatrixForDE(mat, pathdf, pctcells=pctcells, selgenes=selgenes)

    # Remove TH if running allregions + nft/plaq_n/plaq_d:
    if (path %in% c('nft','plaq_n','plaq_d')){
        pathdf = pathdf[pathdf$region != 'TH',]
        mat = mat[,pathdf$barcode]
    }

    # Put together the regression formula:
    # ------------------------------------
    pathdf = runRUVpsbulk(pathdf, mat, path)
    if (length(unique(pathdf$region)) > 1){
        indvar = 'rp'
    } else { indvar = 'projid' }

    # Run NEBULA with the RUV results and other vars:
    # -----------------------------------------------
    if (geneset == 'FA'){
        ll = makeRegFormula(pathdf, path=path, nruv=2)
    } else {
        ll = makeRegFormula(pathdf, path=path, nruv=10)
    }

    chunksize=500
    nchunk = ceiling(nrow(mat) / chunksize)
    fulldf = c()
    offset = log10(pathdf$ncounts)
    t0 = proc.time()
    t1 = t0
    for (chunk in 1:nchunk){
        ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
        print(paste(chunk, length(ind)))
        chunkrda = paste0(regdir, 'nebula_ruv.', prefstr,'.', 
                          chunksize, '_', chunk, '.rda')
        if (!file.exists(chunkrda)){
            submat = mat[ind,]
            re = nebula(submat, as.character(pathdf[[indvar]]), 
                        pred=ll$mdx, offset=offset, model='PMM') 
            rdf = re$summary
            chunkdf = rdf[order(rdf[[ll$peff]]),c('gene', ll$leff, ll$peff)]
            names(chunkdf) = c('gene','logFC','p')
            save(chunkdf, file=chunkrda)
        } else {
            load(chunkrda)
        }
        fulldf = rbind(fulldf, chunkdf)
        t1 = est.endtime(t0, t1, chunk, nchunk)
    }

    # Process the full results:
    # -------------------------
    fulldf = fulldf[order(fulldf$p),]
    fulldf$padj = p.adjust(fulldf$p, 'fdr')
    fulldf$q = qvalue(fulldf$p)$q
    fulldf$col = 1 * (fulldf$q < pcut) * (2 - 1 * (fulldf$logFC < 0))
    fulldf$pc = pctcells[fulldf$gene]
    labdf = fulldf[fulldf$col != 0,]

    # Plot as a volcano plot:
    # -----------------------
    pcols = brewer.pal(12,'Paired')
    gplot = ggplot(fulldf, aes(logFC, -log10(p), color=factor(col))) + 
        geom_point(cex=.25) + 
        geom_text_repel(data=labdf, aes(logFC, -log10(p), label=gene, color=factor(col)), size=2, max.overlaps=20) + 
        scale_color_manual(values=c('grey80',pcols[1],pcols[5])) + 
        scale_y_continuous(expand=c(0,0)) + 
        theme_pubr() + theme(legend.position='none')
    pltprefix = paste0(imgpref, 'volcano_',prefstr,'_nebula_frommat')
    saveGGplot(gplot, pltprefix, w=6, h=6)

    # Write out the results:
    # ----------------------
    # Write out nsig: 
    ctvec = c(celltype=celltype, subtype=subtype, 
              region=region, path=path, geneset=geneset, pcut=pcut)
    nsig = t(c(ctvec, table(fulldf$col)))
    write.table(nsig, nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(outtsv), quote=F, row.names=F, sep="\t")
    save(fulldf, nsig, file=outrda)
} else {
    print("[STATUS] Nebula regression output files already exist")
}

print("[STATUS] Finished script")
