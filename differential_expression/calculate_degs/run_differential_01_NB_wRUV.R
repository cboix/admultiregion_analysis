#!/usr/bin/R
# --------------------------------------------------------------
# Nebula as in minimal example, but using the precomp. matrices:
# Updated: 10/31/23 to add interaction
# --------------------------------------------------------------
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
# runfile = paste0(sdbdir, 'DEG_multiRegion_SI_ACE_runlist.tsv')
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
outtsv = paste0(regdir, 'nebula_ruv.', prefstr,'.tsv.gz')
outrda = paste0(regdir, 'nebula_ruv.', prefstr,'.rda')
nsigfile = paste0(regdir, 'nebula_ruv.', prefstr,'.nsig.tsv')
pcut = 0.05
print(prefstr)


# Check if already computed:
if (!file.exists(outrda)){
    # Run differential expression using Nebula + RUV:
    # -----------------------------------------------
    # Load data and subset matrix:
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

    # Run RUV on the matrix:
    # ----------------------
    pathdf = runRUVpsbulk(pathdf, mat, path)
    if (length(unique(pathdf$region)) > 1){
        indvar = 'rp'
    } else { indvar = 'projid' }

    # Run NEBULA with the RUV results and other vars:
    # -----------------------------------------------
    if (int){ int.var = 'nrad' } else { int.var = NULL }
    ll = makeRegFormula(pathdf, path=path, nruv=10, int.var=int.var)

    chunksize = 500
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
            re = try(nebula(submat, as.character(pathdf[[indvar]]), 
                    pred=ll$mdx, offset=offset, model='PMM'))
            # If the original chunked regression fails (singular), run each gene separately.
            if (class(re) == 'try-error'){
                rdf = c()
                for (j in 1:length(ind)){
                    print(j)
                    x = submat[j,, drop=F]
                    re = try(nebula(x, as.character(pathdf[[indvar]]), 
                            pred=ll$mdx, offset=offset, model='PMM'))
                    while (class(re) == 'try-error'){  # TODO: EVAL WHILE
                        x[x > (max(x) - 5)] = (max(x) - 5)
                        re = try(nebula(x, as.character(pathdf[[indvar]]), 
                                pred=ll$mdx, offset=offset, model='PMM'))
                    }
                    rdf = rbind(rdf, re$summary)
                }
            } else { rdf = re$summary }
            # Process + save the chunk:
            chunkdf = rdf[order(rdf[[ll$peff[1]]]),c('gene', ll$leff, ll$peff)]
            if (length(ll$eff) > 1){
                names(chunkdf) = c('gene','logFC','p')
            } else {
                names(chunkdf) = c('gene',ll$leff, ll$peff)
            }
            save(chunkdf, file=chunkrda)
        } else {
            load(chunkrda)
        }
        names(chunkdf) = paste0("X", 1:ncol(chunkdf))
        fulldf = rbind(fulldf, chunkdf)
        t1 = est.endtime(t0, t1, chunk, nchunk)
    }

    if (length(ll$eff) > 1){
        names(fulldf) = c('gene','logFC','p')
    } else {
        names(fulldf) = c('gene',ll$leff, ll$peff)
    }

    # Process the full results:
    # -------------------------
    # NOTE: if leff is multiple:
    fulldf$pc = pctcells[fulldf$gene]
    if (length(ll$peff) > 1){
        fulldf = fulldf[order(fulldf[[ll$peff[1]]]),]
        for (pvar in ll$peff){
            lvar = sub("^p_", "logFC_", pvar)
            qvar = sub("^p_", "q_", pvar)
            cvar = sub("^p_", "col_", pvar)
            fulldf[[qvar]] = qvalue(fulldf[[pvar]])$q
            fulldf[[cvar]] = 1 * (fulldf[[qvar]] < pcut) * (2 - 1 * (fulldf[[lvar]] < 0))
        }
    } else {
        fulldf = fulldf[order(fulldf$p),]
        fulldf$padj = p.adjust(fulldf$p, 'fdr')
        fulldf$q = qvalue(fulldf$p)$q
        fulldf$col = 1 * (fulldf$q < pcut) * (2 - 1 * (fulldf$logFC < 0))
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
    }

    # Write out the results:
    # ----------------------
    # Write out nsig: 
    if (length(ll$peff) > 1){
        nsig = c()
        for (pvar in ll$peff){
            cvar = sub("^p_", "col_", pvar)
            ndf = table(fulldf[,cvar])
            ndf = as.data.frame(ndf)
            ndf = spread(ndf, Var1, Freq)
            ndf$eff = sub("^p_", "", pvar)
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
    write.table(nsig, nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(outtsv), quote=F, row.names=F, sep="\t")
    save(fulldf, nsig, file=outrda)
} else {
    print("[STATUS] Nebula regression output files already exist")
}

print("[STATUS] Finished script")
