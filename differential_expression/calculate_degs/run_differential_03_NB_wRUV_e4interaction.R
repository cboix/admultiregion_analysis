#!/usr/bin/R
# ------------------------------------------------------
# Run Nebula on e4 interaction with selected genes only:
# - Used for Oli-M7 + Opc-M9 test
# Updated: 06/27/23
# ------------------------------------------------------
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
options(width=170)

# Functions for DE:
cbindir = paste0(sbindir, 'differential_expression/calculate_degs/')
source(paste0(cbindir, 'auxiliary_differential_functions.R'))


# Set the directories:
# --------------------
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)


# Arguments (for now just Oli/Opc M7/M9)
# --------------------------------------
celltype = 'Oli'

modmap = c('Oli'=7, 'Opc'=9)
subtype = celltype
module = modmap[celltype]
region = 'allregions'
# path = 'plaq_n'
path = 'cogdxad'
mtag = paste0(celltype, '-M', module)
run.interaction = FALSE


# Load modules:
# -------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Set the run outputs: 
# ---------------------
prefstr = gsub("[()/]","_", paste(celltype, subtype, region, path, module, sep="_"))
suffstr = ifelse(run.interaction, 'module_e4', 'module_e4noint')
respref = paste0(regdir, 'nebula_ruv.', suffstr, '.', prefstr)
outtsv = paste0(respref, '.tsv.gz')
outrda = paste0(respref, '.rda')
nsigfile = paste0(respref, '.nsig.tsv')
pcut = 0.05
print(prefstr)


# Read the relevant gene set from the modules:
# --------------------------------------------
use.core = FALSE  # Take all associated genes
coremap = cmlist[[celltype]]
genemap = gmlist[[celltype]]
if (use.core){ usemap = coremap } else { usemap = genemap }
selgenes = sort(names(usemap)[usemap == module])

print(paste("[STATUS] Read module", mtag, 'with', length(selgenes), 'genes'))
print(head(selgenes))


# Check if already computed:
if (!file.exists(outrda)){
    # Run differential expression using Nebula + RUV:
    # -----------------------------------------------
    # Load data and subset matrix:
    commandArgs = function(x){ c(celltype, subtype, region)}
    source(paste0(bindir, 'multiRegion/load_difftl_data.R'))

    # Subset matrix to genes in 20%+ of cells:
    mat = subsetMatrixForDE(mat, pathdf, pctcells=pctcells, selgenes=selgenes)
    print(paste("[STATUS] Subset an original", length(selgenes), 'genes to', nrow(mat), 'genes'))

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
    int.type = ifelse(run.interaction, '*', '+')
    ll = makeRegFormula(pathdf, path=paste(path, int.type, 'Apoe_e4'), nruv=10)

    # Redefine pathstr for interactions, etc:
    pathstr = ifelse(path == 'cogdxad', paste0(path,'AD'), path)
    e4str = 'Apoe_e4yes'
    intstr = paste0(pathstr, ':', e4str)
    if (run.interaction){
        kept.eff = c(intstr, e4str, pathstr)
    } else {
        kept.eff = c(e4str, pathstr)
    }
    leff = paste0('logFC_', kept.eff)
    peff = paste0('p_', kept.eff)
    kept.cols = c(leff, peff)

    # Run NEBULA:
    chunksize = ifelse(nrow(mat) > 750, 500, nrow(mat))
    nchunk = ceiling(nrow(mat) / chunksize)
    fulldf = c()
    offset = log10(pathdf$ncounts)
    t0 = proc.time()
    t1 = t0
    for (chunk in 1:nchunk){
        ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
        print(paste(chunk, length(ind)))
        chunkrda = paste0(regdir, 'nebula_ruv.', suffstr, '.', prefstr,'.', 
                          chunksize, '_', chunk, '.rda')
        if (!file.exists(chunkrda)){
            submat = mat[ind,]
            re = nebula(submat, as.character(pathdf[[indvar]]), 
                        pred=ll$mdx, offset=offset, model='PMM') 
            rdf = re$summary
            chunkdf = rdf[order(rdf[[peff[1]]]),c('gene', kept.cols)]
            save(chunkdf, file=chunkrda)
        } else {
            load(chunkrda)
        }
        fulldf = rbind(fulldf, chunkdf)
        t1 = est.endtime(t0, t1, chunk, nchunk)
    }

    # Process the full results:
    # -------------------------
    fulldf$celltype = celltype
    fulldf$subtype = subtype
    fulldf$module = module
    fulldf$path = path
    fulldf$region = region
    fulldf$pc = pctcells[fulldf$gene]

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(outtsv), quote=F, row.names=F, sep="\t")
    save(fulldf, file=outrda)
} else {
    print("[STATUS] Nebula regression output files already exist")
}

print("[STATUS] Finished script")
