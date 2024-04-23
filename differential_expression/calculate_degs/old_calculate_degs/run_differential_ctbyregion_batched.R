#!/usr/bin/R
# -----------------------------------------------------------
# Use MAST + RE to run differential gene expression
# - Per region x per celltype subdivision
# Submit as chunks (2500 -> 8 chunks) per combination.
# TODO: Write separate script to batch and one to combine
# Updated: 01/21/2021
# OLD version of analysis
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)

library(MAST)
library(data.table)

print(version)

celltype = 'Mic_Immune'
subtype = 'Mic'
region = 'EC'
chunksize = 1250
chunk = 1
path = 'nrad'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype subtype region chunksize chunk ascertainment")
} else {        
    celltype = args[1]
    subtype = args[2]
    region = args[3]
    chunksize = as.integer(args[4])
    chunk = as.integer(args[5])
    path = args[6]
}

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# ------------------
# Load the metadata:
# ------------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

# Data directories:
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    # matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')

# Load in data:
ststr = gsub("/","_",gsub(" ","_", subtype))
matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                 celltype,'.',ststr,'.',region)
rdafile = paste0(matpref, '.rda')  # In Matrix format
# Load `mat` from rdafile:
load(rdafile)
print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
barcodes = colnames(mat)
genes = rownames(mat)
ngenes = nrow(mat)

# ------------------------------
# Load the appropriate metadata:
# ------------------------------
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

rownames(cellmeta) = cellmeta$barcode
pathdf = cellmeta[barcodes,]
pathdf = merge(pathdf, metadata[,c('projid','rind','age_death','msex','pmi', 'Apoe_e4', 'cogdx', 'niareagansc')])
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode
pathdf$nrad = 'CTRL'
pathdf$nrad[pathdf$niareagansc %in% c(1,2)] = 'AD'
pathdf$nrad = factor(pathdf$nrad, levels=c('CTRL','AD'))


# -----------------------------------------------------------
# Run the regression on a set of data:
# Req: count matrix, random effects, model matrix, and offset
# -----------------------------------------------------------
# Run on each subtype:
print(paste("[STATUS] Running regression on", subtype, 'in',region, 'on var', path, 'with ncell:', ncol(mat)))
ststr = gsub("/","_",gsub(" ","_", subtype))

# Build model:
pathdf$nGene = marg[pathdf$barcode, 'count']
mdx = model.matrix(asform(c('~',path, '+ age_death + msex + pmi + nGene')), data=pathdf)

# for (chunk in rev(1:16)){
# print(chunk)
# Output files:
regfile = paste0(regdir, prefix, '.mastlmm_reg.', path, '.', region, '.major.', celltype, '.minor.', ststr, '.', chunksize, '.', sprintf("%03d",chunk), '.Rda')
regtsv  = paste0(regdir, prefix, '.mastlmm_reg.', path, '.', region, '.major.', celltype, '.minor.', ststr, '.', chunksize, '.', sprintf("%03d",chunk), '.tsv.gz')
if (!file.exists(regfile)){
    t0 = proc.time()
    cat(subtype, '\t', path,'\t')
    mastdf = c()
    ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, ngenes)) 
    submat = mat[ind, pathdf$barcode]

    # --------------------
    # Prepare MAST object:
    # --------------------
    submat = as.matrix(submat)
    # sum(rownames(mdx) != colnames(submat))  # Check
    fdata = data.frame(primerid=rownames(submat))
    covmat = data.frame(mdx)
    covmat$wellKey = rownames(covmat)
    covmat$indiv = sub(".*-","",rownames(covmat))
    # Normalize to log(TPM + 1)
    norm = sweep(submat, 2, marg[colnames(submat),'count'] / 10000,'/')
    norm = log(norm + 1)
    sca = FromMatrix(exprsArray=norm, cData=covmat, fData=fdata, check_sanity=FALSE)

    # Filter down, only ones with 5%+ expr:
    freq_expressed = 0.05
    expressed_genes <- freq(sca) > freq_expressed
    sca <- sca[expressed_genes,]
    cat("Dim SCA:", dim(sca), '\t')
    if (nrow(sca) > 0){
        # ---------------------
        # Run MAST on raw data:
        # ---------------------
        # Calculate regression, run LRT:
        run.fe = FALSE
        if (run.fe){
            zlm.obj = suppressMessages(zlm(~ nGene + nradAD + msex + pmi, sca))
            summaryCond <- suppressMessages(summary(zlm.obj, doLRT='nradAD'))
            print(summaryCond, n=10)
            summaryDt <- summaryCond$datatable
            fcHurdle <- merge(summaryDt[contrast=='nradAD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                              summaryDt[contrast=='nradAD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
            FCTHRESHOLD = 0.05
            fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
            fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], data.table::as.data.table(mcols(sca)), by='primerid')
            data.table::setorder(fcHurdleSig, fdr)
            print(fcHurdleSig)
        }

        # Random effects model. Each of these steps will take a while:
        zlmCond <- zlm(~ nGene + nradAD + msex + pmi + (1 | indiv),
                       sca, method='glmer',ebayes = F, strictConvergence = FALSE)
        summaryCond <- summary(zlmCond, doLRT='nradAD')
        print(summaryCond, n=10)
        summaryDt <- summaryCond$datatable

        fcHurdle <- merge(summaryDt[summaryDt$contrast=='nradAD'
                          & summaryDt$component=='logFC', c(1,7,5,6,8)],
        summaryDt[summaryDt$contrast=='nradAD'
                  & summaryDt$component=='H', c(1,4)], by='primerid')
        fcHurdle <- merge(summaryDt[contrast=='nradAD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt[contrast=='nradAD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
        FCTHRESHOLD = 0.05
        fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
        fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], data.table::as.data.table(mcols(sca)), by='primerid')
        data.table::setorder(fcHurdle, fdr)
        data.table::setorder(fcHurdleSig, fdr)
        print(fcHurdleSig)

        # Write resuls out:
        regdf = data.frame(fcHurdle)
        save(regdf, summaryDt, file=regfile)
        write.table(regdf, file=gzfile(regtsv), quote=F, row.names=F, sep="\t")
    }
}

print("Finished running regressions.")
