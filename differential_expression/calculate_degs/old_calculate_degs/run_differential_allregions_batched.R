#!/usr/bin/R
# -----------------------------------------------------------
# Use MAST + RE to run differential gene expression
# - ACROSS ALL REGIONS x per celltype subdivision
# Submit as chunks (633 -> 30 chunks) per combination.
# Updated: 01/29/2021
# OLD version of analysis.
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
chunksize = 633
chunk = 1
path = 'nrad'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype subtype region chunksize chunk ascertainment")
} else {        
    celltype = args[1]
    subtype = args[2]
    chunksize = as.integer(args[3])
    chunk = as.integer(args[4])
    path = args[5]
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

# ------------------------------
# Load in data from all regions:
# ------------------------------
keep.reg = regions[regions != 'MB']
amat = c()
barcodes = c()
for (region in keep.reg){
    ststr = gsub("/","_",gsub(" ","_", subtype))
    matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                     celltype,'.',ststr,'.',region)
    rdafile = paste0(matpref, '.rda')  # In Matrix format
    # Load `mat` from rdafile:
    if (file.exists(rdafile)){
        load(rdafile)
        print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
        barcodes = c(barcodes, colnames(mat))
        genes = rownames(mat)
        ngenes = nrow(mat)
        amat = cbind(amat, mat)
    }
}
rm(mat)
gcout = gc()

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
pathdf$cogdxad = 'CTRL'
pathdf$cogdxad[pathdf$cogdx %in% c(4,5)] = 'AD'
pathdf$cogdxad = factor(pathdf$cogdxad, levels=c('CTRL','AD'))
pathdf$nrad = 'CTRL'
pathdf$nrad[pathdf$niareagansc %in% c(1,2)] = 'AD'
pathdf$nrad = factor(pathdf$nrad, levels=c('CTRL','AD'))

if (path %in% c('nft','plaq_d','plaq_n')){
    # Get the pathology mapped to each region:
    regmap = c('AG','HC','PFC','MT','EC')
    names(regmap) = c('ag','hip','mf','mt','ec')
    pqdf = NULL
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    pqdf = slong[,c('rind','value','region')]
    names(pqdf)[2] = path
    pathdf = merge(pathdf, pqdf) # Drops out TH region

    # Subset to only kept cells:
    barcodes = pathdf$barcode
    rownames(pathdf) = pathdf$barcode
    amat = amat[,barcodes]
    gcout = gc()
}

# -----------------------------------------------------------
# Run the regression on a set of data:
# Req: count matrix, random effects, model matrix, and offset
# -----------------------------------------------------------
# Run on each subtype:
print(paste("[STATUS] Running regression on", subtype, 'across all regions, on var', path, 'with ncell:', ncol(amat)))
ststr = gsub("/","_",gsub(" ","_", subtype))

# Build model:
pathdf$nGene = marg[pathdf$barcode, 'count']
mdx = model.matrix(asform(c('~',path, '+ age_death + msex + pmi + region + nGene')), data=pathdf)

# Output files:
# for (chunk in 1:30){
regfile = paste0(regdir, prefix, '.mastlmm_reg.', path, '.allreg.major.', celltype, '.minor.', ststr, '.', chunksize, '.', sprintf("%03d",chunk), '.Rda')
regtsv  = paste0(regdir, prefix, '.mastlmm_reg.', path, '.allreg.major.', celltype, '.minor.', ststr, '.', chunksize, '.', sprintf("%03d",chunk), '.tsv.gz')
if (!file.exists(regfile)){
    t0 = proc.time()
    cat(subtype, '\t', path,'\t')
    mastdf = c()
    ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, ngenes)) 
    submat = amat[ind,pathdf$barcode]
    # rm(amat)
    gcout = gc()

    # --------------------
    # Prepare MAST object:
    # --------------------
    submat = as.matrix(submat)
    # sum(rownames(mdx) != colnames(submat))  # Check
    fdata = data.frame(primerid=rownames(submat))
    covmat = data.frame(mdx)
    covmat$wellKey = rownames(covmat)
    covmat$indiv = sub(".*-","",rownames(covmat))
    covmat$region = sub("_.*","",rownames(covmat))
    # Normalize to log(TPM + 1)
    norm = sweep(submat, 2, marg[colnames(submat),'count'] / 10000,'/')
    norm = log(norm + 1)
    sca = FromMatrix(exprsArray=norm, cData=covmat, fData=fdata, check_sanity=FALSE)

    # Filter down, only ones with 5%+ expr:
    freq_expressed = 0.05
    expressed_genes <- freq(sca) > freq_expressed
    sca <- sca[expressed_genes,]
    cat("Dim SCA:", dim(sca), '\t')

    # Run regressions:
    if (nrow(sca) > 0){
        # ---------------------
        # Run MAST on raw data:
        # ---------------------
        # Calculate regression, run LRT:
        run.fe = FALSE
        if (run.fe){
            zlm.obj = zlm(~ nGene + nradAD + msex + pmi + nradAD * region, sca)
            summaryCond <- summary(zlm.obj, doLRT='nradAD')
            print(summaryCond, n=10)
            summaryDt <- summaryCond$datatable
            fcHurdle <- merge(summaryDt[contrast=='nradAD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                              summaryDt[contrast=='nradAD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
            FCTHRESHOLD = 0.05
            fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
            fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], data.table::as.data.table(mcols(sca)), by='primerid')
            data.table::setorder(fcHurdle, fdr)
            data.table::setorder(fcHurdleSig, fdr)
            print(fcHurdleSig)
        }

        # Random effects model. Each of these steps will take a while:
        # zlmCond <- zlm(~ nGene + nradAD + msex + pmi + (1 | indiv), # Should model region
        # zlmCond <- zlm(~ nGene + nradAD + msex + pmi + nradAD * region + (1 | indiv), # Not appropriate
        if (path == 'nrad'){
            pathstr = 'nradAD'
        } else if (path == 'cogdxad'){
            pathstr = 'cogdxadAD'
        } else {
            pathstr = path
        }
        zlmCond <- zlm(as.formula(paste0('~ nGene +',pathstr, '+ msex + pmi + region + (1 | indiv)')), # Appropriate model
                       sca, method='glmer',ebayes = F, strictConvergence = FALSE)
        summaryCond <- summary(zlmCond, doLRT=pathstr)
        print(summaryCond, n=10)
        summaryDt <- summaryCond$datatable

        fcHurdle <- merge(summaryDt[summaryDt$contrast==pathstr
                          & summaryDt$component=='logFC', c(1,7,5,6,8)],
        summaryDt[summaryDt$contrast==pathstr
                  & summaryDt$component=='H', c(1,4)], by='primerid')
        fcHurdle <- merge(summaryDt[contrast==pathstr & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt[contrast==pathstr & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
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
# }

print("Finished running regressions.")
