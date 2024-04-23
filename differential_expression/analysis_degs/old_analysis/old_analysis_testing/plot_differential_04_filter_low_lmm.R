#!/usr/bin/R
# -----------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Using Liang's Nebula package to run the fast NBLMM
# Basic overall + per-region interactions, for relevant subtypes
# Updated: 10/28/2020
# 
# Run with filters for low expr. genes, etc.
# NOTE: First, try to replicate Mic PFC results.
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(rhdf5)
library(nebula)
library(qvalue)
library(Matrix)

celltype = 'Mic_Immune'
chunk = 1
chunksize=4
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need celltype")
} else {        
    celltype = args[1]
    if (length(args) > 1){
        chunk=as.integer(args[2])
    } else { chunk=NULL }
}
print(celltype)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
pathlist = c('nft','plaq_d','plaq_n')

# Function to estimate end time of a loop:
est.endtime <- function(start, last, done, total){
    curr = proc.time() # Get current time
    # Calculate how much time has passed:
    last.step = (proc.time() - last)[3]
    elapsed = (proc.time() - start)[3]
    cat(paste0(round(last.step,1),'s\t'))
    cat(paste0(round(elapsed,1),'s\t'))
    # Estimate the remaining time:
    # Estimate the end time:
    each.time = elapsed / done
    est.left = (total - done) * each.time
    cat(paste0('Left: ', round(est.left,1),'s\t'))
    fin.time = Sys.time() + est.left
    cat(paste0('Est. ', format(fin.time, "%H:%M:%S"),'\n'))
    return(curr)
}

# ------------------------
# Load pathology measures:
# ------------------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

# Get the pathology mapped to each region:
regmap = c('AG','HC','PFC','MT','EC')
names(regmap) = c('ag','hip','mf','mt','ec')
pqdf = NULL
for (path in pathlist){
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    sdf = slong[,c('rind','value','region')]
    names(sdf)[2] = path
    if (is.null(pqdf)){
        pqdf = sdf 
    } else {
        pqdf = merge(pqdf, sdf)
    }
}
# pqdf = merge(pqdf, metadata[,c('rind','projid')])
# write.table(pqdf, file=paste0(datadir, 'region_pathology_scores.tsv'), quote=F, row.names=F, sep="\t")

# --------------------------------
# Load in the barcodes, pathology:
# --------------------------------
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    h5file = paste0(matdir, rawpref, '.majorcelltype.', celltype, '.hdf5')
    # SWMR file doesnt work with rhdf5:
    # h5file = paste0(matdir, rawpref, '.majorcelltype.', celltype, '.swmr.hdf5')
} else {
    matdir = paste0(datadir, 'matrices/')
    h5file = paste0(matdir, rawpref, '.majorcelltype.', celltype, '.hdf5')
}

# Load margin (for offset term):
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

# Extract metadata from hdf5:
h5f = H5Fopen(h5file)
genes = h5f$genes
barcodes = h5f$barcodes
H5Fclose(h5f)
ngenes = length(genes)

cellpref = paste0('multiRegion/rw_top_imputed_scores/', celltype)
fbc = scan(paste0(cellpref, '_barcodes.txt'), 'c', quiet=T)
impdf = data.frame(barcode=fbc)
impdf$imp_nft = as.numeric(scan(paste0(cellpref, '_imputed_scores_per_region_False_', 'nft','.txt'), 'c', quiet=T))
impdf$imp_plaq_n = as.numeric(scan(paste0(cellpref, '_imputed_scores_per_region_False_', 'plaq_n','.txt'), 'c', quiet=T))
impdf$imp_plaq_d = as.numeric(scan(paste0(cellpref, '_imputed_scores_per_region_False_', 'plaq_d','.txt'), 'c', quiet=T))

# ----------------------
# Make the model matrix:
# ----------------------
rownames(cellmeta) = cellmeta$barcode
submeta = cellmeta[barcodes,]
print(nrow(submeta[submeta$region != 'TH',]))
pathdf = merge(submeta, pqdf,all.x=TRUE)
pathdf = merge(pathdf, impdf, all.x=TRUE) # Imputed scores
print(nrow(pathdf))
pathdf = merge(pathdf, metadata[,c('projid','rind','age_death','msex','pmi', 'Apoe_e4', 'cogdx', 'niareagansc')])
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode
pathdf$cogdxad = 'NCI'
pathdf$cogdxad[pathdf$cogdx %in% c(4,5)] = 'AD'
pathdf$nrad = 'CTRL'
pathdf$nrad[pathdf$niareagansc %in% c(1,2)] = 'AD'
# Order properly
pathdf$nrad = factor(pathdf$nrad, levels=c('CTRL','AD'))
pathdf$cogdxad = factor(pathdf$cogdxad, levels=c('NCI','AD'))

# Split by ct:
split.var = 'hcelltype'
if (celltype %in% c('Exc','Inh', 'Vasc_Epithelia')){ 
    split.var = 'cell_type_high_resolution' 
} else if (celltype == 'Mic_Immune'){
    split.var = 'hcluster'
}
subtypes = unique(pathdf[[split.var]])
if (!is.null(chunk)){
    ind = ((chunk -1)* (chunksize) + 1):min(c(chunk * chunksize,length(subtypes)))
    subtypes = subtypes[ind]
}
print(subtypes)


# -----------------------------------------------------------
# Run the regression on a set of data:
# Req: count matrix, random effects, model matrix, and offset
# -----------------------------------------------------------
# Run by region here:
for (reg in c('PFC','AG','MT','EC','HC')){ #,'TH')){
    reg.pathdf = pathdf[pathdf$region == reg,]
    if (dbdir == '~/data/DEVTRAJ/db/') {
        matrda.file = paste0(matdir, prefix, '.subsetted.', reg, '.', celltype, '.matrix.Rda')
        if (!file.exists(matrda.file)){
            print(paste('[STATUS] Getting full matrix for', celltype, 'in', reg))
            bind = match(reg.pathdf$barcode, barcodes)
            # Pre-load data:
            h5f = H5Fopen(h5file)
            genes = h5f$genes
            bcs = h5f$barcodes
            h5d = h5f&"matrix"
            fullmat = h5d[,bind]
            H5Dclose(h5d)
            H5Fclose(h5f)
            rownames(fullmat) = genes
            colnames(fullmat) = bcs[bind]
            # Order as model matrix:
            fullmat = fullmat[,reg.pathdf$barcode]
            fullmat = Matrix(fullmat)
            gcout = gc()
            print(paste('[STATUS] Saving full matrix for', celltype, 'in', reg, 'in file:', matrda.file))
            save(fullmat, file=matrda.file)
        } else {
            load(matrda.file)
        }
    }
    for (subtype in subtypes){
        print(paste("[STATUS] Running on", subtype))
        for (path in c('nrad',paste0('imp_', pathlist), 'cogdxad',pathlist)) {
            # Select regions:
            print(reg)
            # reg = 'All'
            sub.pathdf = reg.pathdf[reg.pathdf[[split.var]] == subtype,] 
            sub.pathdf = sub.pathdf[!(is.na(sub.pathdf[[path]])),]
            if (reg != 'All') {
                sub.pathdf = sub.pathdf[sub.pathdf$region == reg,]
            } else {
                # NOTE: Need to filt regions for nft, etc.
                if (path %in% pathlist){
                    sub.pathdf = sub.pathdf[sub.pathdf$region == c('PFC','AG','MT','EC','HC'),]
                    sub.pathdf$region = factor(sub.pathdf$region, levels=c('PFC','AG','MT','EC','HC'))
                } else {
                    sub.pathdf$region = factor(sub.pathdf$region, levels=c('PFC','AG','MT','TH','EC','HC'))
                }
            }
            print(paste("Number of cells:", nrow(sub.pathdf)))
            ststr = gsub("/","_",gsub(" ","_", subtype))
            offset = marg[sub.pathdf$barcode, 'count']
            # TODO: Could add the ribo_frac and mt_frac as covariates?
            # Add if there:
            if (split.var != 'cell_type_high_resolution'){
                cvars ='+ age_death + msex + pmi + cell_type_high_resolution'
            } else { cvars = '+ age_death + msex + pmi' }
            if (reg == 'All'){
                mdx = model.matrix(asform(c('~',path, '+', path,'* region', cvars)), data=sub.pathdf)
                # m2 = model.matrix(asform(c('~',path, '+ region', cvars)), data=sub.pathdf)
            } else {
                mdx = model.matrix(asform(c('~',path, cvars)), data=sub.pathdf)
            }
            # Output files:
            regfile = paste0(regdir, prefix, '.nblmm_reg.', path, '.', reg, '.major.', celltype, '.minor.', ststr, '.Rda')
            regtsv  = paste0(regdir, prefix, '.nblmm_reg.', path, '.', reg, '.major.', celltype, '.minor.', ststr, '.tsv.gz')
            pathstr = path 
            if (path %in% c('nrad','cogdxad')){ pathstr = paste0(path,'AD') }
            leff = paste0('logFC_',pathstr)
            peff = paste0('p_',pathstr)

            if (!file.exists(regfile)){
                print(path)
                chunksize = 1000
                # chunksize = 20000
                nchunk = floor(ngenes / chunksize) + 1
                regdf = c()
                regdf2 = c()
                # Subset of locations to extract: 
                if (dbdir == '~/data/DEVTRAJ/db/') {
                    bind = match(sub.pathdf$barcode, colnames(fullmat))
                    subbcs = colnames(fullmat)[bind] 
                } else {
                    bind = match(sub.pathdf$barcode, barcodes)
                    subbcs = barcodes[bind] 
                }
                t0 = proc.time()
                t1 = t0
                for (i in 1:nchunk){
                    cat(i,'\t')
                    ind = (1 + (i-1) * chunksize):min(c(i * chunksize, ngenes)) 
                    # Open handle, extract genes we care about and close:
                    if (dbdir == '~/data/DEVTRAJ/db/') {
                        mat = fullmat[ind,bind]
                    } else {
                        h5f = H5Fopen(h5file)
                        genes = h5f$genes
                        bcs = h5f$barcodes
                        h5d = h5f&"matrix"
                        mat = h5d[ind,bind]
                        H5Dclose(h5d)
                        H5Fclose(h5f)
                        rownames(mat) = genes[ind]
                        colnames(mat) = bcs[bind]
                        # Order as model matrix:
                        mat = mat[,sub.pathdf$barcode]
                        mat = Matrix(mat)
                    }
                    gcout = gc()
                    rmean.mat = rowMeans(mat)
                    rmeanZ.mat = rowMeans(mat > 0)
                    # NOTE: Works fine as both PMM, NBGMM, and NBLMM, choose NBGMM - more conservative
                    re = nebula(mat, as.character(sub.pathdf$projid), pred=mdx, offset=log10(offset), model='NBGMM', cpc=0.05, verbose=F) 
                    rdf = re$summary
                    rdf$cpc = rmean.mat[rdf$gene]
                    rdf$zpc = rmeanZ.mat[rdf$gene]
                    # (As if no individual splits):
                    re2 = nebula(mat, rep(1, nrow(sub.pathdf)), pred=mdx, offset=log10(offset), model='NBGMM',cpc=0.05, verbose=F)
                    rdf2 = re2$summary
                    rdf2$cpc = rmean.mat[rdf2$gene]
                    rdf2$zpc = rmeanZ.mat[rdf2$gene]
                    if (nrow(rdf) > 0){ regdf = rbind(regdf, rdf) }
                    if (nrow(rdf2) > 0){ regdf2 = rbind(regdf2, rdf2) }
                    # Timing:
                    t1 = est.endtime(t0, t1, i, nchunk)
                }

                regdf$path = path
                regdf2$path = path
                regdf$q_path = qvalue(regdf[[peff]])$q
                regdf$log10q = -log10(regdf$q_path)
                regdf = regdf[order(regdf$log10q, decreasing=T),]
                regdf2$q_path = qvalue(regdf2[[peff]])$q
                regdf2$log10q = -log10(regdf2$q_path)
                regdf2 = regdf2[order(regdf2$log10q, decreasing=T),]
                # Save:
                save(regdf, regdf2, file=regfile)
                write.table(regdf, file=gzfile(regtsv), quote=F, row.names=F, sep="\t")
            } else {
                load(regfile)
            }
            print("With RANDOM EFFECTS")
            # print(head(regdf[,c(leff,'log10q','gene','cpc')], 50))
            # regdf[regdf$gene == 'APOE',]
            print(head(regdf[regdf[[leff]] > 0,c(leff,'log10q','gene','cpc')], 20))
            print(head(regdf[regdf[[leff]] < 0,c(leff,'log10q','gene','cpc')], 20))
            print("FIXED EFFECT ONLY")
            print(head(regdf2[regdf2[[leff]] > 0,c(leff,'log10q','gene','cpc')], 20))
            print(head(regdf2[regdf2[[leff]] < 0,c(leff,'log10q','gene','cpc')], 20))
        }
    }
}


print("Finished running regressions.")
