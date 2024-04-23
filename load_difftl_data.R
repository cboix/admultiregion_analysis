#!/usr/bin/R
# -----------------------
# Load differential data:
# -----------------------
library(cbrbase)
library(Matrix)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype subtype region")
} else if (args != 'local'){
    celltype = args[1]
    subtype = args[2]
    region = args[3]
}

# -----------------
# Data directories:
# -----------------
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')

# -------------
# Load in data:
# -------------
load_subtype = function(celltype, st, region){
    ststr = gsub("/","_",gsub(" ","_", st))
    if (region == 'allregions' | region == 'neocortex'){
        amat = NULL
        if (region == 'neocortex'){ 
            selreg = c('AG','MT','PFC')
        } else {
            selreg = reg.nomb
        }
        for (reg in selreg){
            matpref = paste0(mtxdir, rawpref,'.majorcelltype.', 
                             celltype,'.',ststr,'.',reg)
            rdafile = paste0(matpref, '.rda')  # In Matrix format
            if (file.exists(rdafile)){
                load(rdafile)
                amat = cbind(amat, mat)
            } else {
                print(paste0("No file for: ", rdafile))
            }
        }
        mat = amat
        rm(amat)
        gcout = gc()
    } else {
        # Load `mat` from rdafile:
        matpref = paste0(mtxdir, rawpref,'.majorcelltype.', 
                         celltype,'.',ststr,'.',region)
        rdafile = paste0(matpref, '.rda')  # In Matrix format
        if (file.exists(rdafile)){
            load(rdafile)
        } else {
            mat = NULL
            print(paste0("No file for: ", rdafile))
        }
    }
    if (!is.null(mat)){
        print(paste("[STATUS] Loaded", st, 'in',
                    region, 'with', ncol(mat), 'cells'))
    }
    return(mat)
}

# Load data, either one or multiple subtypes:
if (length(subtype) > 1){
    mat = c()
    for (st in subtype){
        print(st)
        stmat = load_subtype(celltype, st, region)
        if (!is.null(stmat)){ mat = cbind(mat, stmat) }
    }
} else {
    mat = load_subtype(celltype, subtype, region)
}
barcodes = colnames(mat)
genes = rownames(mat)
ngenes = nrow(mat)

# NOTE: Already subsetted to protein coding genes:
# pc_genes = anno$symbol[anno$type == 'protein_coding']
# genes = genes[genes %in% pc_genes]

# Properties of cells:
nc = colSums(mat)
ng = colSums(mat > 0)
cpg = nc / ng
nmt = colSums(mat[c(grep("^MT-", genes), grep("^MTRNR", genes)),])
pct_mt = nmt / nc

# Properties of genes:
pctcells = rowSums(mat > 0) / ncol(mat)


# ------------------------------
# Load the appropriate metadata:
# ------------------------------
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode


# -------------------------------------------------------
# Merge in cell, sample, and individual-level covariates:
# -------------------------------------------------------
extsi_tsv = 'Annotation/extended_si_ace_vars.tsv'
ext.si.vars = read.delim(extsi_tsv, header=F)[,1]
advars = c('cogdx', 'nrad', 'cogdxad', 'braaksc.early','braaksc.ad')
covars = c('age_death','msex','pmi', 'Apoe_e4') 
ext.covars = c('soc_net_bl', 'social_isolation_avg', 
    'social_isolation_lv', ext.si.vars)

ext.meta_csv = 'Annotation/ROSMAP_All_Pts_Extended_Basic_Long_652_03-23-2022.csv'
ext.metadata = read.delim(ext.meta_csv, header=T, sep=",")

# Make design dataframe:
rownames(cellmeta) = cellmeta$barcode
cols = c('barcode','rind','region','projid','major.celltype','cell_type_high_resolution')
pathdf = cellmeta[barcodes,cols]

# Merge all individual and sample-level variables:
pathdf = merge(pathdf, metadata[,c('projid','rind', covars, advars)])
pathdf = merge(pathdf, ext.metadata[,c('projid',ext.covars)])
pathdf = merge(pathdf, pqdf, all.x=TRUE)
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode

# Add cell-level covariates:
pathdf$ngenes = ng[pathdf$barcode]
pathdf$ncounts = nc[pathdf$barcode]
pathdf$cpg = cpg[pathdf$barcode]
pathdf$pctMT = pct_mt[pathdf$barcode]

