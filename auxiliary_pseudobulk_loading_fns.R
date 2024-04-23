#!/usr/bin/R
# ---------------------------------------------
# Functions to load data at a pseudobulk level:
# Averaged at level of subtype + individual + region
# Updated 12/01/21
# ---------------------------------------------
library(cbrbase)
library(Matrix)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


# Data directories where RDA matrices are located:
# ------------------------------------------------
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')


# Read in the margins for normalizing the matrix:
# -----------------------------------------------
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
margin = as.numeric(scan(gzfile(margfile), 'c', quiet=TRUE))
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
names(margin) = mbcs


# Functions to load and normalize each matrix:
# -------------------------------------------
load_single_matrix_rda = function(celltype, subtype, region, normalize=FALSE){
    ctstr = gsub("/","_",gsub(" ","_", celltype))
    ststr = gsub("/","_",gsub(" ","_", subtype))
    rdafile = paste0(mtxdir, rawpref,'.majorcelltype.', 
                     ctstr,'.',ststr,'.',region, '.rda')
    if (file.exists(rdafile)){
        load(rdafile)
        if ((nrow(mat) > 1) & (ncol(mat) > 1)){
            # Normalize:
            if (normalize){
                bcs = colnames(mat)
                amarg = margin[bcs]
                fact = amarg / 10000  # Normalize to uniform 10k
                mat@x <- mat@x / rep.int(fact, diff(mat@p))
                gc()
            }
        } else { mat = NULL }
    } else {
        mat = NULL
        print(paste0("No file for: ", rdafile))
    }
    if (!is.null(mat)){
        print(paste("[STATUS] Loaded", subtype, 'in',
                    region, 'with', ncol(mat), 'cells'))
    }
    return(mat)
}


load_pseudobulk_subtype = function(celltype, subtype, region, clscol, byind=TRUE, normalize=TRUE){
    # Load `mat` from rdafile:
    mat = load_single_matrix_rda(celltype, subtype, region, normalize=normalize)

    if (!is.null(mat)){
        # Aggregate to the level of individual x cell type x region:
        submeta = cellmeta[cellmeta$major.celltype == celltype,]
        rownames(submeta) = submeta$barcode
        smeta = submeta[colnames(mat), c('projid', 'barcode','region', clscol)]
        fmt.cthr = gsub("/","_",gsub(" ","_", smeta[[clscol]]))
        if (byind){
            smeta$ptype = with(smeta, paste0(projid, '_', fmt.cthr, '_', region))
        } else {
            smeta$ptype = with(smeta, paste0(fmt.cthr, '_', region))
        }
        ptlist = sort(unique(smeta$ptype))
        # DON'T NORMALIZE IF JUST ADDING UP COUNTS
        tform = make.tform(smeta$ptype, u=sort(unique(smeta$ptype)), norm=normalize)
        pmat = mat %*% tform # norm

        # Return with unique metadata:
        if (byind){
            uqmeta = agg.rename(asform(c('barcode ~ projid + region +',
                                         clscol, '+ ptype')), smeta, length, 'ncell')
        } else {
            uqmeta = agg.rename(asform(c('barcode ~ region +',
                                         clscol, '+ ptype')), smeta, length, 'ncell')
        }
        rownames(uqmeta) = uqmeta$ptype
        uqmeta = uqmeta[ptlist,]
        if (!is.null(mat)){
            print(paste('[STATUS] Reduced to', ncol(pmat), 'batches'))
        }
        gc()
        psdata = list('mat'=pmat, 'meta'=uqmeta)
    } else { psdata = list('mat'=NULL, 'meta'=NULL) }
    return(psdata)
}



# Load data, either one or multiple subtypes, or across regions, etc.
# -------------------------------------------------------------------
load_pseudobulk_dataset = function(celltype, subtype, region, byind=TRUE, normalize=TRUE){
    if (celltype == 'All'){
        celltype = unique(cellmeta$major.celltype)
        if (length(subtype) == 0){
            subtype = c(unique(cellmeta$cell_type_high_resolution), unique(cellmeta$minor.celltype))
        }
    } 
    datasets = expand.grid(ct=celltype, st=subtype, reg=region)
    print(datasets)
    ps.data = list(mat = c(), meta = c())
    clscol = 'cell_type_high_resolution'
    for (i in 1:nrow(datasets)){
        ll = load_pseudobulk_subtype(celltype=datasets$ct[i], 
                                     subtype=datasets$st[i], 
                                     region=datasets$reg[i], 
                                     clscol=clscol, byind=byind,
                                     normalize=normalize)
        if (!is.null(ll$mat)){ 
            ps.data$mat = cbind(ps.data$mat, ll$mat)
            ps.data$meta = rbind(ps.data$meta, ll$meta)
        }
    }
    return(ps.data)
}


# Load data, either one or multiple subtypes, or across regions, etc.
# -------------------------------------------------------------------
load_full_dataset = function(celltype, subtype, region, normalize=FALSE){
    datasets = expand.grid(ct=celltype, st=subtype, reg=region)
    # print(datasets)
    mat = c()
    for (i in 1:nrow(datasets)){
        imat = load_single_matrix_rda(celltype=datasets$ct[i],
                                      subtype=datasets$st[i],
                                      region=datasets$reg[i],
                                      normalize=normalize)
        if (!is.null(imat)){ mat = cbind(mat, imat) }
    }
    rm(imat)
    gcout = gc()
    return(mat)
}


# Aggregate pseudobulk data at the region level:
# ----------------------------------------------
aggregate_psbulk_samplelevel = function(ps.data){
    ps.data$meta$pr = with(ps.data$meta, paste0(projid, '_', region))
    prlist = sort(unique(ps.data$meta$pr))
    tform = make.tform(ps.data$meta$pr, u=prlist, norm=T)
    tform = sweep(tform, 1, ps.data$meta$ncell, '*')
    tform = sweep(tform, 2, apply(tform, 2, sum), '/')
    ps.data$mat = ps.data$mat %*% tform # norm

    # Fix metadata:
    ps.data$meta = aggregate(ncell ~ pr + projid + region, ps.data$meta, sum)
    ps.data$mat = ps.data$mat[, ps.data$meta$pr]
    return(ps.data)
}


# Make individual-level pseudobulk annotation function:
# -----------------------------------------------------
make_ind_annotation = function(indorder, ux=1.5, which='column'){
    require(ComplexHeatmap)
    require(circlize)
    umeta = unique(metadata[metadata$region == 'PFC',
        c('projid','nrad','cogdxad','cogdx', 'niareagansc', 'age_death',
            'msex', 'braaksc', 'Apoe_e4','gpath','tangles','amyloid', 'cogn_global_lv')])
    rownames(umeta) = umeta$projid
    umeta = umeta[indorder,]

    age.col_fun = colorRamp2(range(umeta$age_death), c("white", "slateblue")) 
    pmi.col_fun = colorRamp2(c(2, 15), c("white", "indianred")) 
    gpath.col_fun = colorRamp2(c(0, max(umeta$gpath)), c("white", "indianred")) 
    gcog.col_fun = colorRamp2(c(min(umeta$cogn_global_lv), max(umeta$cogn_global_lv)),
        c("red", 'white')) 
    # mat.col_fun = colorRamp2(c(0, max(pmat, na.rm=T)), c("white", "blue")) 

    # Make metadata annotation:
    ha = HeatmapAnnotation(
        annotation_name_gp = gpar(fontsize=5),
        simple_anno_size = unit(ux, 'mm'),
        gap=0,
        Sex=ifelse(umeta$msex == 0, 'female','male'), 
        # PMI=umeta$pmi,
        Age=umeta$age_death,
        Apoe_e4=umeta$Apoe_e4,
        GPath=umeta$gpath,
        Braak=umeta$braaksc,
        AD=umeta$niareagansc,
        Cognition=umeta$cogdx,
        GlobalCog=umeta$cogn_global_lv,
        col=list(AD=colvals[['niareagansc']],
            Apoe_e4=c('no'='grey95','yes'='grey70'),
            Age=age.col_fun,
            GPath=gpath.col_fun,
            GlobalCog=gcog.col_fun,
            # PMI=pmi.col_fun,
            Braak=colvals[['braaksc']],
            Cognition=colvals[['cogdx']],
            Sex=colvals[['sex']]), which=which)
    return(ha)
}


