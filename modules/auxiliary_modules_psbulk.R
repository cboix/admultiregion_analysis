#!/usr/bin/R
# --------------------------------------------------
# Functions to process to individual / region level:
# --------------------------------------------------
library(cbrbase)


aggregatePsbulkIndRegion = function(plt.data){
    # Aggregate at the projid x region level (should be single-set already):
    plt.data$meta$pr = with(plt.data$meta, paste(projid, region, sep="_"))
    pr.sets = unique(plt.data$meta$pr)
    tform = make.tform(plt.data$meta$pr, u=pr.sets, norm=FALSE)

    # Normalize by number of cells:
    tform = sweep(tform, 1, plt.data$meta$ncell, '*')
    tform = sweep(tform, 2, apply(tform, 2, sum), '/')
    plt.data$mat = plt.data$mat %*% tform

    # Reduce the metadata to this level of granularity:
    plt.data$meta$cell_type_high_resolution = NULL
    plt.data$meta$ptype = NULL
    plt.data$meta$set = NULL
    form = as.formula(paste('ncell ~', paste(colnames(plt.data$meta), collapse='+')))
    plt.data$meta = aggregate(form, plt.data$meta, sum)

    # Order all by metadata:
    rownames(plt.data$meta) = plt.data$meta$pr
    plt.data$mat = plt.data$mat[, plt.data$meta$pr]
    return(plt.data)
}


annotatePsbulkMetadata = function(meta, samplemeta){
    if ('ptype' %in% colnames(meta)){
        ordcol = 'ptype'
    }  else {
        ordcol = 'pr'
    }
    metaorder = meta[[ordcol]]
    advars = c('braaksc','cogdx', 'niareagansc', 'nrad','cogdxad')
    covars = c('msex','age_death','pmi', 'Apoe_e4')
    meta = merge(meta, unique(samplemeta[,
            c('projid','region','rind', advars, covars)]))
    rownames(meta) = meta[[ordcol]]
    meta = meta[metaorder,]
    return(meta)
}

