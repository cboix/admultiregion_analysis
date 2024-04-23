#!/usr/bin/R
# -----------------------------------------------
# Load pseudobulk data for differential analysis:
# Updated: 10/28/23
# -----------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(Matrix)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype subtype")
} else if (args != 'local'){
    celltype = args[1]
    subtype = args[2] 
    region = args[3]
}


source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


# Load major celltype-level dataset:
# ----------------------------------
ct = gsub("_","/", celltype)
psdata.rda = paste0(srdir, 'pseudobulk_data_', celltype, '.counts.rda')
if (!file.exists(psdata.rda)){
    if (celltype %in% c('Mic_Immune', 'Opc', 'Oli')){
        ctcol = 'minor.celltype'
    } else {
        ctcol = 'cell_type_high_resolution'
    }
    subtypes = unique(cellmeta[cellmeta$major.celltype == ct, ctcol])
    ps.data = load_pseudobulk_dataset(
        ct, subtypes, reg.nomb, normalize=FALSE)
    save(ps.data, file=psdata.rda)
} else { 
    load(psdata.rda) 
}


# Specific subsets - subtype and region-specific:
# -----------------------------------------------
if (celltype != subtype){
    if (celltype == 'Mic_Immune'){
        if (subtype == 'CAM'){
            set = c('CAMs')
        } else if (subtype == 'T'){
            set = c('T cells')
        } else {
            subtypes = unique(ps.data$meta$cell_type_high_resolution)
            set = subtypes[grep("^Mic", subtypes)]
        }
        ind = ps.data$meta$cell_type_high_resolution %in% set
    } else {
        submap = unique(ps.data$meta$cell_type_high_resolution)
        names(submap) = gsub(" ", "_", submap)
        ind = ps.data$meta$cell_type_high_resolution == submap[subtype]
    }
    ps.data$meta = ps.data$meta[ind,]
    ps.data = harmonizePsbulkData(ps.data)
}


if (region != 'allregions'){
    ind = ps.data$meta$region == region
    ps.data$meta = ps.data$meta[ind,]
    ps.data = harmonizePsbulkData(ps.data)
}


# Load and add the appropriate metadata:
# --------------------------------------
extsi_tsv = 'Annotation/extended_si_ace_vars.tsv'
ext.si.vars = read.delim(extsi_tsv, header=F)[,1]

advars = c('cogdx', 'nrad', 'cogdxad', 'braaksc.early','braaksc.ad')
covars = c('age_death','msex','pmi', 'Apoe_e4') 
ext.covars = c('soc_net_bl', 'social_isolation_avg', 
    'social_isolation_lv', ext.si.vars)

ext.meta_csv = 'Annotation/ROSMAP_All_Pts_Extended_Basic_Long_652_03-23-2022.csv'
ext.metadata = read.delim(ext.meta_csv, header=T, sep=",")

# Merge all these variables in:
ps.data$meta = merge(ps.data$meta, 
    metadata[,c('projid','region', 'rind', covars, advars)])
ps.data$meta = merge(ps.data$meta, ext.metadata[,c('projid',ext.covars)])
ps.data$meta = merge(ps.data$meta, pqdf, all.x=TRUE)
ps.data$meta = ps.data$meta[order(ps.data$meta$projid),]

# Remove samples with very few cells:
ncell.cutoff = 10
ind = ps.data$meta$ncell >= ncell.cutoff
ps.data$meta = ps.data$meta[ind,]

# Harmonize metadata with mat
ps.data = harmonizePsbulkData(ps.data)


