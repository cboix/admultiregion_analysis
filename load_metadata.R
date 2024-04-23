#!/usr/bin/R
# ------------------------------------------
# Load metadata for the multiregion project:
# ------------------------------------------
library(cbrbase)
library(RColorBrewer)
library(tidyr)
set_proj('DEVTRAJ')

# Arguments:
args = commandArgs()
if (length(args) > 0){
    project = args[1]
} else { 
    project = 'RNA_Regions'
}

# ---------------
# Metadata files:
# ---------------
mfiles = list.files(path='Annotation', pattern='^metadata_[A-Z]*.tsv')
regions = sub(".tsv", "", sub("metadata_", "",mfiles))
reg.nomb = regions[regions != 'MB']
NREGIONS = length(regions)
reg.order = c('All','AG','MT','PFC','EC','HC','TH')

# Load in metadata:
metadatafile = 'Annotation/regions_metadata_list.Rda'
if (!file.exists(metadatafile)){ 
    metalist = list()
    for (i in 1:NREGIONS){
        tab = read.delim(paste0('Annotation/', mfiles[i]), header=T, sep="\t")
        tab$region = regions[i]
        metalist[[regions[i]]] = tab
        if (i == 1){
            kcols = colnames(tab)
        } else {
            kcols = intersect(kcols, colnames(tab))
        }
    }
    kcols = kcols[kcols != 'X']
    redmetalist = list()
    for (i in 1:NREGIONS){
        redmetalist[[regions[i]]] = metalist[[regions[i]]][,kcols]
    }
    metadata = do.call(rbind, redmetalist)
    metadata$rind = rownames(metadata)
    metadata$cogdxad = 'CTRL'
    metadata$cogdxad[metadata$cogdx %in% c(4,5)] = 'AD'
    metadata$cogdxad = factor(metadata$cogdxad, levels=c('CTRL','AD'))
    metadata$nrad = 'CTRL'
    metadata$nrad[metadata$niareagansc %in% c(1,2)] = 'AD'
    metadata$nrad = factor(metadata$nrad, levels=c('CTRL','AD'))
    # Late and early braak stages:
    metadata$braaksc.ad = 'CTRL'
    metadata$braaksc.ad[metadata$braaksc >= 5] = 'AD'
    metadata$braaksc.ad= factor(metadata$braaksc.ad, levels=c('CTRL','AD'))
    metadata$braaksc.early = 'CTRL'
    metadata$braaksc.early[metadata$braaksc >= 3] = 'AD'
    metadata$braaksc.early = factor(metadata$braaksc.early, levels=c('CTRL','AD'))
    save(metadata, metalist, file=metadatafile)
} else {
    load(metadatafile)
}

kept.individuals = scan('Annotation/multiRegion_individuals.txt', 'c')

# -------------------
# Genome Annotations:
# -------------------
gencode_version = 'v28lift37.annotation'
gencoderda = paste0('Annotation/Processed.gene.', gencode_version, '.Rda')
if (!file.exists(gencoderda)){
    geneanno = read.delim(paste0('Annotation/Gene.', gencode_version, '.bed'), header=F, stringsAsFactors=F)
    anno = read.delim(paste0('Annotation/Tss.', gencode_version, '.bed'), header=F, stringsAsFactors=F)
    names(geneanno) <- c('chr','start', 'end', 'strand', 'ENSG','type','symbol')
    names(anno) <- c('chr','tss','ENSG','type','symbol')
    anno$length_35 = geneanno$end - geneanno$start

    # -----------------------
    # Get transcript lengths:
    # -----------------------
    alt.txlenfile = paste('Annotation/ExonMerge.Homo_sapiens.GRCh37.87.gene.totals.bed')
    alt.txdf = read.table(alt.txlenfile, header=F)
    names(alt.txdf) = c('symbol','tx.length','version')
    txlenfile = paste0('gencode.',gencode_version,'.txlen.tsv')
    if (file.exists(txlenfile)){
        print("Loading transcript lengths")
        txdf = read.delim(txlenfile, sep="\t", header=T)
    } else {
        # Import GTF:
        require(rtracklayer)
        print("Reading in GTF to get transcript lengths")
        GTFfile = paste0('Annotation/gencode.',gencode_version,'.gtf')
        GTF <- import.gff(GTFfile, format="gtf", genome="GRCh37", feature.type="exon")
        # Process gtf to get transcript length:
        grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
        reducedGTF <- unlist(grl, use.names=T)
        elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
        elementMetadata(reducedGTF)$widths <- width(reducedGTF)
        calc_length <- function(x) { sum(elementMetadata(x)$widths) }
        txlen <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
        txdf = data.frame(ENSG=colnames(txlen), length=as.numeric(txlen))
        write.table(txdf, txlenfile, sep="\t", quote=F, col.names=T, row.names=F)
    }
    # Update annotation:
    anno = merge(anno, txdf)
    save(anno, file=gencoderda)
} else {
    load(gencoderda)
}


# --------------------
# Colors for plotting:
# --------------------
source(paste0(sbindir, 'set_colors.R'))


# -------------------------------
# Read in the full cell metadata:
# -------------------------------
datadir = 'multiRegion/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy_norm'
lblset = 'leiden_r15_n100'
# Only for sub-clustering:
load(paste0(datadir, prefix, '.', lblset, '.ext.lbl.Rda'))

hlvls = as.character(sort(unique(celldf$lbl)))
lbl.cols = rep(snap.cols,4)[1:length(hlvls)]
names(lbl.cols) = as.character(hlvls)
if (length(grep('hdb', lblset)) > 0){ ctype = 'HDBSCAN' } else { ctype = 'Leiden' }


# --------------------------------------
# Load in the final metadata (cellmeta):
# --------------------------------------
load(file=paste0(datadir, prefix, '.final_noMB.cell_labels.Rda'))

# Colors for full:
typelvls = unique(cellmeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
type.cols = c(type.cols, major.col['Inh'], major.col['Exc'])
tsp.type.cols = sapply(type.cols, tsp.col)

sampmeta.rds = 'Annotation/multiregion_sample_metadata.Rds'
if (!file.exists(sampmeta.rds)){
    uqrind = unique(cellmeta$rind)
    saveRDS(metadata[metadata$rind %in% uqrind,], file=sampmeta.rds)
}

cellmeta.rds = 'Annotation/multiregion_cell_metadata.Rds'
if (!file.exists(cellmeta.rds)){
    sel.cols = c('barcode','rind','region','projid','U1','U2',
                 'major.celltype','minor.celltype','neuronal.layer','inh.subtype',
                 'neuronal.exttype','full.exttype','cell_type_high_resolution')
    saveRDS(cellmeta[,sel.cols], file=cellmeta.rds)
}

full.projids = sort(unique(cellmeta$projid))
projid.cols = snap.cols[1:length(full.projids)]
names(projid.cols) = full.projids


# All multiregion color sets, for visualization:
# ----------------------------------------------
col.rds = 'Annotation/multiregion_color_sets.Rds'
if (!file.exists(col.rds)){
    mrc = list()
    mrc[['major']] = major.col
    mrc[['major.tsp']] = tsp.major.col
    mrc[['covariates']] = colvals
    mrc[['individual']] = ind.cols
    mrc[['lbl']] = lbl.cols
    mrc[['region']] = reg.cols[reg.nomb]
    mrc[['region_label']] = reg.long[reg.nomb]
    mrc[['celltype']] = tcols
    mrc[['celltype.tsp']] = tsp.tcols
    mrc[['Paired']] = col.paired
    mrc[['Reds']] = colr
    mrc[['Blues']] = colb
    mrc[['RdBu']] = colrb
    mrc[['extra']] = snap.cols
    mrc[['projid']] = projid.cols
    saveRDS(mrc, file=col.rds)
}


# Get the pathology mapped to each region:
# ----------------------------------------
pathrda = 'Annotation/multiregion_pathology_df.Rda'
if (!file.exists(pathrda)){
    pqdf = NULL
    for (path in c('nft','plaq_d','plaq_n')){
        regmap = c('AG','HC','PFC','MT','EC')
        names(regmap) = c('ag','hip','mf','mt','ec')
        vars = colnames(metadata)[grep(path, colnames(metadata))]
        vars = vars[vars != path]
        submeta = unique(metadata[,c('projid','region', vars, 'rind')])
        slong = gather(submeta, path, value, -projid, -region, -rind)
        slong$path.region = regmap[sub(".*_","", slong$path)]
        slong = slong[slong$region == slong$path.region,]
        rownames(slong) = slong$rind
        if (is.null(pqdf)){
            pqdf = slong[,c('rind','value','region')]
            names(pqdf)[2] = path
        } else { 
            sub.pqdf = slong[,c('rind','value','region')]
            names(sub.pqdf)[2] = path
            pqdf = merge(pqdf, sub.pqdf)
        }
    }
    save(pqdf, file=pathrda)
} else { load(pathrda) }


# Excitatory neuron to region annotations:
# ----------------------------------------
source(paste0(sbindir, 'set_neuron_subsets.R'))
