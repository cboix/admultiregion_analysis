#!/usr/bin/R
# --------------------------------------------
# Plot the expression patterns for GWAS genes:
# Updated 02/21/2022
# --------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'gwas/')
imgpref = paste0(plotdir, 'gwas_')
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)


source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])


# Load pseudobulk data for all runsets:
# -------------------------------------
reg.ps.rda = paste0(srdir, 'pseudobulk_data_all_regionaverages.rda')
if (!file.exists(reg.ps.rda)){
    umeta = NULL
    pmat = NULL
    runlist = c('Ast', 'Mic_Immune', 'Vasc_Epithelia',
                'Oli','Inh','Opc', 'HCneurons', 'ECneurons',
                'THneurons', 'CTXneurons')

    for (runset in runlist){
        print(runset)
        psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
        load(psdata.rda)
        # Metadata for the pseudobulk matrix:
        ps.data$meta$runset = runset
        umeta = rbind(umeta, ps.data$meta)
        pmat = cbind(pmat, ps.data$mat)
    }


    # Merge the individuals + neurons, keep regional differences:
    # -----------------------------------------------------------
    umeta$runset[umeta$runset %in% names(exc.sets)] = 'Exc'
    umeta$regset = with(umeta, paste0(runset, '-', region))
    regsets = unique(umeta$regset)

    # Average, weight by number of cells
    tform = make.tform(umeta$regset, u=regsets)
    tform = sweep(tform, 1, umeta$ncell, '*')
    totncell = apply(tform, 2, sum)
    avg.mat = pmat %*% tform
    avg.mat = sweep(avg.mat, 2, totncell, '/')
    save(avg.mat, file=reg.ps.rda)
} else {
    load(reg.ps.rda)
}


# Check: scale average matrix to number of counts for each cell type?
# NOTE: No, this is not needed, current normalization is appropriate.
# -------------------------------------------------------------------
# Load ncounts margin data:
oldpref = 'all_brain_regions_filt_preprocessed_scanpy'
margfile = paste0(sdbdir, 'matrices/', oldpref, '_fullmatrix_margin.tsv.gz')
bcfile = paste0(sdbdir, oldpref, '_norm.barcodes.tsv.gz')
marg = as.integer(scan(gzfile(margfile), 'c', quiet=T))
bcs = scan(gzfile(bcfile), 'c', quiet=T)
names(marg) = bcs

# Get ncounts per cell type
cellmeta$ncounts = marg[cellmeta$barcode]
margdf = aggregate(ncounts ~ major.celltype, cellmeta, median)
rownames(margdf) = sub("/", "_", margdf$major.celltype)

sfactor = margdf[sub("-.*", "", colnames(avg.mat)), 'ncounts']


# Plot all of the kept genes (protein-coding):
# --------------------------------------------
kept.genes = gwgenes[gwgenes %in% rownames(avg.mat)]
cmat = as.matrix(avg.mat[kept.genes,])
cmat = log1p(cmat)

# Norm to max:
cmat = sweep(cmat, 1, apply(cmat, 1, max) + 1e-4, '/')

# Reorder:
cmat = reord(cmat, measure='euclidean', method='ward.D')
cts = sub("-.*","", colnames(cmat))
uqcts = c('Exc','Inh','Ast','Mic_Immune','Opc','Oli','Vasc_Epithelia')
tform = make.tform(cts, u=uqcts, norm=T)
amat = cmat %*% tform
amat = sweep(amat, 1, apply(amat, 1, max) + 1e-4, '/')
amat = t(diag.mat2(t(amat), cutoff=.75, ratio=0.8)[[1]])
cmat = cmat[rev(rownames(amat)),]

reg.ordered = c('AG','MT','PFC','HC','EC','TH')
cdf = expand.grid(reg=reg.ordered, ct=uqcts)
pastedash = function(x,y){paste0(x,'-', y)}
colord = outer(uqcts, reg.ordered, FUN="pastedash")
colord = c(t(colord))
cmat = cmat[,colord]

ux = 2
plt = Heatmap(cmat,
              # col=rev(colrb),
              col=viridis(50),
              name='expr',
              use_raster=TRUE,
              column_split=sub("-.*","", colnames(cmat)),
              cluster_columns=FALSE,
              cluster_column_slices=FALSE,
              cluster_rows=FALSE,
              width = ncol(cmat)*unit(ux / 2, "mm"), 
              height = nrow(cmat)*unit(ux, "mm"),
              border_gp = gpar(col="black", lty = 1, lwd=.5),
)

h = 3 + 1 / 15 * nrow(cmat)
w = 5 + 1 / 15 * ncol(cmat)
pltprefix = paste0(imgpref, 'gwas_avgexpr_allgwgenes_heatmap')
saveHeatmap(plt, pltprefix, w=w, h=h)



# Save top expressed:
# -------------------
topdf = data.frame(gene=rownames(amat),
                   ct=colnames(amat)[apply(amat, 1, which.max)])
write.table(topdf, 'Annotation/ADGWAS_topct_multiregion_121721.tsv', quote=F, row.names=F, sep="\t")


# Save matrices and expression attributes:
# ----------------------------------------
gw.ps.rda = paste0(srdir, 'pseudobulk_data_gwasgenes_regionaverages.rda')
norm.exprmat = cmat
full.exprmat = as.matrix(avg.mat[rownames(norm.exprmat),])
save(norm.exprmat, full.exprmat, file=gw.ps.rda)


