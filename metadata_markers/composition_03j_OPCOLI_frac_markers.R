#!/usr/bin/R
# ---------------------------------------------------
# Investigate differences between opc + oli subtypes:
# Updated 03/08/2021 
# TODO: Re-factor as with ast, immune.
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(qvalue)
library(lme4)
library(emmeans)

library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)
library(patchwork)

library(ComplexHeatmap)
library(circlize)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
remove.batches = TRUE
suff = '_subset_final_noMB'

subset = 'OpcOli'
ststr = gsub("/","_", subset)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
subtypes = unique(ctdf$cls)
celltypes = unique(ctdf$major.celltype)


# ------------------------------------------------------------
# Load in the full ast data for these subtypes (not strained):
# ------------------------------------------------------------
# Data directories:
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    # matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')

# Load in Immune data:
fns = list.files(path=mtxdir, pattern=paste0(rawpref,'.majorcelltype.Opc.*rda'))
fns = c(fns, list.files(path=mtxdir, pattern=paste0(rawpref,'.majorcelltype.Oli.*rda')))
amat = c(); lbs = c(); bcs = c()
for (fn in fns){
    rdafile = paste0(mtxdir, fn)
    basename = sub(".rda", "", sub(paste0(".*.majorcelltype.Mic_Immune."),"",fn))
    subtype = gsub("_"," ",sub("\\.[A-Z]*","",basename))
    region = sub(".*\\.","",basename)
    # Load `mat` from rdafile:
    load(rdafile)
    print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
    barcodes = colnames(mat)
    genes = rownames(mat)
    ngenes = nrow(mat)
    amat = cbind(amat, mat)
    bcs = c(bcs, barcodes)
    lbs = c(lbs, rep(region, length(barcodes)))
}
dim(amat)

# Normalized matrix:
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

amarg = marg[bcs,'count']
fact = amarg / median(amarg)
nmat = amat 
nmat@x <- nmat@x / rep.int(fact, diff(nmat@p))
gc()

submeta = cellmeta[cellmeta$major.celltype %in% c('Opc','Oli'),]
submeta = merge(submeta, unique(metadata[,c('projid','cogdx','niareagansc')]))
ctrl.bcs = submeta$barcode[submeta$niareagansc %in% c(3,4)]
bind = bcs %in% ctrl.bcs

# Aggregate to the level of individual x cell type:
rownames(submeta) = submeta$barcode
smeta = submeta[colnames(amat), c('projid','cell_type_high_resolution', 'barcode','region')]
ptype = paste0(smeta$projid, '_', gsub("/","_",gsub(" ","_", smeta$cell_type_high_resolution)), '_', smeta$region)
tform = make.tform(ptype, u=sort(unique(ptype)), norm=T)
apmat = amat %*% tform # raw
pmat = nmat %*% tform # norm
gc()

# -------------------------------------------------------------
# Compare the celltypes to each other at the pseudo-bulk level:
# -------------------------------------------------------------
# Pseudo-bulk metadata:
umeta = agg.rename(barcode ~ projid + cell_type_high_resolution + region, smeta, length, 'count')
umeta$ptype = with(umeta, paste0(projid, '_', gsub("/","_",gsub(" ","_", cell_type_high_resolution)), '_', region))
umeta = merge(umeta, unique(metadata[,c('projid','region','braaksc','cogdx','niareagansc','msex','age_death','pmi', 'Apoe_e4', 'nrad','cogdxad')]))
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(pmat),]
umeta = umeta[umeta$count > 5,] # Remove + will weight

# Filter by avg. expr:
ecut = 0.5
pmean = apply(pmat, 1, mean)
kept.genelist = names(pmean)[pmean > ecut]
avgidf = data.frame(symbol=names(pmean), val=pmean)

# Dummy variable of cell type:
st = topec[3]
full.regdf = c()
for (st in topec){
    print(st)
    umeta$is.st = umeta$cell_type_high_resolution == st

    gene = 'DPP10'
    est.regdf = c()
    for (gene in kept.genelist){
        x = pmat[gene, umeta$ptype]
        umeta$val = x
        # Regression - general effects of covariates on expression:
        fit = glm(val ~ is.st * nrad + Apoe_e4 + age_rescaled + msex + pmi, umeta,
                  weights=log(umeta$count), family='gaussian') # Corrects for cell ct. but not inflated
        # Alt:
        # fit = glm(val ~ cell_type_high_resolution + nrad + Apoe_e4 + age_rescaled + msex + pmi, subdf, 
        #           weights=log(subdf$count), family='gaussian') # Corrects for cell ct. but not inflated
        cfit = coefficients(summary(fit))
        df = data.frame(cfit)
        colnames(df) = c('Est','SE','t','p')
        pval = df['is.stTRUE','p']
        if (!is.na(pval)){
            if (pval < 1e-4){ cat(gene,'\t', sprintf('%0.2e',pval), '\n') }
        }
        df$var = rownames(df)
        rownames(df) = NULL
        df$symbol = gene
        est.regdf = rbind(est.regdf, df)
    }

    # Plot volcano of these assoc. with depletions:
    adf = est.regdf[est.regdf$var == 'is.stTRUE',]
    adf = merge(adf, avgidf)
    adf = adf[!is.na(adf$p),]
    adf = adf[adf$val > 1.5,]
    adf = adf[order(adf$p), ]
    adf$padj = p.adjust(adf$p, 'fdr')
    adf$log10q = -log10(adf$padj)
    pcut = 1e-3
    adf$color = 0
    adf$color[adf$padj < pcut] = 1
    adf$color[adf$padj < pcut & adf$Est > 0] = 2

    labdf = rbind(head(adf[adf$color == 1,],15),
                  head(adf[adf$color == 2,],10))

    pcols = brewer.pal(12, 'Paired')
    gplot = ggplot(adf, aes(Est, log10q, col=factor(color))) + 
        scale_color_manual(values=c('grey85',pcols[1],pcols[5])) + 
        geom_vline(xintercept=0, lwd=.25, lty='dashed') + 
        geom_point(cex=.25, alpha=1) + theme_pubr() + 
        geom_text_repel(data=labdf, aes(Est, log10q, label=symbol), max.overlaps=30, size=2, segment.size=.5) + 
        scale_y_continuous(expand=c(0,0)) + 
        labs(x='Coefficient * Average Expression') + 
        theme(legend.position = 'none')
    w = 7; h=7
    ggsave(paste0(imgpref, 'AST_subtypes_pseudobulk_ct_genes_',sub(" ", "_", st), '_volcano.png'), gplot, dpi=450, units='in', width=w, height=h)
    ggsave(paste0(imgpref, 'AST_subtypes_pseudobulk_ct_genes_',sub(" ", "_", st), '_volcano.pdf'), gplot, dpi=450, units='in', width=w, height=h)

    est.regdf$st = st
    full.regdf = rbind(full.regdf, est.regdf)
}


# -------------------------------------------------
# Plot heatmap of top types (adapt from following):
# -------------------------------------------------
adf = full.regdf[full.regdf$var == 'is.stTRUE',]
adf = adf[order(adf$p),]
adf = adf[adf$Est > 0,]
pgenes = c()
pset = c()
ntop = 7
for (st in topec){
    sdf = adf[adf$st == st,]
    pgenes = c(pgenes, head(unique(sdf$symbol),ntop))
    pset = c(pset, rep(st, ntop))
}

# White matter:
# Literature OPC + Oli genes:
# pgenes = c('C1QA', 'C3', 'CSF1R', 'CD74', 'CD33', pgenes)
# pset = c(rep('Lit.',4), pset)

kmeta = umeta[umeta$count > 5,]
plt.mat = pmat[pgenes,kmeta$ptype]
clsplit = kmeta$cell_type_high_resolution
ha = HeatmapAnnotation(CT=kmeta$cell_type_high_resolution, 
                       Region=kmeta$region,
                       col=list(Braak=colvals[['braaksc']],
                                Region=reg.cols,
                                CT=tcols[topec]))
udsplit = pset
hb = rowAnnotation(Set=pset,
                   col=list(Set=tcols[topec]))

smat = as.matrix(log(plt.mat + 1))
# smat.scaled = t(scale(t(log(plt.mat+1))))
smat.scaled = t(scale(t(smat), center=FALSE))

# png(paste0(imgpref, 'Immune_topDE_heatmap_normalized_individ.png'), res=400, units='in', width=11, height=8)
pdf(paste0(imgpref, 'Immune_topDE_heatmap_normalized_individ.pdf'), width=11, height=7)
Heatmap(smat.scaled, name='scaled\n logcounts', 
        col=viridis(100),
        use_raster=TRUE,
        top_annotation=ha, 
        column_split=clsplit, 
        show_column_names=FALSE,
        row_split=udsplit,
        right_annotation=hb
)
dev.off()








