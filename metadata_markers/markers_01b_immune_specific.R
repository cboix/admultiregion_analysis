#!/usr/bin/R
# ---------------------------------------------------
# Investigate differences between immune subtypes:
# Updated 03/08/2021 
# TODO: Refactor as with astrocytes
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

# -----------------------------------------------
# Plot immune counts across region by subtype:
# -----------------------------------------------
totdf = agg.rename(barcode ~ projid + region, cellmeta, length, 'total')
etotdf = agg.rename(barcode ~ projid + region, cellmeta[cellmeta$major.celltype == 'Mic/Immune',], length, 'exc.total')
submeta = cellmeta[cellmeta$major.celltype == 'Mic/Immune',]

metadata$cogdxad = 'CTRL'
metadata$cogdxad[metadata$cogdx %in% c(4,5)] = 'AD'
metadata$cogdxad = factor(metadata$cogdxad, levels=c('CTRL','AD'))
metadata$nrad = 'CTRL'
metadata$nrad[metadata$niareagansc %in% c(1,2)] = 'AD'
metadata$nrad = factor(metadata$nrad, levels=c('CTRL','AD'))

ctdf = agg.rename(barcode ~ projid + region + cell_type_high_resolution, submeta, length, 'count')
combdf = expand.grid(cell_type_high_resolution=unique(submeta$cell_type_high_resolution), region=unique(submeta$region), projid=unique(submeta$projid))
ctdf = merge(ctdf, combdf, all.y=TRUE)
ctdf$count[is.na(ctdf$count)] = 0
ctdf = merge(ctdf, totdf) # Removes a couple missing projid x region comb.
ctdf = merge(ctdf, etotdf)
ctdf$other = ctdf$total - ctdf$count
names(ctdf)[3] = 'celltype'
ctdf = merge(ctdf, unique(metadata[,c('projid','nft_ec','plaq_d_ec','plaq_n_ec', 'braaksc','cogdx','niareagansc','msex','age_death','pmi', 'Apoe_e4', 'nrad','cogdxad','rind', 'region')]))
ctdf$projid = factor(ctdf$projid)
topec = unique(ctdf$celltype) # For later
ctdf = ctdf[ctdf$celltype %in% topec,]
ctdf$frac = ctdf$count / ctdf$total

# --------------------
# Plot agnostic to AD:
# --------------------
aggdf = aggregate(count ~ celltype + region, ctdf, sum)
totdf = aggregate(count ~ celltype, ctdf, sum)
totdf$region = 'All'
aggdf = rbind(aggdf, totdf[,colnames(aggdf)])
aggdf$region = factor(aggdf$region, levels= c('All','AG','MT','PFC','EC','HC','TH'))

gplot = ggplot(aggdf, aes(region, count, fill=celltype)) + 
    geom_bar(position='fill', stat='identity') + 
    scale_fill_manual(values=tcols) + 
    theme_pubr() + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    labs(x='Region', y='Percentage')
ggsave(paste0(imgpref, 'fractions_immune_hct_major_cols_', lblset, '.pdf'), gplot, dpi=450, units='in', width=5, height=6)

# ---------------
# Plot w.r.t. AD:
# ---------------
aggdf = aggregate(count ~ celltype + region + nrad, ctdf, sum)
totdf = aggregate(count ~ celltype + nrad, ctdf, sum)
totdf$region = 'All'
aggdf = rbind(aggdf, totdf[,colnames(aggdf)])
aggdf$region = factor(aggdf$region, levels= c('All','AG','MT','PFC','EC','HC','TH'))

gplot = ggplot(aggdf, aes(nrad, count, fill=celltype)) + 
    facet_wrap(~region, nrow=1) + 
    geom_bar(position='fill', stat='identity') + 
    scale_fill_manual(values=tcols) + 
    theme_pubr() + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    labs(x='Region', y='Percentage')
ggsave(paste0(imgpref, 'fractions_immune_hct_major_cols_byad_', lblset, '.pdf'), gplot, dpi=450, units='in', width=7, height=6)

# Plot by individual / region:
# ggplot(ctdf, aes(celltype, frac, alpha=nrad, fill=celltype)) + 
gplot = ggplot(ctdf, aes(region, frac, alpha=nrad, fill=celltype)) + 
    facet_wrap(~celltype, nrow=1) + 
    geom_boxplot(outlier.shape=NA) + 
    scale_fill_manual(values=tcols) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
    theme_pubr() + 
    # scale_y_continuous(expand=c(0,0),labels=scales::percent) + 
    stat_compare_means(hide.ns=TRUE, label='p.format') + 
    labs(x='Region', y='Percentage')
ggsave(paste0(imgpref, 'boxplots_immune_hct_major_cols_byad_', lblset, '.pdf'), gplot, dpi=450, units='in', width=8, height=4)

# ---------------------------
ctdf = merge(ctdf, pqdf, all.x=TRUE)
ctdf$lnft = log(ctdf$nft + 1)

pathval = 'nrad'
lmdf = ctdf[!is.na(ctdf[[pathval]]),]
asform = function(x){as.formula(paste(x, collapse=" "))}
gform = asform(c('cbind(count, other) ~ ',pathval, '* celltype + celltype *region * msex'))
fit = glm(gform, lmdf, family='quasibinomial')

emform = asform(c('revpairwise ~', pathval, '|celltype'))
emm1 <- emmeans(fit, specs=emform)
cdf = as.data.frame(rbind(summary(emm1$contrasts, 
                                  infer = TRUE, type = 'response')))
cdf = cdf[order(cdf$odds.ratio),]
cdf$p.adj = p.adjust(cdf$p.value, 'fdr')
cdf$celltype = factor(cdf$celltype, levels=rev(cdf$celltype))

p.cut = 0.05
cdf$col = 1 * (cdf$p.adj < p.cut) + 2 * (cdf$odds.ratio > 1 )
labdf = cdf[cdf$p.adj < p.cut,]
gplot = ggplot(cdf, aes(x=odds.ratio, y=celltype, col=factor(col))) +
    geom_point() +
    geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend =celltype)) +
    geom_text(data=cdf, aes(x=max(labdf$odds.ratio) * 1.1, y=celltype, label=sprintf("OR=%0.2f p=%0.1e",odds.ratio,p.adj)),col='black') +
    # scale_x_log10() +
    scale_color_manual(values =c('0' ='grey75','1'='royalblue','2'='grey75','3'='indianred')) + 
    geom_vline(xintercept=1, lty='dashed') + 
    theme_pubr() + theme(legend.position = 'none') +
    labs(x = cdf$contrast[1], y='Immune Subtype')

h = nrow(labdf) / 30 * 8 + 1
ggsave(paste0(imgpref, 'oddratios_immune_by_', pathval,'.pdf'), gplot, dpi=450, units='in', width=3, height=h)

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
fns = list.files(path=mtxdir, pattern=paste0(rawpref,'.majorcelltype.Mic.*rda'))
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

submeta = cellmeta[cellmeta$major.celltype == 'Mic/Immune',]
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



# --------------------------------------------------
# Plot some of these genes on the Immune -only UMAP:
# --------------------------------------------------
ind = 1:nrow(submeta)
typelvls = unique(submeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
tsp.type.cols = sapply(type.cols, tsp.col)
celltype.loc = aggregate(cbind(U1, U2) ~ cell_type_high_resolution, submeta, mean)
cex = 0.025

xlim = c(min(submeta$U1) + 1.25, min(submeta$U1) + 9.25)
ylim = c(max(submeta$U2)-12.5, max(submeta$U2) - 4)
tsp.tcols = sapply(tcols, tsp.col)
tsp.rcols = sapply(reg.cols, tsp.col)

png(paste0(imgpref, 'astumap_final_highres_cols_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=3, height=3)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(submeta$U1[ind], submeta$U2[ind], col=tsp.tcols[submeta$cell_type_high_resolution[ind]], ylim=ylim, xlim=xlim, pch=19, cex=cex, axes=F)
dev.off()

png(paste0(imgpref, 'astumap_region_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=3, height=3)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=rep(sp, 4))
plot(submeta$U1[ind], submeta$U2[ind], col=tsp.rcols[submeta$region[ind]], ylim=ylim, xlim=xlim, pch=19, cex=cex, axes=F)
dev.off()


palette = viridis(100)
col_fun = function(x, pal=palette){
                bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
                palette[bin]  }

pgenes = c('P2RY12', 'SIGLEC1', 'CENPP','DUSP1','FYN','TPT1',
           'CD81', 'CD74', 'APOE', 'HLA-DRB1'
           )

for (gene in rev(pgenes)){
    print(gene)
    x = log(nmat[gene,submeta$barcode] + 1)
    png(paste0(imgpref, 'astumap_gene_', gene,'_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=3, height=3)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    par(mar=rep(sp, 4))
    plot(submeta$U1[ind], submeta$U2[ind], col=col_fun(x[ind]), ylim=ylim, xlim=xlim, pch=19, cex=cex, axes=F)
    mtext(gene, side=1, line=-1, font=2, cex=1.75)
    dev.off()
}








