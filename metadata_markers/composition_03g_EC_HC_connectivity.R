#!/usr/bin/R
# ----------------------------------------------------------------------------
# Investigate the co-variation of selective vulnerability in EC/HC and others:
# TODO: Neurons alone and neurons + glia
# Updated 05/10/2021 
# -----------------------------------------------------
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
imgpref = paste0(plotdir, 'connect_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# ----------------------------------------
# Plot neuron counts in entorhinal cortex:
# ----------------------------------------
totdf = agg.rename(barcode ~ projid + region, cellmeta, length, 'total')
etotdf = agg.rename(barcode ~ projid + region, cellmeta[cellmeta$major.celltype == 'Exc',], length, 'exc.total')
submeta = cellmeta[cellmeta$major.celltype == 'Exc' & cellmeta$region %in% c('EC', 'HC'),]

ctdf = agg.rename(barcode ~ projid + region + cell_type_high_resolution, submeta, length, 'count')
combdf = expand.grid(cell_type_high_resolution=unique(submeta$cell_type_high_resolution), region=unique(submeta$region), projid=unique(submeta$projid))
ctdf = merge(ctdf, combdf, all.y=TRUE)
ctdf$count[is.na(ctdf$count)] = 0
ctdf = merge(ctdf, totdf) # Removes a couple missing projid x region comb.
ctdf = merge(ctdf, etotdf)
ctdf$other = ctdf$total - ctdf$count
names(ctdf)[3] = 'celltype'
ctdf = merge(ctdf, unique(metadata[,c('projid','nft_ec','plaq_d_ec','plaq_n_ec', 'braaksc','cogdx','niareagansc','msex','age_death','pmi', 'Apoe_e4', 'nrad','cogdxad')]))
gdf = aggregate(gpath ~ projid, metadata, mean)
gdf = gdf[gdf$projid %in% ctdf$projid,] 
gdf$gfact = 0
gdf$gfact[order(gdf$gpath)] = 1:48
ctdf = merge(ctdf, gdf)
ctdf$projid = factor(ctdf$projid)
ctdf$frac = ctdf$count / ctdf$total

# Select neuronal subtypes in EC only:
rdf = agg.rename(barcode ~ region + cell_type_high_resolution, cellmeta[cellmeta$major.celltype == 'Exc',], length, 'count')
rdf = spread(rdf, region, count, fill=0)
rmat = as.matrix(rdf[,-1])
rownames(rmat) = rdf[,1]
topec = apply(rmat,1, max) == rmat[,'EC']
topec = names(topec)[topec]
topec = topec[topec != 'Exc SV2C LINC02137']
tophc = apply(rmat,1, max) == rmat[,'HC']
tophc = names(tophc)[tophc]
# ctdf = ctdf[ctdf$celltype %in% topec * ,]
df1 = ctdf[ctdf$celltype %in% topec & ctdf$region == 'EC',]
df2 = ctdf[ctdf$celltype %in% tophc & ctdf$region == 'HC',]
ctdf = rbind(df1, df2)

ggplot(unique(ctdf), aes(gfact, frac, color=celltype)) + 
    facet_wrap(~region) + 
    geom_point() + 
    geom_smooth(method='lm') + 
    theme_pubr()

# By AD + sex, is there a difference: 
# - No overall sex-difference
# - AD depletion occurs in both M/F in the specified cell types 
ggplot(unique(ctdf), aes(celltype, frac, color=factor(msex), fill=nrad)) + 
    facet_wrap(~region) + 
    geom_boxplot() + 
    theme_pubr() + coord_flip()

# -----------------------------------------------------------
# Load in the full EC data for these subtypes (not strained):
# -----------------------------------------------------------
# Data directories:
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    # matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')

# Load in data matrices:
celltype = 'Exc'
amat = c(); lbs = c(); bcs = c()
for (subtype in c(topec, tophc)){
    if (subtype %in% topec){region='EC'} else {region ='HC'}
    ststr = gsub("/","_",gsub(" ","_", subtype))
    matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                     celltype,'.',ststr,'.',region)
    rdafile = paste0(matpref, '.rda')  # In Matrix format
    # Load `mat` from rdafile:
    load(rdafile)
    print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
    barcodes = colnames(mat)
    genes = rownames(mat)
    ngenes = nrow(mat)
    amat = cbind(amat, mat)
    bcs = c(bcs, barcodes)
    lbs = c(lbs, rep(subtype, length(barcodes)))
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

submeta = cellmeta[cellmeta$major.celltype == 'Exc' & cellmeta$region %in% c('HC','EC'),]
submeta = merge(submeta, unique(metadata[,c('projid','cogdx','niareagansc')]))
ctrl.bcs = submeta$barcode[submeta$niareagansc %in% c(3,4)]
bind = bcs %in% ctrl.bcs

# Aggregate to the level of individual x cell type:
rownames(submeta) = submeta$barcode
smeta = submeta[colnames(amat), c('projid','cell_type_high_resolution', 'barcode')]
cts = gsub(" ","_", smeta$cell_type_high_resolution)
# ptype = paste0(smeta$projid, '_', gsub("/","_",gsub(" ","_", smeta$cell_type_high_resolution)))
# tform = make.tform(ptype, u=sort(unique(ptype)), norm=T)
# apmat = amat[c('CDH20','HS6ST3','PRNP','APP'),] %*% tform #raw
# npmat = nmat[c('CDH20','HS6ST3','PRNP','APP'),] %*% tform #raw
tform = make.tform(cts, u=sort(unique(cts)), norm=T)
acmat = amat[c('CDH20','HS6ST3','PRNP','APP'),] %*% tform #raw
ncmat = nmat[c('CDH20','HS6ST3','PRNP','APP'),] %*% tform #raw

# Stellate: 
# Grid:
# RELN: Layer 2
# https://www.hindawi.com/journals/np/2010/108190/
# Layer II neurons show a variety of molecular alterations in AD, including reductions in muscarinic acetylcholine receptor 1, GABAA receptor delta, and ionotropic glutamate receptor NMDA 1
markers = c('CALB1','CALB2','RELN', 'CHRM1','GABRD', 'GRIN1', 'BDNF', 
            'TOX3','AGBL1','POSTN','DLC1','NTS', 'GPC5','GPC6', 'NEFH',
            'SLC18A3','ETV1','BCL11B','WDR16','FABP5','IGFBP6','KCTD16', 'NOV',
            'FSTL1','NEF3','CUTL2','COL5A1','COL5A2','SEMA3C','KCNG1', 
            'GRP','NTS', 'JUP','NXPH4','COBLL1','THSD7B', 'DCC','CCK',
            'KITL','FEZF2','WFS1','IL1RAPL2','MRG1','COL20A1','VIPR1'
)
adrs = c(scan('Annotation/adrenoreceptors.tsv', 'c'), 'SLC6A2')
markers = c(markers, adrs)

markers = markers[markers %in% rownames(nmat)]
ncmat = nmat[markers,] %*% tform #raw

ncmat = t(scale(t(ncmat),center=FALSE))
ncmat[is.na(ncmat)] = 0
topec2 = gsub("[,)(]","", gsub(" ", "_", c(topec)))
colnames(ncmat) = gsub("[,)(]","", gsub(" ", "_", colnames(ncmat)))

png(paste0(imgpref, 'heatmap_markers_ECHC_adrs.png'), units='in', res=450, width=8, height=12)
Heatmap(as.matrix(ncmat), 
        column_split=ifelse(colnames(ncmat) %in% topec2, 'Entorhinal Ctx','Hippocampus'),
        row_split=ifelse(rownames(ncmat) %in% adrs, 'Adrenoceptors','Other'),
        col=viridis(100))
dev.off()

gcout = gc()

# -----------------------------------------
# Plot a large number of neuronal receptors
# -----------------------------------------
nrdf = read.delim('Annotation/neuronal_receptors.tsv', sep=",")
npdf = read.delim('Annotation/neuropeptide_ligands.tsv', sep=",")

markers = nrdf$Approved.symbol
sets = nrdf$Group
mind = markers %in% rownames(nmat)
markers = markers[mind]
sets = sets[mind]
ncmat = nmat[markers,] %*% tform #raw

ncmat = t(scale(t(ncmat),center=FALSE))
ncmat[is.na(ncmat)] = 0
topec2 = gsub("[,)(]","", gsub(" ", "_", c(topec)))
colnames(ncmat) = gsub("[,)(]","", gsub(" ", "_", colnames(ncmat)))

png(paste0(imgpref, 'heatmap_markers_ECHC_receptors.png'), units='in', res=450, width=25, height=4)
Heatmap(t(as.matrix(ncmat)), 
        row_split=ifelse(colnames(ncmat) %in% topec2, 'Entorhinal Ctx','Hippocampus'),
        column_split=sets,
        row_names_gp=gpar(fontsize=11),
        column_names_gp=gpar(fontsize=10),
        col=viridis(100))
dev.off()

gcout = gc()



# -------------------------------------
# Plot these markers on the EC neurons:
# -------------------------------------
tab = read.delim('multiRegion/metadata_test_mrad_Exc_EC_combat_filthvg.tsv', header=T)
tab = tab[,c('barcode','projid','celltype','U1','U2')]
tab = merge(tab, unique(metadata[,c('projid','nrad')]))

# Shuffle points for plotting:
ind = 1:nrow(tab)
ad.ind = which(tab$nrad == 'AD')
ctrl.ind = which(tab$nrad == 'CTRL')
full.ind = sample(ind, length(ind), replace=FALSE)
ad.ind = sample(ad.ind, length(ad.ind), replace=FALSE)
ctrl.ind = sample(ctrl.ind, length(ctrl.ind), replace=FALSE)
xr = diff(range(tab$U1))
yr = diff(range(tab$U2))

# Plot the umap, with appropriate proportions:
w = 3
h = w * yr / xr
cex = 0.025
sp = 0.1
png(paste0(imgpref, 'EC_vuln_neurons_umap_alone.png'), units='in', res=450, width=w, height=h)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(tab$U1[full.ind], tab$U2[full.ind], 
     col=tcols[tab$celltype[full.ind]], pch=19, cex=cex, axes=F)
dev.off()

palette = c('grey90',viridis(100))
# palette = c(NA,viridis(100))
col_fun = function(x, pal=palette){
    bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
    palette[bin] 
}

# Stellate: 
# Grid:
# RELN: Layer 2
# https://www.hindawi.com/journals/np/2010/108190/
# Layer II neurons show a variety of molecular alterations in AD, including reductions in muscarinic acetylcholine receptor 1, GABAA receptor delta, and ionotropic glutamate receptor NMDA 1

# Set of 20 markers:
plt.genes = c('TOX3','AGBL1','POSTN','DLC1','NTS',
              'CALB1','CALB2','RELN','GPC5','GPC6',
              'BCL11B','ETV1','NXPH4','COL5A1','COL5A2',
              adrs)
scale = .5
w2 = 5 * w * scale
h2 = 5 * h * scale
cex=0.01
png(paste0(imgpref, 'EC_vuln_neurons_umap_marker_genes.png'), units='in', res=450, width=w2, height=h2)
par(xaxs='i',yaxs='i')
layout(matrix(1:25, nrow=5,5, byrow=TRUE))
par(mar=rep(sp,4))
for (gene in plt.genes){
    x = log(nmat[gene, tab$barcode] + 1)
    xind = order(x)
    # plot(tab$U1[full.ind], tab$U2[full.ind], 
    #      col=col_fun(x[full.ind]), pch=19, cex=cex, axes=F)
    plot(tab$U1[xind], tab$U2[xind], 
         col=col_fun(x[xind]), pch=19, cex=cex, axes=F)
    mtext(side=3, gene, line=-1, cex=.75)
}
dev.off()


# ----------------------------------
# Analyze the rates of co-depletion:
# ----------------------------------
pathval = 'nrad'
# TODO: CLR transform?
cdf = ctdf[,c('projid','celltype','frac','nft_ec', 'nrad', 'cogdxad')]
cdf$celltype = gsub("[,)(]","", gsub(" ", "_", cdf$celltype))
ctwide = spread(cdf, celltype,frac)

c1 = 'Exc_AGBL1_GPC5'
c2 = 'Exc_RELN_COL5A2'
c2 = 'Exc_RELN_GPC5'
c2 = 'CA1_pyramidal_cells'
c2 = 'CA2_CA3_pyramidal_cells'

# c2 = 'Exc_TOX3_TTC6'
ggplot(ctwide, aes_string(x=c1, y=c2, color='nrad')) +
    # g0a = ggplot(ctdf, aes(x=celltype, y=count / exc.total, color=nrad)) +
    geom_point()+ 
    # geom_smooth(method='lm', color='black') + 
    geom_smooth(method='lm') + 
    scale_color_manual(values=colvals[['nrad']]) + 
    theme_pubr()

# Heatmap: which neurons predictive of which others (+Region) 
# 1. Overall
# 2. When split into AD / non-AD
# 3. When split by Cognition

# Add in interaction by expression of HS6ST3 and CDH20
# pids = sub("_.*","",colnames(npmat))
# cts = sub("[0-9]*_","",colnames(npmat))
# edf = data.frame(cell_type_high_resolution=cts, val=npmat['HS6ST3',])
# edf = data.frame(cell_type_high_resolution=cts, val=npmat['CDH20',])
# edf = aggregate(val ~ cell_type_high_resolution, edf, mean)

topec2 = gsub("[,)(]","", gsub(" ", "_", c(topec, tophc)))
topec2 = gsub("[,)(]","", gsub(" ", "_", c(topec)))
topec2 = topec2[topec2 != 'Exc_SOX11_NCKAP5']
NC = length(topec2)
cmat = matrix(NA, NC, NC, dimnames=list(topec2, topec2))
pmat = matrix(NA, NC, NC, dimnames=list(topec2, topec2))
advar = 'none'
advar = 'nrad'
for (c1 in topec2){
    for (c2 in topec2){
        if (c1 != c2){
            df = ctwide[,c(c1,c2,'nrad', 'cogdxad')]
            if (advar == 'none'){
                val = c1
                fit = lm(as.formula(paste0(c2, '~', c1)), df)
            } else {
                val = paste0(c1,":",advar,'AD')
                fit = lm(as.formula(paste0(c2, '~', c1, '*',advar)), df)
                # val = c1
                # fit = lm(as.formula(paste0(c2, '~', c1)), df[df[[advar]] == 'AD',])
            }
            cfit = data.frame(coefficients(summary(fit)))
            names(cfit) = c('Est','SE','t','p')
            cmat[c1,c2] = cfit[val,'Est']
            pmat[c1,c2] = cfit[val,'p']
        }
    }
}

diag(pmat) = 1
diag(cmat) = 0
lpmat = -log10(pmat)
lcmat = (pmat < 0.05) * cmat
# Heatmap(lpmat, col=colb,
#         cluster_rows=FALSE,
#         cluster_columns=FALSE)

Heatmap(lcmat,
        cluster_rows=FALSE,
        cluster_columns=FALSE)

cdf = data.frame(cmat)
cdf$c1 = rownames(cmat)
ddf = gather(cdf, c2, val,-c1)
pdf = data.frame(pmat)
pdf$c1 = rownames(pmat)
ddf = merge(ddf, gather(pdf, c2, p,-c1))

nodes = topec2

vuln = c('Exc_RELN_COL5A2',
         'Exc_TOX3_TTC6',
         'Exc_AGBL1_GPC5',
         'Exc_RELN_GPC5')

# Simple network: just the links/points:
pcols = brewer.pal(12, 'Paired')
library(igraph)
pcut = 0.2
sdf = ddf[ddf$val > 0 & ddf$p < pcut,]
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=T) 
vcol = rep(tsp.col('grey75',.5), length(nodes))
vcol[nodes %in% vuln] = tsp.col(pcols[6], .5)
# Node color vs. node border:
expr = ncmat['HS6ST3',]
# expr = ncmat['CDH20',]
names(expr) = gsub("[,)(]","", gsub(" ", "_", names(expr)))
col_fun = colorRamp2(seq(min(expr), max(expr), length.out=100), viridis(100))
fcol = col_fun(expr)[nodes]
# ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
ecol = 'grey'
V(net)$size = 20
V(net)$label = nodes
V(net)$label.cex = .5
V(net)$label.color = 'black'
# V(net)$color = vcol
V(net)$color = fcol
V(net)$frame.color <- vcol
# V(net)$frame.color <- NA
V(net)$pch = 19
E(net)$color = ecol 
elty = rep('dotted', length(sdf$val))
elty[sdf$sim >= .85] = 'dashed'
elty[sdf$sim >= .95] = 'solid'
E(net)$lty = elty
E(net)$width = sqrt(sdf$val)
# E(net)$weight = sdf$sim * .5
E(net)$weight = sqrt(sdf$val) 
# E(net)$weight = sdf$sim * .5
# set.seed(2)
set.seed(8)
l <- layout_with_fr(net, grid='nogrid') # Usually best

# pdf(paste0(imgpref, npref, '.pdf'), width=6, height=6)
png(paste0(imgpref, 'EC_frac_corr_by_',advar,'.png'), res=300, units='in',width=4, height=4)
# pdf(paste0(imgpref, 'EC_frac_corr_by_',advar,'.pdf'), width=4, height=4)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
dev.off()
# plot(net, layout=l, edge.curved=seq(-0.5, 0.5, length = ecount(net)))
# dev.off()









# geom_text(data=labdf, aes(x=celltype, y=max(ctdf$frac), label=padj), color='black') + 
# stat_compare_means(label='p.format', method.args = list(alternative = "greater")) + 
# scale_y_continuous(label=scales::percent, expand=c(0,0.001)) + 





g0a = ggplot(ctdf, aes(x=celltype, y=count / total, color=nrad)) +
    # g0a = ggplot(ctdf, aes(x=celltype, y=count / exc.total, color=nrad)) +
    geom_boxplot() + 
    theme_pubr() + 
    geom_text(data=labdf, aes(x=celltype, y=max(ctdf$frac), label=padj), color='black') + 
    # stat_compare_means(label='p.format', method.args = list(alternative = "greater")) + 
    scale_color_manual(values=colvals[['nrad']]) + 
    scale_y_continuous(label=scales::percent, expand=c(0,0.001)) + 
    coord_flip()

g0b = ggplot(ctdf, aes(x=celltype, y=count / total, color=cogdxad)) +
# g0b = ggplot(ctdf, aes(x=celltype, y=count / exc.total, color=cogdxad)) +
    geom_boxplot() + 
    theme_pubr() + 
    stat_compare_means(label='p.signif') +
    scale_color_manual(values=colvals[['nrad']]) + 
    scale_y_continuous(label=scales::percent, expand=c(0,0.001)) + 
    coord_flip()

garr = g0a + g0b
ggsave(paste0(imgpref, 'ECfrac_againststrat_overall.png'), garr, dpi=450, units='in', width=12, height=6)
ggsave(paste0(imgpref, 'ECfrac_againststrat_overall.pdf'), garr, dpi=450, units='in', width=12, height=6)


g0c = ggplot(ctdf, aes(x=celltype, y=count / total, color=nrad)) +
    # g0a = ggplot(ctdf, aes(x=celltype, y=count / exc.total, color=nrad)) +
    facet_wrap(~Apoe_e4) + 
    geom_boxplot(outlier.shape=NA) + 
    theme_pubr() + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
    # geom_text(data=labdf, aes(x=celltype, y=max(ctdf$frac), label=padj), color='black') + 
    # stat_compare_means(label='p.format', method.args = list(alternative = "greater")) + 
    scale_color_manual(values=colvals[['nrad']]) + 
    scale_y_continuous(label=scales::percent, expand=c(0,0.001)) + 
    coord_flip()

fit = lm(frac ~ Apoe_e4 * celltype, ctdf[ctdf$nrad == 'AD',])

g0c = ggplot(ctdf[ctdf$nrad == 'AD',], aes(x=celltype, y=count / total, color=Apoe_e4)) +
    geom_boxplot(outlier.shape=NA) + 
    theme_pubr() + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
    # geom_text(data=labdf, aes(x=celltype, y=max(ctdf$frac), label=padj), color='black') + 
    stat_compare_means(label='p.format', method.args = list(alternative = "greater")) + 
    scale_color_manual(values=c('grey','black')) + 
    scale_y_continuous(label=scales::percent, expand=c(0,0.001)) + 
    coord_flip()



ctdf$frac = ctdf$count / ctdf$total
g1 = ggplot(ctdf, aes(x=nft_ec, y=frac)) +
    facet_wrap(~celltype, scales='free_y') + 
    geom_point() + geom_smooth(method='lm') + 
    theme_pubr() + 
    stat_fit_glance(method = 'lm', geom = 'text',
                    aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                    label.x = 'center', label.y = 'top', size = 3)

g2 = ggplot(ctdf, aes(x=plaq_n_ec, y=frac)) +
    facet_wrap(~celltype, scales='free_y') + 
    geom_point() + geom_smooth(method='lm') + 
    theme_pubr() + 
    stat_fit_glance(method = 'lm', geom = 'text',
                    aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                    label.x = 'center', label.y = 'top', size = 3)

g3 = ggplot(ctdf, aes(x=plaq_d_ec, y=frac)) +
    facet_wrap(~celltype, scales='free_y') + 
    geom_point() + geom_smooth(method='lm') + 
    theme_pubr() + 
    stat_fit_glance(method = 'lm', geom = 'text',
                    aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                    label.x = 'center', label.y = 'top', size = 3)

garr = g1 + g2 + g3
ggsave(paste0(imgpref, 'ECfrac_againstpath_overall.png'), garr, dpi=450, units='in', width=25, height=8)
ggsave(paste0(imgpref, 'ECfrac_againstpath_overall.pdf'), garr, dpi=450, units='in', width=25, height=8)


# -------------------------------------------------------
# Look at each other neuronal subtype; vs. nft, vs. nrad.
# -------------------------------------------------------
ctdf = agg.rename(barcode ~ projid + region + cell_type_high_resolution + major.celltype, cellmeta, length, 'count')
combdf = expand.grid(cell_type_high_resolution=unique(cellmeta$cell_type_high_resolution), region=unique(cellmeta$region), projid=unique(cellmeta$projid))
combdf = merge(combdf, unique(cellmeta[,c('cell_type_high_resolution','major.celltype')]))
ctdf = merge(ctdf, combdf, all.y=TRUE)
ctdf$count[is.na(ctdf$count)] = 0
ctdf = merge(ctdf, totdf) # Removes a couple missing projid x region comb.
ctdf$other = ctdf$total - ctdf$count
names(ctdf)[3] = 'celltype'
ctdf = merge(ctdf, unique(metadata[,c('projid','rind','region','nft_ec','plaq_d_ec','plaq_n_ec', 'braaksc','cogdx','niareagansc','msex','age_death','pmi', 'Apoe_e4', 'nrad','cogdxad')]))
ctdf$projid = factor(ctdf$projid)

mfracdf = aggregate(cbind(count, total) ~ region + celltype + nrad, ctdf, sum)
mfracdf$tot.frac = mfracdf$count / mfracdf$total

# Are specific celltypes altered?
subdf = ctdf[grep("Mic", ctdf$celltype),]
subdf = merge(subdf, pqdf, all.x=TRUE)
subdf$frac = subdf$count / subdf$total
subdf = merge(subdf, mfracdf[mfracdf$nrad == 'CTRL',c('region','celltype','tot.frac')])

g0a = ggplot(subdf, aes(x=celltype, y=frac, color=nrad)) +
    # facet_wrap(~region) + 
    geom_boxplot() + 
    theme_pubr() + 
    stat_compare_means(label='p.signif') + 
    scale_color_manual(values=colvals[['nrad']]) + 
    coord_flip()

g1 = ggplot(subdf, aes(x=nft, y=frac / tot.frac)) +
    # facet_wrap(region~celltype, scales='free_y') + 
    # facet_grid(region~celltype, scales='free_y') + 
    facet_grid(celltype~region, scales='free_y') + 
    geom_hline(yintercept=1, lty='dashed') + 
    geom_point() + geom_smooth(method='lm') + 
    theme_pubr() + 
    stat_fit_glance(method = 'lm', geom = 'text',
                    aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                    label.x = 'center', label.y = 'top', size = 3)

pathval = 'nrad'
# pathval = 'cogdxad'
asform = function(x){as.formula(paste(x, collapse=" "))}
subdf$cls = subdf$celltype
formula = asform(c('cbind(count, other) ~ cls *',pathval,'+ cls*pmi + cls*msex + cls* age_death + cls*region'))
fit <- glm(formula = formula, family = 'quasibinomial', data = subdf)

emform = asform(c('revpairwise ~', pathval, '|cls '))
emm1 <- emmeans(fit, specs=emform)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results

# head(subdf[order(subdf$frac / subdf$tot.frac, decreasing=T),])


# ---------------------------------------
# Plot the separate UMAPs for EC neurons:
# ---------------------------------------
tab = read.delim('multiRegion/metadata_test_mrad_Exc_EC_combat_filthvg.tsv', header=T)
tab = tab[,c('barcode','projid','celltype','U1','U2')]
tab = merge(tab, unique(metadata[,c('projid','nrad')]))

# Shuffle points for plotting:
ind = 1:nrow(tab)
ad.ind = which(tab$nrad == 'AD')
ctrl.ind = which(tab$nrad == 'CTRL')
full.ind = sample(ind, length(ind), replace=FALSE)
ad.ind = sample(ad.ind, length(ad.ind), replace=FALSE)
ctrl.ind = sample(ctrl.ind, length(ctrl.ind), replace=FALSE)
xr = diff(range(tab$U1))
yr = diff(range(tab$U2))

# Plot the umap, with appropriate proportions:
w = 3
h = w * yr / xr
cex = 0.025
sp = 0.1
png(paste0(imgpref, 'EC_vuln_neurons_umap_alone.png'), units='in', res=450, width=w, height=h)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(tab$U1[full.ind], tab$U2[full.ind], 
     col=tcols[tab$celltype[full.ind]], pch=19, cex=cex, axes=F)
dev.off()

png(paste0(imgpref, 'EC_vuln_neurons_umap_alone_CTRL.png'), units='in', res=450, width=w, height=h)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(tab$U1[ctrl.ind], tab$U2[ctrl.ind], 
     col=tcols[tab$celltype[ctrl.ind]], pch=19, cex=cex, axes=F)
dev.off()

png(paste0(imgpref, 'EC_vuln_neurons_umap_alone_AD.png'), units='in', res=450, width=w, height=h)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(tab$U1[ad.ind], tab$U2[ad.ind], 
     col=tcols[tab$celltype[ad.ind]], pch=19, cex=cex, axes=F)
dev.off()

# Load the TOX3 TTC6 regression results:
rdafile = paste0(mtxdir, 'all_brain_regions_filt_preprocessed_scanpy.majorcelltype.Exc.Exc_TOX3_TTC6.EC.rda')
load(rdafile)
celltype = 'Exc'
ststr = 'Exc_TOX3_TTC6'
path = 'nrad'
region = 'EC'
fpref = paste0(prefix, '.mastlmm_reg.', path, '.', region, '.major.', celltype, '.minor.', ststr)
regdir = paste0(datadir,'dereg/')
fnlist = list.files(pattern=paste0(fpref,'.*.Rda'), path=regdir)
mastdf = c()
for (fn in fnlist){
    print(fn)
    load(paste0(regdir,fn))
    mastdf = rbind(mastdf, regdf)
}
dim(mastdf)

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

submeta = cellmeta[cellmeta$major.celltype == 'Exc' & cellmeta$region == 'EC',]
submeta = merge(submeta, unique(metadata[,c('projid','cogdx','niareagansc')]))
ctrl.bcs = submeta$barcode[submeta$niareagansc %in% c(3,4)]
bind = bcs %in% ctrl.bcs

# Aggregate to the level of individual x cell type:
rownames(submeta) = submeta$barcode
smeta = submeta[colnames(amat), c('projid','cell_type_high_resolution', 'barcode')]
ptype = paste0(smeta$projid, '_', gsub("/","_",gsub(" ","_", smeta$cell_type_high_resolution)))
tform = make.tform(ptype, u=sort(unique(ptype)), norm=T)
apmat = amat %*% tform #raw
pmat = nmat %*% tform #norm
gc()

# One-hot index of selected celltypes:
eid = c(1,0,0,1,1,0,0,0,1)
names(eid) = topec

# --------------------------------------------------------------
# Plot the TOX3 TTC6 differential genes on the average matrices:
# --------------------------------------------------------------
names(mastdf) = c('symbol','p','Est','ci.hi','ci.lo','fdr')
mastdf = mastdf[order(mastdf$p),]
ntop = 10
upgene = head(mastdf[mastdf$Est > 0,'symbol'], ntop)
dwgene = head(mastdf[mastdf$Est < 0,'symbol'], ntop)
# upgene = c('PRNP', 'CHD3','MAP1A','ADNP','SYT11',
#            'HSP90B1','HSP90AB1','ST8SIA3','SV2A','FTL')

smeta = submeta[submeta$cell_type_high_resolution %in% topec,] 
umeta = agg.rename(barcode ~ projid + cell_type_high_resolution + region, smeta, length, 'count')
umeta$ptype = with(umeta, paste0(projid, '_', gsub("/","_",gsub(" ","_", cell_type_high_resolution))))
umeta = umeta[umeta$count > 10,]
umeta = umeta[umeta$cell_type_high_resolution != 'Exc SOX11 NCKAP5',]

# Plot heatmap by braak stage / nft:
pgenes = c(upgene, dwgene)
pset = c(rep('Up',ntop), rep('Down',ntop))
use.norm=TRUE
if (use.norm){
    smat = as.matrix(pmat[pgenes,umeta$ptype])
} else {
    smat = as.matrix(apmat[pgenes,umeta$ptype])
}


# Individual + CT annotations:
indmap = unique(metadata[,c('projid','cogdx','niareagansc', 'braaksc','nft_ec', 'nrad','cogdxad')])
rownames(indmap) = as.character(indmap$projid)
cn = colnames(smat)
pids = sub("_.*","",cn)
pmeta = indmap[as.character(pids),]
cts = gsub("_"," ", sub(".*_Exc","Exc",cn))
ad.diff = ifelse(cts %in% topec[eid==1],'Depleted','No Change')
clsplit = paste0(pmeta$nrad, ", ", ad.diff)
nft.col_fun = colorRamp2(range(pmeta$nft_ec), c("white", "indianred"))
ha = HeatmapAnnotation(AD.Diff=ad.diff, CT=cts, 
                       Braak=pmeta$braaksc,
                       NFT=pmeta$nft_ec,
                       AD=pmeta$niareagansc,
                       col=list(AD=colvals[['niareagansc']],
                                Braak=colvals[['braaksc']],
                                CT=type.cols[topec],
                                NFT=nft.col_fun,
                                AD.Diff=c('Depleted'='slateblue',
                                          'No Change'='grey75')))

udsplit = pset
hb = rowAnnotation(Set=pset,
                   col=list(Set=c('Down'='blue','Up'='indianred')))

if (use.norm){
    # png(paste0(imgpref, 'EC_tox3ttc6_DEgenes_heatmap_normalized_individ.png'), res=400, units='in', width=11, height=5)
    pdf(paste0(imgpref, 'EC_tox3ttc6_DEgenes_heatmap_normalized_individ.pdf'), width=11, height=4)
} else {
    # png(paste0(imgpref, 'EC_tox3ttc6_DEgenes_heatmap_unnormalized_individ.png'), res=400, units='in', width=11, height=5)
    pdf(paste0(imgpref, 'EC_tox3ttc6_DEgenes_heatmap_unnormalized_individ.pdf'), width=11, height=4)
}
smat.scaled = t(scale(t(log(smat+1)), center=TRUE))
Heatmap(smat.scaled, name='scaled\n logcounts', 
        use_raster=TRUE,
        top_annotation=ha, 
        column_split=clsplit, 
        show_column_names=FALSE,
        row_split=udsplit,
        right_annotation=hb)
dev.off()

# -----------------------------------------------------
# Make more systematic calls of what is signif vs. not:
# -----------------------------------------------------
ctdf$frac = ctdf$count / ctdf$total
ctdf = merge(ctdf, pqdf, all.x=TRUE)
ctdf = merge(ctdf, mfracdf[mfracdf$nrad == 'CTRL',c('region','celltype','tot.frac')])
ctdf$rel.frac = ctdf$frac / ctdf$tot.frac
ctdf$rel.frac[ctdf$tot.frac == 0] = 0
subtypes = unique(ctdf$celltype)
regdf = c()  # NOTE: Should we do a quasibinomial model instead?
for (subtype in subtypes){
    cat(subtype)
    subdf = ctdf[ctdf$celltype == subtype,]
    for (region in unique(subdf$region)){
        cat('\t', region)
        fit1 = lm(frac ~ nrad, subdf[subdf$region == region,])
        fit1b = lm(rel.frac ~ nrad, subdf[subdf$region == region,])
        if (region != 'TH'){
            fit2 = lm(frac ~ nft, subdf[subdf$region == region,])
            fit2b = lm(rel.frac ~ nft, subdf[subdf$region == region,])
            cfit2 = coefficients(summary(fit2))
            cfit2b = coefficients(summary(fit2b))
        }
        fit3 = lm(frac ~ cogdxad, subdf[subdf$region == region,])
        fit3b = lm(rel.frac ~ cogdxad, subdf[subdf$region == region,])
        cfit1 = coefficients(summary(fit1))
        cfit1b = coefficients(summary(fit1b))
        cfit3 = coefficients(summary(fit3))
        cfit3b = coefficients(summary(fit3b))
        if (region != 'TH'){
            df = data.frame(rbind(cfit1['nradAD',c(1,4)], cfit1b['nradAD',c(1,4)], 
                                  cfit2['nft',c(1,4)], cfit2b['nft',c(1,4)],
                                  cfit3['cogdxadAD',c(1,4)], cfit3b['cogdxadAD',c(1,4)]))
            df$reg = c('nrad','nrad.rel', 'nft', 'nft.rel','cogdxad','cogdxad.rel')
        } else {
            df = data.frame(rbind(cfit1['nradAD',c(1,4)], cfit1b['nradAD',c(1,4)], 
                                  cfit3['cogdxadAD',c(1,4)], cfit3b['cogdxadAD',c(1,4)]))
            df$reg = c('nrad','nrad.rel', 'cogdxad','cogdxad.rel')
        }
        colnames(df) = c('Est','p','reg')
        df$celltype = subtype
        df$region = region
        regdf = rbind(regdf, df)
    }
    cat('\n')
}

rdf = regdf[regdf$reg == 'nrad',]
rdf = regdf[regdf$reg == 'nft',]
rdf = regdf[regdf$reg == 'nft.rel',]
rdf = rdf[order(-abs(rdf$Est)),]

rdf = regdf[regdf$reg == 'cogdxad',]
rdf = rdf[order(rdf$p),]
rdf = regdf[regdf$reg == 'cogdxad.rel',]
rdf = rdf[order(-abs(rdf$Est)),]


# ----------------------------------------
# Create average matrices; calculate corr:
# ----------------------------------------
tform = make.tform(lbs, u=topec, norm=TRUE)
tform.ctrl = make.tform(lbs[bind], u=topec, norm=TRUE)
tform.ad = make.tform(lbs[!bind], u=topec, norm=TRUE)
# Normalized to median or not:
use.norm=TRUE
if (!use.norm){
    avgmat = as.matrix(amat %*% tform)
    avgmat.ctrl = as.matrix(amat[,bind] %*% tform.ctrl)
    avgmat.ad = as.matrix(amat[,!bind] %*% tform.ad)
} else {
    avgmat = as.matrix(nmat %*% tform)
    avgmat.ctrl = as.matrix(nmat[,bind] %*% tform.ctrl)
    avgmat.ad = as.matrix(nmat[,!bind] %*% tform.ad)
}

# Calculate the correlation between the averages + vulnerable cells
cv = cor(t(avgmat.ctrl), t(t(eid)))
cdf = data.frame(cor=cv, im=apply(avgmat.ctrl[,eid ==1],1, mean),
                 om=apply(avgmat.ctrl[,eid == 0],1, mean), 
                 tm=apply(avgmat.ctrl,1, mean))
cdf = cdf[!is.na(cdf$cor),]
cdf = cdf[order(cdf$cor, decreasing=T),, drop=FALSE]
cdf$symbol = rownames(cdf)
cdf$diff = cdf$im - cdf$om

labdf = rbind(head(cdf[cdf$im > 1,], 20),
              tail(cdf[cdf$om > 1,], 20))

# Plot MA plot - average/difference
g3 = ggplot(cdf, aes(log2(im+om), log2(im/om))) +
    geom_point(col='grey85', cex=.5) + 
    geom_smooth(method='gam') + 
    geom_point(data=labdf, aes(log2(im+om), log2(im/om)), cex=.5) +
    geom_text_repel(data=labdf, aes(log2(im+om), log2(im/om), label=symbol), max.overlaps=20) +
    theme_pubr()
ggsave(paste0(imgpref, 'ECcorr_differences_basic_maplot.png'), g3, dpi=450, units='in', width=8, height=6)
ggsave(paste0(imgpref, 'ECcorr_differences_basic_maplot.pdf'), g3, dpi=450, units='in', width=8, height=6)


# Plot:
lab2df = rbind(head(cdf[cdf$im > 1,], 20),
               tail(cdf[cdf$om > 1,], 20))
smat = avgmat[lab2df$symbol,]
smat = log(smat + 1)
smat = t(scale(t(smat)))

ad.diff = ifelse(eid==1,'Depleted','No Change')
udsplit = ifelse(lab2df$cor > 0 ,'Up','Down')
ha = HeatmapAnnotation(AD.Diff=ad.diff, 
                       col=list(AD.Diff=c('Depleted'='slateblue',
                                          'No Change'='grey75')))

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "indianred"))
# ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
hb = rowAnnotation(Cor=lab2df$cor,
                   col=list(Cor=col_fun))

if (use.norm){
    png(paste0(imgpref, 'EC_topcorr_heatmap_normalized.png'), res=400, units='in', width=8, height=10)
} else {
    png(paste0(imgpref, 'EC_topcorr_heatmap_unnormalized.png'), res=400, units='in', width=8, height=10)
}
Heatmap(smat, name='log counts', column_names_rot=45, 
        # col=viridis(100), 
        top_annotation=ha, column_split=ad.diff, 
        row_split=udsplit,
        right_annotation=hb)
dev.off()

# Plot heatmap by braak stage / nft:
lab2df = rbind(head(cdf[cdf$im > 1,], 50),
               tail(cdf[cdf$om > 1,], 20))
if (use.norm){
    smat = as.matrix(pmat[lab2df$symbol,])
} else {
    smat = as.matrix(apmat[lab2df$symbol,])
}

# Individual + CT annotations:
indmap = unique(metadata[,c('projid','cogdx','niareagansc', 'braaksc','nft_ec', 'nrad','cogdxad')])
rownames(indmap) = as.character(indmap$projid)
cn = colnames(pmat)
pids = sub("_.*","",cn)
pmeta = indmap[as.character(pids),]
cts = gsub("_"," ", sub(".*_Exc","Exc",cn))
ad.diff = ifelse(cts %in% topec[eid==1],'Depleted','No Change')
clsplit = paste0(pmeta$nrad, ", ", ad.diff)
ha = HeatmapAnnotation(AD.Diff=ad.diff, CT=cts, 
                       Braak=pmeta$braaksc,
                       NFT=pmeta$nft_ec,
                       AD=pmeta$niareagansc,
                       col=list(AD=colvals[['niareagansc']],
                                Braak=colvals[['braaksc']],
                                CT=type.cols[topec],
                                AD.Diff=c('Depleted'='slateblue',
                                          'No Change'='grey75')))


col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "indianred"))
udsplit = ifelse(lab2df$cor > 0 ,'Up','Down')
hb = rowAnnotation(Cor=lab2df$cor,
                   col=list(Cor=col_fun))

if (use.norm){
    png(paste0(imgpref, 'EC_topcorr_heatmap_normalized_individ.png'), res=400, units='in', width=11, height=12)
} else {
    png(paste0(imgpref, 'EC_topcorr_heatmap_unnormalized_individ.png'), res=400, units='in', width=11, height=12)
}
smat.scaled = t(scale(t(log(smat+1))))
Heatmap(smat.scaled, name='scaled\n logcounts', 
        top_annotation=ha, 
        column_split=clsplit, 
        show_column_names=FALSE,
        row_split=udsplit,
        right_annotation=hb)
dev.off()


# --------------------------------------------------------
# Plot the heparan sulfate associated genes in particular:
# --------------------------------------------------------
pathway='heparan_sulfate'
# pathway='inositol'
# pathway='oglycan_biosynthesis'
hsgenes = scan(paste0(pathway,'_genes.txt'),'c')
hsdf = read.delim(paste0(pathway, '_ECde.txt'), header=T)
# hsgenes = scan('inositol_genes.txt','c')
rownames(hsdf) = hsdf$symbol
# Plot heatmap by braak stage / nft:
kgenes = hsgenes[hsgenes %in% rownames(pmat)]
if (use.norm){
    smat = as.matrix(pmat[kgenes,])
} else {
    smat = as.matrix(apmat[kgenes,])
}
smat = smat[apply(smat,1,sum) != 0,]
kgenes = rownames(smat)

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "indianred"))
udsplit = ifelse(hsdf[kgenes,'Est'] > 0 ,'Up','Down')
hb = rowAnnotation(Est=hsdf[kgenes,'Est'],
                   log10q=hsdf[kgenes,'log10q'],
                   Signif=ifelse(hsdf[kgenes,'padj'] < 0.001, 'Yes','No'),
                   col=list(Est=col_fun,
                            Signif=c('Yes'='goldenrod1','No'='white')))


if (use.norm){
    imgfile = paste0(imgpref, 'EC_',pathway,'_heatmap_normalized_individ')
} else {
    imgfile = paste0(imgpref, 'EC_',pathway,'_heatmap_unnormalized_individ')
}
png(paste0(imgfile,'.png'), res=400, units='in', width=12, height=10)
# pdf(paste0(imgfile,'.pdf'), width=12, height=10)
smat.scaled = t(scale(t(log(smat+1))))
Heatmap(smat.scaled, name='scaled\n logcounts', 
        top_annotation=ha, 
        column_split=clsplit, 
        row_split=udsplit,
        show_column_names=FALSE,
        right_annotation=hb
)
dev.off()

# -----------------------------------------------------------
# Create a regression score for the heparan sulfate gene set:
# -----------------------------------------------------------
sdf = data.frame(t(as.matrix(pmat[hsgenes,])))
cn = rownames(sdf)
cts = gsub("_"," ", sub(".*_Exc","Exc",cn))
sdf$celltype = cts
rdf = regdf[regdf$region == 'EC' & regdf$reg == 'nrad',]
sdf = merge(sdf, rdf[,c('celltype','Est')])
sdf$celltype=NULL
hs.fit = lm(Est ~ ., sdf)
summary(hs.fit)
# HS6ST3 is most signif, followed by GPC5, SDC2/3/4, NDST3, and HS3ST3A1


# -----------------------------------------------
# Plot boxplots/violins of the genes of interest:
# -----------------------------------------------
pltgenes = c('HS6ST1','HS6ST2','HS6ST3','NDST1','NDST3','NDST4','EXTL1','EXTL2','GPC5','GPC6','SULF1', 'SULF2', 'INPP4B','SLC2A13')
sdf = data.frame(as.matrix(pmat[pltgenes,]))
sdf$symbol = rownames(sdf)
slong = gather(sdf, ptype, val, -symbol)
slong$ptype = sub("^X","", slong$ptype)
hdf = data.frame(ptype=cn,
                 ad=pmeta$niareagansc, ct=cts, 
                 ad.diff=ad.diff)
slong = merge(slong, hdf)

gplot = ggplot(slong, aes(ad.diff, log(val + 1), fill=ad.diff)) + 
    facet_wrap(~symbol, scales='free_y') + 
    scale_fill_manual(values=c('slateblue','grey85')) +
    geom_violin(alpha=.5) + 
    geom_boxplot(width=.25, outlier.shape=NA, alpha=.5) + 
    geom_jitter(cex=.25) + 
    # stat_compare_means() + 
    labs(x='Individuals x Subtypes depleted in AD or not', y='log(Expression + 1)') +
    theme_pubr() + theme(legend.position='none')

if (use.norm){
    imgfile = paste0(imgpref, 'EC_heparan_boxplots_normalized_individ')
} else {
    imgfile = paste0(imgpref, 'EC_heparan_boxplots_unnormalized_individ')
}
ggsave(paste0(imgfile, '.png'), gplot, dpi=450, units='in', width=8.5, height=7)
ggsave(paste0(imgfile, '.pdf'), gplot, dpi=450, units='in', width=8.5, height=7)


# ---------------------------------
# Run regression to find top genes:
# ---------------------------------
# Log2 + 1
pdf = data.frame(as.matrix(pmat))
pdf$symbol = rownames(pdf)
plong = gather(pdf, projid, val, -symbol)
plong$celltype = sub(".*_Exc","Exc",plong$projid)
plong$celltype = gsub("_"," ",plong$celltype)
plong$projid = sub("^X","",plong$projid)
plong$projid = sub("_.*","",plong$projid)
plong$val = log(plong$val + 1)
plong$vuln = 1 * (plong$celltype %in% topec[eid == 1])

nrmeta = unique(ctdf[,c('projid','nrad','cogdxad','nft_ec','plaq_n_ec', 'Apoe_e4','msex','pmi','age_death')])
rownames(nrmeta) = nrmeta$projid
plong$nrad = nrmeta[as.character(plong$projid),'nrad']
plong$nft_ec = nrmeta[as.character(plong$projid),'nft_ec']
plong$cogdxad = nrmeta[as.character(plong$projid),'cogdxad']
plong$pmi = nrmeta[as.character(plong$projid),'pmi']
plong$age_rescaled = nrmeta[as.character(plong$projid),'age_death'] / 100
plong$msex = nrmeta[as.character(plong$projid),'msex']

# Average expr:
avg.pdf = aggregate(val ~ symbol, plong, mean)

# Run the pseudobulk differential expression:
ecvuln.regfile.rda = 'multiRegion/EC_vulnerability_associated_reg_pseudobulk.Rda'
if (!file.exists(ecvuln.regfile.rda)){
    genelist = unique(plong$symbol)
    resdf = c()
    # for (gene in kgenes){
    for (gene in genelist){
        sdf = plong[plong$symbol == gene,]
        fit = glm(vuln ~ val * nrad + val * nft_ec + age_rescaled + msex + pmi,sdf, family='gaussian')
        # print(summary(fit))
        cfit = coefficients(summary(fit))
        df = data.frame(cfit)
        colnames(df) = c('Est','SE','t','p')
        pval = df['val','p']
        if (!is.na(pval)){
            if (pval < 1e-4){ cat(gene,'\t', sprintf('%0.2e',pval), '\n') }
        }
        df$var = rownames(df)
        rownames(df) = NULL
        df$symbol = gene
        resdf = rbind(resdf, df)
    }
    save(resdf, file=ecvuln.regfile.rda)
} else {
    load(ecvuln.regfile.rda)
}

# resdf$var = '(Intercept)'
# resdf$var[grep('vuln', rownames(resdf))] = 'vuln'
# resdf$var[grep('nradAD', rownames(resdf))] = 'nradAD'
# resdf$var[grep('vuln:nradAD', rownames(resdf))] = 'vuln:nradAD'

# Look at the top genes:
adf = resdf[resdf$var == 'val',c('symbol','Est','p')]
adf = merge(adf, avg.pdf[avg.pdf$val > 1,])  # Require some level expr:
adf = adf[!is.na(adf$p),]
adf = adf[order(adf$p), ]
adf$padj = p.adjust(adf$p, 'fdr')
adf$log10q = -log10(adf$padj)
pcut = 0.001
adf$color = 0
adf$color[adf$padj < pcut] = 1
adf$color[adf$padj < pcut & adf$Est > 0] = 2

# write.table(adf, 'heparan_sulfate_ECde.txt',quote=F, row.names=F, sep="\t")
# write.table(adf, 'inositol_ECde.txt',quote=F, row.names=F, sep="\t")

labdf = rbind(head(adf[adf$color == 1,],5),
              head(adf[adf$color == 2,],10))

gplot = ggplot(adf, aes(Est, log10q, col=factor(color))) + 
    scale_color_manual(values=c('grey85','royalblue','indianred')) + 
    geom_vline(xintercept=0, lwd=.25, lty='dashed') + 
    geom_point(cex=.5, alpha=1) + theme_pubr() + 
    labs(x='',y='') +
    scale_y_continuous(expand=c(0,0), lim=c(0,max(labdf$log10q) + 1)) + 
    theme(legend.position = 'none')

gplab = gplot + geom_text_repel(data=labdf, aes(Est, log10q, label=symbol), max.overlaps=30, size=2.2)
w = 2.25; h=2.75
gpnoaxes = gplot + theme(axis.text.x=element_blank(), 
                         axis.ticks.y=element_line(size=.1),
                         axis.ticks.x=element_line(size=.1),
                         axis.line=element_blank(),
                         axis.text.y=element_blank())
ggsave(paste0(imgpref, 'EC_pseudobulk_vuln_genes_volcano.png'), gplab, dpi=450, units='in', width=w, height=h)
ggsave(paste0(imgpref, 'EC_pseudobulk_vuln_genes_volcano_noaxes.png'), gpnoaxes, dpi=450, units='in', width=w, height=h)
ggsave(paste0(imgpref, 'EC_pseudobulk_vuln_genes_volcano.pdf'), gplab, dpi=450, units='in', width=w, height=h)

# NOTE: No interacting genes are significantly differential.

smat = avgmat[head(adf$symbol,30),]
smat = log(smat + 1)
Heatmap(smat, column_names_rot=45, col=viridis(100))


# ---------------------------------------------------------------------
# Differential gene expression on vulnerable + non-vulnerable subtypes:
# ---------------------------------------------------------------------
# Pseudo-bulk metadata:
smeta = submeta[submeta$cell_type_high_resolution %in% topec,] 
umeta = agg.rename(barcode ~ projid + cell_type_high_resolution + region, smeta, length, 'count')
umeta$ptype = with(umeta, paste0(projid, '_', gsub("/","_",gsub(" ","_", cell_type_high_resolution))))
umeta = merge(umeta, unique(metadata[,c('projid','region','braaksc','cogdx','niareagansc','msex','age_death','pmi', 'Apoe_e4', 'nrad','cogdxad')]))
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
umeta$vuln = umeta$cell_type_high_resolution %in% topec[eid==1]
umeta = umeta[colnames(pmat),]
umeta = umeta[umeta$count > 10,] # Remove + will weight

# Filter by avg. expr:
ecut = 0.5
pmean = apply(pmat, 1, mean)
kept.genelist = names(pmean)[pmean > ecut]
avgidf = data.frame(symbol=names(pmean), val=pmean)

gene = 'PRNP'
est.ddf = c()
for (gene in kept.genelist){
    x = pmat[gene, umeta$ptype]
    umeta$val = x
    # Regression - general effects of covariates on expression:
    fit = glm(val ~ nrad + cell_type_high_resolution + Apoe_e4 + age_rescaled + msex + pmi, umeta,
              weights=log(umeta$count), family='gaussian') # Corrects for cell ct. but not inflated
    # Alt:
    cfit = coefficients(summary(fit))
    df = data.frame(cfit)
    colnames(df) = c('Est','SE','t','p')
    pval = df['nradAD','p']
    if (!is.na(pval)){
        if (pval < 1e-2){ cat(gene,'\t', sprintf('%0.2e',pval), '\n') }
    }
    df$var = rownames(df)
    rownames(df) = NULL
    df$symbol = gene
    est.ddf = rbind(est.ddf, df)
}


# Plot volcano of these assoc. with depletions:
adf = est.ddf[est.ddf$var == 'nradAD',]
adf = merge(adf, avgidf)
adf = adf[!is.na(adf$p),]
adf = adf[adf$val > 1.5,]
adf = adf[order(adf$p), ]
adf$padj = p.adjust(adf$p, 'fdr')
adf$log10q = -log10(adf$padj)
adf$Est[adf$Est > 1] = 1
adf$Est[adf$Est < -1] = -1

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
ggsave(paste0(imgpref, 'EC_nrad_subtypes_pseudobulk_ct_genes_',sub(" ", "_", st), '_volcano.png'), gplot, dpi=450, units='in', width=w, height=h)
ggsave(paste0(imgpref, 'EC_nrad_subtypes_pseudobulk_ct_genes_',sub(" ", "_", st), '_volcano.pdf'), gplot, dpi=450, units='in', width=w, height=h)





# --------------------------------------------------------------
# Load in all neuronal cell types to look at heparan metabolism:
# --------------------------------------------------------------
# Load in average tracks for each cell type, regardless of individual:
avg.file.rda = 'multiRegion/neuronal_average_profiles.Rda'
if (!file.exists(avg.file.rda)){
    avgdf = c()
    for (celltype in c('Exc','Inh')){
        fns = list.files(path=mtxdir, paste0(rawpref,'.majorcelltype.', celltype,'..*rda'))
        for (fn in fns[62:length(fns)]){
            rdafile = paste0(mtxdir, fn)
            basename = sub(".rda", "", sub(paste0(".*.majorcelltype.", celltype,'.'),"",fn))
            subtype = gsub("_"," ",sub("\\.[A-Z]*","",basename))
            region = sub(".*\\.","",basename)
            # Load data, normalize and average:
            load(rdafile)
            if (nrow(mat) != 1 && ncol(mat) > 0){
                print(paste("[STATUS] Loaded", subtype, 'in', 
                            region,'with',ncol(mat), 'cells'))
                barcodes = colnames(mat)
                genes = rownames(mat)
                # Normalize:
                amarg = marg[barcodes,'count']
                fact = amarg / 15000  # Rough median for all neurons
                mat@x <- mat@x / rep.int(fact, diff(mat@p))
                # gc()
                # Average:
                avg.profile = rowMeans(mat)
                df = data.frame(major.celltype=celltype,
                                celltype=subtype,
                                region=region,
                                ncell=ncol(mat),
                                symbol=names(avg.profile),
                                val=avg.profile)
                avgdf = rbind(avgdf, df)
                cat(subtype,'\t', region,'\t',avg.profile[c('HS6ST1','HS6ST2','HS6ST3','SLC2A13')],'\n')
                rm(mat)
                gc()
            }
        }
    }
    # Also spread the matrix:
    genelist = unique(avgdf$symbol)
    avgdf$celltype = sub("L5 6 NP","L5/6 NP",avgdf$celltype)
    avgdf$celltype = sub("L5 6 IT","L5/6 IT",avgdf$celltype)
    avgwide = spread(avgdf, symbol, val)
    save(avgdf, avgwide, genelist, file=avg.file.rda)
} else {
    load(avg.file.rda)
}


# -----------------------------------------------------------------
# Plot subtype profiles in the heparan sulfate + inositol pathways:
# -----------------------------------------------------------------
# Heparan
# Inositol
# Differential (Re-run difftl...)
stmat = as.matrix(avgwide[,genelist])
pathway='heparan_sulfate'
# pathway='inositol'
hsgenes = scan(paste0(pathway,'_genes.txt'),'c')
hsdf = read.delim(paste0(pathway, '_ECde.txt'), header=T)
rownames(hsdf) = hsdf$symbol
# Plot heatmap by braak stage / nft:
kgenes = hsgenes[hsgenes %in% genelist]
smat = t(stmat[,kgenes])
smat = smat[apply(smat,1,sum) != 0,]
kgenes = rownames(smat)

colnames(smat) = paste0(avgwide$celltype, '; ', avgwide$region)

# Region + CT annotations:
ha = HeatmapAnnotation(Major=avgwide$major.celltype,
                       CT=avgwide$celltype, 
                       Region=avgwide$region,
                       col=list(Region=reg.cols,
                                CT=type.cols,
                                Major=c('Exc'='slateblue',
                                        'Inh'='grey75')))

# Gene annotations:
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "indianred"))
udsplit = ifelse(hsdf[kgenes,'Est'] > 0 ,'Up','Down')
hb = rowAnnotation(Est=hsdf[kgenes,'Est'],
                   log10q=hsdf[kgenes,'log10q'],
                   Signif=ifelse(hsdf[kgenes,'padj'] < 0.001, 'Yes','No'),
                   col=list(Est=col_fun,
                            Signif=c('Yes'='goldenrod1','No'='white')))

imgfile = paste0(imgpref, 'ECdriven_allsubtypes_',pathway,'_heatmap_normalized')
png(paste0(imgfile, '.png'), res=400, units='in', width=25, height=10)
smat.scaled = t(scale(t(log(smat+1))))
Heatmap(smat.scaled, name='scaled\n logcounts', 
        top_annotation=ha, 
        column_split=avgwide$major.celltype, 
        show_column_names=TRUE,
        column_names_gp = gpar(fontsize=7),
        row_split=udsplit,
        right_annotation=hb)
dev.off()

sdf = avgdf[avgdf$region == 'TH' & avgdf$symbol == 'HS6ST3',]
sdf = sdf[sdf$ncell > 100,]


# --------------------------------
# Plot genes against coefficients:
# --------------------------------
regset = 'nrad'
# regset = 'nrad.rel'
# regset = 'cogdxad'
rdf = regdf[regdf$reg == regset,]
if (regset == 'nrad'){
    ecut = 0.0075
} else if (regset == 'cogdxad.rel'){
    ecut = 0.01
} else if (regset == 'cogdxad'){
    ecut = 0.01
} else { 
    ecut = 0.0005
}

mfracdf = aggregate(cbind(count, total) ~ region + celltype, ctdf, sum)
mfracdf$frac = mfracdf$count / mfracdf$total

ecdriven=FALSE
if (ecdriven){
    genes = c('HS6ST1','HS6ST2','HS6ST3', 'SULF1', 'SULF2', 'GPC5') # Heparan + EC
} else {
    genes = c('HS6ST3','TRPS1','CDH20', 'MAN1A1','SGCD','OXR1', 'PIK3CB', 'HS6ST2') # Discovered across all
    # genes = c('HS6ST3','MAN1A1','SPOCK3', 'FILIP1L') # Discovered across all
}
subdf = avgdf[avgdf$symbol %in% genes,]
subdf = merge(subdf, rdf)
subdf = merge(subdf, mfracdf[,c('region','celltype','frac')])
subdf = subdf[subdf$frac > 0.005,] # Can include/ or not - removes noisy estimates
subdf = subdf[subdf$major.celltype == 'Exc',]
subdf = subdf[subdf$celltype != 'Exc SV2C LINC02137',]

# Add p-values, fit lines:
# regs = c('PFC','AG','MT')
regs = c('EC','HC','TH')
# regs = names(reg.cols)
if (length(grep('rel', regset)) ==  0){ subdf$rel.Est = subdf$Est / subdf$frac } else { subdf$rel.Est = subdf$Est }

lmdf = c()
for (gene in genes){
    sdf = subdf[subdf$symbol == gene & subdf$region %in% regs,]
    fit = lm(Est ~ val,sdf)
    cfit = coefficients(summary(fit))
    df = data.frame(cfit['val',c(1,4), drop=F])
    names(df) = c('lm.est','p')
    df$intercept = cfit['(Intercept)',1]
    df$symbol = gene
    lmdf = rbind(lmdf,df)
}
lmdf$eqn = sprintf("y = %0.2f%% - %0.2f%% * Expression\np = %0.3f",
                   round(lmdf$intercept * 100,2),
                   round(lmdf$lm.est * -100,2), round(lmdf$p,3))

labdf = subdf[abs(subdf$Est) > ecut | subdf$val > 30,]
gplot = ggplot(subdf[subdf$region %in% regs,], aes(val, Est, color=region, size=ncell)) +
    facet_wrap(~symbol, scale='free_x') +
    geom_hline(yintercept=0,lwd=.25, lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point() + 
    geom_text_repel(data=labdf[labdf$region %in% regs,], aes(val, Est, label=paste0(celltype, '; ', region)), max.overlaps=30, size=3) + 
    geom_text(data=lmdf, aes(x=1,y=-.02, label=eqn), size=3, color='black') + 
    scale_color_manual(values=reg.cols) + 
    scale_y_continuous(label=scales::percent)+
    labs(x='Average Expression Level (Normalized)', y='Absolute change vs. all cells') +
    theme_pubr()

if (ecdriven){
    ggsave(paste0(imgpref, 'ECdriven_allsubtypes_est_expr_scatter_comparison_',regset,'.png'), gplot, dpi=450, units='in', width=15, height=10)
    ggsave(paste0(imgpref, 'ECdriven_allsubtypes_est_expr_scatter_comparison_',regset,'.pdf'), gplot, dpi=450, units='in', width=15, height=10)
} else {
    ggsave(paste0(imgpref, 'fromall4_allsubtypes_est_expr_scatter_comparison_',regset,'.png'), gplot, dpi=450, units='in', width=7, height=7)
    ggsave(paste0(imgpref, 'fromall4_allsubtypes_est_expr_scatter_comparison_',regset,'.pdf'), gplot, dpi=450, units='in', width=7, height=7)
}


# ----------------------------------------------------------------
# Perform an EC/HC/TH screen for the best predictors of depletion:
# ----------------------------------------------------------------
# Load in all neuronal cell types in these regions, by individual:
indavg.file.rda = 'multiRegion/neuronal_exc_ECHCTH_ind_average_profiles.Rda'
rownames(cellmeta) = cellmeta$barcode
if (!file.exists(indavg.file.rda)){
    celltype = 'Exc'
    fns = c(list.files(path=mtxdir, paste0(rawpref,'.majorcelltype.', celltype,'..*.EC.rda')),
            list.files(path=mtxdir, paste0(rawpref,'.majorcelltype.', celltype,'..*.HC.rda')),
            list.files(path=mtxdir, paste0(rawpref,'.majorcelltype.', celltype,'..*.TH.rda')))
    indavgdf = c()
    for (fn in fns[1:length(fns)]){
        rdafile = paste0(mtxdir, fn)
        basename = sub(".rda", "", sub(paste0(".*.majorcelltype.", celltype,'.'),"",fn))
        subtype = gsub("_"," ",sub("\\.[A-Z]*","",basename))
        region = sub(".*\\.","",basename)
        # Load data, normalize and average:
        load(rdafile)
        if (nrow(mat) != 1 && ncol(mat) > 0){
            cat(paste("[STATUS] Loaded", subtype, 'in', region,'with',ncol(mat), 'cells\n'))
            barcodes = colnames(mat)
            genes = rownames(mat)
            # Normalize:
            amarg = marg[barcodes,'count']
            fact = amarg / 15000  # Rough median for all neurons
            mat@x <- mat@x / rep.int(fact, diff(mat@p))
            pid = cellmeta[barcodes,'projid']
            # Average across individuals:
            tform = make.tform(as.character(pid), norm=TRUE)
            pid.mat = mat %*% tform
            df = data.frame(as.matrix(pid.mat))
            colnames(df) = sub("^X","", colnames(df))
            df$symbol = rownames(pid.mat)
            df = gather(df, projid, val,-symbol)
            df$celltype = subtype 
            df$region = region
            # df$major.celltype = celltype
            # df$ncell=ncol(mat),
            indavgdf = rbind(indavgdf, df)
            rm(mat, df, pid.mat, tform)
            gc()
        }
    }

    # Also spread the matrix:
    genelist = unique(indavgdf$symbol)
    indavgdf$celltype = sub("L5 6 NP","L5/6 NP",indavgdf$celltype)
    indavgdf$celltype = sub("L5 6 IT","L5/6 IT",indavgdf$celltype)
    indavgwide = spread(indavgdf, symbol, val)
    save(indavgdf, indavgwide, genelist, file=indavg.file.rda)
} else {
    load(indavg.file.rda)
}
rm(indavgwide); gc()

# For filtering genes:
avgidf = aggregate(val ~ symbol, indavgdf, mean)

# ---------------------------------------------
# Merge with the estimates; perform regression:
# ---------------------------------------------
# For weighted regression:
ptdf = agg.rename(barcode ~ cell_type_high_resolution + projid + region, cellmeta, length, 'count')
names(ptdf)[1] = 'celltype'

# Faster indexing:
rdf$cr = paste0(rdf$celltype, '_',rdf$region)
ptdf$pcr = paste0(ptdf$projid,'_',ptdf$celltype, '_',ptdf$region)
rownames(rdf) = rdf$cr
rownames(ptdf) = ptdf$pcr
rownames(pqdf) = pqdf$rind

# gene = 'HS6ST3'
ecut = 0.5
kept.genelist = avgidf[avgidf$val > ecut, 'symbol']
est.regfile.rda = 'multiRegion/gene_fraction_reg_association_table.Rda'
if (!file.exists(est.regfile.rda)){
    est.regdf = c()
    est.expdf = c()
    # How many cells:
    NGENE = length(genelist) # Full genelist, for use in indexing.
    for (gene in kept.genelist){ 
        # Subset to gene:
        gidx = which(genelist == gene)
        gind = ((1:1953) -1) * NGENE + gidx
        subdf = indavgdf[gind,]
        gene = unique(subdf$symbol)
        # Add info:
        subdf$cr = paste0(subdf$celltype, '_',subdf$region)
        subdf$pcr = paste0(subdf$projid,'_',subdf$celltype, '_',subdf$region)
        subdf$pr = paste0(subdf$projid,'_',subdf$region)
        subdf$Est = rdf[subdf$cr, 'Est']
        subdf$count = ptdf[subdf$pcr, 'count']
        # Add attributes:
        nrmeta = unique(ctdf[,c('rind','projid','region','nrad','cogdxad','nft_ec','plaq_n_ec', 'Apoe_e4','msex','pmi','age_death')])
        rownames(nrmeta) = paste0(nrmeta$projid,"_",nrmeta$region)
        subdf$nrad = nrmeta[as.character(subdf$pr),'nrad']
        subdf$nft_ec = nrmeta[as.character(subdf$pr),'nft_ec'] # NOTE: EC NFT proxy for others
        subdf$cogdxad = nrmeta[as.character(subdf$pr),'cogdxad']
        subdf$pmi = nrmeta[as.character(subdf$pr),'pmi']
        subdf$age_rescaled = nrmeta[as.character(subdf$pr),'age_death'] / 100
        subdf$msex = nrmeta[as.character(subdf$pr),'msex']
        subdf$Apoe_e4 = nrmeta[as.character(subdf$pr),'Apoe_e4']
        subdf$rind = nrmeta[as.character(subdf$pr),'rind']
        subdf$nft = pqdf[as.character(subdf$rind),'nft']
        subdf$nft[is.na(subdf$nft)] = subdf$nft_ec[is.na(subdf$nft)] # Proxy for TH
        # subdf$plaq_n = pqdf[as.character(subdf$rind),'plaq_n']
        # Remove very small counts:
        subdf = subdf[subdf$count > 10,]

        # Regression - effect of expression on depletion:
        fit = glm(Est ~ val * nrad + Apoe_e4 + age_rescaled + msex + pmi, subdf, 
                  weights=log(subdf$count), family='gaussian') # Corrects for cell ct. but not inflated
                  # weights=sqrt(subdf$count), family='gaussian') # Corrects for cell ct., still oddly inflated
                  # weights=subdf$count, family='gaussian') # Way inflated
        # Regression - general effects of covariates on expression:
        fit2 = glm(val ~ nft * Est + celltype + Apoe_e4 + age_rescaled + msex + pmi, subdf, 
                  weights=log(subdf$count), family='gaussian') # Corrects for cell ct. but not inflated
        cfit = coefficients(summary(fit))
        cfit2 = coefficients(summary(fit2))
        df = data.frame(cfit)
        df2 = data.frame(cfit2)
        colnames(df) = c('Est','SE','t','p')
        colnames(df2) = c('Est','SE','t','p')
        pval = df['val','p']
        if (!is.na(pval)){
            if (pval < 1e-4){ cat(gene,'\t', sprintf('%0.2e',pval), '\n') }
        }
        df$var = rownames(df)
        df2$var = rownames(df2)
        rownames(df) = NULL
        rownames(df2) = NULL
        df$symbol = gene
        df2$symbol = gene
        est.regdf = rbind(est.regdf, df)
        est.expdf = rbind(est.expdf, df2)
    }
    save(est.regdf, est.expdf, file=est.regfile.rda)
} else { 
    load(est.regfile.rda)
}

est.regdf[est.regdf$symbol == 'OXR1',]
# genes = c('SGCD','SGCE','SGCZ','SGCG','DMD','ANO5', 'DYSF') # SGCE perturbed in other direction
# genes = c('WNT3A','DKK1','LRP6','CTNNB1','LRP5') # WNT assoc.
# est.regdf[est.regdf$symbol %in% genes & est.regdf$var == 'val',]

# Plot volcano of these assoc. with depletions:
adf = est.regdf[est.regdf$var == 'val',]
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
# cat(head(adf[adf$color == 1,'symbol'],50)) # GO enrichment turns up nothing of note
# adf[adf$symbol == 'HS6ST3',]

# adf[adf$symbol == 'MTERF1',]
# adf[adf$Est == min(adf$Est),]

tcols = brewer.pal(12, 'Paired')
gplot = ggplot(adf, aes(Est * val, log10q, col=factor(color))) + 
    scale_color_manual(values=c('grey85',tcols[1],tcols[5])) + 
    geom_vline(xintercept=0, lwd=.25, lty='dashed') + 
    geom_point(cex=.25, alpha=1) + theme_pubr() + 
    geom_text_repel(data=labdf, aes(Est * val, log10q, label=symbol), max.overlaps=30, size=2, segment.size=.5) + 
    scale_y_continuous(expand=c(0,0)) + 
    labs(x='Coefficient * Average Expression') + 
    theme(legend.position = 'none')
w = 7; h=7
ggsave(paste0(imgpref, 'ECHCTH_exc_subtypes_ind_pseudobulk_vuln_genes_volcano.png'), gplot, dpi=450, units='in', width=w, height=h)
ggsave(paste0(imgpref, 'ECHCTH_exc_subtypes_ind_pseudobulk_vuln_genes_volcano.pdf'), gplot, dpi=450, units='in', width=w, height=h)

# ---------------------------------------------------------------
# Plot the effect on expression from AD/non-AD for the top genes:
# ---------------------------------------------------------------
unique(est.expdf$var)
est.expdf[est.expdf$symbol == 'SYT1',c('var','p','Est')]

# vars = paste0("nradAD:celltype",topec[eid == 1])
vars = c('nft', 'Est','age_rescaled','msex','pmi','nft:Est')
edf = est.expdf[est.expdf$var == 'nft',]
# vars = 'age_rescaled'
# edf = est.expdf[est.expdf$var %in% vars,]
edf = merge(edf, avgidf)
edf = edf[!is.na(edf$p),]
edf = edf[edf$val > 1.5,]
edf = edf[order(edf$p), ]

head(edf[edf$Est > 0,], 20)

dgenes = head(adf$symbol,50)
dwdf = edf[edf$symbol %in% dgenes,]

head(dwdf[dwdf$Est > 0,], 20)

# Add in the nrad effect:
ndf = est.expdf[est.expdf$var == 'nradAD',c('Est','p','symbol')]
names(ndf)[1:2] = c('nrad.est', 'nrad.p')

adf = merge(adf, ndf)
adf$Est2 = adf$Est + adf$nrad.est

adf = adf[order(adf$Est2, decreasing=F), ]
head(adf[,c('symbol','Est2','p','nrad.p')], 50)

adf[adf$symbol == 'HS6ST3',]


adf$padj = p.adjust(adf$p, 'fdr')
adf$log10q = -log10(adf$padj)
pcut = 1e-3
adf$color = 0
adf$color[adf$padj < pcut] = 1
adf$color[adf$padj < pcut & adf$Est > 0] = 2

# -------------------------------------
# Plot some of these genes as boxplots:
# -------------------------------------
subdf = indavgdf[indavgdf$symbol %in% c('LRP1B','MT-ND3','NRXN3','HS6ST3','CDH20'),]
subdf = merge(subdf, rdf[,c('Est','region','celltype')])
subdf = merge(subdf, ptdf)
subdf = merge(subdf, mfracdf[,c('region','celltype','frac')])
nrmeta = unique(ctdf[,c('projid','nrad','cogdxad','nft_ec','plaq_n_ec', 'Apoe_e4','msex','pmi','age_death')])
rownames(nrmeta) = nrmeta$projid
subdf$nrad = nrmeta[as.character(subdf$projid),'nrad']
subdf$nft_ec = nrmeta[as.character(subdf$projid),'nft_ec']
subdf$cogdxad = nrmeta[as.character(subdf$projid),'cogdxad']
subdf$pmi = nrmeta[as.character(subdf$projid),'pmi']
subdf$age_rescaled = nrmeta[as.character(subdf$projid),'age_death'] / 100
subdf$msex = nrmeta[as.character(subdf$projid),'msex']
subdf$Apoe_e4 = nrmeta[as.character(subdf$projid),'Apoe_e4']
# Remove very small counts:
subdf = subdf[subdf$count > 10,]

gplot = ggplot(subdf[subdf$celltype %in% topec & subdf$region == 'EC',], aes(nft_ec, val, size=count)) +
    facet_grid(symbol~celltype, scale='free_y') +
    geom_point() + 
    geom_smooth(method='lm') + 
    scale_color_manual(values=colvals[['nrad']]) + 
    scale_y_continuous(expand=c(0,1))+
    theme_pubr() + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))


# -------------------------------------------------------
# Plot all the pseudobulk samples (+ split by AD/non-AD):
# -------------------------------------------------------
regset = 'nrad'
rdf = regdf[regdf$reg == regset,]
if (regset == 'nrad'){
    ecut = 0.0075
} else if (regset == 'cogdxad.rel'){
    ecut = 0.01
} else if (regset == 'cogdxad'){
    ecut = 0.01
} else { 
    ecut = 0.0005
}

mfracdf = aggregate(cbind(count, total) ~ region + celltype, ctdf, sum)
mfracdf$frac = mfracdf$count / mfracdf$total

ecdriven=FALSE
if (ecdriven){
    genes = c('HS6ST1','HS6ST2','HS6ST3', 'SULF1', 'SULF2', 'GPC5') # Heparan + EC
} else {
    genes = c('HS6ST3','TRPS1','CDH20', 'MAN1A1','SGCD','OXR1') # Discovered across all
    # genes = c('HS6ST3','MAN1A1','SPOCK3', 'FILIP1L') # Discovered across all
}
subdf = indavgdf[indavgdf$symbol %in% genes,]
subdf = merge(subdf, rdf[,c('Est','region','celltype')])
subdf = merge(subdf, ptdf)
subdf = merge(subdf, mfracdf[,c('region','celltype','frac')])
nrmeta = unique(ctdf[,c('projid','nrad','cogdxad','nft_ec','plaq_n_ec', 'Apoe_e4','msex','pmi','age_death')])
rownames(nrmeta) = nrmeta$projid
subdf$nrad = nrmeta[as.character(subdf$projid),'nrad']
subdf$nft_ec = nrmeta[as.character(subdf$projid),'nft_ec']
subdf$cogdxad = nrmeta[as.character(subdf$projid),'cogdxad']
subdf$pmi = nrmeta[as.character(subdf$projid),'pmi']
subdf$age_rescaled = nrmeta[as.character(subdf$projid),'age_death'] / 100
subdf$msex = nrmeta[as.character(subdf$projid),'msex']
subdf$Apoe_e4 = nrmeta[as.character(subdf$projid),'Apoe_e4']
# Remove very small counts:
subdf = subdf[subdf$count > 10,]


# subdf = subdf[subdf$frac > 0.005,] # Can include/ or not - removes noisy estimates
# subdf = subdf[subdf$major.celltype == 'Exc',]
# subdf = subdf[subdf$celltype != 'Exc SV2C LINC02137',]

# subdf = indavgdf[indavgdf$symbol %in% genes,]
# # Add info:

# Add p-values, fit lines:
# regs = c('PFC','AG','MT')
regs = c('EC','HC','TH')
# regs = names(reg.cols)
if (length(grep('rel', regset)) ==  0){ subdf$rel.Est = subdf$Est / subdf$frac } else { subdf$rel.Est = subdf$Est }

# lmdf = c()
# for (gene in genes){
#     sdf = subdf[subdf$symbol == gene & subdf$region %in% regs,]
#     fit = lm(Est ~ val,sdf)
#     cfit = coefficients(summary(fit))
#     df = data.frame(cfit['val',c(1,4), drop=F])
#     names(df) = c('lm.est','p')
#     df$intercept = cfit['(Intercept)',1]
#     df$symbol = gene
#     lmdf = rbind(lmdf,df)
# }
# lmdf$eqn = sprintf("y = %0.2f%% - %0.2f%% * Expression\np = %0.3f",
#                    round(lmdf$intercept * 100,2),
#                    round(lmdf$lm.est * -100,2), round(lmdf$p,3))

labdf = subdf[abs(subdf$Est) > ecut | subdf$val > 30,]
gplot = ggplot(subdf[subdf$region %in% regs,], aes(val, Est, color=region, size=count)) +
    facet_wrap(~symbol, scale='free_x') +
    geom_hline(yintercept=0,lwd=.25, lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point() + 
    # geom_text_repel(data=labdf[labdf$region %in% regs,], aes(val, Est, label=paste0(celltype, '; ', region)), max.overlaps=30, size=3) + 
    # geom_text(data=lmdf, aes(x=1,y=-.02, label=eqn), size=3, color='black') + 
    scale_color_manual(values=reg.cols) + 
    scale_y_continuous(label=scales::percent)+
    labs(x='Average Expression Level (Normalized)', y='Absolute change vs. all cells') +
    theme_pubr()



labdf = subdf[abs(subdf$Est) > ecut | subdf$val > 30,]
gplot = ggplot(subdf[subdf$region %in% regs,], aes(val, Est, color=region, size=count, fill=nrad)) +
    facet_wrap(~symbol, scale='free_x') +
    geom_hline(yintercept=0,lwd=.25, lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point(pch=22) + 
    # geom_text_repel(data=labdf[labdf$region %in% regs,], aes(val, Est, label=paste0(celltype, '; ', region)), max.overlaps=30, size=3) + 
    # geom_text(data=lmdf, aes(x=1,y=-.02, label=eqn), size=3, color='black') + 
    scale_color_manual(values=reg.cols) + 
    scale_y_continuous(label=scales::percent)+
    labs(x='Average Expression Level (Normalized)', y='Absolute change vs. all cells') +
    theme_pubr()


# Ask if these genes have difftl expr in AD/CTRL:
dedf = c()
for (gene in genes){
    sdf = subdf[subdf$symbol == gene,]
    fit = glm(val ~ nrad + Apoe_e4 + age_rescaled + msex + pmi + celltype, sdf, 
              weights=log(sdf$count), family='gaussian') # Corrects for cell ct. but not inflated
    cfit = data.frame(coefficients(summary(fit)))
    colnames(cfit) = c('Est','SE','t','p')
    cfit$symbol = gene
    cfit$var = rownames(cfit)
    rownames(cfit) = NULL
    dedf = rbind(dedf, cfit)
}

dedf[dedf$var == 'nradAD',]
dedf[dedf$var == 'age_rescaled',]



















# Plot these genes as well:
genes = labdf$symbol
subdf = avgdf[avgdf$symbol %in% genes,]
subdf = merge(subdf, rdf)
# subdf = subdf[subdf$celltype != 'Exc L2-3 CBLN2 LINC02306',] # Removing this completely = better signal tbh..
labdf = subdf[abs(subdf$Est) > ecut | subdf$val > 40,]

gplot = ggplot(subdf, aes(val, Est, color=region, size=ncell)) +
    facet_wrap(~symbol, scale='free_x') +
    geom_hline(yintercept=0,lwd=.25, lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point() + 
    geom_text_repel(data=labdf, aes(val, Est, label=paste0(celltype, '; ', region)), max.overlaps=30, size=3) + 
    scale_color_manual(values=reg.cols) + 
    theme_pubr()

ggsave(paste0(imgpref, 'fromall_allsubtypes_est_expr_scatter_comparison.png'), gplot, dpi=450, units='in', width=15, height=10)
ggsave(paste0(imgpref, 'fromall_allsubtypes_est_expr_scatter_comparison.pdf'), gplot, dpi=450, units='in', width=15, height=10)


# ----------------------------------------------------------
# Plot a score computed from the heparan sulfate regression:
# ----------------------------------------------------------
hsmat = avgwide[,c('celltype','region','ncell','major.celltype',hsgenes)]
hsmat$pred = predict(hs.fit, hsmat)
subdf = merge(hsmat[,c('celltype','region','ncell','major.celltype','pred')], rdf, all.x=TRUE)
subdf = merge(subdf, mfracdf[mfracdf$nrad == 'CTRL',])
subdf = subdf[subdf$frac > 0.005,]
# subdf = subdf[subdf$celltype != 'Exc L2-3 CBLN2 LINC02306',] # Removing this completely = better signal tbh..
subdf = subdf[!is.na(subdf$Est),]
labdf = subdf[abs(subdf$Est) > ecut | abs(subdf$pred) > ecut,]

gplot = ggplot(subdf, aes(pred, Est, color=region, size=ncell)) +
    facet_wrap(~major.celltype) + 
    geom_hline(yintercept=0,lwd=.25, lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point() + 
    geom_text_repel(data=labdf, aes(pred, Est, label=paste0(celltype, '; ', region)), max.overlaps=30, size=3) + 
    labs(x='Predicted absolute decrease', y='Actual absolute decrease') + 
    scale_color_manual(values=reg.cols) + 
    theme_pubr()

ggsave(paste0(imgpref, 'ECheparanscore_allsubtypes_est_expr_scatter_comparison.png'), gplot, dpi=450, units='in', width=12, height=7)
ggsave(paste0(imgpref, 'ECheparanscore_allsubtypes_est_expr_scatter_comparison.pdf'), gplot, dpi=450, units='in', width=12, height=7)


# --------------------------------------------------------------------------
# UMAP plots of the EC neurons only; colored by vulnerability and HS markers
# --------------------------------------------------------------------------
# TODO: Make separate UMAP for the EC neurons only?
submeta = cellmeta[cellmeta$cell_type_high_resolution %in% topec & cellmeta$region %in% 'EC',]
xr = range(submeta$U1)
yr = range(submeta$U2)
yr = c(-4,9)
xr = c(15,25)
dx = diff(xr)
dy = diff(yr)
dd = max(c(dx, dy)) / 2
xr = c(mean(xr) - dd, mean(xr) + dd)
yr = c(mean(yr) - dd, mean(yr) + dd)
submeta = cellmeta[cellmeta$U1 <= xr[2] & cellmeta$U1 >= xr[1],]
submeta = submeta[submeta$U2 <= yr[2] & submeta$U2 >= yr[1],]
submeta$celltype = as.character(submeta$cell_type_high_resolution)
submeta$celltype[!(submeta$celltype %in% topec)] = 'None'
submeta$celltype[submeta$region != 'EC'] = 'None'
submeta = merge(submeta, unique(metadata[c('projid','nrad')]))

# Shuffle points for plotting:
ind = 1:nrow(submeta)
ad.ind = which(submeta$nrad == 'AD')
ctrl.ind = which(submeta$nrad == 'CTRL')
full.ind = sample(ind, length(ind), replace=FALSE)
ad.ind = sample(ad.ind, length(ad.ind), replace=FALSE)
ctrl.ind = sample(ctrl.ind, length(ctrl.ind), replace=FALSE)

# Labels + colors:
tcols = c('None'='grey97',type.cols[topec])
labdf = aggregate(cbind(U1, U2) ~ cell_type_high_resolution, submeta, mean)
labdf$celltype = labdf$cell_type_high_resolution
labdf$celltype[!(labdf$celltype %in% topec)] = 'None'
labdf = labdf[labdf$celltype != 'None',]

w = 6
cex = 0.02
sp = 0.1
png(paste0(imgpref, 'EC_vuln_neurons_umap.png'), units='in', res=450, width=w, height=w)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(submeta$U1[full.ind], submeta$U2[full.ind], 
     col=tcols[submeta$celltype[full.ind]], pch=19, cex=cex, axes=F)
text(labdf$U1, labdf$U2, labdf$cell_type_high_resolution, font=2, cex=.8,
     col=ifelse(labdf$celltype == 'None', 'grey75','black'), xpd=TRUE)
mtext('UMAP 1', side=1, line=-1, cex=1.25)
mtext('UMAP 2', side=2, line=-1, cex=1.25)
dev.off()

png(paste0(imgpref, 'EC_vuln_neurons_umap_AD.png'), units='in', res=450, width=w, height=w)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(submeta$U1[ad.ind], submeta$U2[ad.ind], 
     col=tcols[submeta$celltype[ad.ind]],
     pch=19, cex=cex, axes=F)
text(labdf$U1, labdf$U2, labdf$cell_type_high_resolution, font=2, cex=.8,
     col=ifelse(labdf$celltype == 'None', 'grey75','black'), xpd=TRUE)
mtext('UMAP 1', side=1, line=-1, cex=1.25)
mtext('UMAP 2', side=2, line=-1, cex=1.25)
mtext('AD Individuals', side=3, line=-1.5, cex=1.5)
dev.off()

png(paste0(imgpref, 'EC_vuln_neurons_umap_CTRL.png'), units='in', res=450, width=w, height=w)
par(xaxs='i',yaxs='i')
par(mar=rep(sp,4))
plot(submeta$U1[ctrl.ind], submeta$U2[ctrl.ind], 
     col=tcols[submeta$celltype[ctrl.ind]],
     pch=19, cex=cex, axes=F)
text(labdf$U1, labdf$U2, labdf$cell_type_high_resolution, font=2, cex=.8,
     col=ifelse(labdf$celltype == 'None', 'grey75','black'), xpd=TRUE)
mtext('UMAP 1', side=1, line=-1, cex=1.25)
mtext('UMAP 2', side=2, line=-1, cex=1.25)
mtext('Non-AD Individuals', side=3, line=-1.5, cex=1.5)
dev.off()

# -------------------------------
# Plot UMAP with gene expression:
# -------------------------------
palette = c('grey90',viridis(100))
col_fun = function(x, pal=palette){
    bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
    palette[bin] 
}

# genes = c('HS6ST3','TRPS1','CDH20', 'MTERF1') # Discovered across all
# for (gene in genes){

gene = 'HS6ST3'
x = rep(0, nrow(submeta))
bind = which(submeta$barcode %in% colnames(nmat))
x[bind] = log(nmat[gene,submeta$barcode[bind]] + 1)

cex = 0.02
w = 6
png(paste0(imgpref, 'EC_vuln_neurons_umap_',gene,'.png'), units='in', res=450, width=w, height=w)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
bsp = 1.5
par(mar=c(sp,sp,sp, sp))
plot(submeta$U1[full.ind], submeta$U2[full.ind], 
     col=col_fun(x[full.ind]),
     pch=19, cex=cex, axes=F)
text(labdf$U1, labdf$U2, labdf$cell_type_high_resolution, font=2, cex=.8,
     col=ifelse(labdf$celltype == 'None', 'grey75','black'), xpd=TRUE)
mtext('UMAP 1', side=1, line=-1, cex=1.25)
mtext('UMAP 2', side=2, line=-1, cex=1.25)
mtext(gene, side=3, font=2, line=-1.5, cex=1.5)
dev.off()

# All four jointly:

cex = 0.01
w = 6
png(paste0(imgpref, 'EC_vuln_neurons_umap_',gene,'_2x2.png'), units='in', res=450, width=w, height=w)
par(xaxs='i', yaxs='i')
layout(matrix(1:4, nrow=2))
par(mar=rep(0,4))
plot(submeta$U1[full.ind], submeta$U2[full.ind], 
     col=tcols[submeta$celltype[full.ind]], pch=19, cex=cex, axes=F)
text(labdf$U1, labdf$U2, labdf$cell_type_high_resolution, font=2, cex=.6,
     col=ifelse(labdf$celltype == 'None', 'grey75','black'), xpd=TRUE)
mtext('UMAP 1', side=1, line=-1, cex=.8)
mtext('UMAP 2', side=2, line=-1, cex=.8)
mtext('All Individuals', side=3, line=-1.25, cex=1)
box(lwd=.5)
# Gene:
plot(submeta$U1[full.ind], submeta$U2[full.ind], 
     col=col_fun(x[full.ind]),
     pch=19, cex=cex, axes=F)
mtext(gene, side=3, line=-1.25, cex=1)
box(lwd=.5)
# AD:
plot(submeta$U1[ad.ind], submeta$U2[ad.ind], 
     col=tcols[submeta$celltype[ad.ind]],
     pch=19, cex=cex, axes=F)
mtext('AD Individuals', side=3, line=-1.25, cex=1)
box(lwd=.5)
# CTRL
plot(submeta$U1[ctrl.ind], submeta$U2[ctrl.ind], 
     col=tcols[submeta$celltype[ctrl.ind]],
     pch=19, cex=cex, axes=F)
mtext('Non-AD Individuals', side=3, line=-1.25, cex=1)
box(lwd=.5)
dev.off()


