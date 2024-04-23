#!/usr/bin/R
# -----------------------
# Nebula minimal example:
# -----------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)

# For nebula: 
library(nebula)
library(DESeq2)
library(RUVSeq)
library(qvalue)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)

# For comparison:
library(lme4)
print(version)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# -------------------------
# Load in these data files:
# -------------------------
prefstr = '_test_deg_Mic_Immune__EC'
matfile = paste0(datadir, 'matrix', prefstr, '.tsv.gz')
varfile = paste0(datadir, 'genes', prefstr, '.tsv')
obsfile = paste0(datadir, 'metadata', prefstr, '.tsv')

mat = read.delim(matfile, header=F)
mat = as.matrix(mat)
obs = read.delim(obsfile, header=T)
var = read.delim(varfile, header=T)
rownames(mat) = obs$obsnames
colnames(mat) = var$varnames

# Normalize internally (approx by colSums)
norm = sweep(mat, 2, colSums(mat) / 10000,'/')
norm = log(norm + 1)

# Split on "bynft":
adind = obs$bynft == 'True'
ctind = obs$bynft == 'False'

# pc_genes = anno$symbol[anno$type == 'protein_coding']
# genes = genes[genes %in% pc_genes]

# ----------------------------------------
# Wilcoxon tests on the normalized matrix:
# ----------------------------------------
# TODO: Try presto for faster runs?
genes = colnames(mat)
wtdf = c()
for (gene in genes){
    x1 = norm[ctind, gene]
    x2 = norm[adind, gene]
    wt = wilcox.test(x1, x2)
    wtdf = rbind(wtdf, 
          data.frame(gene=gene,
                     lp=-log10(wt$p.value),
                     ratio=mean(x2) / mean(x1)))
}
wtdf = wtdf[order(wtdf$lp, decreasing=T),]

# --------------------------
# GLMM with log-normal data:
# --------------------------
genes = colnames(mat)
glmdf = c()
subdf = obs[,c('projid','bynft','nft','nGene','cpg','pctMT')]
for (gene in genes){
    print(gene)
    subdf$x = norm[,gene]
    fit1 = lmer(x ~ nGene + cpg + (1|projid), subdf, REML=F)
    fit2 = lmer(x ~ nft + nGene + cpg + (1|projid), subdf, REML=F)
    av = anova(fit1, fit2)
    p = av['Pr(>Chisq)'][2,]
    chisq = av['Chisq'][2,]
    cfit = data.frame(coefficients(summary(fit2)))
    glmdf = rbind(glmdf,
                  data.frame(gene=gene, lp=-log10(p), chisq=chisq,
                             est=cfit['nft','Estimate']))
}
glmdf = glmdf[order(glmdf$lp, decreasing=T),]

glmdf$p = 10^(-glmdf$lp)
glmdf$padj = p.adjust(glmdf$p, 'fdr')

head(glmdf[glmdf$est > 0,], 40)


# --------------------------------------------
# Run RUV on the data at the individual level:
# --------------------------------------------

# Make the individual-aggregate matrix:
pids = as.character(unique(obs$projid))
tform = make.tform(obs$projid, u=pids)
data_ind = t(mat) %*% tform 

# Make the aggregate design matrix:
uqobs = unique(obs[,c('projid','nft')])
rownames(uqobs) = uqobs$projid
uqobs = uqobs[colnames(data_ind),]
design = model.matrix(~ nft, data=uqobs)

# DESeq2 object
d_e <- DGEList(data_ind, genes=rownames(data_ind))
keep <- rowSums(cpm(d_e)>1) >= 3
d_e <- d_e[keep, , keep.lib.sizes=FALSE]
d_e <- calcNormFactors(d_e, method="TMM")
d_e <- estimateGLMCommonDisp(d_e, design)
d_e <- estimateGLMTagwiseDisp(d_e, design)
fit1 <- glmFit(d_e, design)
res1 <- residuals(fit1, type="deviance")
ruvn <- 10
ruv_cov <- RUVr(round(d_e$counts), 
                as.character(rownames(d_e$counts)), 
                k=ruvn, res1)

# Merge the learned factors back into the data.frame:
uqobs = cbind(uqobs, ruv_cov$W)
obs = merge(obs, uqobs, all.x=TRUE)
rownames(obs) = obs$obsnames
obs = obs[rownames(mat),]

# ------------------------------------------
# Test nebula results for differential expr:
# ------------------------------------------
nc = rowSums(mat)
ng = rowSums(mat > 0)
cpg = nc / ng
nmt = rowSums(mat[,c(grep("^MT-", genes), grep("^MTRNR", genes))])
pct_mt = nmt / nc
obs$nGene = ng
obs$cpg = cpg
obs$pctMT = pct_mt

# Model of "bynft" has decent concordance with GLMM on log1p(norm):
# mdx = model.matrix(~ bynft + cpg + nGene, data=obs) 
# Alternatively, nft + RUV terms:
mdx = model.matrix(~ nft + cpg + nGene + celltype +
                   W_1 + W_2 + W_3 + W_4 + W_5 + 
                   W_6 + W_7 + W_8 + W_9 + W_10, data=obs)

pathstr = 'nft'
lint = paste0('logFC_(Intercept)')
leff = paste0('logFC_', pathstr)
peff = paste0('p_', pathstr)

# For if we test subsampling or not:
oind = 1:nrow(obs)
# scale = 2
# oind = sort(sample(oind, round(length(oind)/scale)))

chunksize=500
fulldf = c()
for (chunk in 1:5){
    print(chunk)
    ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, ncol(mat))) 
    submat = t(mat[oind,ind])
    re = nebula(submat, as.character(obs$projid)[oind], pred=mdx[oind,],offset=log10(nc)[oind], model='PMM') 
    rdf = re$summary
    resdf = rdf[order(rdf[[peff]]),c(lint,leff, peff, 'gene')]
    names(resdf) = c('logFC_int','logFC','p','gene')
    fulldf = rbind(fulldf, resdf)
}
fulldf = fulldf[order(fulldf$p),]
fulldf$padj = p.adjust(fulldf$p, 'fdr')
fulldf$q = qvalue(fulldf$p)$q

pcut = 0.01
fulldf$col = 1 * (fulldf$q < pcut) * (2 - 1 * (fulldf$logFC < 0))
labdf = fulldf[fulldf$col != 0,]

pcols = brewer.pal(12,'Paired')
gplot = ggplot(fulldf, aes(logFC, -log10(p), color=factor(col))) + 
    geom_point(cex=.25) + 
    geom_text_repel(data=labdf, aes(logFC, -log10(p), label=gene, color=factor(col)), size=2, max.overlaps=20) + 
    scale_color_manual(values=c('grey80',pcols[1],pcols[5])) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'volcano_',prefstr,'_nebula_minimal.png'),gplot, units='in', dpi=450, width=6, height=6)
ggsave(paste0(imgpref, 'volcano_',prefstr,'_nebula_minimal.pdf'),gplot, units='in', dpi=450, width=6, height=6)

head(labdf[order(labdf$logFC,decreasing=T),], 40)

# TODO: what happens if we also model the sub-celltype

# ---------------------------------------
# Plot comparison between different runs:
# ---------------------------------------
# GLMM comparison:
gfdf = merge(glmdf[,c('gene','est','lp')], fulldf[,c('gene','logFC','p')])
gfdf$col = 1 * (gfdf$p < 0.05) + 2 * (gfdf$lp > -log10(0.05))
labdf = gfdf[gfdf$col != 0,]

gplot = ggplot(gfdf, aes(est, logFC, color=factor(col))) + 
    geom_point() + 
    geom_smooth(color='black',method='lm') + 
    geom_text_repel(data=gfdf, aes(est, logFC, label=gene, color=factor(col))) + 
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=0) + 
    scale_color_manual(values=c('grey85','indianred','slateblue','goldenrod3')) +
    labs(x='Estimate from GLMM', y='logFC from PMM (Nebula)') + 
    theme_pubr()
ggsave(paste0(imgpref, 'glmcomp_',prefstr,'_nebula_minimal.png'),gplot, units='in', dpi=450, width=6, height=7)

# Wilcoxon comparison:
wfdf = merge(wtdf[,c('gene','ratio','lp')], fulldf[,c('gene','logFC','p')])
wfdf$col = 1 * (wfdf$p < 0.05) + 2 * (wfdf$lp > -log10(0.05))
labdf = wfdf[wfdf$col != 0,]

gplot = ggplot(wfdf, aes(log2(ratio), logFC, color=factor(col))) + 
    geom_point() + 
    geom_smooth(color='black',method='lm') + 
    geom_text_repel(data=wfdf, aes(log2(ratio), logFC, label=gene, color=factor(col))) + 
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=0) + 
    scale_color_manual(values=c('grey85','indianred','slateblue','goldenrod3')) +
    labs(x='logFC from Wilcoxon', y='logFC from PMM (Nebula)') + 
    theme_pubr()
ggsave(paste0(imgpref, 'wxcomp_',prefstr,'_nebula_minimal.png'),gplot, units='in', dpi=450, width=6, height=7)






