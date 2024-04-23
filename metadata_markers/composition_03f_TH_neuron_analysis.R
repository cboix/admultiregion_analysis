#!/usr/bin/R
# --------------------------------------
# Investigate differences in TH neurons:
# Updated 04/26/2021 
# --------------------------------------
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
load('Annotation/multiregion_celltypes_colors.Rda')
tsp.tcols = sapply(tcols, tsp.col)

# Get the pathology mapped to each region:
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

metadata$cogdxad = 'CTRL'
metadata$cogdxad[metadata$cogdx %in% c(4,5)] = 'AD'
metadata$cogdxad = factor(metadata$cogdxad, levels=c('CTRL','AD'))
metadata$nrad = 'CTRL'
metadata$nrad[metadata$niareagansc %in% c(1,2)] = 'AD'
metadata$nrad = factor(metadata$nrad, levels=c('CTRL','AD'))


# ------------------------------------
# Load in average profiles of neurons:
# ------------------------------------
# TODO: Compare MEIS2 FOXP2 to other Inhibitory neurons
# TODO: Compare NXPH1 to other Excitatory neurons
# TODO: Compare each region to others.
indavg.file.rda = 'multiRegion/neuronal_exc_ECHCTH_ind_average_profiles.Rda'
load(indavg.file.rda)
# For filtering genes:
avgidf = aggregate(val ~ symbol, indavgdf, mean)

# For weighted regression:
ptdf = agg.rename(barcode ~ cell_type_high_resolution + projid + region, cellmeta, length, 'count')
names(ptdf)[1] = 'celltype'

# Faster indexing:
ptdf$pcr = paste0(ptdf$projid,'_',ptdf$celltype, '_',ptdf$region)
rownames(ptdf) = ptdf$pcr
rownames(pqdf) = pqdf$rind

# gene = 'HS6ST3'
ecut = 0.5
kept.genelist = avgidf[avgidf$val > ecut, 'symbol']
est.regfile.rda = 'multiRegion/gene_fraction_ECHCTH_byregion_association_table.Rda'
nrmeta = unique(metadata[,c('rind','projid','region','nrad','cogdxad', 'Apoe_e4','msex','pmi','age_death')])
rownames(nrmeta) = paste0(nrmeta$projid,"_",nrmeta$region)

# TODO: Repeat for the inhibitory neurons.
if (!file.exists(est.regfile.rda)){
    full.regdf = c()
    for (reg in c('EC','TH','HC')){
        # How many cells:
        est.regdf = c()
        NGENE = length(genelist) # Full genelist, for use in indexing.
        nrmeta$is.reg = nrmeta$region == reg
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
            subdf$count = ptdf[subdf$pcr, 'count']
            # Add attributes:
            subdf$nrad = nrmeta[as.character(subdf$pr),'nrad']
            subdf$nft_ec = nrmeta[as.character(subdf$pr),'nft_ec'] # NOTE: EC NFT proxy for others
            subdf$cogdxad = nrmeta[as.character(subdf$pr),'cogdxad']
            subdf$pmi = nrmeta[as.character(subdf$pr),'pmi']
            subdf$is.reg = nrmeta[as.character(subdf$pr),'is.reg']
            subdf$age_rescaled = nrmeta[as.character(subdf$pr),'age_death'] / 100
            subdf$msex = nrmeta[as.character(subdf$pr),'msex']
            subdf$Apoe_e4 = nrmeta[as.character(subdf$pr),'Apoe_e4']
            subdf$rind = nrmeta[as.character(subdf$pr),'rind']
            subdf$nft = pqdf[as.character(subdf$rind),'nft']
            # Remove very small counts:
            subdf = subdf[subdf$count > 100,]

            fit = glm(val ~ is.reg + nrad + Apoe_e4 + age_rescaled + msex + pmi, subdf, 
                      weights=log(subdf$count), family='gaussian')
            cfit = coefficients(summary(fit))
            df = data.frame(cfit)
            colnames(df) = c('Est','SE','t','p')
            pval = df['is.regTRUE','p']
            if (!is.na(pval)){
                if (pval < 1e-4){ cat(gene,'\t', sprintf('%0.2e',pval), '\n') }
            }
            df$var = rownames(df)
            rownames(df) = NULL
            df$symbol = gene
            est.regdf = rbind(est.regdf, df)
        }
        est.regdf$region = reg
        full.regdf = rbind(full.regdf, est.regdf)
    }
    save(full.regdf, file=est.regfile.rda)
} else { 
    load(est.regfile.rda)
}

# Plot volcano of these assoc. with depletions:
adf = est.regdf[est.regdf$var == 'is.regTRUE',]
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
ggsave(paste0(imgpref, 'ECHCTH_byregion_differences_ind_pseudobulk_vuln_genes_volcano.png'), gplot, dpi=450, units='in', width=w, height=h)
ggsave(paste0(imgpref, 'ECHCTH_byregion_differences_ind_pseudobulk_vuln_genes_volcano.pdf'), gplot, dpi=450, units='in', width=w, height=h)


# -----------------------------------------------------------------
# Plot heatmap of top genes for each region (adapt from following):
# -----------------------------------------------------------------
adf = full.regdf[full.regdf$var == 'is.regTRUE',]
adf = unique(adf[order(adf$p),])
adf = adf[adf$Est > 0,]
pgenes = c()
pset = c()
ntop = 8
for (reg in c('EC','TH','HC')){
    sdf = adf[adf$region == reg,]
    pgenes = c(pgenes, head(unique(sdf$symbol),ntop))
    pset = c(pset, rep(reg, ntop))
}

# Matrix:
umeta = indavgwide[,c('projid','region','celltype')]
umeta$pcr = with(umeta, paste0(projid,'_',celltype, '_',region))
umeta$count = ptdf[umeta$pcr,'count']
ind = umeta$count > 100
umeta = umeta[ind,]
pmat = t(as.matrix(indavgwide[ind, pgenes]))

# Make matrix and annotations:
# kmeta = umeta[umeta$count > 100,]
# plt.mat = pmat[pgenes,kmeta$ptype]
clsplit = umeta$region
ha = HeatmapAnnotation(CT=umeta$celltype, 
                       Region=umeta$region,
                       col=list(Region=reg.cols,
                                CT=tcols))
udsplit = pset
hb = rowAnnotation(Set=pset, col=list(Set=reg.cols))

smat = as.matrix(log(pmat + 1))
smat.scaled = t(scale(t(smat), center=FALSE))

# png(paste0(imgpref, 'ECHCTH_EXC_topDE_heatmap_normalized_individ.png'), res=400, units='in', width=11, height=4)
pdf(paste0(imgpref, 'ECHCTH_EXC_topDE_heatmap_normalized_individ.pdf'), width=11, height=4)
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

# GABA:
# adf = full.regdf[full.regdf$var == 'is.regTRUE',]
# adf[grep("^GABR",adf$symbol),]


# GO ENR:
adf = full.regdf[full.regdf$var == 'is.regTRUE',]
adf = unique(adf[order(adf$p),])
adf = adf[adf$Est > 0,]
pgenes = c()
pset = c()
ntop = 50
for (reg in c('EC','TH','HC')){
    sdf = adf[adf$region == reg,]
    pgenes = head(unique(sdf$symbol), ntop)
}



