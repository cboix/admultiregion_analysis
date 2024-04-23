#!/usr/bin/R
# ---------------------------------------------------------
# Plot functional and category enrichments for inset plots:
# Updated 12/09/2022
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ggpubr)
library(ggplot2)
library(ggrepel)

library(ComplexHeatmap)
library(circlize)
options(width=170)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'modules/auxiliary_gprofiler_functions.R'))

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/pltres/')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}

imgpref = paste0(plotdir, 'modules_indpt_res_', runset, '_')


# Load in and process data and UMAP coordinates:
# ----------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


commandArgs <- function(trailingOnly=TRUE){c(runset)}
source(paste0(sbindir, 'modules/load_modules_umap_coords.R'))

# Select modules with at least 10 genes:
ctdf = aggregate(gene ~ leiden, nodedf, length)
plt.modules = ctdf$leiden[ctdf$gene >= 10]


# Load in set of top by p-value functional enrichments for each module:
# ---------------------------------------------------------------------
useset = 'coregenes_'
# useset = 'allgenes_'
# useset = 'deonly_'
moduleenr.file = paste0(moddir, 'module_enrichments_',
    useset, fullpref, '.rda')
load(moduleenr.file)

sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
sub.sources = c("REAC","WP","KEGG","CORUM")


# Plots for functional enrichments:
mset = c(17,19)
setname = '2panel'
mset = c(1,2)
setname = 'APOEpanel'
mset = c(0,2)
setname = 'Ast-APOE'
# mset = c(12, 24, 7, 19, 9, 0, 17, 3, 8)
# setname = '8panel'
rownames(mmap) = as.character(mmap$module)
mns = mmap[as.character(mset), 'mname']

genes = names(coremap)[coremap == mset[1]]
cat(genes)


# Heatmap:
NTOP = 20
gpdf = gp2.result$result
gpdf = gpdf[gpdf$source %in% sources,]
gpdf = gpdf[gpdf$term_size < 500,]
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP, keep.sets=mns)
subpmat[subpmat > 5] = 5  # For plotting, cut off very small p-values
# subpmat = t(diag.mat2(t(subpmat))[[1]]

# Plot this enrichment matrix:
pltmat = subpmat
gap = 0; ux = 1.5;
plt = Heatmap(pltmat,
              cluster_columns=FALSE,
              cluster_rows=FALSE,
              use_raster=FALSE,
              border_gp=gpar(color='black'),
              name='-log10(p)',
              row_gap=unit(gap, "mm"),
              column_gap=unit(gap, "mm"),
              width = ncol(pltmat)*unit(ux, "mm"), 
              height = nrow(pltmat)*unit(ux, "mm"),
              col=c('white',colb))

h = 2.25 + 1 / 15 * nrow(pltmat)
w = 5 + 1 / 15 * ncol(pltmat)
pltprefix = paste0(imgpref, 'funcenr_n', NTOP, '_', graphpref, '_', setname,  '_inset')
saveHeatmap(plt, pltprefix, w=w, h=h)


# Plot barplots as well:
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP, keep.sets=mns, thresh=NULL)
sdf = data.frame(subpmat, check.names=F)
sdf$term = rownames(sdf)
sdf = gather(sdf, mname, log10p, -term)
sdf = sdf[sdf$log10p > -log10(0.05),]
sdf$term = factor(sdf$term, levels=rev(unique(sdf$term)))

coldf = unique(nodedf[,c('leiden','col')])
names(coldf)[1] = c('module')
coldf = merge(coldf, mmap)
mcolmap = coldf$col
names(mcolmap) = coldf$mname

gp = ggplot(sdf, aes(term, log10p, fill=mname)) + 
    facet_wrap(~mname, scale='free', ncol=1) + 
    scale_fill_manual(values=mcolmap) + 
    geom_bar(stat='identity') + coord_flip() + 
    theme_pubr() + theme(legend.position='none')
pltprefix = paste0(imgpref, 'funcenr_n', NTOP, '_', graphpref, '_', setname,  '_bar_inset')
saveGGplot(gp, pltprefix, w=8, h=2 * nrow(sdf) / 8 + 1)


# Load in the metadata enrichments for each module:
# -------------------------------------------------
hgann.file = paste0(moddir, 'module_covariate_hgdf_', useset, fullpref, '.tsv')
hgdf = read.delim(hgann.file, sep="\t")
hgdf$cl = paste0(hgdf$covariate, ' ', hgdf$level)
hgdf = merge(hgdf, aggregate(log10p ~ cl + mname, hgdf, max)) # Take max of cls 1/0 pval
# Calculate log2FC in each scenario:
hgdf$lfc[hgdf$cls == 1] = with(hgdf[hgdf$cls == 1,], log2((q/draw)/(m/N)))
hgdf$lfc[hgdf$cls == 0] = with(hgdf[hgdf$cls == 0,], log2(((m-q)/(N-draw))/(m/N)))

# Depletion in case of 0:
# Number with covar in total / number w/out vs. pop ratio
# (N- q) (N - draw)
# q = number with covariate level and z-score > 1 (cls=1) or z-score < 1 (cls=0)
# draws = number with z-score > 1 (cls=1) or z < 1 (cls=0)
# m = number with covariate level
# N = nrow of module (psbulk)

mset = c(17,19)
setname = '2panel'
mset = c(12, 24, 7, 19, 9, 0, 17, 3, 8)
setname = '8panel'
mset = c(0,1)
setname = 'APOEpanel'
rownames(mmap) = as.character(mmap$module)
mns = mmap[as.character(mset), 'mname']


# Plot the enrichment/depletions of covariates:
# ---------------------------------------------
# Reduce to selected names and significant covariates:
sdf = hgdf[hgdf$mname %in% mns,]
rm.covars = c('ceradsc','braaksc','niareagansc', 'apoe_genotype')
sdf = sdf[!(sdf$covariate %in% rm.covars),]
mdf = aggregate(log10p ~ covariate, sdf, max)
sdf = merge(sdf, mdf[mdf$log10p > -log10(0.01),'covariate', drop=F])
sdf$sig = as.character(sdf$log10p > 2)

gp = ggplot(sdf, aes(cl, lfc, fill=mname, alpha=sig)) + 
    facet_wrap(.~ mname, nrow=1) + 
    geom_hline(yintercept=0, lty='dashed', lwd=.5) +
    scale_fill_manual(values=mcolmap) + 
    scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.25)) + 
    # geom_segment(data=sdf, aes(y=0, yend=lfc, x=level, xend=level), lty='dotted', lwd=.5, alpha=1, color='grey') + 
    geom_bar(stat='identity') + coord_flip() +
    # geom_point(pch=21, cex=3) + 
    theme_pubr() + theme(legend.position='none')

pltprefix = paste0(imgpref, 'covenr_', graphpref, '_', setname,  '_bar_inset')
saveGGplot(gp, pltprefix, w=3 + 1 * length(mset), h=2 * length(unique(sdf$cl)) / 10 + .5)


# Plot all braaksc and niareagansc 
# ---------------------------------
setname = 'allmodBraak'
rownames(mmap) = as.character(mmap$module)
mns = mmap[as.character(mset), 'mname']
sdf = hgdf[hgdf$mname %in% mns,]
rm.covars = c('ceradsc','braaksc','niareagansc', 'apoe_genotype')
sdf = sdf[(sdf$covariate %in% rm.covars),]
mdf = aggregate(log10p ~ covariate, sdf, max)
sdf = merge(sdf, mdf[mdf$log10p > -log10(0.01),'covariate', drop=F])
sdf$sig = as.character(sdf$log10p > 2)
mset = plt.modules


sdf$mn = sub(" \\(.*","", sdf$mname)
gp = ggplot(sdf, aes(cl, lfc, fill=mname, alpha=sig)) + 
    facet_wrap(.~ mn, nrow=1) + 
    geom_hline(yintercept=0, lty='dashed', lwd=.5) +
    scale_fill_manual(values=mcolmap) + 
    scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.25)) + 
    # geom_segment(data=sdf, aes(y=0, yend=lfc, x=level, xend=level), lty='dotted', lwd=.5, alpha=1, color='grey') + 
    geom_bar(stat='identity') + coord_flip() +
    # geom_point(pch=21, cex=3) + 
    theme_pubr() + theme(legend.position='none')

pltprefix = paste0(imgpref, 'covenr_', graphpref, '_', setname,  '_bar_inset')
saveGGplot(gp, pltprefix, w=3 + 1 * length(mset), h=2 * length(unique(sdf$cl)) / 10 + .5)


# Plot as heatmap:
sdf$p = 1 - 1 * (sdf$sig == 'TRUE')
cmat = pivot.tomatrix(sdf[,c('mname','cl','lfc'),], 'mname','lfc')
pmat = pivot.tomatrix(sdf[,c('mname','cl','p'),], 'mname','p')
cmat[is.infinite(cmat)] = 0

rsplit = sub(" .*","", rownames(cmat))
mx = max(abs(cmat), na.rm=T)
col_fun = colorRamp2(c(-mx, 0, mx), c('blue', "white", "red")) 

ux = 1.5
ht = Heatmap(cmat,
    use_raster=FALSE,
    name='logFC',
    col=col_fun,
    row_split=rsplit,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(cmat)*unit(ux, "mm"), 
    height = nrow(cmat)*unit(ux, "mm"),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        p = pmat[i,j]
        ann = ifelse(p < 0.05, '*', '')
        grid.text(ann, x, y, gp=gpar(fontsize=gridtxt.fs))}
)

h = 2 + 1 / 15 * nrow(cmat)
w = 2 + 1 / 15 * ncol(cmat)
pltprefix = paste0(imgpref, 'covenr_', graphpref, '_', setname, '_heatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)

