#!/usr/bin/R
# ------------------------------------------
# Print the GWAS genes in the modules / DEGs
# Updated 12/16/2021
# ------------------------------------------
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
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'gwas_')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)


# Set the run arguments:
# ----------------------
# TODO: Load across all?
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Load in the GWAS locus data:
# ----------------------------
anndir = paste0(dbdir, 'Annotation/')
gwdf = read.delim(paste0(anndir, '20210915_ADGENES_CHROM_Tanzi.tsv'), header=T)
gwgenes = unique(gwdf$gene[gwdf$evidence == 'GWAS'])


# Overall statistics for GWAS genes vs. modules:
# ----------------------------------------------
# Number of GWAS genes in core modules vs. total
incore = gwgenes[gwgenes %in% names(coremap)]
infull = gwgenes[gwgenes %in% names(genemap)]

# Module enrichment for GWAS genes:
coregw = coremap[incore]
cdf = data.frame(table(coremap))
names(cdf) = c('module','nmod')
cgwdf = data.frame(table(coregw))
names(cgwdf) = c('module','ngw')

cdf = merge(cgwdf, cdf, all.y=TRUE)
cdf$ngw[is.na(cdf$ngw)] = 0
cdf$ntot = length(coremap)
cdf$ngwtot = length(incore)
cdf$p = apply(cdf[,c('ngw','nmod','ngwtot','ntot')], 1, run.hyper)
cdf = cdf[order(cdf$p),]
cdf$genes = sapply(cdf$module, function(i){paste(sort(names(coregw)[coregw == i]), collapse=',')})

head(cdf, 20)


# Look for genes in DE table:
# ---------------------------
path = 'cogdxad'
gwdedf = dedf[(dedf$gene %in% gwgenes) & (dedf$dkey == path),]
gwdedf$gene[gwdedf$gset == 'Up']
gwdedf$gene[gwdedf$gset == 'Down']
gwdedf$gene[gwdedf$gset == '--']


# Plot GWAS genes DE tested as heatmap:
# -------------------------------------
gwdedf = dedf[(dedf$gene %in% gwgenes),]
gwdedf$lp = ifelse(gwdedf$col_nm > 0, gwdedf$log10p, 0)
cmat = pivot.tomatrix(gwdedf[,c('gene','dkey','logFC_nb')], 'dkey', 'logFC_nb')
pmat = pivot.tomatrix(gwdedf[,c('gene','dkey','lp')], 'dkey', 'lp')
cmat[,c('nft','plaq_n','plaq_d')] = cmat[,c('nft','plaq_n','plaq_d')] * 10
pmat = 10**(-pmat)

col_fun = colorRamp2(c(-.25, 0, .25), c('blue', "white", 'red'))

plt = Heatmap(cmat,
              col=col_fun,
              use_raster=TRUE,
              # column_split=ll$column_split,
              cluster_columns=TRUE,
              cluster_rows=TRUE,
              width = ncol(cmat)*unit(8, "mm"), 
              height = nrow(cmat)*unit(4, "mm"),
              border_gp = gpar(col="black", lty = 1),
              cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                  p = pmat[i,j]
                  ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                  grid.text(ann, x, y)}
)

h = 2.25 + 2.5 / 15 * nrow(cmat)
w = 5 + 2.5 / 15 * ncol(cmat)
pltprefix = paste0(imgpref, 'gwas_desimple_heatmap_', fullpref)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()




