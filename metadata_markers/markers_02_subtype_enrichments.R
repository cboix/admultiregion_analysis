#!/usr/bin/R
# -------------------------------------------------------------------
# Compute pathway enrichments for differences subtypes in a celltype:
# Updated 12/02/21
# -------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(gprofiler2)

library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)
library(patchwork)

library(ComplexHeatmap)
library(circlize)

options(width=175)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')  # TODO: Change dir for these analyses
plotdir = paste0(imgdir, 'metadata/')  # TODO: Change dir for these analyses
imgpref = paste0(plotdir, 'enr_')
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
remove.batches = TRUE
suff = '_subset_final_noMB'

subset = 'Inh'

ststr = gsub("/","_", subset)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
subtypes = unique(ctdf$cls)
celltypes = unique(ctdf$major.celltype)

if (subset == 'Inh'){
    subtypes = c('Thalamus','PAX6','LAMP5',
                 'VIP', 'SST','PVALB')
    subtypecols = snap.cols[j:(j + length(subtypes) - 1)]
    names(subtypecols) = subtypes
} else if (subset == 'Exc'){
} else { subtypecols = c(tcols[subtypes]) }



# Load the scores for differential genes:
# ---------------------------------------
full.regdf = c()
for (st in subtypes){
    print(st)
    sub.ststr = gsub("/","_",gsub(" ","_", st))
    subtype.diff.rda = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_', sub.ststr, '.rda')
    load(subtype.diff.rda)
    full.regdf = rbind(full.regdf, est.regdf)
}

# Subset the regression results to up + is.subtype:
updf = full.regdf[full.regdf$var == 'is.stTRUE',]
updf = updf[order(updf$p),]
updf = updf[updf$Est > 0,]
updf$p.adj = p.adjust(updf$p, 'fdr')



# Get the differentially up genesets for each subtype:
# ----------------------------------------------------
pcutoff = 0.01
genesets = lapply(subtypes, function(x){
                      updf[(updf$st == x) & (updf$p.adj < pcutoff), 'symbol']
})
print(sapply(genesets, length))
names(genesets) = sapply(subtypes, function(x){sub(".* ","",x)})
geneset.map = subtypes
names(geneset.map) = names(genesets)


# Run the pathway enrichments for all subtypes together:
# ------------------------------------------------------
overall.enr.rda = paste0(srdir, 'enrichments_difftl_pseudobulk_markers_', ststr, '.rda')
if (!file.exists(overall.enr.rda)){
    sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
    gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
                                  ordered_query=FALSE, multi_query=TRUE,
                                  sources = sources)
    gp2df = gp2.result$result

    # Process table:
    pmat = t(as.matrix(data.frame(gp2df$p_values)))
    rownames(pmat) = NULL  # Format 
    pvals = c()
    for (j in 1:length(genesets)){
        pstr = paste0('p_', names(genesets)[j])
        gp2df[[pstr]] = pmat[,j]
        pvals = c(pvals, pstr)
    }
    gp2cols = c('term_id',pvals,'source','term_name', 'term_size')
    gp2df = gp2df[,gp2cols]
    save(gp2.result, gp2df, genesets, file=overall.enr.rda)
} else { load(overall.enr.rda) }


# For each subtype, get the unique (or only 2 subtypes):
# ------------------------------------------------------
gp2df$nc = nchar(gp2df$term_name)
subdf = gp2df[gp2df$term_size < 500,]
pvals = colnames(subdf)[grep("^p_", colnames(subdf))]
pmat = subdf[, pvals]
pcutoff = 0.05
max.subtypes = 2
uqterm = apply(pmat < pcutoff, 1, sum)
topdf = subdf[uqterm <= max.subtypes,]


# Aggregate a list of top terms for a heatmap:
# --------------------------------------------
ntop = 10
keep.terms = c()
pset = c()
for (pval in pvals){
    subtopdf = topdf[topdf[[pval]] < pcutoff,]
    subtopdf = subtopdf[subtopdf$nc <= 40,]
    subtopdf = subtopdf[order(subtopdf[[pval]]),]
    keep.terms = c(keep.terms, head(subtopdf$term_id, ntop))
    pset = c(pset, rep(sub("p_", "", pval), ntop))
}


# Plot heatmap of top terms for each type:
# ----------------------------------------
rownames(gp2df) = gp2df$term_id
seldf = gp2df[keep.terms,]
pmat = as.matrix(seldf[,pvals])
pmat = -log10(pmat)
rownames(pmat) = seldf$term_name

hb = rowAnnotation(Set=geneset.map[pset],
                   col=list(Set=subtypecols))
# Cap p-values on top + bottom:
pmat[pmat > 5] = 5
pmat[pmat < -log10(0.05)] = 0

ux = 1.5
plt = Heatmap(pmat, name='-log10p',
              col=c('white', colb),
              use_raster=FALSE,
              width=ncol(pmat) * unit(ux, 'mm'),
              height=nrow(pmat) * unit(ux, 'mm'),
              row_split=pset, 
              right_annotation=hb,
              border_gp=gpar(color='black', lwd=.5))

pltprefix = paste0(imgpref, ststr, '_topDE_subtype_enrichments')
saveHeatmap(plt, pltprefix, w=5, h=6)

