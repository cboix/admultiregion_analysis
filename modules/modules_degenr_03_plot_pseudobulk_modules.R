#!/usr/bin/R
# ----------------------------------------------
# Plot the module scores on the pseudobulk data:
# Updated 11/27/2021
# ----------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ggpubr)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
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


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Load in the full pseudobulk data for these subtypes:
# ----------------------------------------------------
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
if (!file.exists(psdata.rda)){
    ps.data = load_pseudobulk_dataset(celltype, subtypes, region.set)
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }

umeta = ps.data$meta
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
umeta = unique(umeta)
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(ps.data$mat),]
# Remove very low abundance batches + will use for weight
umeta = umeta[umeta$ncell > 10,] 
umeta = unique(umeta)


# Score all modules for (a) all genes and (b) tested DE genes:
# ------------------------------------------------------------
core.only = TRUE
on.degenes = FALSE

# Score coregenes only:
if (core.only){
    useset = 'coregenes_'
    usemap = coremap
} else { 
    useset = ''
    usemap = genemap
}

# For all genes:
modules = sort(unique(usemap))
kept.genes = rownames(ps.data$mat)
kept.genes = kept.genes[kept.genes %in% names(usemap)]

# For DE tested only:
if (on.degenes){ 
    tested.degenes = unique(dedf$gene)
    kept.genes = kept.genes[kept.genes %in% tested.degenes]
    imgsuff = paste0(fullpref, '_degs')
} else {
    imgsuff = fullpref
}

# Conversion matrix:
tform = make.tform(usemap[kept.genes], u=modules, norm=TRUE)
mod.mat = t(tform) %*% ps.data$mat[kept.genes,umeta$ptype]
mod.mat = as.matrix(mod.mat)
rownames(mod.mat) = mmap$mname[1:nrow(mod.mat)]

# Remove non-scored (if using DE genes):
kept.mod = which(!is.na(rowSums(mod.mat)))
mod.mat = mod.mat[kept.mod,]

full.projids = sort(unique(cellmeta$projid))
projid.cols = snap.cols[1:48]
names(projid.cols) = full.projids

# Make annotation from pseudobulk data:
clsplit = umeta$nrad
nft.col_fun = colorRamp2(range(umeta$nft), c("white", "indianred"))
ha = HeatmapAnnotation(CT=umeta$cell_type_high_resolution, 
                       Region=umeta$region,
                       Braak=umeta$braaksc,
                       NFT=umeta$nft,
                       nrad=umeta$nrad,
                       cogdxad=umeta$cogdxad,
                       ncell=umeta$ncell,
                       e4=umeta$Apoe_e4,
                       projid=as.character(umeta$projid),
                       col=list(CT=tcols[plt.subtypes],
                                Region=reg.cols,
                                Braak=colvals[['braaksc']],
                                nrad=colvals[['nrad']],
                                NFT=nft.col_fun,
                                cogdxad=colvals[['cogdxad']],
                                projid=projid.cols,
                                e4=c('no'='grey90','yes'='slateblue')
                                ))

pltmat = log1p(mod.mat)[, umeta$ptype]
# pltmat = pltmat / apply(pltmat, 1, max)
pltmat = t(scale(t(pltmat)))
plt = Heatmap(pltmat, 
              name='Module\nscore\n(scaled)', 
              use_raster=TRUE,
              top_annotation=ha, 
              column_split=clsplit, 
              show_column_names=FALSE,
              )

h = 2.25 + 2.5 / 15 * nrow(pltmat)
w = 5 + 2.5 / 50 * ncol(pltmat)
w = ifelse(w > 25, 25, w)
pltprefix = paste0(imgpref, 'module_pseudobulk_scaled_', useset, imgsuff)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()


# Subset to top modules and plot as boxplots vs. AD variables:
# ------------------------------------------------------------
topnames = unique(statsdf$mname[(statsdf$p.adj < 0.05) & (statsdf$key != '--')])
topnames = topnames[topnames %in% rownames(mod.mat)]
scdf = data.frame(mod.mat[topnames,], check.names=FALSE)
scdf$mname = rownames(scdf)
scdf = gather(scdf, ptype, score, -mname)
scdf$ptype = gsub("^X","",scdf$ptype)
scdf = merge(scdf, umeta)  # Only ones passing cutoff:

# Plot as boxplots vs. cognition:
gplot = ggplot(scdf, aes(cell_type_high_resolution, score, color=cogdxad)) + 
    facet_wrap(~mname, scales='free_y') + 
    geom_boxplot(position=position_dodge(.75), width=.6, fill=NA, outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.5) + 
    scale_color_manual(values=colvals[['cogdxad']]) + 
    stat_compare_means(label='p.signif') + 
    theme_pubr() + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

pltprefix = paste0(imgpref, 'module_topscores_boxplots_', useset, imgsuff)
ggsave(paste0(pltprefix, '.png'), gplot, units='in', dpi=450, width=12, height=8)
ggsave(paste0(pltprefix, '.pdf'), gplot, units='in', dpi=450, width=12, height=8)


# Save the pseudobulk module scores:
# (to use in comparison across cell types)
# ----------------------------------------
scdf = data.frame(mod.mat, check.names=FALSE)
scdf$mname = rownames(scdf)
scdf = gather(scdf, ptype, score, -mname)
scdf$ptype = gsub("^X","",scdf$ptype)
scdf = merge(scdf, umeta[,c('ptype','projid','region',
                            'cell_type_high_resolution','ncell')])

scores.file = paste0(moddir, 'module_pseudobulk_scores_', 
    useset, imgsuff, '.tsv.gz')
write.table(scdf, gzfile(scores.file), quote=F, row.names=F, sep="\t")


