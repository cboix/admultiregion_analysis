#!/usr/bin/R
# ---------------------------------------------------------------------
# Compute and plot differences between subtypes for a certain celltype:
# Updated 12/02/21
# ---------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
options(width=100)


# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')  # TODO: Change dir for these analyses
plotdir = paste0(imgdir, 'metadata/')  # TODO: Change dir for these analyses
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
remove.batches = TRUE
suff = '_subset_final_noMB'

subset = 'ECneurons'
ststr = gsub("/","_", subset)

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){c(subset, remove.batches)}
source(paste0(sbindir, 'metadata_markers/load_proportions_data.R'))
subtypes = unique(ctdf$cls)
celltypes = unique(ctdf$major.celltype)


# Load in the full ast data for these subtypes:
# ---------------------------------------------
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
load(psdata.rda)

pmat = ps.data$mat
umeta = ps.data$meta

# Further annotate the pseudo-bulk metadata:
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(pmat),]


# Remove very low abundance batches:
# ----------------------------------
if (subset == 'Ast'){ 
    umeta = umeta[umeta$ncell > 100,]  # For consistency with figures
} else if (subset == 'Mic/Immune'){ 
    umeta = umeta[umeta$ncell > 0,]  # No filtering cycling!
} else {
    umeta = umeta[umeta$ncell > 25,] 
}
pmat = pmat[,umeta$ptype]
clscol = 'cell_type_high_resolution'


# Compare the celltypes to each other at the pseudo-bulk level:
# score genes for how well they distinguish the subtypes from each other
# -------------------------------------------------------------
# Filter gene list down by avg. expr:
ecut = 0.5
pcut = 1e-3
pmean = apply(pmat, 1, mean)
kept.genelist = names(pmean)[pmean > ecut]
avgidf = data.frame(symbol=names(pmean), val=pmean)
gc()

# Dummy variable of cell type:
full.regdf = c()
for (st in subtypes){
    print(st)
    sub.ststr = gsub("/","_",gsub(" ","_", st))
    subtype.diff.rda = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_', sub.ststr, '.rda')
    subtype.diff.tsv = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_', sub.ststr, '.tsv.gz')
    if (!file.exists(subtype.diff.rda)){
        umeta$is.st = umeta[[clscol]] == st
        est.regdf = c()
        for (gene in kept.genelist){
            x = pmat[gene, umeta$ptype]
            umeta$val = x
            # Regression - general effects of covariates on expression:
            fit = glm(val ~ is.st * nrad + Apoe_e4 + age_rescaled + msex + pmi, umeta,
                      weights=log(umeta$ncell), family='gaussian') # Corrects for cell ct. but not inflated
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
        ggsave(paste0(imgpref, ststr, '_subtypes_pseudobulk_ct_genes_', sub.ststr, '_volcano.png'), 
               gplot, dpi=450, units='in', width=w, height=h)
        ggsave(paste0(imgpref, ststr, '_subtypes_pseudobulk_ct_genes_', sub.ststr, '_volcano.pdf'),
               gplot, dpi=450, units='in', width=w, height=h)

        # Save:
        est.regdf$st = st
        write.table(est.regdf, gzfile(subtype.diff.tsv), quote=F, row.names=F, sep='\t')
        save(est.regdf, file=subtype.diff.rda)
    } else { load(subtype.diff.rda) }

    full.regdf = rbind(full.regdf, est.regdf)
}


# Save top differences:
# ---------------------
adf = full.regdf[full.regdf$var == 'is.stTRUE',]
adf = adf[order(adf$p),]
adf = adf[adf$Est > 0,]

topdf = c()
for (st in subtypes){
    sdf = head(adf[adf$st == st,], 100)
    topdf = rbind(topdf, sdf)
}
top.diff.tsv = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_top100_alltypes.tsv.gz')
write.table(topdf, gzfile(top.diff.tsv), quote=F, row.names=F, sep='\t')





# Plot heatmap of top types (adapt from following):
# -------------------------------------------------
top.diff.tsv = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_top100_alltypes.tsv.gz')
adf = read.delim(gzfile(top.diff.tsv), header=T)
adf = adf[grep("^MT-",adf$symbol, invert=TRUE),]

pgenes = c()
pset = c()
subtypecols = c(tcols[subtypes])


ntop = 10
for (st in subtypes){
    sdf = adf[adf$st == st,]
    selgenes = head(unique(sdf$symbol),ntop)
    if (subset %in% c('Ast', 'Opc','Oli')){
        stgene = sub(".* ","", st)  # Ensure we keep subtype marker in plot
        selgenes = head(unique(c(stgene, selgenes)), ntop)
    } else if (subset == 'Inh'){
        if (st == 'Thalamus'){
            selgenes = head(unique(c('GRM3', selgenes)), ntop)
        }
    }
    pgenes = c(pgenes, selgenes)
    pset = c(pset, rep(st, ntop))
}


kmeta = umeta[umeta$ncell > 10,]
kmeta$region = sapply(kmeta$ptype, function(x){sub(".*_","",x)})
plt.mat = pmat[pgenes,kmeta$ptype]
clsplit = kmeta[[clscol]]
clsplit = sub(" ", "\n", sub("Exc ","", clsplit))
ha = HeatmapAnnotation(
    CT=kmeta[[clscol]], 
    Region=kmeta$region,
    col=list(Braak=colvals[['braaksc']],
        Region=reg.cols,
        CT=subtypecols),
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(ux / 1.25, 'mm'),
    gap = unit(0, "mm"))

udsplit = pset
hb = rowAnnotation(Set=pset,
                   col=list(Set=subtypecols),
                   annotation_name_gp = gpar(fontsize=5),
                   simple_anno_size = unit(ux / 1.25, 'mm'),
                   gap = unit(0, "mm"))

smat = as.matrix(log(plt.mat + 1))
# smat.scaled = t(scale(t(smat), center=FALSE))
smat.scaled = sweep(smat, 1, apply(smat, 1, max), '/')

ux = 1.5
pltmat = smat.scaled
plt = Heatmap(pltmat,
    name='scaled\n logcounts', 
    use_raster=TRUE,
    width = ncol(pltmat)*unit(ux, "mm") / 2, 
    height = nrow(pltmat)*unit(ux, "mm"),
    # column_title='Average expression in each individual',
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    top_annotation=ha, 
    column_split=clsplit, 
    show_column_names=FALSE,
    row_split=udsplit,
    right_annotation=hb
)

pltprefix = paste0(imgpref, ststr, '_topDE_heatmap_n', ntop)
h = 1 + 1 / 15 * nrow(pltmat)
w = 3 + 1 / 15 * ncol(pltmat) / 2
saveHeatmap(plt, pltprefix, w=w, h=h)


# Only selected EC neurons:
# -------------------------
pgenes = c()
pset = c()
ntop = 30
for (st in c('Exc RELN GPC5', 'Exc RELN COL5A2')){
    sdf = adf[adf$st == st,]
    selgenes = head(unique(sdf$symbol),ntop)
    if (subset %in% c('Ast', 'Opc','Oli')){
        stgene = sub(".* ","", st)  # Ensure we keep subtype marker in plot
        selgenes = head(unique(c(stgene, selgenes)), ntop)
    } else if (subset == 'Inh'){
        if (st == 'Thalamus'){
            selgenes = head(unique(c('GRM3', selgenes)), ntop)
        }
    }
    pgenes = c(pgenes, selgenes)
    pset = c(pset, rep(st, ntop))
}


kmeta = umeta[umeta$ncell > 10,]
kmeta$region = sapply(kmeta$ptype, function(x){sub(".*_","",x)})
plt.mat = pmat[pgenes,kmeta$ptype]
clsplit = kmeta[[clscol]]
clsplit = sub(" ", "\n", sub("Exc ","", clsplit))
ha = HeatmapAnnotation(
    CT=kmeta[[clscol]], 
    Region=kmeta$region,
    col=list(Braak=colvals[['braaksc']],
        Region=reg.cols,
        CT=subtypecols),
    annotation_name_gp = gpar(fontsize=5),
    simple_anno_size = unit(ux / 1.25, 'mm'),
    gap = unit(0, "mm"))

udsplit = pset
hb = rowAnnotation(Set=pset,
                   col=list(Set=subtypecols),
                   annotation_name_gp = gpar(fontsize=5),
                   simple_anno_size = unit(ux / 1.25, 'mm'),
                   gap = unit(0, "mm"))

smat = as.matrix(log(plt.mat + 1))
# smat.scaled = t(scale(t(smat), center=FALSE))
smat.scaled = sweep(smat, 1, apply(smat, 1, max), '/')

ux = 1.5
pltmat = smat.scaled
plt = Heatmap(pltmat,
    name='scaled\n logcounts', 
    use_raster=TRUE,
    width = ncol(pltmat)*unit(ux, "mm") / 2, 
    height = nrow(pltmat)*unit(ux, "mm"),
    # column_title='Average expression in each individual',
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    top_annotation=ha, 
    column_split=clsplit, 
    show_column_names=FALSE,
    row_split=udsplit,
    right_annotation=hb
)

pltprefix = paste0(imgpref, ststr, '_topDE_heatmap_RELNsubtypes_n', ntop)
h = 1 + 1 / 15 * nrow(pltmat)
w = 3 + 1 / 15 * ncol(pltmat) / 2
saveHeatmap(plt, pltprefix, w=w, h=h)


# RAW:
ux = 1.5
pltmat = smat
plt = Heatmap(pltmat,
    name='scaled\n logcounts', 
    use_raster=TRUE,
    width = ncol(pltmat)*unit(ux, "mm") / 2, 
    height = nrow(pltmat)*unit(ux, "mm"),
    col=viridis(100),
    # column_title='Average expression in each individual',
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    top_annotation=ha, 
    column_split=clsplit, 
    show_column_names=FALSE,
    row_split=udsplit,
    right_annotation=hb
)

pltprefix = paste0(imgpref, ststr, '_topDE_heatmap_RELNsubtypes_raw_n', ntop)
h = 1 + 1 / 15 * nrow(pltmat)
w = 3 + 1 / 15 * ncol(pltmat) / 2
saveHeatmap(plt, pltprefix, w=w, h=h)


