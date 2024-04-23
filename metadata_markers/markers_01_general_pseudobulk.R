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
library(ggplot2)
library(ggrepel)
library(ggpubr)
# library(ggpmisc)
library(patchwork)

library(ComplexHeatmap)
library(circlize)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


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

subset = 'Oli'
subset = 'ECneurons'
subset = 'Vasc/Epithelia'
subset = 'Mic/Immune'
subset = 'Inh'
subset = 'Ast'
subset = 'Exc'
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
if (!file.exists(psdata.rda)){
    if (subtype == 'OpcOli'){
        # Subtypes (TODO: if mult cell type, load all + subset):
        ps.data = load_pseudobulk_dataset('Opc', 'Opc', reg.nomb)
        ps.data2 = load_pseudobulk_dataset('Oli', 'Oli', reg.nomb)
        # Merge the two datasets:
        ps.data$mat = cbind(ps.data$mat, ps.data2$mat)
        ps.data$meta = rbind(ps.data$meta, ps.data2$meta)
        rm(ps.data2)
    } else {
        ps.data = load_pseudobulk_dataset(celltypes, subtypes, reg.nomb)
    }
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }

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


# For excitatory and inhibitory, reduce to either subsets or regions:
# -------------------------------------------------------------------
if (subset %in% c('Exc', 'Inh')){
    # Load subtype information:
    anndf = read.delim(paste0(sdbdir, subset, '_subclass_annotation.tsv'), header=T)
    major.subclasses = sort(unique(anndf$subclass))

    # Average at that subtype's level:
    umeta = merge(umeta, anndf[,c('cell_type_high_resolution','subclass')], all.x=TRUE)
    umeta$stype = with(umeta, paste(projid, subclass, region, sep="_"))
    tform = make.tform(umeta$stype, u=unique(umeta$stype), norm=FALSE)

    # Normalize relative to ncell and transform:
    tform = sweep(tform, 1, umeta$ncell, '*')
    tform = sweep(tform, 2, apply(tform, 2, sum), '/')
    pmat = pmat[,umeta$ptype] %*% tform

    # Reduce metadata again:
    umeta = aggregate(ncell ~ stype + subclass + nrad + 
                      Apoe_e4 + age_rescaled + msex + pmi, umeta, sum)
    names(umeta)[1] = 'ptype'
    clscol = 'subclass'
    subtypes = major.subclasses
} else { clscol = 'cell_type_high_resolution' }



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
adf = adf[grep("^MT-",adf$symbol, invert=TRUE),]

pgenes = c()
pset = c()
if (subset == 'Ast'){
    pgenes = c('SLC1A2','SLC1A3','GFAP','AQP4')
    pset = c('All','All','All','All')
    subtypecols = c(tcols[subtypes], 'All'='white')
} else if (subset %in% 'Inh'){ 
    j = 14
    subtypecols = snap.cols[j:(j + length(subtypes) - 1)]
    names(subtypecols) = subtypes
} else { subtypecols = c(tcols[subtypes]) }


ntop = 4
for (st in subtypes){
    sdf = adf[adf$st == st,]
    selgenes = head(unique(sdf$symbol),ntop)
    selset = rep(st, ntop)
    if (subset %in% c('Ast', 'Opc','Oli')){
        stgene = sub(".* ","", st)  # Ensure we keep subtype marker in plot
        selgenes = head(unique(c(stgene, selgenes)), ntop)
        if (st == 'Ast LUZP2'){
            selgenes = head(unique(c('LGR6', selgenes)), ntop + 1)
            selset = c(selset, selset[1])
        }
    } else if (subset == 'Inh'){
        if (st == 'Thalamus'){
            selgenes = head(unique(c('GRM3', 'MEIS2', selgenes)), ntop +1)
            selset = c(selset, selset[1])
        }
    }
    pgenes = c(pgenes, selgenes)
    pset = c(pset, selset) 
}


kmeta = umeta[umeta$ncell > 10,]
kmeta$region = sapply(kmeta$ptype, function(x){sub(".*_","",x)})
plt.mat = pmat[pgenes,kmeta$ptype]
clsplit = kmeta[[clscol]]
ha = HeatmapAnnotation(CT=kmeta[[clscol]], 
                       Region=kmeta$region,
                       col=list(Braak=colvals[['braaksc']],
                                Region=reg.cols,
                                CT=subtypecols))
udsplit = pset
hb = rowAnnotation(Set=pset,
                   col=list(Set=subtypecols))

smat = as.matrix(log(plt.mat + 1))
# smat.scaled = t(scale(t(smat), center=FALSE))
smat.scaled = sweep(smat, 1, apply(smat, 1, max), '/')

plt = Heatmap(smat.scaled, 
    name='scaled\n logcounts', 
    use_raster=TRUE,
    top_annotation=ha, 
    column_split=clsplit, 
    show_column_names=FALSE,
    row_split=udsplit,
    right_annotation=hb
)

png(paste0(imgpref, ststr, '_topDE_heatmap_normalized_individ.png'), res=400, units='in', width=11, height=5)
print(plt)
dev.off()

pdf(paste0(imgpref, ststr, '_topDE_heatmap_normalized_individ.pdf'), width=11, height=5)
print(plt)
dev.off()


