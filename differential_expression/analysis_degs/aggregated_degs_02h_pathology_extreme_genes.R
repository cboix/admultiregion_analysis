#!/usr/bin/R
# -----------------------------------------------------
# Plots/stats for genes more of one pathology than other
# - Show overlap with modules, GO, ndegs, etc.
# Updated: 06/15/23
# -----------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(ggrastr)
print(version)
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
setnames = c("Mic_Immune_Mic"='Microglial', "Ast_Ast"='Astrocyte', 
    "Opc_Opc"='OPC', "Oli_Oli"='Oligodendrocyte', 
    'Inh_Inh'='Inhibitory neuron','Exc_Exc'='Excitatory neuron')

pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
degcols = c('NS'='grey90','Down'=col.paired[2],'Up'=col.paired[6])

# Read in genes:
res.file = paste0(regdir, 'allmethods.major.plaq_nft.extremegenes.tsv')
resdf = read.delim(res.file)
dir.ind = apply(resdf[,c('sig.nft.up','sig.plaq.up', 'sig.nft.dw', 'sig.plaq.dw')], 1, which.max)
resdf$sig.dir = ifelse(dir.ind > 2, 'Down', 'Up')


# Explore these genes:
# --------------------
# Shared genes:
shdf = agg.rename(nft ~ dir + sig.dir + gene, resdf, length, 'nhit')
ctshdf = agg.rename(gene ~ nhit + dir + sig.dir, shdf, length, 'freq')

gp = ggplot(ctshdf, aes(nhit, freq, fill=sig.dir, alpha=dir)) +
    geom_bar(stat='identity', position='dodge') + 
    labs(x='Number of cell types', y='Number of genes') +
    # scale_y_continuous(expand=c(0,0)) + 
    # scale_fill_manual(values=c('grey75', 'grey25'), name='DEG for:') + 
    theme_pubr()

meandf = aggregate(cbind(nft, plaq_n, nft.resid) ~ dir + gene, resdf, mean)
shdf = merge(shdf, meandf)
shdf = shdf[order(abs(shdf$nft.resid), decreasing=T),]
shdf = shdf[order(shdf$nhit, decreasing=T),]
shdf = shdf[shdf$nhit > 1,]
head(shdf, 20)

length(unique(shdf$gene[shdf$dir == 'Plaque']))
length(unique(shdf$gene[shdf$dir == 'NFT']))
shared.extr.file = paste0(regdir, 'allmethods.major.plaq_nft.shared.top.extremegenes.tsv')
write.table(shdf, shared.extr.file, quote=F, sep="\t", row.names=F)
# NOTE: average is over tested DEGs only; some may look odd over balance

dir.cols = c('NFT'=col.paired[8], 'Plaque'=col.paired[4]) 
gp = ggplot(shdf, aes(nft, plaq_n, color=dir, size=as.character(nhit), label=gene)) +
    geom_point() + 
    geom_text_repel(box.padding=0.075, cex=3, max.overlaps=20, min.segment.length=0.1) + 
    geom_abline(intercept=0, slope=1, lty='dotted') +
    geom_hline(yintercept=0, lty='dashed') + 
    geom_vline(xintercept=0, lty='dashed') + 
    labs(x='NFT effect size', y='Plaque effect size') +
    scale_size_manual(values=c('2'=.5,'3'=1,'4'=1.5), name='# cell types DEG:') +
    scale_color_manual(values=dir.cols, name='Higher in:') +
    theme_pubr()
pltprefix = paste0(imgpref, 'path_extreme_genes.volcano')
saveGGplot(gp, pltprefix, w=3.5, h=3.5)
saveGGplot(gp, paste0(pltprefix, '_large'), w=6, h=6)

# Look up which cell types:
gene = 'CRYAB'
resdf[resdf$gene == gene,]


# GO enrichment on the shared genes:
# ----------------------------------
shared.enr.file = paste0(regdir, 'allmethods.major.plaq_nft.shared.top.extremegenes.enrichments.rda')
if (!file.exists(shared.enr.file)){
    genesets = list(NFT = unique(shdf$gene[shdf$dir == 'NFT']),
        Plaque = unique(shdf$gene[shdf$dir == 'Plaque']))
    sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
    gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
        ordered_query=FALSE, multi_query=TRUE, sources=sources)

    # Get the genes in the term intersections (slow):
    allgenes = unique(unlist(genesets))
    gp2.res.int = gprofiler2::gost(allgenes, organism='hsapiens',
        ordered_query=FALSE, multi_query=FALSE, user_threshold=1,
        sources = sources, evcodes=TRUE)
    intdf = gp2.res.int$result
    save(gp2.result, intdf, file=shared.enr.file)
} else { load(shared.enr.file) }



suffix = 'path_extreme_genes.enrichments'
NTOP = 10
gpdf = gp2.result$result
gpdf = gpdf[gpdf$term_size < 1000,]
gpdf = pruneWithInt(gpdf, intdf)
pltprefix = paste0(imgpref, suffix)
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)

gpdf = gp2.result$result
gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
gpdf = pruneWithInt(gpdf, intdf)
pltprefix = paste0(imgpref, suffix, '_src')
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)

gpdf = gp2.result$result
gpdf = gpdf[gpdf$term_size < 500,]
gpdf = pruneWithInt(gpdf, intdf, cutoff=0.5)
pltprefix = paste0(imgpref, suffix, '_small')
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)


# Plot these GO terms as a barplot:
# ---------------------------------
# NOTE: using small results (terms < 500 genes) with pruning
df = data.frame(subpmat)
df$term = factor(rownames(subpmat), levels=rev(rownames(subpmat)))
df$dir = ifelse(df$NFT > df$Plaque, 'NFT', 'Plaque')
df$val = ifelse(df$NFT > df$Plaque, df$NFT, df$Plaque)

gp = ggplot(df, aes(term, val, fill=dir)) + 
    facet_wrap(~dir, ncol=1, scales='free_y') +
    geom_bar(stat='identity') + 
    labs(x='GO term', y='-log10p') + 
    scale_fill_manual(values=dir.cols, name='Higher in:') +
    theme_pubr() + 
    coord_flip()
pltprefix = paste0(imgpref, 'path_extreme_genes.shared.enrichments.barplot')
saveGGplot(gp, pltprefix, w=5.5, h=5)


# V2 (up/down) GO enrichment on the shared genes:
# -----------------------------------------------
shared.enrv2.file = paste0(regdir, 'allmethods.major.plaq_nft.shared.top.extremegenes.updown.enrichments.rda')
if (!file.exists(shared.enr.file)){
    genesets = list()
    for (path in unique(shdf$dir)){
        for (sig.dir in unique(shdf$sig.dir)){
            tag = paste0(path, ' (', sig.dir, ')')
            genesets[[tag]] = shdf$gene[shdf$dir == path & shdf$sig.dir == sig.dir]
        } 
    }
    sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
    gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
        ordered_query=FALSE, multi_query=TRUE, sources=sources)

    # Get the genes in the term intersections (slow):
    allgenes = unique(unlist(genesets))
    gp2.res.int = gprofiler2::gost(allgenes, organism='hsapiens',
        ordered_query=FALSE, multi_query=FALSE, user_threshold=1,
        sources = sources, evcodes=TRUE)
    intdf = gp2.res.int$result
    save(gp2.result, intdf, file=shared.enrv2.file)
} else { load(shared.enrv2.file) }

suffix = 'path_extreme_genes.enrichments.updown'
NTOP = 10
gpdf = gp2.result$result
gpdf = gpdf[gpdf$term_size < 1000,]
pltprefix = paste0(imgpref, suffix)
gpdf = pruneWithInt(gpdf, intdf)
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)

gpdf = gp2.result$result
gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
pltprefix = paste0(imgpref, suffix, '_src')
gpdf = pruneWithInt(gpdf, intdf)
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)

gpdf = gp2.result$result
gpdf = gpdf[gpdf$term_size < 500,]
pltprefix = paste0(imgpref, suffix, '_small')
gpdf = pruneWithInt(gpdf, intdf)
subpmat = gpPvalMatrix(gpdf, genesets, ntop=NTOP)
plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)


# Relative numbers of each type:
# ------------------------------
resdf$uq.gene = !(resdf$gene %in% shdf$gene)
# Astrocytes have highest plaque direction:
table(resdf[,c('dir','set', 'uq.gene')])
ctdf = aggregate(gene ~ dir + set + uq.gene, resdf, length)
adf = aggregate(gene ~ dir, shdf, length)
adf$set = 'Shared'
adf$uq.gene = TRUE
ctdf = rbind(ctdf, adf[, colnames(ctdf)])
ctdf$set = factor(ctdf$set, levels=rev(c(sort(keep.sets), 'Shared')))
ctdf$dir = factor(ctdf$dir, levels=c('NFT','Plaque'))

gp = ggplot(ctdf[ctdf$uq.gene,], aes(set, gene, fill=dir)) + 
    geom_bar(stat='identity', position='dodge') + 
    labs(x='Major cell type', y='# genes') + 
    scale_fill_manual(values=dir.cols, name='Higher in:') +
    theme_pubr() + coord_flip()
pltprefix = paste0(imgpref, 'path_extreme_genes.unique.ngene.barplot')
saveGGplot(gp, pltprefix, w=2.5, h=2.5)


# Plots for specific cell types / GO analysis:
# --------------------------------------------
# Plot top genes (?) --> Ast works.
set = 'Ast_Ast'
for (set in keep.sets){
    print(set)
    subdf = resdf[resdf$set == set,]
    dir.cols = c('NFT'=col.paired[8], 'Plaque'=col.paired[4], 'Non-unique'='grey75') 
    subdf$dir = ifelse(subdf$gene %in% shdf$gene, 'Non-unique', subdf$dir)

    gp = ggplot(subdf[subdf$dir != 'Non-unique',], aes(nft.resid, plaq_n, color=dir, label=gene)) +
        geom_point(cex=.5) + 
        geom_text_repel(box.padding=0.075, cex=3, max.overlaps=20, min.segment.length=0.1) + 
        geom_abline(intercept=0, slope=1, lty='dotted') +
        geom_hline(yintercept=0, lty='dashed') + 
        geom_vline(xintercept=0, lty='dashed') + 
        labs(x='Residual NFT effect size', y='Plaque effect size') +
        # scale_size_manual(values=c('2'=.5,'3'=1,'4'=1.5), name='# cell types DEG:') +
        scale_color_manual(values=dir.cols, name='Higher in:') +
        theme_pubr()
    pltprefix = paste0(imgpref, 'path_extreme_genes_resid.', set, '.volcano')
    saveGGplot(gp, pltprefix, w=3.5, h=3.5)
    saveGGplot(gp, paste0(pltprefix, '_large'), w=6, h=6)

    gp = ggplot(subdf[subdf$dir != 'Non-unique',], aes(nft, plaq_n, color=dir, label=gene)) +
        geom_point(cex=.5) + 
        geom_text_repel(box.padding=0.075, cex=3, max.overlaps=20, min.segment.length=0.1) + 
        geom_abline(intercept=0, slope=1, lty='dotted') +
        geom_hline(yintercept=0, lty='dashed') + 
        geom_vline(xintercept=0, lty='dashed') + 
        labs(x='NFT effect size', y='Plaque effect size') +
        # scale_size_manual(values=c('2'=.5,'3'=1,'4'=1.5), name='# cell types DEG:') +
        scale_color_manual(values=dir.cols, name='Higher in:') +
        theme_pubr()
    pltprefix = paste0(imgpref, 'path_extreme_genes.', set, '.volcano')
    saveGGplot(gp, pltprefix, w=3.5, h=3.5)
    saveGGplot(gp, paste0(pltprefix, '_large'), w=6, h=6)


}


# Plots for specific cell types / GO analysis:
# --------------------------------------------
# Plot top genes (?) --> Ast works.
set = 'Ast_Ast'
for (set in keep.sets){
    print(set)
    subdf = resdf[resdf$set == set,]

    genesets = list()
    for (path in unique(subdf$dir)){
        for (sig.dir in unique(subdf$sig.dir)){
            tag = paste0(path, ' (', sig.dir, ')')
            genesets[[tag]] = subdf$gene[subdf$dir == path & subdf$sig.dir == sig.dir]
        } 
    }

    # Get enrichments for shared and specific DEGs:
    # ---------------------------------------------
    fullenr.file = paste0(regdir, 'extremeDEGs_enrichments_', set, '.rda')
    if (!file.exists(fullenr.file)){
        # Remove empty sets:
        dl = c(lapply(genesets, length))
        keep = names(dl)[dl > 0]
        genesets = genesets[keep]

        print('--GO enrichments')
        sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
        gp2.result = gprofiler2::gost(genesets, organism='hsapiens',
            ordered_query=FALSE, multi_query=TRUE,
            sources = sources)
        save(gp2.result, file=fullenr.file)
    } else { load(fullenr.file) }


    # Plot enrichments as heatmaps:
    # -----------------------------
    ntop = 4
    gpdf = gp2.result$result
    gpdf$nc = nchar(gpdf$term_name)
    gpdf = gpdf[gpdf$nc <= 40,]
    gpdf = gpdf[gpdf$term_size < 1000,]
    pltprefix = paste0(imgpref, 'extremeDEGs_enrichments_', set)
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=ntop)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)

    gpdf = gp2.result$result
    gpdf$nc = nchar(gpdf$term_name)
    gpdf = gpdf[gpdf$nc <= 40,]
    gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
    pltprefix = paste0(imgpref, 'extremeDEGs_enrichments_', set, '_src')
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=ntop)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)

    gpdf = gp2.result$result
    gpdf$nc = nchar(gpdf$term_name)
    gpdf = gpdf[gpdf$nc <= 40,]
    gpdf = gpdf[gpdf$term_size < 500,]
    pltprefix = paste0(imgpref, 'extremeDEGs_enrichments_', set, '_small')
    subpmat = gpPvalMatrix(gpdf, genesets, ntop=ntop)
    plt = plotGpPvalMatrix(subpmat, pltprefix, cluster_columns=FALSE, use_raster=FALSE)
}


# Show cell-type specific overlap with modules:
# (to motivate modules analysis)
# ---------------------------------------------
# Load modules:
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Plot top genes (?) --> Ast works.
set = 'Ast_Ast'
use.core = TRUE
for (set in keep.sets){
    setstr = sub("_[A-Za-z]+$","", set)
    cat(set, setstr, '\n')
    coremap = cmlist[[setstr]]
    genemap = gmlist[[setstr]]
    if (use.core){ usemap = coremap } else { usemap = genemap }

    subdf = resdf[resdf$set == set,]
    dir.cols = c('NFT'=col.paired[8], 'Plaque'=col.paired[4], 'Non-unique'='grey75') 
    # subdf$uq.gene = !(subdf$gene %in% shdf$gene)
    # subdf = subdf[subdf$uq.gene,]
    subdf$module = usemap[subdf$gene]
    subdf = subdf[!is.na(subdf$module),]
    kept.genes = unique(subdf$gene)

    # Count # genes per module:
    ndf = data.frame(table(usemap))
    names(ndf) = c('module', 'nmod')
    mod.cutoff = 10 # Modules with at least 10 genes
    kept.modules = ndf$module[ndf$nmod >= mod.cutoff]

    # Enrichment:
    hgdf = agg.rename(gene ~ module + dir + sig.dir, subdf, length, 'nint')
    hgdf = merge(hgdf, agg.rename(gene ~ dir + sig.dir, subdf, length, 'nsig'))
    hgdf = merge(hgdf, ndf)
    hgdf = hgdf[hgdf$nmod >= mod.cutoff,]
    hgdf$ntot = sum(usemap %in% kept.modules)

    # Run hypergeometric tests:
    # -------------------------
    pvdf = hgdf[,c('nint','nmod', 'nsig','ntot')]
    pout <- apply(pvdf, 1, run.hyper)
    hgdf$p = pout
    hgdf$lp = -log10(pout)
    hgdf = hgdf[order(hgdf$p),]
    hgdf$padj = p.adjust(hgdf$p,'BH') # Correct by BH
    hgdf$log2FC = with(hgdf, log2((nint / nsig) / (nmod / ntot)))


    # Plot these modules:
    # -------------------
    plt.modules = hgdf$module[hgdf$p < 0.1]
    pltdf = hgdf[hgdf$module %in% plt.modules,]
    pltdf$lab = paste0(pltdf$dir, '-', pltdf$sig.dir)


    # TOP 3 DEGs
    subdf$ord = 1:nrow(subdf)
    topdf = merge(aggregate(p ~ module , pltdf[pltdf$sig.dir == 'Up',], min),
        pltdf[,c('p','module','sig.dir','dir')])
    topdf = merge(topdf, subdf)
    topdf = topdf[order(topdf$ord),]
    aggregate(gene ~ module, topdf, function(x){
        paste(head(x, 4), collapse=', ') })

    cmat = pivot.tomatrix(pltdf[,c('module','lab','log2FC')], 'lab', 'log2FC')
    pmat = pivot.tomatrix(pltdf[,c('module','lab','p')], 'lab', 'p')
    cmat[is.na(cmat)] = 0
    pmat[is.na(pmat)] = 1
    cn = rev(colnames(cmat))
    cmat = cmat[,cn]
    pmat = pmat[,cn]

    mx = 3
    col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))

    ux = 1.5
    plt = Heatmap(cmat,
        use_raster=FALSE,
        name='log2FC',
        col=col_fun,
        cluster_columns=FALSE,
        cluster_rows=TRUE,
        width = ncol(cmat)*unit(ux * 1.75, "mm"), 
        height = nrow(cmat)*unit(ux, "mm"),
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lwd=.5),
        column_title=paste0('Modules (', setstr, ')'),
        cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
            p = pmat[i,j]
            if (!is.na(p)){
                ann = ifelse(p < 0.1, ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***', '**'), '*'), '.'),'')
                grid.text(ann, x, y,gp=gpar(fontsize=gridtxt.fs))
            }
        })

    pltprefix = paste0(imgpref, 'extremeDEGs_moduleEnr_', set)
    h = 1 + 1 / 15 * nrow(cmat)
    w = 1 + 1 / 15 * ncol(cmat) * 1.5
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

