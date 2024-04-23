#!/usr/bin/R
# ---------------------------------------------------------
# Plot basic marker differences between astrocyte subtypes:
# Updated 11/26/2021
# ---------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)
library(patchwork)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'markers/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, srdir)
system(cmd)


# Run composition analysis for each of these subsets:
# ---------------------------------------------------
remove.batches = TRUE
suff = '_subset_final_noMB'

subset = 'Ast'
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
    ps.data = load_pseudobulk_dataset(subset, subtypes, reg.nomb)
    save(ps.data, file=psdata.rda)
} else { load(psdata.rda) }


# Further annotate the pseudo-bulk metadata:
# ------------------------------------------
pmat = ps.data$mat
umeta = ps.data$meta
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(pmat),]
# Remove very low abundance batches + will use for weight
umeta = umeta[umeta$ncell > 100,] 


# Plot specific genes for astrocytes:
# -----------------------------------
pgenes = c('DCLK1', 'CFAP54', 'CFAP61', 'CFAP43', 'CFAP44', 'SPAG1',
           'CSPP1', 'EFHC2', 'EFCAB2', 'NEK11', 'TMEM232',
           'DNAH6', 'RFX3', 'BBS9', 'CAPS', 'CAPS2','FAP','CD44',
           'BCL6','HILPDA','IRS2',
           'EMP1','GFAP','CXCL14','GLS','SLC38A1','SLC38A2','GLUL','VIM',
           'ITM2B','ITM2C','CLU','PRNP','MT1X','MT1M','MT1E','FTH1')

kmeta = umeta[umeta$ncell > 100,]
plt.mat = pmat[pgenes,kmeta$ptype]
clsplit = kmeta$cell_type_high_resolution

plt.mat = pmat[pgenes,kmeta$ptype]
ha = HeatmapAnnotation(CT=kmeta$cell_type_high_resolution, 
                       Region=kmeta$region,
                       nrad=kmeta$nrad,
                       col=list(Braak=colvals[['braaksc']],
                                Region=reg.cols,
                                nrad=colvals[['nrad']],
                                CT=tcols[subtypes]))

smat = as.matrix(log(plt.mat + 1))
# smat.scaled = t(scale(t(log(plt.mat+1))))
smat.scaled = t(scale(t(smat), center=FALSE))
Heatmap(smat.scaled, name='scaled\n logcounts', 
        # col=viridis(100),
        use_raster=TRUE,
        top_annotation=ha, 
        column_split=kmeta$nrad, 
        show_column_names=FALSE,
)



# Plot some highlighted genes from sub-astrocyte types:
# -----------------------------------------------------
# Metabolic genes:
lgenes = c('GAPDH','HILPDA','IRS2','PFKFB3','PFKP','PLOD2')

x = pmat[lgenes, umeta$ptype]
x = data.frame(t(as.matrix(x)))
x$ptype = rownames(x)
x = merge(x, umeta)
xdf = gather(x[,c('ptype',lgenes,'nrad','Apoe_e4','region','cell_type_high_resolution')], gene, val, -ptype, -nrad,-Apoe_e4,-region, -cell_type_high_resolution)

gplot = ggplot(xdf, aes(gene, val, fill=nrad, alpha=Apoe_e4)) + 
    facet_wrap(~region) +
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=colvals[['nrad']]) + 
    theme_pubr()
ggsave(paste0(imgpref, ststr, '_fullsets_metab_genes_e4_nrad_boxplots.pdf'), gplot, dpi=400, units='in', width=9,height=6)


# Ciliary genes:
cgenes = c('CFAP54','CFAP44','NEK11','EFCAB2','SPAG1','CSPP1')
x = pmat[cgenes, umeta$ptype]
x = data.frame(t(as.matrix(x)))
x$ptype = rownames(x)
x = merge(x, umeta)
xdf = gather(x[,c('ptype',cgenes,'nrad','Apoe_e4','region','cell_type_high_resolution')], gene, val, -ptype, -nrad,-Apoe_e4,-region, -cell_type_high_resolution)

gplot = ggplot(xdf, aes(gene, val, fill=nrad)) + 
    facet_grid(cell_type_high_resolution~region) +
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.75), cex=.8) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=colvals[['nrad']]) + 
    theme_pubr()
ggsave(paste0(imgpref, ststr, '_fullsets_cilia_genes_e4_nrad_boxplots.pdf'), gplot, dpi=400, units='in', width=9,height=6)


# Get the non-pseudobulk data for plotting on a cell-level:
# ---------------------------------------------------------
nmat = load_full_dataset(celltypes, subtypes, reg.nomb)
rownames(cellmeta) = cellmeta$barcode
submeta = cellmeta[colnames(nmat),]


# Plot some of these genes on the Astrocyte-only UMAP:
# ----------------------------------------------------
ind = 1:nrow(submeta)
ind = sample(ind, length(ind), replace=FALSE)
typelvls = unique(submeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
tsp.type.cols = sapply(type.cols, tsp.col)
celltype.loc = aggregate(cbind(U1, U2) ~ cell_type_high_resolution, submeta, mean)
cex = 0.025

cp = c(0, 8.5, -8.5, 0)
xlim = c(min(submeta$U1) + cp[1], min(submeta$U1) + cp[2])
ylim = c(max(submeta$U2) + cp[3], max(submeta$U2) + cp[4])
tsp.tcols = sapply(tcols, tsp.col)
tsp.rcols = sapply(reg.cols, tsp.col)

lrg.fmt = FALSE
if (lrg.fmt){
    lblset = paste0(lblset, '_medium')
    wone = 1
} else {
    wone = 1 / 2.54
}
sp = 0
png(paste0(imgpref, 'astumap_final_highres_cols_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=wone, height=wone)
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
plot(submeta$U1[ind], submeta$U2[ind], 
    col=tsp.tcols[submeta$cell_type_high_resolution[ind]],
    ylim=ylim, xlim=xlim, pch='.', cex=cex, axes=F)
dev.off()

png(paste0(imgpref, 'astumap_region_', lblset, '_wout_doublets_notext.png'), units='in', res=450, width=wone, height=wone)
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
plot(submeta$U1[ind], submeta$U2[ind],
    col=tsp.rcols[submeta$region[ind]],
    ylim=ylim, xlim=xlim, pch='.', cex=cex, axes=F)
dev.off()


palette = viridis(100)
load.colors()
palette = c('grey85', col2) 
col_fun = function(x, pal=palette, mx=NULL){
    if (is.null(mx)){mx=max(x)} else {x[x > mx] = mx}
    bin <- cut(x, seq(0, mx, length.out=length(palette)), include.lowest=T) 
    palette[bin]  }

# Plotting function for each of the module scores:
# ------------------------------------------------
# Set the x and y limits using centering params:

# Plot specific genes for astrocytes:
# -----------------------------------
pgenes = c('SLC1A2', 'GFAP', 'GRM3', 'SLC38A1', 'NRXN3', 'SLC6A11')
ntgenes = c('GRID2', 'GRIA1','GRIK4','GRIK3','NRXN3','NRXN1','NLGN1','RIMS1','SLC6A11','SLC6A1',
    'GABBR2','GABRA2','GABRB1','GRIN2A','RYR1','RYR3','ADRA1A','SLC38A1','SLC38A2')
pltmat = as.matrix(log(nmat[pgenes,submeta$barcode] + 1))

gene = 'SLC1A2'
for (gene in pgenes){
# mx = NULL
# for (gene in ntgenes){
    mx = 4
    print(gene)
    # x = log(nmat[gene,submeta$barcode] + 1)
    x = pltmat[gene,]
    ordind = ind[order(x)]
    pltprefix = paste0(imgpref, 'astumap_gene_', gene,'_', lblset, '_wout_doublets_notext')
    png(paste0(pltprefix, '.png'), units='in', res=450, width=wone, height=wone)
    sp = 0.0
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    plot(submeta$U1[ordind], submeta$U2[ordind], col=col_fun(x[ordind], mx=mx), 
        ylim=ylim, xlim=xlim, pch='.', cex=cex, axes=F)
    # mtext(gene, side=1, line=-1, font=2, cex=1.75)
    dev.off()
}


# Compare the GRM3 astrocyte markers to the MEIS2 Inhibitory markers:
# -------------------------------------------------------------------
st = 'Ast GRM3'
sub.ststr = gsub("/","_",gsub(" ","_", st))
subtype.diff.rda = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_', sub.ststr, '.rda')
load(subtype.diff.rda)

subtype.diff.tsv = paste0(srdir, 'difftl_pseudobulk_markers_', ststr, '_', sub.ststr, '.tsv.gz')
df = read.delim(gzfile(subtype.diff.tsv))

full.regdf = rbind(full.regdf, est.regdf)
}



