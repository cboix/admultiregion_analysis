#!/usr/bin/R
# -------------------------------------------------
# Make a heatmap across all regional shared/specific
# Updated: 11/06/23
# -------------------------------------------------
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
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
remove.shared = TRUE
run.intersections = FALSE
degcols = c('NS'='grey90','Down'=col.paired[2],'Up'=col.paired[6])
denrcols = c('NS'='grey90',
    'Down (1-2)'=col.paired[2],
    'Down (3+)'=col.paired[1],
    'Down (multi-CT)'='slateblue4',
    'Up (1-2)'=col.paired[6],
    'Up (3+)'=col.paired[5],
    'Up (multi-CT)'='brown4')


# Load full results:
# ------------------
path = 'nrad'
mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)

kept.cols = c('gene','col_nm','path','region', 'logFC_nb')
alldf = c()
for (set in keep.sets){
    print(set)
    setdf = setdflist[[set]][, kept.cols]
    setdf = setdf[setdf$col_nm != 0,]
    setdf$set = set
    alldf = rbind(alldf, setdf)
}


# Load in the annotated DEG results:
# ----------------------------------
mstr = paste0('allmethods.regional_', path)
outpref = paste0(regdir, mstr, '.merged.sharedspecific')
allresdf = readRDS(paste0(outpref, '.rds'))
allresdf$sr = paste0(allresdf$set, '@', allresdf$region)

degenes = unique(allresdf[allresdf$denr != 'NS', 'gene'])  # 7.3k genes
df = unique(allresdf[allresdf$gene %in% degenes, c('denr', 'gene')])  # 18k, non-unique assignment, plotting how

regset = c('allregions', 'EC', 'HC', 'TH', 'AG', 'MT', 'PFC')
regset = regset[regset %in% allresdf$region]

# Plot heatmap for sets separately:
set = 'Mic_Immune_Mic'
for (set in keep.sets){

    print(set)
    subdf = allresdf[allresdf$set == set,]
    degenes = unique(subdf[(subdf$denr != 'NS'), 'gene'])
    # df = subdf[(subdf$gene %in% degenes), ]
    df = allresdf[(allresdf$gene %in% degenes) & (allresdf$set %in% keep.sets),]
    cmat = pivot.tomatrix(df[,c('gene', 'sr', 'logFC_nb')], 'sr', 'logFC_nb')

    # Columns:
    cn = expand.grid(keep.sets, regset)
    cn = paste0(cn[,1], '@', cn[,2])
    cmat = cmat[,cn]

    # Reorder the matrix:
    nmat = cmat
    nmat[is.na(nmat)] = 0
    nmat = reord(nmat)
    # nmat = t(diag.mat2(t(nmat))[[1]])
    cmat = cmat[rownames(nmat),]

    # Split rows and columns:
    proggenes = unique(subdf[grep('3+', subdf$denr),'gene'])
    mctgenes = unique(subdf[grep('multi-CT', subdf$denr),'gene'])
    rowsplit = ifelse(rownames(cmat) %in% mctgenes, 'Cross-CT', ifelse(rownames(cmat) %in% proggenes, 'Consistent', 'Region-specific'))
    colsplit = sub("@.*","", colnames(cmat))


    scale = 1 / 20
    ux = 1.5
    ht = Heatmap(cmat,
        use_raster=FALSE,
        name='log2FC',
        # col=col_fun,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        column_split=colsplit,
        show_row_names=FALSE,
        row_split=rowsplit,
        width = ncol(cmat)*unit(ux, "mm"), 
        height = nrow(cmat)*unit(ux, "mm") * scale,
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lwd=.5),
        column_title=paste0('Cell type + region'))

    pltprefix = paste0(imgpref, 'shared_specific_', set)
    h = 1 + 1 / 15 * nrow(cmat) * scale
    w = 1 + 1 / 15 * ncol(cmat)
    saveHeatmap(ht, pltprefix, w=w, h=h)

}



# Plot the top consistent and specific genes for each cell type:
#---------------------------------------------------------------
# Plot heatmap for sets separately:
set = 'Mic_Immune_Mic'
for (set in keep.sets){
    print(set)
    suff = paste0(set, '.', path)
    subdf = allresdf[allresdf$set == set,]

    # Label gene sets:
    degenes = unique(subdf[(subdf$denr != 'NS'), 'gene'])
    proggenes = unique(subdf[grep('3+', subdf$denr),'gene'])
    rsgenes = unique(subdf[grep('1-2', subdf$denr),'gene'])
    mctgenes = unique(subdf[grep('multi-CT', subdf$denr),'gene'])
    gset = ifelse(degenes %in% proggenes, 'Consistent', ifelse(degenes %in% mctgenes, 'Cross-CT', 'Region-specific'))
    names(gset) = degenes
    df = subdf[subdf$gene %in% degenes,]
    df$geneset = gset[df$gene]
    df = df[order(df$log10p_nm, decreasing=T),]
    
    # Select genes to plot:
    reg.sigdf = df[(df$col_nm != 0) & (df$region != 'allregions'),]
    topdf = aggregate(log10p_nm ~ gene + geneset, reg.sigdf, max)
    topdf = merge(topdf, reg.sigdf[,c('log10p_nm', 'gene', 'geneset', 'de', 'region')])
    topdf = unique(topdf)
    topdf = topdf[order(topdf$log10p_nm, decreasing=T),]
    # NOTE: Remove MT-genes, not as interesting:
    topdf = topdf[grep("^MT-", topdf$gene, invert=TRUE),]


    # Plot consistent and specific genes:
    # -----------------------------------
    NTOP = 16
    NREG = 4
    for (plot.uponly in c(TRUE, FALSE)){
        if (plot.uponly){
            updf = topdf[topdf$de == 'Up',]
            cg = head(unique(updf$gene[updf$geneset == 'Consistent'], NTOP))
            # sg = head(updf$gene[updf$geneset == 'Region-specific'], NTOP)
            sg = sapply(regset[2:length(regset)], function(x){
                head(updf$gene[(updf$geneset == 'Region-specific') & (topdf$region == x)], NREG)})
            pltprefix = paste0(imgpref, 'shared_specific_topgenes_heatmap.', suff, '.uponly')
        } else {
            cg = head(unique(topdf$gene[topdf$geneset == 'Consistent']), NTOP)
            # sg = head(topdf$gene[topdf$geneset == 'Region-specific'], NTOP)
            sg = sapply(regset[2:length(regset)], function(x){
                head(topdf$gene[(topdf$geneset == 'Region-specific') & (topdf$region == x)], NREG)})
            sg = c(unique(sg))
            pltprefix = paste0(imgpref, 'shared_specific_topgenes_heatmap.', suff)
        }
        seldf = df[df$gene %in% c(cg, sg),]

        # Make matrices for plotting:
        cmat = pivot.tomatrix(seldf[,c('gene','region','logFC_nb')], 'region', 'logFC_nb')
        pmat = pivot.tomatrix(seldf[,c('gene','region','col_nm')], 'region', 'col_nm')
        cmat = cmat[c(cg, sg), regset]
        pmat = pmat[c(cg, sg), regset]

        # Order in same way:
        nmat = cmat[cg,]
        nmat[is.na(nmat)] = 0
        nmat = reord(nmat)
        # nmat = t(diag.mat2(t(nmat))[[1]])
        cmat = cmat[c(rownames(nmat), sg),]
        pmat = pmat[rownames(cmat),]
        rowsplit = ifelse(rownames(cmat) %in% cg, 'Consistent', 'Region-specific')

        mx = .5
        if (path %in% c('nft','plaq_n', 'plaq_d')){ mx = mx / 20 }
        col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
        ux = 1.5
        ht = Heatmap(cmat, 
            use_raster=FALSE, 
            name='logFC_nb',
            col=col_fun,
            column_title=suff,
            cluster_columns=FALSE, 
            cluster_rows=FALSE,
            cluster_row_slices=FALSE,
            row_split=rowsplit, 
            row_dend_width = unit(.25, "cm"),
            column_dend_height = unit(.25, "cm"),
            row_dend_gp = gpar(lwd=.5),
            column_dend_gp = gpar(lwd=.5),
            border_gp = gpar(col="black", lwd=.5),
            width=ncol(cmat) * unit(ux, 'mm'),
            height=nrow(cmat) * unit(ux, 'mm'),
            cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                p = pmat[i,j]
                if (!is.na(p) & (p != 0)){
                    grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.25))
                }
            })
        h = .5 + 1 / 15 * nrow(cmat)
        w = 1.5 + 1 / 15 * ncol(cmat)
        saveHeatmap(ht, pltprefix, w=w, h=h)
    }
}


# Get the broadly shared (3+ ct) genes as well:
# ---------------------------------------------
mctgenes = unique(allresdf[grep('multi-CT', allresdf$denr),'gene'])
subdf = allresdf[allresdf$gene %in% mctgenes,]
subdf$lp = ifelse(subdf$log10p_nm > 100, 100, subdf$log10p_nm)

# Across all regions:
topdf = merge(aggregate(lp ~ gene, subdf, mean),
    agg.rename(lp ~ gene + de, subdf, length, 'nct'))
topdf = topdf[topdf$de != 'NS',]

# Plot top ordered by average log10p-value with at least 5 signif):
topdf = topdf[order(topdf$lp, decreasing=T),]
topdf = topdf[order(topdf$nct, decreasing=T),]
topdf = topdf[topdf$nct >= 5,] # At least 5 signif
ug = head(topdf$gene[topdf$de == 'Up'], 8)
dg = head(topdf$gene[topdf$de == 'Down'], 8)
top.genes = c(ug, dg)

# Plot these top genes:
seldf = allresdf[allresdf$gene %in% top.genes,]

# Make matrices for plotting:
cmat = pivot.tomatrix(seldf[,c('gene','sr','logFC_nb')], 'sr', 'logFC_nb')
pmat = pivot.tomatrix(seldf[,c('gene','sr','col_nm')], 'sr', 'col_nm')
cn = expand.grid(keep.sets, regset)
cn = paste0(cn[,1], '@', cn[,2])
cmat = cmat[c(ug, dg), cn]
pmat = pmat[c(ug, dg), cn]
colsplit = sub("_.*@.*","", colnames(cmat))
rowsplit = ifelse(rownames(cmat) %in% ug, 'Up', 'Down')

mx = .5
if (path %in% c('nft','plaq_n', 'plaq_d')){ mx = mx / 20 }
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
ux = 1.5
ht = Heatmap(cmat, 
    use_raster=FALSE, 
    name='logFC_nb',
    col=col_fun,
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    column_split=colsplit,
    row_split=rowsplit,
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    width=(ncol(cmat) + (length(unique(colsplit)) - 1)* 3/4) * unit(ux, 'mm'),
    height=(nrow(cmat) + (length(unique(rowsplit)) - 1)* 3/4) * unit(ux, 'mm'),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        p = pmat[i,j]
        if (!is.na(p) & (p != 0)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.25))
        }
    })

h = .5 + 1 / 15 * nrow(cmat)
w = 1.5 + 1 / 15 * ncol(cmat)
pltprefix = paste0(imgpref, 'regionshared_heatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)




# By region:
regdf = merge(aggregate(lp ~ gene, subdf, mean),
    agg.rename(lp ~ gene + de + region, subdf, length, 'nct'))
regdf = regdf[regdf$de != 'NS',]
regdf = regdf[order(regdf$lp, decreasing=T),]
regdf = regdf[order(regdf$nct, decreasing=T),]

print(set)
suff = paste0(set, '.', path)
subdf = allresdf[allresdf$set == set,]
