#!/usr/bin/R
# ------------------------------------------------
# Plot example DEGs from msex*AD interaction model
# Updated: 11/16/23
# ------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)

library(ggplot2)
library(ggrepel)
library(ggrastr)
library(ggpubr)
options(width=170) 
print(version)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir, srdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


# Functions:
# ----------
load_psbulk_matrix = function(set){
    ststr = sub('_[A-Za-z]+$', '', set)
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
    load(psdata.rda)
    # If Mic, remove T, CAMs from psbulk data:
    if (set == 'Mic_Immune_Mic'){
        ind = !(ps.data$meta$cell_type_high_resolution %in% c('T cells', 'CAMs'))
        ps.data$meta = ps.data$meta[ind,]
        ps.data$mat = ps.data$mat[,rownames(ps.data$meta)]
    }
    # Merge to sample-level (individual x region):
    ps.data = aggregate_psbulk_samplelevel(ps.data)
    return(ps.data)
}


# Load the aggregated results:
# ----------------------------
path = 'msex'
fullaggrda = paste0(regdir, 'allmethods.', path, '.merged.rda')
load(fullaggrda)
sets = names(setdflist)
keep.sets = c('Ast_Ast','Oli_Oli','Opc_Opc',
    'Mic_Immune_Mic','Inh_Inh','Exc_Exc')


# Load ABC scores and add metadata:
# ---------------------------------
indmeta_tsv = 'Annotation/metadata_PFC_all_individuals_092520.tsv'
ext.meta = read.delim(indmeta_tsv, header=T)
ext.meta$kept.ind = ifelse(ext.meta$projid %in% kept.individuals, 'Our Cohort\n(48 individ.)', 'ROSMAP\n(Sept. 2020)')
abcdf = read.delim('Annotation/abc_scores_092023.tsv', header=T)
abcdf = merge(abcdf, ext.meta[,c('projid', 'niareagansc', 'cogdx', 'apoe_genotype', 'msex')])


# Plot the composition of the cohort:
# -----------------------------------
abcdf$cogdx_AD_yes = ifelse(abcdf$niareagansc %in% c(3,4), 'nonAD', ifelse(abcdf$cogdx > 3, 'CI', 'Resil'))
abcdf$has_apoe = ifelse(abcdf$apoe_genotype %in% c(34,44), 'yes','no')
abcdf$nrad = ifelse(abcdf$niareagansc > 2, 'CTRL','AD')
table(abcdf[,c('msex','cogdx')])
table(abcdf[,c('msex','cogdx_AD_yes')])
table(abcdf[,c('msex','apoe_genotype')])
table(abcdf[,c('msex','has_apoe')])

df = as.data.frame(table(abcdf[,c('msex','has_apoe', 'nrad')]))
gp = ggplot(df, aes(nrad, Freq, fill=msex, alpha=has_apoe, label=Freq)) + 
    geom_bar(stat='identity', position='dodge') + 
    geom_text() + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() 
pltprefix = paste0(imgpref, 'msex_cohort_split_barplot')
saveGGplot(gp, pltprefix, w=5, h=3.5)


# Load resources for specific set:
# --------------------------------
set = 'Mic_Immune_Mic'
resdf = c()
for (set in sets){
    print(set)
    setdf = setdflist[[set]]
    peff = colnames(setdf)[grep('^p_.*_nb$', colnames(setdf))]
    pathstr = sub("_nb$", "", sub("^p_","", peff))
    cvars = paste0('col_', pathstr, '_nm')

    # Regional DEGs that are supported in "allregions" and 1+ indpdt region(s):
    for (pstr in pathstr){
        cvar = paste0('col_', pstr, '_nm')
        mat = pivot.tomatrix(setdf[, c('gene','region',cvar)], 'region', cvar)
        nmat = sweep(mat, 1, mat[,'allregions'], '/')
        nmat[is.na(nmat)] = 0
        nmarg = apply(nmat == 1, 1, sum)
        repgenes = names(which(nmarg > 1))
        if (length(repgenes) > 0){
            print(mat[repgenes,])
            resdf = rbind(resdf, data.frame(gene=repgenes, 
                    marg=nmarg[repgenes], set=set, eff=cvar))
        }
    }
}

resdf = resdf[order(resdf$marg, decreasing=T),]
head(resdf, 20)


# Plot the logFC for some of these genes:
# ---------------------------------------
ll = lapply(keep.sets, function(x){
    genes = head(unique(resdf[resdf$set == x,'gene']), 3)
    if (length(genes) > 0){
        # df = data.frame(gene=genes, set=x)
        setdf = setdflist[[x]]
        df = setdf[setdf$gene %in% genes,]
        df$set = x
        return(df)
    } 
    })

seldf = do.call(rbind, ll)
peff = colnames(seldf)[grep('^p_.*_nb$', colnames(seldf))]
pathstr = sub("_nb$", "", sub("^p_","", peff))
cvars = paste0('col_', pathstr, '_nm')
lvars = paste0('logFC_', pathstr, '_nb')
seldf = seldf[,c('gene','region',lvars, cvars, 'set')]

ldf = gather(seldf[,c('gene','region', 'set', lvars)], 'eff','logFC', -gene,  -region, -set)
cdf = gather(seldf[,c('gene','region', 'set', cvars)], 'eff','col', -gene,  -region, -set)
ldf$eff = sub("_nb$", "", sub("^logFC_","", ldf$eff))
cdf$eff = sub("_nm$", "", sub("^col_","", cdf$eff))
seldf = merge(ldf, cdf)

# Matrices:
seldf$gs = with(seldf, paste0(gene, '@', set))
seldf$re = with(seldf, paste0(region, '@', eff))
cmat = pivot.tomatrix(seldf[,c('gs','re','logFC')], 're','logFC')
pmat = pivot.tomatrix(seldf[,c('gs','re','col')], 're','col')

# Colsplit + heatmap:
colsplit = sub(".*@","", colnames(cmat))
rowsplit = sub(".*@","", rownames(cmat))

# Plot matrix:
mx = .5
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
    width=ncol(cmat) * unit(ux, 'mm'),
    height=nrow(cmat) * unit(ux, 'mm'),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        p = pmat[i,j]
        if (!is.na(p) & (p != 0)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.25))
        }
    })

pltprefix = paste0(imgpref, 'example_shared_msex_DEGs')
h = 1.5 + 1 / 15 * nrow(cmat)
w = 1.5 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)

