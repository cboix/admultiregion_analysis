#!/usr/bin/R
# -------------------------------------------------------
# Plot energy metabolism for glia and excitatory neurons:
# Updated: 11/07/23
# ------------------------------------------------------
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

pivot.tomatrix.v2 = function(df, index, key, value){
    require(dplyr)
    df = df[,c(index, key, value)]
    wide = spread(df, key, value)
    mat = as.matrix(wide[, -1])
    rownames(mat) = wide[, 1]
    return(mat)
}


# Load in gene sets:
# ------------------
geneset.rds = paste0(sdbdir, 'energy_metab_genesets.Rds')
if (!file.exists(geneset.rds)){
    # All ETC genes, by complex/function:
    df = read.delim(paste0(sdbdir, 'HGNC-group-639.csv'), header=T, sep=",", skip=1, check.names=FALSE)
    names(df) = c('HGNC_id', 'symbol', 'name','prev_symbol', 'aliases', 'chr', 'group')
    # Make sure all genes match:
    df$gene = df$symbol
    ind = !(df$symbol %in% anno$symbol)
    df$gene[ind] = sapply(df$prev_symbol[ind], function(x){
        x = strsplit(x, ',')[[1]]
        out = NA
        for (v in x){ if (v %in% anno$symbol){ out = v }}
        return(out)
        })
    mt.genes = df$gene[df$chr == 'mitochondria']
    df = df[!(df$gene %in% mt.genes),]
    etc.groups = unique(df$group)
    geneset = lapply(etc.groups, function(x){
        df$gene[df$group == x] })
    names(geneset) = etc.groups
    geneset[['MT-ETC']] = mt.genes

    # Specified:
    geneset[['Lactate']] = c('LDHAL6B','LDHB','LDHA','LDHC','SLC16A1','SLC16A3', 'SLC16A7', 'SLC4A4','CA2')
    geneset[['Glutamine']] = c('GS', 'SLC1A2','SLC1A3', 'SLC38A2','SLC38A1', 'SLC38A3','SLC38A5')

    # From WikiPathways:
    df = read.delim(paste0(sdbdir, 'WP78-datanodes.tsv'), header=T)
    tca.genes = df$Label[df$Type == 'GeneProduct']
    tca.genes[!(tca.genes %in% anno$symbol)]
    geneset[['TCA']] = unique(c('PDHA1','DLAT', 'DLD', 'MDH1','MDH2', 'PC', tca.genes))
    # other: 'GOT1', 'GOT2'

    df = read.delim(paste0(sdbdir, 'WP534-datanodes.tsv'), header=T)
    glyc.genes = df$Label[df$Type == 'GeneProduct']
    glyc.genes = c(glyc.genes, 'PKM')  # is PKM1/2
    glyc.genes = glyc.genes[!(glyc.genes %in% c('PGI', 'PKM1', 'PKM2'))]
    glyc.genes[!(glyc.genes %in% anno$symbol)]

    # Remove non-core glycolysis genes from glyc.genes:
    glyc.genes = glyc.genes[!(glyc.genes %in% geneset[['TCA']])]
    geneset[['Glycolysis']] = glyc.genes[!(glyc.genes %in% c('GOT1','GOT2', 'PCK1', geneset[['Lactate']]))]

    saveRDS(geneset, file=geneset.rds)
} else {
    geneset = readRDS(geneset.rds)
}


# Re-label sets (for plotting:
names(geneset) = sub("Mitochondrial complex ", "ETC-", names(geneset))
names(geneset) = sub("ubiqu.*core.*", "core", names(geneset))
names(geneset) = sub("ubiqu.*super.*", "addtl", names(geneset))
names(geneset) = sub(": .*", "", names(geneset))
# NOTE: Non-unique genes = SDH in both TCA and ETC.
geneset[['TCA']] = geneset[['TCA']][!(geneset[['TCA']] %in% geneset[['ETC-II']])]

ll = lapply(names(geneset), function(x){
    data.frame(gene=geneset[[x]], group=x) })
ldf = do.call(rbind, ll)
gmap = ldf$group
names(gmap) = ldf$gene


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
setnames = c("Mic_Immune_Mic"='Microglial', "Ast_Ast"='Astrocyte', 
    "Opc_Opc"='OPC', "Oli_Oli"='Oligodendrocyte', 
    'Inh_Inh'='Inhibitory neuron','Exc_Exc'='Excitatory neuron')
pathlist = c('plaq_d', 'plaq_n', 'nft', 'nrad', 'cogdxad')
regset = c('allregions', 'EC', 'HC', 'TH', 'AG', 'MT', 'PFC')


# Load all DEG sets:
# -----------------
kept.cols = c('gene','col_nm','path','region', 'logFC_nb', 'p_nb')
alldf = c()
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)
    for (set in keep.sets){
        print(set)
        setdf = setdflist[[set]][, kept.cols]
        setdf$set = set
        alldf = rbind(alldf, setdf)
    }
}


# Plot coefficients of these sets for glia + neurons, by AD variable.
# -------------------------------------------------------------------
metab.df = merge(alldf[alldf$region == 'allregions',], ldf)
metab.df$sp = paste0(metab.df$set, '@', metab.df$path)
metab.df$adjFC = ifelse(metab.df$path %in% c('nft','plaq_n', 'plaq_d'),
    metab.df$logFC_nb * 25, metab.df$logFC_nb)
cmat = pivot.tomatrix(metab.df[,c('gene','sp','adjFC')], 'sp','adjFC')
pmat = pivot.tomatrix(metab.df[,c('gene','sp','col_nm')], 'sp','col_nm')

# Order in same way:
cn = expand.grid(keep.sets, pathlist)
cn = paste0(cn[,1], '@', cn[,2])
nmat = cmat
nmat[is.na(nmat)] = 0
nmat = reord(nmat)
# nmat = t(diag.mat2(t(nmat))[[1]])
cmat = cmat[rownames(nmat),cn]
pmat = pmat[rownames(cmat),cn]
rowsplit = gmap[rownames(cmat)]

# Plot matrix:
mx = .5
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
ux = 1.5
ht = Heatmap(cmat, 
    use_raster=FALSE, 
    name='logFC_nb',
    col=col_fun,
    # column_title=',
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    column_split=sub("_.*","", colnames(cmat)),
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

pltprefix = paste0(imgpref, 'energy_metab_allregions_heatmap')
h = .5 + 1 / 15 * nrow(cmat)
w = 1.5 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)



# Repeat for all regions:
# -----------------------
metab.df = merge(alldf, ldf)
metab.df$sp = with(metab.df, paste0(set, '@', path, ':', region))
metab.df$adjFC = ifelse(metab.df$path %in% c('nft','plaq_n', 'plaq_d'),
    metab.df$logFC_nb * 25, metab.df$logFC_nb)
cmat = pivot.tomatrix(metab.df[,c('gene','sp','adjFC')], 'sp','adjFC')
pmat = pivot.tomatrix(metab.df[,c('gene','sp','col_nm')], 'sp','col_nm')

# Order in same way:
cn = expand.grid(keep.sets, pathlist[1:3], regset)
cn = paste0(cn[,1], '@', cn[,2], ':', cn[,3])
cn = cn[cn %in% colnames(cmat)]
nmat = cmat
nmat[is.na(nmat)] = 0
nmat = reord(nmat)
# nmat = t(diag.mat2(t(nmat))[[1]])
cmat = cmat[rownames(nmat),cn]
pmat = pmat[rownames(cmat),cn]
rowsplit = gmap[rownames(cmat)]
colsplit = sub(":.*","",sub("_.*@","-", colnames(cmat)))
colsplit = sub("-plaq_d","\n1-plaq_d", colsplit)
colsplit = sub("-plaq_n","\n2-plaq_n", colsplit)
colsplit = sub("-nft","\n3-nft", colsplit)

# Plot matrix:
mx = .5
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
ux = 1.5
ht = Heatmap(cmat, 
    use_raster=FALSE, 
    name='logFC_nb',
    col=col_fun,
    # column_title=',
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
    cluster_column_slices=FALSE, 
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

pltprefix = paste0(imgpref, 'energy_metab_regional_heatmap')
h = 1 + 1 / 15 * nrow(cmat)
w = 1.5 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)


# Score at the pseudobulk level, for each sample:
# -----------------------------------------------
geneset.ext = geneset
geneset.ext[['ETC']] = c(geneset[['ETC-II']], geneset[['ETC-III']],
    geneset[['ETC-IV']], geneset[['ETC-V']], geneset[['NADH:core']],  geneset[['NADH:addtl']])

psbulk.scores.rds = paste0(sdbdir, 'energy_metab.psbulk_scores.Rds')
glyc.scores.rds = paste0(sdbdir, 'glycolysis.psbulk_scores.Rds')
if (!file.exists(psbulk.scores.rds)){
    scoredf = c()
    glycdf = c()
    for (set in keep.sets){
        print(set)
        ps.data = load_psbulk_matrix(set)
        # Score each set:
        setdf = c()
        for (key in names(geneset.ext)){
            # TODO: sum > log1p? or log1p > mean
            subgenes = geneset.ext[[key]]
            subgenes = intersect(subgenes, rownames(ps.data$mat))
            if (length(subgenes) > 0){
                submat = log1p(ps.data$mat[subgenes,])
                subscore = apply(submat, 2, mean)
                subdf = data.frame(pr=names(subscore), score=subscore, key=key, set=set)
                setdf = rbind(setdf, subdf)
                # Also save glycolysis separately:
                if (key == 'Glycolysis'){
                    gsdf = data.frame(as.matrix(submat), check.names=FALSE)
                    gsdf$gene = rownames(gsdf)
                    gsdf = gather(gsdf, pr, score, -gene)
                    gsdf$set = set
                }
            }
        }
        gsdf = merge(gsdf, ps.data$meta)
        setdf = merge(setdf, ps.data$meta)
        glycdf = rbind(glycdf, gsdf)
        scoredf = rbind(scoredf, setdf)
    }
    saveRDS(scoredf, file=psbulk.scores.rds)
    saveRDS(glycdf, file=glyc.scores.rds)
} else {
    scoredf = readRDS(psbulk.scores.rds)
    glycdf = readRDS(glyc.scores.rds)
}


# Load ABC scores and add metadata:
# ---------------------------------
indmeta_tsv = 'Annotation/metadata_PFC_all_individuals_092520.tsv'
ext.meta = read.delim(indmeta_tsv, header=T)
ext.meta$kept.ind = ifelse(ext.meta$projid %in% kept.individuals, 'Our Cohort\n(48 individ.)', 'ROSMAP\n(Sept. 2020)')
abcdf = read.delim('Annotation/abc_scores_092023.tsv', header=T)
abcdf = merge(abcdf, ext.meta[,c('projid', 'niareagansc', 'cogdx')])

# Add metadata variables:
metadata$braak.group = ifelse(metadata$braaksc %in% c(0,1,2), '0-2', ifelse(metadata$braaksc %in% c(3,4), '3-4', '5-6'))
meta.cols = c('projid', 'region', 'rind', 'braaksc', 'nrad','niareagansc','cogdx', 'Apoe_e4', 'braak.group')
scoredf = merge(scoredf, unique(metadata[,meta.cols]))
scoredf = merge(scoredf, abcdf)
glycdf = merge(glycdf, unique(metadata[,meta.cols]))
glycdf = merge(glycdf, abcdf)

# NOTE: Could weight mean by set size?

# Plot means + SE by Braak group + region:
# ----------------------------------------
use.apoe = FALSE
ncell.cutoff = 10  # At least remove small sets

var = 'nia_aa_sc'
suff = var
formvec = c('score ~ ', var, '+ region + set + key')
if (use.apoe){
    suff = paste0(suff, '-e4')
    formvec = c(formvec, '+Apoe_e4')
}
# var = 'braak.group'
form = asform(formvec)
mn.scoredf = scoredf[scoredf$ncell >= ncell.cutoff,]
mdf = merge(merge(agg.rename(form, mn.scoredf, mean, 'mean'),
    agg.rename(form, mn.scoredf, sd, 'sd')),
    agg.rename(form, mn.scoredf, length, 'nind'))
mdf$se = mdf$sd/ sqrt(mdf$nind)
mdf$var = mdf[[var]]
mdf$fvar = factor(mdf$var)


set = 'Ast_Ast'
for (set in keep.sets){
    if (use.apoe){
        gp = ggplot(mdf[mdf$set == set,], aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region, linetype=Apoe_e4))
    } else {
        gp = ggplot(mdf[mdf$set == set,], aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region))
    }
    gp = gp + facet_wrap(~ key, scales='free') + 
        geom_point() + 
        labs(x=var, title=paste0(set, ' (',var, ')')) + 
        scale_color_manual(values=reg.cols) +
        geom_errorbar(width=.25) + 
        theme_pubr() + theme(legend.position='none')
    if (class(mdf$var) != 'character'){
        gp = gp + geom_line()
    }
    pltprefix = paste0(imgpref, 'energy_genesets_profiles.', suff, '.', set)
    saveGGplot(gp, pltprefix, w=9, h=6.5)
}


key = 'Glycolysis'
subsets = c('Ast_Ast', 'Mic_Immune_Mic', 'Exc_Exc')
sub.mdf = mdf[(mdf$key == key) & (mdf$set %in% subsets),]
sub.mdf$region = factor(sub.mdf$region, levels=regset[regset %in% sub.mdf$region])
if (use.apoe){
    gp = ggplot(sub.mdf, aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region, linetype=Apoe_e4))
} else {
    gp = ggplot(sub.mdf, aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region))
}
gp = gp + facet_grid(set ~ region, scales='free_y') + 
    geom_point() + 
    labs(x=var, title=paste0(set, ' (',var, ')')) + 
    scale_color_manual(values=reg.cols) +
    geom_errorbar(width=.25) + 
    theme_pubr() + theme(legend.position='none')
if (class(mdf$var) != 'character'){
    gp = gp + geom_line()
}
pltprefix = paste0(imgpref, 'energy_genesets_profiles.', suff, '.', key)
saveGGplot(gp, pltprefix, w=7, h=4)


# Against nft, plaque:
msdf = merge(mn.scoredf, pqdf, all.x=TRUE)
set = 'Ast_Ast'
for (set in keep.sets){

    gp = ggplot(msdf[msdf$set == set,], aes(log1p(plaq_n), score, color=region)) + 
        facet_wrap(~ key, scales='free') + 
        geom_point() + 
        geom_smooth(method='lm', alpha=.1) +
        labs(x=var, title=paste0(set, ' (',var, ')')) + 
        scale_color_manual(values=reg.cols) +
        theme_pubr() + theme(legend.position='none')

    pltprefix = paste0(imgpref, 'energy_genesets_profiles.', suff, '.', set)
    saveGGplot(gp, pltprefix, w=9, h=6.5)
}





key = 'Glycolysis'
subsets = c('Ast_Ast', 'Mic_Immune_Mic', 'Exc_Exc')
sub.mdf = mdf[(mdf$key == key) & (mdf$set %in% subsets),]
sub.mdf$region = factor(sub.mdf$region, levels=regset[regset %in% sub.mdf$region])
if (use.apoe){
    gp = ggplot(sub.mdf, aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region, linetype=Apoe_e4))
} else {
    gp = ggplot(sub.mdf, aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region))
}
gp = gp + facet_grid(set ~ region, scales='free_y') + 
    geom_point() + 
    labs(x=var, title=paste0(set, ' (',var, ')')) + 
    scale_color_manual(values=reg.cols) +
    geom_errorbar(width=.25) + 
    theme_pubr() + theme(legend.position='none')
if (class(mdf$var) != 'character'){
    gp = gp + geom_line()
}
pltprefix = paste0(imgpref, 'energy_genesets_profiles.', suff, '.', key)
saveGGplot(gp, pltprefix, w=7, h=4)



# Separate out glycolysis genes:
# ------------------------------
key = 'Glycolysis'
suff = var
formvec = c('score ~ ', var, '+ region + set + gene')
# var = 'braak.group'
form = asform(formvec)
mn.glycdf = glycdf[glycdf$ncell >= ncell.cutoff,]
mdf = merge(merge(agg.rename(form, mn.glycdf, mean, 'mean'),
    agg.rename(form, mn.glycdf, sd, 'sd')),
    agg.rename(form, mn.glycdf, length, 'nind'))
mdf$se = mdf$sd/ sqrt(mdf$nind)
mdf$var = mdf[[var]]
mdf$fvar = factor(mdf$var)

sub.mdf = mdf[(mdf$set %in% subsets),]
sub.mdf$region = factor(sub.mdf$region, levels=regset[regset %in% sub.mdf$region])
gp = ggplot(sub.mdf, aes(var, mean, ymin=mean - 2 * se, ymax=mean + 2 * se, color=region))
gp = gp + facet_grid(gene ~ set, scales='free_y') + 
    geom_point() + 
    labs(x=var, title=paste0(set, ' (',var, ')')) + 
    scale_color_manual(values=reg.cols) +
    geom_errorbar(width=.25) + 
    theme_pubr() + theme(legend.position='none')
if (class(mdf$var) != 'character'){
    gp = gp + geom_line(lwd=1.25)
}
pltprefix = paste0(imgpref, 'glycolysis_genes_profiles.', suff, '.', key)
saveGGplot(gp, pltprefix, w=3.5, h=30)



# Compare Ast, Mic to Exc:
# ------------------------
keys = c('Glycolysis', 'TCA', 'ETC', 'MT-ETC', 'Lactate')
scoredf$sk = paste0(scoredf$set, '@', scoredf$key)
sub.scoredf = scoredf[(scoredf$set %in% subsets) & (scoredf$key %in% keys),]
# sub.scoredf = sub.scoredf[sub.scoredf$ncell >= ncell.cutoff,]
sub.scoredf = sub.scoredf[sub.scoredf$ncell >= 25,]
smat = pivot.tomatrix.v2(sub.scoredf, 'pr','sk', 'score')

cn = colnames(smat)
cmat = cor(smat, use='pairwise.complete.obs')
pmat = cmat * 0 + 1.0
# TODO: cut in half
for (c1 in colnames(smat)){
    s1 = smat[,c1]
    for (c2 in colnames(smat)){
        s2 = smat[,c2]
        ind = !(is.na(s1)) & !(is.na(s2))
        ct = cor.test(s1[ind], s2[ind])
        pmat[c1, c2] = ct$p.value
    }
}

rowsplit = sub("@.*", "", rownames(cmat))
colsplit = sub("@.*", "", colnames(cmat))

mx = 1
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
ux = 1.5
ht = Heatmap(cmat, 
    use_raster=FALSE, 
    name='Correlation',
    col=col_fun,
    cluster_columns=FALSE, 
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    cluster_column_slices=FALSE,
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
        if (!is.na(p) & (p < 0.001)){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.25))
        }
    }
)

pltprefix = paste0(imgpref, 'energy_metab.psbulk_scores.corr_heatmap')
h = 1.5 + 1 / 15 * nrow(cmat)
w = 2 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)

