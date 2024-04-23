#!/usr/bin/R
# --------------------------------------------------
# Part of metabolic modules figures
# - GO enrichments
# - Distr vs. regions/path
# Updated 06/19/2023
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


library(tidyr)
library(viridis)
library(gprofiler2)

library(ggpubr)
library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
resdir = paste0(sdbdir, 'modules/resources/')
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'module_metab_')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
system(cmd)


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))

getGeneSet <- function(module, uselist=gmlist){
    rs = sub("-.*", "", module)
    num = as.numeric(sub(".*-", "", module))
    coremap = uselist[[rs]]
    x = names(coremap)[coremap == num]
    return(x)
}

# Load DE results:
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
load(fullaggrda)


# Get modules that match specific sets:
# -------------------------------------
genesets = list()
genesets[['gly']] = c('PDK1','PFKL','PFKP','LDHA','VEGFA','DDIT4')
genesets[['mt']] = c('MT-ND3','MT-CO3','MT-CYB','MT-ND4')
genesets[['ins']] = c('HIF3A','FKBP5','FOXG1', 'PIK3R1','INSR','FOXO1')
genesets[['oli']] = c('SLC38A2','MID1IP1','SPP1', 'SLC39A11')

uselist = gmlist
dflist = lapply(names(uselist), function(x){
    y = uselist[[x]]
    df = data.frame(gene=names(y), module=y, runset=x)
    return(df)
})

df = do.call(rbind, dflist)

setdf = c()
for (set in names(genesets)){
    gs = genesets[[set]]
    ng = length(gs)
    df$in.set = df$gene %in% gs
    sdf = aggregate(in.set ~ runset + module, df, sum)
    print(head(sdf[order(sdf$in.set, decreasing=T),]))
    sdf = sdf[sdf$in.set == ng,]
    sdf$set = set
    setdf = rbind(setdf, sdf)
}
setdf = setdf[!(setdf$runset %in% names(exc.sets)),]


# Reduce to overarching individual-level scores:
# ----------------------------------------------
runscdf$totscore = with(runscdf, score * ncell)
runscdf$projid = sub("-.*", "", runscdf$pr)
aggdf = aggregate(cbind(totscore, ncell) ~ projid + runset + module + rm, runscdf, sum)
aggdf$score = aggdf$totscore / aggdf$ncell 
aggdf = merge(aggdf, setdf)
aggdf$rms = with(aggdf, paste0(rm, '@', set))
aggdf = merge(aggdf, unique(metadata[,c('projid','Apoe_e4','cogdxad','msex', 'nrad')]))

scdf = merge(runscdf, setdf)
scdf$rms = with(scdf, paste0(rm, '@', set))
scdf = merge(scdf, unique(metadata[,c('projid','Apoe_e4','cogdxad','msex', 'nrad')]))
scdf$region = sub(".*-", "", scdf$pr)


# Sets for plotting:
setcts = c('Ast','Mic_Immune','Opc')
sets = with(setdf[(setdf$set %in% c('gly','mt')) &
    (setdf$runset %in% setcts),], paste0(runset, '-', module))


# Plot boxplots vs. AD + region:
# ------------------------------
gp = ggplot(scdf[scdf$rm %in% sets,], aes(region, log1p(score), fill=nrad)) + 
    facet_wrap(~rm, scales='free', nrow=2) + 
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.25, dodge.width=.8), cex=.75) +
    scale_fill_manual(values=c('CTRL'='grey80', 'AD'='red')) + 
    stat_compare_means(label='p.format') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'boxplots_ad-region')
saveGGplot(gp, pltprefix, w=8, h=8)

# TODO: BOXPLOTS VS STAGE/PROGRESSION


# Plot scores as heatmap for one cell type:
# -----------------------------------------
matlist = list() # Store all
setmap = c("Mic_Immune"="Mic_Immune_Mic", "Ast"="Ast_Ast", 
    "Opc"="Opc_Opc", "Oli"="Oli_Oli")

ct = 'Ast'
# ct = 'Mic_Immune'
use.de = TRUE
mlist = unique(scdf$rm[scdf$runset == ct])
if (use.de){
    deset = setmap[ct]
    dedf = setdflist[[deset]]
    # Keep DE in any cross-region region/indiv-level DEG set
    dedf = dedf[dedf$col_nm != 0,]
    # Load pseudobulk data for scoring:
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ct, '.rda')
    load(psdata.rda)
    tdf = agg.rename(ncell ~ projid + region, ps.data$meta, sum, 'totcell')
    ps.data$meta = merge(ps.data$meta, tdf)
    ps.data$mat = ps.data$mat[, ps.data$meta$ptype]
    ps.data$meta$cellfrac = with(ps.data$meta, ncell / totcell)
}

htlist = list()
ht = NULL
for (module in mlist){
    if (use.de){
        genes = getGeneSet(module)
        aggde = aggregate(logFC_nb ~ gene, dedf[dedf$gene %in% genes,], mean)
        dir = abs(aggde$logFC_nb) # Or just -1,1
        # Weighted average score over cell subtypes:
        ps.data$meta$score = colSums(sweep(ps.data$mat[aggde$gene,], 1, dir, '*'))
        ps.data$meta$score = ps.data$meta$score * ps.data$meta$cellfrac
        df = aggregate(score ~ region + projid, ps.data$meta, sum)
        pmat = pivot.tomatrix(df, 'projid','score')
    } else {
        ind = scdf$rm == module
        pmat = pivot.tomatrix(scdf[ind,c('region','projid','score')], 'projid','score')
    }
    pmat = pmat[reg.order[-1],]
    matlist[[module]] = pmat

    ux = 1.5
    # Heatmap:
    pltmat = log1p(pmat)
    # pltmat = sweep(pltmat, 1, apply(pltmat, 1, max, na.rm=T), '/')
    # Annotation:
    if (module == mlist[1]){
        ha = make_ind_annotation(colnames(pltmat), ux=ux)
    } else { ha = NULL }
    plt = Heatmap(
        pltmat, 
        name=module,
        width=ncol(pmat) * unit(ux, 'mm'),
        height=nrow(pmat) * unit(ux, 'mm'),
        cluster_rows=FALSE,
        use_raster=FALSE,
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp=gpar(color='black', lwd=.5),
        heatmap_legend_param = list(
            legend_height = unit(1.5, "cm"),
            border='black'),
        top_annotation=ha,
        col=viridis(50)
    )
    htlist[[module]] = plt
    ht = ht %v% plt
}

plttype = 'heatmap_indregion_'
if (use.de){ plttype = paste0(plttype, 'deonly_') }
pltprefix = paste0(imgpref, plttype, ct)
w = 3 + ncol(pmat) / 15
h = 3 + nrow(pmat) / 15
saveHeatmap(ht, pltprefix, w=w, h=h)


# Get all module DE scores:
# -------------------------
dscdf = NULL
for (module in sets){
    ct = sub('-.*','', module)
    # Keep DE in any cross-region region/indiv-level DEG set
    deset = setmap[ct]
    dedf = setdflist[[deset]]
    dedf = dedf[dedf$col_nm != 0,]
    # Load pseudobulk data for scoring:
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ct, '.rda')
    load(psdata.rda)
    tdf = agg.rename(ncell ~ projid + region, ps.data$meta, sum, 'totcell')
    ps.data$meta = merge(ps.data$meta, tdf)
    ps.data$mat = ps.data$mat[, ps.data$meta$ptype]
    ps.data$meta$cellfrac = with(ps.data$meta, ncell / totcell)

    # Score with DE genes in module:
    genes = getGeneSet(module)
    aggde = aggregate(logFC_nb ~ gene, dedf[dedf$gene %in% genes,], mean)
    dir = abs(aggde$logFC_nb) # Or just -1,1
    # Weighted average score over cell subtypes:
    ps.data$meta$score = colSums(sweep(ps.data$mat[aggde$gene,], 1, dir, '*'))
    ps.data$meta$score = ps.data$meta$score * ps.data$meta$cellfrac
    df = aggregate(score ~ region + projid, ps.data$meta, sum)
    df$module = module
    df$ct = ct
    dscdf = rbind(dscdf, df)
}

# Merge with pathology data:
dscdf = merge(dscdf, unique(metadata[,c('rind', 'projid','region', 'nrad','cogdxad','braaksc')]))
dscdf = merge(dscdf, pqdf, all.x=TRUE)


# Plot scatters vs. pathology:
# ----------------------------
gp = ggplot(dscdf[dscdf$region != 'TH',], aes(log1p(plaq_d), log1p(score))) + 
    facet_grid(module~region, scales='free') + 
    geom_point(cex=.25) + 
    geom_smooth(method='lm') + 
    stat_cor() + 
    theme_pubr()
pltprefix = paste0(imgpref, 'scatter.plaq_d')
saveGGplot(gp, pltprefix, w=8, h=9)

gp = ggplot(dscdf[dscdf$region != 'TH',], aes(plaq_n + plaq_d, log1p(score))) + 
    facet_grid(module~region, scales='free') + 
    geom_point(cex=.25) + 
    geom_smooth(method='lm') + 
    stat_cor() + 
    theme_pubr()
pltprefix = paste0(imgpref, 'scatter.plaq')
saveGGplot(gp, pltprefix, w=8, h=9)

gp = ggplot(dscdf[dscdf$region != 'TH',], aes(nft, log1p(score))) + 
    facet_grid(module~region, scales='free') + 
    geom_point(cex=.25) + 
    geom_smooth(method='lm') + 
    stat_cor() + 
    theme_pubr()
pltprefix = paste0(imgpref, 'scatter.nft')
saveGGplot(gp, pltprefix, w=8, h=9)


# Make corr. heatmap:
# -------------------
pathlist = c('plaq_d','plaq_n','nft')
regs = reg.order[-1]
regs = regs[-length(regs)]
mx = 0.5
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
htfull = NULL
# Order sets:
sets = sort(sets)
num = as.numeric(sub(".*-","", sets))
sets = sets[order(num > 20)]
for (path in pathlist){
    crmat = matrix(0, nrow=length(sets), ncol=length(regs), dimnames=list(sets, regs))
    cpmat = crmat
    for (module in sets){
        for (region in regs){
            ind = (dscdf$module == module) & (dscdf$region == region)
            x = log1p(dscdf[ind, 'score'])
            ct = cor.test(x, dscdf[ind, path], use='complete')
            crmat[module, region] = ct$estimate
            cpmat[module, region] = ct$p.value
        }
    }

    ux = 1.5
    plt = Heatmap(crmat,
        use_raster=FALSE,
        name='Corr.',
        col=col_fun,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        width = ncol(crmat)*unit(ux, "mm"), 
        height = nrow(crmat)*unit(ux, "mm"),
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lwd=.5),
        cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
            p = cpmat[i,j]
            if (!is.na(p)){
                ann = ifelse(p < 0.1, ifelse(p < 0.05, '*', '.'),'')
                grid.text(ann, x, y,gp=gpar(fontsize=gridtxt.fs))
            }
        })
    htfull = htfull + plt

    pltprefix = paste0(imgpref, 'path_corr_heatmap_', path)
    h = 1 + 1 / 15 * nrow(crmat)
    w = 1 + 1 / 15 * ncol(crmat) * 1.5
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

# NOTE: pvals don't annotate properly?
pltprefix = paste0(imgpref, 'path_corr_heatmap_all')
h = 1 + 1 / 15 * nrow(crmat)
w = 1 + 1 / 15 * ncol(crmat) * 1.5
saveHeatmap(htfull, pltprefix, w=w, h=h)


# GO terms across all of the different glycolysis + MT terms:
# -----------------------------------------------------------
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))

setcts = c('Ast','Mic_Immune','Opc')
sets = with(setdf[(setdf$set %in% c('gly','mt')) &
    (setdf$runset %in% setcts),], paste0(runset, '-', module))

genes = sapply(sets, getGeneSet)
gp2.result = gprofiler2::gost(genes, organism='hsapiens',
    ordered_query=FALSE, multi_query=TRUE,
    sources = sources) 
gpdf = gp2.result$result


# Plot as-is, top functional terms:
# ---------------------------------
gpdf = gp2.result$result
gpdf$nc = nchar(gpdf$term_name)
gpdf = gpdf[gpdf$nc < 40,]
gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
pltprefix = paste0(imgpref, 'enrheatmap_glyMT', '_src')
subpmat = gpPvalMatrix(gpdf, genes, ntop=10)
plt = plotGpPvalMatrix(subpmat, pltprefix)

gpdf = gp2.result$result
gpdf$nc = nchar(gpdf$term_name)
gpdf = gpdf[gpdf$nc < 40,]
gpdf = gpdf[gpdf$term_size < 500,]
pltprefix = paste0(imgpref, 'enrheatmap_glyMT', '_small')
subpmat = gpPvalMatrix(gpdf, genes, ntop=4)
plt = plotGpPvalMatrix(subpmat, pltprefix, use_raster=FALSE)

