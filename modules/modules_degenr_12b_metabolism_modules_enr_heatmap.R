#!/usr/bin/R
# --------------------------------------------------
# Second part of metabolic modules figures
# - Enrichments for modules
# - Module scores heatmap
# - Boxplots vs. e4
# Updated 04/04/2022
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
library(PRROC)
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


# Plot boxplots vs. e4, cognition:
# --------------------------------
e4cols = c('no'='grey80', 'yes'='slateblue')
gp = ggplot(aggdf[aggdf$set != 'mt',], aes(rm, log1p(score), fill=Apoe_e4)) + #, color=Apoe_e4)) + 
    # facet_wrap(~set, scales='free', nrow=1) + 
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.25, dodge.width=.8), cex=.75) +
    scale_fill_manual(values=e4cols) + 
    scale_color_manual(values=e4cols) + 
    stat_compare_means(label='p.format') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'boxplots_e4')
saveGGplot(gp, pltprefix, w=4, h=4)

gp = ggplot(scdf, aes(rm, log1p(score), fill=Apoe_e4)) + 
    facet_wrap(~set, scales='free', nrow=1) + 
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.25, dodge.width=.8), cex=1) +
    scale_fill_manual(values=c('no'='grey80', 'yes'='slateblue')) + 
    stat_compare_means(label='p.format') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'boxplots_e4_sc')
saveGGplot(gp, pltprefix, w=12, h=4)

gp = ggplot(aggdf, aes(rm, score, fill=cogdxad)) + 
    facet_wrap(~set, scales='free', nrow=1) + 
    geom_boxplot() + 
    scale_fill_manual(values=colvals[['cogdxad']]) + 
    stat_compare_means(label='p.signif') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'boxplots_cogdxad')
saveGGplot(gp, pltprefix, w=12, h=4)

gp = ggplot(scdf, aes(rm, score, fill=cogdxad)) + 
    facet_wrap(~set, scales='free', nrow=1) + 
    geom_boxplot() + 
    scale_fill_manual(values=colvals[['cogdxad']]) + 
    stat_compare_means(label='p.signif') + 
    theme_pubr()
pltprefix = paste0(imgpref, 'boxplots_cogdxad_sc')
saveGGplot(gp, pltprefix, w=12, h=4)


# Plot scores as heatmap:
# -----------------------
# ind = aggdf$set != 'mt'
ind = 1:nrow(aggdf)
pmat = pivot.tomatrix(aggdf[ind,c('rms','projid','score')], 'projid','score')
rsplit = sub(".*@","",rownames(pmat))

# Annotation:
ux = 1.5
ha = make_ind_annotation(colnames(pltmat), ux=ux)

# Heatmap:
pltmat = log1p(pmat)
pltmat = sweep(pltmat, 1, apply(pltmat, 1, max), '/')
plt = Heatmap(
    pltmat, 
    row_split=rsplit,
    width=ncol(pmat) * unit(ux, 'mm'),
    height=nrow(pmat) * unit(ux, 'mm'),
    use_raster=FALSE,
    top_annotation=ha,
    col=viridis(50)
    )

pltprefix = paste0(imgpref, 'heatmap_indscores')
w = 3 + ncol(pmat) / 15
h = 3 + nrow(pmat) / 15
saveHeatmap(plt, pltprefix, w=w, h=h)



# Plot gene enrichments for allgenes + coregenes:
# -----------------------------------------------
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))
sets = with(setdf[setdf$set != 'mt',], paste0(runset, '-', module))

# Load DEGs to add top DEGs for each enrichment:
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
load(fullaggrda)

setmap = c("Mic_Immune"="Mic_Immune_Mic", "Ast"="Ast_Ast", 
    "Opc"="Opc_Opc", "Oli"="Oli_Oli")
path = 'cogdxad'

for (set in sets){
    deset = setmap[sub("-.*", "", set)]
    dedf = setdflist[[deset]]
    dedf = dedf[dedf$path == path,]

    genes = getGeneSet(set)
    gp2.result = gprofiler2::gost(genes, organism='hsapiens',
        ordered_query=FALSE, multi_query=FALSE,
        sources = sources, evcodes=TRUE)

    gpdf = gp2.result$result

    # Get the top N genes per term:
    degs = dedf$gene
    gpdf$topgenes = sapply(gpdf$intersection, function(x, ntop=5){
        x = strsplit(x, ',')[[1]]
        x = head(degs[degs %in% x], ntop)
        return(paste0(x, collapse=','))
        })
    # Mapping term ids to genes:
    mapgenes = gpdf$topgenes
    names(mapgenes) = gpdf$term_id

    # Medium size (<1000)
    gpdf = gp2.result$result
    gpdf = gpdf[order(gpdf$p_value),]
    gpdf = gpdf[gpdf$term_size < 1000,]
    gpdf = pruneWithInt(gpdf, gpdf)
    gpdf$intersection = mapgenes[gpdf$term_id]
    pltprefix = paste0(imgpref, 'enrbarplot_lt1000_', set)
    gp = plotGObarplot(gpdf, pltprefix, ntop=10)

    # gpdf = gp2.result$result
    # gpdf = gpdf[order(gpdf$p_value),]
    # gpdf = gpdf[gpdf$source %in% c("REAC","WP","KEGG","CORUM"),]
    # gpdf = pruneWithInt(gpdf, gpdf)
    # gpdf$intersection = mapgenes[gpdf$term_id]
    # pltprefix = paste0(imgpref, 'enrbarplot_subsource_', set)
    # gp = plotGObarplot(gpdf, pltprefix, ntop=10)

    gpdf = gp2.result$result
    gpdf = gpdf[order(gpdf$p_value),]
    gpdf = gpdf[gpdf$term_size < 500,]
    gpdf = pruneWithInt(gpdf, gpdf, cutoff=0.75)
    gpdf$intersection = mapgenes[gpdf$term_id]
    pltprefix = paste0(imgpref, 'enrbarplot_lt500_', set)
    gp = plotGObarplot(gpdf, pltprefix, ntop=10)

}

