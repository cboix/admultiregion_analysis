#!/usr/bin/R
# ------------------------------------
# Plot specific pairs of correlations:
# - score using DEGs
# Updated 06/20/2023
# ------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

# Settings for plots:
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))

library(tidyr)
library(gprofiler2)

library(viridis)
library(ggpubr)
library(ggrastr)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
resdir = paste0(sdbdir, 'modules/resources/')
plotdir = paste0(imgdir, 'modules/')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
system(cmd)




# Functions + variables:
# ----------------------
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))
col_log10p = colorRamp2(c(0, 10), c('white','black'))
colg = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(50)
pcols = brewer.pal(12,'Paired')


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))

getGeneSet <- function(rs, num, uselist=gmlist){
    usemap = uselist[[rs]]
    x = names(usemap)[usemap == num]
    return(x)
}

# Load DE results:
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
load(fullaggrda)


# Set cell type and prep pseudobulk data:
# ---------------------------------------
runset = 'Mic_Immune'
imgpref = paste0(plotdir, 'module_corrplots_', runset,'_')
# subtype.cols = tcols[sort(unique(scdf$cell_type_high_resolution))]
setmap = c("Mic_Immune"="Mic_Immune_Mic", "Ast"="Ast_Ast", 
    "Opc"="Opc_Opc", "Oli"="Oli_Oli")

# Keep DE in any cross-region region/indiv-level DEG set
deset = setmap[runset]
dedf = setdflist[[deset]]
dedf = dedf[dedf$col_nm != 0,]

# Load pseudobulk data for scoring:
psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
load(psdata.rda)

if (runset == 'Mic_Immune'){
    ps.data$meta = ps.data$meta[ps.data$meta$cell_type_high_resolution != 'T cells',]
}
tdf = agg.rename(ncell ~ projid + region, ps.data$meta, sum, 'totcell')
ps.data$meta = merge(ps.data$meta, tdf)
ps.data$mat = ps.data$mat[, ps.data$meta$ptype]
ps.data$meta$cellfrac = with(ps.data$meta, ncell / totcell)


# Set the comparisons of interest:
# --------------------------------
if (runset == 'Ast'){
    mod.comp = list(c(3, 17), c(0, 12), c(6, 27), c(8,15), c(6,13))
} else if (runset == 'Mic_Immune'){
    mod.comp = list(c(1, 20), c(1, 11), c(1, 15),
        c(1, 2), c(11, 2), c(11, 25))
} else {
    mod.comp = list(c(0, 1))
}


nums = unique(unlist(mod.comp))
mdf = NULL
for (num in nums){
    genes = getGeneSet(runset, num)
    aggde = aggregate(logFC_nb ~ gene, dedf[dedf$gene %in% genes,], mean)
    dir = abs(aggde$logFC_nb) # Or just -1,1
    # dir = 1
    # Weighted average score over cell subtypes:
    ps.data$meta$score = colSums(sweep(ps.data$mat[aggde$gene,], 1, dir, '*'))
    ps.data$meta$score = ps.data$meta$score * ps.data$meta$cellfrac
    df = aggregate(score ~ region + projid, ps.data$meta, sum)
    df$score = log1p(df$score)
    df$mod = paste0('M', num)
    mdf = rbind(mdf, df)
}


pltlist = list()
for (i in 1:length(mod.comp)){
    numpair = mod.comp[[i]]
    pair = paste0('M',numpair)
    pstr = paste(pair, collapse='_')
    subdf = mdf[mdf$mod %in% pair,]
    subdf = spread(subdf, mod, score)
    subdf = merge(subdf, unique(metadata[,c('projid','cogdxad')]))

    attr = 'cog'
    if (attr == 'cog'){
        use.attr = 'cogdxad'
        use.cols = c('AD'='red', 'CTRL'='grey75')
    } else {
        use.attr = 'cell_type_high_resolution'
        use.cols = subtype.cols
    }

    gp = ggplot(subdf, aes_string(pair[1], pair[2], color=use.attr)) + 
        rasterize(geom_point(cex=.25), dpi=450) +
        geom_smooth(method='lm', color='black',lwd=.5) + 
        scale_color_manual(values=use.cols) + 
        stat_cor(color='black', cex=3, output.type='text', label.sep='\n', label.y.npc=.95) + 
        theme_pubr() + theme(legend.position='none')
    # gp2 = rasterize(gp, layers='Point', dpi=450)
    pltlist[[i]] = gp
}
garr = ggarrange(plotlist=pltlist)
pltprefix = paste0(imgpref, 'corrplot_', attr)
saveGGplot(garr, pltprefix, w=2*3, h=2*2)


mwide = spread(mdf, mod, score)
mmat = mwide[,3:ncol(mwide)]

