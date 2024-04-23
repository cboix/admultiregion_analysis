#!/usr/bin/R
# --------------------------------------------------
# Part of metabolic modules figures
# - Systematic scan for e4-associated modules
# Updated 06/20/2023
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

# Filter small modules:
ngene = 10
ngdf = unique(runscdf[,c('rm', 'ng')])
mods = sort(unique(ngdf$rm[ngdf$ng >= ngene]))


# Individual-level scores versus e4 status:
# -----------------------------------------
runscdf$totscore = with(runscdf, score * ncell)
runscdf$projid = sub("-.*", "", runscdf$pr)
aggdf = aggregate(cbind(totscore, ncell) ~ projid + runset + module + rm, runscdf, sum)
aggdf$score = aggdf$totscore / aggdf$ncell 
aggdf = merge(aggdf, unique(metadata[,c('projid','Apoe_e4','cogdxad','msex', 'nrad')]))

# Assoc of e4 with module score:
resdf = NULL
for (mod in mods){
    fit = lm(score ~ Apoe_e4, aggdf[aggdf$rm == mod,])
    cfit = coefficients(summary(fit))
    cfit = data.frame(cfit)[2,, drop=FALSE]
    names(cfit) = c('Est','SE','t','p')
    cfit$rm = mod
    resdf = rbind(resdf, cfit)
}
resdf = resdf[order(resdf$p),]
head(resdf[(resdf$p < 0.05) & (resdf$Est > 0),], 20)


# Region-level scores versus e4 status:
# -------------------------------------
runscdf = merge(runscdf, unique(metadata[,c('projid','Apoe_e4','cogdxad','msex', 'nrad')]))
runscdf$region = sub(".*-","", runscdf$pr)

# Assoc of e4 with module score:
resdf = NULL
for (mod in mods){
    # fit = lm(score ~ Apoe_e4 * region, runscdf[runscdf$rm == mod,])
    fit = lm(score ~ Apoe_e4, runscdf[runscdf$rm == mod,])
    cfit = coefficients(summary(fit))
    cfit = data.frame(cfit)[2,, drop=FALSE]
    names(cfit) = c('Est','SE','t','p')
    cfit$rm = mod
    resdf = rbind(resdf, cfit)
}
resdf = resdf[order(resdf$p),]
head(resdf[(resdf$p < 0.05) & (resdf$Est > 0),], 20)

topdf = head(resdf[(resdf$p < 0.05) & (resdf$Est > 0),], 10)
topdf$rm = factor(topdf$rm, levels=rev(topdf$rm))

# Plot these scores
gp = ggplot(topdf, aes(rm, -log10(p))) + #, color=Apoe_e4)) + 
    geom_bar(fill='slateblue', stat='identity') + 
    theme_pubr() + coord_flip()
pltprefix = paste0(imgpref, 'barplot_lm_e4')
saveGGplot(gp, pltprefix, w=2, h=2)



# Plot gene enrichments for allgenes + coregenes:
# - Opc-9 and Oli-7
# -----------------------------------------------
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))
# sets = with(setdf[setdf$set != 'mt',], paste0(runset, '-', module))
sets = c('Opc-9', 'Oli-7')

# Load DEGs to add top DEGs for each enrichment:
fullaggrda = paste0(regdir, 'allmethods.allmajor.merged.rda')
load(fullaggrda)

setmap = c("Mic_Immune"="Mic_Immune_Mic", "Ast"="Ast_Ast", 
    "Opc"="Opc_Opc", "Oli"="Oli_Oli")
path = 'cogdxad'
# path = 'plaq_n'

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
    degs = dedf$gene[dedf$col_nm != 0]
    gpdf$topgenes = sapply(gpdf$intersection, function(x, ntop=5){
        x = strsplit(x, ',')[[1]]
        x = head(degs[degs %in% x], ntop)
        return(paste0(x, collapse=','))
        })
    # Mapping term ids to genes:
    mapgenes = gpdf$topgenes
    names(mapgenes) = gpdf$term_id

    gpdf = gp2.result$result
    gpdf$ngint = sapply(gpdf$intersection, function(x){
       length(strsplit(x, ',')[[1]]) })
    gpdf = gpdf[order(gpdf$p_value),]
    gpdf$nc = nchar(gpdf$term_name)
    gpdf = gpdf[gpdf$nc < 40,]
    gpdf = gpdf[gpdf$ngint > 2,]
    gpdf = gpdf[gpdf$term_size < 500,]
    # gpdf = pruneWithInt(gpdf, gpdf, cutoff=0.75)
    gpdf$intersection = mapgenes[gpdf$term_id]
    pltprefix = paste0(imgpref, 'enrbarplot_lt500_', set)
    gp = plotGObarplot(gpdf, pltprefix, ntop=10)
}

