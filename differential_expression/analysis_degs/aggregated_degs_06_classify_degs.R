#!/usr/bin/R
# ------------------------------------------------
# Classify DEGs as individual vs. region-specific:
# Updated: 03/29/22
# ------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(circlize)
library(ComplexHeatmap)

library(progress)
print(version)
options(width=175)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
srdir = paste0(sdbdir, 'subtype_reg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


processFit = function(fit, gene){
    cfit = coefficients(summary(fit))
    df = data.frame(cfit)
    colnames(df) = c('Est','SE','t','p')
    df$var = rownames(df)
    rownames(df) = NULL
    df$gene = gene
    return(df)
}


# Set run parameters:
# -------------------
subset = 'Mic_Immune_Mic'
ststr = 'Mic_Immune'


# Load the DEGs:
# --------------
full.rds = paste0(regdir, 'aggregated_fullset.', subset, '.rds')
dedf = readRDS(full.rds)

degenes = unique(dedf$gene[dedf$col_nm != 0])
print(length(degenes))


# Load pseudobulk data:
# ---------------------
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
load(psdata.rda)

pmat = ps.data$mat
umeta = unique(ps.data$meta)

# Further annotate the pseudo-bulk metadata:
umeta = merge(umeta, unique(metadata[,c('projid','region', 'rind',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
umeta = merge(umeta, pqdf, all.x=TRUE)
umeta$age_rescaled = umeta$age_death / 100
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(pmat),]

# Remove very low abundance batches:
umeta = umeta[!(umeta$cell_type_high_resolution %in% c('T cells', 'CAMs')), ]
umeta = umeta[umeta$ncell > 10,] 
pmat = pmat[,umeta$ptype]
kept.genes = degenes[degenes %in% rownames(pmat)]

umeta$nft = log1p(umeta$nft)
umeta$plaq_n = log1p(umeta$plaq_n)
umeta$plaq_d = log1p(umeta$plaq_d)


# Run regression at the subtype level:
# ------------------------------------
# TODO: reduce to regions/subset?
ext.covars = '+ Apoe_e4 + age_rescaled + msex + pmi'
reg.var = 'plaq_n'
ind.var = 'cogdxad'
regdf = c()
pb <- progress_bar$new(format="[:bar] :percent eta: :eta",
  total=length(degenes), clear=FALSE, width=60)
for (gene in degenes){
    pb$tick()
    x = pmat[gene, umeta$ptype]
    umeta$expr = x
    # form = asform(c('expr ~ cell_type_high_resolution + region + ',
    #         reg.var, ' * region + region * ', ind.var, ext.covars))
    form = asform(c('expr ~ region + ', reg.var, ' + ', ind.var))
    fit = glm(form, umeta, weights=log(umeta$ncell), 
        family='gaussian')  # Corrects for cell ct. but not inflated
    df = processFit(fit, gene)
    regdf = rbind(regdf, df)
}


# Plot individual vs. region variables:
# -------------------------------------
indstr = paste0(ind.var, 'AD')
varmap = c('p1','p2')
covars = c(reg.var, indstr)
names(varmap) = covars
subdf = regdf[regdf$var %in% covars, ]
ewide = spread(subdf[,c('gene','Est','var')], 'var', 'Est')
names(ewide)[2:3] = varmap[names(ewide[2:3])]
# pwide = spread(subdf[,c('gene','p','var')], 'var', 'p')
# labdf = ewide[abs(ewide[[reg.var]]) > 1 / 25 | abs(ewide[[indstr]]) > 1,]
labdf = ewide[abs(ewide$p1) > 1 /25 | abs(ewide$p2) > 1,]

gp = ggplot(ewide, aes(p1, p2)) + 
    geom_smooth(method='lm') + 
    geom_point(cex=.5, color='grey50') + 
    geom_text_repel(data=labdf, aes(p1, p2, label=gene), cex=2, max.overlaps=20) +
    labs(x=covars[1], y=covars[2]) + 
    geom_hline(yintercept=0, lty='dashed') +
    geom_vline(xintercept=0, lty='dashed') +
    stat_cor() + 
    theme_pubr()

pltprefix = paste0(imgpref, 'ADvar_comparison_scatter_', reg.var, '_', ind.var)
saveGGplot(gp, pltprefix, w=8, h=7)


# Run regression at the subtype level:
# ------------------------------------
gene = 'MT-CO2'
umeta$projid = as.character(umeta$projid)
pids = unique(umeta$projid)
rids = unique(umeta$region)

NR = length(rids)
vec = c()
for (i in 1:NR){
    if (i+1 <= NR){
        for (j in (i+1):NR){
            vec = c(vec, paste0(rids[i], '-', rids[j]))
        }
    }
}

# Turn region x projid matrix > correlation > upper vector 
regionCorVec = function(gene){
    x = pmat[gene, umeta$ptype]
    indices = umeta[,c('region', 'projid')]
    xmat = matrix(0, nr=length(rids), nc=length(pids), dimnames=list(rids, pids))
    xmat[as.matrix(indices)] = x
    xmat = t(scale(t(xmat)))
    xcr = cor(t(xmat), method='spearman')
    xcr = xcr[upper.tri(xcr)]
    return(xcr)
}

xcmat = mclapply(degenes, regionCorVec, mc.cores=20)
xcmat = sapply(xcmat, function(x){x})
xcmat = t(xcmat)
rownames(xcmat) = degenes
colnames(xcmat) = vec


# Cluster the genes by their correlation vectors:
# -----------------------------------------------
ux = 2
plt = Heatmap(xcmat, 
    use_raster=TRUE, 
    col=colb,
    row_km=5,
    width=ncol(xcmat) * unit(ux, 'mm'),
    height=nrow(xcmat) * unit(ux / 10, 'mm'), 
    show_row_names=FALSE)

pltprefix = paste0(imgpref, 'ADvar_comparison_correlation_heatmap_', ststr)
w = 2 + ncol(xcmat) / 15
h = 2 + nrow(xcmat) / 15 / 5
saveHeatmap(plt, pltprefix, w=w, h=h)


pgenes = c('OXR1','MT-CO3','CD74','AOAH','SLC11A1')
xcmat[pgenes,]


regionCorVec(gene)








