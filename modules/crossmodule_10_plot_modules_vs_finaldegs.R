#!/usr/bin/R
# ---------------------------------------------------
# Calculate enrichments for the final aggregated DEGs 
# vs. the final modules (+ save resources for plots).
# Updated 03/18/2021 to add plots for each module
# ---------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
options(width=175)

# Settings for plots:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))

# Directories:
moddir = paste0(sdbdir, 'modules/')
regdir = paste0(sdbdir, 'dereg/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'module_vs_degs_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Load in processed enrichments tables:
# -------------------------------------
sets = c('Ast'='Ast_Ast','Mic_Immune'='Mic_Immune_Mic',
    'Opc'='Opc_Opc','Oli'='Oli_Oli','Inh'='Inh_Inh', 
    'Exc'='Exc_Exc', 'All'='All_All')

fulldf = c()
for (use.core in c(TRUE, FALSE)){
    if (use.core){ 
        dbpref = paste0(moddir, 'modenr_core.on_aggregated_fullset.')
    } else { 
        dbpref = paste0(moddir, 'modenr_all.on_aggregated_fullset.')
    }

    # Load each runset:
    for (i in 1:length(sets)){
        set = sets[i]  # DEGs
        runset = names(sets)[i]  # Modules
        cat(set,'\n')

        # Load the table of enrichments:
        enr.rds = paste0(dbpref, set, '.', runset, '.rds')
        ctdf = readRDS(enr.rds)
        ctdf$runset = runset
        ctdf$use.core = use.core
        fulldf = rbind(fulldf, ctdf)
    }
} 

mndf = unique(fulldf[,c('mname','runset')])
mnmap = mndf$runset
names(mnmap) = mndf$mname
fulldf$abs.lfc = abs(fulldf$log2FC)



# Plot heatmaps of enrichments:
# -----------------------------
use.core = FALSE
imgpref = paste0(plotdir, 'module_vs_degs_')
if (use.core){ 
    imgpref = paste0(imgpref, 'core_')
} else { 
    imgpref = paste0(imgpref, 'all_')
}

ctdf = fulldf[fulldf$use.core == use.core,]
ctdf$pr = paste0(ctdf$path, '_', ctdf$region)

# Aggregate up/down by most significant, breaking ties by log2FC and count:
ctdf = merge(ctdf, aggregate(p.adj ~ pr + mname, ctdf, min))
ctdf = merge(ctdf, aggregate(abs.lfc ~ pr + mname, ctdf, max))
ctdf = merge(ctdf, aggregate(count ~ pr + mname, ctdf, max))
# Flip direction for depletion:
ctdf$log2FC[ctdf$col_nm == 1] = -ctdf$log2FC[ctdf$col_nm == 1]

# Keep only significant modules:
kept.mn = unique(ctdf$mname[ctdf$p.adj < 0.05])
ctdf = ctdf[ctdf$mname %in% kept.mn,]
col_fun = colorRamp2(c(-4, 0, 4), c('blue', "white", 'red'))

# Turn into matrices + plot:
rmat = pivot.tomatrix(ctdf[,c('mname', 'pr', 'log2FC')], 'pr','log2FC')
pmat = pivot.tomatrix(ctdf[,c('mname', 'pr', 'p.adj')], 'pr', 'p.adj')
rmat[is.na(rmat)] = 0
pmat[is.na(pmat)] = 1

# Reorder columns:
pathlist = c('plaq_d','plaq_n','nft','nrad','cogdxad')
reglist = c('allregions', 'AG','MT','PFC', 'TH','HC','EC')
coldf = expand.grid(region=reglist, path=pathlist)
coldf$pr = with(coldf, paste0(path, '_', region))
cn = coldf$pr[coldf$pr %in% colnames(rmat)]
rmat = rmat[,cn]
pmat = pmat[,cn]

row.split = mnmap[rownames(rmat)]
column.split = sub("_[A-Za-z]*$","", colnames(rmat))


ux = 1.5
ht = Heatmap(rmat, 
    name='-log2FC',
    cluster_rows=TRUE,
    cluster_columns=FALSE,
    cluster_column_slices=FALSE,
    cluster_row_slices = TRUE,
    column_split=column.split,
    row_split=row.split,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(rmat) * unit(ux, "mm"), 
    height = nrow(rmat) * unit(ux, "mm"),
    # col=rev(colrb),
    col=col_fun,
    cell_fun = function(j, i, x, y, w, h, col){
        p = pmat[i,j]
        lfc = rmat[i,j]
        if (p < 0.05){
            # ann = ifelse(p < 0.01, ifelse(p < 0.001, '***', '**'), '*')
            ann = '*'
            grid.text(ann, x, y, 
                gp=gpar(col=ifelse(lfc > 6, 'white','black'), 
                    fontsize=gridtxt.fs))
        }
    }
)

h = 2 + 1 / 15 * nrow(rmat)
w = 3 + 1 / 15 * ncol(rmat)
pltprefix = paste0(imgpref, 'full_degenr_heatmap_allruns')
saveHeatmap(ht, pltprefix, w, h)


# Only plot the allregions runs and add top genes by cognition:
# -------------------------------------------------------------
ctdf = fulldf[fulldf$use.core == use.core & fulldf$region == 'allregions',]

# Aggregate up/down by most significant, breaking ties by log2FC and count:
ctdf = merge(ctdf, aggregate(p.adj ~ path + mname, ctdf, min))
ctdf = merge(ctdf, aggregate(abs.lfc ~ path + mname, ctdf, max))
ctdf = merge(ctdf, aggregate(count ~ path + mname, ctdf, max))
# Flip direction for depletion:
ctdf$log2FC[ctdf$col_nm == 1] = -ctdf$log2FC[ctdf$col_nm == 1]

# Keep only modules significant in two + more genes:
ctdf$keep = ctdf$p.adj < 0.05 & ctdf$count >= 5
keepdf = aggregate(keep ~ mname, ctdf, sum)
kept.mn = keepdf$mn[keepdf$keep >= 2]
# kept.mn = unique(ctdf$mname[ctdf$p.adj < 0.05 & ctdf$count >= 5])
ctdf = ctdf[ctdf$mname %in% kept.mn,]
col_fun = colorRamp2(c(-4, 0, 4), c('blue', "white", 'red'))

# Annotate by maximal enrichment:
path = 'cogdxad'
topdf = merge(ctdf, aggregate(p.adj ~ mname, ctdf, min))
topdf = merge(topdf, aggregate(abs.lfc ~ mname, topdf, max))
topdf = merge(topdf, aggregate(count ~ mname, topdf, max))

dirdf = unique(topdf[, c('path','mname','runset','region', 'col_nm', 'use.core')])
dirdf = merge(dirdf, fulldf)
genes.mapping = dirdf$genes
path.mapping = dirdf$path
dircol = dirdf$col_nm
names(path.mapping) = dirdf$mname
names(genes.mapping) = dirdf$mname
names(dircol) = dirdf$mname

# Turn into matrices + plot:
rmat = pivot.tomatrix(ctdf[,c('mname', 'path', 'log2FC')], 'path','log2FC')
pmat = pivot.tomatrix(ctdf[,c('mname', 'path', 'p.adj')], 'path', 'p.adj')
rmat[is.na(rmat)] = 0
pmat[is.na(pmat)] = 1
row.split = mnmap[rownames(rmat)]


hmod = rowAnnotation(mod = anno_text(sub(" .*","", rownames(rmat)), gp=gpar(fontsize=5)))
hpath = rowAnnotation(path = anno_text(path.mapping[rownames(rmat)],
        gp=gpar(fontsize=5, col=ifelse(dircol[rownames(rmat)] == 2, colrb[10],colrb[90]))))
hgenes = rowAnnotation(top.genes = anno_text(genes.mapping[rownames(rmat)],
        gp=gpar(fontsize=5, col=ifelse(dircol[rownames(rmat)] == 2, colrb[10],colrb[90]))))

ux = 1.5
ht = Heatmap(rmat, 
    name='-log2FC',
    cluster_rows=TRUE,
    cluster_columns=FALSE,
    cluster_column_slices = TRUE,
    cluster_row_slices = TRUE,
    row_split=row.split,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(rmat) * unit(ux, "mm"), 
    height = nrow(rmat) * unit(ux, "mm"),
    # col=rev(colrb),
    col=col_fun,
    cell_fun = function(j, i, x, y, w, h, col){
        p = pmat[i,j]
        lfc = rmat[i,j]
        if (p < 0.05){
            # ann = ifelse(p < 0.01, ifelse(p < 0.001, '***', '**'), '*')
            ann = '*'
            grid.text(ann, x, y, 
                gp=gpar(col=ifelse(lfc > 6, 'white','black'), 
                    fontsize=gridtxt.fs))
        }
    }
)

ht = ht + hmod + hpath + hgenes

h = 1 + 1 / 10 * nrow(rmat)
w = 5 + 1 / 10 * ncol(rmat)
pltprefix = paste0(imgpref, 'full_degenr_heatmap')
saveHeatmap(ht, pltprefix, w, h)



