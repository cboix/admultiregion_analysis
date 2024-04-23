#!/usr/bin/R
# ---------------------------------------------------
# Look at the astrocyte + oligodendrocyte modules
# associated with plaque over NFT, plot glyc.
# Updated 11/04/2023
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
srdir = paste0(sdbdir, 'subtype_reg/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = paste0(plotdir, 'module_selection_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the modules:
# --------------------
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Load in DEG data:
# -----------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')

path = 'plaq_n'
mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)


# Calculate the set shared in each region:
# ----------------------------------------
ctlist = list('Ast'=c(6,12), 'Oli'=c(7,9,16))
resdf = c()
kept.cols = c('gene','de','logFC_nb', 'p_nb', 'path','region', 'celltype', 'module', 'uselist')
for (uselist in c('all','core')){
    for (ct in names(ctlist)){
        set = paste0(ct, '_', ct)
        print(set)

        suff = paste0(path, '.', set)
        fulldf = setdflist[[set]]
        fulldf$de = ifelse(fulldf$col_nm == 1, 'Down', ifelse(fulldf$col_nm == 2, 'Up', 'NS'))
        fulldf$pr = paste0(fulldf$region, '_', fulldf$de)

        if (uselist == 'core'){
            usemap = cmlist[[ct]]
        } else {
            usemap = gmlist[[ct]] # Use all
        }

        fulldf = fulldf[fulldf$gene %in% names(usemap),]
        fulldf$module = usemap[fulldf$gene]

        subdf = fulldf[fulldf$module %in% ctlist[[ct]],]
        subdf$celltype = ct
        subdf$uselist = uselist
        for (module in ctlist[[ct]]){
            region = 'allregions'
            moddf = subdf[subdf$module == module,]
            # print(moddf[moddf$region == region,kept.cols])

            cmat = pivot.tomatrix(moddf[,c('gene','region','logFC_nb')], 'region','logFC_nb')
            pmat = pivot.tomatrix(moddf[,c('gene','region','col_nm')], 'region','col_nm')

            mmarg = apply(cmat, 1, mean, na.rm=T)
            roword = order(mmarg, decreasing=T)
            colord = c('allregions', 'EC', 'HC', 'AG', 'MT', 'PFC')

            cmat = cmat[roword, colord]
            pmat = pmat[roword, colord]

            ux = 1.5
            ht = Heatmap(cmat,
                use_raster=FALSE,
                name='logFC',
                cluster_columns=FALSE,
                cluster_rows=FALSE,
                column_title=paste0(ct, ' - M', module, '\n(', uselist, ' genes)'),
                width = ncol(cmat)*unit(ux, "mm"), 
                height = nrow(cmat)*unit(ux, "mm"),
                row_dend_width = unit(.25, "cm"),
                column_dend_height = unit(.25, "cm"),
                row_dend_gp = gpar(lwd=.5),
                column_dend_gp = gpar(lwd=.5),
                border_gp = gpar(col="black", lwd=.5),
                cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                    p = pmat[i,j]
                    if (!is.na(p) & (p > 0)){
                        grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.1))
                    }
                })

            pltprefix = paste0(imgpref, 'plaq_assoc_heatmap.', uselist, '.', ct, '-M', module)
            h = .5 + 1 / 15 * nrow(cmat)
            w = 1 + 1 / 15 * ncol(cmat)
            saveHeatmap(ht, pltprefix, w=w, h=h)
        }
        resdf = rbind(resdf, subdf[, kept.cols])
    }
}

outfile = paste0(regdir, 'astoli_plaque_assoc_modules_genes.tsv')
write.table(resdf, outfile, quote=F, sep="\t", row.names=F)


# Also plot heatmap of glycolysis genes, all regions, diffuse plaque
# ------------------------------------------------------------------
df = read.delim(paste0(sdbdir, 'WP534-datanodes.tsv'), header=T)
genes = df$Label[df$Type == 'GeneProduct']
genes = c(genes, 'PKM')  # is PKM1/2
# NOTE: PGI == GPI!?
genes[!(genes %in% anno$symbol)]
genes = genes[!(genes %in% c('PGI', 'PKM1', 'PKM2'))]

df.glycogen = read.delim(paste0(sdbdir, 'WP500-datanodes.tsv'), header=T)
genes.glycogen = df.glycogen$Label[df.glycogen$Type == 'GeneProduct']
genes.glycogen = c(genes.glycogen, 'PYGM', 'PYGB', 'PYGL', 'GYG1', 'GYG2')
genes.glycogen[!(genes.glycogen %in% anno$symbol)]

df.lipid = read.delim(paste0(sdbdir, 'WP3965-datanodes.tsv'), header=T)
genes.lipid = df.lipid$Label[df.lipid$Type == 'GeneProduct']
genes.lipid[!(genes.lipid %in% anno$symbol)]


# Intersect with DEGs:
# --------------------
plotSetHeatmap = function(df, pltprefix){
    cmat = as.matrix(df$logFC_nb)
    pmat = as.matrix(df$col_nm)
    rownames(cmat) = df$gene
    rownames(pmat) = df$gene
    ux = 1.5
    ht = Heatmap(cmat,
        use_raster=FALSE,
        name='logFC',
        width = ncol(cmat)*unit(ux, "mm"), 
        height = nrow(cmat)*unit(ux, "mm"),
        row_dend_width = unit(.25, "cm"),
        column_dend_height = unit(.25, "cm"),
        row_dend_gp = gpar(lwd=.5),
        column_dend_gp = gpar(lwd=.5),
        border_gp = gpar(col="black", lwd=.5),
        cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
            p = pmat[i,j]
            if (!is.na(p) & (p > 0)){
                grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.1))
            }
        })
    h = .5 + 1 / 15 * nrow(cmat)
    w = 1 + 1 / 15 * ncol(cmat)
    saveHeatmap(ht, pltprefix, w=w, h=h)
}


path = 'plaq_d'
set = 'Ast_Ast'
print(paste(path, set))

mstr = paste0('allmethods.regional_', path)
fullaggrda = paste0(regdir, mstr, '.merged.rda')
load(fullaggrda)
fulldf = setdflist[[set]]

scale_coef = sd(fulldf$coef_mast) / sd(fulldf$logFC_nb)
fulldf$avg_mn = (fulldf$coef_mast / scale_coef + fulldf$logFC_nb)

# Intersection:
keep.cols = c('gene','col_nm', 'col_mast', 'p_mast', 'logFC_nb')
glydf = fulldf[(fulldf$region == 'allregions') & (fulldf$gene %in% genes),]
print(table(glydf$col_nm))
glydf[,keep.cols]

plotSetHeatmap(glydf, paste0(imgpref, 'ast_plaq_d_glycolysis_genes'))


gendf = fulldf[(fulldf$region == 'allregions') &(fulldf$gene %in% genes.glycogen),]
print(table(gendf$col_nm))
gendf[,keep.cols]

lipdf = fulldf[(fulldf$region == 'allregions') &(fulldf$gene %in% genes.lipid),]
print(table(lipdf$col_nm))
lipdf[,keep.cols]

genes.ady = c('PLIN1', 'ADCY8', 'PNPLA2', 'LIPE', 'PRKACA', 'PRKACB', 'PRKACG', 'GNAS', 'GCG', 'PRKAG1', 'PRKAG2', 'PPARG', 'UCP1')

adydf = fulldf[(fulldf$region == 'allregions') &(fulldf$gene %in% genes.ady),]
adydf = fulldf[(fulldf$gene %in% genes.ady),]
print(table(adydf$col_nm))
adydf[,keep.cols]


