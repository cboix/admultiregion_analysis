#!/usr/bin/R
# ------------------------------------------------------------------------
# Vignette 01 - glycolysis associated metabolic changes across cell types:
# Updated 12/20/2021
# ------------------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
regdir = paste0(sdbdir, 'dereg/')
crossdir = paste0(sdbdir, 'vignettes/')
plotdir = paste0(imgdir, 'vignettes/')
imgpref = paste0(plotdir, 'v01_metab_')
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Load the modules in OPC / Mic / Ast:
# ------------------------------------
graph_id = 'boot'
gmlist = list()
for (runset in c('Ast','Mic_Immune','Opc')){
    commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
    source(paste0(sbindir, 'modules/load_modules_degenr.R'))
    print(runset)
    gm = coremap['PDK1']
    gmlist[[runset]] = sort(names(coremap)[coremap == gm])
}

coregenes = intersect(intersect(gmlist[['Opc']], 
                                gmlist[['Ast']]),
                      gmlist[['Mic_Immune']])
allgenes = sort(unique(c(unlist(gmlist))))
cat(sort(coregenes),'\n')


# Load DEG status for these genes (show whether or not DEGs by region, ct):
# -------------------------------------------------------------------------
region = 'allregions'
full.file = paste0(regdir, 'aggregated_allres.', region, '.rda')
load(full.file)

dedf = c()
for (set in names(setdflist)){
    setdf = setdflist[[set]]
    if (!is.null(setdf)){
        setdf$set = set
        dedf = rbind(dedf, setdf)
    }
}

path = 'cogdxad'
dedf = dedf[dedf$gene %in% allgenes,]
subdf = dedf[dedf$path == path,]
subdf$lp = ifelse(subdf$col_nm > 0, subdf$log10p_nm, 0)
cmat = pivot.tomatrix(subdf[,c('gene','set','logFC_nb')], 'set', 'logFC_nb')
pmat = pivot.tomatrix(subdf[,c('gene','set','lp')], 'set', 'lp')
if (path %in% c('nft','plaq_n','plaq_d')){
    col_fun = colorRamp2(c(-.025, 0, .025), c('blue', "white", 'red'))
} else {
    col_fun = colorRamp2(c(-.25, 0, .25), c('blue', "white", 'red'))
}

pmat = 10**(-pmat)
cmat[is.na(cmat)] = 0
pmat[is.na(pmat)] = 1

full.cmat = matrix(0, nrow=length(allgenes), ncol=ncol(cmat), dimnames=list(allgenes, colnames(cmat)))
full.pmat = matrix(1, nrow=length(allgenes), ncol=ncol(pmat), dimnames=list(allgenes, colnames(pmat)))
full.cmat[rownames(cmat),] = cmat
full.pmat[rownames(pmat),] = pmat




# Load pseudobulk expression across cell types:
# ---------------------------------------------
reg.ps.rda = paste0(srdir, 'pseudobulk_data_all_indregion.rda')
if (!file.exists(reg.ps.rda)){
    umeta = NULL
    pmat = NULL
    runlist = c('Ast', 'Mic_Immune', 'Vasc_Epithelia',
                'Oli','Inh','Opc', # 'Exc')
                'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons')
    for (runset in runlist){
        print(runset)
        psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
        load(psdata.rda)
        # Metadata for the pseudobulk matrix:
        ps.data$meta$runset = runset
        umeta = rbind(umeta, ps.data$meta)
        pmat = cbind(pmat, ps.data$mat)
    }

    umeta$runset[umeta$runset %in% names(exc.sets)] = 'Exc'
    save(umeta, pmat, file=reg.ps.rda)
} else {
    load(reg.ps.rda)
}

umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
rownames(umeta) = umeta$ptype
# Remove very low abundance batches + will use for weight
umeta = umeta[umeta$ncell > 25,] 


# Score all modules for (a) all genes and (b) tested DE genes:
# ------------------------------------------------------------
# For all genes:
# kept.genes = rownames(pmat)
# kept.genes = kept.genes[kept.genes %in% names(genemap)]
imgsuff = 'allgenes'

# Conversion matrix:
# tform = make.tform(genemap[kept.genes], u=modules, norm=TRUE)
# kept.mat = t(tform) %*% ps.data$mat[kept.genes,]
# kept.mat = as.matrix(kept.mat)
# rownames(kept.mat) = mmap$mname[1:nrow(kept.mat)]

kept.mat = pmat[allgenes,umeta$ptype]

full.projids = sort(unique(cellmeta$projid))
projid.cols = snap.cols[1:48]
names(projid.cols) = full.projids

# Make annotation from pseudobulk data:
# clsplit = umeta$nrad
clsplit = ifelse(umeta$runset %in% c('Mic_Immune','Opc','Ast'), 'Mic/Opc/Ast', 'All Other')
clsplit = umeta$runset
nft.col_fun = colorRamp2(range(umeta$nft), c("white", "indianred"))
ha = HeatmapAnnotation(CT=umeta$cell_type_high_resolution, 
                       Region=umeta$region,
                       Braak=umeta$braaksc,
                       NFT=umeta$nft,
                       nrad=umeta$nrad,
                       cogdxad=umeta$cogdxad,
                       ncell=umeta$ncell,
                       e4=umeta$Apoe_e4,
                       projid=as.character(umeta$projid),
                       col=list(CT=tcols,
                                Region=reg.cols,
                                Braak=colvals[['braaksc']],
                                nrad=colvals[['nrad']],
                                NFT=nft.col_fun,
                                cogdxad=colvals[['cogdxad']],
                                projid=projid.cols,
                                e4=c('no'='grey90','yes'='slateblue')
                                ))

pltmat = log1p(kept.mat)[, umeta$ptype]
# pltmat = pltmat / apply(pltmat, 1, max)
pltmat = t(scale(t(pltmat)))
plt = Heatmap(pltmat, 
              name='Module\nscore\n(scaled)', 
              use_raster=TRUE,
              top_annotation=ha, 
              column_split=clsplit, 
              show_column_names=FALSE,
              )

plt2 = Heatmap(full.cmat,
               col=col_fun,
               use_raster=TRUE,
               column_split=ifelse(1:ncol(full.cmat) %in% grep("Vasc",colnames(cmat)),'Vasculature','All'),
               cluster_columns=TRUE,
               cluster_rows=TRUE,
               width = ncol(full.cmat)*unit(8, "mm"), 
               height = nrow(full.cmat)*unit(4, "mm"),
               border_gp = gpar(col="black", lty = 1),
               cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
                   p = full.pmat[i,j]
                   ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
                   grid.text(ann, x, y)}
)

plt3 = plt + plt2

h = 2.25 + 2.5 / 15 * nrow(pltmat)
w = 5 + 2.5 / 50 * (ncol(pltmat) + ncol(full.cmat))
w = ifelse(w > 25, 25, w)
pltprefix = paste0(imgpref, 'pseudobulk_heatmap_scaled_', imgsuff)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt3)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt3)
dev.off()




