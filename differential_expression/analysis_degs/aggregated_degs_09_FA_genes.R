#!/usr/bin/R
# -----------------------------------------------------
# Look at DEGs for fatty acid metab genes for LA expts.
# Updated: 09/18/23
# -----------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
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
# source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))



# Arguments for runs:
# -------------------
set = "Ast_Ast"
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
remove.shared = TRUE
run.intersections = TRUE

marker.genes = c('GFAP', 'SLC1A2', 'SLC1A3', 'AQP4')
metab.genes = c('PFKP', 'HILPDA', 'ADCY8', 'PFKL')
genes = c(marker.genes, metab.genes)

path = 'nrad'
alldf = c()
for (path in pathlist){
    print(path)
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)

    # Calculate the set shared in each region:
    # ----------------------------------------
    kept.cols = c('gene','log10p_nm', 'logFC_nb', 'p_nb','col_nm','path','region')
    setdf = setdflist[[set]][, kept.cols]
    alldf = rbind(alldf, setdf[setdf$gene %in% genes,])
}

alldf$pr = with(alldf, paste0(path, "@", region))
alldf$lfc.mod = ifelse(alldf$path %in% c('nrad','cogdxad'), alldf$logFC_nb, alldf$logFC_nb * 25)

cmat = pivot.tomatrix(alldf[,c('pr','gene','lfc.mod')], 'pr','lfc.mod')
colmat = pivot.tomatrix(alldf[,c('pr','gene','col_nm')], 'pr','col_nm')
colmat[is.na(colmat)] = 0
colmat[colmat > 0] = 1
pmat = 1 - colmat

# Reorder columns by name of region:
df = expand.grid(region=c('allregions', 'EC', 'HC', 'TH', 'AG', 'MT', 'PFC'), path=pathlist)
cns = paste0(df$path, '@', df$region)
cns = cns[cns %in% colnames(cmat)]
cmat = cmat[,cns]
pmat = pmat[,cns]

mx = 1
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
colsplit = sub("@.*","", colnames(cmat))
rowsplit = ifelse(rownames(cmat) %in% metab.genes, 'Metab', 'Marker')

ux = 1.5
ht = Heatmap(cmat,
    use_raster=FALSE,
    name='logFC',
    col=col_fun,
    cluster_columns=FALSE,
    cluster_rows=TRUE,
    column_split = colsplit,
    row_split = rowsplit,
    width = ncol(cmat)*unit(ux, "mm"), 
    height = nrow(cmat)*unit(ux, "mm"),
    row_dend_width = unit(.25, "cm"),
    column_dend_height = unit(.25, "cm"),
    row_dend_gp = gpar(lwd=.5),
    column_dend_gp = gpar(lwd=.5),
    border_gp = gpar(col="black", lwd=.5),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        p = pmat[i,j]
        if (p < 0.05){
            grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.1))
        }
    })

pltprefix = paste0(imgpref, 'ast_FA_markers_DE_heatmap')
h = 1 + 1 / 15 * nrow(cmat)
w = 2 + 1 / 15 * ncol(cmat)
saveHeatmap(ht, pltprefix, w=w, h=h)



adf = alldf[order(abs(alldf$logFC_nb), decreasing=T),]
adf = alldf[order(alldf$p_nb),]
adf = adf[adf$gene %in% metab.genes,]
adf = adf[adf$col_nm != 0,]
write.table(adf, paste0(regdir, 'FA_genes_DEGstatus.tsv'), quote=F, sep="\t", row.names=F)


# Load expression data to find individuals for testing:
# -----------------------------------------------------
runset = 'Ast'
psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '.rda')
load(psdata.rda)

umeta = ps.data$meta
umeta = unique(umeta)
umeta = merge(umeta, unique(metadata[,c('projid','region',
                                        'braaksc','cogdx', 'niareagansc',
                                        'msex','age_death','pmi', 
                                        'Apoe_e4', 'nrad','cogdxad')]))
rownames(umeta) = umeta$ptype
umeta = umeta[colnames(ps.data$mat),]
# Remove very low abundance batches
umeta = umeta[umeta$ncell > 10,] 
gmat = ps.data$mat[genes,]

# Aggregate at the individual-level:
mod.data = list('mat'=gmat[,umeta$ptype], 'meta'=umeta)
ind.data = aggregatePsbulkIndRegion(mod.data)

# Aggregate over regions without PFC:
umeta = ind.data$meta[ind.data$meta$region != 'PFC',]
mat = scale(t(ind.data$mat[,umeta$pr]))
tform = make.tform(umeta$projid, u=as.character(unique(umeta$projid)), norm=TRUE)
avg.mat = t(tform) %*% mat

mts = apply(avg.mat[,metab.genes], 1, sum)
mks = apply(avg.mat[,marker.genes], 1, sum)
mts = sort(mts, decreasing=T) # Metab score
mks = mks[names(mts)] # Marker score

c(mean(mks[1:24]), mean(mks[25:48]))
c(mean(mts[1:24]), mean(mts[25:48]))

df = data.frame(projid=names(mts), metab.score=mts, marker.score=mks)
rownames(df) = NULL
df = merge(df, unique(metadata[,c('projid','cogdxad','nrad')]))
df = df[order(df$metab.score, decreasing=T),]

write.table(df, paste0(dbdir, 'FA_genes_metabscore.tsv'), quote=F, sep="\t", row.names=F)


