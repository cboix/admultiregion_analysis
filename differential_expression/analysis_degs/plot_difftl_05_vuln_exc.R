#!/usr/bin/R
# ------------------------------------------------
# Differential expression for excitatory subtypes:
# Updated: 09/10/21
# ------------------------------------------------
# Aggregate number of sign. genes:
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
print(version)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'exc-difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# -----------------------------
# Load in the regression files:
# -----------------------------
# Load in the EC / HC / TH differential results
nbvuln.rda = paste0(regdir, 'nebula_ruv.collated_Exc_vuln_ECHC.Rda')
if (!file.exists(nbvuln.rda)){
    # TODO: Repeat with the neocortex results:
    fnlist = list.files(path=regdir, pattern=paste0('nebula_ruv.Exc_.*', '[HTE][HC].*rda'))
    nsigdf = c(); alldf = c();
    for (fn in fnlist){
        load(paste0(regdir, fn))
        region = nsig[1,'region']
        path = nsig[1,'path']
        ndf = data.frame(nsig)
        if (!('X1' %in% colnames(ndf))){ ndf$X1 = 0 }
        if (!('X2' %in% colnames(ndf))){ ndf$X2 = 0 }
        if (is.null(nsigdf)){
            nsigdf = ndf 
        } else {
            nsigdf = rbind(nsigdf, ndf[,colnames(nsigdf), drop=F])
        }
        fulldf$path = path
        fulldf$region = region
        fulldf$subtype = nsig[1, 'subtype']
        # fulldf = fulldf[order(fulldf$col != 0, abs(fulldf$logFC), decreasing=T),]
        fulldf = fulldf[order(fulldf$p),]
        fulldf$rank = 1:nrow(fulldf)
        alldf = rbind(alldf, fulldf)
    }
    names(nsigdf)[names(nsigdf) == 'X1'] = 'ndown'
    names(nsigdf)[names(nsigdf) == 'X2'] = 'nup'
    nsigdf$nup = as.numeric(nsigdf$nup)
    nsigdf$ndown = as.numeric(nsigdf$ndown)
    nsigdf = nsigdf[nsigdf$subtype != 'Exc',]
    save(alldf, nsigdf, file=nbvuln.rda)
} else {
    load(nbvuln.rda)
}

mdf = aggregate(nup ~ subtype, nsigdf, mean)
mdf = mdf[order(-mdf$nup),]

subdf = nsigdf
subdf$subtype = factor(subdf$subtype, levels=rev(mdf$subtype))
pcols = brewer.pal(12, 'Paired')
gplot = ggplot(subdf, aes(subtype, nup)) + 
    # facet_wrap(~path, ncol=1) + 
    facet_wrap(~path, scales='free_x') + 
    geom_bar(position='dodge',stat='identity', fill=pcols[6]) + 
    geom_bar(data=subdf, aes(subtype, -ndown), position='dodge',stat='identity', fill=pcols[2]) + 
    labs(x='AD variable', y='Number of DEGs') + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, 'nbRUV_comparison_ndeg.excsubtypes.ECHCTH.png'), gplot, dpi=400, units='in', width=9,height=5.5)
ggsave(paste0(imgpref, 'nbRUV_comparison_ndeg.excsubtypes.ECHCTH.pdf'), gplot, width=9,height=5.5)


# ---------------------------------------------------
# Look at shared DEGs for cognition (or any measure):
# ---------------------------------------------------
path = 'cogdxad'
subdf = alldf[(alldf$path == path) & (alldf$subtype != 'Exc'),]

vuln = c('CA1_pyramidal_cells','Exc_TOX3_TTC6','Exc_RELN_GPC5',
         'Exc_AGBL1_GPC5','Exc_RELN_COL5A2')
# vuln = unique(alldf$subtype)
# vuln = c('CA1_pyramidal_cells')

updf = subdf[(subdf$col == 2) & (subdf$subtype %in% vuln),]
dwdf = subdf[(subdf$col == 1) & (subdf$subtype %in% vuln),]
updf$log10q = -log10(updf$q)
updf$count = 1
dwdf$log10q = -log10(dwdf$q)
dwdf$count = 1

ctdf = aggregate(cbind(count, log10q) ~ gene, updf, sum)
ctdf = ctdf[order(-ctdf$count, -ctdf$log10q),]
dtdf = aggregate(cbind(count, log10q) ~ gene, dwdf, sum)
dtdf = dtdf[order(-dtdf$count, -dtdf$log10q),]

# Image top up / down genes across subtypes:
# ------------------------------------------
subdf = alldf[(alldf$path == path),]
subdf$subtype[subdf$subtype == "Exc"] = paste0("Exc_", subdf$region[subdf$subtype=='Exc'])

NTOP = 30
gup = head(ctdf$gene[ctdf$count == length(vuln)], NTOP)
gdw = head(dtdf$gene[dtdf$count == length(vuln)], NTOP)
# gup = ctdf$gene[ctdf$count > 2]
# gdw = dtdf$gene[dtdf$count > 2]

topdf = subdf[subdf$gene %in% c(gup, gdw),]

cwide = spread(topdf[,c('subtype','gene','logFC')], subtype, logFC)
pwide = spread(topdf[,c('subtype','gene','p')], subtype, p)
cmat = as.matrix(cwide[,-1])
pmat = as.matrix(pwide[,-1])
rownames(cmat) = cwide$gene
rownames(pmat) = pwide$gene
cmat[is.na(cmat)] = 0
pmat[is.na(pmat)] = 1
reord = reord(cmat)

ind = 1:nrow(reord)
use = round(ind / 6) %% 2
ha = rowAnnotation(foo = anno_mark(at = ind[use == 1], labels = rownames(reord)[use == 1]))
hb = rowAnnotation(foo = anno_mark(at = ind[use == 0], labels = rownames(reord)[use == 0], side='left'))

colsplit = ifelse(colnames(reord) %in% vuln, 'Vulnerable', 'Stable')
colsplit[grep("^Exc_[ECTH]+$",colnames(reord))] = 'Overall'
png(paste0(imgpref, 'ECHCTH_top',NTOP, '_consistentDEGs_heatmap.png'), res=450, units='in', width=6, height=9)
Heatmap(reord, name='logFC',
        cluster_rows=FALSE,
        use_raster=TRUE,
        clustering_distance_rows='euclidean',
        column_split=colsplit,
        row_split=ifelse(rownames(reord) %in% gup, 'Up-regulated','Down-regulated'),
        show_row_names=FALSE,
        right_annotation=ha,
        left_annotation=hb,
        # cell_fun = function(j, i, x, y, w, h, fill) {
        #     p = pmat[i,j] # Adjusted p-values from regression.
        #     ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
        #     grid.text(ann, x, y) }
)
dev.off()



# Plot as a volcano plot:
# -----------------------
subdf$log10q = -log10(subdf$q)
aggdf = merge(aggregate(log10q ~ gene,subdf[subdf$subtype %in% vuln,], sum),
              aggregate(logFC ~ gene,subdf[subdf$subtype %in% vuln,], mean))

gup = ctdf$gene[ctdf$count == length(vuln)]
gdw = dtdf$gene[dtdf$count == length(vuln)]

aggdf$col = 0 
aggdf$col[aggdf$gene %in% gup] = 2
aggdf$col[aggdf$gene %in% gdw] = 1
labdf = aggdf[aggdf$col != 0,]
labdf = labdf[labdf$log10q > 40,]

gplot = ggplot(aggdf, aes(logFC, log10q, col=factor(col))) + 
    geom_point(cex=.5) + 
    theme_pubr() + 
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=0) + 
    geom_text_repel(data=labdf, aes(logFC, log10q, label=gene, col=factor(col)), max.overlaps=20) + 
    labs(x='Mean logFC in vulnerable subtypes', y='Summed log10q across vulnerable subtypes') + 
    scale_color_manual(values=c('lightgrey',pcols[2],pcols[6]))
ggsave(paste0(imgpref, 'ECHCTH_top',NTOP, '_consistentDEGs_jointvuln_volcano.png'), gplot, dpi=450, units='in', width=8, height=9)
ggsave(paste0(imgpref, 'ECHCTH_top',NTOP, '_consistentDEGs_jointvuln_volcano.pdf'), gplot, dpi=450, units='in', width=8, height=9)


# --------------------------
# Load some expression data:
# TODO: Load all ct:
# --------------------------
celltype = 'Exc'
subtype = 'Exc_TOX3_TTC6'
region = 'EC'
commandArgs = function(x){ c(celltype, subtype, region)}
source(paste0(bindir, 'multiRegion/load_difftl_data.R'))

fact = colSums(mat) / 10000

# Correlation matrix: 
# --------------------
NTOP = 100
gup = head(ctdf$gene[ctdf$count == length(vuln)], NTOP)
gdw = head(dtdf$gene[dtdf$count == length(vuln)], NTOP)
gup = ctdf$gene[ctdf$count == length(vuln)]
# gdw = dtdf$gene[dtdf$count > 2]

# Subset + normalize only genes of interest:
norm = as.matrix(mat[c(gup),])
norm = sweep(norm, 2, fact,'/')
norm = log(norm + 1)

cv = cor(t(norm))
cv = cv - diag(diag(cv))

Heatmap(cv)






topdf = subdf[subdf$gene %in% c(gup, gdw),]

cwide = spread(topdf[,c('subtype','gene','logFC')], subtype, logFC)
pwide = spread(topdf[,c('subtype','gene','p')], subtype, p)
cmat = as.matrix(cwide[,-1])
pmat = as.matrix(pwide[,-1])
rownames(cmat) = cwide$gene
rownames(pmat) = pwide$gene
cmat[is.na(cmat)] = 0
pmat[is.na(pmat)] = 1
reord = reord(cmat)

ind = 1:nrow(reord)
use = round(ind / 6) %% 2
ha = rowAnnotation(foo = anno_mark(at = ind[use == 1], labels = rownames(reord)[use == 1]))
hb = rowAnnotation(foo = anno_mark(at = ind[use == 0], labels = rownames(reord)[use == 0], side='left'))

colsplit = ifelse(colnames(reord) %in% vuln, 'Vulnerable', 'Stable')
colsplit[grep("^Exc_[ECTH]+$",colnames(reord))] = 'Overall'
png(paste0(imgpref, 'ECHCTH_top',NTOP, '_consistentDEGs_heatmap.png'), res=450, units='in', width=6, height=9)
Heatmap(reord, name='logFC',
        cluster_rows=FALSE,
        use_raster=TRUE,
        clustering_distance_rows='euclidean',
        column_split=colsplit,
        row_split=ifelse(rownames(reord) %in% gup, 'Up-regulated','Down-regulated'),
        show_row_names=FALSE,
        right_annotation=ha,
        left_annotation=hb,
        # cell_fun = function(j, i, x, y, w, h, fill) {
        #     p = pmat[i,j] # Adjusted p-values from regression.
        #     ann = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, '***','**'),'*'),'')
        #     grid.text(ann, x, y) }
)
dev.off()






# Co-expression network:







# Localization on the UMAP:







