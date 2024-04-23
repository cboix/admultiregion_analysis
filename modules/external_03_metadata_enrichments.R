#!/usr/bin/R
# ---------------------------------------------------
# Get the metadata enrichments for external datasets:
# Updated 12/01/2023
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
options(width=170)

source(paste0(sbindir, 'modules/auxiliary_gprofiler_functions.R'))

# Set the run arguments:
# ----------------------
# dataset = 'Mathys_Nature2019'
# dataset = 'TabulaSapiens'
# dataset = 'COVID19_Cell2021'
# dataset = 'EasySci_2023'
# graph_id = 'base'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    dataset = args[1]
    graph_id = args[2]
}

# Directories:
extdir = paste0(sdbdir, 'external_datasets/')
dsdir = paste0(extdir, dataset, '/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'external_', dataset, '_')


# Load in the modules data:
# -------------------------
commandArgs <- function(trailingOnly=TRUE){c(dataset, graph_id, FALSE)}
source(paste0(sbindir, 'modules/load_external_modules_data.R'))

# Functional enrichments data:
set = 'allgenes'
toppvals.file = paste0(dsdir, 'module_enrichments_toppvals_',set,'_', graphpref, '_small.tsv')
pvalsdf = read.delim(toppvals.file, header=T)

# Load metadata counts:
hgfile = paste0(dsdir, 'modules_hgdata_', graphpref, '.tsv')
hgdf = read.delim(hgfile, header=T)[,-1]
hgdf$p <- apply(hgdf[,c('nint','nmod','ncat','ntot')], 1, run.hyper)
hgdf$p.adj = p.adjust(hgdf$p, 'BH')
hgdf$log10p = -log10(hgdf$p)
hgdf$logFC = with(hgdf, log((nint / nmod) / (ncat/ ntot)))
hgdf = hgdf[order(hgdf$logFC, decreasing=T),]
hgdf = hgdf[order(hgdf$p),]
head(hgdf, 20)

if (dataset == 'COVID19_Cell2021'){
    kept.vars = c('majorType','Outcome',
        'CoVID-19 severity', 'Sex','SARS-CoV-2')
    ctvar = 'majorType'
} else if (dataset == 'Mathys_Nature2019'){
    kept.vars = c('Subcluster','broad.cell.type')
    ctvar = 'broad.cell.type'
} else if (dataset == 'EasySci_2023'){
    kept.vars = unique(hgdf$var) 
    ctvar = 'Cell_type'
} else if (dataset == 'TabulaSapiens'){
    kept.vars = unique(hgdf$var) 
    ctvar = 'organ_tissue'
} else { 
    kept.vars = unique(hgdf$var) 
} 
sub.hgdf = hgdf[hgdf$var %in% kept.vars,]
ll = lapply(unique(hgdf$module), function(x){ 
    head(sub.hgdf[sub.hgdf$module == x,],1) })
rdf = do.call(rbind, ll)

# Report these results:
sdf = merge(rdf[,c('module','var','cat','logFC','log10p')], 
    pvalsdf[,c('module','p_value','term_name')], all=TRUE)

# min.genes = 10
min.genes = 10
ngdf = data.frame(table(coremap))
names(ngdf) = c('module', 'ngenes')
mmap = merge(mmap, ngdf)
sdf = merge(sdf, mmap)
sdf = sdf[sdf$ngenes >= min.genes,]
sdf = merge(sdf, unique(rdf[,c('module','nmod')]))
sdf = sdf[order(sdf$module),]
rownames(sdf) = sdf$mname


# Covariate heatmap:
# ------------------
hgdf$vc = paste0(hgdf$var, '@', hgdf$cat)
if (dataset %in% c('EasySci_2023', 'COVID19_Cell2021')){
    use.hgdf = hgdf[hgdf$var %in% kept.vars,]
    use.var = 'vc'
} else {
    use.hgdf = hgdf[hgdf$var == ctvar,]
    use.var = 'cat'
}
use.hgdf = merge(mmap, use.hgdf)

cmat = pivot.tomatrix(use.hgdf[, c('mname', use.var, 'logFC')], use.var,'logFC')
pmat = pivot.tomatrix(use.hgdf[, c('mname', use.var, 'p')], use.var,'p')
pmat = -log10(pmat)
cmat = cmat[as.character(sdf$mname),]
pmat = pmat[as.character(sdf$mname),]
mx = max(cmat)
mx = 2
cmat[cmat < -mx] = -mx
# cov.col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", 'red'))
cov.col_fun = colorRamp2(c(0, mx), c("white", 'red'))
cmat = t(reord(t(cmat)))
cmat = reord(cmat)
cmat = t(diag.mat2(t(cmat))[[1]])
cmat = cmat[rev(rownames(cmat)),]
cmat[cmat < 0] = 0
pmat = pmat[rownames(cmat), colnames(cmat)]

sdf = sdf[rownames(cmat),]
topids = rownames(cmat)[apply(cmat, 2, which.max)]
if (use.var == 'vc'){
    ctind = ctvar == sub("@.*", "", colnames(cmat))
    topids = topids[ctind]
}
idsplit = ifelse(sdf$mname %in% topids, '0_Identity', '1_Functional')


# Elements for vis:
# -----------------
ux = 1.75
rn = sdf$mname
# Text:
hg1 = rowAnnotation(top.genes = anno_text(sub(".*\\: ","",rn), gp=gpar(fontsize=5)))
hn1 = rowAnnotation(top.genes = anno_text(sub(" \\(.*","",rn), gp=gpar(fontsize=5)))


# Colors:
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))
col_log10p = colorRamp2(c(0, 10), c('white','black'))
colg = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(50)
pcols = brewer.pal(12,'Paired')

# Number of genes and pct of cells:
ngmat = as.matrix(sdf[, 'ngenes', drop=F])
htng = Heatmap(ngmat, 
    cluster_rows=FALSE,
    name='ngenes',
    width = ncol(ngmat)*unit(ux * 2, "mm"), 
    height = nrow(ngmat)*unit(ux, "mm"),
    row_split=idsplit,
    border_gp=gpar(color='black', lwd=.5),
    col=colg,
    cell_fun = function(j, i, x, y, w, h, col){
        grid.text(ngmat[i,j], x, y,
            gp=gpar(col=ifelse(ngmat[i,j] > .60 * max(ngmat),'white','black'), 
                fontsize=gridtxt.fs))}
)

# Create matrix for plotting enrichment:
# TODO: fill with NA if missing vals
enrmat = as.matrix(-log10(sdf[,'p_value',drop=F]))
enrmat[is.na(enrmat)] = 0
colnames(enrmat) = '-log10p'
rownames(enrmat) = sdf$term_name
th.enrmat = enrmat  # Thresholded, for colors:
th.enrmat[th.enrmat > 10] = 10
htenr = Heatmap(th.enrmat, 
    cluster_rows=FALSE,
    name='-log10p',
    row_split=idsplit,
    width = ncol(enrmat)*unit(ux * 2, "mm"), 
    height = nrow(enrmat)*unit(ux, "mm"),
    border_gp=gpar(color='black', lwd=.5),
    col=colg,
    cell_fun = function(j, i, x, y, w, h, col){
        enr = enrmat[i,j]
        val = sub("^[ ]*","", formatC(round(enr,1), digits=2))
        if (enr != 0){
            grid.text(val, x, y, 
                gp=gpar(col=ifelse(enr < 5, 'black','white'), 
                    fontsize=gridtxt.fs))
        }
    }
)

col_fun = colorRamp2(c(0, max(sdf$logFC)), c("white", 'red'))
lfcmat = as.matrix(sdf[,'logFC',drop=F])
lfcmat[is.na(lfcmat)] = 0
rownames(lfcmat) = sdf$cat # paste0(sdf$var, ' - ', sdf$cat)
henr1 = rowAnnotation(top.genes = anno_text(sdf$cat, gp=gpar(fontsize=5)))
htlfc = Heatmap(lfcmat, 
    cluster_rows=FALSE,
    name='logFC',
    row_split=idsplit,
    width = ncol(enrmat)*unit(ux * 2, "mm"), 
    height = nrow(enrmat)*unit(ux, "mm"),
    border_gp=gpar(color='black', lwd=.5),
    col=col_fun
)


# Covariate enrichment heatmaps:
if (use.var == 'vc'){
    colsplit = sub("@.*","", colnames(cmat))
    colscale = (ncol(cmat) + length(unique(colsplit)) * 1/2) / (ncol(cmat))
} else {
    colsplit=NULL
    colscale=1 
}
htcovar = Heatmap(cmat, 
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    name='logFC',
    row_split=idsplit,
    column_split=colsplit,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(cmat)*unit(ux, "mm") * colscale, 
    height = nrow(cmat)*unit(ux, "mm"),
    col=cov.col_fun,
    # cell_fun = function(j, i, x, y, w, h, col){
    #     # Add the p-value text (placeholder for now)
    #     p = pmat[i,j]
    #     cr = abs(cmat[i,j])
    #     if (p >= 3){
    #         grid.text('*', x, y, 
    #             gp=gpar(col=ifelse(cr > 6, 'white','black'), 
    #                 fontsize=gridtxt.fs))}
    # }
)

# ht = hn1 + hg1 + htng + henr1 + htcovar + htenr
ht = hn1 + hg1 + htng + htcovar + htenr
h = 2 + 1 / 15 * nrow(sdf)
w = 5 + 1 / 10 * (5 + ncol(cmat))
w = ifelse(w > 25, 25, w)
pltprefix = paste0(imgpref, 'representative_panel')
saveHeatmap(ht, pltprefix, w=w, h=h)


