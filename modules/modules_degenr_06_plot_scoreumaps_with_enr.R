#!/usr/bin/R
# ---------------------------------------------------------------
# Plot module umaps and with the functional enrichment resources:
# Updated 03/18/2022 to standardize sizes
# ---------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ggpubr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Match metadata to UMAP coordinates:
# -----------------------------------
if (runset %in% c('Exc', 'Inh', names(exc.sets))){
    submeta = read.delim(ctumap.file, sep="\t")
    submeta = submeta[, c('barcode','C1','C2')]
    names(submeta) = c('barcode','U1','U2')
    rownames(cellmeta) = cellmeta$barcode
    submeta$region = cellmeta[submeta$barcode, 'region']
    submeta$projid = cellmeta[submeta$barcode, 'projid']
    submeta$cell_type_high_resolution = cellmeta[submeta$barcode, 'cell_type_high_resolution']
    rownames(submeta) = submeta$barcode
    submeta = submeta[rownames(scoremat),]
} else if (runset == 'Vasc_Epithelia'){
    submeta = read.delim(ctumap.file)
    rownames(submeta) = submeta$barcode
    submeta = submeta[rownames(scoremat),c('barcode','region','celltype','C1','C2')]
    names(submeta)[3:5] = c('cell_type_high_resolution', 'U1','U2')
} else {
    rownames(cellmeta) = cellmeta$barcode
    submeta = cellmeta[rownames(scoremat),]
}


# Load in set of top by p-value functional enrichments for each module:
# On core genes, with term_size < 500 (in this case: "small")
# ---------------------------------------------------------------------
# if (runset %in% c('Mic_Immune','Vasc_Epithelia', 'Exc', 'Inh', names(exc.sets))){
# useset = 'allgenes_'
# } else { useset = 'deonly_' }
useset = 'coregenes_'
toppvals.file = paste0(moddir, 'module_enrichments_toppvals_',
                       useset, fullpref, '_small.tsv')
pvalsdf = read.delim(toppvals.file, sep="\t")
pvalsdf$nc = nchar(pvalsdf$term) # For filtering:


# Plotting function for each of the module scores:
# ------------------------------------------------
ind = 1:nrow(submeta)
cex = 0.025

# Table of params for each celltype to center + not change aspect ratio:
center.params = list("Ast"=c(0, 8.5, -8.5, 0),
                     "Mic_Immune"=c(1.25, 9.25, -12.5, -4),
                     "Opc"=c(0, 6.5, -8.5, 0),
                     "Oli"=c(-1,11,-14,-2),
                     "HCneurons"=c(14,32,-18,0))

# Set the x and y limits using these:
if (runset %in% names(center.params)){
    cp = center.params[[runset]]
    xlim = c(min(submeta$U1) + cp[1], min(submeta$U1) + cp[2])
    ylim = c(max(submeta$U2) + cp[3], max(submeta$U2) + cp[4])
} else {
    xlim = range(submeta$U1) 
    ylim = range(submeta$U2)
    r = max(c(diff(xlim), diff(ylim))) / 2
    xlim = c(mean(xlim) - r, mean(xlim) + r)
    ylim = c(mean(ylim) - r, mean(ylim) + r)
}

palette = viridis(100)
col_fun = function(x, pal=palette){
    bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
    palette[bin]  }

plot.module.score = function(m, ind, pch=19, textsize=1, cex=0.025){
    module = paste0('M',m)
    x = log(scoremat[,module] + 1)
    sp = 0.1
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    plot(submeta$U1[ind], submeta$U2[ind], type='n', ylim=ylim, xlim=xlim, axes=F)
    points(submeta$U1[ind], submeta$U2[ind], col=col_fun(x[ind]), pch=pch, cex=cex)
    # mtext(module, side=1, line=-1, font=2, cex=textsize)
}


# Select top 30 modules by enrichment p-values:
# Also require that they have 4+ genes in core cluster. 
# -----------------------------------------------------
ctdf = aggregate(gene ~ leiden, nodedf, length)
kept.modules = ctdf$leiden[ctdf$gene >= 4]

pdf = aggregate(p ~ mname, pvalsdf, min)
pdf = merge(pdf, mmap)
pdf = pdf[pdf$module %in% kept.modules,]
pdf = pdf[order(pdf$p),]
plt.modules = head(pdf$module, 30)

# Fill in if needed + if possible:
minmod = min(30, length(kept.modules))
while(length(plt.modules) < minmod){
    notin = kept.modules[!(kept.modules %in% plt.modules)]
    plt.modules = c(plt.modules, head(sort(notin),1))
}
plt.modules = sort(plt.modules)


# Plot each of the module scores + text:
# --------------------------------------
nr = 10
nc = 3
wone = 1 / 2.54  # 1 cm sized
modmat = matrix(plt.modules, nrow=nr, ncol=nc, byrow=F)
sp = 0.1
ntop = 5

# Make resources separately, so we can combine later for full figure:
last = -1
for (icol in 1:nc){
    print(icol)
    colmodules = modmat[,icol]
    png(paste0(imgpref, 'modules_umap_resource_', runset, '_c', icol, '.png'), units='in', res=300, width=wone, height=wone * nr)
    layout(matrix(1:(1 * nr), nr, 1))
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    for (i in colmodules){ 
        if (i > last){
            plot.module.score(m=i, ind=ind, pch='.', textsize=0.75) 
            last = i
        }
    }
    dev.off()
}

# Make a full UMAP panel:
png(paste0(imgpref, 'modules_umap_full_resource_', runset, '.png'), units='in', res=300, width=wone * 4.25 * nc, height=wone * nr)
sp = 0.1
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
layout(matrix(1:(nr*nc * 2), nr, nc * 2, byrow=FALSE), widths=rep(c(1, 3.25), nc))
last = -1
for (icol in 1:nc){
    colmodules = modmat[,icol]
    for (i in colmodules){ 
        if (i > last){
            plot.module.score(m=i, ind=ind, pch='.', textsize=0.75) 
            last = i
        }
    }
    for (i in 1:nr){ plot(1,1, type='n', axes=F) }
}
dev.off()


# Convert enrichments text with a table of common abbreviations:
abdf = read.delim(paste0(dbdir, 'enrichments_abbreviation_table.tsv'), header=T)
abdf$repl = paste0(abdf$abbv, '.')
pvalsdf$abterm = pvalsdf$term
for (i in 1:nrow(abdf)){
    pvalsdf$abterm = sub(abdf$word[i], abdf$repl[i], pvalsdf$abterm)
}


# Make the text enrichment:
sp = 0.1
pdf(paste0(imgpref, 'modules_enrichment_text_resource_', runset, '.pdf'), width=wone * 4.25 * nc, height=wone * nr)
par(xaxs='i', yaxs='i', mar=rep(sp, 4), lheight=0.8, xpd=NA)
layout(matrix(1:(nr*nc), nr, nc, byrow=FALSE))
for (i in plt.modules){
    plot(c(0,1),c(0,1), type='n', axes=F, xlab='', ylab='')
    mn = mmap[i+1,'mname']
    subdf = pvalsdf[(pvalsdf$mname == mn) & (pvalsdf$nc < 40),, drop=F]
    text(1 / 4.25, 1, mn, adj=0, cex=.5, font=2)
    if (nrow(subdf) > 0){
        subdf = head(subdf[order(subdf$p),], ntop)
        pstr = paste0(subdf$abterm, ' (', sprintf('%0.1e', subdf$p),')')
        textstr = paste(pstr, collapse="\n")
        text(1 / 4.25, 0.55, textstr, adj=0, cex=.5)
    }
}
dev.off()



# Plot the region and cell types as resources:
# --------------------------------------------
tsp.tcols = sapply(tcols, tsp.col)
tsp.rcols = sapply(reg.cols, tsp.col)

w = wone
sp = 0.1
png(paste0(imgpref, 'modules_umap_ct_resource_', runset, '.png'), units='in', res=450, width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
plot(submeta$U1[ind], submeta$U2[ind], col=tsp.tcols[submeta$cell_type_high_resolution[ind]], 
    ylim=ylim, xlim=xlim, pch='.', cex=cex, axes=F)
dev.off()

if (runset == 'All'){
    png(paste0(imgpref, 'modules_umap_mct_resource_', runset, '.png'), units='in', res=450, width=w, height=w)
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    plot(submeta$U1[ind], submeta$U2[ind], 
        col=major.col[submeta$major.celltype[ind]], 
        ylim=ylim, xlim=xlim, pch='.', cex=cex, axes=F)
    dev.off()
}


png(paste0(imgpref, 'modules_umap_reg_resource_', runset, '.png'), units='in', res=450, width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
plot(submeta$U1[ind], submeta$U2[ind], col=tsp.rcols[submeta$region[ind]],
    ylim=ylim, xlim=xlim, pch='.', cex=cex, axes=F)
dev.off()


plt.subtypes = unique(submeta$cell_type_high_resolution)
subtypes = plt.subtypes
w = 1.75
pdf(paste0(imgpref, 'resource_', runset, '_legend.pdf'), width=w, height=w)
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
plot(1,1, xlim=c(0,1), ylim=c(0,1), type='n', axes=F)
legend('topright', legend=subtypes, col=tcols[plt.subtypes], pch=15, cex=.5)
dev.off()


