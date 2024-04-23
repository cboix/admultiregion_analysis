#!/usr/bin/R
# ---------------------------------------------------------------
# Plot module umaps and functional enrichment as indpt resources:
# Updated 02/16/2022
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
plotdir = paste0(imgdir, 'modules/pltres/')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

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

imgpref = paste0(plotdir, 'modules_indpt_res_', runset, '_')


# Load in and process data and UMAP coordinates:
# ----------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


commandArgs <- function(trailingOnly=TRUE){c(runset)}
source(paste0(sbindir, 'modules/load_modules_umap_coords.R'))


# Load in set of top by p-value functional enrichments for each module:
# ---------------------------------------------------------------------
set = 'coregenes'
toppvals.file = paste0(moddir, 'module_enrichments_toppvals_',
                       set, '_', fullpref, '_small.tsv')
pvalsdf = read.delim(toppvals.file, sep="\t")
pvalsdf$nc = nchar(pvalsdf$term) # For filtering long terms

# Select modules with at least 10 genes:
ctdf = aggregate(gene ~ leiden, nodedf, length)
plt.modules = ctdf$leiden[ctdf$gene >= 10]


# Plot each of the module scores + text:
# --------------------------------------
palette = viridis(100)
palette = c('grey85',col2)
colgy <- colorRampPalette(brewer.pal(n=7, name="Greys"))(100)
palette = colgy
col_fun = function(x, pal=palette){
    bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
    palette[bin]  }

plot.module.score = function(m, ind, pch=19, textsize=1){
    module = paste0('M',m)
    x = log(scoremat[,module] + 1)
    sp = 0.1
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    plot(submeta$U1[ind], submeta$U2[ind], type='n', ylim=ylim, xlim=xlim, axes=F)
    points(submeta$U1[ind], submeta$U2[ind], col=col_fun(x[ind]), pch=pch, cex=cex)
    if (!is.na(textsize)){
        mtext(module, side=1, line=-1, font=2, cex=textsize)
    }
}

last = -1
wone = 0.6
sp = 0.05
for (i in plt.modules){
    png(paste0(imgpref, 'umap_M', i, '.png'), units='in', res=300, width=wone, height=wone)
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    plot.module.score(m=i, ind=ind, pch='.', textsize=NA) 
    last = i
    dev.off()
}


ntop = 5
for (i in plt.modules){
    pdf(paste0(imgpref, 'enrtext_M', i, '.pdf'), width=1.5, height=.75)
    par(xaxs='i', yaxs='i', mar=rep(sp, 4))
    plot(c(0,1),c(0,1), type='n', axes=F, xlab='', ylab='')
    mn = mmap[i+1,'mname']
    subdf = pvalsdf[(pvalsdf$mname == mn) & (pvalsdf$nc < 40),, drop=F]
    if (nrow(subdf) > 0){
        subdf = head(subdf[order(subdf$p),], ntop)
        pstr = paste0(subdf$term, ' (', sprintf('%0.1e', subdf$p),')')
        textstr = paste(c(mn, pstr), collapse="\n")
        text(0.5, 0.5, textstr, cex=.45)
    }
    dev.off()
}


# Plot specific modules as grid:
# Plot each of the module scores + text:
# --------------------------------------
if (runset == 'Ast'){
    plt.modules = c(12, 24, 7, 19, 9, 17, 3, 8)
    nr = 4
    nc = 2
    modmat = matrix(plt.modules, nrow=nr, ncol=nc, byrow=F)
}

sp = 0.1
ntop = 5

# Make resources separately, so we can combine later for full figure:
last = -1
for (icol in 1:nc){
    print(icol)
    colmodules = modmat[,icol]
    png(paste0(imgpref, 'umap_selmodules.png'), units='in', res=300, width=wone, height=wone * nr)
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


pdf(paste0(imgpref, 'enrtext_selmodules.pdf'), width=1.5 * nc, height=.75*nr)
par(xaxs='i', yaxs='i', mar=rep(sp, 4))
layout(matrix(1:(nr*nc), nr, nc, byrow=FALSE))
for (i in plt.modules){
    plot(c(0,1),c(0,1), type='n', axes=F, xlab='', ylab='')
    mn = mmap[i+1,'mname']
    subdf = pvalsdf[(pvalsdf$mname == mn) & (pvalsdf$nc < 40),, drop=F]
    textstr = mn
    if (nrow(subdf) > 0){
        subdf = head(subdf[order(subdf$p),], ntop)
        pstr = paste0(subdf$term, ' (', sprintf('%0.1e', subdf$p),')')
        textstr = paste(c(textstr, pstr), collapse="\n")
    }
    text(0.5, 0.5, textstr, cex=.45)
}
dev.off()

