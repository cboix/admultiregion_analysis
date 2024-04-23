#!/usr/bin/R
# ---------------------------------------------------
# Plot contours around peaks for modules on the UMAP:
# Updated 01/19/2023 
# ---------------------------------------------------
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
library(soundgen)
library(scattermore) # For raster pts
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'modules/auxiliary_contour_functions.R'))


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
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id)}
source(paste0(sbindir, 'modules/load_modules_coords.R'))


# Calculate the overall number of cells + normalizing mat:
# --------------------------------------------------------
NBREAK = 500
df = data.frame(
    x=submeta$U1[ind],
    y=submeta$U2[ind],
    ind=ind)
df = df[df$x >= xlim[1] & df$x <= xlim[2],]
df = df[df$y >= ylim[1] & df$y <= ylim[2],]
df$bx = as.character(cut(df$x, breaks=NBREAK))
df$by = as.character(cut(df$y, breaks=NBREAK))
binx = sort(unique(df$bx))
biny = sort(unique(df$by))
xmap = sapply(binx, range.mean)
ymap = sapply(biny, range.mean)
df$xavg = xmap[df$bx]
df$yavg = ymap[df$by]
df$ncell = 1

kern = 25
kSD = 1
ncmat = smooth.mat('ncell', df, kern=kern, kSD=kSD)
mx = 0.025
ncmat[ncmat < mx & ncmat > 0] = mx


# Calculate scores for each of the modules separately:
# ----------------------------------------------------
mns = c( # A specific set
    9, 12, 24, 7, 19, # Main
    0, 17, 3, 8,  # Func1 
    30, 11, 2, 16, 28, 6, 13, 27 # Func2
)
mns = 1:30 - 1 
slist = list()
make.plot = FALSE
for (mn in mns){
    module = paste0('M',mn)
    print(module)
    x = log(scoremat[,module] + 1)
    df$z = x[df$ind]
    # Make the smoothed matrix:
    smat = smooth.mat('z', df, kern=kern, kSD=kSD)
    rn = rownames(smat)
    cn = colnames(smat)
    adjmat = smat / ncmat[rn, cn]
    adjmat[is.na(adjmat)] = 0
    adjmat[is.infinite(adjmat)] = 0
    # Add padding to the matrix for contour computation:
    padmat = pad.matrix(adjmat)
    x = as.numeric(rownames(padmat))
    y = as.numeric(colnames(padmat))
    # Get the contour lines:
    slist[[module]] = contourLines(x, y, padmat, levels=seq(0, 1, .1))
    if (make.plot){
        png('~/test.png', res=450, units='in', width=6, height=6)
        contour(x, y, padmat, levels=seq(0, 1, .1), xlim=xlim, ylim=ylim)
        dev.off()
    }
}


# Plot contour lines for multiple modules:
# ----------------------------------------
modcol = unique(nodedf[,c('leiden','col')])
mcolmap = modcol$col
names(mcolmap) = paste0('M', modcol$leiden)

# Aggregate all levels above an arbitrary cutoff:
NM = length(slist)
levdf = c()
for (module in names(slist)){
    cobj = slist[[module]]
    # lvs = sapply(cobj, function(x){x$level})
    # cutoff = max(lvs) - 0.05
    outobj = lapply(1:length(cobj), function(i){
        ll = cobj[[i]]
        if (ll$level >= 0.1){ 
            df = data.frame(x=ll$x, y=ll$y, level=ll$level, i=i)
            return(df) 
        }})
    mdf = do.call(rbind, outobj)
    if (!is.null(mdf)){
        mdf$module = module
        levdf = rbind(levdf, mdf)
    }
}
levdf$col = mcolmap[levdf$module]


# Plot overlaid:
# --------------
# Something weird about contour level cutoffs:
pltdf = levdf[levdf$level >= 0.39 & levdf$level <= 0.41,]
pltdf$mset = paste0(pltdf$module, '-', pltdf$i)
msets = unique(pltdf$mset)
mods = unique(pltdf$module)

# Plot the overlaid contours, with no internal color:
png('~/test.png', res=450, units='in', width=6, height=6)
par(mar=rep(0.1, 4))
plot(submeta$U1, submeta$U2, col='grey90', pch=19, cex=0.1, xlim=xlim, ylim=ylim, axes=F)
for (set in msets){
    subdf = pltdf[pltdf$mset == set,]
    polygon(x=subdf$x, y=subdf$y, col=NA, border=subdf$col[1])
}
box(lwd=.5)
legend('bottomleft', mods, col=mcolmap[mods], pch=15, cex=.75, bty='n', ncol=2)
dev.off()

# Plot the overlaid contours, with no transparent color:
png('~/test_tsp.png', res=450, units='in', width=6, height=6)
par(mar=rep(0.1, 4))
plot(submeta$U1, submeta$U2, col='grey90', pch=19, cex=0.1, xlim=xlim, ylim=ylim, axes=F)
for (set in msets){
    subdf = pltdf[pltdf$mset == set,]
    tcol = tsp.col(subdf$col[1], subdf$level[1])
    polygon(x=subdf$x, y=subdf$y, col=tcol, border=subdf$col[1])
}
box(lwd=.5)
legend('bottomleft', mods, col=mcolmap[mods], pch=15, cex=.75, bty='n', ncol=2)
dev.off()


# Plot selection, overlaid:
# -------------------------
# Something weird about contour level cutoffs:
cutoff = 0.4
pltdf = levdf[levdf$level >= cutoff - 0.1 & levdf$level <= cutoff + 0.01,]
id.modules = c(9,12,24,7, 19)
func.modules = c(0, 6, 13, 27, 16, 17, 3, 8, 28)
selmns = c(id.modules, func.modules)
idmod = paste0('M', id.modules)
pltdf = pltdf[pltdf$module %in% paste0('M', selmns),]
pltdf$mset = paste0(pltdf$module, '-', pltdf$i)
msets = unique(pltdf$mset)
mods = unique(pltdf$module)


# Plot smoothed + overlaid contours, with no transparent color:
# -------------------------------------------------------------
pdf('~/test_sel_smooth.pdf', width=4, height=4)
par(mar=rep(0, 4))
plot(submeta$U1, submeta$U2, xlim=xlim, ylim=ylim, axes=F, type='n')
for (set in msets){
    subdf = pltdf[pltdf$mset == set,]
    smat = smooth.contour(subdf$x, subdf$y,
        xmap, ymap, smoothing=0.1, min.radius=0.25)
    cxy = get.centroid(smat[,1], smat[,2])
    if (abs(cxy[3]) > .25){
        tcol = tsp.col(subdf$col[1], subdf$level[1])
        if (subdf$module[1] %in% idmod){
            lty = 'dashed'
            lwd.adj = .5
        } else {
            lty = 'dotted'
            lwd.adj = 0
        }
        lwd = lwd.adj + ifelse(subdf$level < 0.35, .75, 1.25)
        polygon(x=smat[,1], y=smat[,2], col=NA, border=subdf$col[1], lwd=lwd, lty=lty)
        text(cxy[1], cxy[2], subdf$module[1], col=subdf$col[1], font=2)
    }
}
box()
legend('bottomleft', mods, col=mcolmap[mods], pch=15, cex=.75, bty='n', ncol=2)
dev.off()

# Contours background:
png('~/test_background.png', res=450, units='in', width=4, height=4)
par(mar=rep(0, 4))
tc = tsp.col('grey85')
plot(submeta$U1, submeta$U2, col=tc, pch=19, cex=0.1, xlim=xlim, ylim=ylim, axes=F)
dev.off()



# Grid of modules:
# ----------------
# Plot each in own par:
pltdf = levdf[levdf$level >= 0.29 & levdf$level <= 0.81,]
pltdf$mset = paste0(pltdf$module, '-', pltdf$i)
msets = unique(pltdf$mset)
mods = unique(pltdf$module)

png('~/test_grid.png', res=450, units='in', width=8, height=8)
par(mfrow=c(5,5), mar=rep(0, 4))
for (mod in head(mods, 25)){
    subsets = unique(pltdf$mset[pltdf$module == mod])
    plot(submeta$U1, submeta$U2, col='grey90', pch=19, cex=0.1, xlim=xlim, ylim=ylim, axes=F)
    for (set in subsets){
        subdf = pltdf[pltdf$mset == set,]
        tcol = tsp.col(subdf$col[1], 0.1) # Add step-shade
        polygon(x=subdf$x, y=subdf$y, col=tcol, border=subdf$col[1])
        text(mean(xlim), mean(ylim), mod, font=2, adj=.5)
        box(lwd=.5)
    }
}
dev.off()


# Plot smoothed + overlaid contours, with no transparent color:
# -------------------------------------------------------------
pltdf = levdf[levdf$level >= 0.19 & levdf$level <= 0.81,]
pltdf$mset = paste0(pltdf$module, '-', pltdf$i)
msets = unique(pltdf$mset)
mods = unique(pltdf$module)

pdf('~/test_sel_smooth_grid.pdf', width=8, height=8)
par(mfrow=c(5,5), mar=rep(0, 4))
for (mod in head(mods, 25)){
    subsets = unique(pltdf$mset[pltdf$module == mod])
    plot(submeta$U1, submeta$U2, xlim=xlim, ylim=ylim, axes=F, type='n', ylab='', xlab='')
    for (set in subsets){
        subdf = pltdf[pltdf$mset == set,]
        smat = smooth.contour(subdf$x, subdf$y,
            xmap, ymap, smoothing=0.1, min.radius=0.25)
        cxy = get.centroid(smat[,1], smat[,2])
        if (abs(cxy[3]) > .25){
            tcol = tsp.col(subdf$col[1], 0.1) # Add step-shade
            polygon(x=subdf$x, y=subdf$y, col=tcol, border=subdf$col[1])
        }
    }
    text(mean(xlim), mean(ylim), mod, font=2, adj=.5)
    box(lwd=.5)
}
dev.off()


