#!/usr/bin/R
# ---------------------------------------------------
# Plot selected points on the glial (astrocyte) umap;
# Also plot module expression for these points
# Updated 01/25/2023 
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
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'module_contours_')
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
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id)}
source(paste0(sbindir, 'modules/load_modules_coords.R'))


ptfile = paste0(moddir, runset, '_selpts_ind.tsv')
if (!file.exists(ptfile)){
    # Interactively select points, do once.
    identifyPch <- function(x, y = NULL, n = length(x), plot = FALSE, pch = 19, ...)
    {
        xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
        sel <- rep(FALSE, length(x))
        while(sum(sel) < n) {
            ans <- identify(x[!sel], y[!sel], 
                labels = which(!sel), n = 1, plot = plot, ...)
            if(!length(ans)) break
            ans <- which(!sel)[ans]
            points(x[ans], y[ans], pch = pch)
            sel[ans] <- TRUE
        }
        ## return indices of selected points
        which(sel)
    }

    dev.new()
    if(dev.interactive()) { 
        x = submeta$U1
        y = submeta$U2
        plot(submeta$U1, submeta$U2, col='grey', pch=19, cex=0.1, xlim=xlim, ylim=ylim)
        ind = identifyPch(x,y, n=18) 
    }
    # 18 pts:
    # 3027   4746   7703  11039  18341  28143  34029  43079  56859  60067  74674  78360  95858 102703 114095
    # 86988  96790 110983
    write.table(ind, ptfile, quote=F, row.names=F, col.names=F, sep="\t")
} else {
    ind = as.numeric(scan(ptfile, 'c', quiet=T))
}


# Name the points:
# ----------------
NI = length(ind)
ind = sort(ind)
df = submeta[ind,c('U1', 'U2')]
df$ind = ind
df$ry = df$U2 < 10
df = df[order(df$ry, df$U1),]
df$i = 1:NI

# Push all away from center, easy:
mx = mean(xlim)
my = mean(ylim)
df$x = df$U1 + 0.5 * (df$U1 - mx)
df$y = df$U2 + 0.5 * (df$U2 - my)


# Coords for point choice:
# ------------------------
pdf('~/test_coords.pdf', width=4, height=4)
par(mar=rep(0, 4))
plot(submeta$U1, submeta$U2, xlim=xlim, ylim=ylim, type='n')
points(df$U1, df$U2, col='black', pch=19, cex=1, xlim=xlim, ylim=ylim)
# Add segments + text numbers
text(x=df$x, y=df$y, labels=df$i, srt=0, adj=0, xpd=TRUE, cex=1)
segments(df$x, df$y, df$U1, df$U2, lwd=.5)
box()
dev.off()


# Calculate neighborhoods:
# ------------------------
nbrs = list()
NTOP = 25
for (i in 1:NI){
    # Find the nearest N points to i
    j = df$ind[df$i == i]
    xi = submeta$U1[j]
    yi = submeta$U2[j]
    d = sqrt((submeta$U1 - xi)**2 + (submeta$U2 - yi)**2)
    d = head(order(d), NTOP)
    nbrs[[i]] = d
}


# Calculate scores:
# -----------------
smat = matrix(0, nrow=NI, ncol=ncol(scoremat))
colnames(smat) = colnames(scoremat)
for (i in 1:length(ind)){
    nind = nbrs[[i]]
    smat[i,] = apply(scoremat[nind,], 2, mean)
}

swide = data.frame(smat, pt=1:NI)
sdf = gather(swide, module, value, -pt)
sdf$pt = factor(sdf$pt, levels=1:NI)


# Plot scores barplot:
# --------------------
modcol = unique(nodedf[,c('leiden','col')])
mcolmap = modcol$col
names(mcolmap) = paste0('M', modcol$leiden)

id.modules = c(9, 12, 24, 7, 19)
func.modules = c(0, 13, 16, 17, 3, 8, 28)
selmns = paste0('M', c(id.modules, func.modules))
pltdf = sdf[sdf$module %in% selmns,]
mcolmap = mcolmap[unique(pltdf$module)]
pltdf$module = factor(pltdf$module, levels=selmns)

gp = ggplot(pltdf, aes(pt, value, fill=module)) + 
    facet_grid(module ~ ., scales='free_y') + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=mcolmap) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() +
    theme(legend.position='none') + 
    theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
pltprefix = paste0(imgpref, runset, '_selpts', NTOP, '_barplot')
saveGGplot(gp, pltprefix, w=4, h=6)

gp = ggplot(pltdf, aes(pt, value, fill=module)) + 
    facet_grid(. ~ module, scales='free_x') + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=mcolmap) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr() + coord_flip() + 
    theme(legend.position='none')
pltprefix = paste0(imgpref, runset, '_selpts', NTOP, '_barplot_horiz')
saveGGplot(gp, pltprefix, w=7, h=4)

