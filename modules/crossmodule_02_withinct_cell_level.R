#!/usr/bin/R
# ------------------------------------------------------------------------
# Plot cross-module scores within the same celltype (at pseudobulk level):
# Updated 11/29/2021
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
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Set the run arguments:
# ----------------------
# Default arguments:
runset = 'Ast' 
graph_id = 'boot'

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

# Set a single color scale:
load.colors()
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))

colnames(scoremat) = sub("M","", colnames(scoremat))
rownames(cellmeta) = cellmeta$barcode
submeta = cellmeta[rownames(scoremat),]


# Calculate the subtype-agnostic cross-module correlation:
# TODO: Can also run glm models, account for projid, etc.
# --------------------------------------------------------
# Require that modules have 4+ genes in core cluster. 
ctdf = aggregate(gene ~ leiden, nodedf, length)
kept.modules = ctdf$leiden[ctdf$gene >= 4]

modcor = cor(scoremat[,as.character(kept.modules)])
rownames(modcor) = mmap$mname[as.numeric(rownames(modcor)) + 1]

plt = Heatmap(modcor, 
              col=col_fun,
              use_raster=TRUE,
              width = ncol(modcor)*unit(5, "mm"), 
              height = nrow(modcor)*unit(5, "mm"),
              border_gp = gpar(col="black", lty = 1))


h = 2.25 + 2.5 / 15 * nrow(modcor)
w = 5 + 2.5 / 15 * ncol(modcor)
pltprefix = paste0(imgpref, 'cell_level_module_correlations_', fullpref)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()


# Repeat, but just on the most common subtype:
# --------------------------------------------
ctdf = data.frame(table(submeta[,c('cell_type_high_resolution','region')]))
top.count = ctdf[which.max(ctdf$Freq),]
topct = top.count$cell_type_high_resolution
topreg = top.count$region

bcs = submeta[((submeta$cell_type_high_resolution == topct) & 
               (submeta$region == topreg)), 'barcode']

modcor = cor(scoremat[bcs,as.character(kept.modules)])
rownames(modcor) = mmap$mname[as.numeric(rownames(modcor)) + 1]

plt = Heatmap(modcor, 
              use_raster=TRUE,
              width = ncol(modcor)*unit(5, "mm"), 
              height = nrow(modcor)*unit(5, "mm"),
              border_gp = gpar(col="black", lty = 1))



# Run regression models instead to calculate module-module networks:
# TODO: several models (+/- covars, +/- Apoe_e4)
# Model accounting for celltype 
# Model accounting for celltype + region
# Model accounting for celltype + region + AD of one type (plaque, for example)
# ------------------------------------------------------------------
# Make a dataframe for regressions:
ad.cols = c('msex','pmi','age_death','Apoe_e4')
covar.cols = c('nrad','cogdxad')
regdf = submeta[,c('projid','region','cell_type_high_resolution','barcode')]
regdf = merge(regdf, unique(metadata[,c('projid','region','rind',covar.cols, ad.cols)]))
regdf = merge(regdf, pqdf, all.x=TRUE)
rownames(regdf) = regdf$barcode
regdf = regdf[rownames(scoremat),]


# Run single regression:
run.module.reg = function(i,j, reg.form, regdf){
    regdf$Mi = scoremat[,as.character(i)]
    regdf$Mj = scoremat[,as.character(j)]
    fit = lm(reg.form, regdf)
    cfit = coefficients(summary(fit))
    cfit = data.frame(cfit)
    names(cfit) = c('est','se','t','p')
    return(cfit['Mi',])
}


regresults.rda = paste0(crossdir, 'cell_level_regr_within_', fullpref, '.rda')
if (!file.exists(regresults.rda)){
    # Everything is significant - not very useful - but better than before.
    reg.form = asform(c("Mj ~ Mi", "+ cell_type_high_resolution * region"))
    # NOTE: THIS ONLY WORKS FOR CONSECUTIVE KEPT.MODULES (will break othw.)
    nm = length(kept.modules)
    resdf = c()
    for (i in kept.modules){
        for (j in kept.modules){
            cat(i,'\t', j,'\t')
            cfit = run.module.reg(i, j, reg.form, regdf)
            cfit$i = i
            cfit$j = j
            cat(round(cfit$est,2), '\n')
            resdf = rbind(resdf, cfit)
        }
    }

    # Turn this into a matrix:
    rmat = matrix(0, length(kept.modules), length(kept.modules))
    dfind = as.matrix(resdf[,c('i','j')]) + 1
    rmat[dfind] = resdf$est
    rownames(rmat) = kept.modules
    colnames(rmat) = kept.modules

    # Save:
    save(resdf, rmat, file=regresults.rda)
} else { load(regresults.rda) }


pltmat = rmat
# pltmat = (rmat + t(rmat)) / 2

pltmat = reord(pltmat)
pltmat = pltmat[,rownames(pltmat)]
rownames(pltmat) = mmap$mname[as.numeric(rownames(pltmat)) + 1]


plt = Heatmap(pltmat, 
              col=col_fun,
              use_raster=TRUE,
              cluster_rows=FALSE,
              cluster_columns=FALSE,
              width = ncol(modcor)*unit(5, "mm"), 
              height = nrow(modcor)*unit(5, "mm"),
              border_gp = gpar(col="black", lty = 1))


h = 2.25 + 2.5 / 15 * nrow(modcor)
w = 5 + 2.5 / 15 * ncol(modcor)
pltprefix = paste0(imgpref, 'cell_level_module_regressions_', fullpref)
pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
print(plt)
dev.off()
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=h)
print(plt)
dev.off()


# Plot the regression results as a graphical network:
# ---------------------------------------------------
library(igraph)

# Correlation based:
modcor = cor(scoremat[,as.character(kept.modules)])
modcor = data.frame(modcor)
modcor$i = rownames(modcor)

coldf = unique(nodedf[,c('leiden','col')])
coldf = coldf[order(coldf$leiden),]

# Using regressions:
# cutoff = .75
# edgedf = resdf[,c('i','j','est')]

# Using correlation:
cutoff = 0.4
cutoff = 0.35
edgedf = gather(modcor, j, est, -i)
edgedf$i = as.numeric(edgedf$i)
edgedf$j = sapply(edgedf$j, function(x){as.numeric(sub("X","",x))})

# Format and cutoff:
edgedf = edgedf[edgedf$i != edgedf$j,]
edgedf = edgedf[edgedf$est > cutoff,]
nodes = unique(c(edgedf$i, edgedf$j))
edge.score = edgedf$est

# Simple network: just the links/points:
net <- graph_from_data_frame(d=edgedf, vertices=nodes, directed=T) 
vcol = as.character(coldf[nodes + 1,'col'])
vcol[is.na(vcol)] = 'black'
# ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
ecol = 'grey'
V(net)$size = 7.5
V(net)$label = ""
# V(net)$label = nodes + 1
V(net)$color = vcol
V(net)$frame.color <- 'black' # vcol
V(net)$frame.color <- NA
V(net)$pch = 19
E(net)$color = ecol 
E(net)$arrow.size = .25
elty = rep('dotted', length(edgedf$est))
elty[edgedf$est >= .75] = 'dashed'
elty[edgedf$est >= 1] = 'solid'
E(net)$lty = elty
ewidth = ((edge.score >= .95) * .4 +
          (edge.score >= .85) * 0.4 +
          (edge.score >= 0.75) * 0.4) + 0.4
E(net)$width = edge.score
E(net)$weight = edge.score * .5

set.seed(1)
l <- layout_with_fr(net, grid='nogrid') # Usually best
lrange = apply(l, 2, range)
l2 = l
l2 = sweep(l2, 2, lrange[1,], '-')
l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

lbcex=0.5
labels = mmap[nodes + 1,'mname']
labels = sapply(labels, function(x){sub(": ","\n", x)})

pltprefix = paste0(imgpref, 'cell_level_module_regressions_network', fullpref)
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
par(yaxs='i',xaxs='i', mar=rep(.25,4))
plot(net, layout=l, curved=F)
rdf = cbrbase:::general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         labels=labels, cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex)
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25)
dev.off()






