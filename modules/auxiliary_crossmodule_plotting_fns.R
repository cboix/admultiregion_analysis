#!/usr/bin/R
# ---------------------------------------------------------
# Auxiliary functions for cross-module matrices / networks:
# Updated 11/30/2021
# ---------------------------------------------------------
library(cbrbase)

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(igraph)


# Functions for matrices:
# -----------------------
plotSymMat = function(mat, col_fun=col_fun, ut=3, cluster=TRUE, raster=TRUE, shownames=TRUE, ra=NULL, ta=NULL){
    plt = Heatmap(mat, 
        col=col_fun,
        cluster_columns=cluster,
        cluster_rows=cluster,
        use_raster=raster,
        show_column_names=shownames,
        show_row_names=shownames,
        right_annotation=ra,
        top_annotation=ta,
        width = ncol(mat)*unit(ut, "mm"), 
        height = nrow(mat)*unit(ut, "mm"),
        border_gp = gpar(col="black", lty = 1))
    return(plt)
}


edgesFromMat = function(mat, return.negative=FALSE){
    # Make edgelist from graphical lasso:
    colnames(mat) = sapply(colnames(mat), function(x){gsub(" ","__",x)})
    edgedf = data.frame(mat)
    edgedf$M1 = rownames(edgedf)
    edgedf = gather(edgedf, M2, est, -M1)
    edgedf$M2 = sapply(edgedf$M2, function(x){gsub("[.]","-",x)})
    edgedf$M2 = sapply(edgedf$M2, function(x){gsub("__"," ",x)})
    nodes = unique(c(edgedf$M1, edgedf$M2))
    # TODO: keep only single direction (option)
    edgedf = edgedf[edgedf$M1 != edgedf$M2,]
    edgedf = edgedf[order(edgedf$est, decreasing=T),]
    if (return.negative){
        edgedf = edgedf[edgedf$est != 0,]
    } else {
        edgedf = edgedf[edgedf$est > 0,]
    }
    return(edgedf)
}

col_fun = colorRamp2(c(-1,0,1), c("blue", "white", "red"))


# Network creation and plotting:
# ------------------------------
make.network = function(edgedf, score, mndf, cutoff=NULL, directed=FALSE, label.edge=FALSE, symmetric=TRUE){
    require(igraph)
    # Format and cutoff:
    edgedf = edgedf[edgedf$i != edgedf$j,]
    if (!(directed) & (symmetric)){
        edgedf = edgedf[edgedf$i > edgedf$j,]
    }
    if (!is.null(cutoff)){
        edgedf = edgedf[edgedf[[score]] > cutoff,]
    }
    edge.score = edgedf[[score]]
    nodes = mndf$ind

    # Use mndf as a mapping:
    labs = as.character(mndf[nodes, 'name'])
    vcol = as.character(mndf[nodes,'col'])
    vcol[is.na(vcol)] = 'black'

    # Simple network: just the links/points:
    net <- graph_from_data_frame(
        d=edgedf[,c('i','j',score)], 
        vertices=nodes, directed=directed) 
    net = simplify(net)

    # Main attr:
    V(net)$color = vcol
    ewidth = ((edge.score >= .95) * .4 +
              (edge.score >= .85) * 0.4 +
              (edge.score >= 0.75) * 0.4) + 0.4
    E(net)$width = edge.score * .5 + .5
    E(net)$weight = edge.score * .5

    if (label.edge){
        E(net)$label = round(edge.score, 1)
        E(net)$label.size = 1
        E(net)$label.color = 'grey50'
    }

    return(list(net=net, nodes=nodes,
                labels=labs, edgedf=edgedf))
}


# Use object from before:
set.network.params = function(netlist, seed=1, vcol=NA, use.lty=TRUE){
    require(igraph)
    net = netlist$net
    edgedf = netlist$edgedf

    ecol = 'grey75'
    V(net)$size = 7.5
    V(net)$label = ""
    V(net)$frame.color <- vcol
    V(net)$pch = 19
    E(net)$color = ecol 
    E(net)$arrow.size = .25
    if (use.lty){
        elty = rep('dotted', length(edgedf$est))
        elty[edgedf$est >= .5] = 'dashed'
        elty[edgedf$est >= .75] = 'solid'
        E(net)$lty = elty
    }

    # Get layout:
    set.seed(seed)
    l <- layout_with_fr(net, grid='nogrid') # Usually best
    lrange = apply(l, 2, range)
    l2 = l
    l2 = sweep(l2, 2, lrange[1,], '-')
    l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

    # Return:
    netlist$net = net
    netlist$layout = l
    netlist$layout.norm = l2
    return(netlist)
}


# Plotting the network w/adjusted (or not) labels:
plot.network = function(netlist, lbcex=0.5, spacing=.25, adjust=FALSE, 
                        pie.values=NULL, pie.cols=NULL){
    require(igraph)
    labels = sapply(netlist$labels, function(x){sub(": ","\n", x)})
    par(yaxs='i',xaxs='i', mar=rep(spacing,4))
    if (is.null(pie.values)){
        plot(netlist$net, layout=netlist$layout, curved=F)
    } else {
        plot(netlist$net, 
             layout=netlist$layout, curved=F,
             vertex.shape="pie", 
             vertex.pie=pie.values, 
             vertex.pie.color=list(pie.cols),
             vertex.pie.border='black',
             vertex.pie.lty=1)  #, vertex.size=maxv)
    }
    if (adjust){
        rdf = cbrbase:::general_repel_text(x=netlist$layout.norm[,1], y=netlist$layout.norm[,2], 
                                        xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                                        hjust=.5, vjust=.5, seed=1, max.iter=5000,
                                        labels=labels, cex=lbcex, pt.cex=.25)
        text(x=rdf$x, y=rdf$y, labels=rdf$lab, srt=0, adj=0, xpd=TRUE, cex=lbcex)
        segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25)
    } else {
        text(x=netlist$layout.norm[,1], y=netlist$layout.norm[,2], 
             labels=labels, srt=0, adj=0, xpd=TRUE, cex=lbcex)
    }
}


# Network clustering:
