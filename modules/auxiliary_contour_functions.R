#!/usr/bin/R
# -----------------------------------------------
# Auxiliary functions for contour plots on UMAPs.
# Updated 01/23/2023 
# -----------------------------------------------
library(cbrbase)
library(tidyr)
library(soundgen)
options(width=170)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Calculate contours on a binned version of the surface:
# ------------------------------------------------------
range.mean = function(x){
    x = gsub("\\[|\\]|[()]","",x)
    x = strsplit(x, ',')[[1]]
    x = mean(as.numeric(x))
    return(x)
}

smooth.mat = function(var, df, kern=25, kSD=1, xb=NULL, yb=NULL){
    require(soundgen)
    aggdf = aggregate(asform(c(var, '~ xavg + yavg')), df, mean)
    if (is.null(xb) | is.null(yb)){
        xavg = unique(aggdf$xavg)
        yavg = unique(aggdf$yavg)
    } else {
        xavg = (xb[-1] + xb[-length(xb)]) / 2
        yavg = (yb[-1] + yb[-length(yb)]) / 2
    }
    combdf = expand.grid(xavg=xavg, yavg=yavg)
    aggdf = merge(aggdf, combdf, all.y=TRUE)
    aggdf[is.na(aggdf[[var]]),var] = 0
    aggdf = aggdf[order(aggdf$xavg),]
    aggdf = aggdf[order(aggdf$yavg),]
    # Reshape for contour function:
    aggwide = spread(aggdf, yavg, var)
    zmat = as.matrix(aggwide[,-1])
    na.mat = is.na(zmat)
    zmat[is.na(zmat)] = 0
    smat = gaussianSmooth2D(zmat, kernelSize=kern, kernelSD=kSD)
    rownames(smat) = aggwide$xavg
    return(smat)
}

pad.matrix = function(mat, n=1){
    x = as.numeric(rownames(mat))
    y = as.numeric(colnames(mat))
    dx = median(diff(x))
    dy = median(diff(y))
    x2 = c(x[1] - rev(1:n*dx), x,
        x[length(x)] + 1:n*dx)
    y2 = c(y[1] - rev(1:n*dy), y,
        y[length(y)] + 1:n*dy)
    xc = as.character(x)
    yc = as.character(y)
    xc2 = as.character(x2)
    yc2 = as.character(y2)
    nx = length(xc2)
    ny = length(yc2)
    outmat = matrix(0, nrow=nx, ncol=ny, dimnames=list(xc2, yc2))
    outmat[rownames(mat), colnames(mat)] <- mat
    return(outmat)
}

# Returns centroid and signed area
get.centroid = function(x, y){
    n = length(x)
    x = c(x, x[1])
    y = c(y, y[1])
    cx = 0
    cy = 0
    a = 0
    for (i in 1:n){
        ai = (x[i] * y[i+1] - x[i+1] * y[i])
        a = a + ai
        cx = cx + (x[i] + x[i+1]) * ai
        cy = cy + (y[i] + y[i+1]) * ai
    }
    a = a / 2
    cx = cx / (6 * a)
    cy = cy / (6 * a)
    return(c(cx, cy, a))
}


# Smooth the contours:
# --------------------
# Cribbed from obreschkow/cooltools
# xavg, yavg, points in grid for bins:
smooth.contour = function(x, y, xavg, yavg, smoothing=0.5, min.radius=1){
    # Difference for x, y:
    nx = length(xavg)
    ny = length(yavg)
    dx = (max(xavg) - min(xavg)) / nx
    dy = (max(yavg) - min(yavg)) / ny
    s = 1 - (1 - smoothing)^2
    n = length(x)
    radius = sqrt(sd(x)^2/dx^2 + sd(y)^2/dy^2)
    if (n>4 & radius>min.radius) {
        is.closed.curve = abs(x[1]-x[n])<stats::sd(x) & abs(y[1]-y[n])<stats::sd(y)
        interval = max(sd(x)/dx, sd(y)/dy)
        if (is.closed.curve | interval>min.radius) {
            if (is.closed.curve) {
                df = round(5*((1-s)*n + s*4))
                sx = smooth.spline(seq(5*n), rep(x,5), df=df)$y[c((2*n+1):(3*n),2*n+1)]
                sy = smooth.spline(seq(5*n), rep(y,5), df=df)$y[c((2*n+1):(3*n),2*n+1)]
            } else {
                w = rep(1,n)
                w[1] = w[n] = n
                df = round((1-s)*n+s*4)
                sx = smooth.spline(seq(n),x,w=w,df=df)$y
                sy = smooth.spline(seq(n),y,w=w,df=df)$y
            }
        } else {
            print("Cannot smooth; not closed or interval too small")
            sx = x
            sy = x
        }
    } else {
        print("Cannot smooth; few points or small radius")
        sx = x
        sy = x
    }
    return(cbind(sx, sy))
}

