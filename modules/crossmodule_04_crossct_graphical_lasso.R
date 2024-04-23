#!/usr/bin/R
# -----------------------------------------------------
# Calculate the cross-ct module graphical lasso models:
# Updated 11/30/2021
# -----------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(glasso)
library(cglasso)  # Conditional lasso

# Directories:
moddir = paste0(sdbdir, 'modules/')
crossdir = paste0(sdbdir, 'crossmodule/')
plotdir = paste0(imgdir, 'crossmodule/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, crossdir, moddir)
system(cmd)


# Functions for matrices + networks:
# ----------------------------------
source(paste0(sbindir, 'modules/auxiliary_crossmodule_plotting_fns.R'))


# Load in the cross module pseudobulk-level data:
# -----------------------------------------------
# TODO: Opt to include neurons or not? Or just add later...
# commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
source(paste0(sbindir, 'modules/load_crossmodule_psbulk.R'))


# Turn the score dataframe into a matrix:
# NOTE: At the runset level or at the subtype level:
# --------------------------------------------------
mingenes = 10
modlevel = 'subtype'
modlevel = 'runset'
if (modlevel == 'subtype'){
    mat = pivot.tomatrix(scoredf[scoredf$ng >= mingenes, c('pr','cm','score')], 'pr','score')
    scmat = log(t(mat))
    ut = 1
} else {
    mat = pivot.tomatrix(runscdf[runscdf$ng >= mingenes, c('pr','rm','score')], 'pr','score')
    scmat = log(t(mat))
    ut = 3
}





if (modlevel == 'subtype'){
    # (scoredf[scoredf$ng >= mingenes, c('pr','cm','score')], 'pr','score')
    # subdf = runscdf[runscdf$ng >= mingenes,]
    sregion = sapply(rownames(scmat), function(x){sub(".*-","",x)})
    smetadf = data.frame(region=sregion)
    mdx = model.matrix(~region, smetadf)

    cgmat = scmat
    colnames(cgmat) = NULL
    sdata = datacggm(Y=scmat, X = smetadf)
    out <- cglasso(. ~ ., data = sdata)
    out

    plt = plotSymMat(out$Yipt, col_fun=col_fun)


    lambda.new <- mean(out$lambda)
    rho.new <- mean(out$rho)
    # imputing missing values
    Y.impute <- impute(out, type = "mar", lambda.new = lambda.new, rho.new = rho.new)

    # o2 = fitted(out, lambda.id = 3L, rho.id = 3L, drop = TRUE)






# Model 1: censored glasso estimator (Augugliaro \emph{and other}, 2020a)
# Y ~ N(0, Sigma) and probability of left/right censored values equal to 0.05
n <- 1000L
p <- 3L
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out

# Model 2: conditional censored glasso estimator (Augugliaro \emph{and other}, 2020b)
# Y ~ N(b0 + XB, Sigma) and probability of left/right censored values equal to 0.05
n <- 1000L
p <- 3L
q <- 2L
b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probl = 0.05, probr = 0.05)
out <- cglasso(. ~ ., data = Z)
out

# Graphical Lasso experiments:
# ----------------------------
modcov = cov(scmat, use='pairwise.complete.obs')
modcov[is.na(modcov)] = 0


rho=0.025
fit = glasso(modcov, rho=rho, nobs = nrow(scmat))
rownames(fit$w) = rownames(modcov)
rownames(fit$wi) = rownames(modcov)
colnames(fit$wi) = colnames(modcov)

plt = plotSymMat(fit$wi, col_fun=col_fun)
wdf = edgesFromMat(fit$wi)

wdf = edgesFromMat(fit$wi)
edgedf = wdf[wdf$est > 0,]
ll = prune.edges(edgedf, modlevel)


# Make a network and plot:
# ------------------------
netlist = make.network(ll$edgedf[,c('i','j','est')], 'est', ll$mndf)
netlist = set.network.params(netlist)
V(netlist$net)$size = 3

pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_graphical_lasso_network')
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.2, spacing=.25, adjust=FALSE)
dev.off()


# Plot the cross-subtype interactions only:
# -----------------------------------------
ind = substr(ll$edgedf$M1,1,3) != substr(ll$edgedf$M2,1,3) 
netlist = make.network(ll$edgedf[ind,c('i','j','est')], 'est', ll$mndf)
netlist = set.network.params(netlist, seed=0)
V(netlist$net)$size = 3

pltprefix = paste0(imgpref, 'pseudobulk_crossctonly_', 
                   modlevel, '_graphical_lasso_network')
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.2, spacing=.25, adjust=FALSE)
dev.off()



# Try regressing signal and then correlation: 
# --------------------------------------------
mdf = data.frame(pr=rownames(scmat),
                 projid = sapply(rownames(scmat), function(x){sub("-.*","",x)}),
                 region = sapply(rownames(scmat), function(x){sub(".*-","",x)}))
mdf = merge(mdf, metadata[,c('projid','region','rind','msex','age_death', 'cogdxad','nrad')])
rownames(mdf) = mdf$pr
mdf = mdf[rownames(scmat),]

scres = scmat * 0
for (i in 1:ncol(scmat)){
    mdf$y = scmat[,i]
    fit = glm(y ~ region + msex + age_death, mdf, family='gaussian')
    pred = predict(fit, mdf)
    scres[,i] = scmat[,i] - pred
}

# Calculate and threshold the correlation matrix:
# -----------------------------------------------
modcor = cor(scres, use='pairwise.complete.obs')
modcor[is.na(modcor)] = 0
plt = plotSymMat(modcor, col_fun=col_fun, ut=1.5)
wdf = edgesFromMat(modcor)

# n = 48 # Effective N may be more between nrow(scres) and 48.
n = nrow(scres)
wdf$p = (1 - wdf$est^2)^(n/2 - 2) / beta(1/2, n/2 - 1)
wdf$p.adj = p.adjust(wdf$p,'fdr')
cutoff = min(wdf$est[wdf$p.adj < 0.01])
plt = plotSymMat(modcor * (modcor >= cutoff), col_fun=col_fun, ut=ut)


# Make edgelist and set list:
# ---------------------------
edgedf = wdf[wdf$est > cutoff,]
ll = prune.edges(edgedf, modlevel)


# Make a network and plot:
# ------------------------
if (modlevel == 'subtype'){cutoff = 0.7}
netlist = make.network(ll$edgedf[,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff)
netlist = set.network.params(netlist)
V(netlist$net)$size = 3

pltprefix = paste0(imgpref, 'pseudobulk_crossct_',
                   modlevel, '_correlation_resid_network')
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.2, spacing=.25, adjust=FALSE)
dev.off()


# Cross-subtype only:
ind = substr(ll$edgedf$M1,1,3) != substr(ll$edgedf$M2,1,3) 
netlist = make.network(ll$edgedf[ind,c('i','j','est')], 'est', ll$mndf, cutoff=cutoff)
netlist = set.network.params(netlist, seed=0)
V(netlist$net)$size = 3

pltprefix = paste0(imgpref, 'pseudobulk_crossctonly_', 
                   modlevel, '_correlation_resid_network')
w = 4
png(paste0(pltprefix, '.png'), res=450, units='in', width=w, height=w)
plot.network(netlist, lbcex=0.2, spacing=.25, adjust=FALSE)
dev.off()

