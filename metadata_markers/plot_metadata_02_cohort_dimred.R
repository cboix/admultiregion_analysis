#!/usr/bin/R
# ----------------------------------------------------------
# Plot the full ROSMAP cohort + our data in the context of dim. red.
# Updated 02/18/2021
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(qvalue)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# ---------------------------------
# Process the full cohort metadata:
# ---------------------------------
pids = unique(sort(celldf$projid))
# Load in all (3642 x 126 vars)
emeta = read.delim('Annotation/metadata_PFC_all_individuals_092520.tsv', header=T)
# Variables we care about (pathology + cognition):
pathreg = c('ag','ec','hip','mf','mt')
pvars = c('nft','plaq_n','plaq_d',
          paste0("nft_", pathreg),
          paste0("plaq_n_", pathreg),
          paste0("plaq_d_", pathreg),
          colnames(emeta)[grep('^cogn',colnames(emeta))])

# Keep only ones with pathology measurements (1506)
pmeta = emeta[,c('projid',pvars)]
pmeta = pmeta[!is.na(apply(pmeta,1, sum)),]
rownames(emeta) = as.character(emeta$projid)
sc = emeta[as.character(pmeta$projid),'niareagansc']

# Now add attributes:
nums <- unlist(lapply(emeta, is.numeric))
nums['projid'] = FALSE
pmat = as.matrix(emeta[as.character(pmeta$projid),nums])
# pmat = as.matrix(pmeta[,pvars]) # Only path vars:

# Fill NA with median:
for (var in colnames(pmat)){
    if (sum(is.na(pmat[,var])) > 0){
        pmat[,var][is.na(pmat[,var])] <- median(pmat[,var], na.rm=TRUE)
    }
}

# Perform a preliminary dimensionality reduction:
library(compositions)
library(factoextra)
library(uwot)
res.pca <- prcomp(pmat, scale = TRUE)

scale = 1.5
png(paste0(imgpref, 'pathvars_pca_fullcohort_niareagansc.png'), res=450, units='in', width=6*scale, height=5.5*scale)
fviz_pca_biplot(res.pca, repel = TRUE, pch=19,
                geom="point",
                col.ind = as.character(sc),
                cex = .1,
                palette = sapply(colvals[['niareagansc']], tsp.col),
                col.var = "grey50", # Variables color
                addEllipses = TRUE, # Concentration ellipses
                ellipse.type = "confidence"
)
dev.off()

scale = 1.5
png(paste0(imgpref, 'pathvars_pca_fullcohort_incohort.png'), res=450, units='in', width=6*scale, height=5.5*scale)
sc = pmeta$projid %in% pids
fviz_pca_biplot(res.pca, repel = TRUE, pch=19,
                geom="point",
                col.ind = as.character(sc),
                cex = .1,
                palette = c('grey80','black'),
                col.var = "grey50", # Variables color
                addEllipses = TRUE, # Concentration ellipses
                ellipse.type = "confidence"
)
dev.off()

nn = 25; mdist=0.1; rstr=.25
u = umap(pmat, n_neighbors=nn, min_dist=mdist, verbose=F, repulsion_strength=rstr)
udf <- data.frame(X1=u[,1],X2=u[,2], projid=pmeta$projid)
udf = cbind(udf, emeta[as.character(pmeta$projid),c('niareagansc','cogdx', 'braaksc')])

png(paste0(imgpref, 'pathvars_umap_fullcohort_incohort.png'), res=450, units='in', width=6, height=6)
ptcx = .25
legendloc = 'bottomright'
layout(matrix(1:4, nrow=2))
par(mar=rep(1.1,4))
plot(udf$X1, udf$X2, pch=19, col=colvals[['niareagansc']][udf$niareagansc], 
     axes=F, ylab='',xlab='', cex=ptcx)
mtext('NIA Reagan Score',side=3, line=0)
legend(legendloc,legend=names(colvals[['niareagansc']]), col=colvals[['niareagansc']], 
       pch=19, bty='n')
plot(udf$X1, udf$X2, pch=19, col=colvals[['braaksc']][udf$braaksc], 
     axes=F, ylab='',xlab='', cex=ptcx)
mtext('Braak Stage',side=3, line=0)
legend(legendloc,legend=names(colvals[['braaksc']]), col=colvals[['braaksc']], 
       pch=19, bty='n', ncol=2)
plot(udf$X1, udf$X2, pch=19, col=colvals[['cogdx']][udf$cogdx], 
     axes=F, ylab='',xlab='', cex=ptcx)
mtext('Cognitive Decline',side=3, line=0)
legend(legendloc,legend=names(colvals[['cogdx']]), col=colvals[['cogdx']], 
       pch=19, bty='n',ncol=2)
plot(udf$X1, udf$X2, pch=19, col=ifelse(pmeta$projid %in% pids, 'black',tsp.col('grey85')), 
     axes=F, ylab='',xlab='', cex=ptcx)
mtext('Presence in Cohort',side=3, line=0)
legend(legendloc,legend=c('Yes','No'), col=c('black','grey85'), 
       pch=19, bty='n')
dev.off()


udf = cbind(udf, emeta[as.character(pmeta$projid),c('msex','age_death')])

plot(udf$X1, udf$X2, pch=19, col=ifelse(udf$msex == 1, 'black',tsp.col('grey85')), 
     axes=F, ylab='',xlab='', cex=ptcx)

mtext('Presence in Cohort',side=3, line=0)
legend(legendloc,legend=c('Yes','No'), col=c('black','grey85'), 
       pch=19, bty='n')





