#!/usr/bin/R
# ----------------------------------------
# Collate all GEPs (broad, exc, inh, ast):
# Updated 02/04/2022
# ----------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
options(width=175)

# Directories:
gepdir = paste0(sdbdir, 'GEPs/')
moddir = paste0(sdbdir, 'modules/')
resdir = paste0(moddir, 'resources/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'GEP_comparison_')
cmd = paste('mkdir -p', moddir, resdir)
system(cmd)


# Functions:
# ----------
saveHeatmap = function(ht, pltprefix, w, h, res=400){
    png(paste0(pltprefix, '.png'), res=res, units='in', width=w, height=h)
    print(ht)
    dev.off()
    pdf(paste0(pltprefix, '.pdf'), width=w, height=h)
    print(ht)
    dev.off()
}


# Get list of GEP results to read in:
# -----------------------------------
fns = list.files(pattern='*.csv', path=gepdir)
geplist = list()
for (fn in fns){
    pref = sub("_.*", "", fn)
    df = read.delim(paste0(gepdir, fn), sep=",")
    df = df[,-1]
    ll = lapply(colnames(df), function(x){ c(df[,x]) })
    names(ll) = colnames(df)
    geplist[[pref]] = ll
}

lapply(geplist, length)

# Save GEP results:
saveRDS(geplist, file=paste0(gepdir, 'collated_GEPlists.Rds'))


# Load the core/full modules from resources:
# ------------------------------------------
respref = paste0(resdir, 'modules_resource_')
coremap.list = readRDS(paste0(respref, 'coremap.Rds'))
genemap.list = readRDS(paste0(respref, 'genemap.Rds'))


# Intersect modules / NMF for astrocytes:
# ---------------------------------------
pref = 'Ast'; ccpref = 'Ast'
pref = 'Inh'; ccpref = 'Inh'
pref = 'Exc'; ccpref = 'CTXneurons'
pref = 'Brain'; ccpref = 'All'
ll = geplist[[pref]]
usemap.list = coremap.list

if (ccpref != 'All'){
    usemap = usemap.list[[ccpref]]
    modnums = sort(unique(usemap))
    cc = lapply(modnums, function(x){
        names(usemap)[usemap == x]})
    names(cc) = paste0('M', modnums)
} else {
    cc = list()
    for (name in names(usemap.list)){
        usemap = usemap.list[[name]]
        modnums = sort(unique(usemap))
        subcc = lapply(modnums, function(x){
            names(usemap)[usemap == x]})
        names(subcc) = paste0(name, '-M', modnums)
        cc = c(cc, subcc)
    }
}


NL = length(ll)
NC = length(cc)
imat = matrix(0, nrow=NL, ncol=NC)
for (i in 1:NL){
    gep = ll[[i]]
    for(j in 1:NC){
        mod = cc[[j]]
        imat[i, j] = length(intersect(gep, mod))
    }
}

lct = sapply(ll, length)
cct = sapply(cc, length)
umat = (matrix(rep(cct, NL), byrow=TRUE, nrow=NL, ncol=NC) + 
    matrix(rep(lct, NC), byrow=FALSE, nrow=NL, ncol=NC)) - imat
jmat = imat / umat
colnames(jmat) = names(cc) 
rownames(jmat) = names(ll)
colnames(imat) = names(cc) 
rownames(imat) = names(ll)

# Plot the top modules vs. the GEPs
if (ccpref == 'All'){
    ind = apply(imat, 2, max) > 20
    pltmat = imat[,ind]
} else { pltmat = imat[,1:30] }

pltmat = reord(pltmat)
pltmat = t(reord(t(pltmat)))
diag.list = diag.mat2(t(pltmat), cutoff=10)
pltmat = t(diag.list[[1]])
pltmat = pltmat[,rev(colnames(pltmat))]
ht = Heatmap(pltmat, name='ngenes',
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    use_raster=TRUE, col=colb)

pltprefix = paste0(imgpref, pref, '_ngene_intersection_heatmap')
scale = 0.3
w=scale * ncol(pltmat) + 1; h=scale *nrow(pltmat) + .5
saveHeatmap(ht, pltprefix, w, h)

# Verdict: we perform better on all glia, they perform better on Exc/Inh -> due to how we decorrelate.
# TODO: discuss with Ben why we underperform so much on Exc/Inh


# Get functional enrichment results for each module / each GEP:
# -------------------------------------------------------------
getSingleGost = function(genes, ordered_query=FALSE, term_size=1000){
    sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
    sub.sources = c('KEGG','REAC','WP','CORUM')
    gp2.result = gprofiler2::gost(genes, organism='hsapiens',
        ordered_query=ordq, multi_query=FALSE,
        sources=sources)
    gpdf = gp2.result$result
    if (!is.null(gpdf)){
        gpdf = gpdf[order(gpdf$p_value),]
        gpdf = gpdf[gpdf$term_size < 1000,]
    }
    return(gpdf)
}

lldf = c()
for (i in 1:NL){
    gpdf = getSingleGost(ll[[i]], ordered_query=TRUE, term_size=1000)
    if (!is.null(gpdf)){
        lldf = rbind(lldf, data.frame(p = gpdf$p_value[1],
                term = gpdf$term_name[1], i=names(ll)[i]))
        print(tail(lldf,1))
    }
}
lldf$log10p = -log10(lldf$p)

ccdf = c()
for (i in 1:30){
    gpdf = getSingleGost(cc[[i]], ordered_query=FALSE, term_size=1000)
    if (!is.null(gpdf)){
        ccdf = rbind(ccdf, data.frame(p = gpdf$p_value[1],
                term = gpdf$term_name[1], i=names(cc)[i]))
        print(tail(ccdf,1))
    }
}


# Additionally, run all together:
# -------------------------------
joint.list = c(ll, cc)

gp2.result = gprofiler2::gost(joint.list, organism='hsapiens',
    ordered_query=FALSE, multi_query=TRUE,
    sources=sources)
gpdf = gp2.result$result
gpdf = gpdf[gpdf$term_size < 1000,]


# Plot the p-values:
# ------------------
pmat = matrix(unlist(gpdf$p_value), nrow=nrow(gpdf), byrow=TRUE)
lpmat = -log10(pmat)
lpmat[lpmat > 10] = 10
colnames(lpmat) = names(joint.list)
rownames(lpmat) = gpdf$term_name
lpmat = lpmat[, apply(lpmat, 2, max) > 2]
lpmat = lpmat[apply(lpmat, 1, max) > 2, ]

lpmat = reord(lpmat)
lpmat = t(reord(t(lpmat)))
diag.list = diag.mat2(t(lpmat))
lpmat = t(diag.list[[1]])
rowcut = diag.list[[3]]


# Label only some: check not deleting any rows + add back:
# --------------------------------------------------------
ind = (1:nrow(lpmat) %% 8)
ind = which(ind == 0)
uqr = unique(rowcut)
miss.uqr = uqr[!(uqr %in% unique(rowcut[ind]))]
while (length(miss.uqr) > 0){
    add.uqr = which(rowcut %in% miss.uqr)
    dist.uqr = sapply(add.uqr, function(x){ 
        min(abs(x - ind))})
    ind = c(ind, add.uqr[which.max(dist.uqr)])
    miss.uqr = uqr[!(uqr %in% unique(rowcut[ind]))]
}
rownames(lpmat)[!(1:nrow(lpmat) %in% ind)] = ''


# Make heatmap and plot to files:
# -------------------------------
ht = Heatmap(lpmat, 
    column_split=ifelse(1:ncol(lpmat) %in% grep("^GEP",colnames(lpmat)), 'cNMF','scdemon'),
    row_split = rowcut,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    row_names_gp = gpar(fontsize=8),
    use_raster=TRUE, col=colb)

pltprefix = paste0(imgpref, pref, '_joint_functional_enrichments_heatmap')
w = 10 * ncol(lpmat) / 40
h = 15 * nrow(lpmat) / 1000
saveHeatmap(ht, pltprefix, w, h)





# TODO: Show how we can separate out components of the NMF modules
