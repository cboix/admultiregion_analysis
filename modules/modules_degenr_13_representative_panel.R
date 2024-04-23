#!/usr/bin/R
# ------------------------------------------------------------
# Create a representative figure panel for each set of modules
# Updated 02/09/2022
# ------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))

library(tidyr)
library(viridis)
library(PRROC)

library(ComplexHeatmap)
library(circlize)
options(width=150)

# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
resdir = paste0(sdbdir, 'modules/resources/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'module_panels_')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
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


# Global parameters for heatmaps:
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Colors:
col_fun = colorRamp2(c(-1, 0, 1), c('blue', "white", 'red'))
col_log10p = colorRamp2(c(0, 10), c('white','black'))
colg = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(50)
pcols = brewer.pal(12,'Paired')


# Load in full data:
# ------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))

if (runset == 'All'){
    cls = 'major.celltype'
} else {
    cls = 'cell_type_high_resolution'
}


# # Load in the full pseudobulk data for these subtypes (core genes):
# # NOTE: From modules_degenr_03_plot_pseudobulk_modules.R
# # ------------------------------------------------------
# useset = 'coregenes_'
# scores.file = paste0(moddir, 'module_pseudobulk_scores_', 
#     useset, fullpref, '.tsv.gz')
# scdf = read.delim(gzfile(scores.file), sep="\t")

# TODO: Get the scoremat for only the core genes

rownames(cellmeta) = cellmeta$barcode
submeta = cellmeta[rownames(scoremat),]

# Change scoremat column names:
mmap$mod = paste0('M', mmap$module)
mn.map = mmap$mname
names(mn.map) = mmap$mod
colnames(scoremat) = mn.map[colnames(scoremat)]


# Select top N by genes (> 10 genes):
# -----------------------------------
NGENES = 10
NTOP = 30
ngdf = data.frame(table(coremap))
names(ngdf) = c('module', 'ngenes')
mmap = merge(mmap, ngdf)
mmap = mmap[order(mmap$module),]
mmap = mmap[mmap$ngenes >= NGENES,]
kept.modules = mmap$mod
rownames(mmap) = mmap$mname


# Get hypergeometric enrichment-based covariate scores:
# -----------------------------------------------------
useset = 'coregenes_'
# scores.file = paste0(moddir, 'module_pseudobulk_scores_', fullpref, '.tsv.gz')
# scdf = read.delim(gzfile(scores.file), sep="\t")
# scdf = merge(scdf, metadata[,c('projid','region','rind')])
# TODO: score coregenes only?
hgdf.file = paste0(moddir, 'module_covariate_hgdf_', useset, fullpref, '.tsv')
hgdf = read.delim(hgdf.file, header=T)
hgdf$log10p = -log10(hgdf$p.value)

unique(hgdf$covariate)
hg.covars = c(cls,'region')
hg.covars2 = c('msex','Apoe_e4','nrad','cogdxad','braaksc56')
sub.hgdf = hgdf[hgdf$covariate %in% hg.covars,]
sub.hgdf = merge(sub.hgdf, aggregate(p.value ~ covariate + level + mname, sub.hgdf, min))
sub.hgdf$log10p = (2 * sub.hgdf$cls - 1) * sub.hgdf$log10p

cov.hgdf = hgdf[hgdf$covariate %in% hg.covars2,]
cov.hgdf = cov.hgdf[cov.hgdf$level %in% c('AD','1','yes'),]
cov.hgdf = merge(cov.hgdf, aggregate(p.value ~ covariate + level + mname, cov.hgdf, min))
cov.hgdf$log10p = (2 * cov.hgdf$cls - 1) * cov.hgdf$log10p

core.hgmat = pivot.tomatrix(sub.hgdf[, c('mname','level','log10p')], 'level', 'log10p')
cov.hgmat = pivot.tomatrix(cov.hgdf[, c('mname','covariate','log10p')], 'covariate', 'log10p')
full.hgmat = cbind(core.hgmat, cov.hgmat)


# Load in set of top by p-value functional enrichments for each module:
# On de genes, with term_size < 500 (in this case: "small")
# ---------------------------------------------------------------------
toppvals.file = paste0(moddir, 'module_enrichments_toppvals_',
                       useset, fullpref, '_small.tsv')
pvalsdf = read.delim(toppvals.file, sep="\t")
# Push the very long terms to the back:
pvalsdf$nc = nchar(pvalsdf$term)
pvalsdf$nc.gt = pvalsdf$nc > 40
pvalsdf = pvalsdf[order(pvalsdf$nc.gt),]

# Top two terms per module
topenrdf = lapply(mmap$mname, function(x){
    head(pvalsdf[pvalsdf$mname == x,],1) })
topenrdf = do.call(rbind, topenrdf)
annenrdf = merge(aggregate(p ~ mname, topenrdf, function(x){max(-log10(x))}),
    aggregate(term ~ mname, topenrdf, function(x){ paste(x, collapse='\n') }))
annenrdf = merge(annenrdf, mmap, all.y=TRUE)
annenrdf$p[is.na(annenrdf$p)] = 0
annenrdf$term[is.na(annenrdf$term)] = ''
rownames(annenrdf) = annenrdf$mname


# Get the average score per subtype + correlation matrices:
# ---------------------------------------------------------
subtypes = unique(submeta[[cls]])
tform = make.tform(submeta[[cls]], u=subtypes, norm=TRUE)
reg.tform = make.tform(submeta$region, u=reg.nomb, norm=TRUE)
ct.scoremat = t(scoremat) %*% tform
ct.scoremat = ct.scoremat[mmap$mname,]
reg.scoremat = t(scoremat) %*% reg.tform
reg.scoremat = reg.scoremat[mmap$mname,]

# Average and percent score overall:
pct.cutoff = 0.5
avg.score = t(t(apply(scoremat, 2, mean)))
pct.score = t(t(apply(scoremat > pct.cutoff, 2, mean)))

# Pct score by subtype:
stmat = pct.score
for (st in subtypes){
    ind = (submeta[[cls]] == st)
    stscore = t(t(apply(scoremat[ind,] > pct.cutoff, 2, mean)))
    stmat = cbind(stmat, stscore)
}
# colnames(stmat) = c('All',subtypes)
for (reg in reg.nomb){
    ind = (submeta$region == reg)
    stscore = t(t(apply(scoremat[ind,] > pct.cutoff, 2, mean)))
    stmat = cbind(stmat, stscore)
}
colnames(stmat) = c('All',subtypes, reg.nomb)
regmat = cbind(CTX=apply(stmat[,c('AG','MT','PFC')], 1, mean), stmat[,c('EC','HC','TH')])


# Correlation and order for heatmap plotting:
# -------------------------------------------
cor.mat = cor(scoremat[, mmap$mname])
# NOTE: Find other distance? jaccard scales with sparsity...
# zmat = scoremat[, mmap$mname] >= pct.cutoff
# library(proxy)
# cor.mat = 1 - as.matrix(dist(t(zmat), method='Jaccard'))
# diag(cor.mat) = 0
rmat = reord(cor.mat)
rn = rownames(rmat)

# TODO: Perform regression tests for co-assoc?
# NOTE: All sig cor for cell-level - need to eval at psbulk level 
# cor.test.mat = cor.mat * 0 + 1
# for (i in rn){
#     for (j in rn){
#         ctres = cor.test(scoremat[,i], scoremat[,j])
#         cor.test.mat[i,j] = cr
#     }
# }


# Split rows by identity / region / AD-assoc / pathway / other covariates:
# ------------------------------------------------------------------------
# Annotate programs by covariates and % expression:
kept.subtypes = subtypes[subtypes %in% colnames(full.hgmat)]
attrdf = data.frame(ct = apply(full.hgmat[, kept.subtypes], 1, max),
    reg = apply(full.hgmat[, reg.nomb], 1, max),
    mname = rownames(full.hgmat))
othercov = colnames(full.hgmat)
othercov = othercov[!(othercov %in% c(reg.nomb, kept.subtypes))]
attrdf$other.val = apply(full.hgmat[, othercov], 1, max)
attrdf$other.cov = othercov[apply(full.hgmat[, othercov], 1, which.max)]

# ID core programs (non-subtype, high pct.score)
# TODO: Core programs should be enriched in ~ 100 % of each subtype?
attrdf$pct.all = stmat[attrdf$mname, 1]
attrdf$pct.max = apply(stmat[attrdf$mname,], 1, max)
top.prog = names(which.max(stmat[,'All']))
id.prog = attrdf[attrdf$ct > 3 & attrdf$pct.max > 0.5,'mname']
id.prog = id.prog[id.prog %in% rn]
id.prog = id.prog[id.prog != top.prog]
id.prog = id.prog[grep(" MT-",id.prog, invert=TRUE)]

# Order id programs by region:
regmat = cbind(CTX=apply(full.hgmat[,c('AG','MT','PFC')], 1, mean),
    full.hgmat[,c('EC','HC','TH')])
idmat = regmat[id.prog,] * (regmat[id.prog,] > 0.10)
idmat = sweep(idmat, 1, apply(idmat, 1, sum), '/')
idmat = apply(sweep(idmat, 2, 1:ncol(idmat), '*'), 1, sum)
id.prog = names(sort(idmat))

# Update order:
rn = c(top.prog, id.prog, rn[!(rn %in% c(top.prog, id.prog))])
rn = rn[rn %in% rownames(cor.mat)]
idsplit = ifelse(rn %in% c(top.prog, id.prog), 'Identity','Other') 

# Order subtypes by id program order:
idmat = t(full.hgmat[id.prog, kept.subtypes])
idmat[idmat > 1000] = 1000
idmat[idmat < -1000] = -1000
idmat = idmat * (idmat > -10)
idmat = sweep(idmat, 1, apply(idmat, 1, sum), '/')
idmat = apply(sweep(idmat, 2, 1:ncol(idmat), '*'), 1, sum)
subtype.ord = names(sort(idmat))
reg.ord = c('AG','MT','PFC','EC','HC','TH')


# Reorder for heatmap plotting:
# -------------------------------
pltcor = cor.mat[rn, rn]
pltmat = cbind(ct.scoremat[rn,], reg.scoremat[rn,])
pltmat = sweep(pltmat, 1, apply(pltmat, 1, max), '/')
pltmat.split = ifelse(colnames(pltmat) %in% reg.nomb, 'Region','Subtype')
rownames(pltcor) = sub(".*\\: ","", rownames(pltcor))
colnames(pltcor) = sub(".*\\: ","", colnames(pltcor))

# Make the annotation-style heatmaps:
# -----------------------------------
ux = 1.75
# Text:
hg1 = rowAnnotation(top.genes = anno_text(sub(".*\\: ","",rn), gp=gpar(fontsize=5)))
hn1 = rowAnnotation(top.genes = anno_text(sub(" \\(.*","",rn), gp=gpar(fontsize=5)))

# Number of DE up/down:
path = 'cogdxad'
subdedf = dedf[dedf$dkey == path,]
cdf = data.frame(table(subdedf[,c('module','gset')]))
cdf = spread(cdf,gset, Freq)
cdf = merge(cdf, mmap, all.y=TRUE)
# m = as.matrix(cdf[,c('Up','--','Down')])
m = as.matrix(cdf[,c('Up','Down')])
rownames(m) = cdf$mname
m = t(apply(m, 1, function(x){ x / sum(x)}))
m[is.na(m)] = 0
hde = rowAnnotation(cog.DE = anno_barplot(m[rn,], 
        gp = gpar(fill = c(pcols[6],pcols[2]), col=NA),
        numbers_gp = gpar(fontsize=5), 
        border=FALSE,
        bar_width = 1, width = unit(ux * 2, "mm")), 
    gp=gpar(fontsize=5), annotation_name_gp=gpar(fontsize=5.5))


# Number of genes and pct of cells:
# ---------------------------------
ngmat = as.matrix(mmap[rn, 'ngenes', drop=F])
htng = Heatmap(ngmat, 
    cluster_rows=FALSE,
    name='ngenes',
    width = ncol(ngmat)*unit(ux * 2, "mm"), 
    height = nrow(ngmat)*unit(ux, "mm"),
    row_split=idsplit,
    border_gp=gpar(color='black', lwd=.5),
    col=colg,
    cell_fun = function(j, i, x, y, w, h, col){
        grid.text(ngmat[i,j], x, y,
            gp=gpar(col=ifelse(ngmat[i,j] > .60 * max(ngmat),'white','black'), 
                fontsize=gridtxt.fs))}
)


pctmat = stmat[rn, c('All', subtype.ord) ,drop=F]
htpct = Heatmap(pctmat, 
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    row_split=idsplit,
    column_split=ifelse(colnames(pctmat) == 'All','', 'Subtype'),
    name='pct.score',
    width = ncol(pctmat)*unit(ux * 2, "mm"), 
    height = nrow(pctmat)*unit(ux, "mm"),
    border_gp=gpar(color='black', lwd=.5),
    col=colg,
    cell_fun = function(j, i, x, y, w, h, col){
        pct = 100 * pctmat[i,j]
        font = ifelse(pctmat[i,j] >= max(pctmat[i,]), 2, 1)
        val = sub("^[ ]*","", formatC(round(pct, 1), digits=2))
        val[val == '1e+02'] = '100'
        grid.text(val, x, y, just='center',
            gp=gpar(col=ifelse(pct > 60 ,'white','black'),
                font=font, fontsize=gridtxt.fs))}
)

# Create matrix for plotting enrichment:
# TODO: fill with NA if missing vals
annenrdf = annenrdf[rn,]
enrmat = as.matrix(annenrdf[,'p',drop=F])
colnames(enrmat) = '-log10p'
rownames(enrmat) = annenrdf$term
th.enrmat = enrmat  # Thresholded, for colors:
th.enrmat[th.enrmat > 10] = 10
htenr = Heatmap(th.enrmat, 
    cluster_rows=FALSE,
    name='-log10p',
    row_split=idsplit,
    width = ncol(enrmat)*unit(ux * 2, "mm"), 
    height = nrow(enrmat)*unit(ux, "mm"),
    col=colg,
    cell_fun = function(j, i, x, y, w, h, col){
        enr = enrmat[i,j]
        val = sub("^[ ]*","", formatC(round(enr,1), digits=2))
        if (enr != 0){
            grid.text(val, x, y, 
                gp=gpar(col=ifelse(enr < 5, 'black','white'), 
                    fontsize=gridtxt.fs))
        }
    }
)


# Make the covariate enrichment heatmaps:
# ---------------------------------------
full.plthg = full.hgmat[rn,]
full.plthg[full.plthg > 10] = 10
full.plthg[full.plthg < -10] = -10
full.plthg[abs(full.plthg) < -log10(0.05)] = 0

# Order:
full.plthg = t(reord(t(full.plthg)))
out = diag.mat2(full.plthg, cutoff=3)[2][[1]]
out = c(subtype.ord, reg.ord, out[!(out %in% c(reg.nomb, kept.subtypes))])
full.plthg = full.plthg[,out]
full.hg.split = ifelse(colnames(full.plthg) %in% reg.nomb, 'Region',
    ifelse(colnames(full.plthg) %in% kept.subtypes, 'Subtype', 'Covariate'))
full.hg.split = factor(full.hg.split, levels=c('Subtype','Region','Covariate'))

htcovar = Heatmap(full.plthg, 
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    column_split=full.hg.split,
    cluster_column_slices = FALSE,
    name='-log10\np-value',
    row_split=idsplit,
    border_gp=gpar(color='black', lwd=.5),
    width = ncol(full.plthg)*unit(ux, "mm"), 
    height = nrow(full.plthg)*unit(ux, "mm"),
    col=rev(colrb),
    cell_fun = function(j, i, x, y, w, h, col){
        # Add the p-value text (placeholder for now)
        cr = abs(full.plthg[i,j])
        if (cr >= 3){
            grid.text('*', x, y, 
                gp=gpar(col=ifelse(cr > 6, 'white','black'), 
                    fontsize=gridtxt.fs))}
    }
)

# Make the correlation heatmaps:
# ------------------------------
ht2 = Heatmap(pltcor,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    row_split=idsplit,
    column_split=idsplit,
    name='corr',
    col=col_fun,
    width = ncol(pltcor)*unit(ux, "mm"), 
    height = nrow(pltcor)*unit(ux, "mm"),
    border_gp=gpar(color='black', lwd=.5),
    cell_fun = function(j, i, x, y, w, h, col){
        # Add the p-value text (placeholder for now)
        cr = pltcor[i,j]
        txtcr = ifelse(cr > 0.4, '*','')
        grid.text(txtcr, x, y, gp=gpar(fontsize=gridtxt.fs))}
)


# Join heatmaps + save plots:
# ---------------------------
# ht = hn1 + htng + hg1 + htpct + htcovar + hde + ht2 + htenr
ht = hn1 + htng + hg1 + htpct + htcovar + ht2 + htenr

h = 2 + 1 / 15 * nrow(pltmat)
w = 2 + 1 / 10 * ((ncol(ngmat) + ncol(pctmat) + ncol(enrmat) + 1) * 2 + ncol(full.plthg) + ncol(cor.mat))
w = ifelse(w > 25, 25, w)
pltprefix = paste0(imgpref, 'score_cor_prelim_', fullpref)
saveHeatmap(ht, pltprefix, w=w, h=h)

