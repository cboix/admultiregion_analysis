#!/usr/bin/R
# -----------------------------------------------------
# Starting from all DEG results,
# - Plot the lFC for diff. pathology against each other
# NOTE: mostly preliminary, not finished analysis
# Updated: 05/18/23
# -----------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(ggrastr)
print(version)
options(width=170)

# Directories:
srdir = paste0(sdbdir, 'subtype_reg/')
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
cmd = paste('mkdir -p', plotdir, regdir, enrdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))
source(paste0(sbindir, 'auxiliary_goterm_functions.R'))
source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))


# Function for simple heatmap:
# ----------------------------
simpleHeatmap = function(pltmat, ux=1.5, col=viridis(50), cluster_columns=FALSE){
    ht = Heatmap(pltmat, 
        use_raster=FALSE, 
        cluster_columns=cluster_columns,
        col=col,
        border_gp=gpar(color='black', lwd=.5),
        width=ncol(pltmat) * unit(ux, 'mm'),
        height=nrow(pltmat) * unit(ux, 'mm'))
    return(ht)
}


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
setnames = c("Mic_Immune_Mic"='Microglial', "Ast_Ast"='Astrocyte', 
    "Opc_Opc"='OPC', "Oli_Oli"='Oligodendrocyte', 
    'Inh_Inh'='Inhibitory neuron','Exc_Exc'='Excitatory neuron')

pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
pcols = brewer.pal(12, 'Paired')
sigmap = c('0'= 'grey85',
    '22' = pcols[6], '11' = pcols[2],
    '21' = pcols[10], '12' = pcols[9],
    '20' = pcols[5], '2' = pcols[7],
    '10' = pcols[1], '1' = pcols[3])


# Load all DEG sets:
# -----------------
kept.cols = c('gene','col_nm','path','region', 'logFC_nb', 'p_nb')
alldf = c()
for (path in pathlist){
    mstr = paste0('allmethods.regional_', path)
    fullaggrda = paste0(regdir, mstr, '.merged.rda')
    load(fullaggrda)

    for (set in keep.sets){
        print(set)
        setdf = setdflist[[set]][, kept.cols]
        setdf$set = set
        alldf = rbind(alldf, setdf)
    }
}


# Raw numbers of DEGs:
# --------------------
sigdf = alldf[alldf$col_nm != 0,]
aggdf = agg.rename(gene ~ col_nm + path + set + region, sigdf, length, 'nde')

nmat = pivot.tomatrix(aggdf[(aggdf$path == 'nft') & (aggdf$col_nm == 2),
    c('set', 'region', 'nde')], 'region', 'nde')
pmat = pivot.tomatrix(aggdf[(aggdf$path == 'plaq_n') & (aggdf$col_nm == 2),
    c('set', 'region', 'nde')], 'region', 'nde')

mx = max(c(nmat,pmat)) * 1.01
nx = min(c(nmat,pmat)) * .99
ovl.col_fun = colorRamp2(seq(nx, mx, length.out=100), col2)

ht.n = simpleHeatmap(nmat, col=ovl.col_fun)
ht.p = simpleHeatmap(pmat, col=ovl.col_fun)
ht = ht.p + ht.n

pltprefix = paste0(imgpref, 'regional_nftplaq_ndeg')
w = 2 + c(ncol(nmat) + ncol(pmat)) / 15
h = 2 + nrow(nmat) / 15 
saveHeatmap(ht, pltprefix, w=w, h=h)



# DEG concordance versus NIA-Reagan:
# ----------------------------------
reg.noth = c('allregions', reg.nomb[reg.nomb != 'TH'])
ovldf = c()
for (set in keep.sets){
    setdf = alldf[alldf$set == set,]
    for (region in reg.noth){
        cat(set, '\t', region, '\n')
        fulldf = setdf[setdf$region == region,]
        cmat = pivot.tomatrix(fulldf[,c('gene','path','col_nm')], 'path', 'col_nm')
        c1 = cmat[cmat[,'nrad'] == 1,]
        c2 = cmat[cmat[,'nrad'] == 2,]
        dwdf = data.frame(path=colnames(c1), novl=colSums(c1 == 1), ng=nrow(c1), dir='Down')
        updf = data.frame(path=colnames(c2), novl=colSums(c2 == 2), ng=nrow(c2), dir='Up')
        df = rbind(updf, dwdf)
        df$frac = df$novl / df$ng
        df$set = set
        df$region = region
        ovldf = rbind(ovldf, df)
    }
}
rownames(ovldf) = NULL

aggdf = aggregate(cbind(novl, ng) ~ path + set + region, ovldf, sum)
aggdf$frac = aggdf$novl / aggdf$ng

nmat = pivot.tomatrix(aggdf[aggdf$path == 'nft',
    c('set', 'region', 'frac')], 'region', 'frac')
pmat = pivot.tomatrix(aggdf[aggdf$path == 'plaq_n',
    c('set', 'region', 'frac')], 'region', 'frac')

mx = max(c(nmat,pmat)) * 1.01
nx = min(c(nmat,pmat)) * .99
ovl.col_fun = colorRamp2(seq(nx, mx, length.out=100), col1)

ht.n = simpleHeatmap(nmat, col=ovl.col_fun)
ht.p = simpleHeatmap(pmat, col=ovl.col_fun)
ht = ht.p + ht.n

pltprefix = paste0(imgpref, 'regional_nftplaq_nradOvl')
w = 2 + c(ncol(nmat) + ncol(pmat)) / 15
h = 2 + nrow(nmat) / 15 
saveHeatmap(ht, pltprefix, w=w, h=h)



# Get the R^2 values for all cell types + regions:
# ------------------------------------------------
reg.noth = c('allregions', reg.nomb[reg.nomb != 'TH'])
r2mat = matrix(0, nc=length(reg.noth), nr=length(keep.sets),
    dimnames=list(keep.sets, reg.noth))
for (set in keep.sets){
    setdf = alldf[alldf$set == set,]
    for (region in reg.noth){
        cat(set, '\t', region, '\n')
        fulldf = setdf[setdf$region == region,]
        ldf = spread(fulldf[,c('gene','path','logFC_nb')], path, logFC_nb)
        fit = lm(nft ~ plaq_n, ldf)
        r2mat[set, region] = summary(fit)$adj.r.sq
    }
}

ht = simpleHeatmap(r2mat, ux=1.5, col=viridis(50))
pltprefix = paste0(imgpref, 'regional_nftplaq_r2heatmap')
w = 2 + ncol(r2mat) / 15
h = 2 + nrow(r2mat) / 15 
saveHeatmap(ht, pltprefix, w=w, h=h)


# Score genes w.r.t. involvement in either pathology:
# ---------------------------------------------------
resdf = c()
for (set in keep.sets){
    print(set)
    setdf = alldf[alldf$set == set,]
    setdf = setdf[setdf$region != 'TH',]
    ldf = spread(setdf[,c('gene','path', 'region','logFC_nb')], path, logFC_nb)
    cdf = spread(setdf[,c('gene','path', 'region','col_nm')], path, col_nm)
    names(cdf)[3:ncol(cdf)] = paste0('sig.',names(cdf)[3:ncol(cdf)])
    ldf = merge(ldf, cdf)
    fit = lm(nft ~ plaq_n * region, ldf)
    ldf$nft.pred = predict(fit, ldf)
    ldf$nft.resid = ldf$nft - ldf$nft.pred

    # Average over all tested sets:
    ldf$ntest = 1
    ldf$sig.nft.up = ldf$sig.nft == 2
    ldf$sig.plaq.up = ldf$sig.plaq_n == 2
    ldf$sig.nft.dw = ldf$sig.nft == 1
    ldf$sig.plaq.dw = ldf$sig.plaq_n == 1
    aggdf = aggregate(cbind(sig.nft.up, sig.plaq.up, sig.plaq.dw, sig.nft.dw, ntest) ~ gene, ldf, sum)
    meandf = aggregate(cbind(nft.resid, nft, plaq_n)  ~ gene, ldf, mean)
    aggdf = merge(aggdf, meandf)

    # Plaque associated:
    aggdf = aggdf[order(aggdf$nft.resid),]
    up.ind = (aggdf$sig.plaq.up >= 3) & (aggdf$sig.nft.up <= 2)
    dw.ind = (aggdf$sig.plaq.dw >= 3) & (aggdf$sig.nft.dw <= 2)
    agg.plaq = aggdf[up.ind | dw.ind,]
    agg.plaq$dir = 'Plaque'
    print(head(agg.plaq, 10))

    aggdf = aggdf[order(aggdf$nft.resid, decreasing=T),]
    up.ind = (aggdf$sig.nft.up >= 3) & (aggdf$sig.plaq.up <= 2)
    dw.ind = (aggdf$sig.nft.dw >= 3) & (aggdf$sig.plaq.dw <= 2)
    # agg.nft = aggdf[aggdf$sig.nft.up >= 3,]  
    # agg.nft = agg.nft[agg.nft$sig.plaq.up <= 2,]  
    agg.nft = aggdf[up.ind | dw.ind,]  
    agg.nft$dir = 'NFT'
    print(head(agg.nft, 10))

    # Add to table:
    topdf = rbind(agg.nft, agg.plaq)
    topdf$set = set
    resdf = rbind(resdf, topdf)
}

# Save table (explore in 2h):
res.file = paste0(regdir, 'allmethods.major.plaq_nft.extremegenes.tsv')
write.table(resdf, res.file, quote=F, row.names=F, sep="\t")

#     # Genes higher in NFT (must be up-regulated in NFT):
#     nft.df = ldf[ldf$sig.nft == 2,]
#     nft.df = nft.df[order(nft.df$nft.resid, decreasing=T),]
#     head(nft.df[nft.df$region == 'EC',])
#     # EC: SNX31, HSPA1A, 
#     head(nft.df[nft.df$region == 'HC',])

#     # Genes higher in plaque (must be up-regulated in plaq_n):
#     plaq.df = ldf[ldf$sig.plaq_n == 2,]
#     plaq.df = plaq.df[order(plaq.df$nft.resid, decreasing=F),]
#     head(plaq.df[plaq.df$region == 'EC',])



# Plot the DEGs for a pair + a region:
# ------------------------------------
# Arguments
region = 'HC'
set = 'Mic_Immune_Mic'
pair = c('nft','plaq_n')
for (region in c('allregions', reg.nomb)){
    suff = paste0(set, '.', region, '.', pair[1], '-', pair[2])

    # Build dataframe:
    fulldf = alldf[(alldf$set == set) & (alldf$region == region),]
    ldf = spread(fulldf[,c('gene','path','logFC_nb')], path, logFC_nb)
    pdf = spread(fulldf[,c('gene','path','col_nm')], path, col_nm)
    names(ldf)[-1] = paste0('logFC_', names(ldf)[-1])
    names(pdf)[-1] = paste0('col_', names(pdf)[-1])
    df = merge(pdf, ldf)

    # Signficance + plot:
    lfc = paste0('logFC_', pair)
    sig = paste0('col_', pair)
    df$col = df[[sig[1]]] + 10 * df[[sig[2]]]
    df = df[order(df$col),]
    df$c1 = ifelse(df[[sig[1]]] == 2, paste0(pair[1], ':up'), ifelse(df[[sig[1]]] == 1, paste0(pair[1],':down'), ''))
    df$c2 = ifelse(df[[sig[2]]] == 2, paste0(pair[2], ':up'), ifelse(df[[sig[2]]] == 1, paste0(pair[2],':down'), ''))
    df$sig = paste(df$c1, df$c2)
    labdf = df[df$sig != " ",]

    # Make the set of colors:
    cdf = unique(df[,c('sig','col')])
    cdf$color = sigmap[as.character(cdf$col)]
    sigcols = cdf$color
    names(sigcols) = cdf$sig

    gp = ggplot(df, aes_string(lfc[1], lfc[2], color='sig')) + 
        rasterize(geom_point(cex=.25), dpi=450) +
        scale_color_manual(values=sigcols, name='DE status:') +
        geom_hline(yintercept=0, lty='dotted') + 
        geom_vline(xintercept=0, lty='dotted') + 
        geom_abline(intercept=0, slope=1, lty='dotted', color='grey50') + 
        geom_text_repel(data=labdf, aes_string(lfc[1], lfc[2], label='gene', color='sig'), 
            box.padding=0.075, cex=3, max.overlaps=20, min.segment.length=0.1) + 
        labs(x=paste0('logFC (', pair[1], ')'), y=paste0('logFC (', pair[2], ')'), title=paste0(setnames[set], ' DEGs (', region, ')')) + 
        theme_pubr() + theme(legend.position='right')
    pltprefix = paste0(imgpref, 'regional_comp.', suff)
    saveGGplot(gp, pltprefix, w=8, h=6)
}


# Show how response to nft + plaque agrees across regions:
# --------------------------------------------------------
for (set in keep.sets){
    print(set)
    suff = paste0(set, '.', pair[1], '-', pair[2])

    # Build dataframe:
    fulldf = alldf[(alldf$set == set) & (alldf$region != 'TH'),]
    ldf = spread(fulldf[,c('gene','region', 'path','logFC_nb')], path, logFC_nb)
    pdf = spread(fulldf[,c('gene','region', 'path','col_nm')], path, col_nm)
    names(ldf)[-c(1:2)] = paste0('logFC_', names(ldf)[-c(1:2)])
    names(pdf)[-c(1:2)] = paste0('col_', names(pdf)[-c(1:2)])
    df = merge(pdf, ldf)

    # Signficance + plot:
    lfc = paste0('logFC_', pair)
    sig = paste0('col_', pair)
    df$col = df[[sig[1]]] + 10 * df[[sig[2]]]
    df = df[order(df$col),]
    df$c1 = ifelse(df[[sig[1]]] == 2, paste0(pair[1], ':up'), ifelse(df[[sig[1]]] == 1, paste0(pair[1],':down'), ''))
    df$c2 = ifelse(df[[sig[2]]] == 2, paste0(pair[2], ':up'), ifelse(df[[sig[2]]] == 1, paste0(pair[2],':down'), ''))
    df$sig = paste(df$c1, df$c2)

    ctdf = agg.rename(gene ~ sig + region, df[df$region != 'TH',], length, 'nDEG')
    tmat = pivot.tomatrix(ctdf, 'region', 'nDEG')
    tmat[is.na(tmat)] = 0
    ctdf = ctdf[ctdf$sig != ' ',]

    # Make the set of colors:
    cdf = unique(df[,c('sig','col')])
    cdf$color = sigmap[as.character(cdf$col)]
    sigcols = cdf$color
    names(sigcols) = cdf$sig
    # order:
    cdf$col = factor(cdf$col, levels=c(2, 20, 0, 1, 10, 21, 12, 22, 11))
    cdf = cdf[order(cdf$col),]
    siglvls = cdf$sig
    ctdf$sig = factor(ctdf$sig, levels=siglvls)

    gp = ggplot(ctdf, aes(region, nDEG, fill=sig)) + 
        geom_bar(stat='identity', position='fill') + 
        scale_y_continuous(expand=c(0,0)) + 
        scale_fill_manual(values=sigcols, name='DE status:') +
        labs(x='Region', y='# DEGs', title=paste0(setnames[set], ' DEG agreement across regions (', paste(pair, collapse=' + '), ')')) + 
        theme_pubr() + theme(legend.position='right')
    pltprefix = paste0(imgpref, 'regional_degbarplot.', suff)
    saveGGplot(gp, pltprefix, w=4.5, h=4)
}


# Look at top degs for AG - Ast, Oli:
# -----------------------------------
region = 'AG'
set = 'Opc_Opc'
path = 'plaq_n'
sigdf$signif = ifelse(sigdf$col_nm == 2, 'Up', 'Down')
setdf = sigdf[(sigdf$set == set) & (sigdf$path == path),]
ndf = agg.rename(region ~ gene + col_nm, setdf, length, 'nsig')
setdf = merge(setdf, ndf)

# Genes unique to AG:
subdf = setdf[(setdf$region == region),]
subdf = subdf[order(subdf$p_nb),]

uqdf = subdf[(subdf$nsig <= 1) & (subdf$col_nm == 2),]
outfile = paste0(regdir, 'unique_degs.', region, '.', set, '.', path, '.tsv')
write.table(uqdf, outfile, quote=F, row.names=F, sep="\t")

cat(head(subdf[(subdf$nsig <= 1) & (subdf$col_nm == 2),'gene'], 100))
head(subdf[(subdf$nsig <= 1) & (subdf$col_nm == 2),], 40)

cat(head(subdf[(subdf$nsig >= 5) & (subdf$col_nm == 2),'gene'], 200))
head(subdf[(subdf$nsig >= 5) & (subdf$col_nm == 2),], 40)

