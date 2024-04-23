#!/usr/bin/R
# --------------------------------------------------
# Analysis of metabolic modules
# with each other / across cell types
# both at sample & cell type level
# Updated 02/07/2022
# --------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))
source(paste0(sbindir, 'modules/auxiliary_modules_psbulk.R'))
source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


library(tidyr)
library(viridis)
library(PRROC)
library(gprofiler2)

library(ggpubr)
library(ComplexHeatmap)
library(circlize)
options(width=170)

# Directories:
moddir = paste0(sdbdir, 'modules/')
srdir = paste0(sdbdir, 'subtype_reg/')
resdir = paste0(sdbdir, 'modules/resources/')
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'module_metab_')
cmd = paste('mkdir -p', plotdir, moddir, resdir)
system(cmd)



# Heatmap plotting functions:
# ---------------------------
plotDEgenesHeatmap = function(cmat, pmat, ux, col.split=NULL, row.split=NULL, cluster=TRUE, topann=NULL, raster=FALSE, hscale=1){
    plt = Heatmap(cmat,
        col=col_fun,
        use_raster=raster,
        column_split=col.split,
        row_split=row.split,
        cluster_columns=cluster,
        cluster_rows=cluster,
        cluster_row_slices=cluster,
        top_annotation=topann,
        width = ncol(cmat)*unit(ux, "mm"), 
        height = nrow(cmat)*unit(ux * hscale, "mm"),
        border_gp = gpar(col="black", lty = 1, lwd=.5),
        cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
            p = pmat[i,j]
            if (p < 0.05){
                grid.text('*', x, y, vjust=.75, gp=gpar(fontsize=gridtxt.fs*1.1)) 
            }
        }
    )
    return(plt)
}


# Set the run arguments:
# ----------------------
runset = 'Mic_Immune'
graph_id = 'boot'

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


getMax = function(map, geneset){
    geneset = geneset[geneset %in% names(map)]
    x = map[geneset]
    x = sort(table(x), decreasing=T)
    names(x)[1]
}

gly.genes = c('PDK1','PFKL','PFKP','LDHA','VEGFA','DDIT4', 'ENO1', 'PGK1','BNIP3L')
mt.genes = c('MT-ND3','MT-CO3','MT-CYB','MT-ND4', 'MT-ND2','MT-ND1','MT-ATP8', 'MT-ND4L')
ins.genes = c('HIF3A','FKBP5','FOXG1')


fulldf = c()
glylist = list()
extlist = list()
for (runset in c('Ast','Mic_Immune','Opc')){
    print(runset)
    # NOTE: NEED TO UPDATE OPC RES IN OUTPOST
    # Load in data (scoremat):
    commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
    source(paste0(sbindir, 'modules/load_modules_degenr.R'))

    # Identify metabolism modules:
    gly.mod = getMax(coremap, gly.genes)
    mt.mod = getMax(coremap, mt.genes)
    if (runset == 'Opc'){ ins.mod = getMax(coremap, ins.genes) }

    glylist[[runset]] = names(coremap[coremap == gly.mod])
    extlist[[runset]] = names(genemap[genemap == gly.mod])

    # Compare the modules correlation by subtype:
    rownames(cellmeta) = cellmeta$barcode
    submeta = cellmeta[rownames(scoremat),]
    gly.x = scoremat[,paste0('M', gly.mod)]
    mt.x = scoremat[,paste0('M', mt.mod)]
    if (runset == 'Opc'){ ins.x = scoremat[,paste0('M', ins.mod)]}
    qcutoff = 0.9
    gly.t = quantile(gly.x, qcutoff)
    mt.t = quantile(mt.x, qcutoff)
    gly.bin = 1 * (gly.x >= gly.t)
    mt.bin = 1 * (mt.x >= mt.t)
    cr = cor(gly.x, mt.x)
    cat(round(cr, 3),'\t', runset, '\n')

    sts = unique(submeta$cell_type_high_resolution)
    for (st in sts){
        # By correlation:
        ind = (submeta$cell_type_high_resolution == st)
        cr = cor(gly.x[ind], mt.x[ind])
        cat(round(cr, 3),'\t', st,'\n')
        # By overlap:
        mat = table(gly.bin[ind], mt.bin[ind])
        ft = fisher.test(mat)
        cat(round(ft$estimate, 3),'\t', ft$p.value,'\n')
    }

    # Save at individual-level:
    keep.cols = c('barcode','region','cell_type_high_resolution','projid')
    if (runset == 'Opc'){
        df = data.frame(gly=gly.x, mt=mt.x, ins=ins.x)
        df = cbind(df, submeta[,keep.cols])
        aggdf = aggregate(cbind(gly, mt, ins) ~ region + projid, df, mean)
    } else {
        df = data.frame(gly=gly.x, mt=mt.x)
        df = cbind(df, submeta[,keep.cols])
        if (runset == 'Mic_Immune'){
            df = df[df$cell_type_high_resolution != 'T cells',]
        }
        aggdf = aggregate(cbind(gly, mt) ~ region + projid, df, mean)
    }
    cind = 3:ncol(aggdf)
    names(aggdf)[cind] = paste0(runset, '-', names(aggdf)[cind])
    if (is.null(fulldf)){
        fulldf = aggdf
    } else {
        fulldf = merge(fulldf, aggdf, all=TRUE)
    }
}


# Plot the correlation between all of these different glial modules
# -----------------------------------------------------------------
modmat = as.matrix(fulldf[,3:(ncol(fulldf)-1)])
cr = cor(modmat)
col_fun = colorRamp2(c(-1, 0, 1), c(colrb[100], "white", colrb[1]))

ux = 1.5
ht = Heatmap(cr, 
    name='corr.',
    col=col_fun,
    border_gp=gpar(color='black', lwd=.5),
    use_raster=TRUE,
    width = ncol(cr)*unit(ux, "mm"), 
    height = nrow(cr)*unit(ux, "mm"),
    cell_fun = function(j, i, x, y, w, h, col){ # Add the p-value text
        txtcr = round(cr[i,j], 2) # TODO: Format better?
        txtcr = formatC(txtcr, digits=1)
        grid.text(txtcr, x, y, gp=gpar(fontsize=4))}
)

h = 5 + 1 / 15 * nrow(cr)
w = 5 + 1 / 15 * ncol(cr)
pltprefix = paste0(imgpref, 'crossct_corrheatmap_', fullpref)
saveHeatmap(ht, pltprefix, w=w, h=h)


# Compute and plot enrichments not shared:
# ----------------------------------------
uselist = extlist
core.ext = intersect(intersect(uselist$Ast, uselist$Opc), uselist$Mic_Immune)

uqlist = list()
uqlist[['Core']] = core.ext
for (nam in names(uselist)){
    nam.genes = uselist[[nam]]
    uqlist[[nam]] = nam.genes[!(nam.genes %in% core.ext)]
}

sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
gp2.result = gprofiler2::gost(uqlist, organism='hsapiens',
                              ordered_query=FALSE, multi_query=TRUE,
                              sources = sources)
gpdf = gp2.result$result


# Load DE results and add to this:
# --------------------------------
geneset = 'metabolic_genes'
full.file = paste0(regdir, 'allmethods.allmajor.', geneset, '.merged.rda')
load(full.file)
regdedf = c()
for (set in names(setdflist)){
    setdf = setdflist[[set]]
    if (!is.null(setdf)){
        setdf$set = set
        # setdf$region = region
        regdedf = rbind(regdedf, setdf)
    }
}
fulldedf = regdedf[regdedf$region == 'allregions',]


# gset = c('HK1','HK2','HK3','GPI','PFKP','PFKL','PFKM','PGK1','PGK2','PGAM1','PGAM2','GAPDH','ENO1','ENO2','PKM','LDHB')
# fulldedf[(fulldedf$gene %in% gset) & (fulldedf$set == 'Ast_Ast') & (fulldedf$path == 'plaq_d'),]


# Functions for DE processing:
# ----------------------------
getDEmatrices = function(dedf, pathlist, setlist, genes){
     # Update scale for binary vars:
    pind = dedf$path %in% c('cogdxad','nrad', 'braaksc.early','braaksc.ad')
    dedf$logFC_nb[pind] = dedf$logFC_nb[pind] / 25 # Because scale at 0.04 is convertible
    # Subset:
    dedf = dedf[dedf$set %in% setlist,]
    dedf = dedf[dedf$path %in% pathlist,]
    # Get matrices
    dedf$lp = ifelse(dedf$col_nm > 0, dedf$log10p_nm, 0)
    cmat = pivot.tomatrix(dedf[, c('set2','gene','logFC_nb')], 'gene','logFC_nb')
    pmat = pivot.tomatrix(dedf[, c('set2','gene','lp')], 'gene','lp')
    pmat = 10**(-pmat)
    cmat[is.na(cmat)] = 0
    pmat[is.na(pmat)] = 1
    full.cmat = matrix(0, nrow=nrow(cmat), ncol=length(genes),
                       dimnames=list(rownames(cmat), genes))
    full.pmat = matrix(1, nrow=nrow(cmat), ncol=length(genes),
                       dimnames=list(rownames(cmat), genes))
    full.cmat[,colnames(cmat)] = cmat
    full.pmat[,colnames(pmat)] = pmat
    # Row order:
    sets = sort(unique(dedf$set))
    rdf = expand.grid(path=pathlist, set=sets)
    rn = paste0(rdf$set, '@', rev(rdf$path))
    row.split = sub("@.*","", rownames(full.cmat))
    return(list(cmat=full.cmat[rn,], pmat=full.pmat[rn,],
                row.split=row.split, rn=rn))
}


# Plot in/out of module for these genes:
# --------------------------------------
setlist = c('Ast_Ast','Mic_Immune_Mic','Opc_Opc')
pathlist = c('cogdxad','nrad','nft','plaq_n','plaq_d')

# List of genes that shows DE effects in any celltype 
degenes = fulldedf$gene[(fulldedf$col_nm != 0) & (fulldedf$set %in% setlist) & (fulldedf$path %in% pathlist)]
# core.genes = c('PDK1','PFKL','PFKP','PGK1','TPI1','LDHA','VEGFA','DDIT4', 'BNIP3L')
core.genes = uqlist[['Core']]
fa.genes = c('ANGPTL4', 'HILPDA', 'IRS2','PRDM16','ADCY8')
glycogen.genes = c('GSK3B','PYGL','PPP1R3E', 'AGL', 'GBE1', 'UGP2', 'GYS1')
other.genes = c('GAPDH','HK1', 'HK2', 'SLC2A1','GFPT2', 'GLUL','SLC38A1','SLC38A2')
core.genes = core.genes[!(core.genes %in% c(fa.genes,glycogen.genes,other.genes))]


fdf = data.frame(rbind(cbind(core.genes, 'Core'),
                       cbind(fa.genes, 'FA'),
                       cbind(glycogen.genes, 'Glycogen'),
                       cbind(other.genes, 'Other')))
names(fdf) = c('gene','set')

fdf = fdf[fdf$gene %in% degenes,]
full.genes = fdf$gene
col.split = fdf$set
rownames(fdf) = fdf$gene
fmat = sapply(extlist, function(x){ full.genes %in% x })
rownames(fmat) = full.genes
fmat = t(fmat)

ux = 1.5
ht = Heatmap(fmat * 1, 
    name='inMod',
    col=c('white','royalblue'),
    column_split=col.split,
    border_gp=gpar(color='black', lwd=.5),
    use_raster=FALSE,
    width = ncol(fmat)*unit(ux, "mm"), 
    height = nrow(fmat)*unit(ux, "mm"),
)

h = 2 + 1 / 15 * nrow(fmat)
w = 2 + 1 / 15 * ncol(fmat)
pltprefix = paste0(imgpref, 'crossct_inmodheatmap')
saveHeatmap(ht, pltprefix, w=w, h=h)


# Add the DE results:
# -------------------
dedf = fulldedf[(fulldedf$gene %in% full.genes),]
dedf$set2 = paste0(dedf$set, '@', dedf$path)
ll = getDEmatrices(dedf, pathlist, setlist, genes=full.genes)

rmat = reord(t(ll$cmat))
cn = rownames(rmat)
col.split = fdf[cn, 'set']

ux = 1.5
ht = Heatmap(fmat[,cn] * 1, 
    name='inMod',
    col=c('white','royalblue'),
    column_split=col.split,
    border_gp=gpar(color='black', lwd=.5),
    use_raster=FALSE,
    width = ncol(fmat)*unit(ux, "mm"), 
    height = nrow(fmat)*unit(ux, "mm"),
)

col_fun = colorRamp2(c(-.025, 0, .025), c('blue', "white", 'red'))
htde = plotDEgenesHeatmap(ll$cmat[,cn], ll$pmat[,cn], row.split=ll$row.split,
                          col.split=col.split, ux=ux, cluster=FALSE)

h = 2 + 1 / 15 * nrow(ll$cmat)
w = 2 + 1 / 15 * ncol(ll$cmat)
pltprefix = paste0(imgpref, 'crossct_DEheatmap')
saveHeatmap(ht %v% htde, pltprefix, w=w, h=h)


# Add which region is highest + covar + celltype:
# -----------------------------------------------
dedf = regdedf[(regdedf$region != 'allregions') & (regdedf$gene %in% full.genes),]
dedf = dedf[dedf$path %in% pathlist,]
dedf = dedf[dedf$col_nm != 0,]
dedf = dedf[dedf$set %in% c(setlist),]
pind = dedf$path %in% c('cogdxad','nrad', 'braaksc.early', 'braaksc.ad')
dedf$logFC_nb[pind] = dedf$logFC_nb[pind] / 25 # Convert, more or less
dedf$abseff = abs(dedf$logFC_nb)

tdf = merge(dedf, aggregate(abseff ~ gene, dedf, max))
rownames(tdf) = tdf$gene
tdf = tdf[order(tdf$region),]
tdf = tdf[order(tdf$set),]
tdf = tdf[,c('gene','set', 'path','logFC_nb','region')]
tdf = merge(tdf, data.frame(gene=full.genes), all.y=TRUE)
tdf$logFC_nb[is.na(tdf$logFC_nb)] = 0
tdf[is.na(tdf)] = ''
rownames(tdf) = tdf$gene

setcols = c('red','purple','goldenrod1')
names(setcols) = setlist
pathcols = brewer.pal(length(pathlist), 'Set1')
names(pathcols) = pathlist

htbr = HeatmapAnnotation(set=as.character(tdf[cn,'set']),
                         path=as.character(tdf[cn,'path']),
                         logFC=tdf[cn,'logFC_nb'],
                         region=as.character(tdf[cn,'region']),
                         col=list(region=reg.cols,
                                  logFC=col_fun,
                                  set=setcols,
                                  path=pathcols),
                         gp=gpar(fontsize=5),
                         annotation_name_gp=gpar(fontsize=5),
                         simple_anno_size=unit(ux, 'mm')
)


ux = 1.5
col_fun = colorRamp2(c(-.025, 0, .025), c('blue', "white", 'red'))
htde = plotDEgenesHeatmap(ll$cmat[,cn], ll$pmat[,cn], row.split=ll$row.split,
                          col.split=col.split, ux=ux, cluster=FALSE, topann=htbr, hscale=0.75)

h = 2 + 1 / 15 * nrow(ll$cmat)
w = 5 + 1 / 15 * ncol(ll$cmat)
pltprefix = paste0(imgpref, 'crossct_DEheatmap_top')
saveHeatmap(ht %v% htde, pltprefix, w=w, h=h)



# Add DE results by region:
# -------------------------
path = 'plaq_d'
dedf = regdedf[(regdedf$gene %in% full.genes),]
dedf = dedf[dedf$path == path,]
dedf$path = dedf$region
dedf$set2 = paste0(dedf$set, '@', dedf$path)
setlist = c('Ast_Ast','Mic_Immune_Mic','Opc_Opc')
if (path %in% c('nrad','cogdxad', 'braaksc.early','braaksc.ad')){
    reglist = c('allregions','AG','MT','PFC','EC','HC','TH')
    col_fun = colorRamp2(25 * c(-.025, 0, .025), c('blue', "white", 'red'))
} else { 
    reglist = c('allregions','AG','MT','PFC','EC','HC')
    col_fun = colorRamp2(c(-.025, 0, .025), c('blue', "white", 'red'))
}
ll = getDEmatrices(dedf, reglist, setlist, genes=full.genes)

ux = 1.5
# col_fun = colorRamp2(c(-.025, 0, .025), c('blue', "white", 'red'))
htde = plotDEgenesHeatmap(ll$cmat[,cn], ll$pmat[,cn], row.split=ll$row.split,
                          col.split=col.split, ux=ux, cluster=FALSE, hscale=0.75)

h = 2 + 1 / 15 * nrow(ll$cmat)
w = 2 + 1 / 15 * ncol(ll$cmat)
pltprefix = paste0(imgpref, 'crossct_DEheatmap_regions')
saveHeatmap(ht %v% htde, pltprefix, w=w, h=h)




# DE effect only for genes difftl in diffuse plaque:
# --------------------------------------------------
# df = fulldedf[fulldedf$path == 'plaq_d',]
# degenes = unique(df$gene[df$col_nm != 0])
setlist = c('Ast_Ast','Mic_Immune_Mic','Opc_Opc')
pathlist = c('cogdxad','nrad','nft','plaq_n','plaq_d')

# List of genes that shows DE effects in any celltype 
path = 'plaq_d'
degenes = unique(fulldedf$gene[(fulldedf$col_nm != 0) & (fulldedf$set %in% setlist) & (fulldedf$path == path)])
core.genes = uqlist[['Core']]
fa.genes = c('ANGPTL4', 'HILPDA', 'IRS2','PRDM16','ADCY8')
glycogen.genes = c('GSK3B','PYGL','PPP1R3E', 'AGL', 'GBE1', 'UGP2', 'GYS1')
other.genes = c('GAPDH','HK1', 'HK2', 'SLC2A1','GFPT2', 'GLUL','SLC38A1','SLC38A2', 
    'EGLN3', 'DGKG', 'PLOD2', 'TPI1', 'PDK4','PKM')
core.genes = core.genes[!(core.genes %in% c(fa.genes,glycogen.genes,other.genes))]

fdf = data.frame(rbind(cbind(core.genes, 'Core'),
                       cbind(fa.genes, 'FA'),
                       cbind(glycogen.genes, 'Glycogen'),
                       cbind(other.genes, 'Other')))
names(fdf) = c('gene','set')

fdf = fdf[fdf$gene %in% degenes,]
full.genes = fdf$gene
col.split = fdf$set
rownames(fdf) = fdf$gene
fmat = sapply(extlist, function(x){ full.genes %in% x })
rownames(fmat) = full.genes
fmat = t(fmat)


# Get the DE results for this variable only:
# ------------------------------------------
dedf = fulldedf[(fulldedf$gene %in% full.genes) & (fulldedf$path == path),]
dedf$set2 = paste0(dedf$set, '@', dedf$path)
ll = getDEmatrices(dedf, path, setlist, genes=full.genes)

rmat = reord(t(ll$cmat))
cn = rownames(rmat)
col.split = fdf[cn, 'set']

ux = 1.5
ht = Heatmap(fmat[,cn] * 1, 
    name='inMod',
    col=c('white','royalblue'),
    column_split=col.split,
    border_gp=gpar(color='black', lwd=.5),
    use_raster=FALSE,
    width = ncol(fmat)*unit(ux, "mm"), 
    height = nrow(fmat)*unit(ux, "mm"),
)

mx = 0.01
col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
htde = plotDEgenesHeatmap(ll$cmat[,cn], ll$pmat[,cn], row.split=ll$row.split,
                          col.split=col.split, ux=ux, cluster=FALSE)

h = 2 + 1 / 15 * nrow(ll$cmat)
w = 2 + 1 / 15 * ncol(ll$cmat)
pltprefix = paste0(imgpref, 'crossct_DEheatmap.', path)
saveHeatmap(ht %v% htde, pltprefix, w=w, h=h)



# Plot DE effect + where it is strongest for these genes
# ------------------------------------------------------


# Plot regional differences (enriched where?)


