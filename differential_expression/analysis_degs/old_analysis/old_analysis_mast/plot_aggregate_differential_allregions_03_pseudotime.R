#!/usr/bin/R
# -----------------------------------------------------------
# Aggregate and plot chunked runs of MAST + RE for DGE
# Plot the comparison between regressions on different pathologies.
# Updated: 02/11/2021
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(Matrix)

library(MAST)
library(data.table)

library(viridis)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

celltype = 'Ast'
subtype = 'Ast'
path = 'nrad'

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype subtype region chunksize chunk ascertainment")
} else {        
    celltype = args[1]
    subtype = args[2]
    region = args[3]
    path = args[4]
}

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# ------------------
# Load the metadata:
# ------------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

# Colors for full:
typelvls = unique(cellmeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
type.cols = c(type.cols, major.col['Inh'], major.col['Exc'])
tsp.type.cols = sapply(type.cols, tsp.col)

# Data directories:
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    # matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
    matdir = paste0(datadir,'matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
mtxdir = paste0(matdir, 'mtx/')
ststr = gsub("/","_",gsub(" ","_", subtype))
cellstr = gsub("/","_",gsub(" ","_", celltype))

# --------------------------------------------
# Load from various different regression runs:
# --------------------------------------------
agg.pref = paste0(regdir, prefix, '.mastlmm_reg.allpath.allreg.major.', celltype, '.minor.', ststr)
agg.file = paste0(agg.pref, '.tsv.gz') 
agg.rda = paste0(agg.pref, '.rda') 
if (!file.exists(agg.rda)){
    alldf = c()
    sumdf = c()
    for (path in c('nrad','nft','plaq_d','plaq_n')){
        print(paste("[STATUS] Loading regression on", subtype, 'across regions', 'on var', path))
        fpref = paste0(prefix, '.mastlmm_reg.', path, '.allreg.major.', celltype, '.minor.', ststr)
        # fnlist = list.files(pattern=paste0(fpref,'.*633.*.Rda'), path=regdir)
        fnlist = list.files(pattern=paste0(fpref,'.*.160.*.Rda'), path=regdir)
        if (length(fnlist) > 0){
            print(paste("[STATUS] Loading from",length(fnlist),"files."))
            subdf = c()
            for (fn in fnlist){
                # print(fn)
                load(paste0(regdir,fn))
                regdf$path = path 
                summaryDt$path = path 
                subdf = rbind(subdf, regdf)
                sumdf = rbind(sumdf, summaryDt)
            }
            names(subdf)[2] = 'pvalue'
            subdf = subdf[order(subdf$pvalue),]
            subdf$padj = p.adjust(subdf$pvalue, 'fdr')
            alldf = rbind(alldf, subdf)
        } else {
            print("No files..")
        }
    }
    dim(alldf)
    write.table(alldf, file=gzfile(agg.file), quote=F, sep="\t", col.names=F)
    save(alldf, file=agg.rda)
} else { 
    load(agg.rda)
}

FCTHRESH=0.02 
alldf$col = 1 * (alldf$padj < 0.05) * (1 + 1 * (alldf$coef > 0)) * (abs(alldf$coef) > ifelse(alldf$path == 'nrad', 1,1/50) * FCTHRESH)

# ------------------------------
# Load in data from all regions:
# ------------------------------
deg = unique(alldf$primerid[alldf$col != 0])
keep.reg = regions[regions != 'MB']
amat = c()
barcodes = c()
for (region in keep.reg){
    ststr = gsub("/","_",gsub(" ","_", subtype))
    matpref = paste0(mtxdir, rawpref,'.majorcelltype.',
                     celltype,'.',ststr,'.',region)
    rdafile = paste0(matpref, '.rda')  # In Matrix format
    if (file.exists(rdafile)){
        # Load `mat` from rdafile:
        load(rdafile)
        deg = deg[deg %in% rownames(mat)]
        if (subtype %in% c('Oli','Exc')){
            mat = mat[deg,]
        }
        print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
        barcodes = c(barcodes, colnames(mat))
        genes = rownames(mat)
        ngenes = nrow(mat)
        amat = cbind(amat, mat)
    }
}
rm(mat)
gcout = gc()

# Margin for norm:
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

# ----------------------------------------
# Get the pathology mapped to each region:
# ----------------------------------------
pathlist = c('nft','plaq_d','plaq_n')
regmap = c('AG','HC','PFC','MT','EC')
names(regmap) = c('ag','hip','mf','mt','ec')
pqdf = NULL
for (path in pathlist){
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    sdf = slong[,c('rind','value','region')]
    names(sdf)[2] = paste0(path, "_act")
    if (is.null(pqdf)){
        pqdf = sdf
    } else {
        pqdf = merge(pqdf, sdf)
    }
}
tdf = NULL

# ----------------------------
# Load in the pseudotime data:
# ----------------------------
use.prreg = FALSE
use.meld = TRUE
pathlist = c('nft','plaq_d','plaq_n')
astr = paste0('allreg.major.', celltype, '.minor.', ststr)
if (use.meld){
    astr = paste0(astr, '.meld')
    impdir = 'multiRegion/meld_scores/'
    df = read.delim(paste0(impdir, celltype,'.',subtype,'.MELD_metadata.tsv'), header=T)
    idf = df[,c('barcode',paste0('bin_',pathlist,'__pseudotime'))]
    tdf = df[,c('barcode',paste0('bin_',pathlist,'__n_clusters.5'))]
    names(idf) = c('barcode',pathlist)
    names(tdf) = c('barcode',paste0(pathlist, '.cls'))
    rm(df)
} else { 
    impdir = 'multiRegion/rw_top_imputed_scores/'
    ibc = scan(paste0(impdir, celltype, '_barcodes.txt'),'c')
    idf = data.frame(barcode=ibc)
    if (use.prreg){
        perreg = 'per_region_True_' 
        astr = paste0(astr, '.perreg')
    } else {
        perreg = 'per_region_False_'
        astr = paste0(astr, '.crossreg')
    }
    for (ipath in pathlist){
        isc = as.numeric(scan(paste0(impdir, celltype, '_imputed_scores_',perreg,ipath,'.txt'),'c'))
        # isc = scale(isc)
        idf[[ipath]] = isc
    }
}
gcout = gc()

idf = idf[idf$barcode %in% colnames(amat),]
idf = merge(idf, cellmeta[,c('barcode','region','cell_type_high_resolution','rind','U1','U2')])
idf = merge(idf, pqdf, all.x=TRUE)

NBIN = 75
nftbks = seq(min(idf$nft), max(idf$nft), length.out=NBIN)
pnbks = seq(min(idf$plaq_n), max(idf$plaq_n), length.out=NBIN)
idf$nftcut = cut(idf$nft, breaks=nftbks)
idf$pncut = cut(idf$plaq_n, breaks=pnbks)
idf$nid = as.numeric(idf$nftcut)
idf$pid = as.numeric(idf$pncut)

gplot = ggplot(idf, aes(x=nft, y=plaq_n)) +
# ggplot(idf, aes(x=plaq_d, y=plaq_n)) +
    geom_hex(bins=100) +
    geom_density_2d(color='white', lwd=.1) + 
    scale_fill_continuous(type='viridis')+ 
    theme_bw()
ggsave(paste0(imgpref, 'pseudo2d.',astr,'.fullcontour.png'), gplot, dpi=300, units='in', width=7, height=6)


gplot = ggplot(idf, aes(x=nft, y=plaq_n)) +
    facet_grid(cell_type_high_resolution~region) + 
    geom_hex(bins=100) +
    geom_density_2d(color='white', lwd=.1) + 
    scale_fill_continuous(type='viridis')+ 
    theme_bw()
ggsave(paste0(imgpref, 'pseudo2d.',astr,'.regioncontour.png'), gplot, dpi=300, units='in', width=12, height=8)


bdf = aggregate(barcode ~ nftcut + pncut, idf, length)
gplot = ggplot(bdf, aes(x=nftcut, y=pncut, z=barcode, color=log(barcode))) +
    geom_point(pch=15, cex=6) + 
    scale_color_continuous(type='viridis')+ 
    theme_bw()

if (!is.null(tdf)){
    bdf = merge(idf, tdf)
    gplot = ggplot(bdf, aes(x=nft, y=plaq_n, color=factor(nft.cls))) +
        geom_point() + 
        geom_density_2d(color='white', lwd=.1) + 
        scale_fill_continuous(type='viridis')+ 
        theme_bw()
    ggsave(paste0(imgpref, 'pseudo2d.',astr,'.meld_nftcls.png'), gplot, dpi=300, units='in', width=12, height=8)
}

for (path in pathlist){
    print(path)
    idf$norm = log(idf[[paste0(path, '_act')]] + 1)
    # Aggregate + make matrix:
    bdf = aggregate(norm ~ nid + pid, idf, mean)
    bmat = matrix(NA, nrow=NBIN, ncol=NBIN)
    dfind = as.matrix(bdf[,c('nid','pid')])
    bmat[dfind] = bdf$norm

    # Plot matrix:
    h = 5
    png(paste0(imgpref, 'pseudo2d.',astr,'.',path,'.png'), res=400, units='in', width=h, height=h)
    sp = 0.25
    par(mar=c(2.25,2.25,sp,sp))
    image(bmat, col=viridis(100), useRaster=T, axes=F)
    mtext('NFT Pseudotime',side=1, line=1)
    mtext('Neuritic Plaque Pseudotime',side=2, line=1.4)
    mtext(path, side=3, line=-1)
    axis(1, at=c(0,.5,1), labels=FALSE, tck=-.01)
    axis(2, at=c(0,.5,1), labels=FALSE, tck=-.01)
    text(x=c(0,0.5,1), y=parpos(2,0.035), label=c(0,0.5,1), xpd=TRUE, adj=.5)
    text(y=c(0,0.5,1), x=parpos(1,0.015), label=c(0,0.5,1), xpd=TRUE, adj=1)
    dev.off()
}


# --------------------------------------
# Plot specific genes on the pseudotime:
# --------------------------------------
gene = 'MT2A'
if (subtype == "End"){
    genes = c('LYVE1','CTGF','SDCBP','MT2A','HSP90AA1','GPCPD1','RHOJ','VEGFA','LDHA', 'LINGO1', 'MT-ND3','HLA-C','ANGPT2', 'PLAUR', 'CD34','WNT2B','HSPA1A','NLGN1','PTMS', 'HSPG2','TUSC3','EDN1')
} else if (subtype == "Per"){
    genes = c('APOD','PCBP3','PRELP','SLC6A1','HSP90AA1','DMD','MXI1','DUSP1')
} else if (subtype == "Opc"){
    genes=c('GPC6','MT1X','FKBP5','CEBPD','APOD','CD81','NRCAM','DUSP1','HIF3A','RGS7','VEGFA', 'ANGPT2','FGF14','VCAN')
} else if (subtype == "Ast"){
    genes=c('CDH20','GFAP','CABLES1','HSP90AA1','SLC39A11','FTH1','DCLK1','DPP10','DGKB','RGMA','GPC5', 'SPARCL1','HMGB1','SDC4')
} else if (subtype == "Mic"){
    genes=c('APOE','FKBP5','IFI44L','XAF1','INPP5D','GLDN','HIF1A','HLA-DRB1','MT-ND3','DPYD','CD74', 'CX3CR1','CD81','SERPINE1', 'LRRK1','TLR2', 'UBC', 'NRP2')
} else if (subtype == "CAMs"){
    genes=c('APOE','FKBP5','IFI44L','SPP1','HS3ST4','OXR1','CD74','CD81','SLC11A1','SIGLEC1','MT-ND3', 'SRGN','FTL','HLA-DRB1', 'LIPA','SLC1A3', 'F13A1','LGMN')
} else if (subtype == "Oli"){
    genes=c('QDPR','ANLN','SLC39A11','SLC44A1','SPP1','FKBP5','CPQ','MAN2A1','PDE1C','PDE1A','GLDN', 'DLG2','CDH20','PLCG2', 'ERBIN','MT-ND3', 'CRYAB','ST18','INSR', 'SLC5A11')
}
if (subtype == 'Mic'){
    idf = idf[grep("Mic",idf$cell_type_high_resolution),]
}

for (gene in genes){
    print(gene)
    # Get gene:
    mx = marg[colnames(amat),'count']
    x = amat[gene,]
    nx = log(1e3 * (x / mx) + 1)
    idf$val = x[idf$barcode]
    idf$norm = nx[idf$barcode]
    # Aggregate + make matrix:
    bdf = aggregate(norm ~ nid + pid, idf, mean)
    bmat = matrix(NA, nrow=NBIN, ncol=NBIN)
    dfind = as.matrix(bdf[,c('nid','pid')])
    bmat[dfind] = bdf$norm

    # Plot matrix:
    h = 5
    png(paste0(imgpref, 'pseudo2d.',astr,'.',gene,'.png'), res=400, units='in', width=h, height=h)
    sp = 0.25
    par(mar=c(2.25,2.25,sp,sp))
    image(bmat, col=viridis(100), useRaster=T, axes=F)
    mtext('NFT Pseudotime',side=1, line=1)
    mtext('Neuritic Plaque Pseudotime',side=2, line=1.4)
    mtext(gene,side=3, line=-1)
    axis(1, at=c(0,.5,1), labels=FALSE, tck=-.01)
    axis(2, at=c(0,.5,1), labels=FALSE, tck=-.01)
    text(x=c(0,0.5,1), y=parpos(2,0.035), label=c(0,0.5,1), xpd=TRUE, adj=.5)
    text(y=c(0,0.5,1), x=parpos(1,0.015), label=c(0,0.5,1), xpd=TRUE, adj=1)
    dev.off()

    # 2x2 plot of smoothed fits of actual/imputed; nft/plaqn
    g1 = ggplot(idf, aes(plaq_n,norm, col=cell_type_high_resolution)) + 
        geom_smooth(method='gam') + 
        scale_color_manual(values=type.cols, name='Subtype') + 
        scale_y_continuous(expand=c(0,0)) + 
        labs(x='Neuritic Plaque Pseudotime', y='Norm. Expression, log(x + 1)', title=gene) + 
        theme_pubr()
    g2 = ggplot(idf, aes(log(plaq_n_act +1),norm, col=cell_type_high_resolution)) + 
        geom_smooth(method='gam') + 
        scale_color_manual(values=type.cols, name='Subtype') + 
        scale_y_continuous(expand=c(0,0)) + 
        labs(x='Neuritic Plaque Meas. (log)', y='Norm. Expression, log(x + 1)', title=gene) + 
        theme_pubr()
    g3 = ggplot(idf, aes(nft,norm, col=cell_type_high_resolution)) + 
        geom_smooth(method='gam') + 
        scale_color_manual(values=type.cols, name='Subtype') + 
        scale_y_continuous(expand=c(0,0)) + 
        labs(x='NFT Pseudotime', y='Norm. Expression, log(x + 1)', title=gene) + 
        theme_pubr()
    g4 = ggplot(idf, aes(log(nft_act+1),norm, col=cell_type_high_resolution)) + 
        geom_smooth(method='gam') + 
        scale_color_manual(values=type.cols, name='Subtype') + 
        scale_y_continuous(expand=c(0,0)) + 
        labs(x='NFT Measurement (log)', y='Norm. Expression, log(x + 1)', title=gene) + 
        theme_pubr()
    garr = ggarrange(g1, g2, g3, g4, align='hv',
                     common.legend=TRUE, ncol=2, nrow=2)
    ggsave(paste0(imgpref, 'pseudo2d.',astr,'.',gene, '.gam_1d.png'), garr, dpi=300, units='in', width=8, height=8)
}


# Diff by set:
norm = sweep(amat, 2, marg[colnames(amat),'count'] / 100000,'/')
norm = log(norm + 1)
gcout = gc()

# gene = 'MT2A'
tform = make.tform(paste0('C', tdf$nft.cls), norm=T)
avg.mat = norm[,tdf$barcode] %*% tform
smat = t(scale(t(avg.mat)))

resdf = alldf[alldf$path == 'nft',]
topgenes = head(resdf$primerid, 75)
Heatmap(as.matrix(smat[topgenes,]))



# -------------------------------------------
# Get human PPI data, partially based off of:
# -------------------------------------------
# https://www.bioconductor.org/packages/release/bioc/vignettes/netSmooth/inst/doc/buildingPPIsFromStringDB.html
thresh = 400
version = '11'
string.ppi.rda = paste0('Annotation/string_ppi_v',version,'_thresh',thresh,'_dataframe.Rda')
string.ppi.tsv = paste0('Annotation/string_ppi_v',version,'_thresh',thresh,'_dataframe.tsv.gz')
if (!file.exists(string.ppi.rda)){
    require(STRINGdb)
    require(igraph)
    require(biomaRt)
    string_db <- STRINGdb$new(species=9606, version=version,score_threshold=thresh)
    human_graph <- string_db$get_graph()
    adj_matrix <- as_adjacency_matrix(human_graph)
    pdf = as.data.frame(summary(adj_matrix))
    rdf = data.frame(i = 1:nrow(adj_matrix), 
                     ensembl_peptide_id=sapply(strsplit(rownames(adj_matrix), '\\.'),
                                               function(x) x[2]))
    # Map gene ids to protein ids
    mart=useMart(host = 'grch37.ensembl.org',
                 biomart='ENSEMBL_MART_ENSEMBL',
                 dataset='hsapiens_gene_ensembl')
    mart_results <- getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"),
                          filters = "ensembl_peptide_id", values = rdf$ensembl_peptide_id,
                          mart = mart)
    rdf = merge(rdf, mart_results, all.x=TRUE)
    fdf = unique(rdf[,c('i','hgnc_symbol')])
    tdf = fdf
    names(fdf) = c('i','from')
    names(tdf) = c('j','to')
    pdf = merge(merge(pdf, fdf), tdf)
    pdf = pdf[!is.na(pdf$from) & !is.na(pdf$to),]
    pdf = aggregate(x ~ from + to, pdf, max)
    ppi.df = pdf
    rm(human_graph, adj_matrix, pdf); gc()
    # Save:
    write.table(ppi.df, file=gzfile(string.ppi.tsv), quote=F, sep="\t", col.names=F)
    save(ppi.df, file=string.ppi.rda)
} else { 
    load(string.ppi.rda)
}


# -----------------------------------------------
# Look at these genes in the context of STRINGdb:
# -----------------------------------------------
# All aggregated results:
fns = list.files(path=regdir, pattern=paste0(prefix,'.mastlmm_reg.allpath.allreg.major.*.rda'))
fulldf = c()
for (fn in fns){
    basename = sub("^.*.mastlmm_reg.allpath.allreg.major.","",fn)
    celltype = sub("\\..*$","",basename)
    subtype = sub(".rda","",sub("^.*minor.","",basename))
    load(paste0(regdir, fn))
    if (!is.null(alldf)){
        alldf = alldf[,c('primerid','coef','padj','path')]
        alldf$celltype = celltype
        alldf$subtype = subtype
        fulldf = rbind(fulldf, alldf)
    }
}

# Only one pathology at a time:
path = 'plaq_n'
resdf = fulldf[fulldf$path == path,]

upgene = resdf[resdf$col == 2,'primerid']
downgene = resdf[resdf$col == 1,'primerid']

udf = resdf[resdf$col == 2,c('padj','coef','primerid')] 
names(udf) = c('pvalue','coef','gene')


# Look at specific genes:
ppi.df[ppi.df$to %in% c('SEMA3F') | ppi.df$from %in% c('SEMA3F'),]


### replace protein ids with gene ids
ix <- match(protein_ids, mart_results$ensembl_peptide_id)
ix <- ix[!is.na(ix)]

newnames <- protein_ids
newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <- mart_results[ix, 'hgnc_symbol']
rownames(adj_matrix) <- newnames
colnames(adj_matrix) <- newnames

ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
nullrows <- Matrix::rowSums(ppi)==0
ppi <- ppi[!nullrows,!nullrows]

upi = ppi[up1,]
upi = upi[,colSums(upi) != 0]


up1 = upgene[upgene %in% rownames(ppi)]
dw1 = downgene[downgene %in% rownames(ppi)]

ppi.todf = function(xpi, keep.diag = FALSE){
    xdf = as.data.frame(summary(xpi))
    xdf$from = rownames(xpi)[xdf$i]
    xdf$to = colnames(xpi)[xdf$i]
    if (!keep.diag){ xdf = xdf[xdf$to != xdf$from,] }
    if (nrow(xdf) > 1){ xdf = aggregate(x ~ from + to, xdf, max) }
    if (nrow(xdf) > 1){ xdf = xdf[order(xdf$from),] }
    return(xdf[,c('from','to','x')])
}

upi = 

updf = ppi.todf(ppi[up1,])
dwdf = ppi.todf(ppi[dw1,])











