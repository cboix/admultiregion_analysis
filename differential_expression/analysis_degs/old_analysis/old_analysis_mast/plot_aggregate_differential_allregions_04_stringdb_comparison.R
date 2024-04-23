#!/usr/bin/R
# -----------------------------------------------------------
# Plot/predict ligand interactions with string given 
# the top DE genes in each CT
# Updated: 03/23/2021
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

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'stringcomp_')
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

pqdf2 = merge(pqdf, metadata[,c('rind','projid')])
write.table(pqdf2, 'Annotation/multiregion_path_by_region_projid.tsv', 
            quote=F, sep="\t", row.names=F)

# -------------------------------------------
# Get human PPI data, partially based off of:
# -------------------------------------------
# https://www.bioconductor.org/packages/release/bioc/vignettes/netSmooth/inst/doc/buildingPPIsFromStringDB.html
thresh = 800
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
        cat(celltype, '\t', subtype, '\t', nrow(alldf), '\n')
        alldf = alldf[,c('primerid','coef','padj','path')]
        alldf$celltype = celltype
        alldf$subtype = subtype
        fulldf = rbind(fulldf, alldf)
    }
}
FCTHRESH=0.02 
fulldf$col = 1 * (fulldf$padj < 0.05) * (1 + 1 * (fulldf$coef > 0)) * (abs(fulldf$coef) > ifelse(fulldf$path == 'nrad', 1,1/50) * FCTHRESH)

# Only one pathology at a time:
path = 'plaq_n'
fulldf = fulldf[!is.na(fulldf$primerid),]
resdf = fulldf[fulldf$path == path,]
resdf = resdf[order(resdf$padj),]
sigdf = resdf[resdf$col > 0,]
# nsigdf = table(sigdf[,c('subtype','col')])

# Get all interactions at this level:
sts = unique(fulldf$subtype)
combdf = expand.grid(S1=sts, S2=sts)
intdf = c()
for (i in 1:nrow(combdf)){
    u1 = sigdf$primerid[sigdf$subtype == combdf$S1[i]]
    u2 = sigdf$primerid[sigdf$subtype == combdf$S2[i]]
    u1 = u1[u1 %in% ppi.df$from]
    u2 = u2[u2 %in% ppi.df$to]
    # u1 = head(u1,250)
    # u2 = head(u2,250)
    # Basic stratification of interactions into putative inter/intra-cell
    # TODO: Could improve on this with the scores + directionalities:
    udf = ppi.df[ppi.df$from %in% u1 & ppi.df$to %in% u2,]
    if (nrow(udf) > 0){
        i1df = ppi.df[ppi.df$from %in% u1 & ppi.df$to %in% u1,]
        i2df = ppi.df[ppi.df$from %in% u2 & ppi.df$to %in% u2,]
        idf = unique(rbind(i1df, i2df))
        idf$status = 'Intra'
        udf = merge(udf, idf, all.x=TRUE)
        udf$status[is.na(udf$status)] = 'Inter'
        # table(udf$status)
        udf$S1 = combdf$S1[i]
        udf$S2 = combdf$S2[i]
        udf$N1 = length(u1)
        udf$N2 = length(u2)
        intdf = rbind(intdf, udf)
    }
}
table(intdf[,c('S1','S2','status')])



statdf = merge(data.frame(table(intdf[,c('S1','S2','status')])), 
               unique(intdf[,c('S1','S2','N1','N2')]))

ggplot(statdf[statdf$status == 'Inter',], aes(S1,S2, color=Freq/N2)) + 
    geom_point(pch=15, size=20) + 
    scale_color_viridis() + 
    theme_pubr()

# Find highly enriched genes?? 
gstat = agg.rename(x ~ from, ppi.df, length, 'nconn')
gdf = agg.rename(x ~ from, intdf[intdf$status == 'Inter',], length, 'npull')
gdf = merge(gdf, gstat)
gdf$frac = gdf$npull / gdf$nconn
gdf = gdf[order(gdf$frac, decreasing=T),]


# Cell polarity:
sigdf[sigdf$primerid %in% c('INTU','WDPCP'),]

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


# ----------------------------------------------
# Load in all data from EC, check that it is ok:
# ----------------------------------------------
# Keep markers + gene sets we care about:
m1df = read.delim('Annotation/broad_ct_markers.tsv', header=T, stringsAsFactors=F)
m2df = read.delim('Annotation/known_markers.tsv', header=T, stringsAsFactors=F)
gnr =  c('NRCAM','NFASC','SCN8A','SCN9A','KCNQ3','CNTN1','CNTNAP2','ANK3','ADAM22','DLG1','DLG2')
kept.genes = unique(c(m1df$gene, m2df$symbol, gnr, 'CLU','BIN1','APOE'))

# Load data:
fns = list.files(path=mtxdir, pattern=paste0(rawpref, '.majorcelltype..*.EC.rda'))
region = 'EC'
amat = c()
sts = c(); barcodes = c()
for (fn in fns){
    basename = sub("^.*.majorcelltype.","",fn)
    celltype = sub("\\..*$","",basename)
    subtype = sub("^.*\\.","",sub(".EC.rda","",basename))
    # Load `mat` from rdafile:
    load(paste0(mtxdir, fn))
    mat = mat[kept.genes,]
    print(paste("[STATUS] Loaded", subtype, 'in',region,'with',ncol(mat), 'cells'))
    barcodes = c(barcodes, colnames(mat))
    sts = c(sts, rep(subtype, ncol(mat)))
    genes = rownames(mat)
    ngenes = nrow(mat)
    amat = cbind(amat, mat)
}
rm(mat)
gcout = gc()

# Normalize first:
mx = marg[barcodes,'count']
norm = sweep(amat, 2, mx / 10000, '/')
norm = log(norm + 1)

# Select neuronal subtypes in EC only:
rdf = agg.rename(barcode ~ region + cell_type_high_resolution, cellmeta[cellmeta$major.celltype == 'Exc',], length, 'count')
rdf = spread(rdf, region, count, fill=0)
rmat = as.matrix(rdf[,-1])
rownames(rmat) = rdf[,1]
topec = apply(rmat,1, max) == rmat[,'EC']
topec = names(topec)[topec]
topec = topec[topec != 'Exc SV2C LINC02137']

# Averages, by individual:
cellmeta$cr = paste0(cellmeta$minor.celltype, '_', cellmeta$rind)
cellmeta$chr = paste0(cellmeta$cell_type_high_resolution, '.', cellmeta$rind)
tind = cellmeta[barcodes, 'cell_type_high_resolution'] %in% topec
rownames(cellmeta) = cellmeta$barcode
crlist = cellmeta[barcodes, 'cr']
chrlist = cellmeta[barcodes[tind], 'chr']
ctlist = cellmeta[barcodes, 'minor.celltype']
tform = make.tform(crlist, norm=T)
tform2 = make.tform(ctlist, norm=T)
tform3 = make.tform(chrlist, norm=T)
kgenes = unique(c(m1df$gene, gnr))
avg.mat = as.matrix(norm[kgenes,] %*% tform)
avg.mat2 = as.matrix(norm[kgenes,] %*% tform2)
avg.mat3 = as.matrix(norm[kgenes,tind] %*% tform3)
# TODO: Also look at normalized
# TODO: Diagonalize, ordering celltypes + rows based on 

# Reorder genes:
smat = t(scale(t(avg.mat)))
smat2 = t(scale(t(avg.mat2)))
smat3 = t(scale(t(avg.mat3)))
rmat = smat2[, sort(colnames(smat2))]
rmat = reord(rmat)
rmat = diag.mat2(t(rmat))[[1]]
smat = smat[colnames(rmat),]
smat3 = smat3[colnames(rmat),]
avg.mat3 = avg.mat3[colnames(rmat),]
avg.mat = avg.mat[colnames(rmat),]

# Annotate:
clsplit = sapply(strsplit(colnames(smat), '_'), function(x){x[1]})
rind = sapply(strsplit(colnames(smat), '_'), function(x){x[2]})
ha = HeatmapAnnotation(NFT=metadata[rind,'nft_ec'], 
                       PLQ=metadata[rind,'plaq_n_ec'], 
                       Braak=metadata[rind,'Braak'], 
                       AD=metadata[rind,'niareagansc'], 
                       col=list(AD=colvals[['niareagansc']],
                                Braak=colvals[['braaksc']],
                                AD.Diff=c('Depleted'='slateblue',
                                          'No Change'='grey75')))

png(paste0(imgpref, 'marker_nranvier_genes_ctind_EC_heatmap.png'), res=300, units='in', width=10, height=6)
Heatmap(avg.mat, name='Norm.\nexpr',
        col=viridis(100),
        row_split = ifelse(rownames(smat) %in% gnr, 'NodeOfRanvier','MarkerGenes'),
        column_split = clsplit,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE,
        top_annotation=ha, 
        # right_annotation=hb
)
dev.off()

# Overall for the glia: 
colnames(avg.mat) = gsub(" ","_", colnames(avg.mat) )
adf = data.frame(as.matrix(avg.mat))
adf$gene = rownames(adf)
adf = gather(adf, cr, val, -gene)
adf$ct = sapply(strsplit(adf$cr, '_'), function(x){x[1]})
adf$rind = sapply(strsplit(adf$cr, '_'), function(x){x[2]})
adf = adf[adf$gene %in% gnr,]
adf = merge(adf, unique(metadata[,c('rind','msex','niareagansc','cogdx','braaksc','nft_ec','plaq_n_ec','plaq_d_ec')]))
adf = adf[adf$ct %in% c('Ast','CAM','End','Exc','Inh','Mic','Oli','Opc'),]
# adf$vuln = ifelse(adf$ct %in% c('Exc_TOX3_TTC6','Exc_RELN_COL5A2','Exc_RELN_GPC5','Exc_AGBL1_GPC5'), 'Vuln.','Stable')
adf$nrad = 'AD'
adf$nrad[adf$niareagansc > 2] = 'CTRL'
adf$nrad = factor(adf$nrad, levels=c('CTRL', 'AD'))

ggplot(adf, aes(plaq_n_ec, val, color=ct)) + 
    # facet_grid(gene~ct, scale='free_y') +
    facet_wrap(~gene, scale='free_y') +
    geom_point() + 
    geom_smooth(method='lm') + 
    # scale_color_manual(values=c('grey70','slateblue')) + 
    theme_pubr()

gplot = ggplot(adf, aes(ct, val, col=nrad)) + 
    facet_wrap(~gene, scale='free_y') +
    geom_boxplot(outlier.size=.5) + 
    scale_color_manual(values=colvals[['nrad']]) + 
    labs(x='Celltype',y='Normalized Expr.') +
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(imgpref, 'marker_nranvier_genes_ctind_EC_boxplots.png'), gplot, dpi=300, units='in', width=8, height=6)

# Fit regression across all of these:
library(emmeans)
fit = lm(val ~ gene *(nrad * ct), adf)
emform = asform(c('revpairwise ~ nrad|gene *ct'))
emm1 <- emmeans(fit, specs=emform)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> crdf 
crdf = crdf[order(crdf$p.value),]
head(crdf[,c('gene','ct','estimate','p.value')])

# -----------------------------------
# Annotate the neuronal subtypes too:
# -----------------------------------
clsplit = sapply(strsplit(colnames(avg.mat3), '\\.'), function(x){x[1]})
rind = sapply(strsplit(colnames(avg.mat3), '\\.'), function(x){paste0(x[2],'.',x[3])})
ha = HeatmapAnnotation(NFT=metadata[rind,'nft_ec'], 
                       PLQ=metadata[rind,'plaq_n_ec'], 
                       Braak=metadata[rind,'Braak'], 
                       AD=metadata[rind,'niareagansc'], 
                       col=list(AD=colvals[['niareagansc']],
                                Braak=colvals[['braaksc']],
                                AD.Diff=c('Depleted'='slateblue',
                                          'No Change'='grey75')))

Heatmap(avg.mat3,
        row_split = ifelse(rownames(avg.mat3) %in% gnr, 'NodeOfRanvier','MarkerGene'),
        column_split = clsplit,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE,
        top_annotation=ha, 
        # right_annotation=hb
)


# Overall for the neurons
colnames(avg.mat3) = gsub(" ","_", colnames(avg.mat3) )
adf = data.frame(as.matrix(avg.mat3))
adf$gene = rownames(adf)
adf = gather(adf, cr, val, -gene)
adf$ct = sapply(strsplit(adf$cr, '\\.'), function(x){x[1]})
adf$rind = sapply(strsplit(adf$cr, '\\.'), function(x){paste0(x[2],'.', x[3])})
adf = adf[adf$gene %in% gnr,]
adf = merge(adf, unique(metadata[,c('rind','msex','niareagansc','cogdx','braaksc','nft_ec','plaq_n_ec','plaq_d_ec')]))
adf$vuln = ifelse(adf$ct %in% c('Exc_TOX3_TTC6','Exc_RELN_COL5A2','Exc_RELN_GPC5','Exc_AGBL1_GPC5'), 'Vuln.','Stable')
adf$nrad = 'AD'
adf$nrad[adf$niareagansc > 2] = 'CTRL'
adf$nrad = factor(adf$nrad, levels=c('CTRL', 'AD'))

ggplot(adf, aes(plaq_n_ec, val, color=vuln)) + 
    facet_wrap(~gene, scale='free_y') +
    geom_point() + 
    geom_smooth(method='lm') + 
    scale_color_manual(values=c('grey70','slateblue')) + 
    theme_pubr()

gplot = ggplot(adf, aes(ct, val, col=nrad)) + 
    facet_wrap(~gene, scale='free_y') +
    geom_boxplot(outlier.size=.5) + 
    scale_color_manual(values=colvals[['nrad']]) + 
    labs(x='Celltype',y='Normalized Expr.') +
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(imgpref, 'marker_nranvier_genes_neurind_EC_boxplots.png'), gplot, dpi=300, units='in', width=8, height=6)

fit = lm(val ~ gene *(nrad * vuln), adf)
emform = asform(c('revpairwise ~ nrad|gene *vuln'))
emm1 <- emmeans(fit, specs=emform)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> crdf 
crdf = crdf[order(crdf$p.value),]
head(crdf[,c('gene','ct','estimate','p.value')], 20)
head(crdf[,c('gene','vuln','estimate','p.value')], 20)










# Overall:
adf = data.frame(as.matrix(avg.mat))
adf$gene = rownames(adf)
adf = gather(adf, cr, val, -gene)
adf$ct = sapply(strsplit(adf$cr, '_'), function(x){x[1]})
adf$rind = sapply(strsplit(adf$cr, '_'), function(x){x[2]})







# Try instead the strained data:

if (reg == 'EC'){
    stfile = paste0(datadir, 'EC_strainedCounts_raw.hdf5')
    h5f = H5Fopen(stfile)
    genes = h5f$genes
    bcs = h5f$barcodes
    ind = which(genes %in% kg)
    bind = match(sub.pathdf$barcode, bcs)
    h5d = h5f&"matrix"
    mat2 = t(h5d[ind,bind])
    H5Dclose(h5d)
    H5Fclose(h5f)
    colnames(mat2) = genes[ind]
    rownames(mat2) = bcs[bind]
    # Order as model matrix + transpose matrix:
    mat2 = mat2[sub.pathdf$barcode,]
    mat2 = t(Matrix(mat2))
    gcout = gc()
    # Margin:
    stmfile = paste0(datadir, 'EC_strainedCounts_margin.tsv.gz')
    stmarg = read.delim(gzfile(stmfile), header=F)
    names(stmarg) = 'count'
    stmarg$barcode = bcs
    rownames(stmarg) = stmarg$barcode
    offset2 = stmarg[sub.pathdf$barcode, 'count']
}












