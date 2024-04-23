#!/usr/bin/R
# -----------------------------------------------------------
# Consistency of DEG directionality in the neuronal subtypes:
# Updated: 12/20/22
# --------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ggplot2)
library(ggpubr)
print(version)
options(width=170)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Get a list of all differential runs:
# -------------------------------------------------------
rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_runlist.tsv'), header=T)
rundf = rbind(rundf, read.delim(paste0(sdbdir, 'nebula_wRUV_excitatory_subsets_runlist.tsv'), header=T))
rundf$prefstr = with(rundf, paste(celltype, subtype, region, path, sep="_"))

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
                            length(list.files(path=regdir, pattern=paste0('allmethods.', x, '.merged.rda'))) })
table(rundf$merged)
head(rundf[rundf$merged == 0,])


# Select runs to use (excitatory subtypes only):
# ----------------------------------------------
pathlist = c('nft', 'plaq_n', 'plaq_d','cogdxad','nrad')
selrundf = rundf[rundf$path %in% pathlist,]
selrundf$setid = with(selrundf, paste0(celltype, '_', subtype))
sets = unique(selrundf$setid)
sets = sets[grep("^Exc_", sets)]
sets = sets[sets != 'Exc_Exc']
selrundf = selrundf[selrundf$setid %in% sets,]


# Aggregate the DEGs and results across all regional runs:
# --------------------------------------------------------
mstr = paste0('allmethods.exc_sets.regional')
fullaggrda = paste0(regdir, mstr, '.merged.rda')
if (!file.exists(fullaggrda)){
    kept.cols = c("gene","pc", "col_nm","log10p_nm",
        "path","region",
        "logFC_nb","p_nb","padj_nb","col_nb",
        "coef_mast","p_mast","padj_mast","col_mast")

    setdflist = lapply(sets, function(x){})
    names(setdflist) = sets
    totnsigdf = NULL
    for (i in 1:nrow(selrundf)){
        prefstr = selrundf$prefstr[i]
        setid = selrundf$setid[i]
        aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
        if (file.exists(aggrda)){
            load(aggrda)  # Loads aggdf, nsig
            cat(nsig, '\n')
            # Concatenate results:
            aggdf$path = nsig[1,'path']
            aggdf$region = nsig[1,'region']
            setdflist[[setid]] = rbind(setdflist[[setid]], aggdf[,kept.cols])
            # Pad nsig:
            nsigdf = data.frame(nsig)
            if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
            if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
            totnsigdf = rbind(totnsigdf, nsigdf)
        }
    }
    save(setdflist, totnsigdf, file=fullaggrda)
} else {
    load(fullaggrda)
}

getLFC = function(df, path){
    df = df[df$path == path,]
    lfc = df$logFC_nb
    names(lfc) = df$gene
    return(lfc)
}


for (path in pathlist){
    print(path)
    fulldf = c()
    for (set in sets){
        setdf = setdflist[[set]]
        if (path %in% setdf$path){
        df = setdf[setdf$path == path, c('gene','logFC_nb', 'region')]
        df$set = set
        fulldf = rbind(fulldf, df)
        }
    }
    fulldf$rs = with(fulldf, paste0(region, ':', set))
    mat = pivot.tomatrix(fulldf[,c('gene', 'rs', 'logFC_nb')], 'rs','logFC_nb')

    # Plot heatmap:
    cmat = cor(mat, use='complete.obs')
    ux = 1.5
    plt = Heatmap(cmat, 
        name='Correlation',
        use_raster=FALSE,
        col=rev(colspec),
        border_gp=gpar(color='black', lwd=.5),
        width = ncol(cmat)*unit(ux, "mm"), 
        height = nrow(cmat)*unit(ux, "mm")
    )

    h = 2.25 + 1 / 15 * nrow(cmat)
    w = 5 + 1 / 15 * ncol(cmat)
    pltprefix = paste0(imgpref, 'correlation.exc_sets.regional_', path, '.heatmap')
    saveHeatmap(plt, pltprefix, w=w, h=h)

    submat = mat[apply(is.na(mat), 1, sum) == 0,]
    pca = princomp(submat)

    w = pca$scores
    h = pca$loadings
    for (i in 1:10){
        x = w[,i]
        x = sort(x, decreasing=T)
        cat(i, '\t', names(x)[1:20], '\n')
    }
}



# Path-level similarity:
NP = length(pathlist)
rmat = matrix(0, nrow=NP, ncol=NP, dimnames=list(pathlist, pathlist))
for (p1 in pathlist){
    l1 = getLFC(fulldf, p1)
    for (p2 in pathlist){
        l2 = getLFC(fulldf, p2)
        mn = names(l1)[names(l1) %in% names(l2)]
        rmat[p1, p2] = cor(l1[mn], l2[mn])
    }
}













    # Plot total number of significant genes as a heatmap:
    # ----------------------------------------------------
    totnsigdf$X1 = as.numeric(totnsigdf$X1)
    totnsigdf$X2 = as.numeric(totnsigdf$X2)
    umat = pivot.tomatrix(totnsigdf[,c('region','subtype','X2')], 'region','X2')
    dmat = pivot.tomatrix(totnsigdf[,c('region','subtype','X1')], 'region','X1')
    umat[is.na(umat)] = 0
    dmat[is.na(dmat)] = 0
    colnames(umat) = paste0(colnames(umat), '_up')
    colnames(dmat) = paste0(colnames(dmat), '_down')

    bmat = cbind(umat, -dmat)
    mx = max(abs(bmat))
    csplit = ifelse(1:ncol(bmat) %in% grep("_up", colnames(bmat)), 'Up','Down')

    col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
    plt = Heatmap(bmat, name='Number\nof DEGs',
        use_raster=TRUE,
        col=col_fun,
        column_split=csplit,
        border_gp=gpar(color='black', lwd=.5),
        width = ncol(bmat)*unit(4.5, "mm"), 
        height = nrow(bmat)*unit(2, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
            ann = abs(bmat[i,j])
            grid.text(ann, x, y, gp=gpar(fontsize=5)) })

    h = 2.25 + 1 / 15 * nrow(bmat)
    w = 5 + 1 / 15 * ncol(bmat)
    pltprefix = paste0(imgpref, 'allmethods_ndeg.regional_', path, '.heatmap')
    saveHeatmap(plt, pltprefix, w=w, h=h)

}

# Jaccard overlap for each set: (+ each condition):
# -------------------------------------------------
for (set in sets){
    fulldf = c()
    for (path in pathlist){
        mstr = paste0('allmethods.regional_', path)
        fullaggrda = paste0(regdir, mstr, '.merged.rda')
        load(fullaggrda)
        setdf = setdflist[[set]]
        kept.cols = c('gene','col_nm','path','region')
        fulldf = rbind(fulldf, setdf[, kept.cols])
    }

    # Make sets:
    # ----------
    fulldf = fulldf[fulldf$col_nm != 0,]
    fulldf$pr = with(fulldf, paste0(path, '@', region, '-', col_nm))
    groupnames = unique(fulldf$pr)
    genesets = lapply(groupnames, function(x){
        fulldf$gene[fulldf$pr == x] })
    names(genesets) = groupnames


    # Calculate the intersection matrix:
    # ----------------------------------
    N = length(genesets)
    jmat = matrix(0, nr=N, nc=N)
    for (i in 1:N){
        s1 = genesets[[i]]
        for (j in 1:N){
            s2 = genesets[[j]]
            jmat[i,j] = length(intersect(s1, s2)) / length(union(s1, s2))
        }
    }
    rownames(jmat) = names(genesets)
    colnames(jmat) = names(genesets)

    rmat = reord(jmat)
    rn = rownames(rmat)
    jmat = jmat[rn,rn]

    # Add an annotation based on the names:
    # -------------------------------------
    ux = 1.5
    rann = sub("-.*","", sub(".*@", "", rownames(jmat)))
    pann = sub("@.*", "", rownames(jmat))
    dann = sub(".*-", "", rownames(jmat))
    path.cols = snap.cols[1:length(pathlist)]
    names(path.cols) = rev(pathlist)
    ha = HeatmapAnnotation(
        Region=rann, Path=pann, Dir=dann,
        annotation_name_gp = gpar(fontsize=5),
        simple_anno_size = unit(ux, 'mm'),
        gap = unit(0, "mm"),
        col=list(
            Region=reg.cols,
            Path=path.cols,
            Dir=c('1'='blue','2'='red')
            ), 
        which='row')


    plt = Heatmap(jmat, 
        name='DEG overlap (jaccard)',
        use_raster=FALSE,
        col=rev(colspec),
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        column_split=dann,
        row_split=dann,
        right_annotation=ha,
        border_gp=gpar(color='black', lwd=.5),
        width = ncol(jmat)*unit(ux, "mm"), 
        height = nrow(jmat)*unit(ux, "mm")
    )

    h = 2.25 + 1 / 15 * nrow(jmat)
    w = 5 + 1 / 15 * ncol(jmat)
    pltprefix = paste0(imgpref, 'allmethods_jaccard_ovl.set_', set, '.heatmap')
    saveHeatmap(plt, pltprefix, w=w, h=h)
}

# Top changing and consistent DEGs
