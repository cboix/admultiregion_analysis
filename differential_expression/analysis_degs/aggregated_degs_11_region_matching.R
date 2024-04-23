#!/usr/bin/R
# -------------------------------------------------------
# Region-region stage matching
# Updated: 11/010/23
# ------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(MASS)

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


# Functions:
# ----------
load_psbulk_matrix = function(set){
    ststr = sub('_[A-Za-z]+$', '', set)
    psdata.rda = paste0(srdir, 'pseudobulk_data_', ststr, '.rda')
    load(psdata.rda)
    # If Mic, remove T, CAMs from psbulk data:
    if (set == 'Mic_Immune_Mic'){
        ind = !(ps.data$meta$cell_type_high_resolution %in% c('T cells', 'CAMs'))
        ps.data$meta = ps.data$meta[ind,]
        ps.data$mat = ps.data$mat[,rownames(ps.data$meta)]
    }
    # Merge to sample-level (individual x region):
    ps.data = aggregate_psbulk_samplelevel(ps.data)
    return(ps.data)
}

pivot.tomatrix.v2 = function(df, index, key, value){
    require(dplyr)
    df = df[,c(index, key, value)]
    wide = spread(df, key, value)
    mat = as.matrix(wide[, -1])
    rownames(mat) = wide[, 1]
    return(mat)
}


# Arguments for runs:
# -------------------
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", "Oli_Oli", 'Inh_Inh','Exc_Exc')
setnames = c("Mic_Immune_Mic"='Microglial', "Ast_Ast"='Astrocyte', 
    "Opc_Opc"='OPC', "Oli_Oli"='Oligodendrocyte', 
    'Inh_Inh'='Inhibitory neuron','Exc_Exc'='Excitatory neuron')
pathlist = c('plaq_d', 'plaq_n', 'nft', 'nrad', 'cogdxad')
regset = c('allregions', 'EC', 'HC', 'TH', 'AG', 'MT', 'PFC')


# Load ABC scores and add metadata:
# ---------------------------------
indmeta_tsv = 'Annotation/metadata_PFC_all_individuals_092520.tsv'
ext.meta = read.delim(indmeta_tsv, header=T)
ext.meta$kept.ind = ifelse(ext.meta$projid %in% kept.individuals, 'Our Cohort\n(48 individ.)', 'ROSMAP\n(Sept. 2020)')
abcdf = read.delim('Annotation/abc_scores_092023.tsv', header=T)
abcdf = merge(abcdf, ext.meta[,c('projid', 'niareagansc', 'cogdx')])

# Add metadata variables:
metadata$braak.group = ifelse(metadata$braaksc %in% c(0,1,2), '0-2', ifelse(metadata$braaksc %in% c(3,4), '3-4', '5-6'))
meta.cols = c('projid', 'region', 'rind', 'braaksc', 'nrad','niareagansc','cogdx', 'Apoe_e4', 'braak.group')
metadf = merge(abcdf, unique(metadata[,meta.cols]))


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


# DEGs in at least 3 calculations, low bar:
aggdf = aggregate(path ~ gene + col_nm, alldf, length)

min.ndeg = 3
degs = unique(aggdf$gene[(aggdf$col_nm != 0) & (aggdf$path >= min.ndeg)])
updegs = unique(aggdf$gene[(aggdf$col_nm == 2) & (aggdf$path >= min.ndeg)])

deglist = list()
min.ndeg.ct = 10
for (set in keep.sets){
    aggdf = aggregate(path ~ gene + col_nm, alldf[alldf$set == set,], length)
    deglist[[set]] = sort(unique(aggdf$gene[(aggdf$col_nm != 0) & (aggdf$path >= min.ndeg.ct)]))
}
sapply(deglist, length)

nr.deglist = list()
min.ndeg.nr = 3
for (set in keep.sets){
    aggdf = aggregate(path ~ gene + col_nm, alldf[(alldf$path == 'nrad') & (alldf$set == set),], length)
    nr.deglist[[set]] = sort(unique(aggdf$gene[(aggdf$col_nm != 0) & (aggdf$path >= min.ndeg.nr)]))
}
sapply(nr.deglist, length)



for (set in keep.sets){
    print(set)
    # Load pseudobulk data:
    ps.data = load_psbulk_matrix(set)

    psdf = merge(ps.data$meta, metadf, all.x=TRUE)
    # y = ps.data$mat[degs, psdf$pr]
    y = ps.data$mat[, psdf$pr]

    df = as.data.frame(as.matrix(y), check.names=F)
    df$gene = rownames(df)
    df = gather(df, pr, value, -gene)
    df = merge(df, psdf)

    run.regr = FALSE
    if (run.regr){
        # Regress out region + baseline:
        # TODO: run per gene
        df$pred = 0
        chunksize = 25
        NG = length(degs)
        nchunk = ceiling(NG / chunksize)
        for (i in 1:nchunk){
            print(i)
            genes = degs[((i-1) * chunksize + 1):min(c(NG, i * chunksize))]
            ind = df$gene %in% genes
            fit = glm(value ~ gene * region, data=df[ind,], family='gaussian')
            df$pred[ind] = predict(fit, df[ind,])
        }

        df$resid = df$value - df$pred
        cor(df$value, df$pred) # R=0.95
        cor(df$value, df$resid) # R=0.31
    }

    # NOTE: EVAL FROM RESID OR FROM REGION "BASELINE" e.g. 0?? -- baseline, resid doesn't work.
    # Match regions:
    var = 'nia_aa_sc'
    pred.var = 'value'
    runsuffix = paste0(var, '.', pred.var, '.', set)

    formvec = c(pred.var, '~ ', var, '+ region + gene')
    form = asform(formvec)
    aggdf = agg.rename(form, df, function(x){
            c(mean(x), sd(x), length(x))}, 'out')
    stats = data.frame(aggdf$out)
    names(stats) = c('mean','sd','n')
    aggdf = cbind(aggdf, stats)
    aggdf$out = NULL

    aggdf$group = paste0(aggdf$region, '@', aggdf[[var]])
    mat = pivot.tomatrix.v2(aggdf, 'gene', 'group', 'mean')

    # geneset = deglist[[set]]
    geneset = nr.deglist[[set]]
    geneset = geneset[geneset %in% rownames(mat)]
    centmat = sweep(mat[geneset,], 1, apply(mat[geneset,],1, mean), '-')
    cr = cor(centmat)
    dt = dist(t(centmat))
    dt = as.matrix(dt)

    for (type in c('distance','correlation')){
        if (type == 'correlation'){
            pltmat = cr
            mx = 1
            col_fun = colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
        } else {
            pltmat = dt
            col_fun = colspec
        } 

        colsplit = sub("@.*","", colnames(pltmat))
        rowsplit = sub("@.*","", rownames(pltmat))
        ux = 1.5
        ht = Heatmap(pltmat, 
            use_raster=FALSE, 
            name=capitalize(type),
            col=col_fun,
            cluster_columns=FALSE, 
            cluster_rows=FALSE,
            cluster_row_slices=FALSE,
            cluster_column_slices=FALSE,
            column_split=colsplit,
            row_split=rowsplit, 
            row_dend_width = unit(.25, "cm"),
            column_dend_height = unit(.25, "cm"),
            row_dend_gp = gpar(lwd=.5),
            column_dend_gp = gpar(lwd=.5),
            border_gp = gpar(col="black", lwd=.5),
            width=ncol(pltmat) * unit(ux, 'mm'),
            height=nrow(pltmat) * unit(ux, 'mm')
        )
        pltprefix = paste0(imgpref, 'region_comparison_heatmap.', type, '.degs_', pred.var, '.', set)
        h = 1.5 + 1 / 15 * nrow(pltmat)
        w = 2 + 1 / 15 * ncol(pltmat)
        saveHeatmap(ht, pltprefix, w=w, h=h)
    }

    # Average score:
    ddf = data.frame(dt, check.names=F)
    ddf$G1 = rownames(ddf)
    ddf = gather(ddf, G2, distance, -G1)
    ddf$R1 = sub('@.*', '', ddf$G1)
    ddf$R2 = sub('@.*', '', ddf$G2)
    ddf$V1 = sub('.*@', '', ddf$G1)
    ddf$V2 = sub('.*@', '', ddf$G2)
    # TODO: Map V1, V2 to levels if necessary
    ddf$V1 = as.numeric(ddf$V1)
    ddf$V2 = as.numeric(ddf$V2)

    estdf = c()
    for (R1 in reg.nomb){
        for (lvl in unique(ddf$V1)){
            for (R2 in reg.nomb){
                sub.ddf = ddf[(ddf$V1 == lvl) & (ddf$R1 == R1) & (ddf$R2 == R2),]
                sim = max(sub.ddf$distance) - sub.ddf$distance 
                est = sum(sim * sub.ddf$V2) / sum(sim)
                estdf = rbind(estdf, data.frame(R1=R1, R2=R2, V1=lvl, Est=est))
            }
        }
    }

    gp = ggplot(resdf, aes(V1, Est, color=R1)) + 
        geom_point(cex=.75) + 
        scale_color_manual(values=reg.cols) + 
        geom_smooth(method='lm', color='black') + 
        stat_cor(color='black') + 
        geom_abline(intercept=0, slope=1, lty='dashed') + 
        theme_pubr()
    pltprefix = paste0(imgpref, 'region_matching_scatter.distance.degs_', pred.var, '.', set)
    saveGGplot(gp, pltprefix, w=4, h=4)


    # Estimate by bootstrap:
    # ----------------------
    boot.pref = paste0('bootstraps_region_mapping.', runsuffix)
    boot.ct.rds = paste0(srdir, boot.pref, '.Rds')
    boot.ct.tsv = paste0(srdir, boot.pref, '.tsv.gz')
    if (!file.exists(boot.ct.rds)){
        # For runs:
        formvec = c(pred.var, '~ ', var, '+ region + gene')
        form = asform(formvec)
        geneset = deglist[[set]]
        geneset = geneset[geneset %in% df$gene]

        # Estimate from full sample:
        # --------------------------
        aggdf = agg.rename(form, df, mean, 'mean')
        aggdf$group = paste0(aggdf$region, '@', aggdf[[var]])

        # Get distance matrix:
        mat = pivot.tomatrix.v2(aggdf, 'gene', 'group', 'mean')
        centmat = sweep(mat[geneset,], 1, apply(mat[geneset,],1, mean), '-')
        # dt = dist(t(mat[geneset,]))
        dt = dist(t(centmat))
        dt = as.matrix(dt)

        # Average score:
        ddf = data.frame(dt, check.names=F)
        ddf$G1 = rownames(ddf)
        ddf = gather(ddf, G2, distance, -G1)
        ddf$R1 = sub('@.*', '', ddf$G1)
        ddf$R2 = sub('@.*', '', ddf$G2)
        ddf$V1 = sub('.*@', '', ddf$G1)
        ddf$V2 = sub('.*@', '', ddf$G2)
        # TODO: Map V1, V2 to levels if necessary
        ddf$V1 = as.numeric(ddf$V1)
        ddf$V2 = as.numeric(ddf$V2)

        # Center of mass:
        estdf = c()
        for (R1 in reg.nomb){
            for (lvl in unique(ddf$V1)){
                for (R2 in reg.nomb){
                    sub.ddf = ddf[(ddf$V1 == lvl) & (ddf$R1 == R1) & (ddf$R2 == R2),]
                    sim = max(sub.ddf$distance) - sub.ddf$distance 
                    est = sum(sim * sub.ddf$V2) / sum(sim)
                    estdf = rbind(estdf, data.frame(R1=R1, R2=R2, V1=lvl, Est=est))
                }
            }
        }


        # Run bootstraps
        # --------------
        NBOOT = 100
        samplesize = 36
        bootdf = c()
        for (k in 1:NBOOT){
            cat('Bootstrap', k, '\n')
            set.seed(k)
            boot.individuals = sample(kept.individuals, size=samplesize, replace=FALSE)
            aggdf = agg.rename(form, df[df$projid %in% boot.individuals,],
                mean, 'mean')
            aggdf$group = paste0(aggdf$region, '@', aggdf[[var]])

            # Get distance matrix:
            mat = pivot.tomatrix.v2(aggdf, 'gene', 'group', 'mean')
            centmat = sweep(mat[geneset,], 1, apply(mat[geneset,],1, mean), '-')
            dt = dist(t(mat[geneset,]))
            dt = as.matrix(dt)

            # Average score:
            ddf = data.frame(dt, check.names=F)
            ddf$G1 = rownames(ddf)
            ddf = gather(ddf, G2, distance, -G1)
            ddf$R1 = sub('@.*', '', ddf$G1)
            ddf$R2 = sub('@.*', '', ddf$G2)
            ddf$V1 = sub('.*@', '', ddf$G1)
            ddf$V2 = sub('.*@', '', ddf$G2)
            # TODO: Map V1, V2 to levels if necessary
            ddf$V1 = as.numeric(ddf$V1)
            ddf$V2 = as.numeric(ddf$V2)

            # Center of mass:
            resdf = c()
            for (R1 in reg.nomb){
                for (lvl in unique(ddf$V1)){
                    for (R2 in reg.nomb){
                        sub.ddf = ddf[(ddf$V1 == lvl) & (ddf$R1 == R1) & (ddf$R2 == R2),]
                        sim = max(sub.ddf$distance) - sub.ddf$distance 
                        est = sum(sim * sub.ddf$V2) / sum(sim)
                        resdf = rbind(resdf, data.frame(R1=R1, R2=R2, V1=lvl, Est=est))
                    }
                }
            }
            resdf$bootstrap = k
            bootdf = rbind(bootdf, resdf)
        }
        estdf$bootstrap = 'Original'
        fulldf = rbind(estdf, bootdf)

        # Save bootstraps estimates:
        saveRDS(fulldf, file=boot.ct.rds)
        write.table(fulldf, file=gzfile(boot.ct.tsv), quote=F, row.names=F, sep="\t")
    } else {
        fulldf = readRDS(file=boot.ct.rds)
        estdf = fulldf[fulldf$bootstrap == 'Original',]
        bootdf = fulldf[fulldf$bootstrap != 'Original',]
    }


    # Plot estimates for Empirical and normal bootstrap:
    # --------------------------------------------------
    estdf.v2 = estdf[,c('R1','R2','V1','Est')]
    names(estdf.v2)[4] = c('full.Est')
    bootdf = merge(bootdf, estdf.v2)
    bootdf$delta = bootdf$Est - bootdf$full.Est

    # 90% CI 
    statdf = aggregate(cbind(delta, Est) ~ R1 + R2 + V1 + full.Est, 
        bootdf, function(x){ quantile(x, c(.05, .95)) })
    cidf = cbind(data.frame(statdf$delta), data.frame(statdf$Est))
    names(cidf) = c('delta.lo', 'delta.hi', 'est.lo', 'est.hi')
    statdf = cbind(statdf, cidf)
    statdf = merge(statdf, agg.rename(Est ~ R1 + R2 + V1 + full.Est, bootdf, sd, 'est.sd'))
    statdf$est.se = statdf$est.sd / sqrt(samplesize)

    gp = ggplot(statdf, aes(full.Est, V1, xmin=full.Est + delta.lo, xmax=full.Est + delta.hi, color=R1)) + 
        facet_grid(R1 ~ R2) + 
        geom_point() + 
        scale_color_manual(values=reg.cols) + 
        geom_errorbar(width=.5) + 
        geom_abline(intercept=0, slope=1, lty='dashed') + 
        labs(x='Estimated stage in second region', y='Actual stage in first region') + 
        theme_bw()
    pltprefix = paste0(imgpref, 'region_matching_scatter_bootdelta.distance.', runsuffix)
    saveGGplot(gp, pltprefix, w=8, h=5)

    zcrit = 1.96 # 95%
    # zcrit = 1.645 # 90
    gp = ggplot(statdf, aes(full.Est, V1, xmin=full.Est - zcrit * est.sd, xmax=full.Est + zcrit * est.sd, color=R1)) + 
        facet_grid(R1 ~ R2) + 
        geom_point() + 
        scale_color_manual(values=reg.cols) + 
        geom_errorbar(width=.5) + 
        geom_abline(intercept=0, slope=1, lty='dashed') + 
        labs(x='Estimated stage in second region', y='Actual stage in first region') + 
        theme_bw()
    pltprefix = paste0(imgpref, 'region_matching_scatter_bootnormal.distance.', runsuffix)
    saveGGplot(gp, pltprefix, w=8, h=5)



    # MDS from full sample:
    # ---------------------
    aggdf = agg.rename(form, df, mean, 'mean')
    aggdf$group = paste0(aggdf$region, '@', aggdf[[var]])
    mat = pivot.tomatrix.v2(aggdf, 'gene', 'group', 'mean')
    # ind = grep("TH", colnames(mat), invert=TRUE)
    # dt = dist(t(mat[geneset, ind]))
    dt = dist(t(mat[geneset, ]))

    # TODO: k=1, 2:
    k = 1
    for (k in c(1,2)){
        fit = cmdscale(dt, eig=TRUE, k=k) # Classical
        # fit = isoMDS(dt, k=k) # Nonmetric (MASS)
        mdf = data.frame(fit$points)
        names(mdf) = paste0('M', 1:k)
        mdf$group = rownames(mdf)
        mdf$region = sub("@.*","", mdf$group)
        mdf$var = sub(".*@","", mdf$group)
        if (k==1){
            gp = ggplot(mdf, aes(x=M1, y=region, color=region, label=var)) 
        } else {
            gp = ggplot(mdf, aes(M1, M2, color=region, label=var))
        }
        gp = gp + geom_point() + 
            geom_text_repel() + 
            scale_color_manual(values=reg.cols) + 
            theme_pubr() + theme(legend.position='FALSE')
        pltprefix = paste0(imgpref, 'region_matching.mds_k', k, '.', runsuffix)
        saveGGplot(gp, pltprefix, w=3, h=1 + k)
    }

}


# TODO: Can we also do this for other variables and NFT, plaque load?
# TODO: Can we build full continuum from this --> e.g. time course


