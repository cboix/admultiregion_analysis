#!/usr/bin/R
# ----------------------------------------------------------
# Plot fractions + fractions on AD
# Updated 09/22/2020 
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

# --------------------------------------
# Load in the final metadata (cellmeta):
# --------------------------------------
load(file=paste0(datadir, prefix, '.final_noMB.cell_labels.Rda'))

# Colors for full:
typelvls = unique(cellmeta$cell_type_high_resolution)
type.cols = rep(snap.cols,3)[1:length(typelvls)]
names(type.cols) = as.character(typelvls)
type.cols = c(type.cols, major.col['Inh'], major.col['Exc'])
tsp.type.cols = sapply(type.cols, tsp.col)

# -----------------------
# Plot of overall counts:
# -----------------------
ctdf = aggregate(barcode ~ major.celltype + minor.celltype, cellmeta, length)
g1 = ggplot(ctdf, aes(major.celltype, barcode, fill=minor.celltype)) + 
    geom_bar(stat='identity') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    scale_fill_manual(values=type.cols, name='Minor\nCell-types:') + 
    labs(x='Major Cell-type', y='Number of cells') + 
    theme_pubr()
ggsave(paste0(imgpref, 'metadata_overall_fractions.png'), g1, dpi=450, units='in', width=6.5, height=4.5)


ctdf = aggregate(barcode ~ major.celltype + cell_type_high_resolution, cellmeta, length)
g2 = ggplot(ctdf, aes(major.celltype, barcode, fill=cell_type_high_resolution)) + 
    geom_bar(stat='identity') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    scale_fill_manual(values=type.cols, name='Full\nTypes:') + 
    labs(x='Major Cell-type', y='Number of cells') + 
    theme_pubr()
ggsave(paste0(imgpref, 'metadata_overall_fractions_fulltypes.png'), g2, dpi=450, units='in', width=11.25, height=11)


# Order for full.exttype:
odf = unique(cellmeta[,c('major.celltype','full.exttype')])
odf = odf[order(odf$full.exttype),]
odf$major.celltype = factor(odf$major.celltype, levels=c('Ast','Mic/Immune','Oli','Opc','Vasc/Epithelia', 'Exc','Inh'))
odf = odf[order(odf$major.celltype),]
odf$full.exttype = factor(odf$full.exttype, levels=odf$full.exttype)

# Order for cell_type_high_resolution:
codf = unique(cellmeta[,c('major.celltype','cell_type_high_resolution')])
codf = codf[order(codf$cell_type_high_resolution),]
codf$major.celltype = factor(codf$major.celltype, levels=c('Ast','Mic/Immune','Oli','Opc','Vasc/Epithelia', 'Exc','Inh'))
codf = codf[order(codf$major.celltype),]
codf$cell_type_high_resolution = factor(codf$cell_type_high_resolution, levels=codf$cell_type_high_resolution)


# -------------------------------------
# Make the per-region pathology scores:
# -------------------------------------
# NFT:
metadata$nft_regional = NA
metadata$nft_regional[metadata$region == 'AG'] = metadata$nft_ag[metadata$region == 'AG']
metadata$nft_regional[metadata$region == 'PFC'] = metadata$nft_mf[metadata$region == 'PFC']
metadata$nft_regional[metadata$region == 'MT'] = metadata$nft_mt[metadata$region == 'MT']
metadata$nft_regional[metadata$region == 'HC'] = metadata$nft_hip[metadata$region == 'HC']
metadata$nft_regional[metadata$region == 'EC'] = metadata$nft_ec[metadata$region == 'EC']
# Plaq_d
metadata$plaq_d_regional = NA
metadata$plaq_d_regional[metadata$region == 'AG'] = metadata$plaq_d_ag[metadata$region == 'AG']
metadata$plaq_d_regional[metadata$region == 'PFC'] = metadata$plaq_d_mf[metadata$region == 'PFC']
metadata$plaq_d_regional[metadata$region == 'MT'] = metadata$plaq_d_mt[metadata$region == 'MT']
metadata$plaq_d_regional[metadata$region == 'HC'] = metadata$plaq_d_hip[metadata$region == 'HC']
metadata$plaq_d_regional[metadata$region == 'EC'] = metadata$plaq_d_ec[metadata$region == 'EC']
# Plaq_d
metadata$plaq_n_regional = NA
metadata$plaq_n_regional[metadata$region == 'AG'] = metadata$plaq_n_ag[metadata$region == 'AG']
metadata$plaq_n_regional[metadata$region == 'PFC'] = metadata$plaq_n_mf[metadata$region == 'PFC']
metadata$plaq_n_regional[metadata$region == 'MT'] = metadata$plaq_n_mt[metadata$region == 'MT']
metadata$plaq_n_regional[metadata$region == 'HC'] = metadata$plaq_n_hip[metadata$region == 'HC']
metadata$plaq_n_regional[metadata$region == 'EC'] = metadata$plaq_n_ec[metadata$region == 'EC']

# By region:
ctdf = agg.rename(barcode ~ major.celltype + cell_type_high_resolution + region, cellmeta, length, 'ncell')
g2 = ggplot(ctdf, aes(major.celltype, ncell, fill=cell_type_high_resolution)) + 
    facet_wrap(~region) + 
    geom_bar(stat='identity') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    scale_fill_manual(values=type.cols, name='Full\nTypes:') + 
    labs(x='Major Cell-type', y='Number of cells') + 
    theme_pubr()
ggsave(paste0(imgpref, 'metadata_overall_fractions_fulltypes_byregion.png'), g2, dpi=450, units='in', width=11.25, height=11)

# ectdf = ctdf[ctdf$major.celltype == 'Exc' & ctdf$region == 'PFC',]
# ectdf = ectdf[order(ectdf$ncell, decreasing=T),]

# ectdf = ctdf[ctdf$major.celltype == 'Inh' & ctdf$region == 'PFC',]
# ectdf = ectdf[order(ectdf$ncell, decreasing=T),]


# Fractions vs. region at the basic level 
vars = c('major.celltype', 'minor.celltype','full.exttype', 'cell_type_high_resolution')
for (ct.var in vars) {
    totdf = agg.rename(barcode ~ rind + region, cellmeta, length, 'totcell')
    ctdf = agg.rename(as.formula(paste0('barcode ~ ', ct.var, '+ rind + region')), cellmeta, length, 'ncell')
    ctdf = merge(ctdf, totdf)
    ctdf$frac = ctdf$ncell / ctdf$totcell
    ctdf = merge(ctdf, metadata[,c('rind','plaq_n','plaq_d','nft','cogdx','braaksc','niareagansc', 'msex','age_death','pmi', 'projid')])


    ctdf$projid = factor(ctdf$projid)
    ctdf$nrad = ctdf$niareagansc > 2
    ctdf$cogad = ctdf$cogdx > 3
    ctdf$ltot = log10(ctdf$totcell)
    ctdf$frac = ctdf$ncell / ctdf$totcell

    # Example of interaction term:
    cell = 'Exc NRGN'
    l1 = glm(frac ~ nrad + region + msex + age_death + pmi, ctdf[ctdf$cell_type_high_resolution == cell,], family=quasibinomial)
    l2 = glm(frac ~ nrad * region + msex + age_death + pmi, ctdf[ctdf$cell_type_high_resolution == cell,], family=quasibinomial)
    summary(l2)
    anova(l1, l2)


    coefdf = c()
    for (cell in unique(ctdf$cell_type_high_resolution)){
        print(cell)
        fit = glm(frac ~ nft + region + msex + age_death + pmi, ctdf[ctdf$cell_type_high_resolution == cell,], family=quasibinomial)
        # fit = glm(frac ~ cogad + region + msex + age_death + pmi, ctdf[ctdf$cell_type_high_resolution == cell,], family=quasibinomial)
        # fit = glm(frac ~ nrad + region + msex + age_death + pmi, ctdf[ctdf$cell_type_high_resolution == cell,], family=quasibinomial)
        cdf = data.frame(coefficients(summary(fit)))
        names(cdf) = c('est','se','t','p')
        cdf$var = rownames(cdf)
        cdf$ct = cell
        coefdf = rbind(coefdf, cdf)
    }
    coefdf = coefdf[order(coefdf$p),]
    # cdf = coefdf[coefdf$var == 'nradTRUE',]
    cdf = coefdf[coefdf$var == 'nft',]
    head(cdf,10)
    cdf$q = qvalue(cdf$p)$q

    # Plot these top ones:
    topct = cdf$ct[cdf$q < 0.05]
    scale = 1
    w = 2 + length(topct) / scale
    h = 2 + 6 / scale *2
    gs = ggplot(ctdf[ctdf$cell_type_high_resolution %in% topct,], aes(nft, frac, fill=region)) +
        facet_grid(cell_type_high_resolution~region, scale='free_y') + 
        geom_smooth(method='lm', color='black') + 
        scale_color_manual(values=reg.cols) + 
        scale_y_continuous(labels=scales::percent) + 
        scale_fill_manual(values=reg.cols) + 
        labs(x='NFT density (z-score aggregate)', y='Percent of cells') + 
        geom_point() + 
        theme_pubr() +
        theme(legend.position='none')
    ggsave(paste0(imgpref, 'metadata_AD_fractions_nft_overall_vs_', ct.var ,'.png'), gs, dpi=450, units='in', width=w, height=h)

    # Plot ec specific in EC/HC
    ecneu = c('Exc AGBL1 GPC5', 'Exc RELN GPC5', 'Exc TRPC6 ANO2', 'Exc COBLL1 UST', 'Exc RELN COL5A2',
              'Exc COL25A1 SEMA3D', 'Exc TOX3 POSTN', 'Exc TOX3 INO80D', 'Exc TOX3 TTC6')
    subctdf = ctdf[ctdf$cell_type_high_resolution %in% ecneu & ctdf$region %in% c('EC','HC'),]
    w = 2 + length(ecneu) / scale / 1.5
    h = 2 + 2 / scale * 2
    gs = ggplot(subctdf, aes(nft, frac, fill=region, color=region)) +
        # facet_grid(cell_type_high_resolution~region, scale='free_y') + 
        facet_wrap(cell_type_high_resolution~., scale='free_y') + 
        geom_smooth(method='lm', color='black') + 
        scale_color_manual(values=reg.cols) + 
        scale_y_continuous(labels=scales::percent) + 
        scale_fill_manual(values=reg.cols) + 
        labs(x='NFT density (z-score aggregate)', y='Percent of cells') + 
        geom_point() + 
        theme_pubr() +
        theme(legend.position='none')
    ggsave(paste0(imgpref, 'metadata_AD_fractions_ecneurons_nft_overall_vs_', ct.var ,'.png'), gs, dpi=450, units='in', width=w, height=h)
    ggsave(paste0(imgpref, 'metadata_AD_fractions_ecneurons_nft_overall_vs_', ct.var ,'.pdf'), gs, width=w, height=h)


    # Plot inhibitory (SST cluster) overall:
    sstcluster = c('Inh L6 SST NPY', 'Inh CUX2 MSR1', 'Inh L3-5 SST MAFB', 'Inh ENOX2 SPHKAP',
                   'Inh FBN2 EPB41L4A', 'Inh GPC5 RIT2', 'Inh PVALB SULF1', 'Inh PVALB HTR4')
    subctdf = ctdf[ctdf$cell_type_high_resolution %in% sstcluster,]
    w = 2 + length(sstcluster) / scale / 1.5
    h = 2 + 2 / scale * 2
    gs = ggplot(subctdf, aes(nft, frac, fill=region, color=region)) +
        # facet_grid(cell_type_high_resolution~region, scale='free_y') + 
        facet_wrap(cell_type_high_resolution~., scale='free_y') + 
        geom_smooth(method='lm', color='black') + 
        scale_color_manual(values=reg.cols) + 
        scale_y_continuous(labels=scales::percent) + 
        scale_fill_manual(values=reg.cols) + 
        labs(x='NFT density (z-score aggregate)', y='Percent of cells') + 
        geom_point() + 
        theme_pubr() +
        theme(legend.position='none')
    ggsave(paste0(imgpref, 'metadata_AD_fractions_inhsstcluster_nft_overall_vs_', ct.var ,'.png'), gs, dpi=450, units='in', width=w, height=h)
    ggsave(paste0(imgpref, 'metadata_AD_fractions_inhsstcluster_nft_overall_vs_', ct.var ,'.pdf'), gs, width=w, height=h)


}


# ----------------------------------
# Fractions against region:
# Both overall and per-region-scores
# ----------------------------------
vars = c('major.celltype', 'minor.celltype','full.exttype', 'cell_type_high_resolution')
for (ct.var in vars) {
    totdf = agg.rename(barcode ~ rind + region, cellmeta, length, 'totcell')
    ctdf = agg.rename(as.formula(paste0('barcode ~ ', ct.var, '+ rind + region')), cellmeta, length, 'ncell')
    ctdf = merge(ctdf, totdf)
    ctdf$frac = ctdf$ncell / ctdf$totcell
    ctdf = merge(ctdf, metadata[,c('rind','plaq_n','plaq_d','nft','cogdx','braaksc','niareagansc', 'msex','age_death','pmi')])

    ctdf$nrad = ctdf$niareagansc > 2
    ctdf$cogad = ctdf$cogdx > 3
    ctdf$ltot = log10(ctdf$totcell)
    ctdf$frac = ctdf$ncell / ctdf$totcell

    # Fraction vs. region regression:
    ctdf$logit.frac = logit(ctdf$frac)
    fit = lm(as.formula(paste0('logit.frac ~ ', ct.var, '* region + 1')), ctdf)
    summary(fit)

    # Logit frac vs. region + reponse:
    resp.vars = c('plaq_n','plaq_d','nft','cogdx','braaksc')
    m1 = lm(as.formula(paste0('logit.frac ~ ', ct.var, '* region + 1')), ctdf)
    vdf = c()
    for (var in resp.vars){
        m2 = lm(as.formula(paste('logit.frac ~ ', var, ' * ', ct.var, ' + ', ct.var, '* region + 1')), ctdf)
        pval = anova(m1, m2)$P[2]
        cat(var, '\t', pval,'\n')
        cf = coefficients(summary(m2))
        print(cf[grep(var, rownames(cf)),])
        sdf = as.data.frame(cf[grep(var, rownames(cf)),c(1,4)])
        sdf$var = var
        sdf$rn = sub(var, '', rownames(sdf))
        sdf$rn[sdf$rn == ''] = 'var'
        sdf$log10p = -log10(sdf$'Pr(>|t|)')
        vdf = rbind(vdf, sdf)
        cat('\n')
    }
    rownames(vdf) = NULL
    vdf$rn =sub(paste0(":", ct.var), '', vdf$rn)
    vdf$rn[vdf$rn == 'var'] = 'Intercept'
    if (ct.var == 'full.exttype'){ 
        vdf$rn = factor(vdf$rn, levels=c('Intercept', rev(as.character(odf$full.exttype))))
    }
    pwide = spread(vdf[,c('var','rn','log10p')], rn, log10p)
    rownames(pwide) = pwide[,1]
    rwide = spread(vdf[,c('var','rn','Estimate')], rn, Estimate)
    # TODO: Reorder.
    # Factor for plotting:
    vdf$var = factor(vdf$var)
    vdf$rn = factor(vdf$rn)

    w = 3 + dim(rwide)[1] * 3 / 5
    h = 1.5 + dim(rwide)[2] / 5
    gp = ggplot(vdf, aes(var, rn, size=log10p, color=Estimate)) + 
        # scale_color_continuous(colryb) + 
        geom_point() + 
        theme_pubr() + 
        scale_color_distiller(palette='RdBu',name='Est:') + 
        labs(x='Response Variable', y='Cell type') + 
        theme(axis.text.x=element_text(angle=45, hjust=1)) + 
        theme(legend.position='right')

    ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'.png'), gp, dpi=450, units='in', width=w, height=h)
    ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'.pdf'), gp, width=w, height=h)


    # TODO: automate:
    # Plot some of these top ones:
    if (ct.var == 'full.exttype'){
        keep.vals = c('Oli','NRGN','Granule','CPEC','CAM', 'CA1_Pyr', 'SST','PVALB')
    } else if (ct.var == 'major.celltype'){
        keep.vals = unique(ctdf[[ct.var]]) 
    }
    subdf = ctdf[ctdf[[ct.var]] %in% keep.vals,]

    scale = 4
    h = 2 * scale 
    w = length(keep.vals) / 2 * scale
    # gs = ggplot(subdf, aes(nft, logit.frac, color=region)) + 
    gs = ggplot(subdf, aes(nft, logit.frac, color=region)) + 
        facet_wrap(as.formula(paste0('~', ct.var)), nrow=2) + 
        geom_point() + geom_smooth(method = 'lm') + 
        scale_color_manual(values=reg.cols) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_examples_nft.png'), gs, dpi=450, units='in', width=w, height=h)

    gs = ggplot(subdf, aes(plaq_d, logit.frac, color=region)) + 
        facet_wrap(as.formula(paste0('~', ct.var)), nrow=2) + 
        geom_point() + geom_smooth(method = 'lm') + 
        scale_color_manual(values=reg.cols) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_examples_plaq_d.png'), gs, dpi=450, units='in', width=w, height=h)


    slong = subdf[,c('frac','nft','plaq_d','plaq_n', ct.var,'region', 'rind')]
    swide = gather(slong, var, value, -frac, -ct.var, -region, -rind)

    for (val in keep.vals){
        print(val)
        valstr = sub("/",".", val)
        gs = ggplot(subdf[subdf[[ct.var]] == val,], aes(nft, frac, fill=region)) + 
            facet_wrap(~region, nrow=1) + 
            geom_point() + geom_smooth(method = 'lm') + 
            labs(x='PHFTau Density (nft)', y=paste0(val, ' % of total cells')) + 
            scale_y_continuous(labels=scales::percent) + 
            scale_fill_manual(values=reg.cols) + 
            theme_pubr()
        gs2 = gs + labs(x='PHFTau Density (nft)', y=paste0(val, ' % of total cells\n(log10 scale)')) + 
            scale_y_log10(labels=scales::percent) + 
        ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_nft.png'), gs, dpi=450, units='in', width=w, height=5)
        ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_nft_log10.png'), gs2, dpi=450, units='in', width=w, height=5)

        gs = ggplot(swide[subdf[[ct.var]] == val,], aes(value, frac, fill=region)) + 
            facet_grid(var~region) + 
            geom_point() + geom_smooth(method = 'lm') + 
            labs(x='Pathology Score', y=paste0(val, ' % of total cells')) + 
            scale_y_continuous(labels=scales::percent) + 
            scale_fill_manual(values=reg.cols) + 
            theme_pubr()
        # gs2 = gs + labs(x='Pathology Score', y=paste0(val, ' % of total cells\n(log10 scale)')) + 
        #     scale_y_log10(labels=scales::percent) + 
        ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_3path.png'), gs, dpi=450, units='in', width=w, height=7)
        # ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_nft_log10.png'), gs2, dpi=450, units='in', width=w, height=5)
    }


    # --------------------------------------------
    # Re-run with the per-region pathology scores:
    # --------------------------------------------
    totdf = agg.rename(barcode ~ rind + region, cellmeta, length, 'totcell')
    ctdf = agg.rename(as.formula(paste0('barcode ~ ', ct.var, '+ rind + region')), cellmeta, length, 'ncell')
    ctdf = merge(ctdf, totdf)
    ctdf$frac = ctdf$ncell / ctdf$totcell
    ctdf = merge(ctdf, metadata[,c('rind','plaq_n_regional','plaq_d_regional','nft_regional')])
    ctdf = ctdf[!is.na(ctdf$nft_regional),]
    ctdf$logit.frac = logit(ctdf$frac)

    # Logit frac vs. region + reponse:
    resp.vars = c('plaq_n_regional','plaq_d_regional','nft_regional')
    m1 = lm(as.formula(paste0('logit.frac ~ ', ct.var, '* region + 1')), ctdf)
    vdf = c()
    for (var in resp.vars){
        m2 = lm(as.formula(paste('logit.frac ~ ', var, ' * ', ct.var, ' + ', ct.var, '* region + 1')), ctdf)
        pval = anova(m1, m2)$P[2]
        cat(var, '\t', pval,'\n')
        cf = coefficients(summary(m2))
        print(cf[grep(var, rownames(cf)),])
        sdf = as.data.frame(cf[grep(var, rownames(cf)),c(1,4)])
        sdf$var = var
        sdf$rn = sub(var, '', rownames(sdf))
        sdf$rn[sdf$rn == ''] = 'var'
        sdf$log10p = -log10(sdf$'Pr(>|t|)')
        vdf = rbind(vdf, sdf)
        cat('\n')
    }
    rownames(vdf) = NULL
    vdf$rn =sub(paste0(":", ct.var), '', vdf$rn)
    vdf$rn[vdf$rn == 'var'] = 'Intercept'
    if (ct.var == 'full.exttype'){ 
        vdf$rn = factor(vdf$rn, levels=c('Intercept', rev(as.character(odf$full.exttype))))
    }
    pwide = spread(vdf[,c('var','rn','log10p')], rn, log10p)
    rownames(pwide) = pwide[,1]
    rwide = spread(vdf[,c('var','rn','Estimate')], rn, Estimate)
    # TODO: Reorder.
    # Factor for plotting:
    vdf$var = factor(vdf$var)
    vdf$rn = factor(vdf$rn)

    w = 1.5 + dim(rwide)[1] * 3 / 5
    h = 2 + dim(rwide)[2] / 5
    gp = ggplot(vdf, aes(var, rn, size=log10p, color=Estimate)) + 
        # scale_color_continuous(colryb) + 
        geom_point() + 
        theme_pubr() + 
        scale_color_distiller(palette='RdBu') + 
        labs(x='Response Variable', y='Cell type') + 
        theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'.png'), gp, dpi=450, units='in', width=w, height=h)
    ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'.pdf'), gp, width=w, height=h)

    # Plot some of these top ones:
    if (ct.var == 'full.exttype'){
        keep.vals = c('Ast','Oli','NRGN','Granule','CPEC','CAM', 'CA1_Pyr', 'SST','PVALB')
    } else if (ct.var == 'major.celltype'){
        keep.vals = unique(ctdf[[ct.var]]) 
    }
    subdf = ctdf[ctdf[[ct.var]] %in% keep.vals,]

    scale = 4
    h = 2 * scale 
    w = length(keep.vals) / 2 * scale
    gs = ggplot(subdf, aes(nft_regional, logit.frac, color=region)) + 
        facet_wrap(as.formula(paste0('~', ct.var)), nrow=2) + 
        geom_point() + geom_smooth(method = 'lm') + 
        scale_color_manual(values=reg.cols) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'_examples_nft.png'), gs, dpi=450, units='in', width=w, height=h)

    gs = ggplot(subdf, aes(plaq_d_regional, logit.frac, color=region)) + 
        facet_wrap(as.formula(paste0('~', ct.var)), nrow=2) + 
        geom_point() + geom_smooth(method = 'lm') + 
        scale_color_manual(values=reg.cols) + 
        theme_pubr()
    ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'_examples_plaq_d.png'), gs, dpi=450, units='in', width=w, height=h)


    slong = subdf[,c('frac','nft_regional','plaq_d_regional','plaq_n_regional', ct.var,'region', 'rind')]
    swide = gather(slong, var, value, -frac, -ct.var, -region, -rind)

    for (val in keep.vals){
        print(val)
        valstr = sub("/",".", val)
        gs = ggplot(subdf[subdf[[ct.var]] == val,], aes(nft_regional, frac, fill=region)) + 
            facet_wrap(~region, nrow=1) + 
            geom_point() + geom_smooth(method = 'lm') + 
            labs(x='PHFTau Density (nft)', y=paste0(val, ' % of total cells')) + 
            scale_y_continuous(labels=scales::percent) + 
            scale_fill_manual(values=reg.cols) + 
            theme_pubr()
        gs2 = gs + labs(x='PHFTau Density (nft)', y=paste0(val, ' % of total cells\n(log10 scale)')) + 
            scale_y_log10(labels=scales::percent) + 
        ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'_',valstr, '_nft.png'), gs, dpi=450, units='in', width=w, height=5)
        ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'_',valstr, '_nft_log10.png'), gs2, dpi=450, units='in', width=w, height=5)

        gs = ggplot(swide[subdf[[ct.var]] == val,], aes(value, frac, fill=region)) + 
            facet_grid(var~region) + 
            geom_point() + geom_smooth(method = 'lm') + 
            labs(x='Pathology Score', y=paste0(val, ' % of total cells')) + 
            scale_y_continuous(labels=scales::percent) + 
            scale_fill_manual(values=reg.cols) + 
            theme_pubr()
        # gs2 = gs + labs(x='Pathology Score', y=paste0(val, ' % of total cells\n(log10 scale)')) + 
        #     scale_y_log10(labels=scales::percent) + 
        ggsave(paste0(imgpref, 'metadata_AD_fractions_regional_', ct.var ,'_',valstr, '_3path.png'), gs, dpi=450, units='in', width=w, height=7)
        # ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_nft_log10.png'), gs2, dpi=450, units='in', width=w, height=5)
    }





    print(1 * (vwide > 1))
    # Make heatmaps:
}


# -----------------------------------------------
# Repeat analysis, but with regional coefficient.
# -----------------------------------------------
for (ct.var in vars) {
    totdf = agg.rename(barcode ~ rind + region, cellmeta, length, 'totcell')
    ctdf = agg.rename(as.formula(paste0('barcode ~ ', ct.var, '+ rind + region')), cellmeta, length, 'ncell')
    ctdf = merge(ctdf, totdf)
    ctdf$frac = ctdf$ncell / ctdf$totcell
    ctdf = merge(ctdf, metadata[,c('rind','plaq_n_regional','plaq_d_regional','nft_regional', 'nft','plaq_d','plaq_n', 'msex', 'cogdx', 'age_death')])
    # ctdf = ctdf[!is.na(ctdf$nft_regional),]
    ctdf$logit.frac = logit(ctdf$frac)

    # Run through each independently, see if overall response to var + regional response to var
    ovdf = c()
    intdf = c()
    for (var in c('plaq_n_regional','plaq_d_regional','nft_regional', 'nft','plaq_d','plaq_n')){
        for (val in unique(ctdf[[ct.var]])){
            cat(var, "\t", val, "\t")
            subdf = ctdf[ctdf[[ct.var]] == val,]
            subdf = subdf[!is.na(subdf[[var]]),]
            # m0 = lm(as.formula(paste0('logit.frac ~ region + 1')), subdf)
            # m1 = lm(as.formula(paste0('logit.frac ~ ', var ,'+ region + 1')), subdf)
            # m2 = lm(as.formula(paste0('logit.frac ~ region * ', var, ' + region + 1')), subdf)
            m0 = lm(as.formula(paste0('logit.frac ~ region + age_death + msex + 1')), subdf)
            m1 = lm(as.formula(paste0('logit.frac ~ ', var ,'+ region + age_death + msex + 1')), subdf)
            m2 = lm(as.formula(paste0('logit.frac ~ region * ', var, ' + region + age_death + msex + 1')), subdf)

            m3 = lm(as.formula(paste0('logit.frac ~ msex * ', var, ' + region + age_death + 1')), subdf)
            m4 = lm(as.formula(paste0('logit.frac ~ region * msex * ', var, ' + region + age_death + 1')), subdf)
            av = anova(m0, m1, m2)
            av2 = anova(m1, m3)
            av3 = anova(m2, m4)
            if ((!is.na(av2$Pr[2])) & av2$Pr[2] < 0.1){
                print(av2$Pr[2])
                print(av3$Pr[2])
            }
            pvals = av$Pr[2:3]
            cat(round(pvals[1],3), '\t', round(pvals[2], 3), '\n')
            # Extract coefficients:
            cf = coefficients(summary(m1))
            est1 = cf[var,'Estimate']
            # cf = coefficients(summary(m2))
            # est2 = cf[grep(var, rownames(cf)),'Estimate']
            ovdf = rbind(ovdf, data.frame(var=var, val=val, p=pvals[1], log10p=-log10(pvals[1]), coeff=est1))
            intdf = rbind(intdf, data.frame(var=var, val=val, p=pvals[2], log10p=-log10(pvals[2])))
        }
    }

    ovdf = ovdf[order(ovdf$p),]
    intdf = intdf[order(intdf$p),]

    # Places where disagreement with regional/overall:
    ovdf$short.var = sub('_regional', '', ovdf$var)

    qqplot(ppoints(nrow(ovdf)), ovdf$p)
    qqplot(1:nrow(ovdf)/nrow(ovdf), ovdf$p)
    plot(1:nrow(ovdf)/nrow(ovdf), ovdf$p)
    abline(0,1)

    hist(ovdf[ovdf$var == 'nft','p'])
    head(ovdf[ovdf$var == 'nft',])
    qvalue(ovdf[ovdf$var == 'nft','p'])$qvalue

    osdf = ovdf[ovdf$var == ovdf$short.var,]
    osdf$q = qvalue(osdf$p)$qvalues

    osdf = ovdf[ovdf$short.var == 'nft',]
    osdf$q = qvalue(osdf$p)$qvalues




    # IF SIG:
    #     ggplot(subdf, aes(nft, logit.frac, fill=region)) + 
    #             facet_grid(~region) + 
    #             geom_point() + geom_smooth(method = 'lm') + 
    #             labs(x='Pathology Score', y=paste0(val, ' % of total cells')) + 
    #             scale_y_continuous() + 
    #             scale_fill_manual(values=reg.cols) + 
    #             theme_pubr()

    val = 'SMC'
    valstr = sub("/", ".", val)
    var = 'plaq_n_regional'
    subdf = ctdf[ctdf[[ct.var]] == val,]
    subdf = subdf[!is.na(subdf[[var]]),]

    gs = ggplot(subdf, aes_string(var, 'logit.frac', fill='region')) + 
        facet_grid(msex~region) + 
        geom_point() + geom_smooth(method = 'lm') + 
        labs(x='Pathology Score', y=paste0(val, ' % of total cells')) + 
        scale_y_continuous() + 
        scale_fill_manual(values=reg.cols) + 
        theme_pubr() + theme(legend.position='none')

    ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_',var,'.png'), gs, dpi=450, units='in', width=10, height=3)


    val = 'SST'
    val='Exc RELN GPC5'
    val='Exc TRPC6 ANO2'
    val='End'
    val='Inh CUX2 MSR1'
    val='Inh L3-5 SST MAFB'
    val='T cells'
    val='Oli RASGRF1'
    val='Ast GRM3'
    val='Exc NRGN'

    valstr = gsub(" ", "_", val)
    # var = 'plaq_n'
    var = 'nft'
    subdf = ctdf[ctdf[[ct.var]] == val,]
    subdf = subdf[!is.na(subdf[[var]]),]
    gs = ggplot(subdf, aes_string(var, 'logit.frac', fill='region')) + 
        facet_grid(~region) + 
        geom_point() + geom_smooth(method = 'lm') + 
        labs(x='Pathology Score', y=paste0(val, ' % of total cells')) + 
        scale_y_continuous() + 
        scale_fill_manual(values=reg.cols) + 
        theme_pubr()

    ggsave(paste0(imgpref, 'metadata_AD_fractions_', ct.var ,'_',valstr, '_',var,'.png'), gs, dpi=450, units='in', width=10, height=5)





    # Logit frac vs. region + reponse:
    resp.vars = c('plaq_n_regional','plaq_d_regional','nft_regional')
    m1 = lm(as.formula(paste0('logit.frac ~ ', ct.var, '* region + 1')), ctdf)
    vdf = c()
    for (var in resp.vars){
        # m2 = lm(as.formula(paste('logit.frac ~ ', var, ' * ', ct.var, ' + ', ct.var, '* region + 1')), ctdf)
        m2 = lm(as.formula(paste('logit.frac ~ ', var, ' * ', ct.var, ' + ', ct.var, '* region + 1')), ctdf)
        pval = anova(m1, m2)$P[2]
        cat(var, '\t', pval,'\n')
        cf = coefficients(summary(m2))
        print(cf[grep(var, rownames(cf)),])
        sdf = as.data.frame(cf[grep(var, rownames(cf)),c(1,4)])
        sdf$var = var
        sdf$rn = sub(var, '', rownames(sdf))
        sdf$rn[sdf$rn == ''] = 'var'
        sdf$log10p = -log10(sdf$'Pr(>|t|)')
        vdf = rbind(vdf, sdf)
        cat('\n')
    }

    rownames(vdf) = NULL
    vdf$rn =sub(paste0(":", ct.var), '', vdf$rn)
    vdf$rn[vdf$rn == 'var'] = 'Intercept'
    if (ct.var == 'full.exttype'){ 
        vdf$rn = factor(vdf$rn, levels=c('Intercept', rev(as.character(odf$full.exttype))))
    }
    pwide = spread(vdf[,c('var','rn','log10p')], rn, log10p)
    rownames(pwide) = pwide[,1]
    rwide = spread(vdf[,c('var','rn','Estimate')], rn, Estimate)
    # TODO: Reorder.
    # Factor for plotting:
    vdf$var = factor(vdf$var)
    vdf$rn = factor(vdf$rn)


