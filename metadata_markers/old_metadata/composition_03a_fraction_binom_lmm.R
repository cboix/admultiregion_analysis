#!/usr/bin/R
# -----------------------------------------------------------
# Calculate the fraction differences, by binomial regression:
# As in:
# https://github.com/vals/Blog/tree/master/201127-cell-count-glm
# Updated 01/19/2021 
# Updated 03/15/2021 
# Final version is composition_02_binom_emmeans.R
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(qvalue)
library(lme4)
library(emmeans)
library(nnet) # multinom
library(mlogit) # mlogit

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = plotdir
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

asform = function(x){as.formula(paste(x, collapse=" "))}

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

# Get the pathology mapped to each region:
pqdf = NULL
for (path in c('nft','plaq_d','plaq_n')){
    regmap = c('AG','HC','PFC','MT','EC')
    names(regmap) = c('ag','hip','mf','mt','ec')
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    if (is.null(pqdf)){
        pqdf = slong[,c('rind','value','region')]
        names(pqdf)[2] = path
    } else { 
        sub.pqdf = slong[,c('rind','value','region')]
        names(sub.pqdf)[2] = path
        pqdf = merge(pqdf, sub.pqdf)
    }
}

# Mapping:
nrmeta = unique(metadata[,c('rind','projid','region','Apoe_e4','msex','pmi','age_death', 'niareagansc', 'cogdx')])
nrmeta = merge(nrmeta, pqdf, all.x=TRUE)
rownames(nrmeta) = nrmeta$rind
nrmeta$nrad = 'AD'
nrmeta$nrad[nrmeta$niareagansc > 2] = 'Control'
nrmeta$nrad = factor(nrmeta$nrad, c('Control','AD'))
nrmeta$cogdxad = 'AD'
nrmeta$cogdxad[nrmeta$cogdx %in% c(1:3)] = 'Control'
nrmeta$cogdxad = factor(nrmeta$cogdxad, c('Control','AD'))

# -------------------------------------
# Format data for abundance estimation:
# -------------------------------------
clsval = 'cell_type_high_resolution'
# clsval = 'minor.celltype'
# clsval = 'major.celltype'

cmap = unique(cellmeta[,c('major.celltype',clsval)])
names(cmap)[2] = 'cls'
# Aggregate counts so that we have 0 tot:
combdf = expand.grid(cls=unique(cellmeta[[clsval]]), region=unique(cellmeta$region), projid=unique(cellmeta$projid))
ctdf = agg.rename(asform(c('barcode ~ projid + region +', clsval)), cellmeta, length, 'Count')
names(ctdf)[3] = 'cls'
ctdf = merge(ctdf, cmap)
combdf = merge(combdf, cmap)
ctdf = merge(ctdf, combdf, all.y=TRUE)
ctdf$Count[is.na(ctdf$Count)] = 0

# Remove region and individual hugely inflated for fibroblasts:
ctdf = ctdf[!(ctdf$projid == '50106280' & ctdf$region == 'HC'),]

totdf = agg.rename(Count ~ projid + region, ctdf, sum, 'Total')
ctdf = merge(ctdf, totdf) # Removes a couple missing projid x region comb.
ctdf$Other = ctdf$Total - ctdf$Count
if (clsval == 'cell_type_high_resolution'){
    ctdf = ctdf[ctdf$cls != 'Exc SV2C LINC02137',] # Remove problematic batch
}

# For subsetting; relative to all cells:
if (clsval == 'cell_type_high_resolution'){
    substr = 'Exc'
} else {
    substr = NULL
}
if (!is.null(substr)){
    sts = unique(cellmeta[cellmeta$major.celltype == substr,clsval])
    ctdf = ctdf[ctdf$cls %in% sts,]
    clsstr = paste0(clsval,'_', sub("/","_",substr))
} else { clsstr = clsval }

# Add other attributes:
# ctdf = merge(ctdf, unique(metadata[,c('projid','nft','plaq_d','plaq_n', 'braaksc','cogdx','niareagansc','msex','age_death','pmi')]))
cf = unique(metadata[,c('projid','region','rind','braaksc','cogdx','niareagansc','msex','age_death','pmi')])
ctdf = merge(ctdf, cf)
ctdf = merge(ctdf, pqdf, all.x=TRUE)

# Format attributes:
ctdf$projid = factor(ctdf$projid)
ctdf$nrad = 'AD'
ctdf$nrad[ctdf$niareagansc > 2] = 'CTRL'
ctdf$nrad = factor(ctdf$nrad, c('CTRL','AD'))
ctdf$cogdxad = 'AD'
ctdf$cogdxad[ctdf$cogdx %in% c(1:3)] = 'CTRL'
ctdf$cogdxad = factor(ctdf$cogdxad, c('CTRL','AD'))

# Model with binomial:
pathval = 'nrad'
# pathval = 'cogdxad'
# pathval = 'plaq_n'
formula = asform(c('cbind(Count, Other) ~ cls *(',pathval,'+ pmi + msex + age_death + region)')) # Orig
# formula = asform(c('cbind(Count, Other) ~ cls *(',pathval,'*region + pmi + msex + age_death + region)')) # By region interaction
# formula = asform(c('cbind(Count, Other) ~ cls *(',pathval,' * region * msex + pmi + age_death )')) # Sex 
fit <- glm(formula = formula, family = 'quasibinomial', data = ctdf)
# fit <- glm(formula = formula, family = binomial(link='logit'), data = ctdf)

# Alternative to quasibinom; mixed eff. binom. (bad fit, still overestimates)
# formula = asform(c('cbind(Count, Other) ~ cls *(',pathval,'+ pmi + msex + age_death + region) + (1|projid)')) # Orig
# fit <- glmer(formula = formula, family = binomial, data = ctdf)

emform = asform(c('revpairwise ~', pathval, '|cls '))
# emform = asform(c('revpairwise ~ ',pathval,'|(cls * region)'))
emm1 <- emmeans(fit, specs=emform)

emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results
c_results = c_results[order(c_results$odds.ratio),]

# See if residuals change along prediction:
ctdf$pred = predict(fit, ctdf)
# Residuals in logit space suggests not overdispersed.
gplot = ggplot(ctdf, aes(pred, log(Count / Total - exp(pred)))) + 
    geom_point() + 
    geom_smooth(method='loess') + 
    theme_pubr()

# Plot AD x Region:
p.cut = 0.01
c_results$p.adj = p.adjust(c_results$p.value, 'fdr')
c_results$col = 1 * (c_results$p.adj < p.cut) + 2 * (c_results$odds.ratio > 1 )
c_results = c_results[order(c_results$odds.ratio, decreasing=F),]
c_results$cls = factor(c_results$cls, levels=unique(c_results$cls))
labdf = c_results[c_results$p.adj < p.cut,]
c_results = c_results[c_results$SE < 20,]
c_results = merge(c_results, data.frame(region=unique(c_results$region), jt=seq(.25,.75,length.out=6)))

gplot = ggplot(c_results, aes(x = odds.ratio, y=cls, col=region)) +
    # facet_wrap(~region) + 
    geom_point() +
    geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend =cls)) +
    geom_text(data=labdf, aes(x=max(labdf$odds.ratio) * 1.1, y=cls, label=sprintf("OR=%0.2f\tp=%0.1e",odds.ratio,p.adj)),col='black') +
    scale_x_log10() +
    scale_color_manual(values=reg.cols) + 
    geom_vline(xintercept=1, lty='dashed') + 
    theme_pubr() + theme(legend.position = 'none') + 
    labs(x = c_results$contrast[1], y=clsval)

# Observe diff for OPC / Oli (cross-region):
# gplot = ggplot(ctdf[ctdf$cls %in% c('Opc','Oli'),], aes(region, Count / Total, color=nrad)) + 
#     facet_wrap(~cls, scale='free_y') + 
#     geom_boxplot() + 
#     scale_color_manual(values=colvals[['nrad']]) +
#     theme_pubr()

# ----------------------------------------------------------------
# Alternatively multinomial logit (due to over-pred. at high pct):
# ----------------------------------------------------------------
# ml = ctdf
# ml$cls = factor(ml$cls)
# ml = unique(ml[,c('cls','nrad','pmi','msex','age_death','region','projid','Count')])
# ml$age_rescaled = ml$age_death / 100
# ml$cls2 <- relevel(ml$cls, ref = "Ast")
# test <- multinom(cls2 ~ nrad + pmi + msex + age_rescaled + region , data = ml, summ=2, weights=ml$Count)

# # Wald test p-values:
# z <- summary(test)$coefficients/summary(test)$standard.errors
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# sprintf("%0.2e",p[,'nradAD'])

# -------------------
# Try library mlogit:
# -------------------
# Very annoying indexing structure...
# data("Heating", package = "mlogit")
# H <- dfidx(Heating, choice = "depvar", varying = c(3:12))
# m <- mlogit(depvar ~ ic + oc | 0, H)
# summary(m)
# m <- mlogit(cls ~ nrad + region | 0, ml)

# ----------------
# Plot as volcano:
# ----------------
p.cut = 0.01
c_results$p.adj = p.adjust(c_results$p.value, 'fdr')
c_results$col = 1 * (c_results$p.adj < p.cut) + 2 * (c_results$odds.ratio > 1 )
labdf = c_results[c_results$p.adj < p.cut,]
gplot = ggplot(c_results, aes(x = odds.ratio, y = -log10(p.adj), col=factor(col))) +
    geom_point() +
    geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.adj))) +
    geom_text_repel(data=labdf, aes(x=odds.ratio, y=-log10(p.adj), label=cls, col=factor(col)),max.overlaps=30) +
    scale_x_log10() +
    scale_color_manual(values =c('0' ='grey75','1'='royalblue','2'='grey75','3'='indianred')) + 
    geom_vline(xintercept=1) + 
    scale_y_continuous(expand=c(0,0.1)) + 
    theme_pubr() +
    theme(legend.position = 'none') + 
    labs(x = c_results$contrast[1])
ggsave(paste0(imgpref, 'metadata_AD_fractions_overall_', clsstr, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=6, height=8)
ggsave(paste0(imgpref, 'metadata_AD_fractions_overall_', clsstr, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=6, height=8)


# -----------------------------------------
# Plot the effects in order for legibility:
# -----------------------------------------
p.cut = 0.01
c_results$p.adj = p.adjust(c_results$p.value, 'fdr')
c_results$col = 1 * (c_results$p.adj < p.cut) + 2 * (c_results$odds.ratio > 1 )
c_results = c_results[order(c_results$odds.ratio, decreasing=F),]
c_results$cls = factor(c_results$cls, levels=c_results$cls)
labdf = c_results[c_results$p.adj < p.cut,]
gplot = ggplot(c_results, aes(x = odds.ratio, y=cls, col=factor(col))) +
    geom_point() +
    geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend =cls)) +
    geom_text(data=labdf, aes(x=max(labdf$odds.ratio) * 1.1, y=cls, label=sprintf("OR=%0.2f\tp=%0.1e",odds.ratio,p.adj)),col='black') +
    scale_x_log10() +
    scale_color_manual(values =c('0' ='grey75','1'='royalblue','2'='grey75','3'='indianred')) + 
    geom_vline(xintercept=1, lty='dashed') + 
    theme_pubr() + theme(legend.position = 'none') + 
    labs(x = c_results$contrast[1], y=clsval)
h = nrow(c_results) / 30 * 8
ggsave(paste0(imgpref, 'metadata_AD_oddratios_overall_', clsstr, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=6, height=h)
ggsave(paste0(imgpref, 'metadata_AD_oddratios_overall_', clsstr, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=6, height=h)

# -----------------------------------------
# Plot the effects in order for legibility:
# -----------------------------------------
p.cut = 0.01
labdf = c_results[c_results$p.adj < p.cut,]
gplot = ggplot(labdf, aes(x = odds.ratio, y=cls, col=factor(col))) +
    geom_point() +
    geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend =cls)) +
    geom_text(data=labdf, aes(x=max(labdf$odds.ratio) * 1.1, y=cls, label=sprintf("OR=%0.2f\tp=%0.1e",odds.ratio,p.adj)),col='black') +
    scale_x_log10() +
    scale_color_manual(values =c('0' ='grey75','1'='royalblue','2'='grey75','3'='indianred')) + 
    geom_vline(xintercept=1, lty='dashed') + 
    theme_pubr() + theme(legend.position = 'none') + 
    labs(x = c_results$contrast[1], y=clsval)
h = nrow(labdf) / 30 * 8 + 1
ggsave(paste0(imgpref, 'metadata_AD_oddratios_signif_', clsstr, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=6, height=h)
ggsave(paste0(imgpref, 'metadata_AD_oddratios_signif_', clsstr, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=6, height=h)

# Plot the percentage barplots here:
cdf = aggregate(Count ~ cls + region,ctdf, sum)
cdf$cls = factor(cdf$cls, levels=c_results$cls)
gplot = ggplot(cdf, aes(cls, Count, fill=region)) + 
    geom_bar(stat='identity', color=NA, position='fill') +
    scale_fill_manual(values=reg.cols) + 
    theme_pubr() + theme(legend.position = 'none') + 
    scale_y_continuous(expand=c(0,0)) + 
    coord_flip()
ggsave(paste0(imgpref, 'metadata_AD_proportions_overall_', clsstr, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=4, height=h)
ggsave(paste0(imgpref, 'metadata_AD_proportions_overall_', clsstr, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=4, height=h)

# ------------------------------
# Plot odds ratio vs. abundance: 
# ------------------------------
emm2 <- emmeans(fit, specs =asform(c('~cls')))
emm2 %>% summary(type = 'response') ->  mean_probs
mean_probs = mean_probs[,c('cls','prob')]
m_results = merge(c_results, mean_probs, all.x=TRUE)

submdf = m_results[abs(log(m_results$odds.ratio)) > log(1.5),]
gplot = ggplot(data=m_results, aes(x = prob, y = odds.ratio, color=factor(col))) +
    geom_point() +
    geom_text_repel(data=submdf, aes(label = cls)) +
    scale_x_log10(labels=scales::percent) +
    scale_y_log10() +
    scale_color_manual(values =c('grey75','royalblue','grey75','indianred')) + 
    theme_pubr() + theme(legend.position='none') + 
    labs(y = paste(m_results$contrast[1], '(odds ratio)'), x = 'Average abundance (probability)')
ggsave(paste0(imgpref, 'metadata_AD_fractions_or_vs_abund_', clsstr, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=6, height=6)
ggsave(paste0(imgpref, 'metadata_AD_fractions_or_vs_abund_', clsstr, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=6, height=6)

pdf = m_results[,c('cls','contrast','odds.ratio','asymp.LCL','asymp.UCL','p.value','prob', 'col')]
pdf$region = 'All' 






# -------------------------------------------
# Repeat each of these at the regional level:
# NOTE: Not used.
# -------------------------------------------
for (reg in regions[regions!='MB']){
    print(reg)
    subdf = ctdf[ctdf$region == reg,]
    agdf = aggregate(asform(c('Count ~', clsval)), subdf, sum)
    agdf$frac = agdf$Count / sum(agdf$Count)
    keep.ct = agdf[agdf$frac > 0.0005, clsval] 
    subdf = subdf[subdf[[clsval]] %in% keep.ct,]
    # Prune to only celltypes with significant presence in the region:

    # Model with binomial:
    pathval = 'nrad'
    formula = asform(c('cbind(Count, Other) ~ ',clsval,'*',pathval,'+ ',clsval,'*pmi + ',clsval,'*msex +',clsval,'* age_death'))
    fit <- glm(formula = formula, family = 'binomial', data = subdf)

    emform = asform(c('revpairwise ~', pathval, '|',clsval))
    emm1 <- emmeans(fit, specs=emform)
    emm1$contrasts %>%
        summary(infer = TRUE, type = 'response') %>%
        rbind() %>%
        as.data.frame() -> c_results

    p.cut = 1e-10
    c_results$col = 1 * (c_results$p.value < p.cut) + 2 * (c_results$odds.ratio > 1 )
    gplot = ggplot(c_results, aes(x = odds.ratio, y = -log10(p.value), col=factor(col))) +
        geom_point() +
        geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.value))) +
        geom_text_repel(data=c_results[c_results$p.value < p.cut,], aes(x=odds.ratio, y=-log10(p.value), label=cell_type_high_resolution, col=factor(col))) +
        scale_x_log10() +
        scale_color_manual(values =c('grey75','royalblue','grey75','indianred')) + 
        # geom_hline(yintercept=-log10(p.cut), ) + 
        geom_vline(xintercept=1) + 
        scale_y_continuous(expand=c(0,5)) + 
        theme_pubr() +
        theme(legend.position = 'none') + 
        labs(x = c_results$contrast[1])
    ggsave(paste0(imgpref, 'metadata_AD_fractions_reg_', reg, '_volcano_', clsval, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=6, height=8)
    ggsave(paste0(imgpref, 'metadata_AD_fractions_reg_', reg, '_volcano_', clsval, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=6, height=8)


    # Plot differences: 
    emm2 <- emmeans(fit, specs =asform(c('~',clsval)))
    emm2 %>% summary(type = 'response') ->  mean_probs
    mean_probs = mean_probs[,c(clsval,'prob')]
    m_results = merge(c_results, mean_probs, all.x=TRUE)

    submdf = m_results[abs(log(m_results$odds.ratio)) > log(1.5),]
    gplot = ggplot(data=m_results, aes(x = prob, y = odds.ratio, color=factor(col))) +
        geom_point() +
        geom_text_repel(data=submdf, aes(label = cell_type_high_resolution)) +
        scale_x_log10(labels=scales::percent) +
        scale_y_log10() +
        scale_color_manual(values =c('grey75','royalblue','grey75','indianred')) + 
        theme_pubr() + theme(legend.position='none') + 
        labs(y = paste(m_results$contrast[1], '(odds ratio)'), x = 'Average abundance (probability)')
    ggsave(paste0(imgpref, 'metadata_AD_fractions_reg_',reg,'_or_vs_abund_', clsval, '_vs_',pathval,'.png'), gplot, dpi=450, units='in', width=6, height=6)
    ggsave(paste0(imgpref, 'metadata_AD_fractions_reg_', reg, '_or_vs_abund_', clsval, '_vs_',pathval,'.pdf'), gplot, dpi=450, units='in', width=6, height=6)

    # Write: 
    df = m_results[,c(clsval,'contrast','odds.ratio','asymp.LCL','asymp.UCL','p.value','prob', 'col')]
    df$region = reg
    pdf = rbind(pdf, df)
}


# --------------------------
# Plot as a heatmap as well:
# --------------------------
owide = spread(pdf[,c('cell_type_high_resolution','region','odds.ratio')], cell_type_high_resolution, odds.ratio)
pwide = spread(pdf[,c('cell_type_high_resolution','region','p.value')], cell_type_high_resolution, p.value)
omat = as.matrix(owide[,-1])
rownames(omat) = owide[,1]
omat = log2(omat)

mx = 2
omat[omat > mx] = mx
omat[omat < -mx] = -mx

regs = regions[regions != 'MB']
rmat = omat[regs,]
amat = omat['All',,drop=F]

# TODO: Three panel, ordered + split by ct + in each ct type.

sp = 0.1
png(paste0(imgpref, 'metadata_AD_fractions_reg_binomlmm_compheat.png'), res=400, units='in', width=2.5, height=8)
layout(matrix(1:3,nrow=1), widths=c(8,6,1))
par(mar=c(sp,sp,1, sp))
image(rmat,zlim=c(-mx, mx), col='white', useRaster=T, axes=F)
xat = seq(0,1, length.out=nrow(rmat))
yat = seq(0,1, length.out=ncol(rmat))
text(x=parpos(1,-1), y=yat, colnames(rmat), xpd=TRUE, adj=1, cex=.8)
#
par(mar=c(sp,sp,1, sp))
image(rmat,zlim=c(-mx, mx), col=rev(colrb), useRaster=T, axes=F)
box(lwd=.25)
text(x=xat, y=parpos(2,-1.008), labels=rownames(rmat), xpd=TRUE, adj=.5, cex=.9)
# 
par(mar=c(sp,sp,1, sp))
image(amat,zlim=c(-mx, mx), col=rev(colrb), useRaster=T, axes=F)
text(x=0, y=parpos(2,-1.008), labels='All', xpd=TRUE, adj=.5, cex=.9)
box(lwd=.25)
dev.off()








