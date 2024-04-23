#!/usr/bin/R
# ----------------------------------
# Prelim plots scDRS
# Updated 10/26/2022
# ----------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)
library(viridis)
library(ggplot2)
library(ggpubr)

library(ComplexHeatmap)
library(circlize)
options(width=170)

# Directories:
scddir = paste0(sdbdir, 'scDRS/')
plotdir = paste0(imgdir, 'gwas/')
imgpref = paste0(plotdir, 'scDRS_')
cmd = paste('mkdir -p', plotdir, scddir)
system(cmd)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Load the module scores information for microglia:
# -------------------------------------------------
runset = 'Mic_Immune'
graph_id = 'boot'
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, TRUE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))

mic.bcs = cellmeta$barcode[cellmeta$minor.celltype == 'Mic']
scoremat = scoremat[mic.bcs,]


# Read scDRS per-cell scores for AD only
# --------------------------------------
ad.gwas = "PASS_Alzheimers_Jansen2019"
gsfile = paste0(scddir, 'results/full.scores-', ad.gwas, '.tsv.gz')
mic.gsfile = paste0(scddir, 'results/full.scores-', ad.gwas, '.Mic.tsv.gz')
if (!file.exists(mic.gsfile)){
    df = read.delim(gsfile, header=T)
    df$fdr = p.adjust(df$pval, method='BH') # Adjust in same manner as scDRS
    names(df)[1] = 'barcode'
    rownames(df) = df$barcode
    # Subset to microglia barcodes only
    df = df[mic.bcs,]
    write.table(df, gzfile(mic.gsfile), quote=F, sep="\t", row.names=F)
} else { 
    df = read.delim(mic.gsfile, header=T)
}


# Merge with microglia data:
# --------------------------
df = merge(df, cellmeta[,c('barcode','rind','projid','region','cell_type_high_resolution', 'U1', 'U2')])
rownames(df) = df$barcode  # Match to mic barcodes again
df = df[mic.bcs,]

gp = ggplot(df, aes(cell_type_high_resolution, norm_score, fill=region)) + 
    geom_boxplot() + 
    stat_compare_means() + 
    scale_fill_manual(values=reg.cols) + 
    theme_pubr()

# mx = max(abs(df$norm_score))
# col_fun = colorRamp2(c(-mx, 0, mx), c('blue', "white", 'red'))
# plot(df$U1, df$U2, col=col_fun(df$norm_score), pch='.')


# Score against modules:
# ----------------------
fdr.cutoffs = c(0.01, 0.05, 0.1, 0.2)
statdf = c()
lmdf = c()
for (mod in colnames(scoremat)){
    x = scoremat[,mod]
    z = scale(x)
    ind = z > 2.5
    df$module_score = x
    fit = glm(norm_score ~ module_score, df, family='gaussian')
    cfit = as.data.frame(coefficients(summary(fit)))
    names(cfit) = c('est','se','t','p')
    cfit$mod = mod
    lmdf = rbind(lmdf, cfit[2,,drop=FALSE])
    for (fdr.cutoff in fdr.cutoffs){
        nf = sum(df$fdr[ind] < fdr.cutoff)
        mfdf = data.frame(
            nint = nf, nfdr = sum(df$fdr < fdr.cutoff),
            nmod = sum(ind), mod = mod, fdr = fdr.cutoff)
        statdf = rbind(statdf, mfdf)
    }
}

statdf$module = as.numeric(sub("M","", statdf$mod))
statdf = merge(statdf, mmap)
statdf$ntot = nrow(df)
statdf$frac = statdf$nint / statdf$nmod
statdf = statdf[order(statdf$frac, decreasing=T),]

hgdf = statdf[,c('nint','nmod','nfdr','ntot')]
statdf$hg.p = apply(hgdf, 1, run.hyper)
statdf$hg.padj = p.adjust(statdf$hg.p, 'BH')
sigdf = statdf[(statdf$fdr == 0.05) & (statdf$hg.padj < 0.01),]
print(sigdf)
sigdf$mname = factor(sigdf$mname, levels=rev(sigdf$mname))

gp = ggplot(sigdf, aes(frac, mname)) + 
    geom_bar(fill='grey75', stat='identity') + 
    labs(x='% of immune cells with scDRS FDR < 0.05', y='Microglia/Immune module', title='Microglia/Immune modules significantly enriched for cells with scDRS < 0.05') +
    scale_x_continuous(expand=c(0,0), labels=scales::percent) + 
    theme_pubr()
pltprefix = paste0(imgpref, 'top_modules_barplot')
saveGGplot(gp, pltprefix, w=6.5, h=3.25)




lmdf$module = as.numeric(sub("M","", lmdf$mod))
lmdf = merge(lmdf, mmap)
lmdf = lmdf[order(lmdf$t, decreasing=T),]
head(lmdf)

outfile = paste0(scddir, 'module_enr_microglia-', ad.gwas, '.tsv')
write.table(statdf, outfile, quote=F, row.names=F, sep="\t")

