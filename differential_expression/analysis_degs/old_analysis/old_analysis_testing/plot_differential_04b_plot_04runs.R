#!/usr/bin/R
# ----------------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Using Liang's Nebula package to run the fast NBLMM
# Basic overall + per-region interactions, for relevant subtypes
# Updated: 10/28/2020
# 
# Actually plot the runs for plot_differential_04_filter_low_lmm.R
# ----------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(rhdf5)
library(nebula)
library(qvalue)
library(Matrix)

library(ggplot2)
library(ggpubr)
library(scales)
library(ggrepel)

celltype = 'Mic_Immune'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need celltype")
} else {        
    celltype = args[1]
}
print(celltype)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
pathlist = c('nft','plaq_d','plaq_n')

# Function to estimate end time of a loop:
est.endtime <- function(start, last, done, total){
    curr = proc.time() # Get current time
    # Calculate how much time has passed:
    last.step = (proc.time() - last)[3]
    elapsed = (proc.time() - start)[3]
    cat(paste0(round(last.step,1),'s\t'))
    cat(paste0(round(elapsed,1),'s\t'))
    # Estimate the remaining time:
    # Estimate the end time:
    each.time = elapsed / done
    est.left = (total - done) * each.time
    cat(paste0('Left: ', round(est.left,1),'s\t'))
    fin.time = Sys.time() + est.left
    cat(paste0('Est. ', format(fin.time, "%H:%M:%S"),'\n'))
    return(curr)
}

# -------------
# Load results:
# -------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

fns = list.files(path='multiRegion/dereg/',pattern=paste0(prefix,'.nblmm_reg.*.major.', celltype, '.*.Rda'))
rdf = c()
rdf2 = c()
for (fn in fns){
    base = sub(".Rda","",sub(".*.nblmm_reg.","",fn))
    base = strsplit(base,"\\.")[[1]]
    load(paste0('multiRegion/dereg/',fn))
    # Run variables:
    regdf$path = base[1]
    regdf$region = base[2]
    regdf$subtype = base[6]
    regdf2$path = base[1]
    regdf2$region = base[2]
    regdf2$subtype = base[6]
    cat(paste0(base, collapse='\t'), '\n')
    # Effect var:
    path = base[1]
    pathstr = path
    if (path %in% c('nrad','cogdxad')){ pathstr = paste0(path,'AD') }
    leff = paste0('logFC_',pathstr)
    peff = paste0('p_',pathstr)
    # Add to all runs:
    cnames = c('gene',leff, peff, 'cpc','zpc', 'q_path',
               'log10q','region','subtype','path')
    regdf = regdf[,cnames]
    regdf2 = regdf2[,cnames]
    names(regdf)[2:3] = c('logFC_path','p_path')
    names(regdf2)[2:3] = c('logFC_path','p_path')
    rdf = rbind(rdf, regdf)
    rdf2 = rbind(rdf2, regdf2)
}
rm(regdf, regdf2)
gc()

rdf = rdf[rdf$region %in% regions,]
rdf2 = rdf2[rdf2$region %in% regions,]

# ------------------------
# Plot overall statistics:
# ------------------------
crdf = agg.rename(log10q ~ path + region + subtype, rdf, function(x){sum(x > 2)}, 'nsig')
crdf2 = agg.rename(log10q ~ path + region + subtype, rdf2, function(x){sum(x > 2)}, 'nsig')
crdf$type = 'LMM'
crdf2$type = 'Fixed'
crdf = rbind(crdf, crdf2)

scale = 1.5
w = length(unique(crdf$path)) * scale
h = 2 + length(unique(crdf$subtype)) * length(unique(crdf$region)) / 10 * scale
# h = 1 + length(unique(crdf$subtype)) * length(unique(crdf$region)) / 6 * scale
gp = ggplot(crdf, aes(subtype, nsig, fill=type)) + 
    facet_grid(region~path) + 
    geom_bar(stat='identity', position='dodge') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    labs(x=paste0(celltype, 'subtype'), y='Number of significant genes') +
    scale_fill_manual(values=c('Fixed'='grey80','LMM'='slateblue'), name='Model type:') + 
    theme_pubr() + coord_flip() + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(paste0(imgpref, 'subtype_nsig_counts_',celltype,'.png'), gp, dpi=400, units='in', width=w, height=h)
ggsave(paste0(imgpref, 'subtype_nsig_counts_',celltype,'.pdf'), gp, dpi=400, units='in', width=w, height=h)


if (celltype == 'Inh') {
    sum(cellmeta$cell_type_high_resolution == 'Inh PVALB HTR4' & cellmeta$region == 'PFC')
    idf = rdf[rdf$subtype == 'Inh_PVALB_HTR4' & rdf$path == 'nft',]
    head(idf, 20)
    head(idf[idf$logFC_path > 0,], 20)
}



# Volcano plots:
rdf$type = 'LMM'
rdf2$type = 'Fixed'
path = 'nrad'
subtype = 'Inh_CUX2_MSR1'
# subtype = 'Inh_L3-5_SST_MAFB'
# subtype = 'Inh_L6_SST_NPY'
for (subtype in c('Inh_CUX2_MSR1','Inh_L6_SST_NPY','Inh_L3-5_SST_MAFB')){
    for (path in c('nrad','cogdxad','nft')){
        cat(subtype,
        print(path)

        srdf = rdf[rdf$subtype == subtype & rdf$path == path,]
        srdf2 = rdf2[rdf2$subtype == subtype & rdf2$path == path,]
        subdf = rbind(srdf, srdf2)
        subdf$is.sig = subdf$log10q > 1.3
        srdf = srdf[order(srdf$log10q, decreasing=T),]
        srdf2 = srdf2[order(srdf2$log10q, decreasing=T),]

        ldf = rbind(head(srdf, 20), head(srdf2, 20))
        ldf = ldf[ldf$log10q > 1.5,]
        ldf$is.sig = ldf$log10q > 1.3

        if (nrow(subdf) > 0){

            scale = 1.75
            w = 1 + length(unique(subdf$region)) * scale
            h = 1 + 2 * scale
            gp = ggplot(subdf, aes(logFC_path, log10q, color=is.sig)) + 
                facet_grid(type~region, scales='free_y') + 
                geom_point(cex=.25) + 
                geom_hline(yintercept=2) +
                geom_hline(yintercept=1.3,lty='dashed') +
                geom_text_repel(data=ldf, aes(logFC_path, log10q, label=gene, color=is.sig), size=2, segment.size=.1) + 
                labs(x=paste0('log2FC of ', path), y='-log10q', title=paste0(celltype, ': ',subtype,' (',path, ')')) +
                scale_color_manual(values=c('grey70','indianred'), name='q < 0.01:') + 
                scale_fill_manual(values=c('Fixed'='grey80','LMM'='slateblue'), name='Model type:') + 
                theme_pubr()
            ggsave(paste0(imgpref, 'subtype_volcsummary_',celltype,'_', subtype,'_',path,'.png'), gp, dpi=400, units='in', width=w, height=h)
            ggsave(paste0(imgpref, 'subtype_volcsummary_',celltype,'_', subtype,'_', path, '.pdf'), gp, dpi=400, units='in', width=w, height=h)


            # Select axis to plot heatmaps of top genes on: (region x genes) or (subtype x genes)
            topgenes2 = head(unique(srdf2$gene), 50)
            swide2 = spread(srdf2[srdf2$gene %in% topgenes2, c('gene','logFC_path','region')], region, logFC_path)
            smat2 = as.matrix(swide2[,-1, drop=F])
            rownames(smat2) = swide2[,1]
            smat2[is.na(smat2)] = 0
            mx = max(abs(range(smat2)))
            rsmat2 = reord(smat2)

            scale = .75
            w = 1 + ncol(rsmat2) / 2 * scale
            h = 1 + nrow(rsmat2) / 5 * scale
            png(paste0(imgpref, 'subtype_heatmap_',celltype,'_', subtype,'_',path,'.png'), res=400, units='in', width=w, height=h)
            sp=0.1
            par(mar=c(sp,6,3,sp))
            image(t(rsmat2), zlim=c(-mx, mx), axes=F, useRaster=T, col=rev(colrb))
            mtext(paste0(celltype, ': ',subtype,' against ',path), side=3, line=1)
            text(x=seq(0,1, length.out=ncol(rsmat2)), y=parpos(2,-1.01), colnames(rsmat2), adj=.5, xpd=TRUE, cex=1)
            text(y=seq(0,1, length.out=nrow(rsmat2)), x=parpos(1,.01), rownames(rsmat2), adj=1, xpd=TRUE, cex=1)
            dev.off()

            rsmat2 = t(rsmat2)

            # Flip:
            scale = .75
            w = 1 + ncol(rsmat2) / 5 * scale
            h = 1 + nrow(rsmat2) / 3 * scale
            png(paste0(imgpref, 'subtype_heatmap_horiz_',celltype,'_', subtype,'_',path,'.png'), res=400, units='in', width=w, height=h)
            sp=0.1
            par(mar=c(5,2,1,sp))
            image(t(rsmat2), zlim=c(-mx, mx), axes=F, useRaster=T, col=rev(colrb))
            mtext(paste0(celltype, ': ',subtype,' against ',path), side=3, line=0)
            text(x=seq(0,1, length.out=ncol(rsmat2)), y=parpos(2,.01), colnames(rsmat2), adj=1, srt=90, xpd=TRUE, cex=1)
            text(y=seq(0,1, length.out=nrow(rsmat2)), x=parpos(1,.005), rownames(rsmat2), adj=1, xpd=TRUE, cex=1)
            dev.off()
        }
    }
}









print("Finished plotting regression results.")
