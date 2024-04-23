#!/usr/bin/R
# ---------------------------------------------------------------
# Plot volcano plots for extended data for aggregated DE results:
# Updated: 12/10/21
# ---------------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'aggenr_')
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



# Select runs to use (allregions, main ct):
# -----------------------------------------
region = 'allregions'
# for (region in reg.nomb){
full.file = paste0(regdir, 'aggregated_allres.', region, '.rda')
if (!file.exists(full.file)){
    # Select runs:
    selrundf = rundf[rundf$region == region,]
    selrundf$setid = with(selrundf, paste0(celltype, '_', subtype))
    sets = unique(selrundf$setid)
    print(sets)

    # Aggregate the DEGs and results across all runs:
    # -----------------------------------------------
    kept.cols = c("gene","pc", "col_nm","log10p_nm","path",
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
            setdflist[[setid]] = rbind(setdflist[[setid]], aggdf[,kept.cols])
            # Pad nsig:
            nsigdf = data.frame(nsig)
            if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
            if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
            totnsigdf = rbind(totnsigdf, nsigdf)
        }
    }
    save(setdflist, totnsigdf, file=full.file)
} else {
    load(full.file)
}


# Reduce to a specific AD variable + keep major only, plot a row:
# ---------------------------------------------------------
pathlist = unique(totnsigdf$path)
keep.sets = c("Mic_Immune_Mic", "Ast_Ast", "Opc_Opc", 
              "Oli_Oli", 'Inh_Inh', 'Exc_Exc')

degdf = NULL
labdf = NULL
ntop = 20
for (path in pathlist){
    for (set in keep.sets){
        setdf = setdflist[[set]]
        setdf = setdf[setdf$path == path, ]
        if (!is.null(setdf)){
            setdf$set = set
            degdf = rbind(degdf, setdf)
            labdf = rbind(labdf, 
                          head(setdf[setdf$col_nm == 2,], ntop),
                          head(setdf[setdf$col_nm == 1,], ntop))
        }
    }
}

sdf = aggregate(cbind(logFC_nb, coef_mast) ~ path + set, degdf, 'sd')
names(sdf) = c('path','set', 'nb_sd', 'mast_sd')
degdf = merge(degdf, sdf)
labdf = merge(labdf, sdf)
degdf$zdir = with(degdf,(logFC_nb / nb_sd + coef_mast / mast_sd) / 2)
labdf$zdir = with(labdf,(logFC_nb / nb_sd + coef_mast / mast_sd) / 2)

degdf = degdf[order(degdf$col_nm),]

degcols = c('0'='grey90','1'=col.paired[1],'2'=col.paired[5])
gp = ggplot(degdf, aes(zdir, log10p_nm, color=factor(col_nm))) + 
    facet_grid(set ~ path, scales="free_x") + 
    scale_color_manual(values=degcols) +
    geom_point(cex=.25) +
    scale_y_continuous(expand=c(0,0)) + 
    geom_vline(xintercept=0) + 
    geom_text_repel(data=labdf, aes(zdir, log10p_nm, label=gene, color=factor(col_nm)), cex=2, max.overlaps=100) + 
    theme_pubr() + theme(legend.position = 'none')

gp2 = rasterize(gp, layers='Point', dpi=450)

pltprefix = paste0(imgpref, 'volcano_majorct_DEgrid_', region)
w = 14; h = 9
ggsave(paste0(pltprefix, '.png'), gp2, dpi=400, units='in', width=w, height=h)
ggsave(paste0(pltprefix, '.pdf'), gp2, dpi=400, units='in', width=w, height=h)




# Plot top up/down DEGs as a volcano plot:
# ----------------------------------------
path = 'cogdxad'
for (path in pathlist){
    subdf = topdegdf[topdegdf$path == path,]
    subdf = subdf[order(subdf$col, decreasing=F),]
    subdf$log10p_nm[subdf$log10p_nm > 75] = 75
    degcols = c('0'='grey90','1'=col.paired[1],'2'=col.paired[5])
    labdf = rbind(head(subdf[subdf$col == 2,],35),
                  head(subdf[subdf$col == 1,],10))

    gp = ggplot(subdf, aes(zdir, log10p_nm, color=factor(col))) + 
        scale_color_manual(values=degcols) +
        geom_point(cex=.25) +
        scale_y_continuous(expand=c(0,0)) + 
        geom_vline(xintercept=0) + 
        geom_text_repel(data=labdf, aes(zdir, log10p_nm, label=gene, color=factor(col)), cex=2, max.overlaps=100) + 
        theme_pubr() + theme(legend.position = 'none')
    gp2 = rasterize(gp, layers='Point', dpi=450)

    pltprefix = paste0(imgpref, 'volcano_consistent_degs_', path)
    w = 4; h = 4
    ggsave(paste0(pltprefix, '.png'), gp2, dpi=400, units='in', width=w, height=h)
    ggsave(paste0(pltprefix, '.pdf'), gp2, dpi=400, units='in', width=w, height=h)
}


