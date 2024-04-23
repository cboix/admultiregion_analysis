#!/usr/bin/R
# -----------------------------------------------
# Plot statistics for cohort and choose datasets:
# Preliminary plots for the metadata figures
# -----------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(ggplot2)
library(ggpubr)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/metadata/')
imgpref = paste0(plotdir, 'meta_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)


# ----------------------------------
# Pivot table to look at individuals
# ----------------------------------
pwide = spread(aggregate(library_id ~ projid + region, metadata, length), 
               region, library_id)
pwide[is.na(pwide)] = 0
pmat = as.matrix(pwide[,-1])
rownames(pmat) = pwide[,1]

# Aggregate number:
ind.nregion = apply(pmat, 1, sum)
print(ind.nregion)
nregion = sapply(1:7, function(x){ sum(apply(pmat, 1, sum) >= x)})
names(nregion) = 1:7
print(nregion)
# Keep inidividuals with 5+ regions:
kept.individuals = sort(names(which(ind.nregion > 4)))
print(length(kept.individuals))
write.table(kept.individuals, 'Annotation/multiRegion_individuals.txt', quote=F, row.names=F, col.names=F, sep="\t")

# Kept libraries:
kept.libraries = as.character(metadata$library_id[metadata$projid %in% kept.individuals])
write.table(kept.libraries, 'Annotation/multiRegion_libraries.txt', quote=F, row.names=F, col.names=F, sep="\t")

# Kept rows:
kept.rnames = as.character(rownames(metadata)[metadata$projid %in% kept.individuals])
write.table(kept.rnames, 'Annotation/multiRegion_rows.txt', quote=F, row.names=F, col.names=F, sep="\t")


# ---------------------------------------
# Plot the breakdown of pathology levels:
# ---------------------------------------
metapath = unique(metadata[metadata$projid %in% kept.individuals, c('projid','cogdx','amyloid','tangles','braaksc', 'niareagansc','ceradsc')])
metapath = metapath[order(metapath$projid),]

# NOTE: We will divide AD, CTRL by the NIA + Reagan score:
aggregate(projid ~ niareagansc,metapath, length)
aggregate(projid ~ niareagansc,metapath, length)
aggregate(projid ~ cogdx,metapath, length)
aggregate(projid ~ braaksc,metapath, length)


# --------------------------------------------
# Plot the tangles, plaq_d, and plaq_n levels:
# --------------------------------------------
for (path in c('nft','plaq_d','plaq_n')){
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    submeta = unique(metadata[,c('projid','niareagansc', vars)])
    var.avg = apply(submeta[,vars], 2, mean)
    varlvls = names(sort(var.avg))
    mlong = gather(submeta, var, value, -projid, -niareagansc)
    mlong$var = factor(mlong$var, levels = varlvls)
    lbls = sub(paste0(path,"_"), "", varlvls)
    lbls[lbls == path] = 'avg'
    lbls = toupper(lbls)

    gplot = ggplot(mlong, aes(as.numeric(var), value, color=factor(niareagansc), fill=factor(projid))) + 
        facet_wrap(~niareagansc) + 
        scale_color_manual(values=colvals[['niareagansc']]) + 
        labs(y=paste(capitalize(path), 'levels'), x='Region of measurement') + 
        geom_point() + geom_line() + 
        scale_x_continuous(breaks=1:length(lbls), labels=lbls) + 
        theme_pubr() + theme(legend.position='none')
    ggsave(paste0(imgpref, path, '_vs_niareagansc_perprojid.png'), gplot, dpi=400, units='in', width=6, height=5)
}


