#!/usr/bin/R
# ---------------------------------------------
# Print the DEGs that drive module enrichments:
# Updated 12/11/2021
# ---------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}


library(tidyr)
library(viridis)

library(ComplexHeatmap)
library(circlize)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = paste0(plotdir, 'repr_')
cmd = paste('mkdir -p', plotdir, moddir)
system(cmd)


# Set the run arguments:
# ----------------------
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: celltype graph_id subtype modsuff region")
} else {
    runset = args[1]
    graph_id = args[2]
}


# Load in and process data (saves to matrices):
# ---------------------------------------------
commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE, FALSE)}
source(paste0(sbindir, 'modules/load_modules_degenr.R'))


# Subset to module x path combinations that are up and significant:
# -----------------------------------------------------------------
subdf = statsdf[(statsdf$p.adj < 0.05),] #  & (statsdf$key == 'Up'),]
subdf = subdf[order(subdf$p),]
dedf = dedf[!is.na(dedf$gset),]

for (i in 1:nrow(subdf)){
    module = subdf$module[i]
    path = subdf$dkey[i]
    cat(module, '\t', path, '\n')
    subdedf = dedf[(dedf$dkey == path) & (dedf$module == module),]
    updedf = subdedf[subdedf$gset == 'Up',]
    dwdedf = subdedf[subdedf$gset == 'Down',]
    cat(head(updedf$gene, 15), '\n')
    cat(head(dwdedf$gene, 15), '\n')
    updedf = updedf[!is.na(updedf$gene),]
    if (nrow(updedf) > 0){
        sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
        gp2.result = gprofiler2::gost(updedf$gene, organism='hsapiens',
                                      ordered_query=FALSE, multi_query=FALSE,
                                      sources = sources)
        gpdf = gp2.result$result
        if (!is.null(gpdf)){
            gpdf = gpdf[order(gpdf$p_value),]
            gpdf = gpdf[gpdf$term_size < 500,]
            if (nrow(gpdf) > 0){
                gpdf$pstr = with(gpdf, paste0(term_name, ' (', sprintf('%0.1e', p_value),')'))
                out = paste(head(gpdf$pstr,5), collapse='\n')
                cat(out)
            }
        }
        cat("\n\n")
    }
}


used.modules = unique(subdf$module)
used.modules = c(1, 2,3, 5, 8 , 15)
nmdf = dedf[(dedf$dkey == 'cogdxad') & !(dedf$module %in% used.modules),]
cat(head(nmdf$gene[nmdf$gset == 'Up'], 20), '\n')
cat(head(nmdf$gene[nmdf$gset == 'Down'], 20), '\n')



print(subdedf[order(subdedf$gene),c('gene','gset')])




# Plot module attributions:
library(ggpubr)
gp = ggplot(dedf, aes(gset, fill=factor(module))) + 
    facet_wrap(~dkey) + 
    scale_fill_manual(values=snap.cols) + 
    # geom_bar(stat='count', position='fill') + 
    geom_bar(stat='count') + 
    coord_flip() + 
    theme_pubr()
pltprefix = paste0(imgpref, 'module_attributions_', fullpref)
h = 10; w = 10
ggsave(paste0(pltprefix, '.png'), gp, dpi=450, units='in', height=h, width=w)



