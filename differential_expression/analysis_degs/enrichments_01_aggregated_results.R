#!/usr/bin/R
# ----------------------------------------------------------
# Enrichments on aggregated differential expression results:
# Updated: 11/25/21
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(gprofiler2)
library(tidyr)
print(version)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
enrdir = paste0(sdbdir, 'dereg/enrichments/')
cmd = paste('mkdir -p', regdir, enrdir)
system(cmd)


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


# For each merged differential run, get functional enrichments:
# -------------------------------------------------------------
merged.ind = which(rundf$merged == 1)
for (i in merged.ind){
    prefstr = rundf$prefstr[i]
    aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
    enrfile = paste0(regdir, 'allmethods.', prefstr, '.merged.enrichments.tsv.gz')
    enrrda = paste0(regdir, 'allmethods.', prefstr, '.merged.enrichments.rda')

    if (file.exists(aggrda) & (!file.exists(enrrda))){
        load(aggrda)  # Loads aggdf, nsig
        cat(nsig, '\n')
        degs = list()
        degs$up = aggdf$gene[aggdf$col_nm == 2]
        degs$down = aggdf$gene[aggdf$col_nm == 1]

        # Remove empty sets:
        dl = c(lapply(degs, length))
        keep = names(dl)[dl > 0]
        degs = degs[keep]

        # Initialize (for if any of these fails, we can save NULL):
        gp2.nenr = NULL
        gp2.result = NULL
        gp2df = NULL
        if (length(keep) > 0){
            # Run functional enrichments with gprofiler2:
            sources = c("GO:CC","GO:BP","GO:MF","REAC","WP","KEGG","CORUM")
            gp2.result = gprofiler2::gost(degs, organism='hsapiens',
                                   ordered_query=FALSE, multi_query=TRUE,
                                   sources = sources)
            gp2df = gp2.result$result
            if (!is.null(gp2df)){
                pmat = t(as.matrix(data.frame(gp2df$p_values)))
                rownames(pmat) = NULL  # Format 
                pvals = c()
                gp2df$col = ''
                for (j in 1:length(keep)){
                    pstr = paste0('p_', keep[j])
                    gp2df[[pstr]] = pmat[,j]
                    pvals = c(pvals, pstr)
                    # Significant:
                    sind = gp2df[[pstr]] < 0.05
                    gp2df$col[sind] = paste0(gp2df$col[sind], keep[j])
                }
                gp2cols = c('term_id',pvals,'source','term_name', 'col')
                gp2df = gp2df[,gp2cols]
                # subdf = gp2df[gp2df$source %in% c('KEGG','REAC','WP','CORUM'),]
                gp2.nenr = table(gp2df[,c('col','source')])
                # Write out this aggregated set of results:
                write.table(gp2df, gzfile(enrfile), quote=F, row.names=F, col.names=T, sep="\t")
            } 
        }
        save(gp2df, gp2.result, gp2.nenr, file=enrrda)
    }
}
