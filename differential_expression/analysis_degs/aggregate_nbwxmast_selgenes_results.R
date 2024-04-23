#!/usr/bin/R
# --------------------------------------
# Aggregate the DE results for selgenes:
# Updated: 02/24/22
# --------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
print(version)
options(width=175)

# Directories:
regdir = paste0(sdbdir, 'dereg/')
plotdir = paste0(imgdir, 'difftl/')
imgpref = paste0(plotdir, 'difftl_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)


# Get a list of all differential runs:
# -------------------------------------------------------
rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_selgenes_runlist.tsv'), header=T)
rundf = rbind(rundf, read.delim(paste0(sdbdir, 'nebula_wRUV_excitatory_subsets_runlist.tsv'), header=T))
rundf$prefstr = with(rundf, gsub("[()/]", "_", paste(celltype, subtype, region, path, geneset, sep="_")))

# Check: how many of the runs are completely done:
rundf$complete = sapply(rundf$prefstr, function(x){
                            length(list.files(path=regdir, pattern=paste0('.', x, '.rda'))) })
table(rundf$complete)
# Which aren't done:
head(rundf[rundf$complete < 3,])

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
                            length(list.files(path=regdir, pattern=paste0('allmethods.', x, '.merged.rda'))) })
table(rundf$merged)
head(rundf[rundf$merged == 0,], 50)


# Functions to load each of the different regression methods:
# -----------------------------------------------------------
process_nebula = function(prefstr, pcut=0.05){
    out = load(paste0(regdir, 'nebula_ruv.', prefstr, '.rda'))
    fulldf$col = 1 * (fulldf$padj < pcut) * (2 - 1 * (fulldf$logFC < 0))
    df = fulldf[,c('gene','logFC','p','padj','col','pc')]
    names(df) = c('gene','logFC_nb','p_nb','padj_nb','col_nb','pc')
    return(list(df, nsig))
}

process_mast = function(prefstr, pcut=0.05){
    out = load(paste0(regdir, 'mast.', prefstr, '.rda'))
    fulldf$col = 1 * (fulldf$fdr < pcut) * (2 - 1 * (fulldf$coef < 0))
    df = fulldf[,c('gene','coef','p','fdr','col','val')]
    names(df) = c('gene','coef_mast','p_mast','padj_mast','col_mast','pc')
    return(df)
}


# For each run, aggregate the results from the three methods:
# -----------------------------------------------------------
methodlist = c('nebula_ruv','mast')
complete.ind = which(rundf$complete >= 2)
for (i in complete.ind){
    prefstr = rundf$prefstr[i]
    aggfile = paste0(regdir, 'allmethods.', prefstr, '.merged.tsv.gz')
    aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
    aggnsig = paste0(regdir, 'allmethods.', prefstr,'.merged.nsig.tsv')

    # Check exists all methods for specific string:
    fnlist = paste0(regdir, methodlist[1:2], '.', prefstr, '.rda')
    if ((!file.exists(aggrda)) & (sum(file.exists(fnlist)) == 2)){
        # Load each of the methods + merge DEGlist (tested on same set of genes):
        ll = process_nebula(prefstr)
        ndf = ll[[1]]
        nsig.nb = ll[[2]]
        mdf = process_mast(prefstr)
        aggdf = merge(ndf, mdf)

        # Both MAST and Nebula agree:
        aggdf$col_nm = 0
        aggdf$col_nm[(aggdf$col_nb == 2) & (aggdf$col_mast == 2)] = 2
        aggdf$col_nm[(aggdf$col_nb == 1) & (aggdf$col_mast == 1)] = 1

        # Order based on MAST + Nebula:
        aggdf$log10p_nm = -(log10(aggdf$p_nb) + log10(aggdf$p_mast))
        aggdf = aggdf[order(-aggdf$log10p_nm),]
        aggdf = aggdf[order(aggdf$col_nm == 0),]

        # Write number of significant genes:
        ctvec = nsig.nb[,1:6,drop=F]
        nsig = cbind(ctvec, t(table(aggdf$col_nm)))
        write.table(nsig, aggnsig, quote=F, row.names=F, sep="\t")
        cat(nsig, '\n')
        cat('UP:\t', head(aggdf[aggdf$col_nm == 2,'gene'], 15),'\n')
        cat('DOWN:\t', head(aggdf[aggdf$col_nm == 1,'gene'], 15),'\n')

        # Write out this aggregated set of results:
        write.table(aggdf, gzfile(aggfile), quote=F, row.names=F, col.names=T, sep="\t")
        save(aggdf, nsig, file=aggrda)
    }
}


# Collate all of these results into one dataframe:
# ------------------------------------------------
for (geneset in unique(rundf$geneset)){
    fullaggrda = paste0(regdir, 'allmethods.allmajor.', geneset, '.merged.rda')
    selrundf = rundf[rundf$geneset %in% geneset, ]
    selrundf$setid = with(selrundf, paste0(celltype, '_', subtype))
    sets = unique(selrundf$setid)
    print(sets)
    kept.cols = c("gene","pc", "col_nm","log10p_nm","path", 'region',
                  "logFC_nb","p_nb","padj_nb","col_nb",
                  "coef_mast","p_mast","padj_mast","col_mast")

    setdflist = lapply(sets, function(x){})
    names(setdflist) = sets
    complete.ind = which(selrundf$complete >= 2)
    totnsigdf = NULL
    for (i in complete.ind){
        prefstr = selrundf$prefstr[i]
        setid = selrundf$setid[i]
        aggrda = paste0(regdir, 'allmethods.', prefstr, '.merged.rda')
        # Check exists all methods for specific string:
        if (file.exists(aggrda)){
            load(aggrda)
            cat(nsig, '\n')
            # Concatenate results:
            aggdf$path = nsig[1, 'path']
            aggdf$region = nsig[1, 'region']
            setdflist[[setid]] = rbind(setdflist[[setid]], aggdf[,kept.cols])
            # Pad nsig:
            nsigdf = data.frame(nsig)
            if (!("X1" %in% colnames(nsigdf))){ nsigdf$X1 = 0}
            if (!("X2" %in% colnames(nsigdf))){ nsigdf$X2 = 0}
            totnsigdf = rbind(totnsigdf, nsigdf)
        }
    }
    save(setdflist, totnsigdf, file=fullaggrda)
}

