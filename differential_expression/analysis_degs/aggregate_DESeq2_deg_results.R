#!/usr/bin/R
# --------------------------------------------
# Aggregate the DE results for DESeq2 results:
# Updated: 11/02/23
# --------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

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
rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_joint_runlist.tsv'), header=T)
rundf$prefstr = with(rundf, gsub("[()/]", "_", paste(celltype, subtype, region, path, sep="_")))

# Check: how many of the runs are completely done:
methodlist = c('psbulk_DEseq2', 'nebula_ruv', 'mast', 'wilcoxon')
cruns = t(sapply(rundf$prefstr, function(x){
        fnlist = paste0(regdir, methodlist, '.', x, '.tsv.gz')
        file.exists(fnlist)}))
colnames(cruns) = methodlist
rundf = cbind(rundf, as.data.frame(cruns, check.names=FALSE))
rundf$complete = apply(cruns[,1:3], 1, sum)
table(rundf$complete)
# Which aren't done:
head(rundf[rundf$complete < 3,])
rundf[rundf$complete < 3,]

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
    outfile = paste0(regdir, 'allmethods_ds2.', x, '.merged.rda')
    1 * file.exists(outfile) })
table(rundf$merged)
head(rundf[rundf$merged == 0,], 50)


# Functions to load each of the different regression methods:
# -----------------------------------------------------------
process_DESeq2 = function(prefstr, pcut=0.05){
    fulldf = readRDS(paste0(regdir, 'psbulk_DEseq2.', prefstr, '.rds'))
    print(head(fulldf))
    fulldf$col = 1 * (fulldf$padj < pcut) * (2 - 1 * (fulldf$log2FoldChange < 0))
    df = fulldf[,c('gene','log2FoldChange','pvalue','padj','col','baseMean')]
    names(df) = c('gene','log2FC_ds','p_ds','padj_ds','col_ds','baseMean')
    return(df)
}

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

process_wilcoxon = function(prefstr, pcut=0.05){
    fname = paste0(regdir, 'wilcoxon.', prefstr, '.rda')
    if (file.exists(fname)){
        out = load(fname)
        fulldf$col_cell = 1 * (fulldf$fdr < pcut) * (2 - 1 * (fulldf$logFC < 0))
        fulldf$col_ind = 1 * (fulldf$fdr_ind < pcut) * (2 - 1 * (fulldf$logFC < 0))
        df = fulldf[,c('gene','logFC','p_cell','p_ind','fdr','fdr_ind','col_cell','col_ind')]
        names(df) = c('gene','logFC_mast','p_cell_wx','p_ind_wx','padj_cell_wx','padj_ind_wx','col_cell_wx','col_ind_wx')
    } else { df = NULL }
    return(df)
}


# For each run, aggregate the results from the three methods:
# -----------------------------------------------------------
complete.ind = which(rundf$complete >= 3)
for (i in complete.ind){
    prefstr = rundf$prefstr[i]
    aggfile = paste0(regdir, 'allmethods_ds2.', prefstr, '.merged.tsv.gz')
    aggrda = paste0(regdir, 'allmethods_ds2.', prefstr, '.merged.rda')
    aggnsig = paste0(regdir, 'allmethods_ds2.', prefstr,'.merged.nsig.tsv')

    # Check exists all methods for specific string:
    fnlist = paste0(regdir, methodlist[1:3], '.', prefstr, '.tsv.gz')
    if ((!file.exists(aggrda)) & (sum(file.exists(fnlist)) == 3)){
        # Load each of the methods + merge DEGlist (tested on same set of genes):
        ll = process_nebula(prefstr)
        ndf = ll[[1]]
        nsig.nb = ll[[2]]
        mdf = process_mast(prefstr)
        aggdf = merge(ndf, mdf, all=TRUE)
        dedf = process_DESeq2(prefstr)
        aggdf = merge(aggdf, dedf, all=TRUE)


        # Load wilcoxon if computed:
        # NOTE: computing wilcoxon is not finished in some cases due to memory issues
        wdf = process_wilcoxon(prefstr)
        if (!is.null(wdf)){
            aggdf = merge(aggdf, wdf)

            nup = ((aggdf$col_nb == 2) + (aggdf$col_mast == 2) +
                   (aggdf$col_cell_wx == 2) + (aggdf$col_ind_wx == 2))
            ndown = ((aggdf$col_nb == 1) + (aggdf$col_mast == 1) +
                     (aggdf$col_cell_wx == 1) + (aggdf$col_ind_wx == 1))
            # At least 3 of the method call up:
            aggdf$col_joint = 0
            aggdf$col_joint[nup >= 3] = 2
            aggdf$col_joint[ndown >= 3] = 1
        }

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


