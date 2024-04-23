#!/usr/bin/R
# --------------------------------------------------
# Aggregate the DE results across different methods:
# Methods: Nebula, MAST, Wilcoxon (cell + ind-level)
# Updated: 11/15/23
# --------------------------------------------------
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
# rundf = read.delim(paste0(sdbdir, 'nebula_wRUV_joint_runlist.tsv'), header=T)
rundf = read.delim(paste0(sdbdir, 'DEG_multiRegion_SI_ACE_runlist.tsv'), header=T)
rundf$prefstr = with(rundf, gsub("[()/]", "_", paste(celltype, subtype, region, path, sep="_")))
rundf = rundf[rundf$path == 'msex',]

# Check: how many of the runs are completely done:
methodlist = c('nebula_ruv', 'mast', 'wilcoxon')
cruns = t(sapply(rundf$prefstr, function(x){
        fnlist = paste0(regdir, methodlist, '.', x, '.tsv.gz')
        file.exists(fnlist)}))
colnames(cruns) = methodlist
rundf = cbind(rundf, as.data.frame(cruns, check.names=FALSE))
rundf$complete = apply(cruns[,1:3], 1, sum)
table(rundf$complete)
# Which aren't done:
head(rundf[rundf$complete < 2,])

# Which have final merged outputs:
rundf$merged = sapply(rundf$prefstr, function(x){
    outfile = paste0(regdir, 'allmethods.', x, '.merged.rda')
    1 * file.exists(outfile) })
table(rundf[,c('merged','complete')])
# head(rundf[rundf$merged == 0,], 50)


# Functions to load each of the different regression methods:
# -----------------------------------------------------------
process_nebula = function(prefstr, pcut=0.05){
    out = load(paste0(regdir, 'nebula_ruv.', prefstr, '.rda'))
    if ('p' %in% colnames(fulldf)){
        fulldf$col = 1 * (fulldf$padj < pcut) * (2 - 1 * (fulldf$logFC < 0))
        df = fulldf[,c('gene','logFC','p','padj','col','pc')]
        names(df) = c('gene','logFC_nb','p_nb','padj_nb','col_nb','pc')
    } else {
        # Interaction model:
        peff = colnames(fulldf)[grep('^p_', colnames(fulldf))]
        pathstr = sub("^p_", "", peff)
        leff = paste0('logFC_', pathstr)
        ceff = paste0('col_', pathstr)
        cvars = c(leff, peff, ceff)
        df = fulldf[,c('gene', cvars ,'pc')]
        names(df) = c('gene', paste0(cvars, "_nb"), 'pc')
    }
    return(list(df, nsig))
}

process_mast = function(prefstr, pcut=0.05){
    out = load(paste0(regdir, 'mast.', prefstr, '.rda'))
    if ('p' %in% colnames(fulldf)){
        fulldf$col = 1 * (fulldf$fdr < pcut) * (2 - 1 * (fulldf$coef < 0))
        df = fulldf[,c('gene','coef','p','fdr','col','val')]
        names(df) = c('gene','coef_mast','p_mast','padj_mast','col_mast','pc')
    } else {
        # Interaction model:
        peff = colnames(fulldf)[grep('^p_', colnames(fulldf))]
        pathstr = sub("^p_", "", peff)
        leff = paste0('coef_', pathstr)
        ceff = paste0('col_', pathstr)
        cvars = c(leff, peff, ceff)
        df = fulldf[,c('gene', cvars ,'val')]
        names(df) = c('gene', paste0(cvars, "_mast"), 'pc')
    }
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
methodlist = c('nebula_ruv','mast', 'wilcoxon')
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
        if ('p_nb' %in% colnames(aggdf)){
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
        } else {
            peff = colnames(aggdf)[grep('^p_.*_nb$', colnames(aggdf))]
            pathstr = sub("_nb$", "", sub("^p_","", peff))
            nsig = c()
            for (pstr in pathstr){
                cvars = paste0("col_", pstr, "_", c("nb", "mast", "nm"))
                aggdf[[cvars[3]]] = 0
                aggdf[[cvars[3]]][(aggdf[[cvars[1]]] == 2) & (aggdf[[cvars[2]]] == 2)] = 2
                aggdf[[cvars[3]]][(aggdf[[cvars[1]]] == 1) & (aggdf[[cvars[2]]] == 1)] = 1

                # Order:
                pvars = paste0("p_", pstr, "_", c("nb", "mast"))
                lp = paste0("log10p_", pstr, "_nm")
                aggdf[[lp]] = -(log10(aggdf[[pvars[1]]]) + log10(aggdf[[pvars[2]]]))
                aggdf = aggdf[order(-aggdf[[lp]]),]
                aggdf = aggdf[order(aggdf[[cvars[3]]]== 0),]

                nsdf = table(aggdf[[cvars[3]]])
                nsdf = as.data.frame(nsdf)
                nsdf = spread(nsdf, Var1, Freq)
                nsdf$eff = pstr
                if (!('1' %in% colnames(nsdf))){ nsdf[['1']] = 0 }
                if (!('2' %in% colnames(nsdf))){ nsdf[['2']] = 0 }
                cn = c('0','1','2','eff')
                nsig = rbind(nsig[,cn, drop=F], nsdf[,cn, drop=F])
            }
            ctvec = nsig.nb[1,1:6,drop=F]
            nsig = merge(ctvec, nsig)
        }

        write.table(nsig, aggnsig, quote=F, row.names=F, sep="\t")
        if (nrow(nsig) == 1){
            cat(nsig, '\n')
            cat('UP:\t', head(aggdf[aggdf$col_nm == 2,'gene'], 15),'\n')
            cat('DOWN:\t', head(aggdf[aggdf$col_nm == 1,'gene'], 15),'\n')
        } else {
            print(nsig)
        } 

        # Write out this aggregated set of results:
        write.table(aggdf, gzfile(aggfile), quote=F, row.names=F, col.names=T, sep="\t")
        save(aggdf, nsig, file=aggrda)
    }
}


