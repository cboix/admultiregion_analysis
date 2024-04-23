#!/usr/bin/R
# ----------------------------------------------------------
# DEseq2 run on pseudobulk data for the multiregion project:
# Updated: 10/27/23
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(Matrix)
library(DESeq2)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
print(version)
options(width=170)

# Arguments:
# runfile=paste0(sdbdir, 'nebula_wRUV_joint_runlist.tsv')
args=(commandArgs(TRUE))
if (length(args)==0) {
    print("No arguments supplied: runfile task")
} else {        
    runfile = args[1]
    task = as.numeric(args[2])
}


# Functions for DE:
cbindir = paste0(sbindir, 'differential_expression/calculate_degs/')
source(paste0(cbindir, 'auxiliary_differential_functions_psbulk.R'))


# Read run table and set run options:
# -----------------------------------
rundf = read.delim(runfile, header=T, stringsAsFactors=FALSE)
celltype = as.character(rundf[task, 1])
subtype  = as.character(rundf[task, 2])
region   = as.character(rundf[task,3])
path     = as.character(rundf[task,4])

# Run interaction model for msex:
if (path == 'msex'){ int = TRUE } else { int = FALSE }


# Set the directories:
# --------------------
plotdir = paste0(imgdir, 'difftl/')
regdir = paste0(sdbdir, 'dereg/')
srdir = paste0(sdbdir, 'subtype_reg/')
imgpref = paste0(plotdir, 'difftl_psbulk_')
cmd = paste('mkdir -p', plotdir, regdir)
system(cmd)


# Set the run outputs: 
# ---------------------
prefstr = gsub("[()/]","_", paste(celltype, subtype, region, path, sep="_"))
if (int){ prefstr = paste0(prefstr, '_interact') }
prefstr = paste0('psbulk_DEseq2', '.', prefstr)
outpref = paste0(regdir, prefstr)
outtsv = paste0(outpref, '.tsv.gz')
outrds = paste0(outpref, '.rds')
nsigfile = paste0(outpref, '.nsig.tsv')
pcut = 0.05  # DEG thresholding
pctcut = 0.2  # Gene thresholding
print(outpref)


# Check if already computed:
if (!file.exists(outrds)){
    # Run differential expression using DESeq2 on psbulk:
    # ---------------------------------------------------
    # Load pseudobulk data:
    commandArgs = function(x){ c(celltype, subtype, region)}
    source(paste0(cbindir, 'load_psbulk_data.R'))

    # Subset dataset to genes found in pctcut (20%+) of samples:
    ps.data = subsetPsbulkMatrixForDE(ps.data, pctcut=pctcut, path=path)
    # Remove TH if pathology measurement:
    if (path %in% c('nft','plaq_n','plaq_d')){
        ps.data$meta = ps.data$meta[ps.data$meta$region != 'TH',]
        ps.data = harmonizePsbulkData(ps.data)
    }

    # Set some columns as factors:
    for (fvar in c('msex', 'region', 'cell_type_high_resolution')){
        ps.data$meta[[fvar]] = factor(ps.data$meta[[fvar]])
    }

    # Make the condition string, depdt on is.numeric, interaction
    if (is.numeric(ps.data$meta[[path]])){
        if (int){
            cond = paste0(path, '.nrad', c('AD', 'CTRL'))
        } else { cond = path }
    } else {
        lvls = levels(ps.data$meta[[path]])
        if (int){
            cond = paste0(path, lvls, '.nradAD')
        } else {
            cond = paste0(path, '_', lvls[2], '_vs_', lvls[1])
        }
    }

    # Make the formula:
    ps.data$meta$age_z = scale(ps.data$meta$age_death)
    ps.data$meta$pmi[is.na(ps.data$meta$pmi)] = median(ps.data$meta$pmi, na.rm=T)
    ps.data$meta$pmi_z = scale(ps.data$meta$pmi)
    covars = '+ msex + age_z + pmi_z'
    if (length(unique(ps.data$meta$region)) > 1){ covars = paste0(covars, '+ region') }
    if (celltype == subtype){ covars = paste0(covars, ' + cell_type_high_resolution') }
    if (int){
        form = asform(c('~', path, ':nrad', covars))
    } else {
        form = asform(c('~', path, covars))
    }
    print(form)

    # Make the DESeq2 object, run DESeq2:
    dds = DESeqDataSetFromMatrix(ps.data$mat, 
        colData=ps.data$meta, design=form)
    dds <- DESeq(dds)

    # Get and save results:
    if (length(cond) > 1){
        resdf = c()
        for (cc in cond){
            res1 = results(dds, name=cc)
            subdf <- data.frame(res1[order(res1$pvalue),])
            subdf$gene = rownames(subdf)
            subdf$col = ifelse(subdf$padj >= pcut, 'NS', 
                ifelse(subdf$log2FoldChange > 0, 'Up','Down'))
            subdf$cond = cc
            resdf = rbind(resdf, subdf)
        }
    } else {
        res1 = results(dds, name=cond)
        resdf <- data.frame(res1[order(res1$pvalue),])
        resdf$gene = rownames(resdf)
        resdf$col = ifelse(resdf$padj >= pcut, 'NS', 
            ifelse(resdf$log2FoldChange > 0, 'Up','Down'))
    }

    # Write out the results:
    # ----------------------
    # Write out number of significant genes:
    if (length(cond) > 1){
        nsig = table(resdf[,c('col','cond')])
        nsig = as.data.frame(nsig)
        nsig = spread(nsig, col, Freq)
    } else {
        nsig = table(resdf[,c('col')])
        nsig = as.data.frame(nsig)
        nsig = spread(nsig, Var1, Freq)
    }
    ctvec = data.frame(celltype=celltype, subtype=subtype, 
        region=region, path=path, 
        pctcut=pctcut, pcut=pcut)
    nsig = merge(ctvec, nsig)
    write.table(nsig, nsigfile, quote=F, row.names=F, sep="\t")
    print(nsig)

    # Write out the regression result dataframes:
    write.table(resdf, gzfile(outtsv), quote=F, row.names=F, sep="\t")
    saveRDS(resdf, file=outrds)

    # Plot results:
    if (length(cond) > 1){
        for (cc in cond){
            subdf = resdf[resdf$cond == cc,]
            plotVolcano(subdf, imgpref=imgpref, prefstr=paste0(prefstr, ".", cc))
        }
    } else {
        plotVolcano(resdf, imgpref=imgpref, prefstr=prefstr)
    }
} else {
    print("[STATUS] Pseudobulk DESeq2 regression output files already exist")
}
