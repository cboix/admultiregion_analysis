#!/usr/bin/R
# ----------------------------------------------------------
# Flag modules that may be due to contamination/ambient RNA:
# Updated 12/01/2021
# ----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
if (!exists("cellmeta")){
    source(paste0(sbindir, 'load_metadata.R'))
}
options(width=250)

source(paste0(sbindir, 'auxiliary_pseudobulk_loading_fns.R'))

library(tidyr)

# Directories:
moddir = paste0(sdbdir, 'modules/')
plotdir = paste0(imgdir, 'modules/')
imgpref = plotdir
cmd = paste('mkdir -p', plotdir, moddir, crossdir)
system(cmd)


# Set the run arguments:
# ----------------------
graph_id = 'boot'
runlist = c('Ast','Oli','Inh','Opc','Mic_Immune','Vasc_Epithelia',
            'HCneurons', 'ECneurons', 'THneurons', 'CTXneurons')

agg.psbulk.rda = paste0(moddir, 'aggregated_psbulk_avg_matrices.rda')
if (!file.exists(agg.psbulk.rda)){
    corelist = list()
    assignlist = list()
    subtypelist = list()
    mmapdf = NULL
    ps.full = list('mat'=NULL, 'meta'=NULL)
    for (runset in runlist){
        # Load in and process the modules data:
        commandArgs <- function(trailingOnly=TRUE){c(runset, graph_id, TRUE)}
        source(paste0(sbindir, 'modules/load_modules_degenr.R'))

        # Save direct assignments:
        coremap = nodedf$leiden
        names(coremap) = nodedf$gene
        corelist[[runset]] = coremap

        # Save all assignments: 
        assignlist[[runset]] = genemap
        subtypelist[[runset]] = subtypes
        mmap$runset = runset
        mmapdf = rbind(mmapdf, mmap)

        # Load in the reduced pseudobulk data (averaging individuals):
        psdata.rda = paste0(srdir, 'pseudobulk_data_', runset, '_noind.rda')
        if (!file.exists(psdata.rda)){
            ps.data = load_pseudobulk_dataset(celltype, subtypes, region.set, byind=FALSE)
            save(ps.data, file=psdata.rda)
        } else { load(psdata.rda) }


        # Merge with full pseudobulk dataset:
        if (!is.null(ps.data$mat)){ 
            ps.full$mat = cbind(ps.full$mat, ps.data$mat)
            ps.full$meta = rbind(ps.full$meta, ps.data$meta)
        }
    }
    save(ps.full, corelist, assignlist, subtypelist, mmapdf, file=agg.psbulk.rda)
}

subtypelist[['Mic_Immune']] = c("Mic DUSP1", "Mic MKI67", "Mic P2RY12", 
                                "Mic TPT1", "T cells", "CAMs")
subtypelist[['Oli']] = c("Oli OPALIN", "Oli RASGRF1")
subtypelist[['Opc']] = c("Opc CEP112", "Opc GPC5", "Opc DOCK5")
mmapdf$runmod = paste0(mmapdf$runset, '_', mmapdf$module)

# Make a table mapping all of the genes to modules:
# -------------------------------------------------
df.from.mapping = function(xmap){
    uq = sort(unique(xmap))
    tform = make.tform(xmap, u=uq, norm=FALSE)
    df = data.frame(tform)
    colnames(df) = uq
    df$gene = rownames(df)
    return(df)
}

dflist = sapply(names(corelist), function(x){
                    df = df.from.mapping(corelist[[x]])
                    colnames(df)[-ncol(df)] = paste0(x, '_', colnames(df)[-ncol(df)])
                    return(df) })
fulldf = dflist[[1]]
for (df in dflist){ fulldf = merge(fulldf, df, all=TRUE)}
fulldf[is.na(fulldf)] = 0

# Check recovery of genes:
nrow(fulldf) == length(unique(unlist(lapply(corelist, names))))
mapmat = as.matrix(fulldf[,-1])
rownames(mapmat) = fulldf$gene


# Score all modules in all cell types, turn into dataframe:
# ---------------------------------------------------------
tform = sweep(mapmat, 2, apply(mapmat, 2, sum), '/')
scmat = t(tform) %*% ps.full$mat[rownames(tform),]

# Turn into long dataframe:
scdf = data.frame(as.matrix(scmat), check.names=FALSE)
scdf$runmod = rownames(scdf)
scdf = gather(scdf, ptype, score, -runmod)
scdf = merge(scdf, ps.full$meta, all.x=TRUE)

# Mark if the subtype is in runset or not:
scdf$runset = sapply(scdf$runmod, function(x){sub("_[0-9]*$","",x)})
scdf$module = sapply(scdf$runmod, function(x){sub(".*_","",x)})
# TODO: Fix for Mic_Immune
scdf$inrun = sapply(1:nrow(scdf), function(i){
                        st = scdf$cell_type_high_resolution[i]
                        runset = scdf$runset[i]
                        stlist = subtypelist[[runset]]
                        st %in% stlist })

# Remove low cell counts:
scdf = scdf[scdf$ncell >= 50,]


# Get the top module and score in + out of runset:
# ------------------------------------------------
topdf = aggregate(score ~ runmod + inrun, scdf[scdf$ncell >= 1000,], max)
topdf = merge(topdf, scdf)

# Show in a paired manner:
rdf = data.frame(runmod=unique(topdf$runmod))
for (attr in c('score','ncell','ptype')){
    wide = spread(topdf[,c(attr,'runmod','inrun')], 'inrun', attr)
    names(wide)[2:3] = paste0(attr,'_',names(wide)[2:3])
    rdf = merge(rdf, wide)
}
rdf$ratio = rdf$score_FALSE / rdf$score_TRUE
rdf = merge(rdf, mmapdf[,c('runmod','mname')])

rdf = rdf[order(rdf$ratio, decreasing=T),]
head(rdf, 50)

# To flag, want to see high ncell + high score + low jaccard

# TODO: score it should be:

# Calculate and plot module overlap:
# ----------------------------------
library(proxy)
dt = dist(t(mapmat), 'jaccard')
dt = as.matrix(dt)
diag(dt) = NA
match = apply(dt, 1, min, na.rm=T)

rdf$jacc = match[rdf$runmod]
head(rdf[rdf$jacc < 0.8,], 50)

rdf[grep("^Ast_",rdf$runmod),]

# Score all modules by in/out of appropriate cell type
# ----------------------------------------------------


# Score modules by how many genes are higher expr than in intended celltype:
# --------------------------------------------------------------------------


tform = sweep(mapmat, 2, apply(mapmat, 2, sum), '/')
scmat = t(tform) %*% ps.full$mat[rownames(tform),]

# Turn into long dataframe:
scdf = data.frame(as.matrix(scmat), check.names=FALSE)
scdf$runmod = rownames(scdf)
scdf = gather(scdf, ptype, score, -runmod)
scdf = merge(scdf, ps.full$meta, all.x=TRUE)

# Mark if the subtype is in runset or not:
scdf$runset = sapply(scdf$runmod, function(x){sub("_[0-9]*$","",x)})
scdf$module = sapply(scdf$runmod, function(x){sub(".*_","",x)})
# TODO: Fix for Mic_Immune
scdf$inrun = sapply(1:nrow(scdf), function(i){
                        st = scdf$cell_type_high_resolution[i]
                        runset = scdf$runset[i]
                        stlist = subtypelist[[runset]]
                        st %in% stlist })

# Remove low cell counts:
scdf = scdf[scdf$ncell >= 50,]


# TODO: Count number and score number



# Overall score
# Number of genes higher expr. in other CT
# Number of genes that are markers for other CT
# Flag genes that are in many modules



# Flag leniently, visualize, flag by hand:
# ----------------------------------------



