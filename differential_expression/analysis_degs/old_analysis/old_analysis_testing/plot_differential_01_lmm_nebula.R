#!/usr/bin/R
# -----------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Using Liang's Nebula package to run the fast NBLMM
# Basic overall + per-region interactions, for relevant subtypes
# Updated: 10/14/2020
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(rhdf5)
library(nebula)
library(qvalue)
library(Matrix)

celltype = 'Mic_Immune'
chunksize=4
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need celltype")
} else {        
    celltype = args[1]
    if (length(args) > 1){
        chunk=as.integer(args[2])
    } else { chunk=NULL }
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

# ------------------------
# Load pathology measures:
# ------------------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)

# Get the pathology mapped to each region:
regmap = c('AG','HC','PFC','MT','EC')
names(regmap) = c('ag','hip','mf','mt','ec')
pqdf = NULL
for (path in pathlist){
    vars = colnames(metadata)[grep(path, colnames(metadata))]
    vars = vars[vars != path]
    submeta = unique(metadata[,c('projid','region', vars, 'rind')])
    slong = gather(submeta, path, value, -projid, -region, -rind)
    slong$path.region = regmap[sub(".*_","", slong$path)]
    slong = slong[slong$region == slong$path.region,]
    rownames(slong) = slong$rind
    sdf = slong[,c('rind','value','region')]
    names(sdf)[2] = path
    if (is.null(pqdf)){
        pqdf = sdf 
    } else {
        pqdf = merge(pqdf, sdf)
    }
}
# pqdf = merge(pqdf, metadata[,c('rind','projid')])
# write.table(pqdf, file=paste0(datadir, 'region_pathology_scores.tsv'), quote=F, row.names=F, sep="\t")

# --------------------------------
# Load in the barcodes, pathology:
# --------------------------------
rawpref = 'all_brain_regions_filt_preprocessed_scanpy'
if (dbdir == '~/data/DEVTRAJ/db/') {
    matdir = paste0('/broad/compbio_ce/cboix/multiRegion/matrices/')
} else {
    matdir = paste0(datadir, 'matrices/')
}
h5file = paste0(matdir, rawpref, '.majorcelltype.', celltype, '.hdf5')

# Load margin (for offset term):
margfile = paste0(matdir, rawpref, '_fullmatrix_margin.tsv.gz')
marg = read.delim(gzfile(margfile), header=F)
names(marg) = 'count'
mbcs = scan(paste0(datadir, prefix,'.barcodes.tsv.gz'), 'c', quiet=T)
marg$barcode = mbcs
rownames(marg) = marg$barcode

# Extract metadata from hdf5:
h5f = H5Fopen(h5file)
genes = h5f$genes
barcodes = h5f$barcodes
H5Fclose(h5f)
ngenes = length(genes)

# ----------------------
# Make the model matrix:
# ----------------------
rownames(cellmeta) = cellmeta$barcode
submeta = cellmeta[barcodes,]
print(nrow(submeta[submeta$region != 'TH',]))
pathdf = merge(submeta, pqdf)
print(nrow(pathdf))
pathdf = merge(pathdf, metadata[,c('projid','rind','age_death','msex','pmi', 'Apoe_e4', 'cogdx')])
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode
pathdf$cogdxad = 'NCI'
pathdf$cogdxad[pathdf$cogdx %in% c(4,5)] = 'AD'

# Split by ct:
split.var = 'hcelltype'
if (celltype %in% c('Exc','Inh', 'Vasc_Epithelia')){ 
    split.var = 'cell_type_high_resolution' 
} else if (celltype == 'Mic_Immune'){
    split.var = 'hcluster'
}
subtypes = unique(pathdf[[split.var]])
if (!is.null(chunk)){
    ind = ((chunk -1)* (chunksize) + 1):min(c(chunk * chunksize,length(subtypes)))
    subtypes = subtypes[ind]
}
print(subtypes)

# -----------------------------------------------------------
# Run the regression on a set of data:
# Req: count matrix, random effects, model matrix, and offset
# -----------------------------------------------------------
# Run on each subtype:
for (subtype in subtypes){
    print(paste("[STATUS] Running on", subtype))
    sub.pathdf = pathdf[pathdf[[split.var]] == subtype,] 
    print(paste("Number of cells:", nrow(sub.pathdf)))
    ststr = gsub("/","_",gsub(" ","_", subtype))
    offset = marg[sub.pathdf$barcode, 'count']
    for (path in c(pathlist,'cogdxad')) {
        # TODO: For now, run without interaction terms, later add. 
        # mdx = model.matrix(asform(c('~',path, '* region + age_death + msex + pmi')), data=sub.pathdf)
        mdx = model.matrix(asform(c('~',path, '+ region + age_death + msex + pmi')), data=sub.pathdf)
        # Output files:
        regfile = paste0(regdir, prefix, '.nblmm_reg.', path, '.major.', celltype, '.minor.', ststr, '.Rda')
        regtsv  = paste0(regdir, prefix, '.nblmm_reg.', path, '.major.', celltype, '.minor.', ststr, '.tsv.gz')
        if (!file.exists(regfile)){
            print(path)
            chunksize = 500
            nchunk = floor(ngenes / chunksize) + 1
            regdf = c()
            # Subset of locations to extract: 
            bind = match(sub.pathdf$barcode, barcodes)
            subbcs = barcodes[bind] 
            t0 = proc.time()
            for (i in 1:nchunk){
                t1 = proc.time()
                print(i)
                ind = (1 + (i-1) * chunksize):min(c(i * chunksize, ngenes)) 
                # Open handle, extract genes we care about and close:
                h5f = H5Fopen(h5file)
                genes = h5f$genes
                bcs = h5f$barcodes
                h5d = h5f&"matrix"
                mat = t(h5d[ind,bind])
                H5Dclose(h5d)
                H5Fclose(h5f)
                colnames(mat) = genes[ind]
                rownames(mat) = bcs[bind]
                # Order as model matrix + transpose matrix:
                mat = mat[sub.pathdf$barcode,]
                mat = t(Matrix(mat))
                gcout = gc()
                # NOTE: Works fine as both PMM and as NBLMM, choose NBLMM - more conservative
                re = nebula(mat, as.character(sub.pathdf$projid), pred=mdx, offset=offset) 
                rdf = re$summary
                if (nrow(rdf) > 0){ regdf = rbind(regdf, rdf) }
                # Timing:
                t2 = (proc.time() - t1)[3]
                ttot = (proc.time() - t0)[3]
                names(t2) = NULL
                names(ttot) = NULL
                cat(paste0(round(t2,1),'s\t'))
                cat(paste0(round(ttot,1),'s\n'))
            }
            regdf$path = path
            if (path == 'cogdxad'){
                pval = paste0('p_', path,'NCI')
            } else {
                pval = paste0('p_',path)
            }
            regdf$q_path = qvalue(regdf[[pval]])$q
            regdf$log10q = -log10(regdf$q_path)
            regdf = regdf[order(regdf$log10q, decreasing=T),]
            # Save:
            save(regdf, file=regfile)
            write.table(regdf, file=gzfile(regtsv), quote=F, row.names=F, sep="\t")
        } else {
            load(regfile)
        }
        if (path == 'cogdxad'){
            eff = paste0('logFC_',path, 'NCI')
        } else {
            eff = paste0('logFC_',path)
        }
        print(head(regdf[regdf[[eff]] > 0,c(eff,'log10q','gene')], 50))
        print(head(regdf[,c(eff,'log10q','gene')], 50))
    }
}


# If local, plot some of these:
if (dbdir != '~/data/DEVTRAJ/db/') {
    library(ggplot2)
    library(ggpubr)
    library(ggrepel)

    for (subtype in subtypes){
        print(paste("[STATUS] Running on", subtype))
        sub.pathdf = pathdf[pathdf[[split.var]] == subtype,] 
        print(paste("Number of cells:", nrow(sub.pathdf)))
        cellstr = gsub("/","_",gsub(" ","_", celltype))
        ststr = gsub("/","_",gsub(" ","_", subtype))
        offset = marg[sub.pathdf$barcode, 'count']
        for (path in pathlist) {
            # Output files:
            regfile = paste0(regdir, prefix, '.nblmm_reg.', path, '.major.', celltype, '.minor.', ststr, '.Rda')
            load(regfile)
            topgenes = c(head(regdf$gene, 25), 'APOE')
            eff = paste0('logFC_',path)
            # print(head(regdf[regdf[[eff]] > 0,c(eff,'log10q','gene')], 50))
            print(head(regdf[,c(eff,'log10q','gene')], 25))
            regdf$signif = regdf$log10q > 2
            regdf$rank = 1:nrow(regdf)

            gp = ggplot(regdf, aes_string(eff, 'log10q', color='signif')) + 
                geom_point() + theme_pubr() + 
                scale_color_manual(values=c('black','red')) + 
                geom_text_repel(data=regdf[regdf$signif & regdf$rank <= 50,], aes_string(eff, 'log10q',label='gene'), color='grey50', size=3) + 
                labs(x=paste('logFC on',path), y='log10 q-value') + 
                geom_vline(xintercept=0, lty='dashed') + 
                geom_hline(yintercept=2, lty='dotted') +
                theme(legend.position='none')
            ggsave(paste0(imgpref, 'volcano_',cellstr,'_sub_',ststr,'_',path,'_top50.png'),gp, units='in', dpi=450, width=6, height=9)

            # Load in the data matrix:
            bind = match(sub.pathdf$barcode, barcodes)
            ind = match(topgenes, genes)
            h5f = H5Fopen(h5file)
            genes = h5f$genes
            bcs = h5f$barcodes
            h5d = h5f&"matrix"
            mat = t(h5d[ind,])
            H5Dclose(h5d)
            H5Fclose(h5f)
            colnames(mat) = genes[ind]
            rownames(mat) = bcs
            # Order as model matrix + transpose matrix:
            mat = mat[sub.pathdf$barcode,]
            mat = t(Matrix(mat))
            gcout = gc()


            mt = data.frame(as.matrix(mat))
            mt$gene = rownames(mat)
            mtdf = gather(mt, barcode, rcount, -gene)
            mtdf$barcode = sub('\\.','-', mtdf$barcode)
            medmarg = median(marg$count)
            rownames(marg) = marg$barcode
            mtdf$count = marg[mtdf$barcode, 'count']
            mtdf$region = pathdf[mtdf$barcode, 'region']
            mtdf[[path]] = pathdf[mtdf$barcode, path]
            mtdf$norm = with(mtdf, log(rcount / count * medmarg + 1))
            mtdf$ptsplit = mtdf[[path]] > 2
            # TODO: Normalize --> by margin

            ggplot(mtdf[mtdf$gene == 'APOE',], aes_string(path, 'norm')) +
                facet_wrap(~region) + 
                geom_smooth(method='lm') + 
                geom_point(alpha=.1, cex=.5) + 
                # scale_y_log10() + 
                theme_pubr()


            ggplot(mtdf[mtdf$region == 'EC',], aes(gene, rcount, fill=factor(ptsplit))) +
                # facet_wrap(~region) + 
                geom_violin() + 
                scale_y_log10() + 
                theme_pubr()

        }
    }

}

print("Finished running regressions.")
