#!/usr/bin/R
# -----------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Using Liang's Nebula package to run the fast NBLMM
# Basic overall + per-region interactions, for relevant subtypes
# Updated: 10/28/2020
# 
# Compare multiple models (in a per-region sense)
# Compare if PLP1 in Mic remains or not?
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
pathdf = merge(pathdf, metadata[,c('projid','rind','age_death','msex','pmi', 'Apoe_e4', 'cogdx', 'niareagansc')])
pathdf = pathdf[order(pathdf$projid),]
rownames(pathdf) = pathdf$barcode
pathdf$cogdxad = 'NCI'
pathdf$cogdxad[pathdf$cogdx %in% c(4,5)] = 'AD'
pathdf$nrad = 'CTRL'
pathdf$nrad[pathdf$niareagansc %in% c(1,2)] = 'AD'

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

cellstr = gsub("/","_",gsub(" ","_", celltype))
path = 'nrad'

# -----------------------------------------------------------
# Run the regression on a set of data:
# Req: count matrix, random effects, model matrix, and offset
# -----------------------------------------------------------
# Run on each subtype:
for (subtype in subtypes){
    print(paste("[STATUS] Running on", subtype))
    sub.pathdf = pathdf[pathdf[[split.var]] == subtype,] 
    # TODO: Look at 
    reg = 'EC'
    # reg = 'PFC'
    sub.pathdf = sub.pathdf[sub.pathdf$region == reg,]
    sub.pathdf$age.z = scale(sub.pathdf$age_death)

    print(paste("Number of cells:", nrow(sub.pathdf)))
    ststr = gsub("/","_",gsub(" ","_", subtype))
    offset = marg[sub.pathdf$barcode, 'count']
    for (path in c(pathlist,'cogdxad', 'nrad')) {
        # TODO: For now, run without interaction terms, later add. 
        # mdx = model.matrix(asform(c('~',path, '* region + age_death + msex + pmi')), data=sub.pathdf)
        # mdx = model.matrix(asform(c('~',path, '+ region + age_death + msex + pmi')), data=sub.pathdf)
        # mdx = model.matrix(asform(c('~',path, '+ age_death + msex + pmi')), data=sub.pathdf)
        mdx = model.matrix(asform(c('~',path, '+ age.z + msex + pmi')), data=sub.pathdf)
        # md2 = model.matrix(asform(c('~',path, '')), data=sub.pathdf)

        # Output files:
        regfile = paste0(regdir, prefix, '.nblmm_reg.', path, '.', reg, '.major.', celltype, '.minor.', ststr, '.Rda')
        regtsv  = paste0(regdir, prefix, '.nblmm_reg.', path, '.', reg, '.major.', celltype, '.minor.', ststr, '.tsv.gz')

        if (!file.exists(regfile)){
            print(path)
            chunksize = 500
            nchunk = floor(ngenes / chunksize) + 1
            regdf = c()
            regdf2 = c()
            # Subset of locations to extract: 
            bind = match(sub.pathdf$barcode, barcodes)
            subbcs = barcodes[bind] 
            t0 = proc.time()
            for (i in 1:nchunk){
                t1 = proc.time()
                print(i)
                ind = (1 + (i-1) * chunksize):min(c(i * chunksize, ngenes)) 
                kg = c('SLC39A1','LINGO1','SLC26A3','APOE','MT-ND3','SIGLEC1','APP','CLU','PLP1','GFAP','DSCAM','PLIN2', 'POMC')
                ind = which(genes %in% kg)
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
                # NOTE: Works fine as both PMM, NBGMM, and NBLMM, choose NBGMM - more conservative
                # TODO: TEST GLM.NB, NEBULA, BASIC + With 
                # NOTE: RUNNING PMM FOR SPEED

                re = nebula(mat, as.character(sub.pathdf$projid), pred=mdx, offset=log(offset), model='PMM') 
                rdf = re$summary
                rdf[order(rdf$p_nradCTRL),c('logFC_nradCTRL','p_nradCTRL', 'gene')]

                library(ggplot2)
                library(ggpubr)
                sub.pathdf$offset = offset
                med.off = median(offset)
                for (gene in kg){
                    print(gene)
                    sub.pathdf$expr = mat[gene,]
                    g1 = ggplot(sub.pathdf, aes(nrad, expr, fill=nrad)) +
                        geom_violin(scale = 'width') + 
                        labs(x='AD Status', title=paste(gene, 'Expr.')) + 
                        scale_y_continuous(expand=c(0,0)) + 
                        theme_pubr()
                    g2 = ggplot(sub.pathdf, aes(nrad, expr / offset * med.off, fill=nrad)) +
                        geom_violin(scale = 'width') + 
                        labs(x='AD Status', title=paste(gene, 'Expr. (norm)')) + 
                        scale_y_continuous(expand=c(0,0)) + 
                        theme_pubr()
                    garr = ggarrange(g1,g2, ncol=2, common.legend=TRUE)
                    ggsave(paste0(imgpref, 'genedist_',gene, '_comparison_',reg,'_', cellstr,'_sub_',ststr,'_',path,'_violin.png'),garr, units='in', dpi=450, width=5, height=6)

                    sdf = sub.pathdf[,c('expr','offset','projid','nrad')]
                    sdf = sdf[order(sdf$nrad),]
                    sdf$norm = sdf$expr / sdf$offset * med.off
                    ss = aggregate(norm ~ projid, sdf, mean)
                    ss = ss[order(ss$norm),]
                    sdf$projid = factor(sdf$projid, levels=unique(ss$projid))
                    g3 = ggplot(sdf, aes(factor(projid), norm, fill=nrad)) +
                        geom_violin(scale = 'width', alpha=0.25) + 
                        geom_boxplot(width=.25) + 
                        labs(x='Individual', title=paste(gene, 'Expr. (norm)')) + 
                        scale_y_continuous(expand=c(0,0)) + 
                        theme_pubr() +
                        theme(axis.text.x=element_text(angle=90, hjust=1))
                    ggsave(paste0(imgpref, 'genedist_ind_',gene, '_comparison_',reg,'_', cellstr,'_sub_',ststr,'_',path,'_violin.png'),g3, units='in', dpi=450, width=8, height=6)
                }

                    # stat_compare_means(label='p.signif')



                # (As if no individual splits):
                re2 = nebula(mat, rep(1, nrow(sub.pathdf)), pred=mdx, offset=log10(offset), model='PMM') 
                rdf2 = re2$summary
                rdf2[order(rdf2$p_nradCTRL),c('logFC_nradCTRL','p_nradCTRL', 'gene')]

                # # Removing pmi etc. doesnt change
                # re = nebula(mat, as.character(sub.pathdf$projid), pred=md2, offset=offset) 
                # rdf = re$summary
                # rdf[order(rdf$p_nradCTRL),c('logFC_nradCTRL','p_nradCTRL', 'gene')]
                # # With poisson mixed model; only slightly more sign.
                # re = nebula(mat, as.character(sub.pathdf$projid), pred=mdx, offset=offset, model='PMM') 
                # rdf = re$summary
                # rdf[order(rdf$p_nradCTRL),c('logFC_nradCTRL','p_nradCTRL', 'p_pmi','gene')]
                # # With log10 offset:
                # re = nebula(mat, as.character(sub.pathdf$projid), pred=mdx, offset=log10(offset), model='PMM') 
                # rdf = re$summary
                # rdf[order(rdf$p_nradCTRL),c('logFC_nradCTRL','p_nradCTRL', 'p_pmi','gene')]
                # re = nebula(mat, rep(1, nrow(sub.pathdf)), pred=mdx, offset=log10(offset), model='PMM') 
                # rdf = re$summary
                # rdf[order(rdf$p_nradCTRL),c('logFC_nradCTRL','p_nradCTRL', 'p_pmi','gene')]

                # # With glm.nb:
                # # TODO: Compare EC here with EC strained
                # library(MASS)
                # sdf = sub.pathdf
                # sdf$marg = marg[sdf$barcode, 'count']
                # acdf = c()
                # for(gene in kg){
                #     print(gene)
                #     sdf$expr = mat[gene,sdf$barcode]
                #     # fit = glmer.nb(expr ~ nrad + age_death + msex + pmi + offset(log(marg)) + (1|projid), sdf)
                #     fit = glm.nb(expr ~ nrad + age_death + msex + pmi + offset(log(marg)), sdf)
                #     cdf = as.data.frame(t(unlist(coefficients(summary(fit))['nradCTRL',])))
                #     cdf$gene = gene
                #     acdf = rbind(acdf, cdf)
                # }
                # cd = as.data.frame(acdf)
                # colnames(cd) = c('est','se','z','p','gene')
                # cd[order(cd$p),]

                if (length(ind) = length(kg)){
                    # With glmer:
                    library(lme4)
                    sdf = sub.pathdf
                    sdf$marg = marg[sdf$barcode, 'count']
                    # zscore the vars
                    sdf$age.z = scale(sdf$age_death)
                    acdf = c()
                    for(gene in kg){
                        print(gene)
                        sdf$expr = mat[gene,sdf$barcode]
                        # fit = glmer.nb(expr ~ nrad + age_death + msex + pmi + offset(log(marg)) + (1|projid), sdf)
                        fit = glmer.nb(expr ~ nrad + age.z + msex + pmi + offset(log(marg)) + (1|projid), sdf)
                        cdf = as.data.frame(t(unlist(coefficients(summary(fit))['nradCTRL',])))
                        cdf$gene = gene
                        acdf = rbind(acdf, cdf)
                    }
                    # Plot comparison to the nebula results:
                    acdf$p = acdf$'Pr(>|z|)'
                    acdf[order(acdf$p),c('Estimate','p','gene')]
                    df = merge(acdf[,c('p','gene')], rdf[,c('p_nradCTRL','gene')])
                    df$p.nebula = -log10(df$p_nradCTRL)
                    df$p.glmer.nb = -log10(df$p)
                    # Compare these in EC:
                    df$not.in = factor(df$gene %in% c('EYA2','GFAP','PLP1','NEXMIF','TMEM109', 'SIGLEC1'))

                    library(ggrepel)

                    g1 = ggplot(df, aes(p.nebula, p.glmer.nb, label=gene, color=not.in)) + 
                        geom_smooth(method='lm', alpha=.25, color='blue') + 
                        geom_point() + 
                        scale_color_manual(values=c('black','red'),name='Not in Mic.') + 
                        geom_text_repel() + 
                        labs(x='log10p for NBGMM', y='log10p for GLMER.NB') + 
                        theme_pubr()
                    ggsave(paste0(imgpref, 'neb_glmer_comparison_',reg,'_', cellstr,'_sub_',ststr,'_',path,'_lp_scatter.png'),g1, units='in', dpi=450, width=5, height=5)


                }


                if (nrow(rdf) > 0){ regdf = rbind(regdf, rdf) }
                if (nrow(rdf2) > 0){ regdf2 = rbind(regdf2, rdf2) }
                print(dim(regdf))
                print(dim(regdf2))
                # Timing:
                t2 = (proc.time() - t1)[3]
                ttot = (proc.time() - t0)[3]
                names(t2) = NULL
                names(ttot) = NULL
                cat(paste0(round(t2,1),'s\t'))
                cat(paste0(round(ttot,1),'s\n'))
            }
            regdf$path = path
            regdf2$path = path
            if (path == 'cogdxad'){
                pval = paste0('p_', path,'NCI')
            } else if (path == 'nrad'){
                pval = paste0('p_', path,'CTRL')
            } else {
                pval = paste0('p_',path)
            }
            regdf$q_path = qvalue(regdf[[pval]])$q
            regdf$log10q = -log10(regdf$q_path)
            regdf = regdf[order(regdf$log10q, decreasing=T),]
            regdf2$q_path = qvalue(regdf2[[pval]])$q
            regdf2$log10q = -log10(regdf2$q_path)
            regdf2 = regdf2[order(regdf2$log10q, decreasing=T),]
            # Save:
            save(regdf, regdf2, file=regfile)
            write.table(regdf, file=gzfile(regtsv), quote=F, row.names=F, sep="\t")
        } else {
            load(regfile)
        }


        if (path == 'cogdxad'){
            eff = paste0('logFC_',path, 'NCI')
        } else if (path == 'nrad'){
            eff = paste0('logFC_', path,'CTRL')
        } else {
            eff = paste0('logFC_',path)
        }
        print(head(regdf[regdf[[eff]] > 0,c(eff,'log10q','gene')], 50))
        print(head(regdf[regdf[[eff]] < 0,c(eff,'log10q','gene')], 50))
        print(head(regdf[,c(eff,'log10q','gene')], 50))

        print(head(regdf2[regdf2[[eff]] > 0,c(eff,'log10q','gene')], 50))
        print(head(regdf2[,c(eff,'log10q','gene')], 50))

    }
}




# ---------------------------------------------------------
# Plot the differences between with/without random-effects:
# + plot summaries of the results for each.
# ---------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(ggrepel)
rd = regdf[,c('gene', 'logFC_nradCTRL', 'p_nradCTRL')]
rd2 = regdf2[,c('gene', 'logFC_nradCTRL', 'p_nradCTRL')]
colnames(rd) = c('gene','est.lmm','p.lmm')
colnames(rd2) = c('gene','est.fixed','p.fixed')
rdf = merge(rd,rd2)
rdf$lp.fixed = -log10(rdf$p.fixed)
rdf$lp.lmm = -log10(rdf$p.lmm)

gp = ggplot(rdf[abs(rdf$est.fixed) < 4,], aes(est.fixed, est.lmm)) + 
    geom_abline(color='grey50') + geom_point(alpha=0.25, cex=.5) +
    theme_pubr() + 
    labs(x='Effect size (model without mixed eff.)', y='Effect size (model with mixed eff.)')
ggsave(paste0(imgpref, 'model_comparison_lmm_',reg,'_', cellstr,'_sub_',ststr,'_',path,'_eff_scatter.png'),gp, units='in', dpi=450, width=4, height=4)

rdf$col = (rdf$lp.fixed > 2) + (rdf$lp.lmm > 2)

gp = ggplot(rdf, aes(lp.fixed, lp.lmm, color=factor(col))) + 
    scale_x_log10() + 
    scale_y_log10() + 
    scale_color_manual(values=c('grey75','grey35','red')) +
    geom_hline(yintercept=2, lwd=.25) + 
    geom_vline(xintercept=2, lwd=.25) + 
    geom_point(pch=16, alpha=.5, cex=.25) + theme_pubr() + 
    theme(legend.position = 'none') + 
    labs(x='-log10p (model without mixed eff.)', y='-log10p (model with mixed eff.)')
ggsave(paste0(imgpref, 'model_comparison_lmm_',reg,'_', cellstr,'_sub_',ststr,'_',path,'_lp_scatter.png'),gp, units='in', dpi=450, width=4, height=4)



gp = ggplot(rdf[rdf$est.fixed > 0,], aes(lp.fixed, lp.lmm, color=factor(col))) + 
    scale_x_log10() + 
    scale_y_log10() + 
    scale_color_manual(values=c('grey75','grey35','red')) +
    geom_hline(yintercept=2, lwd=.25) + 
    geom_vline(xintercept=2, lwd=.25) + 
    geom_point(pch=16, alpha=.5, cex=.5) + theme_pubr() + 
    theme(legend.position = 'none') + 
    labs(x='-log10p (model without mixed eff.)', y='-log10p (model with mixed eff.)')
ggsave(paste0(imgpref, 'model_comparison_lmm_',reg,'_', cellstr,'_sub_',ststr,'_',path,'_lp_scatter_neg.png'),gp, units='in', dpi=450, width=4, height=4)


# Which are high in lp.fixed but not signif in lp.lmm:
# rdf[(rdf$lp.fixed > 100) & (rdf$lp.lmm < 2),]
# ggplot(rdf, aes(p.fixed, p.lmm)) + 
#     geom_point() + theme_pubr()

head(rdf[order(rdf$p.fixed),], 25)
head(rdf[order(rdf$p.lmm),], 25)

hist(rdf$p.lmm)
hist(rdf$p.fixed)

sum(rdf$p.lmm < 0.05)
sum(rdf$p.fixed < 0.05)

# rdf$p.adj = p.adjust(rdf$p.lmm, 'fdr')
rdf$p.adj = p.adjust(rdf$p.lmm, 'fdr')
sum(rdf$p.adj < 0.05)
head(rdf[order(rdf$p.adj),c('gene','est.lmm','lp.fixed','lp.lmm','p.adj')], 33)

rdf$p.adj2 = p.adjust(rdf$p.fixed, 'fdr')
sum(rdf$p.adj2 < 0.05)
head(rdf[order(rdf$p.adj2),c('gene','est.fixed','lp.fixed','lp.lmm','p.adj2')], 33)

# Look at pathways:
rd = rdf[order(rdf$p.adj),]
g1 = head(rd[rd$est.lmm < 0,'gene'],100)
g1d = head(rd[rd$est.lmm > 0,'gene'],100)

rd2 = rdf[order(rdf$p.adj2),]
g2 = head(rd2[rd2$est.fixed < 0,'gene'],100)
g2d = head(rd2[rd2$est.fixed > 0,'gene'],100)

library(gprofiler2)
go1 = gost(g1, ordered_query=TRUE)$result
go1d = gost(g1d, ordered_query=TRUE)$result
go2 = gost(g2, ordered_query=TRUE)$result
go2d = gost(g2d, ordered_query=TRUE)$result


print(go1[go1$source == 'REAC',c('p_value','term_name')])
print(go2[go2$source == 'REAC',c('p_value','term_name')])

print(go1[go1$source == 'KEGG',c('p_value','term_name')])
print(go2[go2$source == 'KEGG',c('p_value','term_name')])

print(go1d[go1d$source == 'REAC',c('p_value','term_name')])
print(go2d[go2d$source == 'REAC',c('p_value','term_name')])

print(go1d[go1d$source == 'KEGG',c('p_value','term_name')])
print(go2d[go2d$source == 'KEGG',c('p_value','term_name')])

# rd[rd$gene == 'SLC38A2',]
# rd2[rd2$gene == 'SLC38A2',]

hist(go1$p_value)
hist(go1d$p_value)

rdf[rdf$est.fixed < -2,]


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
