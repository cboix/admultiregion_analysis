#!/usr/bin/R
# -----------------------------------------------------------
# Use the pathology to run some tests on linear mixed models:
# Using Liang's Nebula package to run the fast NBLMM
# TEST on the PFC data:
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
library(MASS)

library(ggplot2)
library(ggpubr)
library(ggrepel)

load('10x_processed/filtered_matrix_attr.Rda')

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'pfcrerun_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }
pathlist = c('nft','plaq_d','plaq_n')

# ----------------
# Add in metadata:
# ----------------
meta = read.delim('Annotation/metadata_PFC_all_individuals_092520.tsv')
dim(rna.df)
rna.df = merge(rna.df, meta[,c('projid','age_death','pmi','msex',
                               'nft','plaq_d','plaq_n','niareagansc','cogdx','braaksc')])
rna.df$barcode = paste0('RNA_',rna.df$TAG)
rna.df$cogdxad = 'NCI'
rna.df$cogdxad[rna.df$cogdx %in% c(4,5)] = 'AD'
rna.df$nrad = 'CTRL'
rna.df$nrad[rna.df$niareagansc %in% c(1,2)] = 'AD'
rna.df$nrad = factor(rna.df$nrad, levels=c('CTRL','AD'))
pqdf = unique(rna.df[, c('projid','age_death','pmi','msex',
                         'nft','plaq_d','plaq_n','nrad','cogdxad')])


# ----------------------------------------
# Run regressions on a specific cell type:
# ----------------------------------------
celltype = 'Mic'
sub.pathdf = rna.df[rna.df$broad.cell.type == celltype,]
path = 'nrad'

# For offset:
marg = colSums(rna.mat)

# -----------------------------------------------------------
# Run the regression on a set of data:
# Req: count matrix, random effects, model matrix, and offset
# -----------------------------------------------------------
# Run on each subtype:
celltypes = as.character(unique(rna.df$broad.cell.type))
for (ct in celltypes){
    print(paste("[STATUS] Running on", ct))
    sub.pathdf = rna.df[rna.df$broad.cell.type == ct,]
    print(paste("Number of cells:", nrow(sub.pathdf)))
    ststr = gsub("/","_",gsub(" ","_", ct))
    offset = marg[sub.pathdf$barcode]
    # for (path in c('nrad',pathlist,'cogdxad')) {
    for (path in c('nrad')) {
        # TODO: For now, run without interaction terms, later add. 
        # mdx = model.matrix(asform(c('~',path, '* region + age_death + msex + pmi')), data=sub.pathdf)
        mdx = model.matrix(asform(c('~',path, '+ age_death + msex + pmi')), data=sub.pathdf)
        # Output files:
        regfile = paste0(regdir, 'pfcrerun', '.nblmm_reg.', path, '.major.', ststr, '.Rda')
        regtsv  = paste0(regdir, 'pfcrerun', '.nblmm_reg.', path, '.major.', ststr, '.tsv.gz')

        if (path == 'cogdxad'){
            pathstr = paste0(path,'NCI')
        } else if (path == 'nrad'){
            pathstr = paste0(path,'AD')
        } else { pathstr = path }
        pval = paste0('p_', pathstr)
        eff = paste0('logFC_',pathstr)
        if (!file.exists(regfile)){
            print(path)
            chunksize = 1000
            ngenes = nrow(rna.mat)
            nchunk = floor(ngenes / chunksize) + 1
            regdf = c(); regdf2 = c();
            # Subset of locations to extract: 
            bind = match(sub.pathdf$barcode, colnames(rna.mat))
            subbcs = colnames(rna.mat)[bind] 
            t0 = proc.time()
            for (i in 1:nchunk){
                t1 = proc.time()
                print(i)
                ind = (1 + (i-1) * chunksize):min(c(i * chunksize, ngenes)) 
                mat = rna.mat[ind, bind]
                mat = mat[,sub.pathdf$barcode]
                gcout = gc()
                # NOTE: Works fine as both PMM and as NBLMM, choose NBLMM - more conservative
                re = nebula(mat, as.character(sub.pathdf$projid), pred=mdx, offset=log10(offset))
                rdf = re$summary
                re2 = nebula(mat, rep(1,nrow(sub.pathdf)), pred=mdx, offset=log10(offset))
                rdf2 = re2$summary
                if (nrow(rdf) > 0){ regdf = rbind(regdf, rdf) }
                if (nrow(rdf2) > 0){ regdf2 = rbind(regdf2, rdf2) }
                # Timing:
                t2 = (proc.time() - t1)[3]
                ttot = (proc.time() - t0)[3]
                names(t2) = NULL
                names(ttot) = NULL
                cat(paste0(round(t2,1),'s\t'))
                cat(paste0(round(ttot,1),'s\n'))
            }


            # --------------------------------------
            # Also run the pseudo-bulk glm.nb model:
            # --------------------------------------
            # Run standard glm.nb on this matrix:
            print("[STATUS] Running pseudo-bulk")
            tform = make.tform(as.character(sub.pathdf$projid))
            tmat = as.matrix(t(rna.mat[,sub.pathdf$barcode] %*% tform))
            gcout = gc()
            tdf = data.frame(tmat)
            tdf$projid = rownames(tdf)
            tdf = gather(tdf, gene, expr, -projid)
            tdf = merge(tdf, pqdf)
            marg2 = agg.rename(expr ~ projid, tdf, sum, 'count')
            tdf = merge(tdf, marg2)
            # tdf$cpm = log2((tdf$expr +1)/ tdf$count * 1e6)

            # Run regression:
            acdf = c()
            genes = unique(tdf$gene)
            for (i in 1:nrow(regdf)){
                gene = regdf$gene[i]
                if (i %% 100 == 0){ print(i) }
                sdf = tdf[tdf$gene == gene,]
                fit = try(glm.nb(asform(c('expr ~ ',path,'+ age_death + msex + pmi + offset(log10(count))')), sdf))
                if (class(fit) != 'try-error'){
                    cdf = as.data.frame(t(unlist(coefficients(summary(fit))[pathstr,])))
                    cdf$gene = gene
                    acdf = rbind(acdf, cdf)
                }
            }
            pseudf = as.data.frame(acdf)
            colnames(pseudf) = c('est','se','z','p','gene')
            pseudf$gene = sub("\\.","-", pseudf$gene)
            pseudf = pseudf[order(pseudf$p),]
            pseudf$q_path = qvalue(pseudf$p)$q
            pseudf$log10q = -log10(pseudf$q_path)
            pseudf$path = path 
            hist(pseudf$p)

            regdf$path = path
            regdf2$path = path
            regdf$q_path = qvalue(regdf[[pval]])$q
            regdf$log10q = -log10(regdf$q_path)
            regdf = regdf[order(regdf$log10q, decreasing=T),]
            regdf2$q_path = qvalue(regdf2[[pval]])$q
            regdf2$log10q = -log10(regdf2$q_path)
            regdf2 = regdf2[order(regdf2$log10q, decreasing=T),]
            # Save:
            save(regdf, regdf2, pseudf, file=regfile)
            write.table(regdf, file=gzfile(regtsv), quote=F, row.names=F, sep="\t")
        } else {
            load(regfile)
        }
        print(head(regdf[regdf[[eff]] > 0,c(eff,'log10q','gene')], 50))
        print(head(regdf[,c(eff,'log10q','gene')], 50))
        print(head(regdf2[regdf2[[eff]] > 0,c(eff,'log10q','gene')], 50))
        print(head(regdf2[,c(eff,'log10q','gene')], 50))

        # ---------------------------------
        # Plot a comparison of the results:
        # ---------------------------------
        rd = regdf[,c('gene', eff, 'log10q')]
        rd2 = regdf2[,c('gene', eff, 'log10q')]
        cd = pseudf[,c('gene', 'est', 'log10q')]
        colnames(rd) = c('gene','est.lmm','lq.lmm')
        colnames(rd2) = c('gene','est.fixed','lq.fixed')
        colnames(cd) = c('gene','est.pbulk','lq.pbulk')
        rdf = merge(merge(rd,rd2), cd, all.x=TRUE)

        gp = ggplot(rdf[abs(rdf$est.fixed) < 4,], aes(est.fixed, est.lmm)) + 
            geom_abline(color='grey50') + geom_point(alpha=0.25, cex=.5) +
            theme_pubr() + 
            labs(x='Effect size (model without mixed eff.)', y='Effect size (model with mixed eff.)')
        ggsave(paste0(imgpref, 'model_comparison_lmm_',ststr,'_',path,'_eff_scatter.png'),gp, units='in', dpi=450, width=3.5, height=3)

        rdf$col = (rdf$lq.fixed > 2) + (rdf$lq.lmm > 2)
        labdf = rdf[rdf$col > 0,]
        gp = ggplot(rdf, aes(lq.fixed, lq.lmm, color=factor(col))) + 
            scale_x_log10(expand=c(0,0)) + 
            scale_y_log10(expand=c(0,0)) + 
            scale_color_manual(values=c('grey50','blue','red')) +
            geom_hline(yintercept=2, lwd=.25) + 
            geom_vline(xintercept=2, lwd=.25) + 
            geom_point(pch=19, alpha=1, cex=.25) + theme_pubr() + 
            geom_text_repel(data=labdf, aes(lq.fixed,lq.lmm, label=gene, color=factor(col)), cex=2, segment.size=.25) +
            theme_pubr() + 
            theme(legend.position = 'none') + 
            labs(x='-log10q (model without mixed eff.)', y='-log10q (model with mixed eff.)', title=ct)
        ggsave(paste0(imgpref, 'model_comparison_lmm_',ststr,'_',path,'_lp_scatter.png'),gp, units='in', dpi=450, width=3.5, height=3.5)

        png(paste0(imgpref, 'model_comparison_lmm_',ststr,'_',path,'_paired_scatter.png'),units='in', res=450, width=5, height=5)
        plot(rdf[,c('lq.lmm','lq.pbulk','lq.fixed')], pch=19, cex=.5, gap=.2)
        dev.off()

    }
}



print("Finished running regressions.")
