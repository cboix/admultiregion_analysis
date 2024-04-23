#!/usr/bin/R
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ')
source(paste0(bindir, 'multiRegion/load_metadata.R'))

library(tidyr)
library(rhdf5)
library(ggplot2)
library(viridis)

load('Annotation/multiregion_celltypes_colors.Rda')
palette = viridis(100)
col_fun = function(x, pal=palette){
    bin <- cut(x, seq(0, max(x), length.out=length(palette)), include.lowest=T) 
    palette[bin] 
}

gene = 'NRGN'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need a gene or genelist")
} else {        
    gene = args[1]
    if (length(args) > 1){
        chunk=as.integer(args[2])
    } else { chunk=NULL }
}
print(gene)

# Directories:
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/geneplots/')
imgpref = paste0(plotdir, 'geneplot_')
cmd = paste('mkdir -p', topimgdir, plotdir)
system(cmd)

# Building functions for regression:
# asform = function(x){ as.formula(paste0(x, collapse='')) }
# pathlist = c('nft','plaq_d','plaq_n')
# genelist = c('MSH3', 'MLH1','PMS1','PMS2','FAN1','LIG1','BCAN','NRGN','PVALB','SYT1', 'MT-ND3','APOE','HSPA1A')
genelist = c('MT-ND3','APOE','HSPA1A', 'EYA2','NEXMIF','PLP1','GFAP', 'SLC39A1','SLC26A3','SPP1')

# -------------------------
# Load metadata and h5attr:
# -------------------------
final.rdafile = paste0(datadir, prefix, '.final_noMB.cell_labels.Rda')
load(final.rdafile)
rm(celldf)
rownames(cellmeta) = cellmeta$barcode

matdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
prefix = 'all_brain_regions_filt_preprocessed_scanpy'
# h5file = paste0(matdir, prefix, "_norm", '_bygene_fullmatrix.swmr.hdf5')

cts = sub("/","_",unique(cellmeta$hcelltype))
# Load in these from full data:
gfile = paste0(matdir, prefix, "_norm.", ".tsv.gz")
# if (!file.exists(gfile)){ 
    adf = c()
    for (ct in cts){
        print(ct)
        h5file = paste0(matdir, prefix, "_norm.majorcelltype.", ct,".hdf5")

        # Extract metadata from hdf5 (is instant - why is python so slow?)
        print(h5ls(h5file))
        h5f = H5Fopen(h5file)
        genes = h5f$genes
        barcodes = h5f$barcodes
        H5Fclose(h5f)

        # Indexes for genes:
        gid = which(genes %in% genelist)
        # Extract data:
        h5f = H5Fopen(h5file)
        h5d = h5f&"matrix"
        mat = h5d[gid,]
        H5Dclose(h5d)
        H5Fclose(h5f)
        mat = t(mat)
        colnames(mat) = genes[gid]
        df = data.frame(mat)
        df$barcode = barcodes
        adf = rbind(adf, df)
    }
    # Write the raw data out:
    write.table(adf, gzfile(gfile), quote=F, sep="\t", row.names=F)
# } else { 
#     adf = read.delim(gzfile(gfile))
# }

# Match barcodes: 
kbc = adf$barcode 
cellmeta = cellmeta[kbc,]

NCELL = nrow(cellmeta)
ind = sample(1:NCELL,NCELL, replace=FALSE)


# Plot diagnostics for each:
# 1. UMAP
# 2. Violin
# 3. ID celltype?
for (gene in genelist){
    print(gene)
    genestr = sub("-",".", gene)

    x = c(adf[[genestr]])

    cex = 0.0125
    # 1. Plot it on a UMAP
    png(paste0(imgpref, 'umap_',gene, '.png'), units='in', res=450, width=8, height=8)
    par(xaxs='i')
    par(yaxs='i')
    sp = 0.1
    bsp = 1.5
    par(mar=c(bsp,sp,2,sp))
    plot(cellmeta$U1[ind], cellmeta$U2[ind], col=col_fun(x[ind]), 
         pch=19, cex=cex, axes=F)
    # legend('topleft', legend=names(cmap), col=cmap, 
    #        pt.cex=2.5, cex=1.25, pch=19, bty='n', ncol=2)
    rect(xleft=par()$usr[1], xright=par()$usr[2],
         ybottom=par()$usr[4] + 0.001 * diff(par()$usr[3:4]),
         ytop=par()$usr[4] + 0.0725 * diff(par()$usr[3:4]), 
         col='grey85', border=NA, lwd=.5, xpd=TRUE)
    mtext(gene, side=3, cex=1.5, col='grey25', font=2, line=0.25)
    mtext('UMAP 1', side=1, line=0.25, cex=1.25)
    mtext('UMAP 2', side=2, line=-1, cex=1.25)
    dev.off()

    # Violin in major and minor celltypes:
    df = data.frame(expr=x, region=cellmeta$region, ct=cellmeta$hcelltype,
                    hct=cellmeta$cell_type_high_resolution, rind=cellmeta$rind)
    df$niareagansc = metadata[df$rind,'niareagansc']
    df$nrad = 'CTRL'
    df$nrad[df$niareagansc <= 2] = 'AD'
    df$nrad = factor(df$nrad, levels=c('CTRL','AD'))

    # 2. Plot % of cells with it for each cell type (VIOLIN)
    gp = ggplot(df, aes(ct, expr, fill=ct)) + 
        geom_violin(scale='width', alpha=.5) + 
        labs(x='Major celltype', y='Expression', title=gene) + 
        scale_fill_manual(values=major.col) + 
        theme_minimal() + 
        theme(legend.position='none')
    ggsave(paste0(imgpref, 'violin_',gene, '.png'), gp, units='in', dpi=450, width=6, height=3.5)

    # 2b. Cross region by cell type?
    gp2 = gp + facet_wrap(~region, ncol=1)
    ggsave(paste0(imgpref, 'violin_byregion_',gene, '.png'), gp2, units='in', dpi=450, width=6, height=8)

    # 2c. By AD status:
    gp = ggplot(df, aes(ct, expr, fill=nrad)) + 
        geom_violin(scale='width', alpha=.5) + 
        labs(x='Major celltype', y='Expression', title=gene) + 
        scale_fill_manual(values=c('AD'='indianred','CTRL'='royalblue'), name='Status') + 
        theme_minimal()
        # theme(legend.position=c(.9,.9))
    ggsave(paste0(imgpref, 'violin_byad_',gene, '.png'), gp, units='in', dpi=450, width=6, height=3.5)

    # 3. Show if it has specificity to any neuronal subtype (Exc or Inh)
    sdf = df[df$ct == 'Exc',]
    gp = ggplot(sdf, aes(hct, expr, fill=hct, alpha=0.5)) + 
        geom_violin(scale='width') + 
        labs(x='Minor celltype', y='Expression', title=gene) + 
        scale_fill_manual(values=tcols) + 
        theme_minimal() + 
         theme(legend.position='none') + 
         theme(axis.text.x =element_text(angle=90, hjust=1))
    ggsave(paste0(imgpref, 'violin_exc_',gene, '.png'), gp, units='in', dpi=450, width=6, height=5)

    # 3b. Inh specificity
    sdf = df[df$ct == 'Inh',]
    gp = ggplot(sdf, aes(hct, expr, fill=hct, alpha=0.5)) + 
        geom_violin(scale='width') + 
        labs(x='Minor celltype', y='Expression', title=gene) + 
        scale_fill_manual(values=tcols) + 
        theme_minimal() + 
         theme(legend.position='none') + 
         theme(axis.text.x =element_text(angle=90, hjust=1))
    ggsave(paste0(imgpref, 'violin_inh_',gene, '.png'), gp, units='in', dpi=450, width=6, height=5)

    # By AD status: 

    # TODO: Enrichment in 
    # 4. If found subtype, show if it has specificity to any region + subtype?
    # 5. Output statistics for plotting elsewhere

}




