#!/usr/bin/R
# ------------------------------------------
# Plot results for Nebula on e4 interaction:
# - Used for Oli-M7 + Opc-M9 test
# Updated: 06/27/23
# ------------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))

library(tidyr)
library(Matrix)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
print(version)
options(width=170)

source(paste0(sbindir, 'auxiliary_plotting_settings.R'))


# Set the directories:
# --------------------
topimgdir = paste0(img, 'multiRegion/')
plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0(datadir,'dereg/')
imgpref = paste0(plotdir, 'difftl_e4int_')
cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
system(cmd)


# Arguments (for now just Oli/Opc M7/M9)
# --------------------------------------
celltype = 'Opc'

modmap = c('Oli'=7, 'Opc'=9)
subtype = celltype
module = modmap[celltype]
region = 'allregions'
# path = 'plaq_n'
path = 'cogdxad'
mtag = paste0(celltype, '-M', module)
run.interaction = FALSE


# Read run output files:
# ----------------------
prefstr = gsub("[()/]","_", paste(celltype, subtype, region, path, module, sep="_"))
suffstr = ifelse(run.interaction, 'module_e4', 'module_e4noint')
respref = paste0(regdir, 'nebula_ruv.', suffstr, '.', prefstr)
outtsv = paste0(respref, '.tsv.gz')
resdf = read.delim(gzfile(outtsv), header=T, check.names=FALSE)


# Basic stats on # of DE
# ----------------------
pathstr = ifelse(path == 'cogdxad', paste0(path,'AD'), path)
e4str = 'Apoe_e4yes'
intstr = paste0(pathstr, ':', e4str)
if (run.interaction){
    kept.eff = c(intstr, pathstr, e4str)
} else {
    kept.eff = c(pathstr, e4str)
}

pcut = 0.05
for (eff in kept.eff){
    leff = paste0('logFC_', eff)
    peff = paste0('p_', eff)
    qeff = paste0('q_', eff)
    ceff = paste0('col_', eff)

    resdf[[qeff]] = p.adjust(resdf[[peff]], 'BH')
    resdf[[ceff]] = ifelse(resdf[[qeff]] < pcut, ifelse(resdf[[leff]] < 0, 'Down', 'Up'), 'NS')
    print(eff)
    print(table(resdf[[ceff]]))
    head(resdf[resdf[[ceff]] != 0,], 20)
}


# Plot pairs of effects:
# ----------------------
if (run.interaction){
    comp.pairs = list(
        c(pathstr, e4str),
        c(pathstr, intstr), 
        c(e4str, intstr))
} else {
    comp.pairs = list(c(pathstr, e4str))
}

pair.cols = brewer.pal(12, 'Paired')
dir.cols = c(
    'NS/NS'='grey90',
    'Down/NS' = pair.cols[1],
    'Up/NS' = pair.cols[5],
    'Down/Down' = pair.cols[2],
    'Up/Up' = pair.cols[6],
    'Up/Down' = pair.cols[3],
    'Down/Up' = pair.cols[7],
    'NS/Down' = pair.cols[4],
    'NS/Up' = pair.cols[8])


for (pair in comp.pairs){
    print(pair)
    leff = paste0('logFC_', pair)
    ceff = paste0('col_', pair)

    resdf$joint_col = paste0(resdf[[ceff[1]]], '/', resdf[[ceff[2]]])
    resdf$lab = ifelse(resdf$joint_col == 'NS/NS', '', resdf$gene)

    gp = ggplot(resdf, aes(.data[[leff[1]]], .data[[leff[2]]], label=lab, col=factor(joint_col))) + 
        rasterize(geom_point(cex=.25), dpi=450) + 
        geom_text_repel(cex=2.5, max.overlaps=10, force=5, 
            min.segment.length=0.04, point.padding=0.05, box.padding=0.05) + 
        scale_color_manual(values=dir.cols, name='Signif.') + 
        geom_hline(yintercept=0, lty='dashed') + 
        geom_vline(xintercept=0, lty='dashed') + 
        labs(x=paste0('logFC (', pair[1], ')'), y=paste0('logFC (', pair[2], ')')) + 
        theme_pubr() + theme(legend.position='right')
    pltpref = paste0(imgpref, 'scattercomp_', mtag, '_', paste0(pair, collapse='-'))
    saveGGplot(gp, pltpref, w=4.75, h=3.25)
    saveGGplot(gp, paste0(pltpref, '_large'), w=8, h=7)
}


table(resdf$joint_col)
resdf[resdf$joint_col %in% c('Up/Up', 'Down/Down'),]

