#!/usr/bin/R
# --------------------------------------
# Set up a matrix for DE runs for SI/ACE
# Updated: 10/05/22
# --------------------------------------
library(cbrbase)
set_proj('DEVTRAJ', 'multiRegion')
source(paste0(sbindir, 'load_metadata.R'))


# ----------------------------------------------------
# List all of the cell type contexts that we will use:
# ----------------------------------------------------
df = data.frame(celltype='Mic_Immune',subtype=c('Mic','CAM','T'))

vasc.sts = unique(cellmeta$cell_type_high_resolution[cellmeta$major.celltype == 'Vasc/Epithelia'])
df = rbind(df, data.frame(celltype='Vasc_Epithelia', subtype=vasc.sts))

cts = c('Ast','Opc','Oli','Inh','Exc')
df = rbind(df, data.frame(celltype=cts, subtype=cts))


# -----------------------------------------------
# Merge the regions + AD ascertainment variables:
# -----------------------------------------------
extsi_tsv = 'Annotation/extended_si_ace_vars.tsv'
ext.si.vars = read.delim(extsi_tsv, header=F)[,1]
pathlist = c('msex', 'soc_net_bl', 'social_isolation_avg',
    'social_isolation_lv', ext.si.vars)

reglist = c('allregions', reg.nomb)
pdf = expand.grid(region=reglist, path=pathlist)

# Merge all:
testdf = merge(df, pdf)

# Write out:
write.table(testdf, paste0(sdbdir, 'DEG_multiRegion_SI_ACE_runlist.tsv'), 
    quote=F, row.names=F, sep="\t")
