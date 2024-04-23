#!usr/bin/python
"""scDRS on each cell type alone."""
# wget https://figshare.com/ndownloader/articles/19312583/versions/1
# ----------------------------------------------
# scDRS (v1.0.2) scores on the multiregion data:
# Updated: 06/21/23
# ----------------------------------------------
import os
import scdrs
import scanpy as sc
import pandas as pd
import numpy as np


# Load in all of the multiregion data at once:
# --------------------------------------------
# NOTE: size is 1.35M cells by 18k genes
dbdir = '/home/cboix/data/DEVTRAJ/db/'
datadir = dbdir + 'multiRegion/'
h5ad_file = datadir + 'preprocessed_formodules_All_log1p.h5ad'
adata = sc.read(h5ad_file)


# Load UMAP coords and metadata:
# ------------------------------
prefix = 'all_brain_regions_filt_preprocessed_scanpy_norm'
metafile = datadir + prefix + '.final_noMB.cell_labels.tsv.gz'
metadf = pd.read_csv(metafile, sep="\t")
metadf.index = metadf.barcode
metadf = metadf[['U1','U2', 'major.celltype', 'cell_type_high_resolution']]

# Merge into adata:
mdf2 = pd.merge(adata.obs, metadf, how='left', left_index=True, right_index=True)
adata.obs['U1'] = mdf2['U1']
adata.obs['U2'] = mdf2['U2']
adata.obs['major.celltype'] = mdf2['major.celltype'].astype('category')
adata.obs['celltype'] = mdf2['cell_type_high_resolution'].astype('category')

# Add columns for MCTxRegion and CTHRxREGION
zmct = zip(adata.obs['major.celltype'].tolist(), adata.obs['region'].tolist())
zcthr = zip(adata.obs['celltype'].tolist(), adata.obs['region'].tolist())
adata.obs['mct_region'] = [x + "_" + y for (x,y) in zmct]
adata.obs['cthr_region'] = [x + "_" + y for (x,y) in zcthr]
adata.obs['mct_region'] = adata.obs['mct_region'].astype('category')
adata.obs['cthr_region'] = adata.obs['cthr_region'].astype('category')


# Pre-process for scDRS (takes the longest):
# ------------------------------------------
# Note: Add covariates if need to regress out AD/non-AD, M/F
scdrs.pp.preprocess(adata)

# Required for grouping
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)


# TODO: Save processed scDRS file


# Load scDRS files:
# -----------------
gsdir = datadir + 'scDRS/gs_file/'
outdir = datadir + 'scDRS/results/'
gsfile = gsdir + 'magma_10kb_top1000_zscore.74_traits.rv1.gs'
df_gs = scdrs.util.load_gs(gsfile)

gslist = df_gs.keys()
gslist = [x for x in gslist]


# Score each of the GWAS:
# -----------------------
def get_cell_scores(gs, df_gs, outdir=outdir, adata=adata):
    gsoutfile = outdir + 'full.scores-' + gs + '.tsv.gz'
    if os.path.exists(gsoutfile):
        df_res = pd.read_csv(gsoutfile, sep="\t", index_col=0)
    else:
        gene_list = df_gs[gs][0]
        gene_weight = df_gs[gs][1]
        df_res = scdrs.score_cell(
            adata, gene_list, gene_weight=gene_weight, n_ctrl=20)
        df_res.to_csv(gsoutfile, sep="\t")
    return(df_res)


def get_group_scores(gs, df_res, ctcol, adata=adata):
    df_stats = scdrs.method.downstream_group_analysis(
        adata=adata, df_full_score=df_res, group_cols=[ctcol])[ctcol]
    order = np.argsort(df_stats['assoc_mcp'].to_numpy())
    df_stats = df_stats.iloc[order]
    df_stats = df_stats.reset_index()
    df_stats['gwas'] = gs
    df_stats['ctcol'] = ctcol
    return(df_stats)


# Run on celltypes + major cell types:
# ------------------------------------
fulldf = None
ctcols = ['celltype', 'major.celltype', 'cthr_region', 'mct_region']
for i, gs in enumerate(gslist):
    print(i, gs)
    df_res = get_cell_scores(gs, df_gs)
    # Group level + aggregate
    gwas_stats = None
    for ctcol in ctcols:
        df_stats = get_group_scores(gs, df_res, ctcol)
        gwas_stats = pd.concat([gwas_stats, df_stats])
    fulldf = pd.concat([fulldf, gwas_stats])
    if gs == 'PASS_Alzheimers_Jansen2019':
        ad_outfile = outdir + 'group.scores-ADGWAS.tsv.gz'
        gwas_stats.to_csv(ad_outfile, sep="\t")

overall_outfile = outdir + 'group.scores-allGWAS.tsv.gz'
fulldf.to_csv(overall_outfile, sep="\t")


# TODO: Plot AD scores on UMAP
# ----------------------------

