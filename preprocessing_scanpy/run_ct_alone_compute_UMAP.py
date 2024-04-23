#!usr/bin/python
"""Compute a UMAP for each major celltype alone for visualization."""
# --------------------------------------------------------------
# Compute a UMAP for each major celltype alone for visualization
# Updated: 09/03/21
# --------------------------------------------------------------
import scanpy as sc
import pandas as pd
from adata_handling import adata_loader

# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# NOTE: Inh recipe: normalize/log1p, no-combat, hvg5k, k=50 + remove MT
# NOTE: Exc recipe: normalize/log1p, combat, hvg1k, k=50 + remove MT
ct = 'Glial'
st = None
ngenes = 1000 if ct == 'Exc' else 5000
usecombat = False if ct == 'Inh' else True
# Same run across all:
ngenes = 1000
usecombat = True
covariates = ['region'] if ct == 'Exc' else None
scdata = adata_loader(ct, st,
                      normalize=True, log1p=True,
                      usecombat=usecombat,
                      usehvg=True,
                      n_top_genes=ngenes,
                      remove_mt=True, filter_TH=False,
                      covariates=covariates)
scdata.load()
prefstr = '_majorctUMAP_' + scdata.csuff

# Compute UMAP from PCA with k=50 and NN=15:
# ------------------------------------------
scdata.compute_representation(n_neighbors=15, n_comps=50, seed=0)

# Save after computing these:
scdata.save_adata()

# Plot the UMAP for the major celltype:
# -------------------------------------
sc.pl.umap(scdata.adata, color='celltype', frameon=False,
           save=prefstr + '_ctalone.png')
sc.pl.umap(scdata.adata, color='projid', frameon=False,
           save=prefstr + '_ctalone_projid.png')
sc.pl.umap(scdata.adata, color='region', frameon=False,
           save=prefstr + '_ctalone_region.png')


# Write UMAP coordinates to file:
# -------------------------------
scdata.adata.obs['C1'] = scdata.adata.obsm['X_umap'][:, 0]
scdata.adata.obs['C2'] = scdata.adata.obsm['X_umap'][:, 1]
adf = pd.DataFrame(scdata.adata.obs)[['barcode', 'C1', 'C2']]
adf.to_csv('multiregion_majorctUMAP_' + scdata.csuff + '.tsv', sep="\t")

# Specific genes for each cell type:
# ----------------------------------
if ct == 'Inh':
    genes = ['PVALB', 'SST', 'LAMP5', 'VIP', 'GPC5', 'RELN']
    sc.pl.umap(scdata.adata, color=genes, frameon=False,
               save=prefstr + '_ctalone_genes.png')

# Leiden + new markers:
# ---------------------
sc.pl.umap(scdata.adata, color='leiden', frameon=False,
           save=prefstr + '_ctalone_leiden.png')
sc.pl.umap(scdata.adata, color='leiden', frameon=False, legend_loc='on data',
           save=prefstr + '_ctalone_leiden2.png')

NG = 25
sc.tl.rank_genes_groups(scdata.adata, 'leiden')
sc.pl.rank_genes_groups(scdata.adata, n_genes=NG, sharey=False,
                        frameon=False, save=prefstr + '_ctalone_leiden.png')
sc.pl.rank_genes_groups_dotplot(scdata.adata, n_genes=5,
                                save=prefstr + '_ctalone_leiden.png')
sc.pl.rank_genes_groups_matrixplot(scdata.adata, n_genes=5,
                                   save=prefstr + '_ctalone_leiden.png')
