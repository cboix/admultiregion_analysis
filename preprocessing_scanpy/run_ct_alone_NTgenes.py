#!usr/bin/python
"""Compute representations using just neurotransmitter-related genes."""
# --------------------------------------------------------------
# Compute a UMAP for each major celltype alone for visualization
# Updated: 09/03/21
# --------------------------------------------------------------
import scanpy as sc
import numpy as np
import pandas as pd
from adata_handling import adata_loader

# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
useset = 'NT'
ct = 'Glial'
ct = 'Exc'
st = None
covariates = ['region'] if ct == 'Exc' else None
covariates = None
usecombat = False
scdata = adata_loader(ct, st,
                      normalize=True, log1p=True,
                      usecombat=usecombat,
                      usehvg=False,
                      # n_top_genes=ngenes,
                      remove_mt=False, filter_TH=False,
                      covariates=covariates)
scdata.load()
prefstr = '_majorctUMAP_' + scdata.csuff + "_" + useset


# Read in the set of neurotransmitter genes:
# ------------------------------------------
anndir = scdata.dbdir + '/Annotation/'
regdir = scdata.dbdir + '/multiRegion/dereg/'
if useset == 'NT':
    ntfile = anndir + 'neurotransmitter_related_genes_merged.txt'
    selgenes = pd.read_csv(ntfile, header=None)[0].to_numpy()
    # elif useset == 'DEGpartitions':
    # ptnfile = regdir + 'topDEGs_bypartition.tsv'
    # ptndf = pd.read_csv(ptnfile, sep="\t")
    # selgenes = ptndf.gene.to_numpy()
else:
    pass


# Subset data to these genes:
agenes = scdata.adata.var_names.tolist()
filtgenes = [x for x in agenes if x in selgenes]
scdata.adata = scdata.adata[:, filtgenes]


# Compute UMAP from PCA with k=50 and NN=15:
# ------------------------------------------
scdata.compute_representation(n_neighbors=15, n_comps=50, seed=0)

# Save after computing these:
# scdata.save_adata()

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
adf.to_csv('multiregion' + prefstr + '_subsetUMAP.tsv', sep="\t")


# Specific genes for each cell type:
# ----------------------------------
if ct == 'Ast':
    genes = [
        'GRM3', 'GABRB1', 'GABRA2', 'GABBR2', 'SLC1A2', 'SLC1A3',
        'GRIA1', 'GRID2', 'GRIK4', 'GRIN1', 'GRIK3', 'GRIN2A',
        'GRIN2B', 'RYR1', 'RYR2', 'RYR3', 'CADPS', 'CADPS2', 'ADRA1A',
        'SLC38A1', 'SLC38A2', 'SLC6A11', 'SLC6A1', 'SLC29A4', 'NRXN1',
        'NRXN2', 'NRXN3', 'NLGN1', 'ADM', 'NF1',
        'RIMS1', 'RIMS2', 'RIMS3']
elif ct == 'Opc':
    genes = [
        'GRIA1', 'GRID1', 'GRID2',
        'GRIK1', 'GRIK2', 'GRIK3', 'GRIK4',
        'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIA4',
        'SLC22A3', 'NRXN1', 'NRXN3', 'VIPR2',
        'SYN2', 'SYN3', 'RIMS2', 'RIMS1', 'RIMS3',
        'CADPS', 'CADPS2', 'SNAP23', 'ERC2', 'SYT1',
        'SLC1A3', 'ADRA1A', 'SLC38A1',
        'GABRB1', 'GABRA2', 'GABBR2', 'GABRB3',
        'SLC1A1', 'NOS1', 'SNCA', 'PCLO', 'NF1', 'SNAP25',
        'P2RX7', 'GAD1']
elif ct == 'Mic_Immune':
    genes = [
        'P2RY12', 'GRID2', 'ERC2', 'RYR1', 'CTBP2', 'SLC1A3', 'NF1',
        'SNCA']
else:
    genes = ['GABRB1', 'GABRA2', 'GABBR2', 'SLC1A2', 'SLC1A3',
                'GRIA1', 'GRID1', 'GRID2',
                'GRIK1', 'GRIK2', 'GRIK3', 'GRIK4',
                'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIA4',
                'RYR1', 'RYR2', 'RYR3', 'CADPS', 'CADPS2', 'ADRA1A',
                'SLC38A1', 'SLC38A2', 'SLC6A11', 'SLC6A1', 'SLC29A4',
                'NRXN1', 'NRXN2', 'NRXN3', 'RIMS1', 'RIMS2', 'RIMS3',
                'NLGN1', 'GRM3', 'SLC1A3', 'P2RY12'
                ]

sc.pl.umap(scdata.adata, color=genes, frameon=False,
           save=prefstr + '_ctalone_genes.png', vmax=2.5)

sc.pl.umap(scdata.adata, color=genes, frameon=False,
           save=prefstr + '_ctalone_genes_full.png')

sc.pl.matrixplot(scdata.adata, genes, groupby='region',
                 save=prefstr + '_ctalone_region.png')
sc.pl.matrixplot(scdata.adata, genes, groupby='celltype',
                 save=prefstr + '_ctalone_celltype.png')

# Leiden + new markers:
# ---------------------
sc.tl.leiden(scdata.adata, resolution=2)
sc.pl.umap(scdata.adata, color='leiden', frameon=False, legend_loc='on data',
           save=prefstr + '_ctalone_leiden.png')

NG = 25
sc.tl.dendrogram(scdata.adata, 'leiden')
sc.tl.rank_genes_groups(scdata.adata, 'leiden')
sc.pl.rank_genes_groups(scdata.adata, n_genes=NG, sharey=False,
                        frameon=False, save=prefstr + '_ctalone_leiden.png')
sc.pl.rank_genes_groups_dotplot(scdata.adata, n_genes=5,
                                save=prefstr + '_ctalone_leiden.png')
sc.pl.rank_genes_groups_matrixplot(scdata.adata, n_genes=5,
                                   save=prefstr + '_ctalone_leiden.png')


# Scaled heatmaps of ranked genes:
scdata.adata.layers['scaled'] = sc.pp.scale(scdata.adata, copy=True).X
sc.pl.rank_genes_groups_matrixplot(
    scdata.adata, n_genes=10, save=prefstr + '_ctalone_leiden.png',
    colorbar_title='mean z-score', layer='scaled',
    vmin=-2, vmax=2, cmap='RdBu_r')
