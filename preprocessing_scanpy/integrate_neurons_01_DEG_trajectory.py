#!usr/bin/python
"""Compute integrated neuronal representation using neuronal DEGs."""
# --------------------------------------------------------------
# Compute an integrated UMAP across all neuron subtypes:
# Updated: 03/01/22
# --------------------------------------------------------------
import re
import scanpy as sc
import numpy as np
import pandas as pd
from adata_handling import adata_loader
from adata_handling import plot_umap_category_grid
import matplotlib.pyplot as plt

from pygam import PoissonGAM, LinearGAM
from pygam import GAM, s, f


# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
useset = 'DEGpartitions'
useset = 'allDEGs'
useset = 'testedgenes'
ct = 'Glial'
ct = 'Exc'
st = None
covariates = None
usecombat = False
scdata = adata_loader(ct, st,
                      normalize=True, log1p=True,
                      usecombat=usecombat, usehvg=False,
                      remove_mt=False, filter_TH=False)
scdata.load()
prefstr = '_majorctUMAP_' + scdata.csuff + "_" + useset
# scdata.save_adata()


# Read in the set of neuronal tested genes or DEGs:
# -------------------------------------------------
anndir = scdata.dbdir + '/Annotation/'
regdir = scdata.dbdir + '/multiRegion/dereg/'
if useset == 'DEGpartitions':
    ptnfile = regdir + 'topDEGs_bypartition.tsv'
    ptndf = pd.read_csv(ptnfile, sep="\t")
    # ptndf = ptndf.loc[ptndf.dir == 'Up', :] # Up-regulated genes only
    # ptndf = ptndf.loc[ptndf.process != 'Middle', :]
    selgenes = ptndf.gene.to_numpy()
elif useset == 'allDEGs':
    ptnfile = regdir + 'allDEGs_neuronal_genes.tsv'
    selgenes = pd.read_csv(ptnfile, sep="\t", header=None)[0].to_numpy()
else:
    ptnfile = regdir + 'DEtested_neuronal_genes.tsv'
    selgenes = pd.read_csv(ptnfile, sep="\t", header=None)[0].to_numpy()


# Remove the chrY genes (othw. splits M/F in some clusters):
# ----------------------------------------------------------
gencode_version = 'v28lift37.annotation'
adf = pd.read_csv(anndir + 'Gene.' + gencode_version + '.bed',
                  sep="\t", header=None)
adf.columns = ['chr', 'start', 'end', 'strand', 'gene', 'type', 'symbol']
chrYgenes = adf.loc[adf.chr == 'Y', 'symbol'].tolist()
selgenes = [x for x in selgenes if x not in chrYgenes]


# Subset data to these genes:
# ---------------------------
agenes = scdata.adata.var_names.tolist()
filtgenes = [x for x in agenes if x in selgenes]
scdata.adata = scdata.adata[:, filtgenes]


# Normalize by combat after subsetting genes:
# -------------------------------------------
usecombat = True
if usecombat:
    prefstr = prefstr + "_combat"
    sc.pp.combat(scdata.adata, key='celltype')

pltpref = 'figures/umap' + prefstr + '_ctalone'


def savelocal(fig, suffix, prefix=pltpref, dpi=350):
    """Save plot with current local prefix."""
    plotname = prefix + "_" + suffix + ".png"
    plt.tight_layout()
    plt.savefig(plotname, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(plotname)


# Compute UMAP from PCA, using BB k-NN for batch correction:
# ----------------------------------------------------------
sc.tl.pca(scdata.adata)
sc.external.pp.bbknn(scdata.adata, batch_key='celltype')
sc.tl.umap(scdata.adata, maxiter=None, random_state=1)


# Plot the UMAP for the major celltype:
# -------------------------------------
sc.pl.umap(scdata.adata, color='celltype', frameon=False,
           save=prefstr + '_ctalone.png')
sc.pl.umap(scdata.adata, color='projid', frameon=False,
           save=prefstr + '_ctalone_projid.png')
sc.pl.umap(scdata.adata, color='region', frameon=False,
           save=prefstr + '_ctalone_region.png')

fig = plot_umap_category_grid(scdata.adata, 'region',
                              scale=1, c=3, fs=20, s=.1)
savelocal(fig, 'byregion')

fig = plot_umap_category_grid(scdata.adata, 'celltype',
                              scale=1, c=8, fs=8, s=.1)
savelocal(fig, 'byct')


# Plot AD ascertainment variables:
# --------------------------------
sc.pl.umap(scdata.adata, color='nft', frameon=False, color_map='viridis',
           save=prefstr + '_ctalone_nft.png')
sc.pl.umap(scdata.adata, color='plaq_n', frameon=False, color_map='viridis',
           save=prefstr + '_ctalone_plaq_n.png')
sc.pl.umap(scdata.adata, color='plaq_d', frameon=False, color_map='viridis',
           save=prefstr + '_ctalone_plaq_d.png')
scdata.adata.uns['cogdxad_colors'] = ['red', 'lightgrey']
scdata.adata.uns['nrad_colors'] = ['red', 'lightgrey']
sc.pl.umap(scdata.adata, color='nrad', frameon=False,
           save=prefstr + '_ctalone_nrad.png')
sc.pl.umap(scdata.adata, color='cogdxad', frameon=False,
           save=prefstr + '_ctalone_cogdxad.png')

# Plot by density (takes a long time, better in R):
# sc.tl.embedding_density(scdata.adata, basis='umap', groupby='cogdxad')
# sc.pl.embedding_density(
#     scdata.adata, basis='umap', key='umap_density_cogdxad',
#     group='AD', save=prefstr + '_ctalone_cogdxad_density.png')


# Write UMAP coordinates to file:
# -------------------------------
scdata.adata.obs['C1'] = scdata.adata.obsm['X_umap'][:, 0]
scdata.adata.obs['C2'] = scdata.adata.obsm['X_umap'][:, 1]
adf = pd.DataFrame(scdata.adata.obs)
adf.to_csv('multiregion' + prefstr + '_subsetUMAP.tsv', sep="\t")


# Specific genes for each cell type:
# ----------------------------------
genes = ['CLTB', 'AP2M1', 'NRGN', 'PRNP', 'PKM', 'CALR', 'FAIM2', 'DHFR',
         'SLC27A6', 'YPEL1', 'CRK', 'NCAM2', 'ABCA7', 'UBTF', 'MT-CO3',
         'MT-ND3', 'PLCG2', 'CLU', 'NEUROD2', 'NOXA1', 'KAZN', 'PSEN1', 'APP']
genes = [x for x in genes if x in filtgenes]

sc.pl.umap(scdata.adata, color=genes, frameon=False, color_map='viridis',
           save=prefstr + '_ctalone_genes.png', vmax=2.5)

sc.pl.umap(scdata.adata, color=genes, frameon=False, color_map='viridis',
           save=prefstr + '_ctalone_genes_full.png')

sc.pl.matrixplot(scdata.adata, genes, groupby='region',
                 save=prefstr + '_ctalone_region.png')
sc.pl.matrixplot(scdata.adata, genes, groupby='celltype',
                 save=prefstr + '_ctalone_celltype.png')


# Leiden and markers:
# -------------------
sc.tl.leiden(scdata.adata)
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


# Use leiden + PAGA to estimate a diffusion pseudotime:
# -----------------------------------------------------
sc.tl.paga(scdata.adata, groups='leiden')
sc.pl.paga(scdata.adata, color='leiden',
           save=prefstr + '_ctalone_paga.png')

scdata.adata.uns['iroot'] = np.flatnonzero(
    scdata.adata.obs['leiden'] == '4')[0]

sc.tl.dpt(scdata.adata, n_branchings=1)

scdata.adata.obs['dpg'] = scdata.adata.obs.dpt_groups.astype('category')
scdata.adata.uns['dpg_colors'] = scdata.adata.uns['leiden_colors'][0:4]
sc.pl.umap(scdata.adata, color=['leiden', 'dpt_pseudotime', 'dpg'],
           frameon=False, legend_loc='on data',
           save=prefstr + '_ctalone_dpt.png')


# Write UMAP coordinates to file:
# -------------------------------
scdata.adata.obs['C1'] = scdata.adata.obsm['X_umap'][:, 0]
scdata.adata.obs['C2'] = scdata.adata.obsm['X_umap'][:, 1]
adf = pd.DataFrame(scdata.adata.obs)
adf.to_csv('multiregion' + prefstr + '_subsetUMAP.tsv', sep="\t")


# # TODO ADAPT:
# _, axs = pl.subplots(
#     ncols=3, figsize=(
#         6, 2.5), gridspec_kw={
#             'wspace': 0.05, 'left': 0.12})
# pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
# for ipath, (descr, path) in enumerate(paths):
#     _, data = sc.pl.paga_path(
#         adata, path, gene_names,
#         show_node_names=False,
#         ax=axs[ipath],
#         ytick_fontsize=12,
#         left_margin=0.15,
#         n_avg=50,
#         annotations=['distance'],
#         show_yticks=True if ipath == 0 else False,
#         show_colorbar=False,
#         color_map='Greys',
#         groups_key='clusters',
#         color_maps_annotations={'distance': 'viridis'},
#         title='{} path'.format(descr),
#         return_data=True,
#         show=False)
#     data.to_csv('./write/paga_path_{}.csv'.format(descr))

# pl.savefig('./figures/paga_path_paul15.pdf')
# pl.show()


# Rough split of trajectory (particular to "testedgenes")
# -------------------------------------------------------
# Remove points in single ind:
cmat = pd.crosstab(scdata.adata.obs['projid'], scdata.adata.obs['leiden'])
cmarg = cmat.max(0) / cmat.sum(0)
cmarg = cmarg.reset_index()
cmarg.columns = ['leiden', 'pct']
cutoff = 0.15
keepleiden = cmarg.loc[cmarg.pct < cutoff, 'leiden'].to_numpy()
keepleiden = keepleiden[np.isin(keepleiden, '5', invert=True)]
flagind = np.isin(scdata.adata.obs.leiden, keepleiden, invert=True)
keepind = np.isin(scdata.adata.obs.leiden, keepleiden)

# Split the groups:
C1 = scdata.adata.obs['C1']
C2 = scdata.adata.obs['C2']
C1med = np.median(C1[keepind])
C2y = (C1 - np.median(C1)) * 4 + np.median(C2)
# scdata.adata.obs['group'] = ['top' if C2[i] >= C2y[i]
scdata.adata.obs['group'] = ['top' if C1[i] <= C1med
                             else 'bottom' for i in range(len(C2))]
scdata.adata.obs.loc[flagind, 'group'] = ''

# Plot these groups:
sc.pl.umap(scdata.adata, color=['C2', 'group'],
           frameon=False, save=prefstr + '_ctalone_group.png')


# Save adata for GAM fits in separate analysis:
h5ad_file = re.sub("formodules", "fordetraj", scdata.h5ad_file)
scdata.adata.write(h5ad_file)


# Compute fits for each side of the trajectory:
# ----------------------------------------------
dgenes = ['KCNIP4','NRXN1', 'MT-ND3', 'MT-CO3', 'CALM1', 'CALM3']
X = scdata.adata.obs['C2'].to_numpy()
X = 1 - X / np.max(X)

# indset =
NG = len(dgenes)
plt.figure()
fig, axs = plt.subplots(2, NG, figsize=(4 * NG, 4))
for i, axvec in enumerate(axs):
    group = ['top', 'bottom'][i]
    ind = scdata.adata.obs['group'] == group
    subX = X[ind]
    for j, ax in enumerate(axvec):
        gene = dgenes[j]
        y = scdata.adata[ind, gene].X.toarray().T[0]
        gam = LinearGAM(s(0))
        gam.gridsearch(subX[:, np.newaxis], y)
        # plotting
        XX = gam.generate_X_grid(term=0)
        ax.scatter(subX, y, alpha=0.05, edgecolors='none', marker='.')
        ax.plot(XX[:, 0], gam.predict(XX))
        ax.plot(XX[:, 0], gam.prediction_intervals(XX, width=.95),
                c='r', ls='--')
        ax.set_title(gene)
        if j == 0:
            ax.set_ylabel(group)

savelocal(fig, 'selgenes_linearGAM')

