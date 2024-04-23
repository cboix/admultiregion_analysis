#!usr/bin/python
"""Compute co-expression modules for Tabula Sapiens using scdemon."""
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# on the Tabula Sapiens dataset
# Created: 11/30/23
# -------------------------------------------------------------
import os
import re
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import scdemon as sm
from scdemon.auxiliary import recipe_full

import socket
domain = socket.getfqdn()
if 'broadinstitute.org' in domain or 'private' in domain:
    topdir = '/home/cboix/data/'
else:
    topdir = '/home/cboix/'


# Arguments:
# dataset = 'COVID19_Cell2021'
dataset = 'TabulaSapiens'
dataset = 'EasySci_2023'
# dataset = 'Mathys_Nature2019'
z = 4.5
filt = 0.05
res = 2.5

if dataset == 'TabulaSapiens':
    filt = 0.1 # Try for lower genes

# Directories:
sdbdir = topdir + 'DEVTRAJ/db/multiRegion/'
extdir = sdbdir + 'external_datasets/'
datadir = extdir + dataset + '/'
h5ad_file = datadir + dataset + '.h5ad'
imgdir = topdir + 'DEVTRAJ/img/multiRegion/modules/'
imgpref = imgdir + 'external_' + dataset + '_'

logging.basicConfig(level=logging.DEBUG)


# Load in dataset:
# ----------------
if not os.path.exists(h5ad_file):
    if dataset == 'EasySci_2023':
        gsm = datadir + 'GSM6657986'
        adata = sc.read_mtx(gsm + '_gene_count.mtx.gz')
        odf = pd.read_csv(gsm + '_cell_annotation.csv.gz')
        vdf = pd.read_csv(gsm + '_gene_annotation.csv.gz')
        d = {'X': adata.X.T, 'obs': odf, 'var': vdf}
        adata = anndata.AnnData(**d)
        adata.var_names = adata.var.gene_short_name
        adata.obs_names = adata.obs.Cell_ID
    if dataset == 'Mathys_Nature2019':
        adata = sc.read_mtx(datadir + 'matrix.mtx.gz')
        odf = pd.read_csv(datadir + 'barcodes.tsv.gz', sep='\t')
        vdf = pd.read_csv(datadir + 'features.tsv.gz', header=None)
        d = {'X': adata.X.T, 'obs': odf, 'var': vdf}
        adata = anndata.AnnData(**d)
        adata.var.columns = ['gene']
        adata.var_names = adata.var['gene']
        adata.obs_names = adata.obs.TAG
        for key in ['projid','pre.cluster']:
            adata.obs[key] = adata.obs[key].astype('category')
    # Basic preprocessing needed for these datasets:
    adata.obs.index.name = 'index'
    adata.var.index.name = 'index'
    adata.var_names_make_unique()
    recipe_full(adata, preprocess=True, annotate=True)
    adata.write_h5ad(h5ad_file, compression='gzip')


adata = sc.read(h5ad_file)
# TODO: need to norm?

# To plot these on a grid of UMAPs later:
if dataset == 'COVID19_Cell2021':
    mod.adata.obsm['X_umap'] = mod.adata.obsm['X_tsne']


# Set up the modules computation object:
# --------------------------------------
max_k = 120
csuff = dataset + "_z" + str(z)
del(adata.obsm['X_pca']) # Let fbpca calculate the decomposition
mod = sm.scdemon(adata, h5ad_file=h5ad_file,
                  csuff=csuff, imgdir=imgdir, seed=1,
                  svd_k=max_k, filter_expr=filt,
                  z=z, calc_raw=False)
mod.setup()
kept_genes = mod.genes  # For bootstraps


# Make graphs from base and merged graph:
# ---------------------------------------
power_list = list(np.arange(0, 1, .25)) + [1]
top_k = 100
graph_id = 'base'
k_ind = np.arange(1, top_k)
mod.make_svd_subset_graph(graph_id, k_ind=k_ind, resolution=res)


graph_id = 'merge'
mod.make_merged_graph(graph_id, power_list=power_list, resolution=res)
# mod.plot_graph(graph_id, attr="leiden", show_labels=True, frac_labels=.5, width=16)


# Plots for these graphs:
# -----------------------
for graph_id in ['base']: # , 'merge']: # , 'boot']:
    # Get the modules/print out:
    mlist = mod.get_modules(graph_id, print_modules=False)
    mod.save_modules(graph_id, filedir=datadir)
    mod.plot_graph(graph_id, attr="leiden", show_labels=False, width=16)
    mod.plot_gene_umap(graph_id, width=16)
    if 'X_umap' in mod.adata.obsm:
        mod.plot_umap_grid(graph_id)
    # Table of module assignments:
    suffix = mod.csuff + '_' + graph_id
    mdf = mod.get_module_assignment(graph_id)
    mdf.to_csv(datadir + 'allgene_assignments_' + suffix + '.tsv',
               sep="\t", index=False)

# mod.plot_graph('base', attr="leiden", show_labels=True, width=16)

# Write representations and scores to file:
# -----------------------------------------
for graph_id in ['base']: #, 'merge']: # , 'boot']:
    g = mod.graphs[graph_id]
    layout = np.array(g.layout.coords)
    gdf = pd.DataFrame({'gene': g.graph.vs['name'],
                        'leiden': g.assign['leiden'],
                        'col': g.colors['leiden'],
                        'L1': layout[:, 0], 'L2': layout[:, 1],
                        'U1': g.umat[:, 0], 'U2': g.umat[:, 1]})
    edges = np.array(g.graph.get_edgelist())
    gdf.to_csv(datadir + 'graph_nodes_info_' +
               mod.csuff + '_' + graph_id + '.tsv', sep="\t")
    np.savetxt(datadir + 'graph_edgelist_' +
               mod.csuff + '_' + graph_id + '.tsv', X=edges)
    # Save the scores as well (can be quite large)
    scdf = pd.DataFrame({'bc': mod.adata.obs_names.to_numpy()})
    for i in range(g.scores['leiden'].shape[1]):
        scdf[f'M{i}'] = g.scores['leiden'][:, i]
    scdf.to_csv(datadir + 'module_cell_scores_' +
                mod.csuff + '_' + graph_id + '.tsv.gz', sep="\t")
    names = mod.graphs[graph_id].mnames['leiden']
    namefile = datadir + 'module_names_' + mod.csuff + '_' + graph_id + '.tsv'
    with open(namefile, 'w') as f:
        for name in names:
            f.write(name + '\n')

# Other plots for graphs:
# -----------------------
cvlist = ['projid', 'celltype', 'region',
          'braaksc', 'cogdx', 'Apoe_e4', 'msex']

if dataset == 'TabulaSapiens':
    cvlist = ['organ_tissue', 'method', 'donor', 'anatomical_information',
            'cell_ontology_class', 'free_annotation', 'compartment','gender']
elif dataset == 'EasySci_2023':
    cvlist = ['Region', 'Cell_type', 'Condition']
elif dataset == 'Mathys_Nature2019':
    cvlist = ['projid', 'broad.cell.type', 'Subcluster']
elif dataset == 'COVID19_Cell2021':
    cvlist = ['majorType', 'celltype', 'City',  'Sex', 'Sample type',
              'CoVID-19 severity', 'Sample time',
              'SARS-CoV-2',
              'Single cell sequencing platform',
              'BCR single cell sequencing', 'TCR single cell sequencing',
              'Outcome', 'Comorbidities',
              'COVID-19-related medication and anti-microbials']
              # 'Sampling day (Days after symptom onset)', 'Age',
              # 'Leukocytes [G/L]', 'Neutrophils [G/L]', 'Lymphocytes [G/L]']


for graph_id in ['base', 'merge']:
    mod.plot_heatmap_avgexpr(graph_id, cvlist=cvlist)




# Get counts to calculate enrichments for covariates:
# ---------------------------------------------------
# Average expression for base graph:
graph_id = 'base'
smat = mod.graphs[graph_id].scores['leiden']

msd = np.std(smat, axis=0)
zmat = 1 * ((smat / msd[np.newaxis, :]) > 2.5)
nmod = np.sum(zmat, axis=0)

# NB: Aggregate over individuals + cell types?

hgdf = None
for cvar in cvlist:
    print(cvar)
    for cat in adata.obs[cvar].cat.categories:
        print('--', cat)
        ind = (adata.obs[cvar] == cat)
        cdf = pd.DataFrame({
            'var': cvar,
            'cat': cat,
            'module': range(len(nmod)),
            'ncat': np.sum(ind),
            'nmod': nmod,
            'nint': np.sum(zmat[ind,:], axis=0),
            'ntot': zmat.shape[0]})
        if hgdf is None:
            hgdf = cdf
        else:
            hgdf = pd.concat([hgdf, cdf])

# Write data out:
hg_data_csv = datadir + 'modules_hgdata_' + mod.csuff + '_' + graph_id + '.tsv'
hgdf.to_csv(hg_data_csv, sep="\t")


# mod.calc_svd_corr(cvlist)
# mod.plot_svd_corr(cvlist)
