#!usr/bin/python
"""Compute co-expression modules for multiRegion data using scdemon."""
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# Using data loader + modules + graph handlers
# Created: 06/16/21 // Last updated: 11/27/21
# python $SBINDIR/modules/modules_02_basic_graphs.py --ct Opc
# -------------------------------------------------------------
import logging
import os
import re
# import scanpy as sc
import pandas as pd
import numpy as np
import scdemon as sm
from adata_handling import adata_loader
import argparse

import socket
domain = socket.getfqdn()
if 'broadinstitute.org' in domain or 'private' in domain:
    topdir = '/home/cboix/data/'
else:
    topdir = '/home/cboix/'

dbdir = topdir + 'DEVTRAJ/db/'
datadir = dbdir + 'multiRegion/'
resultsdir = dbdir + 'multiRegion/modules/'
bootdir = resultsdir + 'bootstrap/'
dedir = datadir + 'dereg/'
imgdir = topdir + 'DEVTRAJ/img/multiRegion/modules/'


# Build argparser so that we can run this from CLI / on batched:
# --------------------------------------------------------------
parser = argparse.ArgumentParser(description='Calculate scdemon modules.')
parser.add_argument('--ct', metavar='ct', type=str,
                    help='celltype')
parser.add_argument('--st', metavar='st', type=str,
                    help='subtype', default=None)
parser.add_argument('--z', metavar='z', type=float,
                    help='zscore', default=4.5)
parser.add_argument('--filt', metavar='filt', type=float,
                    help='filter low expressed genes', default=0.05)
parser.add_argument('--res', metavar='res', type=float,
                    help='leiden clustering resolution', default=2.5)
args = parser.parse_args()

# Set logging level (TODO: set by argparser)
logging.basicConfig(level=logging.INFO)
logging.info(str(args))


# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
scdata = adata_loader(celltype=args.ct, subtype=args.st,
                      dbdir=dbdir,
                      normalize=True,
                      usecombat=False, usehvg=False,
                      remove_mt=False, filter_TH=False)
scdata.load()
# scdata.compute_representation()  # For later
# scdata.save_adata()


# Set up the modules computation object:
# --------------------------------------
max_k = 200
csuff = scdata.csuff + "_z" + str(args.z)
mod = sm.scdemon(scdata.adata, h5ad_file=scdata.h5ad_file,
                  csuff=csuff, imgdir=imgdir, seed=1,
                  svd_k=max_k, filter_expr=args.filt,
                  z=args.z, calc_raw=False)
mod.setup()
kept_genes = mod.genes  # For bootstraps


# Load best hyperparameter for k:
# -------------------------------
paramfile = resultsdir + 'paramsearch_' + mod.csuff + '.tsv'
paramdf = pd.read_csv(paramfile, sep="\t")
logging.info(paramdf.head())
ng = paramdf.ng.to_numpy()
nm = paramdf.nm.to_numpy()
ki = np.argmax(paramdf.score.to_numpy())
top_k = paramdf.k.to_numpy()[ki]
logging.info(f"Setting k={top_k}, ng={ng[ki]}, nm={nm[ki]}")


# Make graphs from precomputed bootstrapped zscore matrix estimates:
# ------------------------------------------------------------------
nperm = 10
power_list = list(np.arange(0, 1, .25)) + [1]

for power in power_list:
    ztotfile = bootdir + f'ztot_n{nperm}_p{power}' + mod.csuff + '.npz'
    ztot = np.load(ztotfile, allow_pickle=True)['ztot']
    graph_id = f'p{power}'
    # Giving the zscores so using a raw cutoff here
    mod.make_graph(graph_id, corr=ztot, cutoff=mod.z,
                   use_zscore=False, full_graph_only=True)


# Make graphs from the bootstraps and from base:
# ----------------------------------------------
graph_id = 'base'
k_ind = np.arange(1, top_k)
mod.make_svd_subset_graph(graph_id, k_ind=k_ind, resolution=args.res)
# mod.make_merged_graph(graph_id, power_list=power_list, resolution=args.res)
# mod.plot_graph(graph_id, attr="leiden", show_labels=True, frac_labels=.5, width=16)


# TODO: specify graphs so that we can do merged on both base and bootstrapped
graph_id = 'boot'
mod.make_merged_graph(graph_id, power_list=power_list, resolution=args.res)
mlist = mod.get_modules(graph_id, print_modules=True)


# Plots for these graphs:
# -----------------------
for graph_id in ['base', 'boot']:
    # Get the modules/print out:
    mlist = mod.get_modules(graph_id, print_modules=False)
    mod.save_modules(graph_id, filedir=resultsdir)
    mod.plot_graph(graph_id, attr="leiden", show_labels=False, width=16)
    mod.plot_gene_umap(graph_id, width=16)
    if 'X_umap' in mod.adata.obsm:
        mod.plot_umap_grid(graph_id)
    # Table of module assignments:
    suffix = mod.csuff + '_' + graph_id
    mdf = mod.get_module_assignment(graph_id)
    mdf.to_csv(resultsdir + 'allgene_assignments_' + suffix + '.tsv',
               sep="\t", index=False)


# Write representations and scores to file:
# -----------------------------------------
for graph_id in ['base', 'boot']:
    g = mod.graphs[graph_id]
    layout = np.array(g.layout.coords)
    gdf = pd.DataFrame({'gene': g.graph.vs['name'],
                        'leiden': g.assign['leiden'],
                        'col': g.colors['leiden'],
                        'L1': layout[:, 0], 'L2': layout[:, 1],
                        'U1': g.umat[:, 0], 'U2': g.umat[:, 1]})
    edges = np.array(g.graph.get_edgelist())
    gdf.to_csv(resultsdir + 'graph_nodes_info_' +
               mod.csuff + '_' + graph_id + '.tsv', sep="\t")
    np.savetxt(resultsdir + 'graph_edgelist_' +
               mod.csuff + '_' + graph_id + '.tsv', X=edges)
    # Save the scores as well (can be quite large)
    scdf = pd.DataFrame({'bc': mod.adata.obs_names.to_numpy()})
    for i in range(g.scores['leiden'].shape[1]):
        scdf[f'M{i}'] = g.scores['leiden'][:, i]
    scdf.to_csv(resultsdir + 'module_cell_scores_' +
                mod.csuff + '_' + graph_id + '.tsv.gz', sep="\t")


# Module enrichment in any ranked gene list / DEG list:
# ----------------------------------------------------
def get_dflist(depref, pathlist):
    """Get the DE datasets for given prefix and test variables."""
    dflist = {}
    for path in pathlist:
        # Load and format data:
        defile = dedir + depref + path + ".merged.tsv.gz"
        if os.path.exists(defile):
            df = pd.read_csv(defile, sep="\t")
            logging.info(str(defile))
            df['gset'] = '--'
            df.loc[df.col_nm == 2, 'gset'] = 'Up'
            df.loc[df.col_nm == 1, 'gset'] = 'Down'
            dflist[path] = df
    return(dflist)


def get_de_prefix(ct, st):
    """Get the DE prefixes for given cell + sub types."""
    method = 'allmethods'
    reg = 'allregions'
    ctpref = ct + "_"
    if ct == 'Vasc_Epithelia':
        st = list(pd.unique(mod.adata.obs.celltype))
    elif ct == 'Mic_Immune':
        st = ['Mic', 'CAM', 'T']
    elif st in ['ECneurons', 'HCneurons', 'THneurons']:
        reg = st[0:2]
        st = [re.sub(" ", "_", x) for x in pd.unique(mod.adata.obs.celltype)]
    elif st == 'CTXneurons':
        st = [re.sub(" ", "_", x) for x in pd.unique(mod.adata.obs.celltype)]
        reg = "neocortex"
    else:
        st = st
    desuff = "_" + reg + "_"
    depref = method + "."
    if st is None:
        depref = [depref + ctpref + ct + desuff]
    elif isinstance(st, list):
        depref = [depref + ctpref + s + desuff for s in st]
    else:
        depref = [depref + ctpref + st + desuff]
    return(depref)


def run_regs(dflist, prefstr, graph_id):
    """Run and save comparisons with the regressions."""
    suffix = 'allpath_' + prefstr + mod.csuff + "_" + graph_id
    # TODO: Don't recalculate - not a big deal in terms of speed though
    statsdf = mod.plot_df_enrichment(
        dflist, 'gset', graph_id=graph_id, suffix=suffix)
    statsdf = mod.plot_df_enrichment(
        dflist, 'gset', graph_id, suffix=suffix, ext='pdf')
    # Save table of enrichment statistics, with short module names:
    mnames = np.array(mod.graphs[graph_id].mnames['leiden'])
    statsdf['mname'] = mnames[statsdf.module.to_numpy()]
    statsdf.to_csv(resultsdir + 'enrstats_' + suffix + '.tsv',
                   sep="\t", index=False)


depref = get_de_prefix(args.ct, args.st)
for prefstr in depref:
    pathlist = ['nrad', 'nft', 'plaq_n', 'plaq_d', 'cogdxad']
    dflist = get_dflist(prefstr, pathlist)
    if dflist != {}:
        # Save dflist:
        alldegdf = pd.concat(dflist)
        alldegdf = alldegdf.reset_index()
        alldegdf.to_csv(resultsdir + 'allpath_' + prefstr + "usedDEGs.tsv",
                        sep="\t")
        for graph_id in ['base', 'boot']:
            run_regs(dflist, prefstr, graph_id)


# Other plots for graphs:
# -----------------------
cvlist = ['projid', 'celltype', 'region',
          'braaksc', 'cogdx', 'Apoe_e4', 'msex']
mod.calc_svd_corr(cvlist)
mod.plot_svd_corr(cvlist)

# TODO: FIX category:
scdata.adata.obs['braaksc'] = scdata.adata.obs['braaksc'].astype('category')
scdata.adata.obs['cogdx'] = scdata.adata.obs['cogdx'].astype('category')
for graph_id in ['boot', 'base']:
    mod.plot_heatmap_avgexpr(graph_id, cvlist=cvlist)
