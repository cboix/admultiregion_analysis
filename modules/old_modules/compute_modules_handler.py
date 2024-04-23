# --------------------------------------------
# Class for computing modules on adata object:
# Updated: 06/13/21
# OUTDATED, use SCMODULE package instead
# --------------------------------------------
import glob
import h5py
import re
import gzip
import numpy as np
import pandas as pd
import time
import gc
import os
import sys

# For graphs
import umap
import fbpca
import igraph
from scipy import sparse
import leidenalg as la
from adjustText import adjust_text

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import Normalize, rgb2hex
import seaborn as sns

# For PCA:
import scanpy as sc
# For GO annotations (may not work on luria private node)
from gprofiler import GProfiler

# For the correlation cutoffs:
from scipy.interpolate import SmoothBivariateSpline
import scipy.stats as st

# Vectorized correlation, matrix to vector:
def vcorrcoef(X,y):
    Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm)*(y-ym),axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)*np.sum((y-ym)**2))
    r = r_num/r_den
    return r


# External / auxiliary functions for graph operations:
def prune_scale(aw, cutoff=0.90, rcut=1):
    # TODO: handle 0s for max:
    awm1 = np.max(aw, axis=1).T.toarray()[0]
    awm2 = np.max(aw, axis=0).toarray()[0]
    awm = np.vstack((awm1, awm2))
    awm = np.max(awm, axis=0)
    aw = aw.tocoo()
    scale = np.vstack((awm[aw.row], awm[aw.col]))
    scale = np.max(scale, axis=0)
    pct_data = aw.data / scale
    kind = pct_data > cutoff
    aw = sparse.coo_matrix((aw.data[kind], (aw.row[kind], aw.col[kind])), shape=aw.shape)
    aw = aw.tocsr()
    # Remove nodes with no/very few links:
    rind = np.array((np.sum(aw, axis=1) > rcut).T)[0]
    cind = np.array((np.sum(aw, axis=0) > rcut))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return(aw, ind)


def prune_knn(aw, k=50, twodir=True, row_only=False):
    # Prune to maximum k for each node:
    aw = aw.tocoo()
    ind = np.argsort(-aw.data) # Sort correlations
    mk = np.zeros(aw.shape[0], int) # Kept margin
    keepind = np.zeros(len(ind), int)
    # Add indices in order of strength:
    for i in range(len(ind)):
        r = aw.row[i]
        c = aw.col[i]
        if twodir:
            cond = mk[r] < k or mk[c] < k
        else:
            cond = mk[r] < k and mk[c] < k
        if cond:
            mk[r] += 1
            if twodir or not row_only:
                mk[c] += 1
            keepind[i] = 1
    # Remove edges:
    kind = ind[keepind == 1]
    aw = sparse.coo_matrix((aw.data[kind],
                            (aw.row[kind], aw.col[kind])),
                           shape=aw.shape)
    aw = aw.tocsr()
    # Remove nodes with no links left (shouldn't be any)
    rind = np.array((np.sum(aw, axis=1) > 0).T)[0]
    cind = np.array((np.sum(aw, axis=0) > 0))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return(aw, ind)


def adj_from_bivariate_cutoff(corr, gmarg, z=None, min_sd=0.01, rcut=1):
    print("Converting to i,j,v df")
    corr = corr - np.diag(np.diag(corr))
    # Set a z-threshold:
    if z is None:
        pcut = 0.05 # TODO: get proper BY testing correction equation
        z = - st.norm.ppf(pcut / corr.shape[0]) # If z is none.
    # Bin genes by their margin:
    lmarg = np.log10(np.array(gmarg))
    lmarg[lmarg < -3] = -3 # TODO: can lower, but not much.
    xs = np.linspace(np.min([np.min(lmarg), -3]), 0, 50)
    dmarg = np.digitize(lmarg, xs)
    # Calculate the means/sdevs:
    u = np.unique(dmarg)
    NS = len(u)
    smean = np.zeros((NS, NS))
    ssd = np.zeros((NS, NS))
    for i,xi in enumerate(u):
        rind = (dmarg == xi)
        ci = corr[rind,:]
        for j,xj in enumerate(u):
            cind = (dmarg == xj)
            cij = ci[:,cind]
            smean[i,j] = np.mean(cij)
            ssd[i,j] = np.std(cij)
    smean = sparse.coo_matrix(smean)
    ssd = sparse.coo_matrix(ssd)
    # Fit the spline to the bivariate data:
    print("Fitting + predicting spline")
    splmean = SmoothBivariateSpline(x=u[smean.row], y=u[smean.col], z=smean.data, kx=2, ky=2)
    splsd = SmoothBivariateSpline(x=u[ssd.row], y=u[ssd.col], z=ssd.data, kx=2, ky=2)
    xloc = np.arange(len(xs))
    smean = splmean(xloc, xloc)
    ssd = splsd(xloc, xloc)
    ssd[ssd < min_sd] = min_sd
    # Calculate the z-scored cutoff based on spline predictions:
    zcut = smean + ssd * z
    zcoo = sparse.coo_matrix(zcut)
    zdf = pd.DataFrame({'x':xloc[zcoo.row], 'y':xloc[zcoo.col], 'cutoff':zcoo.data})
    # Merge cutoffs against pre-cut:
    print("Thresholding adjacency matrix")
    min_cut = np.min(zcoo.data)
    zcorr = corr * (corr >= min_cut)
    aw = sparse.coo_matrix(zcorr) # TODO: Bypass this creation
    aind = aw.data >= min_cut
    awcut_df = pd.DataFrame({'x':dmarg[aw.row], 'y':dmarg[aw.col],
        'r': aw.row, 'c': aw.col, 'dat':aw.data})
    awcut_df = awcut_df.merge(zdf)
    awcut_df = awcut_df[awcut_df.dat >= awcut_df.cutoff]
    # Create thresholded adjacency matrix:
    aw2 = sparse.coo_matrix((awcut_df.dat, (awcut_df.r, awcut_df.c)), shape=aw.shape)
    aw = aw2.tocsr()
    rind = np.array((np.sum(aw2, axis=1) > rcut).T)[0]
    cind = np.array((np.sum(aw2, axis=0) > rcut))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return(aw, ind)


# Visual style dictionary for plotting graphs:
# TODO: Put into graph object
visual_style = {}
visual_style["vertex_size"] = 10
visual_style["vertex_color"] = 'slateblue'
visual_style["vertex_frame_color"] = 'slateblue'
visual_style["vertex_label"] = ''
visual_style["edge_color"] = 'lightgrey'
visual_style["bbox"] = (2400, 2400)
visual_style["margin"] = 10
visual_style["edge_width"] = 0.5


class gene_graph(object):
    def __init__(self, corr, cutoff, suffix, imgdir=None,
                 edge_weight=None, labels=None, z=None, gmarg=None,
                 layout_method='fr', k=None, scale=None):
        self.corr = corr
        self.cutoff = cutoff
        self.suffix = suffix
        self.labels = labels
        self.z = z # For adjacency cutoff
        self.gmarg = gmarg # Margin of X
        self.k = k
        self.scale = scale
        self.edge_weight = edge_weight
        self.layout_method = layout_method
        if imgdir is not None:
            self.imgdir = imgdir
        else:
            self.imgdir = './'
        # For graph colors:
        self.dbdir = '/home/cboix/data/DEVTRAJ/db'
        self.anndir = self.dbdir + '/Annotation/'
        self.snapcols = pd.read_csv(
            self.anndir + '/snap_colors.tsv', header=None)
        # For covariate annotations
        self.assign = {}
        self.colors = {}

    def build_graph(self, resolution=None):
        self.compute_adjacency()
        if self.scale is not None:
            self.prune_by_scale()
        if self.k is not None:
            self.prune_by_knn()
        self.make_graph_object()
        self.layout_graph()
        self.detect_modules(resolution=resolution)

    # TODO: Find better param than rcut
    def compute_adjacency(self, rcut=2):
        if self.z is None:
            self.aw = self.corr * (self.corr > self.cutoff)
            self.aw = self.aw - np.diag(np.diag(self.aw))
        else:
            self.aw, ind = adj_from_bivariate_cutoff(
                    self.corr, self.gmarg, z=self.z)
            self.labels = self.labels[ind]
        self.aw = sparse.csr_matrix(self.aw)
        rind = np.array((np.sum(self.aw, axis=0) > rcut))[0]
        # Remove nodes with no links:
        self.aw = self.aw[rind, :]
        self.aw = self.aw[:, rind]
        self.labels = self.labels[rind]

    def prune_by_scale(self):
        print("Pruning by relative cutoff, scale=" + str(self.scale))
        self.aw, keptind = prune_scale(self.aw, cutoff=self.scale)
        if self.labels is not None:
            self.labels = self.labels[keptind]

    def prune_by_knn(self):
        print("Pruning by KNN, K=" + str(self.k))
        self.aw, keptind = prune_knn(self.aw, k=self.k)
        if self.labels is not None:
            self.labels = self.labels[keptind]

    def make_graph_object(self, del_nodes=True):
        print("Making graph object")
        self.aw = self.aw.tocoo()
        ind = self.aw.row < self.aw.col
        edg_list = list(tuple(zip(self.aw.row[ind], self.aw.col[ind])))
        edg_wt = self.aw.data[ind]
        self.gw = igraph.Graph(edg_list, directed=False)
        self.gw.vs['name'] = self.labels
        # Delete after labeling:
        if del_nodes:
            todel = [i for i,x in enumerate(self.gw.degree()) if x < 1]
            self.gw.delete_vertices(todel)
        if self.edge_weight is not None:
            self.gw.es['weight'] = self.edge_weight

    def layout_graph(self):
        print("Laying out graph, method=" + str(self.layout_method))
        if self.layout_method == 'fr':
            self.layout = self.gw.layout_fruchterman_reingold(niter=500, grid=False)
        else:
            self.layout = self.gw.layout(self.layout_method)
        # TODO: Fix the really distant cliques (?)
        print("Kept", len(self.layout), "nodes, out of", self.corr.shape[0])

    def detect_modules(self, method='leiden',
            resolution=None, partition_type=None,
            use_weights=False, n_iterations=-1, random_state=1):
        print("Running module detection")
        if method == 'leiden':
            partition_kwargs = {}
            if partition_type is None:
                if resolution is None:
                    partition_type = la.ModularityVertexPartition
                else:
                    partition_type = la.RBConfigurationVertexPartition
                    partition_kwargs['resolution_parameter'] = resolution
            if use_weights:
                partition_kwargs['weights'] = np.array(
                        self.gw.es['weight']).astype(np.float64)
            partition_kwargs['n_iterations'] = n_iterations
            partition_kwargs['seed'] = random_state
            # Run leidenalg:
            ptn = la.find_partition(self.gw, partition_type, **partition_kwargs)
        else:
            # TODO: implement louvain, resolution seems more tunable?
            print("Louvain, other methods not ready yet")
        # Colors:
        ptlist = np.array([np.array(x) for x in list(ptn)])
        nc = len(ptlist)
        print("Found " + str(nc) + " clusters")
        j = 1
        lcols = self.snapcols[j:j+nc].to_numpy().T[0]
        self.assign['leiden'] = (np.zeros(len(self.gw.vs)) - 1).astype(int)
        for i in range(nc):
            self.assign['leiden'][ptlist[i]] = i
        self.colors['leiden'] = lcols[self.assign['leiden']]

    # TODO: Improve method by connectivity / (cross-connect leiden)
    def select_labels(self, labels, frac_labels):
        frac_labels = 0 if frac_labels < 0 else frac_labels
        # Subset labels:
        if frac_labels < 1.0:
            la = self.assign['leiden']
            u,c = np.unique(la, return_counts=True)
            tot = np.round(np.sum(c) * frac_labels,0)
            # At least one per leiden + rest divided up:
            cr = c - 1
            cr[cr < 0] = 0
            alloc = c - cr
            # Remaining:
            r1 = np.sum(alloc)
            if tot - r1 > 0:
                r2 = tot - r1
                alloc += np.round(
                        cr * (r2/np.sum(cr))).astype(int)
            # Allocate:
            keeplabs = []
            for ll in u:
                ind = np.random.choice(np.where(la == ll)[0], alloc[ll])
                keeplabs = keeplabs + ind.tolist()
            subset_labels = [labels[i] if i in keeplabs else ''
                    for i in range(len(labels))]
            labels = subset_labels
        return(labels)

    # TODO: allow extra suffix for graph (e.g. if params are different)
    def plot_graph(self, attr=None, color=None,
            show_labels=True, frac_labels=1.0,
            adj_txt=False, w=24):
        # Set color + plotname:
        plotname = self.imgdir + "graph_"
        if attr is None:
            visual_style["vertex_color"] = 'slateblue'
            visual_style["vertex_frame_color"] = 'slateblue'
        else:
            if color is not None:
                self.colors[attr] = color
            visual_style["vertex_color"] = self.colors[attr]
            visual_style["vertex_frame_color"] = self.colors[attr]
            plotname += str(attr) + '_'
        plotname += self.suffix + '.png'
        # Set label specifications:
        labels = self.select_labels(self.gw.vs['name'], frac_labels)
        if adj_txt or not show_labels:
            visual_style["vertex_label"] = ''
        else:
            visual_style["vertex_label"] = labels
        # Plot graph on matplotlib axes:
        fig = plt.figure(figsize=(w,w))
        ax = plt.gca()
        igraph.plot(self.gw, layout=self.layout, target=ax, **visual_style)
        # Adjusted text:
        if adj_txt and show_labels:
            keptind = np.where(np.array(labels) != '')[0].tolist()
            if len(keptind) > 0:
                llist = list(self.layout)
                text = [plt.text(llist[i][0], llist[i][1], labels[i],
                    ha='center', va='center', fontsize=12)
                    for i in keptind]
                adjust_text(text, lim=25)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(plotname)
        plt.close()
        print("Plotted graph to " + plotname)


# Pseudotime method handler:
class handler_ZCA_module(object):
    def __init__(self, adata, csuff, filter_expr=None, z=None,
            get_raw=False, imgdir=None, use_fbpca=True, use_v=False):
        # Arguments:
        self.adata = adata
        self.genes = self.adata.var_names # labels
        self.csuff = csuff
        self.z = z
        self.filter_expr = filter_expr
        self.use_fbpca = use_fbpca
        self.use_v = use_v
        # Directories/Files:
        self.dbdir = '/home/cboix/data/DEVTRAJ/db'
        self.datadir = self.dbdir + '/multiRegion/'
        self.anndir = self.dbdir + '/Annotation/'
        if imgdir is not None:
            self.imgdir = imgdir
        else:
            self.imgdir = './'
        # Additional metadata / files:
        self.snapcols = pd.read_csv(
            self.anndir + '/snap_colors.tsv', header=None)
        self.snapcols = self.snapcols.to_numpy().T[0][1:].tolist()
        self.graphs = {} # Store computed graphs
        self.get_raw = get_raw
        self.gp = GProfiler(return_dataframe=True)
        self.modules = {}

    def setup(self):
        self.filter_adata()
        self.compute_pca()
        self.compute_decorr()
        self.compute_corr()

    def filter_adata(self):
        # Filter to top expr:
        if self.filter_expr is not None:
            gmarg = np.mean(self.adata.X > 0, axis=0)
            if sparse.issparse(self.adata.X):
                gmarg = np.array(gmarg)[0]
            filtgenes = self.adata.var_names[gmarg >= self.filter_expr]
            self.adata = self.adata[:, filtgenes]
            self.genes = self.adata.var_names
            print(self.adata.shape)
        self.gmarg = np.mean(self.adata.X > 0, axis=0) # For filt later
        if sparse.issparse(self.adata.X):
            self.gmarg = np.array(self.gmarg)[0]

    def compute_pca(self, k=100):
        if self.use_fbpca:
            print("Computing PCA with fbpca")
            self.U, self.s, self.V = fbpca.pca(self.adata.X, k=k, raw=False)
        else:
            if 'X_pca' not in self.adata.obsm.keys():
                print("Computing PCA through scanpy")
                sc.tl.pca(self.adata, n_comps=k)
                # TODO: Fix U in scanpy PCA --> U * S
            self.U = self.adata.obsm['X_pca'] # TODO: FIX
            self.s = self.adata.uns['pca']['variance']
            self.V = self.adata.varm['PCs'].T

    def compute_decorr(self):
        sinv = np.diag(1/self.s)
        if not self.use_v:
            XtU = self.adata.X.T.dot(self.U)
            self.Xw = XtU.dot(sinv.dot(self.U.T)) / self.V.shape[1]
        else:
            pass
            # XV = self.adata.X.dot(self.V.T)
            # self.Xw = (XV.dot(sinv.dot(self.V)) / self.U.shape[0]).T

    def compute_corr(self):
        if self.get_raw:
            self.corr_raw = np.corrcoef(self.adata.X.T)
        if self.use_v:
            # Decorrelated correlation directly from SVD:
            cv = self.V.T.dot(self.V) / self.U.shape[0]
            sd = np.sqrt(np.diag(cv))
            cv = cv / sd[:,np.newaxis]
            self.corr_w = cv / sd[np.newaxis,:]
        else:
            # Compute decorr:
            self.corr_w = np.corrcoef(self.Xw)

    def compute_umap(self):
        uw = umap.UMAP()
        self.umat = uw.fit_transform(self.corr_w)

    # Build the graph object and store in dictionary:
    def make_graph(self, suff, cutoff, corr=None, edge_weight=None,
            k=None, scale=None, z=None, layout_method='fr', resolution=None):
        if z is not None:
            self.z = z
        if corr is None:
            if suff == 'raw':
                usecorr = self.corr_raw
            else:
                usecorr = self.corr_w
        else:
            usecorr = corr
        self.graphs[suff] = gene_graph(usecorr, cutoff=cutoff,
                                       edge_weight=edge_weight,
                                       suffix=self.csuff + "_" + suff,
                                       layout_method=layout_method,
                                       z=self.z, gmarg=self.gmarg,
                                       labels=self.genes, k=k, scale=scale,
                                       imgdir=self.imgdir)
        self.graphs[suff].build_graph(resolution=resolution)

    def plot_graph(self, suff, attr=None, show_labels=None,
            frac_labels=1.0, adj_txt=False, w=20):
        if suff not in self.graphs:
            # TODO: Default to quantiles.
            self.make_graph(suff, cutoff=.4) # Default to this
        # Plot graph object, both ways:
        self.graphs[suff].plot_graph(
            attr=attr, show_labels=show_labels,
            frac_labels=frac_labels,
            adj_txt=adj_txt, w=w)

    def plot_gene_logfc(self, suff, attr='celltype', fc_cut=2, **kwargs):
        print("Running rank genes on", attr)
        sc.tl.rank_genes_groups(self.adata, attr)
        iscat = (self.adata.obs[attr].dtype.name == 'category')
        # TODO: How to deal with non-categorical
        # Color map:
        cmap = plt.cm.RdYlBu # TODO: allow other options
        norm = Normalize(vmin=-fc_cut, vmax=fc_cut)
        # Turn into assignment:
        narr = self.adata.uns['rank_genes_groups']['names']
        larr = self.adata.uns['rank_genes_groups']['logfoldchanges']
        parr = self.adata.uns['rank_genes_groups']['pvals_adj']
        names = list(narr.dtype.names)
        for name in names:
            namestr = attr + "_" + re.sub("[ /]", "_", name)
            print("Plotting for", namestr)
            # TODO: Pass the list of names + a score to the gene graph instead
            n1 = np.array(narr[name])
            xind = np.array([np.where(n1 == x)[0][0]
                             for x in self.graphs[suff].gw.vs['name']])
            n1 = n1[xind]
            l1 = np.array(larr[name])[xind]
            p1 = np.array(parr[name])[xind]
            l1 = l1 * (p1 < 0.05)
            l1[l1 < -fc_cut] = -fc_cut
            l1[l1 > fc_cut] = fc_cut
            # Color mapping:
            vcols = [rgb2hex(c) for c in cmap(norm(-l1))]
            # Plot graph:
            self.graphs[suff].plot_graph(attr=namestr, color=vcols, **kwargs)

    # TODO: add plotting scale parameters:
    def plot_heatmap_avgexpr(self, suff, cvlist=None,
            attr='leiden', cbar=False):
        self.modules[suff] = {}
        la = self.graphs[suff].assign[attr]
        lc = self.graphs[suff].colors[attr]
        nam = np.array(self.graphs[suff].gw.vs['name'])
        las = np.unique(la)
        lmat = np.zeros((self.U.shape[0], np.max(las)+1))
        issp = sparse.issparse(self.adata.X)
        for ll in las:
            lnam = np.sort(nam[np.where(la == ll)[0]])
            avgexpr = np.mean(self.adata[:,lnam].X, axis=1)
            if issp:
                avgexpr = np.array(avgexpr).T[0]
            lmat[:,ll] = avgexpr
            self.modules[suff][ll] = lnam
            # print(ll)
            # print(" ".join(lnam))
        if cvlist is None:
            cvlist = ['celltype', 'region', 'niareagansc',
                    'cogdx', 'Apoe_e4', 'msex']
        hratio = [1] + [len(pd.unique(self.adata.obs[x])) for x in cvlist]
        hfull = np.sum(hratio)
        # For subgraph colors heatmap:
        u, c = np.unique(la, return_index=True)
        cmat = u[np.newaxis,:]
        ccol = lc[c]
        ccmap = sns.color_palette(ccol)
        # Make heatmap figure:
        w = 7 * lmat.shape[1] / 20 + 1.2 + 0.8 * cbar
        h = 2 * hfull / 6 + 0.1 * len(cvlist)
        fig, axs = plt.subplots(len(cvlist) + 1, 1,
            figsize=(w, h),
            gridspec_kw={'height_ratios': hratio, 'hspace': 0.05,
                'left': 1.2 / w, 'right': 1 - 0.8 / w,
                'top':1 - 0.1 / h, 'bottom': .4 / h})
        sns.heatmap(cmat, yticklabels=False, xticklabels=False,
                cmap=ccmap, ax=axs[0], cbar=False)
        for i, covar in enumerate(cvlist):
            covar_dummy = pd.get_dummies(self.adata.obs[covar])
            covar_labels = covar_dummy.columns.tolist()
            covar_dummy = covar_dummy.to_numpy()
            covar_cols = np.sum(covar_dummy, axis=0)
            tform = covar_dummy / covar_cols[np.newaxis, :]
            avg_lmat = tform.T.dot(lmat)
            scaled_lmat = avg_lmat / np.sum(avg_lmat, axis=0)[np.newaxis,:]
            ht = sns.heatmap(scaled_lmat, yticklabels=covar_labels,
                    xticklabels=(i == len(cvlist)-1),
                    cmap='Blues', ax=axs[i+1], cbar=cbar)
            txt = axs[i+1].set_ylabel(covar, fontsize=14)
        # Save figure
        plt.tight_layout()
        plotname = self.imgdir + "heatmap_" + self.csuff + \
                "_" + attr + "_" + suff + ".png"
        plt.savefig(plotname)
        plt.close()
        print("Plotted graph to " + plotname)

    def calc_svd_corr(self, cvlist):
        print("Calculating covariate correlations with SVD")
        self.ut = self.U.T
        # TODO: ADD covar matrices, generate cvhratio for specific graphs
        # So that we don't have to re-calc each time
        self.cvhratio = []
        self.cvmats = {}
        self.cvtcks = {}
        for covar in cvlist:
            cvcol = self.adata.obs[covar]
            if cvcol.dtype.name == 'category':
                covar_dummy = pd.get_dummies(cvcol)
                self.cvhratio.append(covar_dummy.shape[1])
                self.cvmats[covar] = np.zeros((self.U.shape[1],covar_dummy.shape[1]))
                self.cvtcks[covar] = covar_dummy.columns.tolist()
                for i, cv_lvl in enumerate(covar_dummy.columns):
                    self.cvmats[covar][:,i] = vcorrcoef(self.ut,
                            covar_dummy[cv_lvl].to_numpy().T)
            else:
                self.cvhratio.append(1)
                cvcol = cvcol.to_numpy()
                if np.sum(np.isnan(cvcol)) > 0:
                    ind = np.where((1 - np.isnan(cvcol)) == 1)[0]
                    self.cvmats[covar] = vcorrcoef(self.ut[:,ind], cvcol[ind,np.newaxis].T)[:,np.newaxis]
                else:
                    self.cvmats[covar] = vcorrcoef(self.ut, cvcol[:,np.newaxis].T)[:,np.newaxis]
                self.cvtcks[covar] = False

    def plot_svd_corr(self, cvlist, cbar=False):
        self.calc_svd_corr(cvlist)
        hfull = np.sum(self.cvhratio)
        w = 2.5 * self.U.shape[1] / 20 + 1.5 + 0.8 * cbar
        h = 1 * hfull / 6 + 0.1 * len(cvlist)
        fig, axs = plt.subplots(len(self.cvmats), 1,
            figsize=(w,h),
            gridspec_kw={'height_ratios': self.cvhratio,
                'left': 1.5 / w, 'right': .99 - 0.8 / w * cbar,
                'top':1 - 0.1 / h, 'bottom': .4 / h})
        # Add heatmaps:
        sfact = np.max(self.s) / self.s
        for i, covar in enumerate(cvlist):
            cmat = self.cvmats[covar].T
            cmat = cmat * sfact[np.newaxis,:]
            ht = sns.heatmap(cmat, yticklabels=self.cvtcks[covar],
                    xticklabels=(i == len(cvlist)-1),
                    cmap='RdBu', ax=axs[i], cbar=cbar, center=0)
            txt = axs[i].set_ylabel(covar, fontsize=12,
                    rotation=0, ha='right', va='center')
        # Save figure:
        plt.tight_layout()
        plotname = self.imgdir + "svd_corr_" + self.csuff + ".png"
        plt.savefig(plotname)
        plt.close()
        print("Plotted graph to " + plotname)

    def make_subset_graph(self, suff, covar, cutoff=0.4, plot=True, invert=False):
        if not hasattr(self, 'cvmats') or covar not in self.cvmats.keys():
            self.calc_svd_corr([covar])
        cmat = self.cvmats[covar].T
        sfact = np.max(self.s) / self.s
        cmat = cmat * sfact[np.newaxis,:]
        cmx = np.max(np.abs(cmat),axis=0)
        if invert:
            ind = np.where(cmx < cutoff)[0]
        else:
            ind = np.where(cmx >= cutoff)[0]
        print("Keeping",len(ind),"components")
        # Make corr matrix:
        self.vvt_cv = self.V.T.dot(self.V) / self.U.shape[0]
        self.vvt_sd = np.sqrt(np.diag(self.vvt_cv))
        sub_cv = self.V[ind,:].T.dot(self.V[ind,:]) / self.U.shape[0]
        sub_cv = sub_cv / self.vvt_sd[:,np.newaxis]
        corr_subset = sub_cv / self.vvt_sd[np.newaxis,:]
        suff = suff + '_' + covar
        print("Making new graph under name:", suff)
        self.make_graph(suff, cutoff=.1, corr=corr_subset,
                k=None, scale=None, layout_method='fr')
        if plot:
            print("Plotting leiden for covariate graph under name:", suff)
            self.plot_graph(suff, attr='leiden',
                    show_labels=True, adj_txt=False, w=18)

    def get_modules(self, suff, attr='leiden', print_modules=False):
        # Construct object:
        if suff not in self.modules.keys():
            self.modules[suff] = {}
            la = self.graphs[suff].assign[attr]
            lc = self.graphs[suff].colors[attr]
            nam = np.array(self.graphs[suff].gw.vs['name'])
            las = np.unique(la)
            for ll in las:
                lnam = np.sort(nam[np.where(la == ll)[0]])
                self.modules[suff][ll] = lnam
        # Print:
        if print_modules:
            for ll in self.modules[suff].keys():
                print(ll)
                print(" ".join(self.modules[suff][ll]))
        return(self.modules[suff])

    # TODO: Save enrichments for each graph/config
    def get_goterms(self, suff, attr='leiden'):
        # NOTE: Can't run from private node - internet issue:
        mlist = self.get_modules(suff, print_modules=False)
        self.gpres = {}
        for ll in mlist.keys():
            print(ll)
            testlist = mlist[ll].tolist()
            print(" ".join(testlist))
            self.gpres[ll] = self.gp.profile(organism='hsapiens', query=testlist)
        return(self.gpres)

    # DONE: Fixed graph layout (nogrid works)
    # DONE: Enrichments on graph / other attributes on graph
    # DONE: P-values on corr for cutoff
    # DONE: Quantiles for corr. cutoff
    # DONE: Expr percent dpdt. cutoffs
    # TODO: Saving enrichments so we don't recompute them
    # TODO: Heatmap for the enrichments (make profiles, pass to scanpy dotplots)
    # TODO: Compute gene properties
    # TODO: Plot umap with any coloring (gene properties)


