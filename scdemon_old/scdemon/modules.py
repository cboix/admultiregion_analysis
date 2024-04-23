#!/usr/bin/python
"""Class for computing modules on adata object."""

# Internal:
from .graph import gene_graph
from .correlation import correlation
from .auxiliary import calculate_svd_covar_corr
from .plotting import plot_umap_grid, plot_svd_corr
from .utils_multigraph import make_graphlist, partition_graphlist

import os
import re
import gc
import logging
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
import scanpy as sc  # TODO: Don't import if not using adata formulation
from scipy import sparse
import igraph
import contextlib

# Plotting:
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, rgb2hex
import seaborn as sns

# For GO annotations (needs connection):
from gprofiler import GProfiler


class scdemon(object):
    """
    Calculate gene modules from anndata or stand-alone matrix.
    """

    def __init__(self,
                 adata,
                 csuff,
                 h5ad_file=None,  # File corresponding to adata (for saving)
                 imgdir=None,
                 estimate_sd=False,
                 seed=1,
                 filter_expr=None,
                 z=4, svd_k=100,
                 calc_raw=False,
                 use_fbpca=True):
        """Initialize gene modules object."""
        # TODO: implement stand-alone matrix
        # Arguments:
        self.adata = adata
        self.h5ad_file = h5ad_file
        self.genes = self.adata.var_names  # labels
        self.csuff = csuff
        self.z = z
        self.svd_k = svd_k
        self.filter_expr = filter_expr
        self.use_fbpca = use_fbpca
        self.estimate_sd = estimate_sd
        self.seed = seed
        self.calc_raw = calc_raw
        # Directories:
        if imgdir is not None:
            self.imgdir = imgdir
        else:
            self.imgdir = "./"
        # Additional metadata / files:
        self.graphs = {}  # Store computed graphs
        self.cv_mats = {}
        self.gp = GProfiler(return_dataframe=True)

    # TODO: Enable this to work with / without an adata object
    def setup(self):
        """Set up standard graph."""
        self._filter_data()
        self._calculate_margin()
        self._calculate_correlation()

    def _filter_data(self):
        # Filter to top expr:
        if self.filter_expr is not None:
            margin = np.mean(self.adata.X > 0, axis=0)
            if sparse.issparse(self.adata.X):
                margin = np.array(margin)[0]
            filtgenes = self.adata.var_names[margin >= self.filter_expr]
            self.adata = self.adata[:, filtgenes]
            self.genes = self.adata.var_names
            print(self.adata.shape)

    def _calculate_margin(self):
        # Calculate margin for filtering later.
        self.margin = np.mean(self.adata.X > 0, axis=0).copy()
        if sparse.issparse(self.adata.X):
            self.margin = np.array(self.margin)[0]

    def _calculate_correlation(self):
        np.random.seed(self.seed)
        if 'X_pca' in self.adata.obsm:
            U = self.adata.obsm["X_pca"]  # TODO: FIX U for scanpy
            s = self.adata.uns["pca"]["variance"]
            V = self.adata.varm["PCs"].T
            self.cobj = correlation(X=self.adata.X,
                                    margin=self.margin,
                                    U=U, s=s, V=V,
                                    k=self.svd_k,
                                    calc_raw=self.calc_raw)

        elif self.use_fbpca:
            # Allow handler to compute PCA with FBPCA library by not feeding in
            # any matrices for U, s, or V
            # TODO: Allow centering if under a certain size / not sparse
            self.cobj = correlation(X=self.adata.X,
                                    margin=self.margin,
                                    k=self.svd_k,
                                    calc_raw=self.calc_raw,
                                    center=False)
        else:
            # Compute PCA with scanpy options if not using FBPCA
            if "X_pca" not in self.adata.obsm.keys():
                logging.debug("Computing PCA through scanpy")
                sc.tl.pca(self.adata, n_comps=self.k)

            U = self.adata.obsm["X_pca"]  # TODO: FIX U for scanpy
            s = self.adata.uns["pca"]["variance"]
            V = self.adata.varm["PCs"].T
            self.cobj = correlation(X=self.adata.X,
                                    margin=self.margin,
                                    U=U, s=s, V=V,
                                    k=self.svd_k,
                                    calc_raw=self.calc_raw)

        # Setup and get the correlation objects
        self.cobj.setup()  # Primarily to calculate PCA if not done yet
        self.corr_est, self.corr_raw = self.cobj.get_correlation()
        # Estimate the std. deviation of the transformed correlation estimate
        if self.estimate_sd:
            self.corr_mean, self.corr_sd = self.cobj.estimate_corr_sd(nperm=50)
        else:
            self.corr_mean, self.corr_sd = None, None

    def _make_graph_object(self,
                           graph_id,
                           corr,
                           graph=None,
                           corr_sd=None,
                           cutoff=None,
                           edge_weight=None,
                           knn_k=None,
                           scale=None,
                           z=None,
                           degree_cutoff=0,
                           min_size=4,
                           use_zscore=True,
                           layout_method="fr"):
        if z is not None:
            self.z = z
        # TODO: Simplify inheritance of kwargs params:
        self.graphs[graph_id] = gene_graph(
            corr,
            genes=self.genes,
            graph=graph,
            corr_sd=corr_sd,
            cutoff=cutoff,
            use_zscore=use_zscore,
            edge_weight=edge_weight,
            margin=self.margin,
            layout_method=layout_method,
            knn_k=knn_k,
            z=self.z,
            scale=scale,
            min_size=min_size,
            degree_cutoff=degree_cutoff)

    # Make a merged graph from multiple decorrelation resolutions:
    # NOTE: Performs better than each graph singly at calling modules
    # TODO: Figure out how to merge...
    # TODO: Assign modules later should be using z-scores
    # TODO: If given new z-score, re-threshold all graphs
    def make_merged_graph(self,
                          graph_id,
                          power_list=[0, .25, .5, .75, 1],
                          resolution=2,
                          keep_all_z=False,
                          **kwargs):
        """Make a merged graph from a list of many parameters."""
        # Make the full list of graphs at each power:
        # TODO: Move make_graphlist into modules?
        graphlist, graphs = make_graphlist(self, plist=power_list,
                                           keep_all_z=keep_all_z)
        # Multiplex cluster the graphs:
        membership = partition_graphlist(graphlist, resolution=resolution)
        # Calculate the average correlation:
        # TODO: Corr is correlation if not loaded as zcutoff.....!
        # NOTE: Extended assignment is always better with zcutoff
        # TODO: Standardize graph: use cutoff to put matrix -> zmat directly
        corr = self.graphs[graphs[0]].corr
        for graph in graphs[1:]:
            # Scale graphs appropriately:
            maxval = np.max(np.abs(self.graphs[graph].corr))
            corr += (self.graphs[graph].corr / maxval)
        corr = corr / (len(graphs) * 1.0)

        # Make the graph object with the pre-computed graph:
        # NOTE: Not good! but otherwise igraph prints a huge message...
        with contextlib.redirect_stdout(None):
            graph = igraph.union(graphlist)
        graph = graph.simplify(combine_edges="max")  # Merge edges
        self._make_graph_object(
            graph_id, corr=corr, graph=graph, **kwargs)
        self.graphs[graph_id].kept_genes = graph.vs['name']

        # Turn multiplex partition into modules (reorder due to merge):
        gn = np.array(graphlist[0].vs['name'])
        reord = np.array([np.where(gn == x)[0][0] for x in graph.vs['name']])
        # gn[reord] == graph.vs['name']  # If we want to check correct.
        membership = np.array(membership)[reord]
        ptns = np.unique(membership)
        partition = [np.where(membership == x)[0] for x in ptns]
        self.graphs[graph_id].get_modules_from_partition(partition, 'leiden')

        # Layout the merged graph:
        self.graphs[graph_id].layout_graph()
        # Populate modules??
        self.graphs[graph_id].populate_modules(self.adata, attr='leiden')
        self.graphs[graph_id].match_genes_to_modules(attr='leiden')

    # Build the graph object and store in dictionary:
    # TODO: Clarify the hierarchy of options - merge use_zscore and z and
    # cutoff if necessary.
    def make_graph(self,
                   graph_id,
                   corr=None,
                   corr_sd=None,
                   power=0,
                   resolution=2,
                   layout=True,
                   adjacency_only=False,
                   full_graph_only=False,
                   keep_all_z=True,
                   **kwargs):
        # For debugging purposes (testing power, if recomputing, etc.)
        self.corr_est, self.corr_raw = self.cobj.get_correlation(power=power)
        # Set parameters + correlations:
        # TODO: Check kwargs passed properly (corr, corr_sd?)
        if corr is None:
            if graph_id == "raw":
                if self.corr_raw is None:
                    self.corr_raw = self.cobj.get_raw_correlation(force=True)
                usecorr = self.corr_raw
                usesd = None
            else:
                usecorr = self.corr_est
                usesd = self.corr_sd
        else:
            usecorr = corr
            usesd = corr_sd
        # TODO: Simplify inheritance of kwargs params:
        # Make the actual graph object:
        self._make_graph_object(
            graph_id, corr=usecorr, corr_sd=usesd, **kwargs)
        # Process graph:
        if adjacency_only:
            # Only compute adjacency, for bootstraps
            _, _, _ = self.graphs[graph_id].adj.get_adjacency()
        elif full_graph_only:
            # Compute the full, unaltered, graph for multiplexing
            self.graphs[graph_id].construct_full_graph(keep_all_z=keep_all_z)
        else:
            # Build the graph, get modules, and annotate genes:
            # TODO: simplify where adata is called (for populate modules):
            self.graphs[graph_id].construct_graph(
                resolution=resolution, layout=layout)
            self.graphs[graph_id].populate_modules(self.adata, attr='leiden')
            self.graphs[graph_id].match_genes_to_modules(attr='leiden')

    def make_svd_subset_graph(self, graph_id, k_ind, power=0, **kwargs):
        """Make a graph from a subset of SVD dims."""
        # Make corr matrix:
        corr_subset = self.cobj.get_correlation_subset(k_ind, power=power)
        self.make_graph(graph_id, corr=corr_subset, **kwargs)

    def make_subset_graph(self,
                          graph_id,
                          covar,
                          cv_cutoff=0.4,
                          invert=False,
                          **kwargs):
        """Make a graph from a covariate-correlated subset of SVD dims."""
        if not hasattr(self, "cv_mats") or covar not in self.cv_mats.keys():
            self.calc_svd_corr([covar])
        # TODO: allow multiple components (e.g. celltype + region)
        cmat = self.cv_mats[covar].T
        sfact = np.max(self.cobj.s) / self.cobj.s
        cmat = cmat * sfact[np.newaxis, :]
        cmx = np.max(np.abs(cmat), axis=0)
        if invert:
            k_ind = np.where(cmx < cv_cutoff)[0]
        else:
            k_ind = np.where(cmx >= cv_cutoff)[0]
        nk = len(k_ind)
        logging.info(f"Subsetting to {nk} components for covariate: {covar}")
        self.make_svd_subset_graph(graph_id, k_ind=k_ind, **kwargs)

    def get_k_stats(self, k_list, power=0, resolution=None, **kwargs):
        """Get statistics on # genes and # modules for each k setting."""
        ngenes = []
        nmodules = []
        for k in k_list:
            ind = np.arange(k)
            graph_id = 'test_k' + str(k)
            # Build the graph - don't layout, but get modules:
            corr_subset = self.cobj.get_correlation_subset(ind, power=power)
            try:
                self._make_graph_object(graph_id, corr=corr_subset, **kwargs)
                self.graphs[graph_id].construct_graph(
                    resolution=resolution, layout=False)
                # Store number of genes and modules:
                nmodules.append(np.max(
                    self.graphs[graph_id].assign['leiden'] + 1))
                ngenes.append(len(self.graphs[graph_id].graph.vs))
                # Delete graph after we have metrics:
                del(self.graphs[graph_id])
                gc.collect()
            except BaseException:
                ngenes.append(0)
                nmodules.append(0)
            logging.info(f"{k}: ngenes={ngenes} and nmodules={nmodules}")
        return(ngenes, nmodules)

    def recluster_graph(self, graph_id, resolution=None):
        """Re-cluster graph with different resolution."""
        # TODO: Allow multiple modules at different resolution
        self.graphs[graph_id].calculate_gene_modules(resolution=resolution)

    def get_modules(self, graph_id, attr="leiden", print_modules=False):
        """Get list of modules from graph and clustering."""
        modules = self.graphs[graph_id].get_modules(
            attr, print_modules=print_modules)
        return modules

    def get_module_assignment(self, graph_id, attr="leiden"):
        """Get module assignment for each gene as a pandas DataFrame."""
        mdf = self.graphs[graph_id].get_module_assignment(attr=attr)
        return mdf

    def find_gene(self, graph_id, gene, return_genes=True, print_genes=True):
        """Find module corresponding to a gene."""
        out = self.graphs[graph_id].find_gene(gene, print_genes=print_genes)
        if return_genes:
            return out

    def save_modules(self, graph_id, attr="leiden", as_df=True,
                     filedir="./", filename=None):
        """Save module list for a specific graph as txt or tsv."""
        if filename is None:
            filename = filedir + "module_df_" + self.csuff + "_" + \
                attr + "_" + graph_id
            filename += ".tsv" if as_df else ".txt"
        self.graphs[graph_id].save_modules(
            attr=attr, as_df=as_df, filename=filename)

    def plot_graph(self, graph_id, attr=None, show_labels=None,
                   frac_labels=1.0, adjust_labels=False, width=20):
        if graph_id not in self.graphs:
            # TODO: Default to quantiles.
            self.make_graph(graph_id)  # Default to this
        # Plot graph object, both ways:
        self.graphs[graph_id].plot(basis='graph',
                                   attr=attr,
                                   show_labels=show_labels,
                                   frac_labels=frac_labels,
                                   adjust_labels=adjust_labels,
                                   width=width,
                                   suffix=self.csuff + "_" + graph_id,
                                   imgdir=self.imgdir)

    # Calculate UMAP representation and plot the genes on it:
    def plot_gene_umap(self, graph_id, attr=None, width=12):
        self.graphs[graph_id].plot(basis='umap', attr=attr, width=width,
                                   suffix=self.csuff + "_" + graph_id,
                                   imgdir=self.imgdir)

    # TODO: Add plot of arbitrary vectors (e.g. margin)
    def plot_gene_logfc(self, graph_id, basis='graph', attr="celltype",
                        fc_cut=2, **kwargs):
        logging.info("Running rank genes on", attr)
        sc.tl.rank_genes_groups(self.adata, attr)
        iscat = self.adata.obs[attr].dtype.name == "category"
        # TODO: Figure out to deal with non-categorical
        # Color map:
        cmap = plt.cm.RdYlBu  # TODO: allow other options
        norm = Normalize(vmin=-fc_cut, vmax=fc_cut)
        # Turn into assignment:
        narr = self.adata.uns["rank_genes_groups"]["names"]
        larr = self.adata.uns["rank_genes_groups"]["logfoldchanges"]
        parr = self.adata.uns["rank_genes_groups"]["pvals_adj"]
        names = list(narr.dtype.names)
        for name in names:
            namestr = attr + "_" + re.sub("[ /]", "_", name)
            logging.info("Plotting for" + str(namestr))
            # TODO: Pass the list of names + a score to the gene graph instead
            n1 = np.array(narr[name])
            xind = np.array([np.where(n1 == x)[0][0]
                             for x in self.graphs[graph_id].graph.vs["name"]])
            n1 = n1[xind]
            l1 = np.array(larr[name])[xind]
            p1 = np.array(parr[name])[xind]
            l1 = l1 * (p1 < 0.05)
            l1[l1 < -fc_cut] = -fc_cut
            l1[l1 > fc_cut] = fc_cut
            # Color mapping:
            vcols = [rgb2hex(c) for c in cmap(norm(-l1))]
            # Plot graph:
            self.graphs[graph_id].plot(basis=basis, attr=namestr, color=vcols,
                                       suffix=self.csuff + "_" + graph_id,
                                       imgdir=self.imgdir, **kwargs)

    def _calculate_covariate_lengths(self):
        if not hasattr(self, 'cv_lengths'):
            self.cv_lengths = {}
            self.cv_ticks = {}
        cvlist = self.adata.obs.columns
        for covar in cvlist:
            if covar not in self.cv_lengths.keys():
                cvcol = self.adata.obs[covar]
                if is_categorical_dtype(cvcol):
                    covar_dummy = pd.get_dummies(cvcol)  # To match cv_mats
                    self.cv_lengths[covar] = covar_dummy.shape[1]
                    self.cv_ticks[covar] = covar_dummy.columns.tolist()
                else:
                    self.cv_lengths[covar] = 1
                    self.cv_ticks[covar] = False

    # TODO: Keep refactoring suff to graph_id below here:
    # TODO: add plotting scale parameters:
    # TODO: Clean up - use pre-computed scores!
    def plot_heatmap_avgexpr(self, graph_id, cvlist=None,
                             attr="leiden", cbar=False):
        """Plot heatmap of module average expression in each covariate."""
        plotname = self.imgdir + "heatmap_" + \
            self.csuff + "_" + graph_id + "_" + attr + ".png"
        if cvlist is None:
            cvlist = ["celltype", "region", "niareagansc",
                      "cogdx", "Apoe_e4", "msex"]
        self._calculate_covariate_lengths()
        hratio = [1] + [self.cv_lengths[covar] for covar in cvlist]
        hfull = np.sum(hratio)
        # For subgraph colors heatmap:
        partition = self.graphs[graph_id].assign[attr]
        part_cols = self.graphs[graph_id].colors[attr]
        u, c = np.unique(partition, return_index=True)
        col_mat = u[np.newaxis, :]
        col_cmap = sns.color_palette(part_cols[c])
        # Make heatmap figure:
        scores = self.graphs[graph_id].scores[attr]
        w = 7 * scores.shape[1] / 20 + 1.2 + 0.8 * cbar
        h = 2 * hfull / 6 + 0.1 * len(cvlist)
        fig, axs = plt.subplots(len(cvlist) + 1, 1, figsize=(w, h),
                                gridspec_kw={"height_ratios": hratio,
                                             "hspace": 0.05,
                                             "left": 1.2 / w,
                                             "right": 1 - 0.8 / w,
                                             "top": 1 - 0.1 / h,
                                             "bottom": 0.4 / h})
        # Plot the module colors:
        sns.heatmap(col_mat, yticklabels=False, xticklabels=False,
                    cmap=col_cmap, ax=axs[0], cbar=False)
        # NOTE: only works for categorical covariates:
        for i, covar in enumerate(cvlist):
            covar_dummy = pd.get_dummies(self.adata.obs[covar])
            covar_dummy = covar_dummy.to_numpy()
            covar_cols = np.sum(covar_dummy, axis=0)
            tform = covar_dummy / covar_cols[np.newaxis, :]
            avg_lmat = tform.T.dot(scores)
            scaled_lmat = avg_lmat / np.sum(avg_lmat, axis=0)[np.newaxis, :]
            sns.heatmap(scaled_lmat,
                        yticklabels=self.cv_ticks[covar],
                        xticklabels=(i == len(cvlist) - 1),
                        cmap="Blues", ax=axs[i + 1], cbar=cbar)
            axs[i + 1].set_ylabel(covar, fontsize=14)
        # Save heatmap figure:
        plt.tight_layout()
        plt.savefig(plotname)
        plt.close()
        print("Plotted graph to " + plotname)

    def calc_svd_corr(self, cvlist):
        self.cv_mats = calculate_svd_covar_corr(
            self.cobj.U.T, self.adata.obs, cvlist, cv_mats=self.cv_mats)

    def plot_svd_corr(self, cvlist, cbar=False):
        self.calc_svd_corr(cvlist)
        plotname = self.imgdir + "svd_corr_" + self.csuff + ".png"
        self._calculate_covariate_lengths()
        cv_hratio = [self.cv_lengths[covar] for covar in cvlist]
        plot_svd_corr(cvlist=cvlist,
                      cv_hratio=cv_hratio,
                      cv_mats=self.cv_mats,
                      cv_ticks=self.cv_ticks,
                      svec=self.cobj.s,
                      plotname=plotname,
                      cbar=cbar)

    # TODO: Plot a single module, several, or all modules:
    def plot_umap_grid(self, graph_id, attr='leiden', plotname=None,
                       ind=None, sel=None, width=2, s=0.5):
        if plotname is None:
            plotname = self.imgdir + "module_umap_grid_" + \
                self.csuff + "_" + graph_id + ".png"
        # Plot umap grid:
        plot_umap_grid(umat=self.adata.obsm["X_umap"],
                       scores=self.graphs[graph_id].scores[attr],
                       titles=self.graphs[graph_id].mnames[attr],
                       plotname=plotname,
                       ind=ind, sel=sel, width=width, s=s)

    # TODO: Save enrichments for each graph/config
    def get_goterms(self, graph_id, attr="leiden", organism="hsapiens",
                    sources=["GO:CC", "GO:BP", "GO:MF",
                             "REAC", "WP", "KEGG", "CORUM"]):
        # NOTE: Can't run from private node - internet issue:
        mlist = self.get_modules(graph_id)
        self.gpres = {}
        for ll in mlist.keys():
            testlist = mlist[ll].tolist()
            logging.info(str(ll) + ": " + " ".join(testlist))
            self.gpres[ll] = self.gp.profile(
                organism=organism, query=testlist, sources=sources)
            logging.info(self.gpres[ll].head())
        return self.gpres

    def plot_df_enrichment(self, df, col, graph_id, suffix,
                           attr='leiden', title=None, ext="png"):
        plotname = self.imgdir + "modules_degenr_heatmap_" + \
            suffix + "." + ext
        statsdf = self.graphs[graph_id].plot_df_enrichment(
            df=df, col=col, attr=attr,
            plotname=plotname, title=title)
        return(statsdf)

    def save_adata(self):
        """Save adata for saving object."""
        # TODO: Currently only saves when file doesn't exist
        # come up with a better choice of when to save adata:
        if not os.path.exists(self.h5ad_file):
            self.adata.write(self.h5ad_file)

    # Functions for properly pickling these objects:
    def __getstate__(self):
        """Get state for object, removing adata, for pickling."""
        # Save the adata object first
        self.save_adata()
        # Copy the object's state from self.__dict__
        state = self.__dict__.copy()
        # Remove the unpicklable entries (adata may be huge):
        del state['adata']
        # TODO: Should we also remove cobj (U, s, V?)
        return state

    def __setstate__(self, state):
        """Load object properly, with adata separately."""
        # Restore instance attributes (i.e., filename and lineno).
        self.__dict__.update(state)
        # Restore the adata with the h5ad_file handle:
        self.adata = sc.read(self.h5ad_file)

    # DONE: Fixed graph layout (nogrid works)
    # DONE: Enrichments on graph / other attributes on graph
    # DONE: P-values on corr for cutoff
    # DONE: Quantiles for corr. cutoff
    # DONE: Expr percent dpdt. cutoffs
    # DONE: Writing module list
    # TODO: Saving enrichments so we don't recompute them
    # TODO: Heatmap for the enrichments (make profiles, pass to scanpy dotplot)
    # TODO: Compute gene properties
    # TODO: Plot umap with any coloring (gene properties)
    # TODO: Handle both X and adata
    # TODO: Plot multiple graphs from the same adata
    # TODO: Better docstrings
    # TODO: Plot multiple graphs on the same layout (for subsetting!)
    # TODO: CELL STATE DISCOVERY + ANNOTATION FROM GENE EXPRESSION MODULES.
