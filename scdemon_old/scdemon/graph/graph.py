#!/usr/bin/python
"""Gene graph handling."""
# --------------------------------
# Class for handling a gene graph:
# Updated: 06/28/21
# --------------------------------
from .adjacency_matrix import adjacency_matrix
from .utils_graph_plotting import plot_gene_graph, plot_gene_umap
from .utils_enrichment import plot_df_enrichment
from ..auxiliary import vcorrcoef
from ..data import snap_colors

import logging
import numpy as np
import pandas as pd
from scipy import sparse

# For graphs
import umap
import igraph
import leidenalg as la


class gene_graph(object):
    def __init__(self,
                 corr,
                 genes,
                 graph=None,
                 corr_sd=None,
                 cutoff=None,
                 use_zscore=True,
                 edge_weight=None,
                 margin=None,
                 layout_method="fr",
                 knn_k=None,
                 min_size=4,
                 z=None,
                 scale=None,
                 degree_cutoff=0):
        self.corr = corr
        self.genes = genes
        self.graph = graph  # For giving pre-computed graph

        # Adjacency matrix processing object:
        self.adj = adjacency_matrix(
            self.corr,
            corr_sd=corr_sd,
            cutoff=cutoff,
            use_zscore=use_zscore,
            labels=self.genes,
            margin=margin,
            z=z,
            knn_k=knn_k,
            scale=scale,
            degree_cutoff=degree_cutoff,
        )

        self.min_size = min_size
        self.edge_weight = edge_weight
        self.layout_method = layout_method
        self.umat = None

        # Load graph colors resources:
        self.snapcols = snap_colors()
        self.snapcols = list(self.snapcols.to_numpy().T[0])
        self.ncols = len(self.snapcols)

        # For covariate annotations
        self.assign = {}
        self.colors = {}
        # For tracking modules:
        self.modules = {}
        self.mnames = {}
        self.scores = {}
        self.module_match = {}

    # TODO: Rename method (main processing not just build graph)
    def construct_graph(self, resolution=2, layout=True, method='leiden'):
        """Construct the graph object and find gene modules."""
        if self.graph is None:
            self.adjacency, self.kept_genes, self.ind = \
                self.adj.get_adjacency()
            self._make_graph_object()
            self.calculate_gene_modules(resolution=resolution)
        if layout:
            self.layout_graph(layout_method=self.layout_method)

    def construct_full_graph(self, keep_all_z=True):
        """Construct the full adjacency graph for multiplexing purposes."""
        if self.graph is None:
            # Process adjacency matrix:
            _, _, _ = self.adj.get_adjacency(keep_all_z=keep_all_z)
            # Use the full, unaltered set of genes:
            self.adjacency = self.adj.full_adjacency
            self.kept_genes = self.adj.full_labels
            # Create a graph with no pruning at all:
            self._make_graph_object(prune_nodes=False, prune_components=False)

    # TODO: Decide if arguments feed in or keep everything in self.
    def _make_graph_object(self, prune_nodes=True, prune_components=True):
        """Make graph object from adjacency."""
        logging.info("Making graph object from adjacency.")
        # Get edge list from adjacency matrix:
        self.adjacency = self.adjacency.tocoo()
        edge_list = list(tuple(zip(self.adjacency.row, self.adjacency.col)))
        edge_wt_data = self.adjacency.data

        # Make graph from edge list:
        self.graph = igraph.Graph(edge_list, directed=False,
                                  vertex_attrs={'name': self.kept_genes})
        self.graph.simplify(combine_edges="max")  # Merge edges
        # self.graph.vs["name"] = self.kept_genes
        if self.edge_weight is not None:
            self.graph.es["weight"] = self.edge_weight
        else:
            self.graph.es["weight"] = edge_wt_data

        # Clean up very small components (1-2 nodes):
        if prune_components:
            self._prune_graph_components(min_size=self.min_size)

        # Clean up orphan nodes with degree 0:
        if prune_nodes:
            self._prune_graph_nodes(degree=0)

    def _prune_graph_nodes(self, degree=0):
        """Prune igraph graph nodes by degree."""
        logging.info("Removing igraph nodes with degree <= " + str(degree))
        todel = [i for i, x in enumerate(self.graph.degree()) if x <= degree]
        self.graph.delete_vertices(todel)

    def _prune_graph_components(self, min_size=2):
        """Prune igraph graph components by size."""
        logging.info("Removing graph components smaller than " + str(min_size))
        memb = self.graph.clusters().membership  # Connected components
        u, c = np.unique(memb, return_counts=True)  # Sizes of components
        comp = set(u[c < min_size])  # Components to remove
        to_delete = np.array([i for i, x in enumerate(memb) if x in comp])
        self.graph.delete_vertices(to_delete)  # Remove nodes in components

    def layout_graph(self, layout_method='fr'):
        """Compute graph layout."""
        logging.info("Laying out graph, method=" + str(layout_method))
        if layout_method == "fr":
            self.layout = self.graph.layout_fruchterman_reingold(
                niter=500, grid=False)
        else:
            self.layout = self.graph.layout(layout_method)

        logging.info("Kept " + str(len(self.layout)) +
                     " nodes, out of " + str(self.corr.shape[0]))

    def _compute_leiden_partition(self,
                                  resolution=None,
                                  partition_type=None,
                                  use_weights=False,
                                  n_iterations=-1,
                                  random_state=1):
        # Set up the leidenalg arguments:
        partition_kwargs = {}
        partition_kwargs["n_iterations"] = n_iterations
        partition_kwargs["seed"] = random_state

        if partition_type is None:
            if resolution is None:
                partition_type = la.ModularityVertexPartition
            else:
                partition_type = la.RBConfigurationVertexPartition
                partition_kwargs["resolution_parameter"] = resolution

        if use_weights:
            partition_kwargs["weights"] = np.array(
                self.graph.es["weight"]).astype(np.float64)

        # Run leidenalg to get the partition:
        partition = la.find_partition(
            self.graph, partition_type, **partition_kwargs)
        return(partition)

    def get_modules_from_partition(self, partition, method):
        # Initialize and assign each node to a cluster:
        self.assign[method] = (np.zeros(len(self.graph.vs)) - 1).astype(int)
        part_list = [np.array(x) for x in list(partition)]
        nclust = len(part_list)
        for i in range(nclust):
            self.assign[method][part_list[i]] = i

        # Assign colors to each cluster:
        rep_col = int(np.ceil((nclust + 2) / (self.ncols * 1.0)))
        part_cols = np.array((self.snapcols * rep_col)[1:(1 + nclust)])
        self.colors[method] = part_cols[self.assign[method]]
        logging.info("Found " + str(nclust) + " clusters")

    def calculate_gene_modules(self,
                               method="leiden",
                               resolution=None,
                               partition_type=None,
                               use_weights=False,
                               n_iterations=-1,
                               random_state=1):
        """Calculate modules from gene-gene graph using graph clustering."""
        logging.info("Running module detection using method=" + str(method))
        if method == "leiden":
            partition = self._compute_leiden_partition(
                resolution=resolution,
                partition_type=partition_type,
                use_weights=use_weights,
                n_iterations=n_iterations,
                random_state=random_state)
        else:
            # TODO: implement louvain, resolution seems more tunable?
            logging.warning("Louvain, other methods not implemented yet.")
        # Get modules once the partition is calculated:
        self.get_modules_from_partition(partition, method)

    # TODO: Allow compute on adjusted adjacency
    def _compute_umap(self):
        """Calculate UMAP for correlation estimate underlying graph."""
        logging.info("Calculating UMAP for graph adjacency.")
        graph_names = self.graph.vs["name"]
        if self.umat is None or self.umat.shape[1] != len(graph_names):
            # Subset to the appropriate genes for the graph:
            gind = np.array([np.where(self.genes == x)[0][0]
                             for x in graph_names])
            corr_gene_subset = self.corr[gind[:, np.newaxis], gind]
        # Fit the UMAP on the subsetted correlation:
        uw = umap.UMAP()
        self.umat = uw.fit_transform(corr_gene_subset)

    def plot(self, basis='graph', attr=None, color=None,
             imgdir='./', suffix=None, ext='png', **kwargs):
        """Plot the genes under a specific basis (umap/graph)."""
        # Set color from attribute and specify plotname:
        # TODO: clean up how colors are handled
        plotname = imgdir + basis + "_"
        if attr is None:
            col = None
        else:
            if color is not None:
                self.colors[attr] = color
            col = self.colors[attr]
            plotname += str(attr) + "_"
        if suffix is not None:
            plotname += suffix + "." + ext

        # Plot the genes on a specific basis (graph/umap):
        if basis == 'graph':
            self.plot_graph(col=col, plotname=plotname, attr=attr, **kwargs)
        elif basis == 'umap':
            if self.umat is None:
                self._compute_umap()
            self.plot_umap(col=col, plotname=plotname, **kwargs)

    def plot_umap(self, col=None, plotname=None,
                  ax=None, title=None, width=12):
        """Plot the gene UMAP plot with a specified attribute."""
        plot_gene_umap(self.umat, col=col,
                       width=width, title=title,
                       plotname=plotname, ax=ax)

    def plot_graph(self, col=None, width=24, plotname=None, ax=None,
                   title=None, attr='leiden',
                   show_labels=True, frac_labels=1.0, adjust_labels=False):
        """Plot the gene-gene graph with a specified attribute."""
        if not hasattr(self, 'layout'):
            self.layout_graph(layout_method=self.layout_method)
        plot_gene_graph(graph=self.graph, layout=self.layout,
                        assign=self.assign['leiden'],
                        plotname=plotname, ax=ax, title=title,
                        col=col, width=width,
                        show_labels=show_labels,
                        frac_labels=frac_labels,
                        adjust_labels=adjust_labels)

    def get_modules(self, attr="leiden", print_modules=False):
        """Get list of modules from graph and clustering."""
        if attr not in self.modules.keys():
            # Construct object if necessary:
            self.populate_modules(self.adata, attr)
        modules = self.modules[attr]
        if print_modules:
            for ll in modules.keys():
                print(ll, " ".join(modules[ll]))
        return modules

    def get_module_assignment(self, attr="leiden"):
        """Get module assignment for each gene as a pandas DataFrame."""
        if attr not in self.modules.keys():
            # Construct object if necessary:
            self.populate_modules(self.adata, attr)
        mdf = pd.DataFrame({'gene': self.genes,
                            'module': self.module_match[attr]})
        return mdf

    def find_gene(self, gene, attr='leiden', print_genes=False):
        """Find module corresponding to a gene."""
        mkey = None
        for key in self.modules[attr].keys():
            if gene in self.modules[attr][key]:
                mkey = key
        if mkey is None:
            out = None
            if print_genes:
                print("Gene (" + gene + ") is not present in any modules")
        else:
            out = self.modules[attr][mkey]
            if print_genes:
                print(mkey, " ".join(out))
        return out

    # TODO: Speed this up (possibly avg. by tform multiplication)
    # Modules lists functions (requires external adata or X):
    def populate_modules(self, adata, attr='leiden'):
        """Populate modules data."""
        logging.info("Populating modules data")
        # TODO: make orderedDict for the module names?
        partition = self.assign[attr]
        nam = np.array(self.graph.vs["name"])
        modules = np.unique(partition)
        issp = sparse.issparse(adata.X)  # TODO: get this into graph class
        # Initialize:
        self.modules[attr] = {}
        self.mnames[attr] = [""] * (np.max(modules) + 1)
        self.scores[attr] = np.zeros((len(adata), np.max(modules) + 1))
        for i in modules:
            mgenes = np.sort(nam[np.where(partition == i)[0]])
            subadata = adata[:, mgenes]
            avgexpr = np.mean(subadata.X, axis=1)
            if issp:
                avgexpr = np.array(avgexpr).T[0]
            else:
                avgexpr = avgexpr.T
            self.modules[attr][i] = mgenes
            self.scores[attr][:, i] = avgexpr  # For plotting on cell umap

            # Find the top 3 genes (for naming module):
            cv = vcorrcoef(subadata.X.T, avgexpr)
            topgenes = ", ".join(mgenes[np.argsort(-cv)][0:3])
            ngstr = str(len(mgenes))
            ngstr = ngstr + " gene" if ngstr == "1" else ngstr + " genes"
            self.mnames[attr][i] = "M" + str(i) + \
                " (" + ngstr + "): " + topgenes

    # TODO: Speed this up - sort step can be slow
    def match_genes_to_modules(self, attr='leiden'):
        """Match genes to the closest to module by its correlation."""
        logging.info("Matching genes to modules")
        # Nearest module by top correlations across the full dataset:
        if attr not in self.module_match.keys():
            self.module_match[attr] = np.zeros(
                (self.corr.shape[0], len(self.modules[attr])))
            mxval = np.max(np.abs(self.corr)) * 10  # Ensure core stays
            for key in self.modules[attr].keys():
                # Subset correlation to all genes vs. module genes:
                genes = self.modules[attr][key]
                lind = np.in1d(self.genes, genes)
                subcorr = self.corr[:, lind]
                # Score a module by the average corr. of top 5 module genes
                self.module_match[attr][:, key] = np.mean(
                    np.sort(subcorr, 1)[:, -5:], 1)
                # Ensure core genes in module stay in it in the full assign:
                self.module_match[attr][lind, key] = mxval
            # Set each gene's module to the top average correlation:
            self.module_match[attr] = np.argmax(self.module_match[attr], 1)

    # Add save functions here (graph and full object as well):
    def save_modules(self, attr="leiden", as_df=True, filename=None):
        """Save module list as txt or tsv."""
        if as_df:
            mdf = pd.DataFrame({"gene": np.array(self.graph.vs["name"]),
                                "leiden": self.assign[attr],
                                "color": self.colors[attr]})
            mdf.to_csv(filename, sep="\t")
        else:
            mod_list = self.get_modules(attr=attr, print_modules=False)
            with open(filename, "w") as f:
                for key in mod_list.keys():
                    line = str(key) + ": " + " ".join(mod_list[key]) + "\n"
                    f.write(line)
        logging.info("Wrote modules to:" + str(filename))

    def plot_df_enrichment(self, df, col, plotname, attr='leiden', title=None):
        statsdf = plot_df_enrichment(df=df, col=col, genes=self.genes,
                                     module_match=self.module_match[attr],
                                     mnames=self.mnames[attr],
                                     plotname=plotname, title=title)
        return(statsdf)
