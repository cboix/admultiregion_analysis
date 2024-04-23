#!usr/bin/python
"""Example code for using scdemon library."""
# ---------------------------------------------------------------------
# Example - compute co-expression modules using a scanpy anndata object
# Updated: 06/20/22
# ---------------------------------------------------------------------
import logging
import os
import scdemon as sm
import scanpy as sc
import numpy as np
from scdemon.auxiliary import recipe_full

# Set logging level:
logging.basicConfig(level=logging.INFO)


# Load and prep datasets:
# -----------------------
# Or one of scanpy's datasets:
tag = "pbmc_example"
annfile = tag + "_ann.h5ad"
if os.path.exists(annfile):
    adata = sc.read(annfile)
else:
    adata = sc.datasets.pbmc3k()
    recipe_full(adata, preprocess=True, annotate=True)
    adata.write(annfile)

logging.info("Loaded example dataset")

# Set arguments for plotting outputs:
imgdir = "./"
# Make the scdemon object:
max_k = 100
mod = sm.scdemon(adata,
                  # Where files+plots go / what they are called:
                  csuff=tag, imgdir=imgdir,
                  # Options for graph creation:
                  estimate_sd=False,
                  svd_k=max_k, filter_expr=0.05, z=4.5,
                  # Overall usage/correlation computation options:
                  calc_raw=False)
mod.setup()  # Setup the object


# Make graph using the selected parameters for basic analysis:
# ------------------------------------------------------------
graph_id = "base"
mod.make_graph(graph_id, resolution=2.5)
mod.plot_graph(graph_id, attr="leiden", show_labels=True, width=16)
mod.plot_gene_umap(graph_id, attr="leiden", width=16)

# Modules on cell UMAP:
mod.plot_umap_grid(graph_id)

# Get the modules/print out:
mlist = mod.get_modules(graph_id, print_modules=False)
mod.save_modules(graph_id)

# Get functional enrichments for the modules:
gpres = mod.get_goterms(graph_id)


# We can plot additional plots on the gene modules graph or UMAP:
# ---------------------------------------------------------------
# Plot the logFC for a specific covariate on the graph:
graph_id = "base"
covariate = "leiden"
mod.plot_gene_logfc(graph_id, attr=covariate,
                    show_labels=False, width=16, fc_cut=2)

# Plot correlation of SVD components with covariates:
cvlist = ["n_genes", "leiden"]
mod.plot_svd_corr(cvlist)

# Plot the average expression of leiden modules on covariates:
mod.plot_heatmap_avgexpr(graph_id, cvlist=cvlist, attr="leiden")

# Make a graph from only specific SVD covariate-correlated components
# Such as components correlated with celltypes:
mod.make_subset_graph(graph_id, "leiden", cv_cutoff=2.0)
mod.plot_graph(graph_id + "_cls", attr="leiden", show_labels=True, width=16)
mod.plot_heatmap_avgexpr(graph_id, cvlist=['leiden'], attr="leiden")
mod.plot_umap_grid(graph_id)
