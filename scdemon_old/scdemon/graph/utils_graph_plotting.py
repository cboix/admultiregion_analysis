#!/usr/bin/python
"""Functions for plotting gene-gene graphs."""

import logging
import numpy as np

import igraph
from adjustText import adjust_text
from matplotlib import pyplot as plt


# Set a standard visual style:
standard_visual_style = {}
standard_visual_style["vertex_size"] = 10
standard_visual_style["vertex_color"] = "slateblue"
standard_visual_style["vertex_frame_color"] = "slateblue"
standard_visual_style["vertex_label"] = ""
standard_visual_style["edge_color"] = "lightgrey"
standard_visual_style["bbox"] = (2400, 2400)
standard_visual_style["margin"] = 10
standard_visual_style["edge_width"] = 0.5


# TODO: Improve method by connectivity / (cross-connect leiden)
def select_label_subset(assign, labels, frac_labels):
    """Select a subset of labels to plot."""
    frac_labels = 0 if frac_labels < 0 else frac_labels
    # Subset labels:
    if frac_labels < 1.0:
        u, c = np.unique(assign, return_counts=True)
        tot = np.round(np.sum(c) * frac_labels, 0)
        # At least one per leiden + rest divided up:
        cr = c - 1
        cr[cr < 0] = 0
        alloc = c - cr
        # Remaining:
        r1 = np.sum(alloc)
        if tot - r1 > 0:
            r2 = tot - r1
            alloc += np.round(cr * (r2 / np.sum(cr))).astype(int)
        # Allocate:
        keeplabs = []
        for ll in u:
            ind = np.random.choice(np.where(assign == ll)[0], alloc[ll])
            keeplabs = keeplabs + ind.tolist()
        subset_labels = [labels[i] if i in keeplabs else ""
                         for i in range(len(labels))]
        labels = subset_labels
    return labels


def plot_adjusted_labels(labels, layout, ax, adjust_labels=True, fontsize=12):
    """Adjust and plot non-empty labels given labels and their positions."""
    keptind = np.where(np.array(labels) != "")[0].tolist()
    if len(keptind) > 0:
        layout_list = list(layout)
        text = [ax.text(layout_list[i][0],
                        layout_list[i][1],
                        labels[i],
                        ha="center",
                        va="center",
                        fontsize=fontsize)
                for i in keptind]
        # TODO: Check that non-adjusted works.
        if adjust_labels:
            adjust_text(text, lim=25)


# TODO: allow extra suffix for graph (e.g. if params are different)
def plot_gene_graph(graph, layout, assign, plotname=None, col=None, width=24,
                    ax=None, title=None,
                    # Label specifications:
                    show_labels=True, frac_labels=1.0, adjust_labels=False):
    """Plot a gene-gene graph."""
    # Parameters for saving / plotting subplots:
    standalone = (ax is None)
    save_plot = (standalone) and (plotname is not None)
    # Set up the visual style for the igraph-plotted graph
    plot_style = standard_visual_style.copy()
    plot_style["vertex_color"] = "slateblue" if col is None else col
    plot_style["vertex_frame_color"] = "slateblue" if col is None else col
    plot_style["bbox"] = (width * 100, width * 100)

    # Plot graph on matplotlib axes:
    if standalone:
        plt.figure(figsize=(width, width))
        ax = plt.gca()
    igraph.plot(graph, layout=layout, target=ax, **plot_style)
    if title is not None:
        ax.set_title(title, fontdict={'fontsize': 16})

    # Set label specifications and plot labels if required:
    if show_labels:
        labels = select_label_subset(assign, graph.vs["name"], frac_labels)
        plot_adjusted_labels(labels, layout, ax=ax,
                             adjust_labels=adjust_labels, fontsize=12)

    # Clean up and save plot if not subplot, etc.
    if save_plot:
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(plotname)
        plt.close()
        logging.info("Plotted graph to " + plotname)


def plot_gene_umap(umat, col=None, plotname=None,
                   ax=None, title=None, width=12):
    # Parameters for saving / plotting subplots:
    standalone = (ax is None)
    save_plot = (standalone) and (plotname is not None)
    if col is None:
        col = 'slateblue'
    if standalone:
        plt.figure(figsize=(width, width))
        ax = plt.gca()

    # Plot UMAP as scatterplot:
    ax.scatter(umat[:, 0], umat[:, 1], color=col, s=8)
    ax.set_aspect("equal", "datalim")
    if title is not None:
        ax.set_title(title, fontdict={'fontsize': 16})

    # Clean up and save plot if not subplot, etc.
    if save_plot:
        plt.tight_layout()
        plt.savefig(plotname)
        plt.close()
        logging.info("Plotted graph to " + plotname)
