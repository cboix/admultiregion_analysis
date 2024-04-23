#!/usr/bin/python
# -----------------------------------------------------------------
# Plot an overview for a single gene:
# TODO LIST:
# 1. Plot it on a UMAP
# 2. Plot % of cells with it for each cell type
# 3. Show if it has specificity to any neuronal subtype (Exc or Inh)
# 4. Show if it has specificity to any region
# 5. Output statistics for plotting elsewhere (R, etc)
# -----------------------------------------------------------------
import glob
import h5py
from scipy import sparse
import re
from sklearn.utils import sparsefuncs
import gzip
import pickle
import numpy as np
import pandas as pd
import time
import gc
import sys
import fire

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

# No scanpy should be required here:
# import scanpy.api as sc
# import anndata

# # Scanpy plotting settings:
# sc.settings.verbosity = 2  # show logging output
# sc.settings.autosave = True  # save figures, do not show them
# sc.settings.set_figure_params(dpi=300)  # set sufficiently high res for saving

# TODO: Add pdfpages
# TODO: Plot the norm version.

class plot_gene_overview(object):
    def __init__(self, gene, outprefix):
        self.gene = gene # NOTE: Could be list
        self.outprefix = str(outprefix)
        # Data file:
        self.matdir = '/broad/compbio_ce/cboix/multiRegion/matrices/'
        self.prefix = 'all_brain_regions_filt_preprocessed_scanpy'
        self.h5file = self.matdir + self.prefix + \
            "_norm" + '_bygene_fullmatrix.swmr.hdf5'
        # Metadata file:
        self.mdir = '/broad/compbio/cboix/DEVTRAJ/db/multiRegion/metadata/'
        self.mfile = self.mdir + self.prefix + '_norm.final_noMB.cell_labels.tsv.gz'
        self.X = {}
        self.bcfile = self.matdir + self.prefix + '_norm.barcodes.tsv.gz'

    # Main workflow:
    def main(self):
        self.get_metadata()
        self.get_genedata(self.gene)
        # Plotting:
        # TODO: Can we put all plots onto a report?
        self.plot_umap(self.gene)
        # self.plot_fractions(self.gene)
        # self.plot_regional_specificity(self.gene)
        # self.plot_neuronal_specificity(self.gene)
        # self.write_stats(self.gene)

    def get_metadata(self):
        print("[STATUS] Reading metadata")
        self.metadf = pd.read_csv(self.mfile, sep="\t")
        self.mbcs = self.metadf.barcode.to_numpy()
        self.metadf.index = self.metadf['barcode']
        self.metadf.index.names = ['index']
        print("[STATUS] Read in metadata, with shape:", self.metadf.shape)

    def get_genedata(self, gene):
        print("[STATUS] Getting gene expression for gene:", gene)
        # Get the gene expression for this gene:
        hf = h5py.File(self.h5file, 'r')
        # Read in genes:
        if not hasattr(self, 'dgenes'):
            print("Genes")
            encgenes = hf['genes']
            self.dgenes = [n.decode() for n in encgenes]
        # Read in barcodes:
        if not hasattr(self, 'dbcs'):
            print("Barcodes")
            encbcs = hf['barcodes']
            self.dbcs = [n.decode() for n in encbcs] # Takes a while.
            self.dbcs = np.array(self.dbcs)
            if 0 == 1:
                # self.dbcs = pd.read_csv(self.bcfile, header=None, sep="\t")[0].to_numpy()
                print('DB', self.dbcs[0:5], '...', self.dbcs[-5:])
                print('MB', self.mbcs[0:5], '...', self.mbcs[-5:])
                dord = np.argsort(self.dbcs)
                print(dord[0:5])
                # self.bind = dord[np.searchsorted(self.dbcs[dord], self.mbcs)]
                self.bind = np.nonzero(np.in1d(self.dbcs, self.mbcs)) # Subset bcs
                self.sbcs = self.dbcs[self.bind]
                print('SB', self.sbcs[0:5], '...', self.sbcs[-5:])
                print(np.sum(self.sbcs == self.mbcs))
                # self.bind = np.searchsorted(self.dbcs, self.mbcs) # Subset bcs
                print('BIND', self.bind[0:5], '...', self.bind[-5:])
                # Rearrange metadata:
                rmeta = self.metadf.loc[self.sbcs]
                print(self.metadf.shape, rmeta.shape)
                self.metadf = rmeta
        # Add gene expr to dictionary, in proper order
        gid = np.searchsorted(self.dgenes, gene)
        print(gid, gene)
        hd = hf['matrix']
        # self.X[gene] = hd[:,gid][self.bind]
        self.X[gene] = hd[:,gid]
        hd.refresh()
        hf.close()
        with open(self.outprefix + "_x.tsv", 'w') as f:
            preformatted_write(self.X[gene], f)
        np.savetxt(self.outprefix + "_db.tsv", self.dbcs)
        print("[STATUS] Extracted gene data for gene:", gene)

    def plot_umap(self, gene, title=None):
        print("[STATUS] Plotting UMAP for gene:", gene)
        print("Quantiles are:", np.quantile(self.X[gene], [0,.25,.5,.75,1]))
        fname = self.outprefix + "_" + gene + "_img_umap.png"
        sns.set(font_scale=1.1)
        ssize=1
        # Plot current trace:
        fig = plt.figure(figsize=(8, 8))
        ax = plt.gca()
        # Plot % change:
        if title is not None:
            ax.set_title(title)
        ax.set_facecolor('white')
        plt.scatter(self.metadf.U1, self.metadf.U2,
                    s=ssize, c=self.X[gene], marker='.',edgecolors='none',
                    cmap=plt.get_cmap('viridis'))
        plt.ylabel('UMAP 1')
        plt.xlabel('UMAP 2')
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(fname, dpi=350, bbox_inches='tight')


# Faster than for loop or np.savetxt:
def preformatted_write(mat, f, fmtstring=None, encode=False):
    if fmtstring is None:
        if len(mat.shape) == 1:
            fmtstring = '%g'
        else:
            fmtstring = '\t'.join(['%g']*mat.shape[1])
    fmt = '\n'.join([fmtstring]*mat.shape[0])
    data = fmt % tuple(mat.ravel())
    if encode:
        data = data.encode('utf-8')
    f.write(data)

if __name__ == "__main__":
    fire.Fire(plot_gene_overview)
