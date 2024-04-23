# ------------------------------------------------
# Processing + plotting of specific cell datasets:
# Created: 04/04/21
# ------------------------------------------------
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

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import rcParams
import seaborn as sns

# Single-cell processing:
import scanpy as sc
import anndata
import meld

# -----------------------------------------------------
# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# Datafiles:
ct = 'Mic/Immune'
st = ''
ctstr = re.sub("_","/", ct)
ststr = re.sub("_","/", st)
datadir = '/home/cboix/data/DEVTRAJ/db/multiRegion/'
pref = 'all_brain_regions_filt_preprocessed_scanpy'
fh5ad = datadir + 'matrices/' + pref + '.majorcelltype.' + ct + '.hdf5'
metafile = datadir + pref + '_norm.final_noMB.cell_labels.tsv.gz'

# Metadata:
mdf = pd.read_csv(metafile, sep='\t')
mdf = mdf.loc[mdf.hcelltype == ctstr,:]
if ct == 'Vasc_Epithelia' and (ststr is not ''):
    mdf = mdf.loc[mdf.cell_type_high_resolution == ststr,:]

# Only the EC specific subtypes:
if ct == 'Exc' and st == 'EC':
    ecs = ['Exc AGBL1 GPC5', 'Exc COL25A1 SEMA3D', 'Exc DLC1 SNTG2',
           'Exc RELN COL5A2', 'Exc RELN GPC5', 'Exc SOX11 NCKAP5',
           'Exc TOX3 INO80D', 'Exc TOX3 POSTN', 'Exc TOX3 TTC6']
    mdf = mdf.loc[mdf.region == 'EC',:]
    mdf = mdf.loc[mdf.cell_type_high_resolution.isin(ecs),:]
elif ct == 'Exc' and st == 'ECall':
    mdf = mdf.loc[mdf.region == 'EC',:]
    mdf = mdf.loc[mdf['major.celltype'] == 'Exc',:]

kbc = mdf.barcode.tolist()

# Read in data from hdf5 file:
hf = h5py.File(fh5ad, 'r')
mat = hf['matrix'][:]
bcs = hf['barcodes'][:]
genes = hf['genes'][:]
hf.close()

bcs2 = [x.decode() for x in bcs]
genes2 = [x.decode() for x in genes]

# Make adata object:
d = {}
d['obs'] = bcs2
d['var'] = genes2
d['X'] = mat
adata = anndata.AnnData(**d)
adata.obs_names = bcs2
adata.var_names = genes2
adata.var.columns = ['varnames']
adata.obs.columns = ['obsnames']
gc.collect()

# Annotate with metadata (projid in particular):
adata = adata[kbc,:]
adata.obs['barcode'] = adata.obs_names
mdf.index = mdf.barcode
mdf2 = pd.merge(adata.obs, mdf, how='left', left_index=True, right_index=True)
adata.obs['projid'] = mdf2['projid'].astype('category')
adata.obs['celltype'] = mdf2['cell_type_high_resolution'].astype('category')
reg = [re.sub("_.*","", x) for x in adata.obs.barcode.tolist()]
adata.obs['region'] = reg
pathlist = ['nft','plaq_d','plaq_n']

# From MELD/RW
if ct == 'Mic_Immune':
    lldf = pd.read_csv(datadir + 'Mic_Immune__metadata.tsv', sep='\t')
    lldf.index = lldf.barcode
    lldf2 = pd.merge(adata.obs, lldf, how='left', left_index=True, right_index=True)
    for path in pathlist:
        adata.obs[path + '_pst'] = lldf2['bin_' + path + '__pseudotime']

# -----------------------------------
# Process this cell type with scanpy:
# -----------------------------------
prefstr = '_test_mrad_' + ct + "_" + st

# Filter cells:
sc.pp.filter_genes(adata, min_cells=3)
adata.obs['cell'] = adata.obs_names
sc.pp.filter_cells(adata, min_genes=100)

# Normalize, log1p:
countsafter = None
sc.pp.normalize_per_cell(adata, counts_per_cell_after=countsafter, copy=False)
sc.pp.log1p(adata)

# Can choose to normalize here:
usecombat = True
if usecombat:
    prefstr = prefstr + "_combat"
    sc.pp.combat(adata, key='projid')

# Here, can filter down to highly variable genes if we want:
usehvg = False
n_top_genes = 5000
if usehvg:
    prefstr = prefstr + "_filthvg"
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=True)
    filtgenes = adata.var_names[filter_result.gene_subset]
    adata = adata[:, filtgenes]

# PCA, Neighbors, UMAP:
sc.tl.pca(adata) # N = 50, default
NN = 13
sc.pp.neighbors(adata, n_neighbors=NN)
sc.tl.umap(adata, maxiter=None, random_state=1)

# Plot + save:
sc.pl.umap(adata, color='celltype',
           frameon=False, save=prefstr + '_umap_ctalone.png')
sc.pl.umap(adata, color='projid',
           frameon=False, save=prefstr + '_umap_ctalone_projid.png')
sc.pl.umap(adata, color='region',
           frameon=False, save=prefstr + '_umap_ctalone_region.png')

# Write UMAP coordinates:
adata.obs['U1'] = adata.obsm['X_umap'][:,0]
adata.obs['U2'] = adata.obsm['X_umap'][:,1]

adf = pd.DataFrame(adata.obs)
adf.to_csv('metadata' + prefstr + '.tsv', sep="\t")

# ------------------------------------------------------
# Incorporate the disease metadata and plot on the umap:
# ------------------------------------------------------
adata = adata[adata.obs.region != 'TH',:] # Remove TH
# mmetafile = datadir + 'regions_metadata_list.tsv.gz'
rmetafile = datadir + 'multiregion_path_by_region_projid.tsv'
rdf = pd.read_csv(rmetafile, sep='\t')
rdf2 = pd.merge(adata.obs, rdf, how='left', left_index=False, right_index=False)
rdf2.index = rdf2.barcode
adata.obs['nft'] = np.log1p(rdf2['nft_act'].astype('float'))
adata.obs['plaq_d'] = np.log1p(rdf2['plaq_d_act'].astype('float'))
adata.obs['plaq_n'] = np.log1p(rdf2['plaq_n_act'].astype('float'))

# Plot + save:
sc.pl.umap(adata, color=['nft','plaq_d','plaq_n'],  #size=4,
           frameon=False, save=prefstr + '_umap_ctalone_path.png')

# Plot + save:
dgenes = ['APOE','XAF1','GLDN', 'CD74','HLA-DRB1','CD81', 'MT-ND3','CX3CR1']

# if ct == 'Mic_Immune':
if ct == 'Vasc_Epithelia':
    dgenes = ['MT2A','MT-ND3','HSPA1A', 'ANGPT2','PLAUR','XAF1', 'SDCBP','WNT2B']
    #           ,'MT-ND3','DPYD','CD74', 'CX3CR1','CD81','SERPINE1', 'LRRK1','TLR2', 'UBC', 'NRP2')

sc.pl.umap(adata, color=dgenes, # size=4,
        frameon=False, save=prefstr + '_umap_ctalone_dgenes.png')

# Remove CAMs and T-cells for the rest of the analysis:
if ct == 'Mic_Immune':
    adata = adata[adata.obs.celltype != 'CAMs',:]
    adata = adata[adata.obs.celltype != 'T cells',:]

# -----------------------------------------------
# Run MELD on the matrix from the scanpy dataset:
# -----------------------------------------------
# path = 'plaq_n'
for path in pathlist:
    print(path)
    # Estimate density of each sample over the graph
    pl = [['CTRL','AD'][x] for x in (1 * (adata.obs[path] > 0)).tolist()]
    path_labels = np.array(pl)
    BETA = 150 # NN=13; BETA=150 for Astrocytes
    meld_op = meld.MELD(beta=BETA)
    path_densities = meld_op.fit_transform(adata.X, path_labels)
    path_densities.index = adata.obs.index
    # Normalize densities to calculate sample likelihoods
    path_likelihoods = meld.utils.normalize_densities(path_densities)
    ad_ll = path_likelihoods['AD']
    # Add and plot the likelihoods:
    pathstr = path + '_ll'
    adata.obs[pathstr] = ad_ll
    # Plot + save:
    sc.pl.umap(adata, color=pathstr, # size=8,
            frameon=False, save=prefstr + '_umap_ctalone_' + pathstr + '.png')

sc.pl.umap(adata, color=[x + "_ll" for x in pathlist], # size=8,
        frameon=False, save=prefstr + '_umap_ctalone_path_ll.png')


# --------------------------------------------
# Alternatively, run pseudotime (random walk):
# --------------------------------------------
# Compute a very basic random-walk run:
sys.path.append('/home/cboix/data/DEVTRAJ/bin/multiRegion')
from random_walk_040521_frz import *

# Arguments for running on the full matrix, non-embedded:
# Alternatively, could feed in X_pca, X_lsi, etc.
path = 'plaq_n'
for path in pathlist:
    imputation_path_arr = adata.obs[path].to_numpy()
    # input_mat = adata.X.toarray()
    input_mat = adata.obsm['X_pca']
    # Hyper-parameters:
    nn = 10; beta = 0.90
    # Compute the transition matrix first:
    sparse_L = similarity_matrix(embeddings=input_mat, nn=nn, vectorize=True)
    M = transition_matrix(L=sparse_L)
    sparse_M = csr_matrix(M, dtype=np.float32)
    del(sparse_L, M)
    gc.collect()
    # Compute random walk with a given beta:
    r = random_walk(sparse_M, imputation_path_arr, beta=beta, max_iter=10000)
    norm_r = (r - np.mean(r)) / np.std(r) # z-score imputed score
    z_path_scores = imputation_path_arr
    # Plot this random walk result on the umap:
    pathstr = path + '_rw'
    adata.obs[pathstr] = r
    sc.pl.umap(adata, color=pathstr, size=4, frameon=False,
            save=prefstr + '_umap_ctalone_' + pathstr + '.png')

sc.pl.umap(adata, color=[x + "_rw" for x in pathlist], # size=8,
        frameon=False, save=prefstr + '_umap_ctalone_path_rw.png')

# FOR PLOTTING WITH PSEUDOTIME ONLY (re-does UMAP)
# # ext = '_rw'
# # ext = '_ll'
# # ext = '_pst'
# # ext = ''
# pathvars = [x + ext for x in pathlist]
# scaleup = 50
# X = adata.obs[pathvars].to_numpy()
# X = X / np.max(X)
# X = X * scaleup
# adata.obsm['X_joint'] = np.hstack((X, adata.obsm['X_pca']))
# sc.pp.neighbors(adata, use_rep='X_joint')
# sc.tl.umap(adata)
# sc.pl.umap(adata, legend_loc='on data', color=dgenes + pathvars, size=50,
#                 frameon=False, save=prefstr + 'umap_ctalone_rerun_with' + ext + '.png')

# ----------------------
# Fit GAM for LL and RW:
# ----------------------
from pygam import LinearGAM, s, f

# ext = '_rw'
# ext = '_pst'
ext = '_ll'
# ext = ''
pathvars = [x + ext for x in pathlist]
subadata = adata[adata.obs.celltype == 'Ast GRM3']
# subadata = adata
X = subadata.obs[pathvars].to_numpy()
if ext == '_rw':
    X = X / np.max(X)

dfit = {}
dx = {}
for title in pathvars:
    dx[title] = {}
    dfit[title] = {}

# dgenes = ['APOE','HIF1A','GLDN', 'CD74','HSP90B1','CD81', 'C3','BIN1']
# dgenes = ['APOE','XAF1','GLDN', 'CD74','HLA-DRB1','CD81', 'MT-ND3','CX3CR1']
dgenes = ['APOE','XAF1','GLDN', 'CD74','HLA-DRB1','CD81', 'MT-ND3','CX3CR1',
          'HIF1A','HSP90B1','C1QC','UBC','INPP5D','C3','BIN1']

dgenes = ['GFAP', 'DPP10', 'CD44', 'GPC5', 'RGMA',
    'DGKB', 'NAV2', 'DGKB', 'SLC39A11', 'CRYAB',
    'HSP90AA1', 'HSPH1', 'FTH1', 'HMGB1', 'SPARCL1', 'CDH20']

for gene in dgenes:
    y = subadata[:,gene].X.toarray().T[0]
    ## plotting
    plt.figure();
    fig, axs = plt.subplots(1,3, figsize=(12,4));
    titles = pathvars
    for i, ax in enumerate(axs):
        print(i)
        gam = LinearGAM(s(0))
        gam.gridsearch(X[:,i][:, np.newaxis], y)
        XX = gam.generate_X_grid(term=0)
        fity = gam.partial_dependence(term=0, X=XX)
        a1 = ax.plot(XX[:, 0], fity)
        a2 = ax.plot(XX[:, 0], gam.partial_dependence(term=0, X=XX, width=.95)[1], c='r', ls='--')
        if i == 0:
            lab =ax.set_ylabel(gene)
        title = ax.set_title(titles[i]);
        dx[titles[i]][gene] = XX[:,0]
        dfit[titles[i]][gene] = fity
    plt.savefig('test_smooth_' + ct + "." + gene + '.allpath' + ext + '.png')
    plt.close()

for title in titles:
    ## plotting
    plt.figure(figsize=(10,5));
    ax = plt.gca()
    for gene in dgenes:
        line = ax.plot(dx[title][gene],dfit[title][gene], label=gene)
    tp = ax.set_title(title);
    plt.legend()
    plt.savefig('test_smooth_' + ct + ".allgenes." + title + '.png')
    plt.close()

# -----------------------------------------------------------------------------
# Separate out and fit a model for AD+ and AD- (per plaq?) - very rough for now
# -----------------------------------------------------------------------------
ad_ind = np.where((adata.obs.plaq_n > 0).to_numpy())[0]
ct_ind = np.where((adata.obs.plaq_n == 0).to_numpy())[0]

ext = '_ll'
# ext = '_rw'
ext = ''
pathvars = [x + ext for x in pathlist]
X = adata.obs[pathvars].to_numpy()
if ext == '_rw':
    X = X / np.max(X)

for gene in dgenes:
    y = adata[:,gene].X.toarray().T[0]
    gam1 = LinearGAM(s(0) + s(1) + s(2))
    gam2 = LinearGAM(s(0) + s(1) + s(2))
    gam1.gridsearch(X[ad_ind,:], y[ad_ind])
    gam2.gridsearch(X[ct_ind,:], y[ct_ind])
    ## plotting
    plt.figure();
    fig, axs = plt.subplots(1,3, figsize=(12,4));
    titles = pathvars
    for i, ax in enumerate(axs):
        XX1 = gam1.generate_X_grid(term=i)
        XX2 = gam2.generate_X_grid(term=i)
        ax.plot(XX1[:, i], gam1.partial_dependence(term=i, X=XX1))
        ax.plot(XX1[:, i], gam1.partial_dependence(term=i, X=XX1, width=.95)[1],
                c='r', ls='--')
        ax.plot(XX2[:, i], gam2.partial_dependence(term=i, X=XX2))
        ax.plot(XX2[:, i], gam2.partial_dependence(term=i, X=XX2, width=.95)[1],
                c='g', ls='--')
        if i == 0:
            ax.set_ylabel(gene)
        ax.set_title(titles[i]);
    plt.savefig('test_smooth_' + ct + "." + gene + '.byadplaq.allpath' + ext + '.png')
    plt.close()

# ----------------------------------
# Fit and plot a 3D Surface as well:
# ----------------------------------
from mpl_toolkits import mplot3d
from pygam import s, te

# Fit the model:
ext = '_rw'
ext = '_ll'
pathvars = [x + ext for x in pathlist]
X = adata.obs[pathvars].to_numpy()
if ext == '_rw':
    X = X / np.max(X)

gene = 'MT2A'
# gene = 'MT-ND3'
# gene = 'HSP90AB1'
gene = 'HSPA1A'
y = adata[:,gene].X.toarray().T[0]
# gam = LinearGAM(te(0, 2))
gam = LinearGAM(te(0, 2))
gam.fit(X, y)

# Plot the 3D contour:
plt.ion()
plt.rcParams['figure.figsize'] = (12, 8)
XX = gam.generate_X_grid(term=0, meshgrid=True)
Z = gam.partial_dependence(term=0, X=XX, meshgrid=True)
ax = plt.axes(projection='3d')
ax.plot_surface(XX[0], XX[1], Z, cmap='viridis')
ax.set_xlabel('NFT')
ax.set_ylabel('Neuritic Plaque')
plt.title(gene + " (" + re.sub("_","",ext) + ")")
plt.savefig('test_smooth_' + ct + "." + gene + '.3d_nft_plaqn' + ext + '.png')
plt.close()




from scipy.interpolate import make_interp_spline, BSpline

# 300 represents number of points to make between T.min and T.max
gene = 'APOE'
t = ad_ll.to_numpy()
ind = np.argsort(t)
t = t[ind]
tnew = np.linspace(t.min(), t.max(), 300)

x = adata[:,gene].X.toarray().T[0][ind]
spl = make_interp_spline(t, x, k=1)  # type: BSpline
x_smooth = spl(tnew)

plt.plot(tnew, x_smooth)
plt.title(gene)
plt.savefig('test_smooth_' + gene + '.' + pathstr + '.png')
plt.close()


geneset = ['APOE','CD74','MT-ND3','CD81','TLR2','HLA-DRB1']
for gene in geneset:
    x = adata[:,gene].X.toarray().T[0][ind]
    spl = make_interp_spline(t, x, k=1)  # type: BSpline
    x_smooth = spl(tnew)
    # Concatenate





# 300 represents number of points to make between T.min and T.max
aind = (adata.obs.region != 'TH').to_numpy() # remove TH
gene = 'CD74'
t = adata.obs[pathstr].to_numpy()[aind]
ind = np.argsort(t)
x = adata[:,gene].X.toarray().T[0][aind]
t = t[ind]
x = x[ind]
tnew = np.linspace(t.min(), t.max(), 300)
spl = make_interp_spline(t, x, k=1)  # type: BSpline
x_smooth = spl(tnew)

plt.plot(tnew, x_smooth)
plt.title(gene)
plt.savefig('test_smooth_' + gene + '.' + pathstr + '.png')
plt.close()


# ----------------------------------------------------
# Extract graph, turn into CSR Graph, run nodevectors:
# ----------------------------------------------------




