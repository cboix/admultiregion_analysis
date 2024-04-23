# Look at multi-region data - ncells etc:
import glob
import h5py
from scipy import sparse
import re
from sklearn.utils import sparsefuncs
from sklearn.decomposition import IncrementalPCA
import gzip
import pickle
from inc_pca import IncPCA
import numpy as np
import pandas as pd
import time
import gc
import fbpca

# For KNN + conn
import umap
from umap.umap_ import nearest_neighbors
from sklearn.utils import check_random_state
from types import MappingProxyType
from umap.umap_ import fuzzy_simplicial_set

# For clustering:
import igraph as ig
import leidenalg
import hdbscan

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

# ----------
# Functions:
# ----------
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

def gzipped_write(matrix, filename):
    with gzip.open(filename, 'wb') as f:
        preformatted_write(matrix, f, encode=True)

def plot_dimred(mat, filename, redtype=None, title=None,
                values=None, cmap=None, size=.1):
    sns.set(font_scale=1.1)
    # Plot current trace:
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    # Plot % change:
    if title is not None:
        ax.set_title(title)
    ax.set_facecolor('white')
    if values is None:
        plt.scatter(mat[:,0], mat[:,1], s=size)
    else:
        plt.scatter(mat[:,0], mat[:,1], s=size,
                    c=values, cmap=cmap)
    if redtype is not None:
        plt.ylabel(redtype + ' 1')
        plt.xlabel(redtype + ' 2')
    plt.tight_layout()
    fig = plt.gcf()
    fig.savefig(filename, dpi=350, bbox_inches='tight')

def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return(matrix)

# For labels write:
def aggregate_labels_over_matrix(matrix, labels, hdb=False):
    NLAB = np.max(labels) + 1 + 1 * hdb
    avgmat = np.zeros((NLAB, matrix.shape[1]))
    for i in range(NLAB):
        print(i)
        ind = np.where(labels == (i - 1 * hdb))[0]
        avgmat[i,:] = np.array(np.mean(matrix[ind,:], axis=0)[0])[0]
    return(avgmat)

def calc_write_avgmatrix(lblset, hdb=False):
    groupfile = prefix + "." + lblset + '.tsv.gz'
    avgfile = prefix + "." + lblset + '.avg.tsv.gz'
    labels = pd.read_csv(groupfile, header=None, sep="\t").to_numpy().T[0]
    avgmat = aggregate_labels_over_matrix(fullmat, labels, hdb=hdb)
    gzipped_write(avgmat, avgfile)

# ------------------
# Load the datasets:
# ------------------
keptbc = {}
with open('../Annotation/multiRegion_rows.txt','r') as f:
    for line in f:
        line = line.rstrip().split(".")
        if line[0] in keptbc.keys():
            keptbc[line[0]].append(line[1])
        else:
            keptbc[line[0]] = [line[1]]

suffix = "filtered_feature_bc_matrix.h5"
fnames = np.sort(glob.glob("*" + suffix))

# Load in the protein coding genes to subset the data:
anno = pd.read_csv('../Annotation/Gene.v28lift37.annotation.bed',
                   header=None, sep="\t",
                   names=['chr','start','end','strand','ENSG','type','symbol'])
keep_type = ['protein_coding']
pc_symbols = list(anno.symbol[[tp in keep_type for tp in anno['type']]])
pc_ensg = list(anno.ENSG[[tp in keep_type for tp in anno['type']]])
#pc_genes = list(anno.symbol)

# load_data = False
load_data = True

ncells = 0
acell = []
agene = []
h5file = fnames[0]
fullmat = None
sind = None
allbc = []

if load_data:
	for h5file in fnames:
		region = re.sub("_" + suffix,"", h5file)
		print(region, h5file)
		hf = h5py.File(h5file, 'r')
		hd = hf['matrix']
		shape = hd['shape'][:]
		acell.append(shape[1])
		agene.append(shape[0])
		ncells = ncells + shape[1]
			# Get matrix:
		print("Loading in information")
		bc = np.array([b.decode() for b in hd['barcodes'][:]])
		keepind = []
		for kb in keptbc[region]:
			keepind = keepind + [i for i,v in enumerate(bc)
								if re.search('-' + kb + '$',v)]
		keepind = np.sort(np.array(keepind))
		# Make the data matrix:
		data = hd['data'][:]
		indices = hd['indices'][:]
		indptr = hd['indptr'][:]
		shape = hd['shape'][:]
		print("Making matrix")
		matrix = sparse.csc_matrix((data, indices, indptr), shape=shape, dtype=float)
		del(data, indices, indptr)
		# Subset the matrix to kept indices
		matrix = matrix[:, keepind]
		bc = bc[keepind]
		# Filtering genes:
		hfeat = hd['features']
		genes = np.array([g.decode() for g in hfeat['id'][:]])
		symbols = np.array([s.decode() for s in hfeat['name'][:]])
		# Could remove:
		ribo_genes = [name for name in list(symbols) if re.match('^RP[0-9]+-', str(name)) or re.match('^RP[SL]', name)]
		mito_genes = [name for name in list(symbols) if re.match('^MT-', name) or re.match('^MTRNR', name)]
		# Calculate mito and ribo pct:
		rind = [np.where(symbols == n)[0][0] for n in ribo_genes]
		mind = [np.where(symbols == n)[0][0] for n in mito_genes]
		gmarg = np.array(np.sum(matrix, axis=0)[0])[0]
		rmarg = np.array(np.sum(matrix[rind,:], axis=0)[0])[0]
		mmarg = np.array(np.sum(matrix[mind,:], axis=0)[0])[0]
		# Percentage:
		rfrac = rmarg / gmarg
		mfrac = mmarg / gmarg
		# Remove cells with greater than 5% ribo and 20% mito:
		mrkeep = np.where((mfrac < .2) * (rfrac < .05))[0]
		matrix = matrix[:, mrkeep]
		# Final barcodes:
		bc = bc[mrkeep]
		allbc = allbc + [region + "_" + b for b in bc]
		# Filter to only the protein coding genes:
		if sind is None:
			sind = np.array([i for i,s in enumerate(symbols) if s in pc_symbols])
		symbols = symbols[sind]
		genes = genes[sind]
		matrix = matrix[sind,:]
		# Transpose so samples are rows:
		matrix = matrix.T
		print(matrix.shape)
		print("Adding matrix")
		if fullmat is None:
			fullmat = matrix
		else:
			fullmat = sparse.vstack([fullmat, matrix])
		del(matrix)
		print(fullmat.shape)
		hf.close()

gc.collect()
gc.collect()


if load_data:
    # Write out barcodes and symbols + genes
    with open('all_brain_regions_filt_preprocessed_barcodes.txt', 'w') as f:
        for line in allbc:
            f.write(line + '\n')

    with open('all_brain_regions_filt_preprocessed_genes.txt', 'w') as f:
        for i in range(len(symbols)):
            f.write(symbols[i] + '\t' + genes[i] + '\n')


# --------------------------
# Pre-processing before PCA:
# --------------------------
prefix = 'all_brain_regions_filt_preprocessed'
norm_mat = True
if norm_mat:
    prefix = prefix + "_norm"  # update prefix
    if load_data:
        # Normalize per cell to the median # counts:
        gmarg = np.array(np.sum(fullmat, axis=1).T[0])[0]
        cafter = np.median(gmarg)
        gmarg = gmarg / cafter
        sparsefuncs.inplace_row_scale(fullmat, 1/gmarg)
        newmarg = np.array(np.sum(fullmat, axis=1).T[0])[0]
        print(newmarg[0:5])
        # As seurat, log1p transform the data (inplace)
        np.log1p(fullmat.data, out=fullmat.data)

print(prefix)


# ----------------------
# Run PCA on the matrix:
# ----------------------
n_iter = 8
if n_iter != 4:
    prefix = prefix + "_" + str(n_iter)

ufile = prefix + '.u.tsv.gz'
dfile = prefix + '.d.tsv.gz'
vfile = prefix + '.v.tsv.gz'

k = 50
raw = False
t1 = time.time()
print("[STATUS] Performing pca with fbpca.pca and k=" + str(k))
(U, s, Va) = fbpca.pca(fullmat, k=k, raw=raw, n_iter=n_iter, l=None)
print('[STATUS] Finished in ' + str(round(time.time() - t1, 2)) + "s")

print("[STATUS] Writing out results")
gzipped_write(U, ufile)
gzipped_write(s, dfile)
gzipped_write(Va, vfile)

# Read back and combine:
U = pd.read_csv(ufile, header=None, sep="\t").to_numpy()
s = pd.read_csv(dfile, header=None).to_numpy()
sdiag = np.diag(s.T[0])
US = U.dot(sdiag)

# ------------------------------------
# Run the UMAP decomp (approx 2.5 hrs)
# ------------------------------------
NNEIGHBORS = 50

t1 = time.time()
print("[STATUS] Performing UMAP")
fit = umap.UMAP(n_neighbors=NNEIGHBORS)
umap_U = fit.fit_transform(US)
print('[STATUS] Finished in ' + str(round(time.time() - t1, 2))+"s")

print("[STATUS] Writing out UMAP results")
umapfile = prefix + '.umap.tsv.gz'
gzipped_write(umap_U, umapfile)

# Read back in and plot:
umap_U = pd.read_csv(umapfile, header=None, sep="\t").to_numpy()

umapimg = prefix + '.umap.png'
plot_dimred(umap_U, filename=umapimg, redtype='UMAP')

# -----------------------------------
# LEARN THE KNN on the decomposition:
# -----------------------------------
random_state = check_random_state(1)
# First compute neighbors (slowest)
NNEIGHBORS = 50

knn_indices, knn_dists, forest = nearest_neighbors(
    X=US, n_neighbors=NNEIGHBORS, random_state=random_state,
    metric='euclidean', metric_kwds= MappingProxyType({}),
    angular=False, verbose=True)

# Save KNN indices + dists:
kifile = prefix + '.knn_indices.tsv.gz'
kdfile = prefix + '.knn_dists.tsv.gz'

print("[STATUS] Writing out KNN results")
gzipped_write(knn_indices, kifile)
gzipped_write(knn_dists, kdfile)

# Read back in:
knn_indices = pd.read_csv(kifile, header=None, sep="\t").to_numpy()
knn_dists = pd.read_csv(kdfile, header=None, sep="\t").to_numpy()

# ----------------------------------
# Then compute connectivities (fast)
# ----------------------------------
NNEIGHBORS = 50
X = sparse.coo_matrix(([], ([], [])), shape=(U.shape[0], 1))
connectivities = fuzzy_simplicial_set(
    X, NNEIGHBORS,
    None, None,
    knn_indices=knn_indices,
    knn_dists=knn_dists,
    set_op_mix_ratio=1.0,
    local_connectivity=1.0)

if isinstance(connectivities, tuple):
    connectivities = connectivities[0]

conn = connectivities.tocsr()

# ---------------------------------------
# Use connectivities object to do leiden:
# ---------------------------------------
# NOTE: Takes a while, but need more than 2 iterations.
# Running 1 8 and 2 6 and 3 10
resolution = 5
directed = True
use_weights = True
random_state = 1
# 20 took about 3.5 hrs with resolution = 5
n_iterations = 40

groupfile = prefix + '.groups_r' +  str(resolution) + '.tsv.gz'
init = pd.read_csv(groupfile, header=None, sep="\t").to_numpy().T[0]

partition_kwargs = dict()

def get_igraph_from_adjacency(adjacency, directed=directed):
	sources, targets = adjacency.nonzero()
	weights = adjacency[sources, targets]
	if isinstance(weights, np.matrix):
		weights = weights.A1
	g = ig.Graph(directed=directed)
	g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
	g.add_edges(list(zip(sources, targets)))
	try:
		g.es['weight'] = weights
	except:
		pass
	if g.vcount() != adjacency.shape[0]:
		print("error")
	return(g)

g = get_igraph_from_adjacency(conn, directed=directed)
partition_type = leidenalg.RBConfigurationVertexPartition

if use_weights:
	partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)

partition_kwargs['n_iterations'] = n_iterations
partition_kwargs['seed'] = random_state
if resolution is not None:
	partition_kwargs['resolution_parameter'] = resolution

# Perform clustering:
t1 = time.time()
part = leidenalg.find_partition(g, partition_type,
                               initial_membership=init + 1,
                                **partition_kwargs)
print(time.time() - t1)

# Get groups:
groups = np.array(part.membership)

# Save with resolution parameter:
groupfile = prefix + '.groups_r' + str(resolution) + '.tsv.gz'
with gzip.open(groupfile, 'wb') as f:
    preformatted_write(groups, f, encode=True)

# Aggregate:
lblset = 'groups_r' + str(resolution)
calc_write_avgmatrix(lblset, hdb=False)

# resolution = 1
# groupfile = prefix + '.groups_r' + str(resolution) + '.tsv.gz'
# groups = pd.read_csv(groupfile, header=None, sep="\t").to_numpy().T[0]

# NLAB = np.max(groups) + 1
# avgmat = np.zeros((NLAB, fullmat.shape[1]))
# for i in range(NLAB):
#     print(i)
#     ind = np.where(groups == i)[0]
#     avgmat[i,:] = np.array(np.mean(fullmat[ind,:], axis=0)[0])[0]

# # Write out avg:
# avgfile = prefix + '.groups_r' + str(resolution) + '.avg.tsv.gz'
# with gzip.open(avgfile, 'wb') as f:
#     preformatted_write(avgmat, f, encode=True)

# ---------------------------------------
# Alternatively, hdbscan groups (faster):
# ---------------------------------------
umapfile = prefix + '.umap.tsv.gz'
umap_U = pd.read_csv(umapfile, header=None, sep="\t").to_numpy()

# 1.
labels = hdbscan.HDBSCAN(min_samples=10,
                         min_cluster_size=500).fit_predict(umap_U)
np.sum(labels < 0)
groupfile = prefix + '.groups_hdb.tsv.gz'
gzipped_write(labels, groupfile)

# 2.
labels_lg = hdbscan.HDBSCAN(min_samples=2,
                            min_cluster_size=500).fit_predict(umap_U)
np.sum(labels_lg < 0)
groupfile = prefix + '.groups_hdb_lg.tsv.gz'
gzipped_write(labels_lg, groupfile)

# 3. Small
labels_sm = hdbscan.HDBSCAN(min_samples=1,
                            min_cluster_size=250).fit_predict(umap_U)
np.sum(labels_sm < 0)
groupfile = prefix + '.groups_hdb_sm.tsv.gz'
gzipped_write(labels_sm, groupfile)

# -------------------------------------------
# Aggregate the matrix according to clusters:
# -------------------------------------------
# Write out avg:
# lblset = 'groups_hdb'
lblset = 'groups_hdb_lg'
lblset = 'groups_hdb_lg_m1'
lblset = 'groups_hdb_sm'
calc_write_avgmatrix(lblset, hdb=True)

# ---------------------------------------------
# Finally, compute umap from the connectivites:
# ---------------------------------------------
from umap.umap_ import find_ab_params, simplicial_set_embedding
min_dist = 0.5
spread = 1.0
# n_components= 2
n_epochs = 2
alpha = 1.0
gamma = 1.0
negative_sample_rate = 5

a, b = find_ab_params(spread, min_dist)
random_state = check_random_state(1)

conn = connectivities.tocoo()
del(g, knn_indices, knn_dists)

# Umap from connectivities:
umap_C = simplicial_set_embedding(
    US, conn,
    2, 1.0, a, b, 1.0,
    negative_sample_rate,
    n_epochs,
    init='spectral',
    random_state=random_state,
    metric='euclidean',
    metric_kwds= MappingProxyType({}),
    verbose=True)


# Plot conn with U*S and nn=50, should be better?:
umapfile = prefix + '.umap_C.tsv.gz'
with gzip.open(umapfile, 'wb') as f:
    preformatted_write(umap_C, f, encode=True)

# umapfile = prefix + '.umap_C.tsv.gz'
# umap_C = pd.read_csv(umapfile, header=None, sep="\t").to_numpy()
# umapimg = prefix + '.umap_C.png'
# plot_dimred(umap_C, filename=umapimg, redtype='UMAP')

# -------------------------------
# Plot genes or other attributes:
# -------------------------------
marker_dict={'Ast': ['AQP4', 'GFAP', 'SLC1A2'],
             'Exc': ['CAMK2A', 'NRGN', 'SLC17A7'],
             'In': ['GAD1', 'GAD2', 'PVALB', 'SST', 'VIP'],
             'Mic': ['C3', 'CD74', 'CSF1R'],
             'Neu': ['CELF4', 'GRIN1', 'SNAP2', 'SYT1', 'THY1', 'VSNL1'],
             'Oli': ['MBP', 'MOBP', 'OLIG1', 'PLP1'],
             'OPC': ['OLIG2', 'PDGFRA', 'VCAN','BCAN'],
             'End': ['FLT1','CLDN5']}

inh_markers = ['ADRA2A', 'CALB1', 'CALB2', 'CCK', 'CHRNA7', 'CNR1', 'CRHBP', 'CXCL14', 'LAMP5', 'LHX6', 'MAFB',
               'NDNF', 'NOS1', 'NPY', 'NR2F2', 'PDE1A', 'PDYN', 'PNOC', 'PVALB', 'RELN', 'SATB1', 'SOX6', 'SP8',
               'SST', 'SULF1', 'SV2C', 'SYNPR', 'TAC1', 'TAC3', 'TH', 'TNFAIP8L', 'TOX', 'VIP']

exc_markers = ['ADRA2A', 'CARTPT', 'CUX2', 'ETV1', 'FOXP2', 'GLRA3', 'GRIK4', 'LAMP5', 'NR4A2',
               'NTNG2', 'OPRK1', 'PRSS12', 'RORB', 'RPRM', 'RXFP1', 'TLE4', 'TOX', 'TPBG']

# NOTE: Slicing is slow because CSR:
# Plot genes:
genelist = [i for k,v in marker_dict.items() for i in v] + inh_markers + exc_markers

for gene in genelist:
    print(gene)
    gind = np.where(symbols == gene)[0]
    if len(gind) > 0:
        vals = fullmat[:,gind]
        vals = vals.toarray()
        geneimg = prefix + ".geneplot_" + gene + ".png"
        plot_dimred(umap_U, filename=geneimg, redtype='UMAP',
                    values=vals, cmap='viridis')

# Plot the batches:
batch = [re.sub("_.*", "",b) for b in allbc]
regions = pd.unique(batch)
nbatch = [np.where(regions == i)[0][0] for i in batch]
geneimg = prefix + ".batchplot.png"
plot_dimred(umap_U, filename=geneimg, redtype='UMAP',
            values=nbatch, cmap='Accent', size=0.02)


# TODO: Get average gene profiles per groups

# TODO: Run wilcox of gene profiles -- one v all? or v diff?





# allfile = 'all_brain_regions_filtered.hdf5'
# hf = h5py.File(allfile, 'w')
# hf.create_dataset(mtype, data=fullmat,
#                     compression='gzip', compression_opts=9)
# hf.close()
# ds.close()
print(agene)


# ---------------------------
# Preprocess the full matrix:
# ---------------------------
ncell = np.array(np.sum(matrix > 0, axis=0)[0])[0]
cind = np.where(ngenes >= 50)[0] # At least 50 cells
genes = genes[cind]
symbols = symbols[cind]
matrix = matrix[:,cind]

# allbc = allbc[gind]

ncells = np.array(np.sum(matrix, axis=1).T[0])[0]
cind = np.where(ncells >= 10)[0] # At least 10 cells
matrix = matrix[cind,:]


# Could remove:
ribo_genes = [name for name in list(symbols) if re.match('^RP[0-9]+-', str(name)) or re.match('^RP[SL]', name)]
mito_genes = [name for name in list(symbols) if re.match('^MT-', name) or re.match('^MTRNR', name)]

# Calculate mito and ribo pct:
rind = [np.where(symbols == n)[0][0] for n in ribo_genes]
mind = [np.where(symbols == n)[0][0] for n in mito_genes]
gmarg = np.array(np.sum(matrix, axis=0)[0])[0]
rmarg = np.array(np.sum(matrix[rind,:], axis=0)[0])[0]
mmarg = np.array(np.sum(matrix[mind,:], axis=0)[0])[0]
# Percentage:
rfrac = rmarg / gmarg
mfrac = mmarg / gmarg

# Remove cells with greater than 5% ribo and 20% mito:
mrkeep = np.where((mfrac < .2) * (rfrac < .05))[0]
matrix = matrix[:, mrkeep]
barcodes = barcodes[mrkeep]

# --------------------------
# Pre-processing before PCA:
# --------------------------
# Normalize per cell to the median # counts:
gmarg = np.array(np.sum(matrix, axis=0)[0])[0]
cafter = np.median(gmarg)
gmarg = gmarg / cafter
sparsefuncs.inplace_column_scale(matrix, 1/gmarg)
newmarg = np.array(np.sum(matrix, axis=0)[0])[0]
print(newmarg)

# As seurat, log1p transform the data (inplace)
np.log1p(matrix.data, out=matrix.data)







# See if we can compute covariance in one sample (then we can run 7x7 jobs):
h5file = fnames[1]
hf = h5py.File(h5file, 'r')
hd = hf['matrix']
bc = hd['barcodes']

# Get matrix:
barcodes = bc[:]
data = hd['data'][:]
indices = hd['indices'][:]
indptr = hd['indptr'][:]
shape = hd['shape'][:]
matrix = sparse.csc_matrix((data, indices, indptr), shape=shape, dtype=float)

# Get genes:
hfeat = hd['features']
genes = np.array([g.decode() for g in hfeat['id'][:]])
symbols = np.array([s.decode() for s in hfeat['name'][:]])

hf.close()

# ---------------------------
# Filter out genes and cells:
# ---------------------------
ncells = np.array(np.sum(matrix, axis=1).T[0])[0]
cind = np.where(ncells >= 10)[0] # At least 10 cells
matrix = matrix[cind,:]
genes = genes[cind]
symbols = symbols[cind]

ngenes = np.array(np.sum(matrix > 0, axis=0)[0])[0]
gind = np.where(ngenes >= 500)[0] # At least 10 cells
matrix = matrix[:,gind]
barcodes = barcodes[gind]

# Could remove:
ribo_genes = [name for name in list(symbols) if re.match('^RP[0-9]+-', str(name)) or re.match('^RP[SL]', name)]
mito_genes = [name for name in list(symbols) if re.match('^MT-', name) or re.match('^MTRNR', name)]

# Calculate mito and ribo pct:
rind = [np.where(symbols == n)[0][0] for n in ribo_genes]
mind = [np.where(symbols == n)[0][0] for n in mito_genes]
gmarg = np.array(np.sum(matrix, axis=0)[0])[0]
rmarg = np.array(np.sum(matrix[rind,:], axis=0)[0])[0]
mmarg = np.array(np.sum(matrix[mind,:], axis=0)[0])[0]
# Percentage:
rfrac = rmarg / gmarg
mfrac = mmarg / gmarg

# Remove cells with greater than 5% ribo and 20% mito:
mrkeep = np.where((mfrac < .2) * (rfrac < .05))[0]
matrix = matrix[:, mrkeep]
barcodes = barcodes[mrkeep]

# --------------------------
# Pre-processing before PCA:
# --------------------------
# Normalize per cell to the median # counts:
gmarg = np.array(np.sum(matrix, axis=0)[0])[0]
cafter = np.median(gmarg)
gmarg = gmarg / cafter
sparsefuncs.inplace_column_scale(matrix, 1/gmarg)
newmarg = np.array(np.sum(matrix, axis=0)[0])[0]
print(newmarg)

# As seurat, log1p transform the data (inplace)
np.log1p(matrix.data, out=matrix.data)

# Finally, scale data and set mean 0 and sd 1:
# sc.pp.scale(adata, max_value=10)
matrix = matrix.T

# Save the matrix as cpickle
def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return(matrix)

region = re.sub("_" + suffix,"", h5file)
# Preprocessed as cp.gz:
ppfile = region + '_preprocessed_' + re.sub("h5","", suffix) + "cp.gz"






# with gzip.open(ppfile, 'wb') as outfile:
#     pickle.dump(matrix, outfile, protocol=2)

ngene = matrix.shape[0]
n_comps = 50
chunksize = 1e4
ncell = matrix.shape[0]
ngene = matrix.shape[1]
pca_ = IncrementalPCA(n_components=n_comps)
X_pca = np.zeros((ncell, n_comps), matrix.dtype)
nchunk = int(np.floor(ncell / chunksize) + 1)
fullcov = np.zeros((ngene, ngene), matrix.dtype)

# print(i)
# start = int(i * chunksize)
# end = int(np.min([(i+1) * chunksize, ncell]))
# chunk = matrix[start:end,:]
# print("Converting to dense")
# chunk = chunk.toarray() if sparse.issparse(chunk) else chunk

ipca = IncPCA(2, 1.0)
m = 2  # number of new points


print(i)
start = int(i * chunksize)
end = int(np.min([(i+1) * chunksize, ncell]))
chunk = matrix[start:end,:]
print("Converting to dense")
chunk = chunk.toarray() if sparse.issparse(chunk) else chunk
print("PCA step")

import time
t1 = time.time()
cv = np.cov(chunk.T)
print(time.time() - t1)


# As dense takes:
t1 = time.time()
xtx = chunk.T.dot(chunk)
print(time.time() - t1)


# As sparse takes:
chunk = matrix[start:end,:]
t1 = time.time()
xtx = chunk.T.dot(chunk)
print(time.time() - t1)

# Full as sparse takes:
t1 = time.time()
xtx = matrix.T.dot(matrix)
print(time.time() - t1)


ipca.partial_fit(chunk)


# process 20 x m new points
for i in range(20):
    ipca.partial_fit(chunk)

Y_a = ipca.transform(chunk)

for i in range(nchunk):
    print(i)
    start = int(i * chunksize)
    end = int(np.min([(i+1) * chunksize, ncell]))
    chunk = matrix[start:end,:]
    print("Converting to dense")
    chunk = chunk.toarray() if sparse.issparse(chunk) else chunk
    print("PCA step")
    ipca.partial_fit(chunk)







for i in range(nchunk):
    print(i)
    start = int(i * chunksize)
    end = int(np.min([(i+1) * chunksize, ncell]))
    chunk = matrix[start:end,:]
    print("Converting to dense")
    chunk = chunk.toarray() if sparse.issparse(chunk) else chunk
    print("PCA step")
    pca_.partial_fit(chunk)


for i in range(nchunk):
    print(i)
    start = int(i * chunksize)
    end = int(np.min([(i+1) * chunksize, ncell]))
    chunk = matrix[start:end,:]
    chunk = chunk.toarray() if sparse.issparse(chunk) else chunk
    X_pca[start:end] = pca_.transform(chunk)



# for chunk, _, _ in adata_comp.chunked_X(chunk_size):
#     chunk = chunk.toarray() if issparse(chunk) else chunk
#     pca_.partial_fit(chunk)

# for chunk, start, end in adata_comp.chunked_X(chunk_size):
#     chunk = chunk.toarray() if issparse(chunk) else chunk
#     X_pca[start:end] = pca_.transform(chunk)









# yes | conda create -n loom_env -c conda-forge loompy h5py numpy
# conda activate loom_env
import loompy
import numpy as np
import h5py

loomfile = '~/DEVTRAJ/db/variants/ALZ/counts/merged_velocyto/ALZ_full_dataset_6k_smartseq2.loom'
loomfile = 'ALZ_full_dataset_6k_smartseq2.loom'


ds = loompy.connect(loomfile, validate=False)
types = ds.layers.keys()

print(types)

h5file = 'ALZ_full_dataset_6k_smartseq2.hdf5'
hf = h5py.File(h5file, 'w')

for mtype in types:
    print(mtype)
    sp = ds[mtype]
    if mtype == '':
        mtype = '_'
    print(sp.shape)
    cells = sp.ds.col_attrs['CellID']
    genes = sp.ds.row_attrs['Gene']
    mat = sp[:,:]
    hf.create_dataset(mtype, data=mat,
                      compression='gzip', compression_opts=9)
    del(mat)

hf.create_dataset('genes', data=genes,
                  compression='gzip', compression_opts=9)
hf.create_dataset('cells', data=cells,
                  compression='gzip', compression_opts=9)
hf.close()
ds.close()

gs = np.array([str(c) for c in genes])
ct = np.array([str(c) for c in cells])


with open('ALZ_full_dataset_6k_smartseq2_genes.txt', 'w') as f:
    for line in gs:
        f.write(str(line) + "\n")


with open('ALZ_full_dataset_6k_smartseq2_cells.txt', 'w') as f:
    for line in ct:
        f.write(str(line) + "\n")


# --------------------------------------------
import pandas as pd
import pickle
import numpy as np
import gzip

dffile = "pairs_ENH_ovl_df.cp.gz"

def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return(matrix)

# Write out each of the chromosomes separately
df = load_pickle_gzip_latin(dffile)
df['ind'] = df.index + 1

chromlist = pd.unique(df.chrom)

for chrom in chromlist:
    print(chrom)
    sdf = df[df.chrom == chrom]
    sdf = sdf[['gene', 'tss','mind','ind']]
    chromfile = chrom + '_pairs_ENH_ovl_df.tsv.gz'
    sdf.to_csv(chromfile, index=False, sep="\t")




