#!usr/bin/python
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# Using data loader + modules + graph handlers
# Created: 06/16/21 // Updated: 08/10/21
# -------------------------------------------------------------
import sys
bindir = '/home/cboix/data/DEVTRAJ/bin/'
sys.path.append(bindir + 'multiRegion/pseudotime')
from importlib import reload
from adata_loader import *
import scdemon as sm


# -----------------------------------------------------
# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# Only if you don't already have an anndata object you want to use for this.
ct = 'Ast'
# ct = 'Mic_Immune'
# ct = 'Vasc_Epithelia'
# st = 'Ast GRM3'
st = None
scdata = scDataLoader(ct, st,
                      usecombat=False, usehvg=False,
                      remove_mt=False, filter_TH=False)
h5ad_file = scdata.datadir + 'preprocessed_formodules_' + scdata.csuff + '.h5ad'
if os.path.exists(h5ad_file):
    scdata.adata = sc.read_h5ad(h5ad_file)
else:
    scdata.load()
    # Prep observed covariates (project specific):
    scdata.adata.obs['nrad'] = 'AD'
    scdata.adata.obs.loc[scdata.adata.obs.niareagansc > 2, 'nrad'] = 'CTRL'
    scdata.adata.obs['cogdxad'] = 'AD'
    scdata.adata.obs.loc[scdata.adata.obs.cogdx < 4, 'cogdxad'] = 'CTRL'
    scdata.adata.obs['cogdxad'] = scdata.adata.obs['cogdxad'].astype('category')
    scdata.adata.obs['region'] = scdata.adata.obs['region'].astype('category')
    scdata.adata.obs['Apoe_e4'] = scdata.adata.obs['Apoe_e4'].astype('category')
    scdata.adata.obs['projid'] = scdata.adata.obs['projid'].astype('category')
    scdata.adata.obs['msex'] = scdata.adata.obs['msex'].astype('category')
    sc.tl.pca(scdata.adata, n_comps=100)
    # Save:
    scdata.adata.write(h5ad_file)

# All glial:
if ct == 'Glial':
    prefstr = scdata.datadir + 'preprocessed_formodules_'
    ast_adata = sc.read_h5ad(prefstr + 'Ast.h5ad')
    mic_adata = sc.read_h5ad(prefstr + 'Mic_Immune.h5ad')
    opc_adata = sc.read_h5ad(prefstr + 'Opc.h5ad')
    oli_adata = sc.read_h5ad(prefstr + 'Oli.h5ad')
    vce_adata = sc.read_h5ad(prefstr + 'Vasc_Epithelia.h5ad')
    # Concatenate all:
    adata = anndata.concat((ast_adata, mic_adata, opc_adata,
                            oli_adata, vce_adata))
    del(ast_adata, mic_adata, opc_adata, oli_adata, vce_adata)
    gc.collect()

outdir = scdata.datadir + 'mlp_gp_pred/'
imgdir = '/home/cboix/data/DEVTRAJ/img/multiRegion/modules/'

# -------------------
# Module computation:
# -------------------
# ct = 'Per'
# adata = scdata.adata[scdata.adata.obs.celltype == ct]
# mod = sm.scdemon(adata, csuff=ct, imgdir=imgdir,
csuff = 'Glial' if ct == 'Glial' else scdata.csuff
adata = adata if ct == 'Glial' else scdata.adata
mod = sm.scdemon(adata, csuff=csuff, imgdir=imgdir,
                  # estimate_sd=True,
                  svd_k=100, filter_expr=0.1, calc_raw=False)
# TODO: Ensure mod.adata.X doesn't turn into ArrayView !!
mod.setup() # Filter - PCA - decorr - correlation

import scipy.stats as st
pcut = 0.01 # TODO: get proper BY testing correction equation
mod.z = -st.norm.ppf(pcut / mod.corr_w.shape[0]) # If z is none.

suff = 'vvt'
mod.make_graph(suff, resolution=2.5, rcut=1)
mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)
mod.plot_graph(suff, attr='leiden', show_labels=False, w=16)
# mod.plot_graph(suff, attr='leiden', show_labels=False, w=16)

def find_module(gene, mlist):
    for key in mlist.keys():
        if gene in mlist[key]:
            print(key, " ".join(mlist[key]))

find_module('HSF1',mlist)


mod.z = -1
suff = 'adj'
mod.corr_w, t_hat = mod.cobj.get_adj_correlation()

# As if corr:
used = s_ij
X_sd = np.sqrt(np.diag(used))
cv = used / X_sd[:,np.newaxis]
cv2 = cv / X_sd[np.newaxis,:]

mod.corr_w = cv2
mod.z = 4.5

mod.make_graph(suff, resolution=2, rcut=2)# , cutoff=t_hat)
# mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)
mod.plot_graph(suff, attr=None, show_labels=True, w=16)


# Copied:
from scipy import sparse
from scipy.interpolate import UnivariateSpline, SmoothBivariateSpline
from sklearn.linear_model import LinearRegression
from scipy.stats import norm

def triu_mask(A):
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:,None] < r
    return A[mask]

# Baseline:
# ---------
X = mod.adata.X.copy()
N = X.shape[0]
mu = np.mean(X, axis=0)
X_adj = X - mu[None,:]
# s_ij = X.T @ X / N
s_ij = X_adj.T @ X_adj / N

X_adj2 = X_adj ** 2
E1 = X_adj2.T.dot(X_adj2) / N
E2 = (s_ij)**2
theta = E1 - E2
t_ij = s_ij / np.sqrt(theta)
t_ij = t_ij - np.diag(np.diag(t_ij))
t_ij = t_ij / np.std(t_ij)

# mod.cobj.t_ij_adj = t_ij / np.std(t_ij)

# zc = sparse.coo_matrix(np.triu(mod.cobj.t_ij))
# zc = sparse.coo_matrix(np.triu(t_ij))
zc = sparse.coo_matrix(np.triu(t_ij_adj))
sdf = pd.DataFrame({'p': mod.cobj.gmarg[zc.row], 'q': mod.cobj.gmarg[zc.col],
                    'r': zc.row, 'c': zc.col, 'dat': zc.data})
sdf['pq'] = np.sqrt(sdf.p * sdf.q)

Z = mod.adata.X.copy()
Z[Z > 0] = 1
N = Z.shape[0]
pij = Z.T.dot(Z) / N
pij = pij - np.diag(np.diag(pij))
pij = sparse.coo_matrix(np.triu(pij))

corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
zc = sparse.coo_matrix(corr_w)

pdf = pd.DataFrame({'r': pij.row, 'c': pij.col, 'pij': pij.data})
sdf['pij'] = pdf.pij

# Investigate top-tail
tdf = sdf[(sdf.pq > 0.92)]
tdf['g1'] = mod.genes[tdf.r]
tdf['g2'] = mod.genes[tdf.c]
tdf = tdf.reset_index()
tdf = tdf.loc[np.argsort(-tdf.dat),:]

# Bin and adjust:
# xs = np.linspace(np.min(sdf.pq), np.max(sdf.pq), 50)
# sdf['d_pq'] = np.digitize(sdf.pq, xs)
# u, c = np.unique(sdf['d_pq'], return_counts=True)
xs = np.linspace(np.min(sdf.pij), np.max(sdf.pij), 50)
sdf['d_pij'] = np.digitize(sdf.pij, xs)
u, c = np.unique(sdf['d_pij'], return_counts=True)
NS = len(u)
# Calculate the means/sdevs:
smean = np.zeros(NS)
ssd = np.zeros(NS)
for i, ui in enumerate(u):
    ci = c[i]
    # rind = (sdf['d_pq'] == ui)
    rind = (sdf['d_pij'] == ui)
    vec = sdf.dat[rind]
    smean[i] = np.mean(vec)
    ssd[i] = np.std(vec)


# Bin2d and adjust:
sqmarg = np.sqrt(mod.gmarg)
xs = np.linspace(np.min(sqmarg), np.max(sqmarg), 50)
sdf['d_p'] = np.digitize(np.sqrt(sdf.p), xs)
sdf['d_q'] = np.digitize(np.sqrt(sdf.q), xs)
u, c = np.unique(sdf['d_p'], return_counts=True)
NS = len(u)

# Calculate the means/sdevs:
scount = np.zeros((NS, NS))
smean = np.zeros((NS, NS))
ssd = np.zeros((NS,NS))
for i, ui in enumerate(u):
    idf = sdf[sdf.d_p == ui]
    for j, uj in enumerate(u):
        jdf = idf[idf.d_q == uj]
        scount[i,j] = len(jdf)
        smean[i,j] = np.mean(jdf.dat)
        ssd[i,j] = np.std(jdf.dat)

ssd[np.isnan(ssd)] = 1
smean[np.isnan(smean)] = 0

# Past halfway, set mean to flat behavior (othw skewed by diff. alt)
# ind = np.where(xs > 0.5)[0][0]
# smean[-ind:,:] = 0
# smean[:,-ind:] = 0

# TODO: Log ?? sqrt??

smean = sparse.coo_matrix(smean)
scount = sparse.coo_matrix(scount + 1)
ssd = sparse.coo_matrix(ssd)
# Fit the spline to the bivariate data:
print("Fitting + predicting spline")
overall_mean = np.mean(smean.data)
overall_sd = np.mean(ssd.data)
print('mean and sd:', overall_mean, overall_sd)
ks = 2
splmean = SmoothBivariateSpline(x=u[smean.row], y=u[smean.col],
                                z=smean.data - overall_mean, kx=ks, ky=ks)
splsd = SmoothBivariateSpline(x=u[ssd.row], y=u[ssd.col],
                              z=ssd.data - overall_sd, kx=1, ky=1)
xloc = np.arange(len(xs))
pmean = splmean(u, u)
pmean = pmean + overall_mean
psd = splsd(xloc, xloc)
psd = psd + overall_sd
pvar = psd **2

# Linear regression for the variance (due to null model):
regr = LinearRegression()
mids = (np.array([0] + list(xs[:-1])) + xs) / 2
X = np.vstack([ssd.row,ssd.col]).T
y = ssd.data
regr.fit(X, y, np.sqrt(scount.data))
sdf['pred_var'] = regr.predict(sdf[['d_p','d_q']]) ** 2
# sdf['pred_var'] = 1 / np.sqrt(sdf.pq) * np.sqrt(NK / NG)
# sdf['pred_var'] = 1 / np.sqrt(sdf.p * sdf.q + np.sqrt(sdf.p * sdf.q * (1 - sdf.p) * (1 - sdf.q)))* np.sqrt(NK / NG)
# sdf['pred_var'] = 1 / np.sqrt(sdf.p * sdf.q + np.sqrt(sdf.p * sdf.q * (1 - sdf.p) * (1 - sdf.q)))
# pq - pq^2 + p^2 * q + p * q^2 + 2 * p * q * pq
val = sdf.pij - sdf.pij**2 + sdf.p**2 * sdf.q + sdf.p * sdf.q**2 + 2 * sdf.p * sdf.q * sdf.pij

# sdf['pred_var'] = sdf.pq * -5.21 + 5.26

sdf['pred_mean'] = pmean[sdf.d_p-1, sdf.d_q-1]
# sdf['pred_var'] = pvar[sdf.d_p-1, sdf.d_q-1]
# sdf['dat2'] = (sdf['dat'] - sdf['pred_mean']) / sdf['pred_var']
sdf['dat2'] = (sdf['dat'] - sdf['pred_mean'])


i = 0
j = 3
# sns.histplot(X_z[:,i] * X_z[:,j])
x = X[:,i] * X[:,j]
sns.histplot(x[x != 0])
plt.tight_layout()
plt.savefig('test.png')
plt.close()

i = 7167
j = 7169
x = X[:,i] * X[:,j]
sns.histplot(x)
plt.tight_layout()
plt.savefig('test2.png')
plt.close()

i = 0
j = 4
X = mod.adata.X.copy()
# x = X_z[:,i] * X_z[:,j]
x = X[:,i] * X[:,j]
sns.histplot(x[x != 0])
plt.tight_layout()
plt.savefig('test3.png')
plt.close()

print(np.std(x))
print(np.mean(x))

# Binomial testing:
# -----------------
Z = mod.adata.X.copy()
Z[Z > 0] = 1
N = Z.shape[0]
pij = Z.T.dot(Z) / N
p = np.mean(Z, axis=0)
pq = p[:, None] * p[None, :]
# del(Z)
gc.collect()

from scipy.stats import binom

i = 7167
print(i, p[i])
t = binom.ppf(.99999, N, pq[i,:])
np.sum(pij[i,:] * N > t)


mod.z = - 1
suff = 'binom'
mod.corr_w = (pij * N) / t
t_hat = 1.0

mod.make_graph(suff, resolution=2, rcut=2, cutoff=t_hat)
# mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)
mod.plot_graph(suff, attr=None, show_labels=True, w=16)








# ---------------------------------------------------
# Past halfway, set mean to flat behavior (othw skewed by diff. alt)
# TODO: Find better way to do this
ind = np.where(xs > 0.5)[0][0]
smean[-ind:] = smean[ind]

# Linear regression for the variance (due to null model):
regr = LinearRegression()
mids = (np.array([0] + list(xs[:-1])) + xs) / 2
regr.fit(mids[:,None], ssd, np.sqrt(c))
sdf['pred_var'] = regr.predict(sdf.pq[:,None]) ** 2

# Univariate spline for the mean (inflated values):
spl = UnivariateSpline(mids[:,None], smean, w=np.sqrt(c), k=2)
sdf['pred_mean'] = spl(sdf.pq[:, None])

# Adjust - pretty close to what we want.
sdf['dat2'] = (sdf['dat'] - sdf['pred_mean']) / sdf['pred_var']
sdf['dat2'] = sdf['dat2'] / np.std(sdf['dat2'])
# Turn back into matrix-like
cu = sparse.coo_matrix((sdf.dat2, (sdf.r, sdf.c)), shape=mod.cobj.t_ij.shape)
cu = cu.toarray()
cu = cu + cu.T
mod.cobj.t_ij_adj = cu # Temporary fix
# mod.cobj.t_ij_adj = t_ij_adj # cu # Temporary fix

# TODO move this to graph.py ??
alpha = 0.01
# Valid search space (for which t_ij ~ N(0,1))
ap = 2 * np.log(np.log(mod.cobj.P))
bp = np.sqrt(4 * np.log(mod.cobj.P) - ap)

# Data vector:
vec = triu_mask(mod.cobj.t_ij_adj)
# Search broadly:
t = np.array(list(np.arange(0, bp,.1)) + [bp])
est_null = (2 - 2 * norm.cdf(t)) * len(vec)
sig = np.array([np.sum(vec > x) for x in t])
efdr = est_null / sig
tind = np.where(efdr < alpha)[0][0]
# TODO: Set as 2 sqrt(log(p)) if none.

# Search second round.
t2 = np.arange(t[tind-1], t[tind], 0.01)
est_null = (2 - 2 * norm.cdf(t2)) * len(sdf)
sig = np.array([np.sum(vec > x) for x in t2])
efdr = est_null / sig
t2ind = np.where(efdr < alpha)[0][0]
mod.cobj.t_hat = t2[t2ind]
print(mod.cobj.t_hat)

mod.corr_w = mod.cobj.t_ij_adj
t_hat = mod.cobj.t_hat

mod.z = - 1
suff = 'adj'
# mod.corr_w, t_hat = mod.cobj.get_adj_correlation()

mod.make_graph(suff, resolution=2, rcut=2, cutoff=6)
# mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)
mod.plot_graph(suff, attr=None, show_labels=True, w=16)




# Average correlation / raw, because why not.
# mod.corr_w = (1 * corr_w + 1 * mod.corr_raw) / 2

suff = 'vvt'
# mod.z=5
mod.make_graph(suff, resolution=2, rcut=1)
mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)

genes = np.array(['GRM3','SLC1A3','SLC1A2','GFAP','APOE','HSPA1A','USP9Y', 'ZFY', 'UTY'])
gind = mod.genes.isin(genes)
mod.genes[gind]
np.argmax(mod.cobj.V[:,gind], 0)

np.sum(mod.cobj.V[:,gind], 0)

U = mod.cobj.U
V = mod.cobj.V
xc = V[0,:][:, np.newaxis].dot(V[0,:][np.newaxis,:]) * mod.cobj.s[0]
X_sd = np.sqrt(np.diag(xc))
cv = xc / X_sd[:,np.newaxis]
cv2 = cv / X_sd[np.newaxis,:]

# Test establishing variance of estimates:
U = mod.cobj.U
V = mod.cobj.V
data = mod.adata.X[mod.adata.X > 0]

# Semi-simulate data:
def gen_data(pi):
    hi = 1.0 * (np.random.uniform(size=mod.adata.shape[0]) <= pi)
    # zi = np.random.choice(data[:100000], size=int(np.sum(hi)))
    # zi = np.random.normal(size=int(np.sum(hi)))
    zi = np.random.gamma(1, size=int(np.sum(hi)))
    hi[hi == 1] = zi
    return(hi)

def upper_tri_masking(A):
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:,None] < r
    return A[mask]

s_m2 = np.diag(1 / mod.cobj.s **2)

N = 1000
prob = .99
arr = np.array([gen_data(p) for p in [prob]*N])

XtUS = arr.dot(U).dot(s_m2)
UtX = U.T.dot(arr.T)
VVt = XtUS.dot(UtX) * V.shape[1] / U.shape[0]
# VVt = arr.dot(arr.T) * V.shape[1] / U.shape[0]

X_sd = np.sqrt(np.diag(VVt))
cv = VVt / X_sd[:,np.newaxis]
cv2 = cv / X_sd[np.newaxis,:]
cv2 = cv2 - np.diag(np.diag(cv2))

vals = upper_tri_masking(cv2)
print(np.mean(vals), np.std(vals))

cv3 = np.corrcoef(arr)
vals = upper_tri_masking(cv3)
print(np.mean(vals), np.std(vals))

# 0.05
# 0.02 and 0.045

# 0.05
# 0.73 and 0.016

gmarg = mod.gmarg
corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
cm = np.max(corr_w, 1)
cm = np.std(corr_w, 1)

import seaborn as sns
from matplotlib import pyplot as plt

fig = plt.figure(figsize=(6,6))
plt.scatter(np.sqrt(mod.gmarg), cm, s=3)
plt.tight_layout()
plt.savefig('test.png')
plt.close()

import seaborn as sns
from matplotlib import pyplot as plt
from scipy import sparse
zc = sparse.coo_matrix(corr_w)

zdf = pd.DataFrame({'p': gmarg[zc.row], 'q': gmarg[zc.col],
                         'r': zc.row, 'c': zc.col, 'dat': zc.data})

ind = np.random.randint(len(zdf), size=1000000)

fig = plt.figure(figsize=(15,15))
ax = plt.gca()
# plt.scatter(np.sqrt(zdf.p[ind] * zdf.q[ind]), zdf.dat[ind], s=1)
plt.scatter(np.sqrt(zdf.p * zdf.q), zdf.dat, s=1)
ax.set_xlabel('sqrt(sparsity of i * sparsity of j)')
ax.set_ylabel('Estimated correlation')
plt.tight_layout()
plt.savefig('corr_vs.png')
plt.close()

# Instead (int / n)/ (sqrt(pq))
import gc
from scipy import sparse
Z = mod.adata.X.copy()
Z[Z > 0] = 1
Zint = Z.T.dot(Z)
del(Z)
gc.collect()

pij = Zint / mod.adata.shape[0]
pij = pij / np.sqrt(mod.gmarg[:, np.newaxis])
pij = pij / np.sqrt(mod.gmarg[np.newaxis, :])
pij = pij - np.diag(np.diag(pij))
pij = sparse.coo_matrix(pij)

corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
zc = sparse.coo_matrix(corr_w)

gmarg = mod.gmarg
zdf = pd.DataFrame({'p': gmarg[zc.row], 'q': gmarg[zc.col],
                         'r': zc.row, 'c': zc.col, 'dat': zc.data})
pdf = pd.DataFrame({'r': pij.row, 'c': pij.col, 'pij': pij.data})

# TODO: Fix to be robust:
zdf['pij'] = pdf.pij

import seaborn as sns
from matplotlib import pyplot as plt

fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(zdf.pij, zdf.dat, s=1)
ax.set_xlabel('p_ij / sqrt(p_i * p_j)')
ax.set_ylabel('Estimated correlation')
plt.tight_layout()
plt.savefig('corr_vs_pij_int.png')
plt.close()

# Correlation of the intersection
zdf['ij_int'] = zdf.pij * np.sqrt(zdf.p * zdf.q)
zdf['pq_cv'] = (zdf.ij_int - zdf.p * zdf.q) / np.sqrt(zdf.q * zdf.p * (1 - zdf.q) * (1 - zdf.p))

fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(zdf.pq_cv, zdf.dat, s=1)
ax.set_xlabel('pq_cv')
ax.set_ylabel('Estimated correlation')
plt.tight_layout()
plt.savefig('corr_vs_pij_diff.png')
plt.close()

sdf = zdf[zdf.pij > .9]
sdf['g1'] = mod.genes[sdf.r]
sdf['g2'] = mod.genes[sdf.c]
sdf = sdf.reset_index()
sdf = sdf.loc[np.argsort(-sdf.dat),:]

ind = np.where(mod.genes == 'APOE')[0][0]
sdf = zdf[zdf.r == ind]
sdf['g1'] = mod.genes[sdf.r]
sdf['g2'] = mod.genes[sdf.c]
sdf = sdf.reset_index()
sdf = sdf.loc[np.argsort(-sdf.dat),:]
sdf.head(30)

# Bin and fit a linear model through variance; plot cutoffs:
# ----------------------------------------------------------
xs = np.linspace(np.min(zdf.pij), np.max(zdf.pij), 25)
# xs = np.linspace(np.min(zdf.pij), .9, 25)
zdf.d_pij = np.digitize(zdf.pij, xs)

# Calculate the means/sdevs:
u, c = np.unique(zdf.d_pij, return_counts=True)
NS = len(u)
smean = np.zeros(NS)
ssd = np.zeros(NS)
sq = np.zeros(NS)
for i,xi in enumerate(u):
    rind = (zdf.d_pij == xi)
    ci = zdf.dat[rind]
    smean[i] = np.mean(ci)
    ssd[i] = np.std(ci)
    sq[i] = np.quantile(ci, 1-0.0005)

fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(zdf.pij, zdf.dat, s=1)
# plt.plot(xs, ssd * 4.5, 'r+-')
plt.plot(xs, sq, 'r+-')
plt.plot(xs, c**.2 / 100, 'y*-')
ax.set_xlabel('p_ij / sqrt(p_i * p_j)')
ax.set_ylabel('Estimated correlation')
plt.tight_layout()
plt.savefig('corr_vs_pij_int_bins.png')
plt.close()




smean = sparse.coo_matrix(smean)
ssd = sparse.coo_matrix(ssd)
# Fit the spline to the bivariate data:
print("Fitting + predicting spline")
overall_mean = np.mean(smean.data)
overall_sd = np.mean(ssd.data)
print('mean and sd:', overall_mean, overall_sd)
ks = 2
splmean = SmoothBivariateSpline(x=u[smean.row], y=u[smean.col], z=smean.data - overall_mean, kx=ks, ky=ks)
splsd = SmoothBivariateSpline(x=u[ssd.row], y=u[ssd.col], z=ssd.data - overall_sd, kx=ks, ky=ks)
xloc = np.arange(len(xs))
# TODO: Improve the absolute numerical shrinkage of mean and SD towards mean:
smean = splmean(xloc, xloc)
smean = smean + overall_mean # - (np.sign(smean) * 0.1 * overall_mean)
ssd = splsd(xloc, xloc)
ssd = ssd + overall_sd # - (np.sign(ssd) * 0.1 * overall_sd)


# Alternatively, calculate the variance of the covariance estimator:
# ------------------------------------------------------------------
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import sparse

N = mod.cobj.U.shape[0]
NK = mod.cobj.U.shape[1]
NG = mod.cobj.V.shape[1]

X_z = mod.cobj.U.dot(mod.cobj.V)
# s = np.diag(mod.cobj.s)
# X_z = mod.cobj.U.dot(s).dot(mod.cobj.V)
# mu_z = np.mean(X_z, axis=1) # W/ or without reduction
# X_z = (X_z - mu_z[:, None])
X_z2 = X_z ** 2
s_ij = X_z.T.dot(X_z) / N
E1 = X_z2.T.dot(X_z2) / N
E2 = (s_ij)**2
N = X_z.shape[0]
theta = E1 - E2
# NOTE: Without S, must add this scaling factor (as before)
t_ij = s_ij / np.sqrt(theta) * np.sqrt(NK) # NG or NK? Need to double check...
# \--> NG for sparsity? NK for internal sum
# t_ij = s_ij / np.sqrt(theta)
t_ij = t_ij - np.diag(np.diag(t_ij))

Z = mod.adata.X.copy()
Z[Z > 0] = 1
N = Z.shape[0]
pij = Z.T.dot(Z) / N
pij = pij - np.diag(np.diag(pij))
pij = sparse.coo_matrix(np.triu(pij))
pq = p[:, None] * p[None, :]
pij = Z.T.dot(Z) / N
p = np.mean(Z, axis=0)
p2 = p**2

stat = pij - pij**2 + \
    p2[:,None] * p[None,:] + \
    p[:,None] * p2[None,:] + \
    2 * p[:,None] * p[None,:] * pij
# stat2 = (pij - p[:, None] * p[None,:])

# Test adj:
t_ij_adj = t_ij * stat / pij


np.mean(t_ij)
np.std(t_ij)

ind = np.unravel_index(np.argmax(t_ij), t_ij.shape)
mod.cobj.corr_w[ind]
mod.genes[ind[0]] # Makes sense... so not losing info
mod.genes[ind[1]]

tri = np.triu(t_ij)
tri = tri - np.diag(np.diag(tri))
vec = tri[tri != 0]
np.mean(vec)
np.std(vec)

# Example: Get one cutoff at FDR = 5%:
from scipy.stats import norm
NT = (NG**2 - NG) / 2 # Number of tests

t = np.arange(2, 5.25, .25)
est_null = (2 - 2 * norm.cdf(t)) * NT
sig = np.array([np.sum(vec > i) for i in t])
efdr = est_null / sig * 100
tind = np.where(efdr < 5)[0][0]
print('cutoff:', t[tind])
print('nsig:', sig[tind])

# Plotting this statistic:
# ------------------------
zc = sparse.coo_matrix(np.triu(t_ij))
zc = sparse.coo_matrix(np.triu(t_ij_adj))
sdf = pd.DataFrame({'p': mod.gmarg[zc.row], 'q': mod.gmarg[zc.col],
                    'r': zc.row, 'c': zc.col, 'dat': zc.data})
sdf['pq'] = np.sqrt(sdf.p * sdf.q)

fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(t_ij, corr_w, s=1)
ax.set_xlabel('Test statistic')
ax.set_ylabel('Estimated correlation')
plt.tight_layout()
plt.savefig('corr_vs_tij.png')
plt.close()

fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(np.sqrt(sdf.p * sdf.q), sdf.dat, s=1)
ax.set_xlabel('sqrt(sparsity of i * sparsity of j)')
ax.set_ylabel('Test statistic')
plt.tight_layout()
plt.savefig('tij_vs_sparsity.png')
plt.close()

# Adjust the test statistic across the bins:
# TODO: adjustment
# ------------------------------------------
from scipy.stats import norm
# Reduce to tested stats:
sdf['pq'] = np.sqrt(sdf.p * sdf.q)

NBIN = 50
xs = np.linspace(np.min(sdf.pq), np.max(sdf.pq), NBIN)
sdf['d_pq'] = np.digitize(sdf.pq, xs)

# Calculate the means/sdevs:
u, c = np.unique(sdf['d_pq'], return_counts=True)
NS = len(u)

smean = np.zeros(NS)
ssd = np.zeros(NS)
for i, ui in enumerate(u):
    ci = c[i]
    rind = (sdf['d_pq'] == ui)
    vec = sdf.dat[rind]
    smean[i] = np.mean(vec)
    ssd[i] = np.std(vec)
    # TODO properly change vals:

# Past halfway, set mean to flat behavior (othw skewed by
ind = np.where(xs > 0.5)[0][0]
smean[-ind:] = smean[ind]

from scipy.interpolate import UnivariateSpline
from sklearn.linear_model import LinearRegression

# Linear regression for the variance (due to null model):
regr = LinearRegression()
mids = (np.array([0] + list(xs[:-1])) + xs) / 2
regr.fit(mids[:,None], ssd, np.sqrt(c))
sdf['pred_var'] = regr.predict(sdf.pq[:,None]) ** 2

sdf[sdf.pred_var > 1.3]

# Univariate spline for the mean (inflated values):
spl = UnivariateSpline(mids[:,None], smean, w=np.sqrt(c), k=2)
sdf['pred_mean'] = spl(sdf.pq[:, None])
# regr.fit(mids[:,None], smean, np.sqrt(c))
# sdf['pred_mean'] = regr.predict(sdf.pq[:,None]) ** 2

# Adjust - pretty close to what we want.
sdf['dat2'] = (sdf['dat'] - sdf['pred_mean']) / sdf['pred_var']

# Search space:
ap = 2 * np.log(np.log(NG))
bp = np.sqrt(4 * np.log(NG) - ap)

# FDR:
alpha = 0.01
t = np.array(list(np.arange(2, bp,.1)) + [bp])
est_null = (2 - 2 * norm.cdf(t)) * len(sdf)
sig = np.array([np.sum(sdf['dat2'] > x) for x in t])
efdr = est_null / sig
tind = np.where(efdr < alpha)[0][0]

# Two rounds, but local FDR is probably best.
t2 = np.arange(t[tind-1], t[tind], 0.01)
est_null = (2 - 2 * norm.cdf(t2)) * len(sdf)
sig = np.array([np.sum(sdf['dat2'] > x) for x in t2])
efdr = est_null / sig
t2ind = np.where(efdr < alpha)[0][0]
t_hat = t2[t2ind]
print(t_hat)

# Plot adjusted:
fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(sdf.pq, sdf.dat, s=1)
t_hat = 0.4
plt.plot([0,1], [t_hat, t_hat], 'r-', lw=2)
plt.plot([0,1], [-t_hat, -t_hat], 'r-', lw=2)
ax.set_xlabel('sqrt(sparsity of i * sparsity of j)')
ax.set_ylabel('Test statistic')
plt.tight_layout()
plt.savefig('tij_vs_sparsity.png')
plt.close()

# Plot adjusted:
fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(sdf.pq, sdf.dat2, s=1)
plt.plot([0,1], [t_hat, t_hat], 'r-', lw=2)
plt.plot([0,1], [-t_hat, -t_hat], 'r-', lw=2)
ax.set_xlabel('sqrt(sparsity of i * sparsity of j)')
ax.set_ylabel('Test statistic')
plt.tight_layout()
plt.savefig('tij_adj_vs_sparsity.png')
plt.close()

# Plot adjusted:
fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(np.sqrt(sdf.pij), sdf.dat, s=1)
plt.plot([0,1], [t_hat, t_hat], 'r-', lw=2)
plt.plot([0,1], [-t_hat, -t_hat], 'r-', lw=2)
ax.set_xlabel('sqrt(p_ij)')
ax.set_ylabel('Test statistic')
plt.tight_layout()
plt.savefig('tij_vs_pij.png')
plt.close()

# Plot adjusted:
fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(np.sqrt(sdf.pij), sdf.dat2, s=1)
plt.plot([0,1], [t_hat, t_hat], 'r-', lw=2)
plt.plot([0,1], [-t_hat, -t_hat], 'r-', lw=2)
ax.set_xlabel('sqrt(p_ij)')
ax.set_ylabel('Test statistic')
plt.tight_layout()
plt.savefig('tij_adj_vs_pij.png')
plt.close()



zdf = sdf[(sdf.p < 0.2) & (sdf.q < 0.2)]
# zdf = sdf[(sdf.p < 0.2)]
t_hat = mod.cobj.t_hat

# Plot adjusted:
fig = plt.figure(figsize=(15,15))
ax = plt.gca()
plt.scatter(zdf.pq, zdf.dat2, s=1)
plt.plot([0,1], [t_hat, t_hat], 'r-', lw=2)
plt.plot([0,1], [-t_hat, -t_hat], 'r-', lw=2)
ax.set_xlabel('sqrt(sparsity of i * sparsity of j)')
ax.set_ylabel('Test statistic')
plt.tight_layout()
plt.savefig('tij_adj_vs_sparsity_lowlow.png')
plt.close()


# Look at top genes:
topind = np.argsort(-sdf.dat2)[0:200]
sdf = sdf.reset_index()
zdf = sdf.loc[np.array(topind)]
zdf['g1'] = mod.genes[zdf.r]
zdf['g2'] = mod.genes[zdf.c]

gene = 'LSS'
gind = np.where(mod.genes == gene)[0][0]
gdf = sdf[(sdf.r == gind) | (sdf.c == gind)]
gdf = gdf[np.abs(gdf.dat2) > 4.2]
gdf['g1'] = mod.genes[gdf.r]
gdf['g2'] = mod.genes[gdf.c]

# TODO: Turn back into correlation matrix + plot as network.
cu = sparse.coo_matrix((sdf.dat2, (sdf.r, sdf.c)), shape=(NG, NG))
cu = cu.toarray()
cu = cu + cu.T
mod.corr_w = cu # Temporary fix

mod.z = - 1
suff = 'adj'
mod.make_graph(suff, resolution=2, rcut=1, cutoff=t_hat)
# mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)
mod.plot_graph(suff, attr=None, show_labels=True, w=16)





# Plotting estimated SD and other data:
# TODO: Clean up
# -------------------------------------
mod.adata.var['gmarg'] = np.mean(mod.adata.X > 0, axis=0)
gmarg = mod.adata.var['gmarg'][mod.genes].to_numpy()

fig = plt.figure(figsize=(5,5))
ax = plt.gca()
plt.scatter(np.log10(gmarg), np.log10(X_sd), s=4)
ax.set_xlabel('log10(% cells express gene)')
ax.set_ylabel('Estimated standard deviation')
plt.tight_layout()
plt.savefig('margin_vs_est_sd.png')
plt.close()

X_sd[gind]

V = mod.cobj.V
ind = np.random.choice(np.arange(V.shape[0]), size=V.shape[0], replace=False)

V = mod.cobj.V
xc = V.T.dot(V)
X_sd = np.sqrt(np.diag(xc))

xc = V[ind,:].T.dot(V)
cv = xc / X_sd[:,np.newaxis]
cv2 = cv / X_sd[np.newaxis,:]

np.argmax(X_sd)



# Varying s_red changes:
# s_pow = 0.01 # s_pow at 2 is unbalanced / wrong
for s_pow in [0.01, 0.1, 0.25, 0.5, 1, 1.5, 2]:
    s_red = np.diag(mod.cobj.s**s_pow) # For comparison with Vs^2V.T (approx cov.)
    X_cov = mod.cobj.V.T.dot(s_red).dot(mod.cobj.V) / mod.cobj.U.shape[0] * mod.cobj.V.shape[1]
    X_sd = np.sqrt(np.diag(X_cov))
    cv = X_cov / X_sd[:,np.newaxis]
    mod.corr_w = cv / X_sd[np.newaxis,:]
    suff = 's_' + str(s_pow)
    # Make graph:
    mod.make_graph(suff, resolution=2, rcut=3)
    mod.plot_graph(suff, attr='leiden', show_labels=True, w=16)


# --------------------------------
# Calculate base graphs + modules:
# --------------------------------
import scipy.stats as st
pcut = 0.05 # TODO: get proper BY testing correction equation
mod.z = -st.norm.ppf(pcut / mod.corr_w.shape[0]) # If z is none.
print(round(mod.z, 2))

suff = 'vvt'
# Make graph:
mod.make_graph(suff, resolution=2, rcut=3)
# mod.plot_graph(suff, attr=None, show_labels=True, w=18)
mod.plot_graph(suff, attr='leiden', show_labels=True, w=18)

# mod.plot_graph(suff, attr='leiden', show_labels=True,
#         frac_labels=0.25, adj_txt=True, w=18)

# Get the modules/print out:
mlist = mod.get_modules(suff, print_modules=False)
mod.save_modules(suff)
mod.plot_umap(suff, attr='leiden')

suff = 'raw'
# Make graph:
mod.make_graph(suff, resolution=2)
mod.plot_graph(suff, attr=None, show_labels=True, w=18)
mod.plot_graph(suff, attr='leiden', show_labels=True,
        frac_labels=0.25, adj_txt=True, w=18)

# Get the modules/print out:
mlist = mod.get_modules(suff, print_modules=False)
mod.save_modules(suff)

mod.plot_umap(suff, attr='leiden')

# TODO: SE as covar (fixed, Xvar and Yvar (diff), and n=k?
# V <- ((n-1)^2/n^3)*(Xvar*Yvar - covar^2) + ((n-1)/n^3)*(covar^2 - Xvar*Yvar) # NOT CORRECT

# ----------------------------------------------------------
# Look at the genes with high corr in raw but low in decorr:
# ----------------------------------------------------------
mod.adata.var['gmarg'] = np.mean(mod.adata.X > 0, axis=0)
genes = np.array(['GRM3','SLC1A3','SLC1A2','GFAP','VIM','EMP1', 'CLIC4', 'CD44','IFITM3','C3', 'APOE','CST3','HLA-A'])
mod.adata.var.loc[genes,:]
# mod.adata.var.loc['GBP2',:]
mod.adata.var.loc[['GRM3','CHRDL1','WIF1','ETV5'],:]

keptgenes = mod.graphs[suff].labels
np.array([[x, x in keptgenes] for x in genes])

corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
corr_z = corr_w / mod.corr_sd
ind = np.unravel_index(np.argsort(-corr_z, axis=None)[:300], corr_z.shape)
arr_ind = np.array(ind).T
mod.genes[arr_ind]
np.sort(-corr_z, axis=None)[:300]

gind = mod.genes.isin(genes)
cz = corr_z[gind,:]
np.max(cz, axis=1)
mod.genes[gind]
mod.genes[np.argmax(cz, axis=1)]


gind = mod.genes.isin(genes)
corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
np.max(corr_w[gind,:], axis=1)



gw = mod.graphs[suff].gw
names = np.array(gw.vs['name'])
gene = 'APOE'
[names[i] for i in gw.neighbors(np.where(names == gene)[0][0])]
# mod.adata.var.loc[keptgenes,:]

ind = np.array([np.where(mod.genes == x)[0][0] for x in genes])

[list(mod.genes[np.argsort(-mod.corr_w[ind[i],:])][0:10]) for i in range(len(genes))]
[[genes[i]] + list(np.round(np.sort(-mod.corr_w[ind[i],:]),2)[0:10]) for i in range(len(genes))]

[list(mod.genes[np.argsort(-mod.corr_raw[ind[i],:])][0:10]) for i in range(len(genes))]
[[genes[i]] + list(np.round(np.sort(-mod.corr_raw[ind[i],:]),2)[0:10]) for i in range(len(genes))]

mod.corr_w[ind[:, np.newaxis],ind]

gmarg = np.array(np.mean(mod.adata.X > 0, axis=0).T)
gmarg_cut = 0.5
gmarg[gmarg > gmarg_cut] = gmarg_cut
lmarg = np.log10(np.array(gmarg))
lmarg[lmarg < -3] = -3
xs = np.linspace(np.min([np.min(lmarg), -3]), np.max(lmarg), 50)
dmarg = np.digitize(lmarg, xs)
u,c = np.unique(dmarg, return_counts=True)

import seaborn as sns
from matplotlib import pyplot as plt
sns.histplot(mod.adata.var.loc[keptgenes,:].gmarg)
plt.tight_layout()
plt.savefig('test.png')
plt.close()

sns.histplot(mod.adata.var.gmarg)
plt.tight_layout()
plt.savefig('test.png')
plt.close()

import seaborn as sns
from matplotlib import pyplot as plt

gmarg = np.array(np.mean(mod.adata.X > 0, axis=0).T)
gmarg_cut = 0.75 # For stability
gmarg[gmarg > gmarg_cut] = gmarg_cut # Remove if no thresh.
lmarg = np.log10(np.array(gmarg))
lmarg[lmarg < -3] = -3
# xs = np.linspace(np.min([np.min(lmarg), -3]), 0, 50)
# xs = np.linspace(np.min([np.min(lmarg), -3]), np.max(lmarg), 50)
xs = np.linspace(np.min(lmarg), np.max(lmarg), 25)
dmarg = np.digitize(lmarg, xs)
u,c = np.unique(dmarg, return_counts=True)


zcut = mod.graphs[suff].zcut
fig = plt.figure(figsize=(16,16))
sns.heatmap(np.round(zcut,2), annot=True)
plt.tight_layout()
plt.savefig('test.png')
plt.close()


evs = mod.cobj.s
plt.scatter(np.arange(len(evs)), evs / np.sum(evs))
plt.tight_layout()
plt.savefig('test.png')
plt.close()




# --------------------------------
# Plot UMAP with leiden from each:
# --------------------------------
import umap
from matplotlib import pyplot as plt
from matplotlib import rcParams

suff = 'vvt'

gname = mod.graphs[suff].gw.vs['name']
gind = np.array([np.where(mod.genes == x)[0][0] for x in gname])
if suff == 'raw':
    subcorr = mod.corr_raw[gind[:,np.newaxis], gind]
else:
    subcorr = mod.corr_w[gind[:,np.newaxis], gind]

uw = umap.UMAP()
umat = uw.fit_transform(subcorr)
vscols = mod.graphs[suff].colors['leiden']

# Plot each:
plt.figure(figsize=(12,12))
plt.scatter(umat[:, 0], umat[:, 1], color=vscols, s=8)
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP of whitened, centered', fontsize=16)
plt.tight_layout()
plt.savefig(mod.imgdir + 'umap_leiden_' + mod.csuff + '_' + suff + '.png')
plt.close()



# --------------------------------
# Plots for slides:
suff = 'vvt'
mod.plot_gene_logfc(suff, attr='celltype', show_labels=False, adj_txt=False, w=16, fc_cut=2)
mod.plot_gene_logfc(suff, attr='nrad', show_labels=False, adj_txt=False, w=16, fc_cut=.5)
# mod.plot_gene_logfc(suff, attr='msex', show_labels=False, adj_txt=False, w=16, fc_cut=2)

# Plot correlation of SVD components with covariates:
cvlist = ['celltype', 'region', 'niareagansc', 'cogdx', 'cogdxad', 'nft','plaq_n','plaq_d','Apoe_e4', 'msex','pmi','age_death','n_genes','n_counts']
mod.plot_svd_corr(cvlist)

# Plot the average expression of leiden modules on covariates:
cvlist = ['celltype', 'region', 'niareagansc', 'cogdx', 'Apoe_e4', 'msex']
mod.plot_heatmap_avgexpr(suff, cvlist=cvlist, attr='leiden')

mod.make_subset_graph(suff, 'msex', cutoff=0.5, plot=True)

mod.plot_graph('vvt_msex', attr='leiden', show_labels=True, adj_txt=False, w=18)

# Test on AD / non-AD:
mod.plot_gene_logfc('vvt', attr='nrad', show_labels=True, adj_txt=False, w=20, fc_cut=.5)


# Also run gprofiler:
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)
# Alternatively, use: GSEAPY with enrichr()

# NOTE: Can't run from internal...
mlist = mod.get_modules(suff, print_modules=False)
gpres = {}
for ll in mlist.keys(): #las:
    print(ll)
    testlist = mlist[ll].tolist()
    print(" ".join(testlist))
    gpres[ll] = gp.profile(organism='hsapiens', query=testlist)
    # lmat[:,ll] = np.mean(mod.adata[:,lnam].X, axis=1)

