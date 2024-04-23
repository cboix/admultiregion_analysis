#!usr/bin/python
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# Using data loader + modules + graph handlers
# Created: 06/16/21
# -------------------------------------------------------------
import sys
bindir = '/home/cboix/data/DEVTRAJ/bin/'
sys.path.append(bindir + 'multiRegion/modules')
sys.path.append(bindir + 'multiRegion/pseudotime')
from adata_loader import *
from compute_modules_handler import *

# -----------------------------------------------------
# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# Only if you don't already have an anndata object you want to use for this.
# ct = 'Ast'
ct = 'Mic_Immune'
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
    sc.tl.pca(scdata.adata, n_comps=100)
    # Save:
    scdata.adata.write(h5ad_file)

# mod.adata.obs['msex'] = mod.adata.obs['msex'].astype('category')

outdir = scdata.datadir + 'mlp_gp_pred/'
imgdir = '/home/cboix/data/DEVTRAJ/img/multiRegion/modules/'

# -------------------
# Module computation:
# -------------------
mod = handler_ZCA_module(scdata.adata, csuff=scdata.csuff,
        filter_expr=0.05, z=4.5,
        use_v=True,
        imgdir=imgdir, use_fbpca=False)
mod.setup() # Filter - PCA - decorr - correlation
# mod.compute_umap()
suff = 'vvt'

suff = 'vvt_50'
ind = np.arange(50)
cv = mod.V[ind,:].T.dot(mod.V[ind,:]) / mod.U.shape[0]
sd = np.sqrt(np.diag(cv))
cv = cv / sd[:,np.newaxis]
mod.corr_w = cv / sd[np.newaxis,:]

# Make graph:
mod.make_graph(suff, cutoff=.5, resolution=2)
# mod.plot_graph(suff, attr='leiden', show_labels=True, w=18)
# mod.plot_graph(suff, attr='leiden', show_labels=False, w=18)
mod.plot_graph(suff, attr='leiden', show_labels=True,
        frac_labels=0.25, adj_txt=True, w=18)



# Plots for slides:
mod.plot_gene_logfc('vvt', attr='celltype', show_labels=False, adj_txt=False, w=16, fc_cut=2)
mod.plot_gene_logfc('vvt', attr='nrad', show_labels=False, adj_txt=False, w=16, fc_cut=.5)
mod.plot_gene_logfc('vvt', attr='msex', show_labels=False, adj_txt=False, w=16, fc_cut=2)

# Get the modules/print out:
mlist = mod.get_modules(suff, print_modules=True)


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

# GSEAPY
# enrichr

# NOTE: Can't run from internal...
mlist = mod.get_modules(suff, print_modules=False)
gpres = {}
for ll in mlist.keys(): #las:
    print(ll)
    testlist = mlist[ll].tolist()
    print(" ".join(testlist))
    gpres[ll] = gp.profile(organism='hsapiens', query=testlist)
    # lmat[:,ll] = np.mean(mod.adata[:,lnam].X, axis=1)


# 

suff = 'vvt_' + covar
mod.plot_gene_logfc(suff, attr=covar, show_labels=True, adj_txt=False, w=20, fc_cut=1)

# --------------------------------
# Design a rotation matrix for SVD 
# to contain a specific covariate:
# --------------------------------
covar = 'niareagansc'
cvcol = mod.adata.obs[covar]
cvcol = cvcol.to_numpy()

# Coeff fit for U s.t. it best fits column:
vcorrcoef(mod.U.T, cvcol)
coeff = mod.U.T.dot(cvcol) / mod.s

r = mod.U.dot(coeff * mod.s)
np.corrcoef(cvcol, r)

np.sum(mod.U, axis=0)

mod.U






# ------------------------------------

mod.make_graph('u', cutoff=.3, k=None, scale=None, layout_method='fr')
# mod.make_graph('u', cutoff=.3, k=50, scale=None, layout_method='fr')
mod.plot_graph('u', covariate=None, show_labels=True, adj_txt=False, w=20)
mod.plot_graph('u', covariate='leiden', show_labels=True, adj_txt=False, w=18)
# mod.plot_graph('u', covariate='leiden', show_labels=True, adj_txt=True, w=18)

# Test on AD / non-AD:
mod.adata.obs['nrad'] = 'AD'
mod.adata.obs.loc[mod.adata.obs.niareagansc > 2, 'nrad'] = 'CTRL'
mod.plot_gene_logfc('u', attr='nrad', show_labels=True, adj_txt=False, w=20, fc_cut=.5)

mod.plot_gene_logfc('u', attr='celltype', show_labels=True, adj_txt=False, w=20)


corr = mod.graphs['u'].corr
gmarg = mod.graphs['u'].gmarg
aw, ind = adj_from_bivariate_cutoff(corr, gmarg)

aw= mod.graphs['u'].aw


# -------------------------------------------------------
# Finding appropriate cutoffs for the correlation matrix:
# -------------------------------------------------------
# p-values (highly inflated, need n-effective)
from scipy.stats import beta
n = mod.Xw.shape[1]
n = 100
dist = beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
# The p-value returned by pearsonr is a two-sided p-value.
# For a given sample with correlation coefficient r, the p-value is the probability that abs(r’) of a random sample x’ and y’ drawn from the population with zero correlation would be greater than or equal to abs(r). In terms of the object dist shown above, the p-value for a given r and length n can be computed as:
r = 0.45
# r = 0.2
p = 2*dist.cdf(-abs(r))


# ---------------------------------------------
# Plotting statistics on the univariate margin:
# ---------------------------------------------
# Look at distr here:
cw = mod.corr_w
cw = cw - np.diag(np.diag(cw))
# Margins
gmarg = np.mean(mod.adata.X > 0, axis=0)
gvar = np.var(mod.adata.X, axis=0)
cvar = np.var(cw, axis=0)
cmax = np.max(cw, axis=0)
cmed = np.median(cw, axis=0)
cmean = np.mean(cw, axis=0)
cq = np.quantile(cw, .99, axis=0)

from scipy.interpolate import UnivariateSpline
rng = np.random.default_rng()
# Use the default value for the smoothing parameter:
cord = np.argsort(gmarg)
xs = np.linspace(0.01, 1, 100)
splmean = UnivariateSpline(gmarg[cord], cmean[cord], k=2)
splvar = UnivariateSpline(gmarg[cord], cvar[cord], k=2)

# Plot corr by marg:
fig, axs = plt.subplots(1, 4, figsize=(24,6))
axs[0].scatter(gmarg, cmax, s=5)
axs[0].set_title('Margin of genes vs. max corr', fontsize=24)
axs[1].scatter(gmarg, cvar, s=5)
axs[1].plot(xs, splvar(xs), 'g', lw=3)
axs[1].set_title('Margin of genes vs. var corr', fontsize=24)
axs[2].scatter(gmarg, cmean, s=5)
axs[2].set_title('Margin of genes vs. mean corr', fontsize=24)
axs[2].plot(xs, splmean(xs), 'g', lw=3)
axs[3].scatter(gmarg, cmean, s=5)
axs[3].set_title('Margin of genes vs. 99% quantile', fontsize=24)
plt.tight_layout()
plt.savefig("marg_vs_corr.png")
plt.close()


# Cutoff by zscore:
ymean = splmean(xs)
ysd = np.sqrt(splvar(xs))

ymean = splmean(gmarg)
ysd = np.sqrt(splvar(gmarg))

z = 4
ycut = ymean + ysd * z # Appropriate cutoffs


# Plot corr by marg:
plt.figure(figsize=(12,12))
fig, axs = plt.subplots(1, 3, figsize=(24,8))
axs[0].scatter(gvar * gmarg, cmax, s=5)
axs[0].set_title('Var of genes vs. max corr', fontsize=24)
axs[1].scatter(gvar * gmarg, cmed, s=5)
axs[1].set_title('Var of genes vs. median corr', fontsize=24)
axs[2].scatter(np.log10(gvar * gmarg), cq, s=5)
axs[2].set_title('Var of genes vs. 99% quantile', fontsize=24)
plt.tight_layout()
plt.savefig("varmarg_vs_corr.png")
plt.close()


cord = np.argsort(gmarg)
spl = UnivariateSpline(np.log10(gmarg[cord]), cmax[cord], k=1)
xs = np.log10(np.linspace(0.01, 1, 100))
ys = spl(xs)

# Plot corr by marg:
plt.figure(figsize=(12,12))
fig, axs = plt.subplots(1, 3, figsize=(24,8))
axs[0].scatter(np.log10(gmarg), cmax, s=5)
axs[0].plot(xs, spl(xs), 'g', lw=3)
axs[0].set_title('log(margin) of genes vs. max corr', fontsize=24)
axs[1].scatter(np.log10(gmarg), cmed, s=5)
axs[1].set_title('log(margin) of genes vs. median corr', fontsize=24)
axs[2].scatter(np.log10(gmarg), cq, s=5)
axs[2].set_title('log(margin) of genes vs. 99% quantile', fontsize=24)
plt.tight_layout()
plt.savefig("logmarg_vs_corr.png")
plt.close()


# -----------------------------------------------------
# Plotting statistics on the bivariate pair of margins:
# -----------------------------------------------------
from scipy.interpolate import SmoothBivariateSpline
from mpl_toolkits.mplot3d import Axes3D

# Correlation matrix:
cw = mod.corr_w
cw = cw - np.diag(np.diag(cw)) # Probably should just not eval on diag (set NA?)
aw = sparse.coo_matrix(cw)

lmarg = np.log10(gmarg)

# TODO: Increased computational performance:
# Bin genes by their margin:
xs = np.linspace(-3, 0, 50)
rinds = np.digitize(lmarg[aw.row], xs)
cinds = np.digitize(lmarg[aw.col], xs)
awdf = pd.DataFrame({'x':rinds, 'y':cinds, 'dat':aw.data})

# Aggregate, compute the mean and variances in bins:
mdf = awdf.groupby(["x","y"]).mean().reset_index()
sdf = awdf.groupby(["x","y"]).std().reset_index()

# Fit the spline to the bivariate data:
splmean = SmoothBivariateSpline(x=mdf.x, y=mdf.y, z=mdf.dat, kx=2, ky=2)
splsd = SmoothBivariateSpline(x=sdf.x, y=sdf.y, z=sdf.dat, kx=2, ky=2)

# Plot z-score splines

# Predict on full data, cutoff given z-scores
xloc = np.sort(np.unique(mdf.x.to_numpy()))
smean = splmean(xloc, xloc)
ssd = splsd(xloc, xloc)
se = 0.01
ssd[ssd < se] = se

# Plot mean and variance splines:
xd = np.diff(xs)
xc = 10**(xs[:-1] + xd)[xloc-1]
xa, ya = np.meshgrid(xc, xc)

import scipy.stats as st
pcut = 0.05 # Set bonferroni testing per gene:
z = - st.norm.ppf(pcut / aw.shape[0]) # If z is none.

zcut = smean + ssd * z

fig, axs = plt.subplots(nrows=1, ncols=3, subplot_kw={'projection': '3d'},
        figsize=(20,6))
axs[0].plot_wireframe(xa, ya, smean, color='k')
axs[0].set_title('Margin of genes vs. mean corr', fontsize=24)
axs[1].plot_wireframe(xa, ya, ssd, color='k')
axs[1].set_title('Margin of genes vs. sd corr', fontsize=24)
axs[2].plot_wireframe(xa, ya, zcut, color='k')
axs[2].set_title('Margin of genes vs. cutoff, z=' + str(round(z,1)), fontsize=24)
plt.tight_layout()
plt.savefig("marg_2d_vs_corr.png")
plt.close()

# Fix: 
zcoo = sparse.coo_matrix(zcut)
zdf = pd.DataFrame({'x':xloc[zcoo.row], 'y':xloc[zcoo.col], 'cutoff':zcoo.data})
# TO multi-index

# Merge cutoffs against pre-cut:
min_cut = np.min(zcoo.data)
aind = aw.data >= min_cut
awcut_df = pd.DataFrame({'x':rinds[aind], 'y':cinds[aind],
    'r': aw.row[aind], 'c': aw.col[aind], 'dat':aw.data[aind]})
awcut_df = awcut_df.merge(zdf)
awcut_df = awcut_df[awcut_df.dat >= awcut_df.cutoff]

aw2 = sparse.coo_matrix((awcut_df.dat, (awcut_df.r, awcut_df.c)), shape=aw.shape)
aw = aw2.tocsr()

rcut = 1
rind = np.array((np.sum(aw2, axis=1) > rcut).T)[0]
cind = np.array((np.sum(aw2, axis=0) > rcut))[0]
ind = rind + cind
aw = aw[ind, :]
aw = aw[:, ind]



tuples = list(zip(zdf['x'].to_numpy(),
    zdf['y'].to_numpy()))

zdf.index = pd.MultiIndex.from_tuples(tuples, names=["x", "y"])


# TODO: Increased computational performance:
# Bin genes by their margin:


# Re-create aw adjacency:


# Reduce to binned, calc ymarg, y



cord = np.argsort(gmarg[aw.row])
xdat = gmarg[aw.row][cord]
ydat = gmarg[aw.col][cord]
zdat = aw.data[cord]


# Fit the spline to the bivariate data:
interp_spline = SmoothBivariateSpline(
        x=xdat, y=ydat, z=zdat)


# Margins
gmarg = np.mean(mod.adata.X > 0, axis=0)
gvar = np.var(mod.adata.X, axis=0)
cvar = np.var(cw, axis=0)
cmax = np.max(cw, axis=0)
cmed = np.median(cw, axis=0)
cmean = np.mean(cw, axis=0)
cq = np.quantile(cw, .99, axis=0)

from scipy.interpolate import UnivariateSpline
rng = np.random.default_rng()
# Use the default value for the smoothing parameter:
cord = np.argsort(gmarg)
xs = np.linspace(0.01, 1, 100)
splmean = UnivariateSpline(gmarg[cord], cmean[cord], k=2)
splvar = UnivariateSpline(gmarg[cord], cvar[cord], k=2)



# Regularly-spaced, coarse grid
xs = np.linspace(0.01, 1, 100)
ys = np.linspace(0.01, 1, 100)

X, Y = np.meshgrid(x, y)
Z = np.exp(-(2*X)**2 - (Y/2)**2)

interp_spline = RectBivariateSpline(y, x, Z)



# Regularly-spaced, fine grid
dx2, dy2 = 0.16, 0.16
x2 = np.arange(-xmax, xmax, dx2)
y2 = np.arange(-ymax, ymax, dy2)
X2, Y2 = np.meshgrid(x2,y2)
Z2 = interp_spline(y2, x2)

fig, ax = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': '3d'})
ax[0].plot_wireframe(X, Y, Z, color='k')

ax[1].plot_wireframe(X2, Y2, Z2, color='k')
for axes in ax:
    axes.set_zlim(-0.2,1)
    axes.set_axis_off()

fig.tight_layout()
plt.show()





self.aw = self.corr * (self.corr > self.cutoff)
self.aw = self.aw - np.diag(np.diag(self.aw))
self.aw = sparse.csr_matrix(self.aw)
rind = np.array((np.sum(self.aw, axis=0) > rcut))[0]
