#!usr/bin/python
# -------------------------------------------------------------
# Compute co-expression modules within python scanpy framework,
# Using data loader + modules + graph handlers
# Created: 06/16/21 // Updated: 08/10/21
# -------------------------------------------------------------
import sys
import re
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype, is_numeric_dtype
from scipy import sparse
from scipy.stats import mannwhitneyu
import scdemon as sm
from adata_handling import adata_loader
from gprofiler import GProfiler

# For plotting:
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import rcParams
from statannot import add_stat_annotation
import statsmodels.api as sma
from statsmodels.formula.api import ols
import umap
import textwrap

# Directories:
maindir = '/home/cboix/data/DEVTRAJ/'
imgdir = maindir + 'img/multiRegion/modules/'
datadir = maindir + 'db/multiRegion/'
pref = 'all_brain_regions_filt_preprocessed_scanpy'

# -----------------------------------------------------
# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
# TODO: UPDATE - adata_loader is updated
ct = 'Mic_Immune'
st = None
ctstr = re.sub("_", "/", ct)
ststr = None if st is None else re.sub("_", "/", st)
scdata = adata_loader(ct, st,
                      usecombat=False, usehvg=False,
                      remove_mt=False, filter_TH=False)
h5ad_file = scdata.datadir + 'preprocessed_formodules_' + \
    scdata.csuff + '.h5ad'
adata = sc.read_h5ad(h5ad_file)
csuff = scdata.csuff

# TODO: Subset neurons to EC, TH, HC (by max repr.) and all CORTICAL separately
subset_EC = True
if subset_EC and ct == 'Exc':
    csuff = csuff + "_ECneurons"
    ecneu = ['Exc AGBL1 GPC5', 'Exc RELN GPC5',
             # 'Exc TRPC6 ANO2', 'Exc COBLL1 UST',
             'Exc RELN COL5A2', 'Exc COL25A1 SEMA3D',
             'Exc DLC1 SNTG2', 'Exc SOX11 NCKAP5',
             'Exc TOX3 POSTN', 'Exc TOX3 INO80D', 'Exc TOX3 TTC6']
    ind = np.isin(adata.obs.celltype, ecneu)
    adata = adata[ind, :]
    # Load UMAP specific to EC neurons:
    metafile = datadir + 'metadata_test_mrad_Exc_EC_combat_filthvg.tsv'
    mdf = pd.read_csv(metafile, sep='\t')
else:
    # Load the UMAP and filtering metadata:
    metafile = datadir + pref + '_norm.final_noMB.cell_labels.tsv.gz'
    mdf = pd.read_csv(metafile, sep='\t')
    mdf = mdf.loc[mdf['major.celltype'] == ctstr, :]
    if st is not None:
        mdf = mdf.loc[mdf.cell_type_high_resolution == ststr, :]

# Annotate adata object with UMAP coords (TODO: add other metadata?):
mdf.index = mdf.barcode
mdf2 = pd.merge(adata.obs, mdf, how='left',
                left_index=True, right_index=True)
adata.obsm['X_umap'] = mdf2[['U1', 'U2']].to_numpy()

# Make sure to keep only the ones on the UMAP (for EC in particular)
kind = np.where(-np.isnan(mdf2.U2))[0]
adata = adata[kind, :]

# ----------------------------------------
# Compute modules for the anndata dataset:
# ----------------------------------------
# TODO: Ensure mod.adata.X doesn't turn into ArrayView !!
mod = sm.scdemon(adata, csuff=csuff, imgdir=imgdir,
                  seed=1, svd_k=100, filter_expr=0.1, calc_raw=False)
mod.setup()


cor = np.corrcoef(mod.cobj.X.T)
mod.corr_raw = cor

# Create graph + plot:
suff = 'vvt'
# suff = 'raw'; mod.z = 5
mod.make_graph(suff, resolution=2.5, rcut=1)
mod.plot_graph(suff, show_labels=True, w=16)
mod.plot_graph(suff, attr='leiden', show_labels=False, w=16)

# Get modules:
mlist = mod.get_modules(suff, print_modules=False)
mod.save_modules(suff)
# mod.plot_umap(suff, attr='leiden')


def find_module(gene, mlist):
    for key in mlist.keys():
        if gene in mlist[key]:
            print(key, " ".join(mlist[key]))
            return(key)


k = find_module('NCDN', mlist)
k = find_module('H1FX', mlist)
k = find_module('CDK5R1', mlist)
k = find_module('PRNP', mlist)

# Score the cells by the modules and plot:
mod.plot_umap_grid(suff)

# -----------------------------------------------------
# Module enrichment in any ranked gene list / DEG list:
# ----------------------------------------------------
method = 'nebula_ruv'
# method = 'mast'
reg = 'allregions'
dedir = datadir + "dereg/"
ctpref = ct + "_"
st = 'Mic' if ct == 'Mic_Immune' else st
# st = 'CPEC'
ctpref = ctpref + ct if st is None else ctpref + st
ctpref = ctpref + "_" + reg
depref = dedir + method + "." + ctpref + "_"
# depref = dedir + method + "." + 'Exc_Exc_TOX3_TTC6_EC' + "_"

pathlist = ['nrad', 'nft', 'plaq_n', 'plaq_d', 'cogdxad']
dflist = {}
for path in pathlist:
    # Load and format data:
    defile = depref + path + ".tsv.gz"
    if os.path.exists(defile):
        df = pd.read_csv(defile, sep="\t")
        df['gset'] = '--'
        df.loc[df.col == 2, 'gset'] = 'Up'
        df.loc[df.col == 1, 'gset'] = 'Down'
        dflist[path] = df
        # Plot:
        suffix = method + "_" + path + "_" + reg + "_" + mod.csuff + "_" + suff
        title = ctpref + " (" + path + ")"
        mod.plot_df_enrichment(df, 'gset', suff=suff,
                               suffix=suffix, title=title)

suffix = method + "_allpath_" + reg + "_" + mod.csuff + "_" + suff
mod.plot_df_enrichment(dflist, 'gset', suff=suff, suffix=suffix)
mod.plot_df_enrichment(dflist, 'gset', suff=suff, suffix=suffix, ext='pdf')

k = find_module('THY1', mlist)
df[np.isin(df.gene, mlist[k])].head(50)

ndf = dflist['nrad'][['gene', 'logFC', 'gset']]
cdf = dflist['cogdxad'][['gene', 'logFC', 'gset']]
ndf.columns = ['gene', 'logFC_path', 'gset_path']
cdf.columns = ['gene', 'logFC_ci', 'gset_ci']
cdf = pd.merge(cdf, ndf)
cdf['in_m0'] = np.isin(cdf.gene, mlist[0])
subdf = cdf[cdf.in_m0]

pd.crosstab(cdf.gset_ci, cdf.gset_path)
pd.crosstab(subdf.gset_ci, subdf.gset_path)

# ndf[np.isin(df.gene, mlist[k])]
x = subdf.logFC_path.to_numpy()
y = subdf.logFC_ci.to_numpy()
lbl = subdf.gene.to_numpy()

fig = plt.figure(figsize=(8, 8))
sns.scatterplot(data=cdf, x='logFC_path', y='logFC_ci', hue='in_m0', s=5)
# add annotations one by one with a loop
for i in range(len(x)):
    # plt.text(df.x[line]+0.2, df.y[line], df.group[line], horizontalalignment='left', size='medium', color='black', weight='semibold')
    txt = plt.text(x[i], y[i], lbl[i])

ax = plt.gca()
ax.axhline(0, ls='--')
ax.axvline(0, ls='--')
plt.tight_layout()
plt.savefig('M0_nrad_cogdxad_comparison.png')
plt.close()


# Compare overall vs. within M0


# ----------------------------------------------------------
# Compare module activities across individuals / covariates:
# ----------------------------------------------------------
def make_tform(x, norm=True):
    unique, inverse = np.unique(x, return_inverse=True)
    onehot = np.eye(unique.shape[0])[inverse]
    if norm:
        onehot = onehot / np.sum(onehot, 0)[None, :]
    return(onehot, unique)


p = adata.obs.projid.to_numpy().astype(str)
r = adata.obs.region.to_numpy().astype(str)
adata.obs['pr'] = [p[i] + "_" + r[i] for i in range(len(r))]
tform, uq = make_tform(adata.obs.pr)

mapdf = adata.obs[['projid', 'region', 'pr', 'nrad', 'nft',
                   'plaq_n', 'plaq_d', 'cogdxad', 'Apoe_e4', 'cogdx']]
mapdf = mapdf.drop_duplicates()
mapdf = mapdf.reset_index()

# Average scores to dataframe:
avg_smat = smat.T @ tform
avg_smat = sparse.coo_matrix(avg_smat)
sdf = pd.DataFrame({'m': avg_smat.row,
                    'pr': uq[avg_smat.col],
                    'score': avg_smat.data})
sdf = pd.merge(sdf, mapdf)

# sns.histplot(mod.adata.var.loc[keptgenes,:].gmarg)
# sns.histplot(data=sdf[sdf.m == 0], hue='nrad', x='score')
w = 3
h = 2.5
nc, nr = 5, 5
fig, axs = plt.subplots(nrows=nr, ncols=nc,
                        gridspec_kw={'hspace': 0.1,
                                     'wspace': 0.05},
                        figsize=(w * nc, h * nr))
for i in range(5):
    for j in range(5):
        k = i * 5 + j
        ax = axs[i, j]
        ht = sns.scatterplot(
            data=sdf[sdf.m == k],
            x='nft', y='score', hue='nrad', ax=ax)
        ht = ax.set_title(k)

# plt.tight_layout()
plt.savefig('test.png')
plt.close()


sdf[(sdf.m == 0) & (sdf.score > 0.7)]

genes = ['GFAP', 'OSMR', 'CLIC4', 'CD44', 'ETV6']
score = np.mean(adata[:, genes].X, 1)
avg_score = score @ tform

adf = pd.DataFrame({'pr': uq, 'score': avg_score})
adf = pd.merge(adf, mapdf)

pd.crosstab(
    index=sdf.cogdxad,
    columns=sdf.nrad,
    values=sdf.score,
    aggfunc=np.mean)
pd.crosstab(index=sdf.cogdxad, columns=sdf.nrad, values=sdf.score, aggfunc=len)

bdf = sdf[sdf.m == 0]

hue_order = ['AD', 'CTRL']
box_pairs = [(('AD', 'AD'), ('AD', 'CTRL')),
             (('CTRL', 'AD'), ('CTRL', 'CTRL'))]

pal = {'AD': 'indianred', 'CTRL': 'royalblue'}
fig, axs = plt.subplots(nrows=1, ncols=2,
                        gridspec_kw={'hspace': 0.1, 'wspace': 0.25},
                        figsize=(10, 5))
sns.violinplot(data=sdf[sdf.m == 0], hue='cogdxad', x='nrad', y='score',
               ax=axs[0], palette=pal)
axs[0].set_title('M0 Score')
test_results = add_stat_annotation(
    axs[0], data=sdf[sdf.m == 0], x='nrad', y='score', hue='cogdxad',
    box_pairs=box_pairs, test='Mann-Whitney', text_format='star',
    loc='inside', verbose=2)
sns.violinplot(data=adf, hue='cogdxad', x='nrad', y='score', ax=axs[1],
               palette=pal)
axs[1].set_title('GFAP/OSMR/CLIC4/CD44/ETV6')
test_results = add_stat_annotation(
    axs[1], data=adf, x='nrad', y='score', hue='cogdxad',
    box_pairs=box_pairs, test='Mann-Whitney', text_format='star',
    loc='inside', verbose=2)
plt.tight_layout()
plt.savefig('M0_reactive_resilience_comparison.png')
plt.close()


adf[adf.score > 1.5]
pd.crosstab(adf.nrad, adf.cogdxad)
pd.crosstab(adf[adf.score > 1.5].nrad, adf[adf.score > 1.5].cogdxad)

# Fraction of cells with these genes:
# -----------------------------------
genes = ['GFAP', 'OSMR', 'CD44']
genes = ['GFAP', 'OSMR', 'CD44', 'SLC1A2']  # Kobayashi paper
score = (np.mean(adata[:, genes].X > 0, 1) == 1)
frac_cells = score @ tform
adf = pd.DataFrame({'pr': uq, 'frac_cells': frac_cells})
adf = pd.merge(adf, mapdf)
adf['path_ci'] = [adf.nrad[i] + "_" + adf.cogdxad[i] for i in range(len(adf))]

# box_pairs = [(('AD_AD','AD_CTRL'),('CTRL_CTRL','AD_CTRL')), (('CTRL_CTRL','AD_AD'),('CTRL','CTRL'))]
box_pairs = []
regs = pd.unique(adf.region)
for r in regs:
    box_pairs += [((r, 'AD_AD'), (r, 'AD_CTRL')),
                  ((r, 'CTRL_CTRL'), (r, 'AD_CTRL')),
                  ((r, 'CTRL_CTRL'), (r, 'AD_AD'))]

fig = plt.figure(figsize=(18, 7))
ax = plt.gca()
sns.violinplot(data=adf, hue='path_ci', x='region', y='frac_cells')
plt.tight_layout()
test_results = add_stat_annotation(
    ax, data=adf, x='region', y='frac_cells', hue='path_ci',
    box_pairs=box_pairs, test='Mann-Whitney', text_format='star',
    loc='inside', verbose=2)
plt.savefig('M0_reactive_fraction_cells.png')
plt.close()

model = ols('frac_cells ~ region + cogdxad', data=adf)
fitted_model = model.fit()
fitted_model.summary()


# ------------------------------------------------
# Compare module activities to dataset covariates:
# ------------------------------------------------

score_mat = mod.scores[suff]
covariates = adata.obs.columns
for covariate in covariates:
    iscat = is_categorical_dtype(mod.adata.obs[covariate])
    isnum = is_numeric_dtype(mod.adata.obs[covariate])
    print(covariate, iscat, isnum)


def zscore(x):
    return((x - np.mean(x)) / np.std(x))


covariate = ['nft', 'plaq_n', 'plaq_d']
covariate = 'plaq_n'
i = 0
NM = score_mat.shape[1]
coef = None
pdf = pd.DataFrame({'pmi': mod.adata.obs['pmi'],
                    'msex': mod.adata.obs['msex'],
                    'age_death': zscore(mod.adata.obs['age_death']),
                    'celltype': mod.adata.obs['celltype'],
                    'region': mod.adata.obs['region']})

if isinstance(covariate, str):
    pdf[covariate] = mod.adata.obs[covariate]
    lmstr = covariate
    pltstr = covariate
else:
    lmstr = " + ".join(covariate)
    pltstr = "_".join(covariate)
    for covar in covariate:
        pdf[covar] = mod.adata.obs[covar]


for i in range(NM):
    pdf['module'] = score_mat[:, i]
    model = ols('module ~ celltype + region + msex + pmi + age_death + ' +
                lmstr, data=pdf).fit()
    if is_categorical_dtype(pdf[covariate]):
        aov_table = sma.stats.anova_lm(model, typ=2)
        # TODO Log these:
        eff = aov_table['sum_sq'][5] / (np.sum(aov_table['sum_sq']))
        p = aov_table['PR(>F)'][5]
    else:
        eff = model.rsquared
        p = model.pvalues[covariate]
    if coef is None:
        coef = np.zeros((NM, len(model.params)))
        coef_names = model.params.index
    coef[i, :] = model.params
    print(i, np.round(eff, 4), np.round(-np.log10(p), 2))

fname = mod.imgdir + "modules_act_vs_" + pltstr + "_" + mod.csuff + "_" + suff + ".png"
fig = plt.figure(figsize=(1 + .75 * (len(coef_names) - 1), 12))
ht = sns.heatmap(coef[:, 1:], xticklabels=coef_names[1:], annot=True,
                 cmap='RdBu_r', cbar=True, center=0)
plt.tight_layout()
plt.savefig(fname)
plt.close()

# Correlation with each?
# Eta-squared for categorical?

# Remove highly correlated with celltype - but this won't remove all genes...
rmmod = (np.max(np.abs(coef[:, 1:4]), axis=1) > .2)
keep = np.where(rmmod == False)[0]

# Use only remaining genes:
genes = []
for k in keep:
    genes += mlist[k].tolist()

pstadata = mod.adata[:, genes]

sc.pp.combat(pstadata, key='projid')

sc.tl.pca(pstadata)
sc.pp.neighbors(pstadata)
sc.tl.umap(pstadata)

# Good start / benchmark for the pseudotime in astrocytes
sc.pl.umap(pstadata, color=['celltype', 'region', 'nrad'],
           save='_pst_umap_covars.png')


# TODO: Function to scdemon
# Vectorized correlation, matrix to vector:
def vcorrcoef(X, y):
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm)*(y-ym), axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2, axis=1)*np.sum((y-ym)**2))
    r = r_num/r_den
    return r


score_mat = np.zeros((len(adata), len(mlist)))
topgenes = [''] * len(mlist)
titles = [''] * len(mlist)
for i, key in enumerate(mlist.keys()):
    genes = mlist[key]
    subadata = adata[:, genes]
    score_mat[:, i] = np.mean(subadata.X, axis=1)
    cv = vcorrcoef(subadata.X.T, score_mat[:, i])
    topgenes[i] = ", ".join(genes[np.argsort(-cv)][0:3])
    titles[i] = "M" + str(key) + " (" + str(len(genes)) + \
        " genes): " + topgenes[i]
    print(i, topgenes[i])

# umap plotting function


def plot_umap(x, c, fname=None, ax=None, title=None, s=1, w=8, axlab=False):
    if ax is None:
        # Define plotting range:
        mx = np.max(x, axis=0)
        mn = np.min(x, axis=0)
        rn = mx - mn
        h = w * rn[1] / rn[0]
        fig = plt.figure(figsize=(w, h))
        ax = plt.gca()
    if title is not None:
        tw = textwrap.fill(title, 18)
        ax.set_title(tw, fontdict={'fontsize': 7})
    ax.set_facecolor('white')
    ax.scatter(x[:, 0], x[:, 1], s=s, c=c,
               marker='.', edgecolors='none',
               cmap=plt.get_cmap('viridis'))
    # Remove tick labels:
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    if axlab:
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
    if fname is not None:
        plt.tight_layout()
        plt.savefig(fname, dpi=350, bbox_inches='tight')
        plt.close()


def plot_umap_grid(x, cmat, nr, nc, fname, sel=None, titles=None, w=2, s=.5):
    # Define plotting range:
    mx = np.max(x, axis=0)
    mn = np.min(x, axis=0)
    rn = mx - mn
    h = w * rn[1] / rn[0]
    # Make subplots:
    fig, axs = plt.subplots(nrows=nr, ncols=nc,
                            gridspec_kw={'hspace': 0.01,
                                         'wspace': 0.01},
                            figsize=(w * nc, h * nr))
    nplot = cmat.shape[1] if sel is None else len(sel)
    # TODO: Auto determine grid NR/NC
    for i in range(nr):
        for j in range(nc):
            k = i * nc + j
            ax = axs[i, j]
            if k < nplot:
                k = k if sel is None else sel[k]
                title = None if titles is None else titles[k]
                plot_umap(x=x, c=cmat[:, k], title=title, s=s, ax=ax)
            else:
                ax.set_facecolor('white')
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.savefig(fname, dpi=350, bbox_inches='tight')
    plt.close()
    print(fname)


x = adata.obsm['X_umap']
# x = pstadata.obsm['X_umap']
if ct == 'Ast':
    ind = (x[:, 0] < 3) & (x[:, 1] > 5)  # For astrocytes, for now
else:
    ind = x[:, 0] < 100  # Default to all

fname = mod.imgdir + "module_umap_grid_" + mod.csuff + ".png"
# fname = mod.imgdir + "module_umap_pst_grid_" + mod.csuff + ".png"
plot_umap_grid(x=x[ind], cmat=score_mat[ind, :], fname=fname,
               nr=5, nc=8, w=2, titles=titles)

# Plot a selection:
sel = [15, 1, 7, 2, 4, 18, 9, 11, 14, 23, 0]
fname = mod.imgdir + "module_umap_grid_sel_" + mod.csuff + ".png"
plot_umap_grid(x=x[ind], cmat=score_mat[ind, :], sel=sel, fname=fname,
               nr=2, nc=6, w=2, titles=titles)


# ----------------------------
# Other base graphs + modules:
# ----------------------------
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
# V <- ((n-1)^2/n^3)*(Xvar*Yvar - covar^2) + ((n-1)/n^3)*(covar^2 -
# Xvar*Yvar) # NOT CORRECT

# ----------------------------------------------------------
# Look at the genes with high corr in raw but low in decorr:
# ----------------------------------------------------------
mod.adata.var['gmarg'] = np.mean(mod.adata.X > 0, axis=0)
genes = np.array(['GRM3',
                  'SLC1A3',
                  'SLC1A2',
                  'GFAP',
                  'VIM',
                  'EMP1',
                  'CLIC4',
                  'CD44',
                  'IFITM3',
                  'C3',
                  'APOE',
                  'CST3',
                  'HLA-A'])
mod.adata.var.loc[genes, :]
# mod.adata.var.loc['GBP2',:]
mod.adata.var.loc[['GRM3', 'CHRDL1', 'WIF1', 'ETV5'], :]

keptgenes = mod.graphs[suff].labels
np.array([[x, x in keptgenes] for x in genes])

corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
corr_z = corr_w / mod.corr_sd
ind = np.unravel_index(np.argsort(-corr_z, axis=None)[:300], corr_z.shape)
arr_ind = np.array(ind).T
mod.genes[arr_ind]
np.sort(-corr_z, axis=None)[:300]

gind = mod.genes.isin(genes)
cz = corr_z[gind, :]
np.max(cz, axis=1)
mod.genes[gind]
mod.genes[np.argmax(cz, axis=1)]


gind = mod.genes.isin(genes)
corr_w = mod.corr_w - np.diag(np.diag(mod.corr_w))
np.max(corr_w[gind, :], axis=1)


gw = mod.graphs[suff].gw
names = np.array(gw.vs['name'])
gene = 'APOE'
[names[i] for i in gw.neighbors(np.where(names == gene)[0][0])]
# mod.adata.var.loc[keptgenes,:]

ind = np.array([np.where(mod.genes == x)[0][0] for x in genes])

[list(mod.genes[np.argsort(-mod.corr_w[ind[i], :])][0:10])
 for i in range(len(genes))]
[[genes[i]] + list(np.round(np.sort(-mod.corr_w[ind[i], :]), 2)[0:10])
 for i in range(len(genes))]

[list(mod.genes[np.argsort(-mod.corr_raw[ind[i], :])][0:10])
 for i in range(len(genes))]
[[genes[i]] + list(np.round(np.sort(-mod.corr_raw[ind[i], :]), 2)[0:10])
 for i in range(len(genes))]

mod.corr_w[ind[:, np.newaxis], ind]

gmarg = np.array(np.mean(mod.adata.X > 0, axis=0).T)
gmarg_cut = 0.5
gmarg[gmarg > gmarg_cut] = gmarg_cut
lmarg = np.log10(np.array(gmarg))
lmarg[lmarg < -3] = -3
xs = np.linspace(np.min([np.min(lmarg), -3]), np.max(lmarg), 50)
dmarg = np.digitize(lmarg, xs)
u, c = np.unique(dmarg, return_counts=True)


sns.histplot(mod.adata.var.loc[keptgenes, :].gmarg)
plt.tight_layout()
plt.savefig('test.png')
plt.close()

sns.histplot(mod.adata.var.gmarg)
plt.tight_layout()
plt.savefig('test.png')
plt.close()


gmarg = np.array(np.mean(mod.adata.X > 0, axis=0).T)
gmarg_cut = 0.75  # For stability
gmarg[gmarg > gmarg_cut] = gmarg_cut  # Remove if no thresh.
lmarg = np.log10(np.array(gmarg))
lmarg[lmarg < -3] = -3
# xs = np.linspace(np.min([np.min(lmarg), -3]), 0, 50)
# xs = np.linspace(np.min([np.min(lmarg), -3]), np.max(lmarg), 50)
xs = np.linspace(np.min(lmarg), np.max(lmarg), 25)
dmarg = np.digitize(lmarg, xs)
u, c = np.unique(dmarg, return_counts=True)


zcut = mod.graphs[suff].zcut
fig = plt.figure(figsize=(16, 16))
sns.heatmap(np.round(zcut, 2), annot=True)
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

suff = 'vvt'

gname = mod.graphs[suff].gw.vs['name']
gind = np.array([np.where(mod.genes == x)[0][0] for x in gname])
if suff == 'raw':
    subcorr = mod.corr_raw[gind[:, np.newaxis], gind]
else:
    subcorr = mod.corr_w[gind[:, np.newaxis], gind]

uw = umap.UMAP()
umat = uw.fit_transform(subcorr)
vscols = mod.graphs[suff].colors['leiden']

# Plot each:
plt.figure(figsize=(12, 12))
plt.scatter(umat[:, 0], umat[:, 1], color=vscols, s=8)
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP of whitened, centered', fontsize=16)
plt.tight_layout()
plt.savefig(mod.imgdir + 'umap_leiden_' + mod.csuff + '_' + suff + '.png')
plt.close()


# --------------------------------
# Plots for slides:
suff = 'vvt'
mod.plot_gene_logfc(
    suff,
    attr='celltype',
    show_labels=False,
    adj_txt=False,
    w=16,
    fc_cut=2)
mod.plot_gene_logfc(
    suff,
    attr='nrad',
    show_labels=False,
    adj_txt=False,
    w=16,
    fc_cut=.5)
# mod.plot_gene_logfc(suff, attr='msex', show_labels=False, adj_txt=False, w=16, fc_cut=2)

# Plot correlation of SVD components with covariates:
cvlist = [
    'celltype',
    'region',
    'niareagansc',
    'cogdx',
    'cogdxad',
    'nft',
    'plaq_n',
    'plaq_d',
    'Apoe_e4',
    'msex',
    'pmi',
    'age_death',
    'n_genes',
    'n_counts']
mod.plot_svd_corr(cvlist)

# Plot the average expression of leiden modules on covariates:
cvlist = ['celltype', 'region', 'niareagansc', 'cogdx', 'Apoe_e4', 'msex']
mod.plot_heatmap_avgexpr(suff, cvlist=cvlist, attr='leiden')

mod.make_subset_graph(suff, 'msex', cutoff=0.5, plot=True)

mod.plot_graph(
    'vvt_msex',
    attr='leiden',
    show_labels=True,
    adj_txt=False,
    w=18)

# Test on AD / non-AD:
mod.plot_gene_logfc(
    'vvt',
    attr='nrad',
    show_labels=True,
    adj_txt=False,
    w=20,
    fc_cut=.5)


# Also run gprofiler:
gp = GProfiler(return_dataframe=True)
# Alternatively, use: GSEAPY with enrichr()

# NOTE: Can't run from internal...
mlist = mod.get_modules(suff, print_modules=False)
gpres = {}
for ll in mlist.keys():  # las:
    print(ll)
    testlist = mlist[ll].tolist()
    print(" ".join(testlist))
    gpres[ll] = gp.profile(organism='hsapiens', query=testlist)
    # lmat[:,ll] = np.mean(mod.adata[:,lnam].X, axis=1)

