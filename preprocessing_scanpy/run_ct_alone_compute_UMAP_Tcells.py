#!usr/bin/python
# --------------------------------------------------------------
# Compute a UMAP for each major celltype alone for visualization
# Updated: 09/03/21
# --------------------------------------------------------------
import sys
bindir = '/home/cboix/data/DEVTRAJ/bin/'
sys.path.append(bindir + 'multiRegion/pseudotime')
from importlib import reload
from adata_loader import *

# Load in all of the relevant, cell type-specific data:
# -----------------------------------------------------
ct = 'Mic_Immune'
st = 'T cells'
scdata = scDataLoader(ct, st,
                      usecombat=True, usehvg=False,
                      remove_mt=False, filter_TH=False)
scdata.load()

# Compute UMAP from PCA with k=50 and NN=15:
# ------------------------------------------
NN = 15
KSVD = 50
if st == 'T cells':
    NN = 8

sc.tl.pca(scdata.adata, n_comps=KSVD)
sc.pp.neighbors(scdata.adata, n_neighbors=NN)
sc.tl.umap(scdata.adata, maxiter=None, random_state=0)

# Plot the UMAP for the major celltype:
# -------------------------------------
prefstr = '_majorctUMAP_' + scdata.csuff
sc.pl.umap(scdata.adata, color='celltype', frameon=False, save=prefstr + '_ctalone.png')
sc.pl.umap(scdata.adata, color='projid', frameon=False, save=prefstr + '_ctalone_projid.png')
sc.pl.umap(scdata.adata, color='region', frameon=False, save=prefstr + '_ctalone_region.png')

# Write UMAP coordinates to file:
# -------------------------------
scdata.adata.obs['C1'] = scdata.adata.obsm['X_umap'][:,0]
scdata.adata.obs['C2'] = scdata.adata.obsm['X_umap'][:,1]
adf = pd.DataFrame(scdata.adata.obs)[['barcode','C1','C2']]
adf.to_csv('multiregion_majorctUMAP_' + scdata.csuff + '.tsv', sep="\t")

# Leiden + new markers:
# ---------------------
sc.tl.leiden(scdata.adata, resolution=2)
sc.pl.umap(scdata.adata, color='leiden', frameon=False,
           save=prefstr + '_ctalone_leiden.png')
sc.pl.umap(scdata.adata, color='leiden', frameon=False, legend_loc='on data',
           save=prefstr + '_ctalone_leiden2.png')

NG = 25
sc.tl.rank_genes_groups(scdata.adata, 'leiden')
sc.pl.rank_genes_groups(scdata.adata, n_genes=NG, sharey=False,
                        frameon=False, save=prefstr + '_ctalone_leiden.png')
sc.pl.rank_genes_groups_dotplot(scdata.adata, n_genes=5,
                                save=prefstr + '_ctalone_leiden.png')
sc.pl.rank_genes_groups_matrixplot(scdata.adata, n_genes=5,
                                   save=prefstr + '_ctalone_leiden.png')


# Annotate T-cells:
# -----------------
tdf = pd.read_csv('../Annotation/cdt_annotation.tsv', sep="\t", header=None)
tdf.columns = ['set','gene']

sets = np.unique(tdf.set)
mjr_dict = {}
for subset in sets:
    mjr_dict[subset] = tdf['gene'][tdf.set == subset].tolist()

## set up sets for heatmaps:
def dict_to_heatmapattr(mjr_dict, adata):
    i = 0; tplist = []; lblist = []; pltmarkers = []
    for key in mjr_dict.keys():
        gl = [name for name in mjr_dict[key] if name in adata.var_names]
        n = len(gl)
        if (n > 0):
            pltmarkers += gl
            lblist.append(key)
            tplist.append((i,i+n-1))
            i += n
    return(tplist, lblist, pltmarkers)

tplist, lblist, pltmarkers = dict_to_heatmapattr(mjr_dict, scdata.adata)

# Make Heatmap plot (plain + reord)
# for group in ['auto.major','auto.cthr']:
group = 'leiden'
# Make heatmap plot:
fig = sc.pl.matrixplot(scdata.adata, pltmarkers, groupby=group,
                       var_group_rotation=90, var_group_positions=tplist,
                       var_group_labels=lblist, log=True, cmap='Blues',
                       save=prefstr + "_" + group + '_markers.png')

fig = sc.pl.heatmap(scdata.adata, pltmarkers, groupby=group,
                    var_group_rotation=90, var_group_positions=tplist,
                    var_group_labels=lblist, log=True, cmap='Blues',
                    save=prefstr + "_" + group + '_markers.png')

genes = ['STAT1','IL4R','TCF7','THEMIS','GZMB','CCL4','CCL3',
         'CCL5','NKG7','ITGA1','SKAP1', 'MAF']
sc.pl.umap(scdata.adata, color=genes, frameon=False,
           save=prefstr + '_ctalone_markers.png')

genes = ['CD8A','CD8B','CD4','CD96','CD69','CD247','CD3E',
         'CD68','HBB','MS4A1', 'CD74', 'CD163']
sc.pl.umap(scdata.adata, color=genes, frameon=False,
           save=prefstr + '_ctalone_markers2.png')

fig = sc.pl.matrixplot(scdata.adata, genes, groupby=group,
                       var_group_rotation=90, # var_group_positions=tplist,
                       # var_group_labels=lblist,
                       log=True, cmap='Blues',
                       save=prefstr + "_" + group + '_markers2.png')

genes = ['CD69','CPA3','TPSB2', 'IL7R','RIPOR2','PSAP','ITGA1','ITGA4',
         'IL4R','STAT1','TNF']
sc.pl.umap(scdata.adata, color=genes, frameon=False,
           save=prefstr + '_ctalone_markers3.png')

# For T cells, re-annotate:
# -------------------------
scdata.adata.obs['subtype'] = ''
tmap = {'CD8': ['5','9','7','12','1','18','0','10','16'],
        'CD4': ['8'],
        'NKT': ['3','4'],
        'CD8_act': ['2'],
        'B cells': ['13'],
        'CD69+': ['14'],
        'T': ['15'],
        'PSAP+_1': ['11'],
        'PSAP+_2': ['19'],
        'Batch1': ['6'],
        'Batch2': ['17'],
        'Cancer_Unkw': ['20']}

for subset in tmap.keys():
    ind = np.where(np.isin(scdata.adata.obs['leiden'], tmap[subset]))[0]
    scdata.adata.obs['subtype'][ind] = subset


sc.pl.umap(scdata.adata, color='subtype', frameon=False, legend_loc='on data',
           save=prefstr + '_ctalone_subtype.png')

df = scdata.adata.obs[['region','projid','subtype', 'leiden','niareagansc']]



import matplotlib as mpl
from matplotlib import pyplot as plt

aggdf = pd.crosstab(scdata.adata.obs['subtype'],
                    scdata.adata.obs['region'])
normdf = aggdf / aggdf.sum(1)[:,None]

res = normdf.plot.barh(stacked=True, figsize=(8,3), width=.8)
fname = 'figures/barplot' + prefstr + '_regsubtype.png'
plt.tight_layout()
plt.savefig(fname, dpi=350, bbox_inches='tight')
plt.close()


aggdf = pd.crosstab(scdata.adata.obs['subtype'],
                    scdata.adata.obs['projid'])
normdf = aggdf / aggdf.sum(1)[:,None]

res = normdf.plot.barh(stacked=True, figsize=(8,3), width=.8, legend=False)
fname = 'figures/barplot' + prefstr + '_pidsubtype.png'
plt.tight_layout()
plt.savefig(fname, dpi=350, bbox_inches='tight')
plt.close()


aggdf = pd.crosstab(scdata.adata.obs['subtype'],
                    scdata.adata.obs['niareagansc'])
normdf = aggdf / aggdf.sum(1)[:,None]

res = normdf.plot.barh(stacked=True, figsize=(8,3), width=.8)
fname = 'figures/barplot' + prefstr + '_nrsubtype.png'
plt.tight_layout()
plt.savefig(fname, dpi=350, bbox_inches='tight')
plt.close()










