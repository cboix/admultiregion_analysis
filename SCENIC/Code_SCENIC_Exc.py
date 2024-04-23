cd /mnt/d/Brain_region_project/SCENIC/
conda activate scenic_protocol
python

import scanpy as sc
import numpy as np
import pandas as pd
adata=sc.read_csv("Exc_counts_matrix_1000.csv")
sc.pp.filter_genes(adata,min_cells=100)
adata.raw=adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=15000)
adata_raw=sc.AnnData(X=adata.raw.X,obs=adata.obs,var=adata.var)
adata_raw=adata_raw[:,adata_raw.var.highly_variable]
del adata_raw.var['highly_variable']


import loompy as lp
row_attrs = {
"Gene": np.array(adata_raw.var_names),
}
col_attrs = {
"CellID": np.array(adata_raw.obs_names),
"nGene": np.array(np.sum(adata_raw.X.transpose()>0, axis=0)).flatten(),
"nUMI": np.array(np.sum(adata_raw.X.transpose(),axis=0)).flatten(),
}
lp.create("Exc_filtered.loom",adata_raw.X.transpose(),row_attrs,
col_attrs)

##########################################################################

cd /mnt/d/Brain_region_project/SCENIC/
conda activate scenic_protocol

pyscenic grn \
--num_workers 20 \
--output Exc.tsv \
--method grnboost2 \
Exc_filtered.loom \
hs_hgnc_tfs.txt


pyscenic ctx \
Exc.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname Exc_filtered.loom \
--mode "dask_multiprocessing" \
--output Exc_reg.csv \
--num_workers 20 \
--mask_dropouts


pyscenic aucell \
Exc_filtered.loom \
Exc_reg.csv \
--output Exc_SCENIC.loom \
--num_workers 20

##########################################################################

conda activate scenic_protocol
cd /mnt/d/Brain_region_project/SCENIC
python
import scanpy as sc

adata=sc.read_csv("Exc_counts_matrix_1000.csv")
adata.obs.to_csv("Exc_counts_matrix_1000_metadata.csv")

import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
adata.raw = adata
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

tsne = TSNE( n_jobs=20 )
adata.obsm['X_tsne'] = tsne.fit_transform( adata.X )
sc.tl.louvain(adata,resolution=0.375)
import json
import zlib
import base64

lf = lp.connect( "Exc_SCENIC.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
import umap

runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "Exc_scenic_umap.txt", sep='\t')

tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "Exc_scenic_tsne.txt", sep='\t')

lf = lp.connect("Exc_SCENIC.loom", mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( './Exc_scenic_umap.txt', sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv( './Exc_scenic_tsne.txt', sep='\t', header=0, index_col=0 )
auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
rt = meta['regulonThresholds']

for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )

for i,x in enumerate( rt ):
        x['motifData'] = 'NA.png'

tsneDF = pd.DataFrame(adata.obsm['X_umap'], columns=['_X', '_Y'])

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
                pd.DataFrame(adata.obsm['X_tsne'],index=adata.obs.index)[0] ,
                pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
                dr_tsne['X'] ,
                dr_umap['X']
            ], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
                pd.DataFrame(adata.obsm['X_tsne'],index=adata.obs.index)[1] ,
                pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
                dr_tsne['Y'] ,
                dr_umap['Y']
            ], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3','4']

garnett = pd.read_csv( 'Exc_cell_type_1000.txt', sep='\t', index_col=0)

metaJson = {}

metaJson['embeddings'] = [
    {
        "id": -1,
        "name": "highly variable genes UMAP"
    },
    {
        "id": 1,
        "name": "highly variable genes t-SNE"
    },
    {
        "id": 2,
        "name": "Scanpy PC1/PC2"
    },
    {
        "id": 3,
        "name": "SCENIC AUC t-SNE"
    },
    {
        "id": 4,
        "name": "SCENIC AUC UMAP"
    },
]

metaJson["clusterings"] = [{
            "id": 0,
            "group": "Scanpy",
            "name": "Scanpy louvain default resolution",
            "clusters": [],
        }]

metaJson["metrics"] = [
        {
            "name": "nUMI"
        }, {
            "name": "nGene"
        }, {
            "name": "Percent_mito"
        }
]

metaJson["annotations"] = [
    {
        "name": "Louvain_clusters_Scanpy",
        "values": list(set( adata.obs['louvain'].astype(np.str) ))
    },
    {
        "name": "Celltype",
        "values": list(set( garnett['cell_type_high_resolution'].astype(np.str) ))
    },
]


metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
    clustDict = {}
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i}'
    metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()
clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)


def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr

col_attrs = {
        "CellID": np.array(adata.obs.index),

        "Celltype": np.array( garnett['cell_type_high_resolution'].values ),

        "Embedding": dfToNamedMatrix(tsneDF),
        "Embeddings_X": dfToNamedMatrix(Embeddings_X),
        "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
        "RegulonsAUC": dfToNamedMatrix(auc_mtx),
        "Clusterings": dfToNamedMatrix(clusterings),
        "ClusterID": np.array(adata.obs['louvain'].values)
}

row_attrs = {
    "Gene": lf.ra.Gene,
    "Regulons": regulons,
}

attrs = {
    "title": "sampleTitle",
    "MetaData": json.dumps(metaJson),
    "Genome": 'hg38',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

f_final_loom = 'Exc_tracks.loom'
lp.create(
    filename = f_final_loom ,
    layers=lf[:,:],
    row_attrs=row_attrs,
    col_attrs=col_attrs,
    file_attrs=attrs
)
lf.close()

##########################################################################

conda activate scenic_protocol
cd /mnt/d/Brain_region_project/SCENIC
python

import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
from scanpy.plotting._tools.scatterplots import plot_scatter
import seaborn as sns

f_final_loom = 'Exc_tracks.loom'

# scenic output
lf = lp.connect( f_final_loom, mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
        [
            pd.DataFrame( lf.ca.Celltype, index=lf.ca.CellID )
        ],
        axis=1
    )
cellAnnot.columns = [
     'Celltype'
     ]
lf.close()

from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

rss_cellType = regulon_specificity_scores( auc_mtx, cellAnnot['Celltype'] )

auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

#Export results

auc_mtx_Z.to_csv('Exc_AUC_mtx_z_1000.csv',index=True)

rss_cellType.to_csv('Exc_rss_Celltype_1000.csv', index = True)

(pd.DataFrame.from_dict(data=regulons, orient='index').to_csv('Exc_regulons_1000.csv', header=False))

auc_mtx.to_csv('Exc_AUC_mtx_1000.csv',index=True)

cellAnnot.to_csv('Exc_cellAnnot_1000.csv',index=True)
