#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_BM_400k_RNA/muon_env/bin/python3.8

import muon as mu
from muon import MuData
import numpy as np
import scanpy as sc
import pandas as pd
import os,sys
from scipy.io import mmwrite


import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# read intergrated RNA and ADT matrix
adata_rna = sc.read('to_muon_rna.h5ad')
adata_adt = sc.read('to_muon_adt.h5ad')

# filter cellbarcodes to the 50k post QC
df = pd.read_csv('Final-clean-titration.txt',sep='\t',index_col=0)
valid_barcodes = df.index.tolist()
adata_rna = adata_rna[valid_barcodes,:]
adata_adt = adata_adt[valid_barcodes,:]

# re-derive neighbors
sc.pp.neighbors(adata_rna,use_rep='X_pca_harmony')
sc.pp.neighbors(adata_adt,use_rep='X_pca')


# create mdata and do wnn
mdata = MuData({'rna':adata_rna,'adt':adata_adt})
mu.pp.neighbors(mdata)
mdata.write('mudata_subset.h5mu')

mdata = mu.read('mudata_subset.h5mu')
mu.tl.umap(mdata)
data = mdata.obsm['X_umap']
f = pd.DataFrame(data=data,index=mdata.obs_names,columns=['umap_x','umap_y'])
f.to_csv('new_umap.txt',sep='\t')


wnn_kernel = mdata.obsp['distances']
os.mkdir('./build_cite_ref_new')
mmwrite('./build_cite_ref_new/wnn_kernel.mtx',wnn_kernel) #this is the neighbor file used in Azimuth
mdata.obs_names.to_series().to_csv('./build_cite_ref_new/barcodes.tsv',sep='\t',header=None,index=None)

# create wnn at different resolutions and save
for r in [0.5,1,2,3,4,5,6]:
    sc.tl.leiden(mdata, resolution=r,key_added='leiden_wnn_{}'.format(r))
mdata.write('mudata_umap_leiden.h5mu')

mdata = mu.read('mudata_umap_leiden.h5mu')
for key in ['leiden_wnn_{}'.format(r) for r in [0.5,1,2,3,4,5,6]]:
    mu.pl.umap(mdata, color=[key],legend_loc='on data',legend_fontsize='xx-small')
    plt.savefig('muon_wnn_umap_{}.pdf'.format(key.replace(':','_')),bbox_inches='tight')
    plt.close()

#save wnn clustering results metadata
pd.DataFrame(data=mdata.obsm['X_umap'],index=mdata.obs_names,columns=['umap_x','umap_y']).to_csv('muon_wnn_umap.txt',sep='\t')
mdata.obs.to_csv('muon_wnn_metadata.txt',sep='\t')

#check modality weight
mu.pl.umap(mdata, color=['rna:mod_weight', 'adt:mod_weight'])
plt.savefig('muon_wnn_umap_weights.pdf',bbox_inches='tight')
plt.close()
for key in ['rna:ori_label','rna:sample','rna:donor']:
    mu.pl.umap(mdata, color=[key],legend_loc='on data',legend_fontsize='xx-small')
    plt.savefig('muon_wnn_umap_{}.pdf'.format(key.replace(':','_')),bbox_inches='tight')
    plt.close()





