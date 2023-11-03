#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import os,sys
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


'''
CD34-1/RNA/CD34-1/outs/SoupX/soupX-contamination-fraction-0.15  (66508, 32738)
CD34-2/RNA/CD34-2/outs/SoupX/soupX-contamination-fraction-0.15  (49252, 32738)
CD34-3/RNA/CD34-3/outs/SoupX/soupX-contamination-fraction-0.15  (53736, 32738)
CD34-4/RNA/CD34-4/outs/SoupX/soupX-contamination-fraction-0.15  (41607, 32738)

CD271-1/RNA/CD271-1/outs/SoupX/soupX-contamination-fraction-0.15  (33031, 32738)
CD271-2/RNA/CD271-2/outs/SoupX/soupX-contamination-fraction-0.15  (35063, 32738)
CD271-3/RNA/CD271-3/outs/SoupX/soupX-contamination-fraction-0.15  (45686, 32738)
CD271-4/RNA/CD271-4/outs/SoupX/soupX-contamination-fraction-0.15  (32132, 32738)

TNC-1-2-1/RNA/TNC-1-2-1/outs/SoupX/soupX-contamination-fraction-0.15   (8107, 32738)
TNC-1-2-2/RNA/TNC-1-2-2/outs/SoupX/soupX-contamination-fraction-0.15   (8509, 32738)
TNC-3-4-1/RNA/TNC-3-4-1/outs/SoupX/soupX-contamination-fraction-0.15   (10200, 32738)
TNC-3-4-2/RNA/TNC-3-4-2/outs/SoupX/soupX-contamination-fraction-0.15   (9917, 32738)
'''

# root_path = '/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/HBM_Titration_Xuan/220114_Grimes_GSL-PY-2582_Run1'
# adata_list = []
# for id_ in ['CD34-1','CD34-2','CD34-3','CD34-4','CD271-1','CD271-2','CD271-3','CD271-4','TNC-1-2-1','TNC-1-2-2','TNC-3-4-1','TNC-3-4-2']:
#     ind_path = '{}/RNA/{}/outs/SoupX/soupX-contamination-fraction-0.15'.format(id_,id_)
#     whole_path = os.path.join(root_path,ind_path)
#     print('reading {} from {}'.format(id_,whole_path))
#     adata = mtx_to_adata(int_folder=whole_path,gene_is_index=True,feature='genes.tsv',feature_col='index',barcode='barcodes.tsv',barcode_col='index',matrix='matrix.mtx')
#     adata.obs_names = ['.'.join([cell,id_]) for cell in adata.obs_names]
#     print(adata.shape)
#     adata_list.append(adata)

# common_gene = set(adata_list[0].var_names)
# for adata in adata_list[1:]:
#     genes = set(adata.var_names)
#     common_gene = common_gene.intersection(genes)
# common_gene = list(common_gene)
# print(len(common_gene)) . # 32738

# new_adata_list = []
# for adata in adata_list:
#     adata_new = adata[:,common_gene]
#     new_adata_list.append(adata_new)

# adata_combined_rna = ad.concat(adatas=new_adata_list,axis=0,label='sample',keys=['CD34-1','CD34-2','CD34-3','CD34-4','CD271-1','CD271-2','CD271-3','CD271-4','TNC-1-2-1','TNC-1-2-2','TNC-3-4-1','TNC-3-4-2'])
# print(adata_combined_rna)
# make_sure_adata_writable(adata_combined_rna)
# adata_combined_rna.write('adata_combined_rna.h5ad')  # 393748 × 32738, obs: sample


adata_combined_rna = sc.read('adata_combined_rna.h5ad')
# adata_combined_rna_subset = sc.pp.subsample(adata_combined_rna,n_obs=1000,copy=True)
# adata_combined_rna_subset.write('subset_rna_for_test.h5ad')

# add_annotations(adata_combined_rna,'Clusters/Azimuth-Marrow.txt',['predicted.celltype.l2'],0,['Azimuth_Bone_Marrow'],'\t','disk')
# add_annotations(adata_combined_rna,'Clusters/Haas.txt',['Predicted Class'],0,['Haas'],'\t','disk')
# add_annotations(adata_combined_rna,'Clusters/Haniffa.txt',['Predicted Class'],0,['Haniffa'],'\t','disk')
# add_annotations(adata_combined_rna,'Clusters/HCA.txt',['Predicted Class'],0,['HCA'],'\t','disk')
# add_annotations(adata_combined_rna,'Clusters/ICGS2-subclustering-R1.txt',['Predicted Class'],0,['ICGS2'],'\t','disk')
# # scanpy_recipe(adata_combined_rna,'human',False,[0.5,1,2,3,4,5,6],'rna',True,pca_n_comps=50,n_top_genes=3000)
# adata_combined_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
# umap_coord = pd.read_csv('UMAP/UMAP_scores.txt',sep='\t',index_col=0,header=None)
# umap_coord.columns = ['umap_x','umap_y']
# valid_cells = list(set(adata_combined_rna.obs_names).intersection(set(umap_coord.index)))
# adata_combined_rna = adata_combined_rna[valid_cells,:]   # 303199
# make_sure_adata_writable(adata_combined_rna)
# adata_combined_rna.write('../Hs_BM_400k_RNA/adata_combined_300k_rna.h5ad')

adata_combined_rna_post_qc = sc.read('../Hs_BM_400k_RNA/adata_combined_300k_rna.h5ad')
valid_cells = adata_combined_rna_post_qc.obs_names
adata_to_sean_mapping = adata_combined_rna[valid_cells,:]
print(adata_to_sean_mapping)
adata_to_sean_mapping.write('adata_to_sean_mapping.h5ad')
sys.exit('stop')


# add_umap(adata_combined_rna,umap_coord,'pandas_memory',['umap_x','umap_y'],0,'X_umap')
# # umap_dual_view_save(adata_combined_rna,cols=['Azimuth_Bone_Marrow','Haas','Haniffa','HCA','ICGS2'])
# # umap_dual_view_save(adata_combined_rna,cols=['sctri_rna_leiden_{}'.format(r) for r in [0.5,1,2,3,4,5,6]])

# just_log_norm(adata_combined_rna)
# # sctri = ScTriangulate(dir='output_two_5rna',adata=adata_combined_rna,query=['Azimuth_Bone_Marrow','Haas','Haniffa','HCA','ICGS2'])
# # sctri.lazy_run(compute_metrics_parallel=False,scale_sccaf=False)

# # sctri = ScTriangulate(dir='output_one_5rna',adata=adata_combined_rna,query=['Azimuth_Bone_Marrow','Haas','Haniffa','HCA','ICGS2'],add_metrics={})
# # sctri.lazy_run(compute_metrics_parallel=False,scale_sccaf=False)

# # sctri = ScTriangulate(dir='output_three_5rna',adata=adata_combined_rna,query=['Azimuth_Bone_Marrow','Haas','Haniffa','HCA','ICGS2'],
# #                       add_metrics={'tfidf5':tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
# # sctri.lazy_run(compute_metrics_parallel=False,scale_sccaf=False,added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}])


sctri = ScTriangulate.deserialize('output_three_5rna/after_pruned_assess.p')
# sctri.adata.obs.to_csv('output_three_5rna/barcode2cellmetadata.txt',sep='\t')
# pd.DataFrame(data=sctri.adata.obsm['X_umap'],index=sctri.adata.obs_names,columns=['umap_x','umap_y']).to_csv('output_three_5rna/umap_coords.txt',sep='\t')
# sctri.plot_umap(col='confidence',kind='continuous',format='png')
# sctri.plot_umap(col='final_annotation',kind='category',format='png')
# sctri.plot_umap(col='pruned',kind='category',format='png')

# how certain populations are being subdivided
# obs = sctri.adata.obs.copy()
'''1. HCA stromal population'''
key = 'HCA'
cluster = 'Stromal'
col = 'pruned'
series = obs.loc[obs[key]==cluster,:].groupby(col).apply(lambda x:x.shape[0])
subset = series.loc[series>10].index.tolist()
adata_s = sctri.plot_heterogeneity(key=key,cluster=cluster,style='build',col=col,subset=subset)
'''2. ERP, MKP, MEP population delineation'''
# dummy_col = []
# query_set = set(['HCA@CD34+_ERP','Haas@Early_MPPs','Haas@Metaphase_MPPs','HCA@CD34+_MKP'])
# for item in obs['pruned']:
#     if item in query_set:
#         dummy_col.append('MEP_MKP_early_ERP')
#     else:
#         dummy_col.append('others')
# obs['dummy'] = dummy_col
# sctri.adata.obs = obs
# sctri.plot_heterogeneity(key='dummy',cluster='MEP_MKP_early_ERP',style='build',col='pruned')
'''3. HSC and HSCP population'''
# dummy_col = []
# for item in obs['pruned']:
#     if 'HSC' in item:
#         dummy_col.append('HSC_HSCP')
#     else:
#         dummy_col.append('others')
# obs['dummy'] = dummy_col
# sctri.adata.obs = obs
# sctri.plot_heterogeneity(key='dummy',cluster='HSC_HSCP',style='build',col='pruned')

# # change umap coordiante
# umap_coord = pd.read_csv('UMAP/UMAP_scores-revised.txt',sep='\t',index_col=0,header=None)
# umap_coord.columns = ['umap_x','umap_y']
# common_barcodes = list(set(sctri.adata.obs_names).intersection(set(umap_coord.index)))
# sctri.adata = sctri.adata[common_barcodes,:]  # 303195 × 32738, 4 cells being removed
# add_umap(sctri.adata,umap_coord,'pandas_memory',['umap_x','umap_y'],0,'X_umap')
# sctri.plot_umap(col='confidence',kind='continuous',format='pdf',umap_cmap='viridis')
# sctri.plot_umap(col='final_annotation',kind='category',format='pdf')
# sctri.plot_umap(col='pruned',kind='category',format='pdf')





# # assess antibody concentration
# group = pd.read_csv('groups.hybrid.txt',sep='\t',index_col=0,header=None)
# concentrations = [float(concentration) for _,concentration in group[1].str.split('--')]
# group[3] = concentrations

# '''
# 0.25    56356
# 0.50    56160
# 1.00    55387
# 4.00    55225
# 2.00    54442
# 1.50    38184
# '''

# # adata_combined_adt = small_txt_to_adata(int_file='TotalVI-ADTs/denoised_adt_values_2021_03_09_400_epochs-transposed.txt',gene_is_index=True)
# # make_sure_adata_writable(adata_combined_adt)
# # adata_combined_adt.write('adata_combined_adt.h5ad')  # 315754 × 275

# adata_combined_adt = sc.read('adata_combined_adt.h5ad')
# umap_coord = pd.read_csv('UMAP/UMAP_scores.txt',sep='\t',index_col=0,header=None)
# umap_coord.columns = ['umap_x','umap_y']
# valid_cells = list(set(adata_combined_adt.obs_names).intersection(set(umap_coord.index)))
# adata_combined_adt = adata_combined_adt[valid_cells,:]   
# add_umap(adata_combined_adt,umap_coord,'pandas_memory',['umap_x','umap_y'],0,'X_umap')
# add_annotations(adata_combined_adt,group,[3],0,['concentration'],'\t','memory')
# add_annotations(adata_combined_adt,'output_two_5rna/barcode2cellmetadata.txt',['pruned'],0,['rna_two_tfidf_pruned'],'\t','disk')

# for c in [0.25,0.5,1.00,2.00,4.00,1.50]:  # 54098 × 275, 53877 × 275, 53191 × 275, 52334 × 275, 53151 × 275, 36548 × 275
#     adata_subset = adata_combined_adt[adata_combined_adt.obs['concentration']==str(c),:]
#     print(adata_subset)
#     adata_subset.write('adata_adt_{}.h5ad'.format(c))
#     adata_subset.obs.rename(columns={'rna_two_tfidf_pruned':'rna_two_tfidf_pruned_adt_{}'.format(c)},inplace=True)
#     umap_dual_view_save(adata_subset,cols=['rna_two_tfidf_pruned_adt_{}'.format(c)])


# c = 0.5
# adata_subset = sc.read('adata_adt_{}.h5ad'.format(c))
# adata_subset.X = Normalization.GMM_normalization(make_sure_mat_dense(adata_subset.X),True)
# sctri_subset = ScTriangulate(dir='output_adt_GMM_non_negative_{}'.format(c),adata=adata_subset,query=['rna_two_tfidf_pruned'],add_metrics={})
# sctri_subset.run_single_key_assessment(key='rna_two_tfidf_pruned',scale_sccaf=False,layer=None,added_metrics_kwargs={})
# sctri_subset.serialize('sctri_subset.p')

# for c in [0.25,0.5,1.00,2.00,4.00,1.50]:
#     root_dir = 'output_adt_GMM_non_negative_{}'.format(c)
#     pickle_file = 'sctri_subset.p'
#     full_path = os.path.join(root_dir,pickle_file)
#     sctri = ScTriangulate.deserialize(full_path)
#     sctri.gene_to_df(mode='exclusive_genes',key='rna_two_tfidf_pruned',n=275)
#     sctri.gene_to_df(mode='marker_genes',key='rna_two_tfidf_pruned',col='purify')
#     sctri.plot_confusion(name='confusion_reassign',key='rna_two_tfidf_pruned',labelsize=1,xticklabels=1,yticklabels=1)
#     sctri.plot_confusion(name='confusion_sccaf',key='rna_two_tfidf_pruned',labelsize=1,xticklabels=1,yticklabels=1)
#     sctri.viewer_cluster_feature_html()
#     sctri.viewer_cluster_feature_figure(select_keys=['rna_two_tfidf_pruned'])

# for c in [0.25,0.5,1.00,2.00,4.00,1.50]:
#     root_dir = 'output_adt_GMM_non_negative_{}'.format(c)
#     pickle_file = 'sctri_subset.p'
#     full_path = os.path.join(root_dir,pickle_file)
#     sctri = ScTriangulate.deserialize(full_path)
#     sctri.adata.to_df().to_csv(os.path.join(root_dir,'scaled_matrix.txt'),sep='\t')



# c = 2.0
# adata = sc.read('adata_adt_{}.h5ad'.format(c))
# adata.X = Normalization.GMM_normalization(make_sure_mat_dense(adata.X),True)
# ncols = 5
# nrows = adata.shape[1] // ncols + 1
# fig, axes = plt.subplots(nrows=nrows,ncols=ncols,gridspec_kw={'wspace':0.8,'hspace':0.5},figsize=(20,80))
# axes = axes.flatten()
# for i,var in tqdm(enumerate(adata.var_names),total=adata.shape[1]):
#     sns.histplot(data=make_sure_mat_dense(adata[:,var].X),bins=100,ax=axes[i])
#     axes[i].set_title(var,fontsize=3)
# plt.savefig('adt_{}_histplot_GMM_non_negative.pdf'.format(c),bbox_inches='tight')
# plt.close()











