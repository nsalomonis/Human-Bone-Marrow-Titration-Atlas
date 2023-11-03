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






# large_txt_to_mtx(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/Merged_validation4/Merged_4validation.txt',
#                  out_folder='./100k_data',
#                  gene_is_index=True,
#                  type_convert_to='float32',
#                  n_lines=32739,
#                  sep='\t')
# adata_combined_100k_rna = mtx_to_adata(int_folder='./100k_data',
#                                        gene_is_index=True,
#                                        feature='genes.tsv',
#                                        feature_col='index',
#                                        barcode='barcodes.tsv',
#                                        barcode_col='index',
#                                        matrix='matrix.mtx')
# adata_combined_100k_rna.obs_names.name = None
# adata_combined_100k_rna.var_names.name = None
# adata_combined_100k_rna.write('adata_combined_100k_rna.h5ad')  # 116341 × 32738

# adata_combined_300k_rna = sc.read('adata_combined_300k_rna.h5ad') # 303199 × 32738
# adata_combined_300k_rna = ad.AnnData(X=adata_combined_300k_rna.X,
#                                      obs=adata_combined_300k_rna.obs['sample'].to_frame(),
#                                      var=pd.DataFrame(index=adata_combined_300k_rna.var_names))
# sc.pp.normalize_total(adata_combined_300k_rna,target_sum=1e4)
# sc.pp.log1p(adata_combined_300k_rna,base=2)  # log2 instead of natural logarithm
# anno = pd.read_csv('QueryGroups.cellHarmony.300k.txt',sep='\t',index_col=0,header=None)
# anno.columns = ['label','repeat']
# add_annotations(adata_combined_300k_rna,inputs=anno,cols_input=['label'],index_col=0,cols_output=['ori_label'],sep='\t',kind='memory')

adata_combined_100k_rna = sc.read('adata_combined_100k_rna.h5ad')
anno = pd.read_csv('QueryGroups.cellHarmony.100k.txt',sep='\t',index_col=0,header=None)
anno.columns = ['label','repeat']
adata_combined_100k_rna = adata_combined_100k_rna[anno.index,:]   # 104660 × 32738
denoised_adt = pd.read_csv('/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/total_VI_2/totalvi_denoised_adt_values_RC2_15.txt',sep='\t',index_col=0) # 97162
# convert from BM27_CD34_AAACCCAAGCGTCAGA-1-0 to AAACCCAAGATTCGCT-1.BF21_032123_CD34
adata_combined_100k_rna.obs['sample'] = [item.split('.')[1] for item in adata_combined_100k_rna.obs_names]
vc = adata_combined_100k_rna.obs['sample'].value_counts()
dic = {'_'.join([item.split('_')[0],item.split('_')[2]]):item for item in vc.index}  # {WF26_CD34: WF26_031423_CD34}
col = []
for item in denoised_adt.index:
    sample_name = '_'.join(item.split('_')[:2])
    barcode = '-'.join(item.split('_')[2].split('-')[:-1])
    new = barcode + '.' + dic[sample_name]
    col.append(new)
denoised_adt.index = col
# continue subset
common = list(set(adata_combined_100k_rna.obs_names).intersection(set(denoised_adt.index)))  # 91910
adata_combined_100k_rna = adata_combined_100k_rna[common,:]  # 91910 × 32738
# inspect
adata_combined_100k_rna.obs['sample'] = [item.split('.')[1] for item in adata_combined_100k_rna.obs_names]
vc = adata_combined_100k_rna.obs['sample'].value_counts()

# '''
# original:
# WF26_031423_CD34     11515
# WF83_021523_CD271     9885
# BF21_032123_CD271     9548
# BM27_120522_CD271     8580
# WF83_021523_TNC       8545
# WF26_031423_CD271     8395
# BM27_120522_CD34      8129
# WM34_120522_CD34      7974
# BM27_120522_TNC       7656
# WF83_021523_CD34      7022
# WM34_120522_CD271     6718
# BF21_032123_CD34      6537
# WM34_120522_TNC       5989
# BF21_032123_TNC       5916
# WF26_031423_TNC       3932

# now:
# WF26_031423_CD34     9546
# WF83_021523_CD271    7903
# BF21_032123_CD271    7789
# BM27_120522_CD271    7340
# BM27_120522_CD34     6989
# WF26_031423_CD271    6875
# WM34_120522_CD34     6796
# WF83_021523_TNC      6083
# WF83_021523_CD34     5745
# BF21_032123_CD34     5390
# WM34_120522_CD271    5356
# BM27_120522_TNC      5329
# BF21_032123_TNC      4333
# WM34_120522_TNC      3958
# WF26_031423_TNC      2478
# '''

# # add annotation
# add_annotations(adata_combined_100k_rna,inputs=anno,cols_input=['label'],index_col=0,cols_output=['ori_label'],sep='\t',kind='memory')
# adata_combined_400k_rna = ad.concat([adata_combined_300k_rna,adata_combined_100k_rna],axis=0,join='inner',merge='first',label='batch',keys=['old','new'])
# adata_combined_400k_rna.write('adata_combined_400k_rna.h5ad')  # 395109 × 32738


# adata_combined_400k_rna = sc.read('adata_combined_400k_rna.h5ad')
# adata_combined_400k_rna = scanpy_recipe(adata_combined_400k_rna,is_log=True,resolutions=[0.5,1,2,3,4,5,6],modality='rna',umap=True,save=True,pca_n_comps=50,n_top_genes=3000)

# adata_combined_400k_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
# umap_dual_view_save(adata_combined_400k_rna,cols=['sctri_rna_leiden_{}'.format(str(r)) for r in [0.5,1,2,3,4,5,6]]+['sample','ori_label','batch'])


'''test if there are batch effect on denoised ADT data'''
# build adt anndata
adata_combined_100k_adt = ad.AnnData(X=denoised_adt.values,var=pd.DataFrame(index=denoised_adt.columns),obs=pd.DataFrame(index=denoised_adt.index))
adata_combined_100k_adt.obs['sample'] = [item.split('.')[1] for item in adata_combined_100k_adt.obs_names]
anno = pd.read_csv('QueryGroups.cellHarmony.100k.txt',sep='\t',index_col=0,header=None)
anno.columns = ['label','repeat']
add_annotations(adata_combined_100k_adt,inputs=anno,cols_input=['label'],index_col=0,cols_output=['ori_label'],sep='\t',kind='memory')
col = []
for row in adata_combined_100k_adt.obs.itertuples():
    if 'WF26' in row.sample:
        col.append('WF26')
    elif 'WF83' in row.sample:
        col.append('WF83')
    elif 'BF21' in row.sample:
        col.append('BF21')
    elif 'BM27' in row.sample:
        col.append('BM27')
    else:
        col.append('WM34')
adata_combined_100k_adt.obs['donor'] = col

'''below we remove the WF83'''
adata_combined_100k_adt = adata_combined_100k_adt[adata_combined_100k_adt.obs['donor']!='WF83',:] # 97162 × 142 to 75979 × 142

# adata_combined_100k_adt = scanpy_recipe(adata_combined_100k_adt,is_log=False,resolutions=[0.5,1,2,3,4,5,6],modality='adt',pca_n_comps=15)
# umap_dual_view_save(adata_combined_100k_adt,cols=['sctri_adt_leiden_{}'.format(str(r)) for r in [0.5,1,2,3,4,5,6]]+['sample','ori_label','donor'])


'''start running scTriangulate'''

adata_combined_400k_rna_harmony = sc.read('adata_combined_400k_rna_harmony.h5ad')
adata_combined_100k_rna_harmony = adata_combined_400k_rna_harmony[adata_combined_400k_rna_harmony.obs['batch']=='new',:]  # 91910 × 32738

'''below we remove the WF83'''
adata_combined_100k_rna_harmony = adata_combined_100k_rna_harmony[adata_combined_100k_rna_harmony.obs['donor']!='WF83',:]  # 72179 × 32738

# umap_dual_view_save(adata_combined_100k_rna_harmony,cols=['sample','ori_label','donor'])

adata_combined_100k_adt = sc.read('adata_after_scanpy_recipe_adt_0.5_1_2_3_4_5_6_umap_True_new.h5ad')  
adata_combined_100k_adt = adata_combined_100k_adt[adata_combined_100k_rna_harmony.obs_names,:] # 91910 × 142   # 72179 × 142

adata_final = concat_rna_and_other(adata_combined_100k_rna_harmony,adata_combined_100k_adt,umap='rna',umap_key='X_umap',name='adt',prefix='AB_')
del adata_final.obsm['X_pca']
del adata_final.obsm['X_pca_harmony']
# umap_dual_view_save(adata_final,cols=['ori_label','sample','donor']+['sctri_adt_leiden_{}'.format(str(r)) for r in [0.5,1,2,3,4,5,6]])

## build the raw layer
root_dir = '/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/ADT'

BF21_CD34 = small_txt_to_adata(int_file=os.path.join(root_dir,'BF21-CD34_ADT_clean.txt'),gene_is_index=True,sep='\t')
BF21_CD34.obs_names = [item + '.BF21_032123_CD34' for item in BF21_CD34.obs_names]
BF21_CD271 = small_txt_to_adata(int_file=os.path.join(root_dir,'BF21-CD271_ADT_clean.txt'),gene_is_index=True,sep='\t')
BF21_CD271.obs_names = [item + '.BF21_032123_CD271' for item in BF21_CD271.obs_names]
BF21_TNC = small_txt_to_adata(int_file=os.path.join(root_dir,'BF21-TNC_ADT_clean.txt'),gene_is_index=True,sep='\t')
BF21_TNC.obs_names = [item + '.BF21_032123_TNC' for item in BF21_TNC.obs_names]

BM27_CD34 = small_txt_to_adata(int_file=os.path.join(root_dir,'BM27-CD34_ADT_clean.txt'),gene_is_index=True,sep='\t')
BM27_CD34.obs_names = [item + '.BM27_120522_CD34' for item in BM27_CD34.obs_names]
BM27_CD271 = small_txt_to_adata(int_file=os.path.join(root_dir,'BM27_CD271_ADT_clean.txt'),gene_is_index=True,sep='\t')
BM27_CD271.obs_names = [item + '.BM27_120522_CD271' for item in BM27_CD271.obs_names]
BM27_TNC = small_txt_to_adata(int_file=os.path.join(root_dir,'BM27-TNC_ADT_clean.txt'),gene_is_index=True,sep='\t')
BM27_TNC.obs_names = [item + '.BM27_120522_TNC' for item in BM27_TNC.obs_names]

WF26_CD34 = small_txt_to_adata(int_file=os.path.join(root_dir,'WF26_CD34_ADT_clean.txt'),gene_is_index=True,sep='\t')
WF26_CD34.obs_names = [item + '.WF26_031423_CD34' for item in WF26_CD34.obs_names]
WF26_CD271 = small_txt_to_adata(int_file=os.path.join(root_dir,'WF26_CD271_ADT_clean.txt'),gene_is_index=True,sep='\t')
WF26_CD271.obs_names = [item + '.WF26_031423_CD271' for item in WF26_CD271.obs_names]
WF26_TNC = small_txt_to_adata(int_file=os.path.join(root_dir,'WF26_TNC_ADT_clean.txt'),gene_is_index=True,sep='\t')
WF26_TNC.obs_names = [item + '.WF26_031423_TNC' for item in WF26_TNC.obs_names]

WF83_CD34 = small_txt_to_adata(int_file=os.path.join(root_dir,'WF83-CD34-hg19_ADT_clean.txt'),gene_is_index=True,sep='\t')
WF83_CD34.obs_names = [item + '.WF83_021523_CD34' for item in WF83_CD34.obs_names]
WF83_CD271 = small_txt_to_adata(int_file=os.path.join(root_dir,'WF83-CD271-hg19_ADT_clean.txt'),gene_is_index=True,sep='\t')
WF83_CD271.obs_names = [item + '.WF83_021523_CD271' for item in WF83_CD271.obs_names]
WF83_TNC = small_txt_to_adata(int_file=os.path.join(root_dir,'WF83-TNC-hg19_ADT_clean.txt'),gene_is_index=True,sep='\t')
WF83_TNC.obs_names = [item + '.WF83_021523_TNC' for item in WF83_TNC.obs_names]

WM34_CD34 = small_txt_to_adata(int_file=os.path.join(root_dir,'WM34-CD34_ADT_clean.txt'),gene_is_index=True,sep='\t')
WM34_CD34.obs_names = [item + '.WM34_120522_CD34' for item in WM34_CD34.obs_names]
WM34_CD271 = small_txt_to_adata(int_file=os.path.join(root_dir,'WM34-CD271_ADT_clean.txt'),gene_is_index=True,sep='\t')
WM34_CD271.obs_names = [item + '.WM34_120522_CD271' for item in WM34_CD271.obs_names]
WM34_TNC = small_txt_to_adata(int_file=os.path.join(root_dir,'WM34_TNC_ADT_clean.txt'),gene_is_index=True,sep='\t')
WM34_TNC.obs_names = [item + '.WM34_120522_TNC' for item in WM34_TNC.obs_names]


'''make sure to exclude WF83'''

adata_raw_adt = ad.concat(adatas=[BF21_CD34,BF21_CD271,BF21_TNC,BM27_CD34,BM27_CD271,BM27_TNC,WF26_CD34,WF26_CD271,WF26_TNC,WM34_CD34,WM34_CD271,WM34_TNC],
                          axis=0,join='inner',merge='first',label='sample',
                          keys=['BF21_CD34','BF21_CD271','BF21_TNC','BM27_CD34','BM27_CD271','BM27_TNC','WF26_CD34','WF26_CD271','WF26_TNC','WM34_CD34','WM34_CD271','WM34_TNC'])
adata_raw_adt = adata_raw_adt[adata_final.obs_names,:]  # 91910 × 142 to 72179 × 142
adata_raw_rna = adata_combined_100k_rna_harmony
adata_raw_combined = concat_rna_and_other(adata_raw_rna,adata_raw_adt,umap='rna',umap_key='X_umap',name='adt',prefix='AB_')
adata_raw_combined = adata_raw_combined[adata_final.obs_names,adata_final.var_names]

adata_final.layers['raw'] = adata_raw_combined.X



'''create two anndata for muon'''
# adata_combined_100k_rna_harmony.write('to_muon_rna.h5ad')
# adata_combined_100k_adt.write('to_muon_adt.h5ad')

# now wnn and clusters are ready, do a dry run of scTriangulate
# replace X_umap and wnn_leiden clusters
add_umap(adata_final,inputs='muon_wnn_umap.txt',mode='pandas_disk',cols=['umap_x','umap_y'],index_col=0,key='X_umap')
add_annotations(adata_final,inputs='muon_wnn_metadata.txt',
                cols_input=['rna:mod_weight','adt:mod_weight'] + ['leiden_wnn_{}'.format(r) for r in [0.5,1,2,3,4,5,6]],
                cols_output=['rna_mod_weight','adt_mod_weight'] + ['leiden_wnn_{}'.format(r) for r in [0.5,1,2,3,4,5,6]])


'''Nathan asked to run scTriangulate without ADT'''
adata_final_minus_adt = adata_final.copy()
del adata_final_minus_adt.layers['raw']
adata_final_minus_adt = adata_final_minus_adt[adata_final.obs_names,adata_combined_100k_rna_harmony.var_names]  # 72179 × 32738

# here, we want to re-run and replace ori_label with pruned cluster 
anno = pd.read_csv('groups.TFIDF3-WNN3.txt',sep='\t',index_col=0,header=None)
anno.columns = ['label','repeat']
add_annotations(adata_final_minus_adt,inputs=anno,cols_input=['label'],cols_output=['ori_label'],kind='memory')

sctri = ScTriangulate(dir='test_run_output_tfidf3_ori_1_2_3_minus_adt_Frank',
                      adata=adata_final_minus_adt,
                      query=['ori_label'] + ['leiden_wnn_{}'.format(r) for r in [1,2,3]],
                      add_metrics={'tfidf5': tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
sctri.lazy_run(layer=None,compute_metrics_parallel=False,
               added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}],
               assess_pruned=True,viewer_cluster=True,viewer_heterogeneity=True,viewer_heterogeneity_keys=['ori_label'])
sys.exit('stop')


sctri = ScTriangulate(dir='test_run_output_tfidf3_ori_1_2_3',
                      adata=adata_final,
                      query=['ori_label'] + ['leiden_wnn_{}'.format(r) for r in [1,2,3]],
                      add_metrics={'tfidf5': tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
# sctri.lazy_run(layer='raw',compute_metrics_parallel=False,
#                added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}],
#                assess_pruned=True,viewer_cluster=True,viewer_heterogeneity=True,viewer_heterogeneity_keys=['ori_label'])

# sctri = ScTriangulate.deserialize('dryrun_muon_output_tfidf3/after_rank_pruning.p')
# for k in ['leiden_wnn_{}'.format(r) for r in [0.5,1,2,3,4,5,6]]:
#     sctri.adata.obs['pruned'] = sctri.adata.obs[k]
#     sctri.viewer_heterogeneity_html(key='ori_label')
#     sctri.viewer_heterogeneity_figure(key='ori_label')
#     os.rename('dryrun_muon_output_tfidf3/figure4viewer','dryrun_muon_output_tfidf3/figure4viewer_hetero_{}'.format(k))

'''the test_run, per Xuan request, generate marker plots for all pruned cluster'''
sctri = ScTriangulate.deserialize('test_run_output_tfidf3_ori_1_2_3/after_pruned_assess.p')
# sctri.modality_contributions(mode='marker_genes',key='pruned',tops=20)
# for col in ['adt_contribution','rna_contribution']:
#     sctri.plot_umap(col,'continuous',umap_cmap='viridis')

# xuan_curation = pd.read_csv('Xuan_curation.csv',sep=',',index_col=0)
# xuan_curation['pruned'] = [item.split('__')[1] for item in xuan_curation.index]
# dic = pd.Series(index=xuan_curation['pruned'].values,data=xuan_curation['Curated annotation by ADT'].values).to_dict()
# os.mkdir(os.path.join(sctri.dir,'feature_plots'))
# for c in np.unique(sctri.adata.obs['pruned'].values):
#     cc = dic[c]
#     sctri.plot_multi_modal_feature_rank(cluster=c)
#     ori_file = os.path.join(sctri.dir,'sctri_multi_modal_feature_rank_{}_{}_{}_{}.{}'.format('marker_genes','pruned',c,20,'.pdf'))
#     new_file = os.path.join(sctri.dir,'feature_plots','sctri_multi_modal_feature_rank_{}_{}.pdf'.format(c,cc))
#     os.rename(ori_file,new_file)

# merge = [('leiden_wnn_1@0','leiden_wnn_1@1')]
# marker_gene_dict = {
#     'leiden_wnn_1@0+leiden_wnn_1@1':['SPINK2','AB_Hu.CD117'],
#     'leiden_wnn_1@8':['ZBTB20'],
#     'leiden_wnn_2@25':['AB_Hu.CD71','AB_Hu.CD82','AB_Hu.CD102']
# }
# sctri.plot_heterogeneity(key='ori_label',cluster='HSCP-4',style='heatmap+umap',
#                          subset=['leiden_wnn_1@0','leiden_wnn_1@1','leiden_wnn_1@8','leiden_wnn_2@25'],
#                          merge=merge,marker_gene_dict=marker_gene_dict)
# sctri.plot_heterogeneity(key='ori_label',cluster='HSCP-4',style='multi_gene',multi_gene=['SPINK2','ZBTB20','AB_Hu.CD71'])
# for gene in ['SPINK2','AB_Hu.CD117','ZBTB20','AB_Hu.CD71','AB_Hu.CD82','AB_Hu.CD102']:
#     sctri.plot_heterogeneity(key='ori_label',cluster='HSCP-4',style='single_gene',single_gene=gene,cmap='viridis')

curation = pd.read_csv('xuan_wnn.csv',index_col=0)
curation.index = [item.split('__')[1] for item in curation.index]
dic = curation['Note'].to_dict()
sctri.adata.obs['Note'] = sctri.adata.obs['pruned'].map(dic).astype('category').values
umap_dual_view_save(sctri.adata,cols=['Note'])
sys.exit('stop')



# # just want to test w/o adt
# adata_final_only_rna = adata_final[:,adata_final.var['modality']=='rna']
# sctri = ScTriangulate(dir='dryrun_muon_output_tfidf3_only_rna',
#                       adata=adata_final_only_rna,
#                       query=['ori_label'],
#                       add_metrics={'tfidf5': tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
# sctri.run_single_key_assessment(key='ori_label',scale_sccaf=False,layer=None,added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}])
# sctri.serialize('ori_label_stability_rna_only.p')



'''now if I want to build the soupx RNA count
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/BF21-CD34/BF21_CD34/outs/soupX-contamination-fraction-0.15/BF21_032123_CD34.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/BF21-CD271/BF21_CD271/outs/soupX-contamination-fraction-0.15/BF21_032123_CD271.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/BF21-TNC/BF21_TNC/outs/soupX-contamination-fraction-0.15/BF21_032123_TNC.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/WF26-CD34/WF26_CD34/outs/soupX-contamination-fraction-0.15/WF26_031423_CD34.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/WF26-CD271/WF26_CD271/outs/soupX-contamination-fraction-0.15/WF26_031423_CD271.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/WF26-TNC/WF26_TNC/outs/soupX-contamination-fraction-0.15/WF26_031423_TNC.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/BM27/BM27_CD34/BM27_CD34/outs/soupX-contamination-fraction-0.15/BM27_120522_CD34.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/BM27/BM27_CD271/BM27_CD271/outs/soupX-contamination-fraction-0.15/BM27_120522_CD271.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/BM27/BM27_TNC/BM27_TNC/outs/soupX-contamination-fraction-0.15/BM27_120522_TNC.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/WM34/WM34_CD34/WM34_CD34/outs/soupX-contamination-fraction-0.15/WM34_120522_CD34.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/WM34/WM34_CD271/WM34_CD271/outs/soupX-contamination-fraction-0.15/WM34_120522_CD271.txt
/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/WM34/WM34_TNC/WM34_TNC/outs/soupX-contamination-fraction-0.15/WM34_120522_TNC.txt
'''
# BF21_CD34 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/BF21-CD34/BF21_CD34/outs/soupX-contamination-fraction-0.15/BF21_032123_CD34.txt',gene_is_index=True,sep='\t')
# BF21_CD34.obs_names = [item + '.BF21_032123_CD34' for item in BF21_CD34.obs_names]
# print(BF21_CD34)
# BF21_CD271 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/BF21-CD271/BF21_CD271/outs/soupX-contamination-fraction-0.15/BF21_032123_CD271.txt',gene_is_index=True,sep='\t')
# BF21_CD271.obs_names = [item + '.BF21_032123_CD271' for item in BF21_CD271.obs_names]
# print(BF21_CD271)
# BF21_TNC = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/BF21-TNC/BF21_TNC/outs/soupX-contamination-fraction-0.15/BF21_032123_TNC.txt',gene_is_index=True,sep='\t')
# BF21_TNC.obs_names = [item + '.BF21_032123_TNC' for item in BF21_TNC.obs_names]
# print(BF21_TNC)

# BM27_CD34 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/BM27/BM27_CD34/BM27_CD34/outs/soupX-contamination-fraction-0.15/BM27_120522_CD34.txt',gene_is_index=True,sep='\t')
# BM27_CD34.obs_names = [item + '.BM27_120522_CD34' for item in BM27_CD34.obs_names]
# print(BM27_CD34)
# BM27_CD271 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/BM27/BM27_CD271/BM27_CD271/outs/soupX-contamination-fraction-0.15/BM27_120522_CD271.txt',gene_is_index=True,sep='\t')
# BM27_CD271.obs_names = [item + '.BM27_120522_CD271' for item in BM27_CD271.obs_names]
# print(BM27_CD271)
# BM27_TNC = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/BM27/BM27_TNC/BM27_TNC/outs/soupX-contamination-fraction-0.15/BM27_120522_TNC.txt',gene_is_index=True,sep='\t')
# BM27_TNC.obs_names = [item + '.BM27_120522_TNC' for item in BM27_TNC.obs_names]
# print(BM27_TNC)

# WF26_CD34 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/WF26-CD34/WF26_CD34/outs/soupX-contamination-fraction-0.15/WF26_031423_CD34.txt',gene_is_index=True,sep='\t')
# WF26_CD34.obs_names = [item + '.WF26_031423_CD34' for item in WF26_CD34.obs_names]
# print(WF26_CD34)
# WF26_CD271 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/WF26-CD271/WF26_CD271/outs/soupX-contamination-fraction-0.15/WF26_031423_CD271.txt',gene_is_index=True,sep='\t')
# WF26_CD271.obs_names = [item + '.WF26_031423_CD271' for item in WF26_CD271.obs_names]
# print(WF26_CD271)
# WF26_TNC = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/WF26-TNC/WF26_TNC/outs/soupX-contamination-fraction-0.15/WF26_031423_TNC.txt',gene_is_index=True,sep='\t')
# WF26_TNC.obs_names = [item + '.WF26_031423_TNC' for item in WF26_TNC.obs_names]
# print(WF26_TNC)

# WM34_CD34 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/WM34/WM34_CD34/WM34_CD34/outs/soupX-contamination-fraction-0.15/WM34_120522_CD34.txt',gene_is_index=True,sep='\t')
# WM34_CD34.obs_names = [item + '.WM34_120522_CD34' for item in WM34_CD34.obs_names]
# print(WM34_CD34)
# WM34_CD271 = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/WM34/WM34_CD271/WM34_CD271/outs/soupX-contamination-fraction-0.15/WM34_120522_CD271.txt',gene_is_index=True,sep='\t')
# WM34_CD271.obs_names = [item + '.WM34_120522_CD271' for item in WM34_CD271.obs_names]
# print(WM34_CD271)
# WM34_TNC = small_txt_to_adata(int_file='/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/WM34/WM34_TNC/WM34_TNC/outs/soupX-contamination-fraction-0.15/WM34_120522_TNC.txt',gene_is_index=True,sep='\t')
# WM34_TNC.obs_names = [item + '.WM34_120522_TNC' for item in WM34_TNC.obs_names]
# print(WM34_TNC)

# adata_raw_rna = ad.concat(adatas=[BF21_CD34,BF21_CD271,BF21_TNC,BM27_CD34,BM27_CD271,BM27_TNC,WF26_CD34,WF26_CD271,WF26_TNC,WM34_CD34,WM34_CD271,WM34_TNC],
#                           axis=0,join='inner',merge='first',label='sample',
#                           keys=['BF21_CD34','BF21_CD271','BF21_TNC','BM27_CD34','BM27_CD271','BM27_TNC','WF26_CD34','WF26_CD271','WF26_TNC','WM34_CD34','WM34_CD271','WM34_TNC'])
# adata_raw_rna = adata_raw_rna[adata_final.obs_names,:]  
# adata_raw_rna.write('to_seurat_raw_rna.h5ad')
# adata_raw_adt.write('to_seurat_raw_adt.h5ad')




'''
A bit test on Seurat results on 100k
'''
# add_umap(adata_final,inputs='wnn_umap.txt',mode='pandas_disk',cols=['wnnUMAP_1','wnnUMAP_2'],index_col=0,key='X_umap')
# add_annotations(adata_final,inputs='wnn_metadata.txt',cols_input=['RNA.weight','ADT.weight'] + ['wsnn_res.{}'.format(r) for r in [0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75]])
# # umap_dual_view_save(adata_final,cols=['ori_label','sample','donor'])
# adata_final.obs['RNA.weight'] = adata_final.obs['RNA.weight'].astype('float32')
# adata_final.obs['ADT.weight'] = adata_final.obs['ADT.weight'].astype('float32')
# for key in ['RNA.weight','ADT.weight']:
#     sc.pl.umap(adata_final,color=key)
#     plt.savefig('{}_umap.pdf'.format(key),bbox_inches='tight')
#     plt.close()
# sys.exit('stop')

# # running using 1 or 2 or 3 TFIDF scores
# sctri = ScTriangulate(dir='newrun_output_three_tfidf_ori0512',
#                       adata=adata_final,
#                       query=['ori_label','sctri_adt_leiden_0.5','sctri_adt_leiden_1','sctri_adt_leiden_2'],
#                       add_metrics={'tfidf5': tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
# sctri.lazy_run(layer='raw',compute_metrics_parallel=False,
#                added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}],
#                assess_pruned=True,viewer_cluster=True,viewer_heterogeneity=True,viewer_heterogeneity_keys=['ori_label'])


    
# downstream
# adata = sc.read('output_three_tfidf/sctriangulate.h5ad')
# umap_color_exceed_102(adata,key='pruned',outdir='output_three_tfidf',dot_size=0.5**2,edgecolor='none')
# ScTriangulate.salvage_run(step_to_start='assess_pruned',last_step_file='output_three_tfidf/after_rank_pruning.p',scale_sccaf=False,
#                           layer='raw',added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}],
#                           viewer_heterogeneity_keys=['ori_label'])


# ## pruned cluster, donor distribution
# sctri = ScTriangulate.deserialize('output_three_tfidf/after_pruned_assess.p')
# df = sctri.adata.obs.loc[:,['donor','pruned']]
# def custom_value_count(x):
#     x = x['donor']
#     result = []
#     vc = x.value_counts().to_dict()
#     # WF26
#     result.append(vc.get('WF26',0))
#     # WF83
#     result.append(vc.get('WF83',0))
#     # BF21
#     result.append(vc.get('BF21',0))
#     # BM27
#     result.append(vc.get('BM27',0))
#     # WM34
#     result.append(vc.get('WM34',0))
#     return result[0],result[1],result[2],result[3],result[4]
# result = df.groupby(by='pruned').apply(func=custom_value_count)
# result = pd.DataFrame(data=result.to_dict(),index=['WF26','WF83','BF21','BM27','WM34']).T
# result.to_csv('donor_origin_pruned_abs.txt',sep='\t')
# result = result / result.sum(axis=1).values
# print(result)
# result.to_csv('donor_origin_pruned_rel.txt',sep='\t')
# fig,ax = plt.subplots(figsize=(20,4.8))
# result.plot.bar(rot=90,fontsize=1,ax=ax,stacked=True)
# plt.savefig('donor_origin_pruned.pdf',bbox_inches='tight')
# plt.close()


# ## no 1. does ADT makes the rna-define cluster more stable?
# adata_only_rna = adata_combined_100k_rna_harmony
# sctri_rna = ScTriangulate(dir='output_three_tfidf_only_rna',adata=adata_only_rna,query=['ori_label'],add_metrics={'tfidf5':tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
# sctri_rna.run_single_key_assessment(key='ori_label',scale_sccaf=False,layer=None,added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}])
# sctri_rna.serialize('ori_label_rna_only_stability.p')


## no 2. does ADT makes the pruned cluster more stable?
# sctri = ScTriangulate.deserialize('output_three_tfidf/after_pruned_assess.p')
# adata_only_rna = sctri.adata[:,sctri.adata.var['modality']=='rna']
# sctri_rna = ScTriangulate(dir='output_three_tfidf_only_rna_pruned',adata=adata_only_rna,query=['pruned'],add_metrics={'tfidf5':tf_idf5_for_cluster,'tfidf1':tf_idf1_for_cluster})
# sctri_rna.run_single_key_assessment(key='pruned',scale_sccaf=False,layer=None,added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None},{'species': 'human', 'criterion': 2, 'layer': None}])
# sctri_rna.serialize('ori_label_rna_only_stability.p')


# sct_rna = ScTriangulate.deserialize('output_three_tfidf_only_rna/ori_label_rna_only_stability.p')
# sctri = ScTriangulate.deserialize('output_three_tfidf/after_pruned_assess.p')
# key = 'ori_label'
# for score in ['cluster_to_reassign','cluster_to_tfidf10','cluster_to_tfidf1','cluster_to_tfidf5','cluster_to_SCCAF']:
#     a = sct_rna.score[key][score]
#     b = sctri.score[key][score]
#     # b = {k.replace('@','_'):v for k,v in b.items()}
#     df = pd.Series(a,name='rna_only').to_frame()
#     df['rna_adt'] = df.index.map(b).values
#     df['delta'] = df['rna_adt'] - df['rna_only']
#     df.sort_values(by='delta',ascending=False,inplace=True)
#     df.to_csv('adt_gain_{}_{}.txt'.format(key,score),sep='\t')
#     df.drop(columns='delta',inplace=True)
#     fig,ax = plt.subplots(figsize=(20,4.8))
#     df.plot.bar(rot=90,fontsize=1,ax=ax)
#     # plt.savefig('adt_gain_{}_{}.pdf'.format(key,score),bbox_inches='tight')
#     # plt.close()
#     from scipy.stats import ttest_ind,ttest_rel
#     print(score)
#     print(ttest_ind(df['rna_only'].values,df['rna_adt'].values,alternative='less'))
#     print(ttest_rel(df['rna_only'].values,df['rna_adt'].values,alternative='less'))







sctri = ScTriangulate.deserialize('test_run_output_tfidf3_ori_1_2_3/after_pruned_assess.p')
# df = pd.read_csv('output_three_tfidf_only_rna/adt_gain_ori_label_cluster_to_SCCAF.txt',sep='\t',index_col=0)
xuan = pd.read_csv('xuan_wnn.csv',sep=',',index_col=0)
xuan = xuan.drop_duplicates(subset='Enrichment assigned')
xuan = xuan.loc[xuan['Enrichment assigned'].notna(),:]
xuan = xuan.loc[:,['Enrichment assigned','Note']]
xuan['Enrichment assigned'] = [item.replace(' ','_') for item in xuan['Enrichment assigned']]

# umap1
# c2d = df['delta'].to_dict()
# sctri.adata.obs['delta'] = sctri.adata.obs['ori_label'].map(c2d)
# col = []
# for item in sctri.adata.obs['delta']:
#     if item > 0:
#         col.append('increase')
#     elif item == 0:
#         col.append('stay_the_same')
#     else:
#         col.append('negligible_descrease')
# sctri.adata.obs['is_increase'] = col
# umap_dual_view_save(sctri.adata,cols=['is_increase'])

# umap2
xuan.set_index(keys='Enrichment assigned',inplace=True)
dic = xuan['Note'].to_dict()
col = []
for item in sctri.adata.obs['ori_label']:
    identity = dic.get(item,'unknown')
    if identity == 'subcluster':
        col.append('subcluster')
    else:
        col.append('no_subcluster')
sctri.adata.obs['subcluster'] = col
umap_dual_view_save(sctri.adata,cols=['subcluster'])

# umap3, we need all cluster, but grey all the one that are equivalant to RNA
xuan = pd.read_csv('xuan.csv',sep=',',index_col=0)
xuan = xuan.loc[xuan['Enrichment assigned'].notna(),:]
dic = xuan['Note'].to_dict()
def umap_color_exceed_102_update(adata,key,dot_size=None,legend_fontsize=6,outdir='.',name=None,**kwargs):
    fig,ax = plt.subplots()
    mapping = colors_for_set(adata.obs[key].unique().tolist())
    mapping = {k:'lightgrey' if dic[k] == 'equivalent' else v for k,v in mapping.items()}
    color = adata.obs[key].map(mapping).values
    if dot_size is None:
        dot_size = 120000/adata.shape[0]
    ax.scatter(adata.obsm['X_umap'][:,0],adata.obsm['X_umap'][:,1],c=color,s=dot_size,**kwargs)
    # add text
    adata.obs['umap_x'] = adata.obsm['X_umap'][:,0]
    adata.obs['umap_y'] = adata.obsm['X_umap'][:,1]
    tmp_df = adata.obs.loc[:,['pruned','umap_x','umap_y']].copy()
    for k,v in mapping.items():
        if v != 'lightgrey':
            mat = tmp_df.loc[tmp_df['pruned']==k,['umap_x','umap_y']]
            centroid = list(np.mean(mat,axis=0))
            ax.text(x=centroid[0],y=centroid[1],s=k,fontsize=1,ha='center',va='center')
    import matplotlib.lines as mlines
    ax.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in mapping.values()],
            labels=[i for i in mapping.keys()],loc='upper left',bbox_to_anchor=(1,1),ncol=3,frameon=False,prop={'size':6})
    if name is None:
        name = 'umap_{}_exceed_102.pdf'.format(key)
    plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
    plt.close()

umap_color_exceed_102_update(sctri.adata,key='pruned',outdir='output_three_tfidf',dot_size=0.5**2,edgecolor='none')



