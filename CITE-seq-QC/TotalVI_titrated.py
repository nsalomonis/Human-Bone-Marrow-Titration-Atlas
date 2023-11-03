#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 13:08:27 2022

@author: khauv3
"""

import pandas as pd
import numpy as np
import sys
import os
import csv
import gzip
import scipy.io
import anndata
import scvi
import scanpy as sc
import pickle
import torch


sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='pdf')

def load_gene_expression_matrix(matrix_dir):
    features_path = os.path.join(matrix_dir, "features.tsv.gz")
    feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
    feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]
    tmp_adata = anndata.read_mtx(os.path.join(matrix_dir, "matrix.mtx.gz")).T
    tmp_adata.obs = pd.DataFrame(index=barcodes)
    tmp_adata.var = pd.DataFrame({"gene": gene_names,
			    "id": feature_ids,
			    "feature_type": feature_types},
			    index=feature_ids)
    return(tmp_adata.copy())

#### 

def load_adt_matrix(path_adt_matrix):
    tmp_adt = pd.read_table(path_adt_matrix, index_col=0)
    tmp_adt.columns = pd.Series(tmp_adt.columns)
    return(tmp_adt.copy())

def load_adata_for_cite_seq(matrix_dir, path_adt_matrix, sample_id):
    print("Loading RNA data for {}...".format(sample_id))
    tmp_rna = load_gene_expression_matrix(matrix_dir)
    print("Loading ADT data for {}...".format(sample_id))
    tmp_adt = load_adt_matrix(path_adt_matrix)
    tmp_rna.obs.index = sample_id + "_" + pd.Series(tmp_rna.obs.index)
    tmp_adt.columns = sample_id + "_" + pd.Series(tmp_adt.columns)
    shared_cells = pd.Series(tmp_rna.obs.index)[list(pd.Series(tmp_rna.obs.index).isin(tmp_adt.columns))]
    filtered_adata = tmp_rna[shared_cells,:].copy()
    filtered_adata.obsm['protein_expression'] = tmp_adt.T.loc[shared_cells]
    print("Done.")
    print()
    return(filtered_adata.copy())
    
def combine_adata_objects(list_of_adata_objs):
  adts = pd.Series(list_of_adata_objs[0].obsm._data['protein_expression'].columns)
  for x in list_of_adata_objs:
    adts = adts[adts.isin(pd.Series(x.obsm._data['protein_expression'].columns))]
  for x in list_of_adata_objs:
    x.obsm._data['protein_expression'] = x.obsm._data['protein_expression'][adts]
  adata = list_of_adata_objs[0].concatenate(list_of_adata_objs[1:])
  return(adata.copy())




#### ---RNA Mtx directory paths --- ####  '/Volumes/salomonis2/Grimes/RNA/scRNASeq/10X-Genomics'
base_path             =   '/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/'
ext_BM27_CD34      =   'LGCHMC53-17GEX/BM27/BM27_CD34/BM27_CD34/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_BM27_CD271    =    'LGCHMC53-17GEX/BM27/BM27_CD271/BM27_CD271/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_BM27_TNC       =    'LGCHMC53-17GEX/BM27/BM27_TNC/BM27_TNC/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WM34_CD34      =    'LGCHMC53-17GEX/WM34/WM34_CD34/WM34_CD34/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WM34_CD271       =  'LGCHMC53-17GEX/WM34/WM34_CD271/WM34_CD271/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WM34_TNC       =    'LGCHMC53-17GEX/WM34/WM34_TNC/WM34_TNC/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WF83_TNC =   'LGCHMC60-17GEX-6ADT/WF83-TNC/WF83-TNC-hg19/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WF83_CD271 = 'LGCHMC60-17GEX-6ADT/WF83-CD271/WF83-CD271-hg19/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WF83_CD34  = 'LGCHMC60-17GEX-6ADT/WF83-CD34/WF83-CD34-hg19/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WF26_CD34  = 'LGCHMC60-17GEX-6ADT/WF26-CD34/WF26_CD34/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WF26_CD271  = 'LGCHMC60-17GEX-6ADT/WF26-CD271/WF26_CD271/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_WF26_TNC  = 'LGCHMC60-17GEX-6ADT/WF26-TNC/WF26_TNC/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_BF21_TNC = 'LGCHMC60-17GEX-6ADT/BF21-TNC/BF21_TNC/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_BF21_CD34 = 'LGCHMC60-17GEX-6ADT/BF21-CD34/BF21_CD34/outs/soupX-contamination-fraction-0.15/RNA_v3'
ext_BF21_CD271= 'LGCHMC60-17GEX-6ADT/BF21-CD271/BF21_CD271/outs/soupX-contamination-fraction-0.15/RNA_v3'




### --- # Example of loading corresponding ADT file --- ####
# Example of loading corresponding ADT file
path_adt_files         =    '/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/ADT'

adata_BM27_CD34     = load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_BM27_CD34),
					      path_adt_matrix = os.path.join(path_adt_files, 'BM27-CD34_ADT_clean.txt'),
					      sample_id = 'BM27_CD34')


adata_BM27_CD271   = load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_BM27_CD271 ),
					     path_adt_matrix = os.path.join(path_adt_files, 'BM27_CD271_ADT_clean.txt'),
					     sample_id = 'BM27_CD271')


adata_BM27_TNC     = load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_BM27_TNC),
					    path_adt_matrix = os.path.join(path_adt_files, 'BM27-TNC_ADT_clean.txt'),
					    sample_id = 'BM27_TNC')

adata_WM34_CD271   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WM34_CD271),
					      path_adt_matrix = os.path.join(path_adt_files, 'WM34-CD271_ADT_clean.txt'),
					      sample_id = 'WM34_CD271')



adata_WM34_CD34   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WM34_CD34),
					      path_adt_matrix = os.path.join(path_adt_files, 'WM34-CD34_ADT_clean.txt'),
					      sample_id = 'WM34_CD34')

adata_WM34_TNC   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WM34_TNC),
					      path_adt_matrix = os.path.join(path_adt_files, 'WM34_TNC_ADT_clean.txt'),
					      sample_id = 'WM34_TNC')
           
           
adata_WF83_CD271   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WF83_CD271),
                          path_adt_matrix = os.path.join(path_adt_files, 'WF83-CD271-hg19_ADT_clean.txt'),
                          sample_id = 'WF83_CD271')
                          
adata_WF83_CD34   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WF83_CD34),
                          path_adt_matrix = os.path.join(path_adt_files, 'WF83-CD34-hg19_ADT_clean.txt'),
                          sample_id = 'WF83_CD34')

adata_WF83_TNC   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WF83_TNC),
                          path_adt_matrix = os.path.join(path_adt_files, 'WF83-TNC-hg19_ADT_clean.txt'),
                          sample_id = 'WF83_TNC')



adata_WF26_CD34   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WF26_CD34),
                          path_adt_matrix = os.path.join(path_adt_files, 'WF26_CD34_ADT_clean.txt'),
                          sample_id = 'WF26_CD34')

adata_WF26_TNC   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WF26_TNC),
                          path_adt_matrix = os.path.join(path_adt_files, 'WF26_TNC_ADT_clean.txt'),
                          sample_id = 'WF26_TNC')
                          
                                     
adata_WF26_CD271   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_WF26_CD271),
                          path_adt_matrix = os.path.join(path_adt_files, 'WF26_CD271_ADT_clean.txt'),
                          sample_id = 'WF26_CD271')

adata_BF21_CD34   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_BF21_CD34),
                          path_adt_matrix = os.path.join(path_adt_files, 'BF21-CD34_ADT_clean.txt'),
                          sample_id = 'BF21_CD34')

adata_BF21_CD271   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_BF21_CD271),
                          path_adt_matrix = os.path.join(path_adt_files, 'BF21-CD271_ADT_clean.txt'),
                          sample_id = 'BF21_CD271')

adata_BF21_TNC   =  load_adata_for_cite_seq(matrix_dir = os.path.join(base_path, ext_BF21_TNC),
                          path_adt_matrix = os.path.join(path_adt_files, 'BF21-TNC_ADT_clean.txt'),
                          sample_id = 'BF21_TNC')


adata_WF83_CD34.obs['sample']="WF83_CD34" 
adata_WF83_CD271.obs['sample']="WF83_CD271"
adata_WF83_TNC.obs['sample']="WF83_TNC" 
adata_BF21_CD271.obs['sample']="BF21_CD271" 
adata_BF21_CD34.obs['sample']="BF21_CD34" 
adata_BF21_TNC.obs['sample']="BF21_TNC" 
adata_BM27_CD34.obs['sample']="BM27_CD34"  
adata_BM27_CD271.obs['sample']="BM27_CD271"   
adata_BM27_TNC.obs['sample']="BM27_TNC"   
adata_WM34_CD271.obs['sample']="WM34_CD271"  
adata_WM34_CD34.obs['sample']="WM34_CD34" 
adata_WM34_TNC.obs['sample']="WM34_TNC" 
adata_WF26_CD271.obs['sample']="WF26_CD271" 
adata_WF26_CD34.obs['sample']="WF26_CD34" 
adata_WF26_TNC.obs['sample']="WF26_TNC" 




list_adta_objs =[adata_BM27_CD34, adata_BM27_CD271, adata_BM27_TNC ,  adata_WM34_CD271 , adata_WM34_CD34,
                 adata_WM34_TNC,adata_WF83_CD34, adata_WF83_CD271, adata_WF83_TNC ,  adata_WF26_CD271 , adata_WF26_CD34,
                 adata_WF26_TNC, adata_BF21_CD271 , adata_BF21_CD34,
                 adata_BF21_TNC ]



print("Combining adata objects...")
adata = combine_adata_objects(list_adta_objs)
adata.write('adata.h5ad', compression="gzip")

print(adata.obs['sample'].value_counts())

print("Normalizing adata object...")
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

print("Finding highly variable genes...")
sc.pp.highly_variable_genes(
  adata,
  n_top_genes=4000,
  flavor="seurat_v3",
  batch_key="batch",
  subset=True,
  layer="counts"
)

print(adata.obs['sample'].value_counts())

print("Setting up anndata object for TOTALVI...")
scvi.model.TOTALVI.setup_anndata(
  adata,
  layer="counts",
  batch_key="batch",
  protein_expression_obsm_key="protein_expression"
)

model = scvi.model.TOTALVI(adata)
model.view_anndata_setup(adata)

print("Running TOTALVI...")
os.chdir('/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/total_VI_2')
vae = scvi.model.TOTALVI(adata,  latent_distribution="normal")
vae.train(max_epochs =400, train_size=0.8)

print("Adding batch corrected values to adata object...")
adata.obsm["X_totalVI"] = vae.get_latent_representation()
rna, protein = vae.get_normalized_expression(
  n_samples=20,
  return_mean=True,
  transform_batch=["0", "1", "2", "3", "4", "5","6","7","8","9","10","11","12","13","14"]
)

adata.layers["denoised_rna"], adata.obsm["denoised_protein"] = rna, protein

print("Adding protein foreground probability values...")
adata.obsm["protein_foreground_prob"] = vae.get_protein_foreground_probability( 
  n_samples=20,
  return_mean=True,
  transform_batch=["0", "1", "2", "3", "4", "5","6","7","8","9","10","11","12","13","14"]
)

print("Saving adata object...")
pickle_out = open("RC2_cite-seq_batch_correction_from_gpu.pickle", "wb")
pickle.dump(adata, pickle_out, protocol=4)
pickle_out.close()
				
# To load the data back into python
os.chdir('/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/total_VI_2')
with open("RC2_cite-seq_batch_correction_from_gpu.pickle", "rb") as f:
  adata = pickle.load(f)

sc.pp.neighbors(adata, use_rep="X_totalVI")
sc.tl.umap(adata, min_dist=0.4)
sc.tl.leiden(adata, key_added="leiden_totalVI")
sc.pl.umap(
  adata,
  color=["leiden_totalVI", "batch"],
  frameon=False,
  ncols=1,
  save="umap_plot_totalvi_signal_RC2-15_combined.pdf"
)


cell_anno = pd.concat([adata.obs, pd.DataFrame(adata.obsm["X_umap"],
					       columns=["umap_x", "umap_y"],
					       index=adata.obs.index)], axis=1)
cell_anno.to_csv("totalvi_signal_RC2_15__combined_cell_annotation.txt", sep="\t", header=True, index=True,
		 index_label="UID")

# Write out the original adt values
adata.obsm["protein_expression"].loc[adata.obs.index].to_csv("totalvi_raw_adt_values_RC2_15.txt", sep="\t",
							     header=True, index=True, index_label="UID")

# Write out the denoised adt values
adata.obsm["denoised_protein"].loc[adata.obs.index].to_csv("totalvi_denoised_adt_values_RC2_15.txt", sep="\t",
							     header=True, index=True, index_label="UID")




#### adding batch  information  as sample  name for plotting 
sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3) 
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png') 

cell_anno['Batch'] = cell_anno['batch']
cell_anno['Batch'].replace({'0':'BM27_CD34', '1':'BM27_CD127', '2':'BM27_TNC',  '3' : 'WM34_CD34', '4':'WM34_CD127',
                            '5': 'WM34_TNC','6':'WF83_CD271', '7':'WF83_TNC', '8':'WF83_CD34',  '9' : 'WF26_CD34',          
                            '10':'WF26_TNC','11': 'WF26_CD271','12':'BF21_CD34', '13':'BF21_CD127', '14':'BF21_TNC'},
                            inplace = True)


cell_anno.to_csv("totalvi_signal_RC2_15_cell_annotation.txt", sep="\t", header=True, index=True,
		 index_label="UID")

adata.obs['Batch']= cell_anno['Batch']
sc.pl.umap(
  adata,
  color=["leiden_totalVI", "Batch"],
  frameon=False,
  ncols=1,
  save="umap_plot_totalvi_signal_with-dataset-name.pdf"
)

### save as individual plots 
sc.pl.umap(adata, color="leiden_totalVI",legend_loc = 'on data', legend_fontsize =4, save = "clussters-umap.png")

sc.pl.umap(adata, color="Batch", palette='Set2', frameon=False, save = "Batch-corrected-umap-by-dataset-name.png")













