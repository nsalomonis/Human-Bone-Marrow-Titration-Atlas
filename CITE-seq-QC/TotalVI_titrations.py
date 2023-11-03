# bsub -q gpu-v100 -W 2:00 -n 8 -M 32000 -gpu "num=1" -Is bash
# module load cuda/10.1
# module load python3

#### Setup
print("Praparing packages...")
import csv
import gzip
import os
import scipy.io
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import collections
import scipy.sparse as sp_sparse
import tables
import anndata
import scvi
import scanpy as sc
import pickle



sample_names = ["TNC-3-4-2",
                "TNC-3-4-1",
                "TNC-1-2-2",
                "TNC-1-2-1",
                "CD34-4",
                "CD34-3",
                "CD34-2",
                "CD34-1",
                "CD271-4",
                "CD271-3",
                "CD271-2",
                "CD271-1"]


def make_list_items_unique(input_list):
  existing_items = {}
  output_list = []
  for tmp_item in input_list:
    if tmp_item in existing_items:
      output_list.append(tmp_item + "." + str(existing_items[tmp_item]))
      existing_items[tmp_item] += 1
    else:
      output_list.append(tmp_item)
      existing_items[tmp_item] = 1
  return(np.array(output_list))


def make_anndata_sprase_gex_from_gz_filtered_feature_directory(path_to_filtered_feature_dir, tmp_sample_name):
  print("Working on sample {}".format(tmp_sample_name))
  print("\tLoading h5 file...")
  sparse_mtx = scipy.io.mmread(gzip.open(os.path.join(path_to_filtered_feature_dir, "matrix.mtx.gz"))).tocsc()
  tmp_fanno = pd.read_table(gzip.open(os.path.join(path_to_filtered_feature_dir, "features.tsv.gz")), sep="\t", header=None)
  tmp_fanno.columns = ["id", "name", "feature_type"]
  tmp_gene_indices = np.where(tmp_fanno["feature_type"] == "Gene Expression")[0]
  tmp_adt_indices = np.where(tmp_fanno["feature_type"] == "Antibody Capture")[0]
  tmp_barcodes = pd.read_table(gzip.open(os.path.join(path_to_filtered_feature_dir, "barcodes.tsv.gz")), 
			       sep="\t", header=None).iloc[:,0].values
  tmp_barcode_df = pd.DataFrame({"cell": tmp_barcodes,
				"sample": tmp_sample_name},
				index=(pd.Series(tmp_barcodes) + "." + tmp_sample_name).values)
  tmp_feature_df = pd.DataFrame(tmp_fanno.iloc[tmp_gene_indices].values,
				index=make_list_items_unique(tmp_fanno.iloc[tmp_gene_indices]["name"].values),
				columns=tmp_fanno.columns.values)
  tmp_protein_df = pd.DataFrame(sparse_mtx[tmp_adt_indices,:].toarray().T,
				index=tmp_barcode_df.index.values,
				columns=make_list_items_unique(tmp_fanno.iloc[tmp_adt_indices]["name"].values))
  print("\tCreating AnnData object...")
  tmp_anndata = anndata.AnnData(X=pd.DataFrame.sparse.from_spmatrix(sparse_mtx[tmp_gene_indices,:].T,
					      index=tmp_barcode_df.index.values,
					      columns=tmp_feature_df.index.values),
				obs=tmp_barcode_df,
				var=tmp_feature_df)
  tmp_anndata.obsm["protein_expression"] = tmp_protein_df
  print("Done.")
  print("")
  return(tmp_anndata)




path_data_to_format = "/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/HBM_Titration_Xuan/220114_Grimes_GSL-PY-2582_Run1/{}/RNA/{}/outs/SoupXall/filtered_feature_bc_matrix/"
list_adata = [make_anndata_sprase_gex_from_gz_filtered_feature_directory(os.path.join(path_data_to_format.format(tmp_sample_name, tmp_sample_name)), 
                                                                          tmp_sample_name) for tmp_sample_name in sample_names]

combo = anndata.concat(list_adata, merge="same")

del list_adata

import gc 
gc.collect()


print("Filtering barcodes...")
# path_barcode_filtering_table = "/data/salomonis2/LabFiles/Kyle/Analysis/2022_02_12_run_totalvi_human_gsl_py_2582_run1/input/Titration-Barcodes.v1.txt"
# barcode_filtering_table = pd.read_table(path_barcode_filtering_table)
path_cell_metadata = "/data/salomonis2/LabFiles/Kyle/Analysis/2022_02_12_run_totalvi_human_gsl_py_2582_run1/input/Cells-Metadata.txt"

cell_metadata_df = pd.read_table(path_cell_metadata, header=0, index_col=0)


# Implementing custom cell filtering based on Seurat pipeline
# filtered_barcodes = np.intersect1d(barcode_filtering_table.loc[barcode_filtering_table["Filter.Seurat"] == 1]["uid"].values, combo.obs.index.values)
filtered_barcodes = np.intersect1d(cell_metadata_df.index.values, combo.obs.index.values)


combo = combo[filtered_barcodes]


# Add batch field
print("Preprocessing RNA content...")
combo.obs["Batch"] = cell_metadata_df.loc[combo.obs.index.values, "HTO"].str.replace("--ADT", "")
print("Adding umi to counts layers...")
combo.layers["counts"] = combo.X.copy()
print("Normalizing features...")
sc.pp.normalize_total(combo, target_sum=1e4)
print("Log transforming values...")
sc.pp.log1p(combo)


print("Finding Variable features...")
sc.pp.highly_variable_genes(combo, n_top_genes=4000, flavor="seurat_v3", batch_key="Batch", subset=True, layer="counts")
print("Setting up anndata object")
scvi.data.setup_anndata(combo, layer="counts", batch_key="Batch", protein_expression_obsm_key="protein_expression") 



#### Prepare and run model
print("Starting TotalVI training...")
vae = scvi.model.TOTALVI(combo, use_cuda=True, latent_distribution="ln")
vae.train(n_epochs=400, train_size=0.9, lr=0.004)

print("Writing out final object...")
with open("tmp_totalvi_vae_object_2022_03_10_ln_400_epochs_soupx_all_samples.pickle", "wb") as tmp_file:
  pickle.dump(vae, tmp_file, protocol=4)


print("Writing out temp anndata object...")
with open("tmp_totalvi_anndata_object_2022_03_10_ln_400_epochs_soupx_all_samples.pickle", "wb") as tmp_file:
  pickle.dump(combo, tmp_file, protocol=4)
  
  
#### Analyze outputs
print("Attaching denoised values...")
rna, protein = vae.get_normalized_expression()

print("Writing out denoised adt...")
protein.to_csv("denoised_adt_values_2021_03_09_ln_400_epochs_soupx_all_samples.txt", sep="\t", index=True, header=True, index_label="UID")
