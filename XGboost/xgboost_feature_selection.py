# module load python3/3.7.1

import os
import pandas as pd
import numpy as np
from xgboost import XGBClassifier

path_analysis = "/data/salomonis2/LabFiles/Kyle/Analysis/2022_09_01_human_cite_titration_nominate_markers/"
path_adts = "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs-BoneMarrow-Titration_Altas/TotalVI-ADTs/denoised_adt_values_2021_03_09_400_epochs-transposed.txt"
path_adts_feather = "input/human_cite_titration_denoised_adts.fea"
path_groups = "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs-BoneMarrow-Titration_Altas/output_two_5rna/barcode2cellmetadata.txt"
path_titration_annotation = "input/titration_groups.hybrid.txt"

os.chdir(path_analysis)

# Move denoised adt matrix to a feather file so it can be loaded faster
# adt = pd.read_table(path_adts, index_col=0)
# adt.T.reset_index().to_feather(path_adts_feather)

# Load the adt matrix
adt = pd.read_feather(path_adts_feather).set_index("index")

# Load the groups file
groups = pd.read_table(path_groups, index_col=0)

# Read in the annotation information for the titration
titration_anno = pd.read_table(path_titration_annotation, header=None, index_col=0)

# Intersect the cell annotation information
shared_cells = np.intersect1d(np.intersect1d(adt.index.values, 
					     groups.index.values), 
				titration_anno.index.values)

cell_anno = pd.DataFrame({"sctri": groups.loc[shared_cells, 'pruned'].values,
			  "titration": [item.split("-")[-1] for item in titration_anno.loc[shared_cells].iloc[:,1].values]},
			 index=shared_cells)

cells_0_25 = cell_anno.loc[cell_anno["titration"] == '0.25'].index.values
cells_0_5 = cell_anno.loc[cell_anno["titration"] == '0.5'].index.values
cells_1 = cell_anno.loc[cell_anno["titration"] == '1'].index.values
cells_2 = cell_anno.loc[cell_anno["titration"] == '2'].index.values
cells_4 = cell_anno.loc[cell_anno["titration"] == '4'].index.values

def build_xgboost_model_feature_importance_df(xgb_model, input_pandas_df):
    init_series = "f" + pd.Series(list(range(input_pandas_df.shape[1])), index=input_pandas_df.columns).astype(str)
    tmp_weight_dict = xgb_model._Booster.get_score(importance_type="weight")
    tmp_gain_dict = xgb_model._Booster.get_score(importance_type="gain")
    tmp_cover_dict = xgb_model._Booster.get_score(importance_type="cover")
    tmp_totalgain_dict = xgb_model._Booster.get_score(importance_type="total_gain")
    tmp_totalcover_dict = xgb_model._Booster.get_score(importance_type="total_cover")
    return(pd.DataFrame({"weight": init_series.apply(lambda x: tmp_weight_dict[x]),
                         "gain": init_series.apply(lambda x: tmp_gain_dict[x]),
                         "cover": init_series.apply(lambda x: tmp_cover_dict[x]),
                         "total_gain": init_series.apply(lambda x: tmp_totalgain_dict[x]),
                         "total_cover": init_series.apply(lambda x: tmp_totalcover_dict[x])}, index=init_series.index))

def process_model(input_df, input_classes, path_save_feature_importance_df):
    tmp_model = XGBClassifier(njobs=10, verbosity=2)
    tmp_model.fit(input_df.values, input_classes)
    tmp_fimps = build_xgboost_model_feature_importance_df(tmp_model, input_df)
    tmp_fimps.sort_values(by="gain", ascending=False).to_csv(path_save_feature_importance_df,
							     sep="\t", header=True)

process_model(input_df = adt.loc[cells_0_25,:], 
	      input_classes = cell_anno.loc[cells_0_25, 'sctri'].values, 
	      path_save_feature_importance_df = "output/feature_importance_df_titration_0_25.tsv")

process_model(input_df = adt.loc[cells_0_5,:], 
	      input_classes = cell_anno.loc[cells_0_5, 'sctri'].values, 
	      path_save_feature_importance_df = "output/feature_importance_df_titration_0_5.tsv")

process_model(input_df = adt.loc[cells_1,:], 
	      input_classes = cell_anno.loc[cells_1, 'sctri'].values, 
	      path_save_feature_importance_df = "output/feature_importance_df_titration_1.tsv")

process_model(input_df = adt.loc[cells_2,:], 
	      input_classes = cell_anno.loc[cells_2, 'sctri'].values, 
	      path_save_feature_importance_df = "output/feature_importance_df_titration_2.tsv")

process_model(input_df = adt.loc[cells_4,:], 
	      input_classes = cell_anno.loc[cells_4, 'sctri'].values, 
	      path_save_feature_importance_df = "output/feature_importance_df_titration_4.tsv")












