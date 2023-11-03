library(Seurat)

setwd("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/HBM_Titration_Xuan/220114_Grimes_GSL-PY-2582_Run1/TNC-3-4-2/RNA/TNC-3-4-2/outs/")

input.data= Read10X_h5("./filtered_feature_bc_matrix.h5")
RNA.counts= input.data$`Gene Expression`
data.seurat= CreateSeuratObject(counts = RNA.counts)
data.seurat[["percent.mt"]] <- PercentageFeatureSet(data.seurat, pattern = "^MT-")
data.seurat= subset(x= data.seurat, subset = nFeature_RNA>500&nCount_RNA>1000&percent.mt<25)
write.table(colnames(data.seurat@assays$RNA@counts), file = "seurat_filteredBC.txt", row.names = F, col.names = F)

