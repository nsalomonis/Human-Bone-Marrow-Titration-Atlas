setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_BM_400k_RNA')

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#convert h5ad to h5seurat and load
Convert('./to_seurat_raw_rna.h5ad','to_seurat_raw_rna.h5seurat',assay='RNA')
Convert('./to_seurat_raw_adt.h5ad','to_seurat_raw_adt.h5seurat',assay='ADT')

seurat_rna <- LoadH5Seurat('./to_seurat_raw_rna.h5seurat')
new_seurat_rna <- CreateSeuratObject(counts=seurat_rna@assays[['RNA']]@counts, project='Xuan', min.cells=3, min.features=2)

seurat_adt <- LoadH5Seurat('./to_seurat_raw_adt.h5seurat')
seurat_cite <- new_seurat_rna
seurat_cite@assays$ADT = seurat_adt@assays$ADT
DefaultAssay(seurat_cite) <- 'RNA'

#QC filtering in Seurat
seurat_cite <- PercentageFeatureSet(seurat_cite,pattern='^MT-',col.name = 'precent.mt')
seurat_cite <- SCTransform(seurat_cite,method='glmGamPoi',vars.to.regress = 'precent.mt',verbose=TRUE)
seurat_cite <- RunPCA(seurat_cite,verbose=TRUE)

DefaultAssay(seurat_cite) <- 'ADT'
VariableFeatures(seurat_cite) <- rownames(seurat_cite[["ADT"]])
seurat_cite <- NormalizeData(seurat_cite, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

# Find neighbors
seurat_cite <- FindMultiModalNeighbors(seurat_cite,reduction.list=list('pca','apca'),
                        dims.list=list(1:50,1:18))
seurat_cite <- RunSPCA(seurat_cite,assay='SCT',graph=seurat_cite@graphs[['wsnn']])
saveRDS(seurat_cite,'./build_cite_ref/file.rds')

seurat_cite <- readRDS('./build_cite_ref/file.rds')

# for Azumith to keep similar UMAP projection, modify weighted.nn and supply to RunsPCA, then re-run spca
# https://rdrr.io/github/mojaveazure/seurat-object/man/as.Graph.html
library(Matrix)
library(Azimuth)
library(presto)

wnn_kernel <- Matrix::readMM('./build_cite_ref/wnn_kernel.mtx')
rownames(wnn_kernel) <- read.table('./build_cite_ref/barcodes.tsv',sep='\t',header=F)$V1
colnames(wnn_kernel) <- read.table('./build_cite_ref/barcodes.tsv',sep='\t',header=F)$V1
g <- as.Graph(x = wnn_kernel)


neighbor <- as.Neighbor(g) # only 16 dims
  
seurat_cite@graphs[['wsnn']] = g
seurat_cite@neighbors[['weighted.nn']] = neighbor

# make sure the cell type is of the column name Clusters, maybe here you can add more
meta <- read.table('test_run_output_tfidf3_ori_1_2_3/sctri_barcode2cellmetadata-curated.txt',header=T,row.names=1,sep='\t')  # provide cellmetadata file
cells <- colnames(seurat_cite@assays[['RNA']]@counts)
meta_need <- meta[cells,c('Level1','Level2','Level3R','Level3M')]  # modify if more levels
new_meta <- cbind(seurat_cite@meta.data,meta_need)
seurat_cite@meta.data <- new_meta

seurat_cite <- RunSPCA(seurat_cite,assay='SCT',graph=seurat_cite@graphs[['wsnn']])
seurat_cite <- RunUMAP(seurat_cite,nn.name='weighted.nn', 
                       reduction.name='umap',reduction.key='UMAP_',umap.method="uwot",return.model=TRUE)

#check UMAP
DimPlot(seurat_cite, reduction = 'umap', group.by = 'Level3M', label = T, 
        repel = TRUE, label.size = 0.5) + NoLegend()

#create Azimuth
ref <- AzimuthReference(
  object = seurat_cite,
  refUMAP = "umap",
  refDR = "spca",
  refAssay = "SCT",
  assays = c("ADT"),
  metadata = c('Level1','Level2','Level3R','Level3M'),    # modify if more levels
  dims = 1:50,
  k.param = 31
)

ref.dir = "./build_cite_ref/ref_adt_xuan_4levels"   # modify
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))





