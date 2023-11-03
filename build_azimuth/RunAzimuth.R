setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_BM_400k_RNA')

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(SeuratDisk)


library(Matrix)
library(Azimuth)
library(presto)

sessionInfo()

#terminal login to internet
#proxy_on
#username
#password
#curl -I https://linuxhint.com


# run azimuth
query <- RunAzimuth(query='/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs-BoneMarrow-Titration_Altas/adata_300k_raw.h5ad',
                    reference='/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_BM_400k_RNA/build_cite_ref_new',
                    do.adt = TRUE,
                    annotation.levels='Level3M') #could use other level c()

saveRDS(query,'./build_cite_ref_new/query.rds')

# visualize and save result
query <- readRDS('./build_cite_ref_new/query.rds')
DimPlot(query, reduction = 'ref.umap', group.by = 'predicted.Level3M', label = T, 
         repel = TRUE, label.size = 0.5) + NoLegend()
write.table(query@meta.data,'level3M_mapping_results.txt',sep='\t',col.names=NA)

# get ADT information
counts= query@assays$impADT@data
counts <- counts[!(grepl("Isotype", rownames(counts)) | grepl("Hu.Galectin.9", rownames(counts))), ]
# change to correct cellbarcodes
uid= colnames(counts)
uid <- gsub("\\.", "-", uid)
uid <- gsub("-1-", "-1.", uid)
# change to correct uid dash
ab= rownames(counts)
ab= gsub("-", "_", ab)
rownames(counts)= ab
# Create a new matrix with column names as the first row
counts_dense <- as.matrix(counts)
counts_dense <- rbind(uid, counts_dense)
write.table(counts_dense,file='/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_BM_400k_RNA/build_cite_ref/impADT300k.txt',sep='\t',row.names = T,col.names=F,quote = F)

#export matrix
a = query@assays[['impADT']]
writeMM(a@data,file='impADT300k.mtx')
write.table(a@data@Dimnames[[1]],file='impADT300k_features.tsv',sep='\t',row.names = F,col.names=F,quote = F)
write.table(a@data@Dimnames[[2]],file='impADT300k_barcodes.tsv',sep='\t',row.names = F,col.names=F,quote = F)
write.table(a@data@Dimnames[[2]],file='impADT300k_barcodes.tsv',sep='\t',row.names = F,col.names=F,quote = F)

