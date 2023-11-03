#the workflow reads in CITE-seq then flow data, please modify as needed
library(cyCombine)
library(tidyverse)
library(Seurat)
library(emdist)
library(cowplot)
library(CATALYST)
library(gridExtra)


# working directory
data_dir <- "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/merged"
# Set a seed to use throughout
seed = 840

# Get some nice colors for plotting
color_clusters <- c(RColorBrewer::brewer.pal(12, 'Paired'), RColorBrewer::brewer.pal(8, 'Dark2'))

expanded_colors <- rep(color_clusters, length.out = 64)

#For WNN clustered cells, read in txt and make into tibble(), bypassing Seurat object
BF21_032123_WNN100 <- read.delim("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/TFID3_ADT_centroids/WNN-new/inputfile_AS/exp.BF21-100.txt", header = TRUE, stringsAsFactors = FALSE)

BM27_120522_WNN100 <- read.delim("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/TFID3_ADT_centroids/WNN-new/inputfile_AS/exp.BM27-100.txt", header = TRUE, stringsAsFactors = FALSE)

WF26_031423_WNN100 <- read.delim("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/TFID3_ADT_centroids/WNN-new/inputfile_AS/exp.WF26-100.txt", header = TRUE, stringsAsFactors = FALSE)

WM34_120522_WNN100 <- read.delim("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC60-17GEX-6ADT/Total_VI_RC2_files/Total_vi_output/TFID3_ADT_centroids/WNN-new/inputfile_AS/exp.WM34-100.txt", header = TRUE, stringsAsFactors = FALSE)

#For CD34_ comparisons, read in Seurat object of CD34 filtered cells for each donor
BF21_CD34_Cite <- readRDS("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/Merged_validation4/input_files_forAS/new_input/BF21_CD34-filtered.rds")
BM27_CD34_Cite <- readRDS("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/Merged_validation4/input_files_forAS/new_input/BM27_CD34-filtered.rds")
WF26_CD34_Cite <- readRDS("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/Merged_validation4/input_files_forAS/new_input/WF26_CD34-filtered.rds")
WM34_CD34_Cite <- readRDS("/data/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/LGCHMC53-17GEX/Merged_validation4/input_files_forAS/new_input/WM34_CD34-filtered.rds")


rownames(BF21_032123_WNN100) <- BF21_032123_WNN100$uid
BF21 <-BF21_032123_WNN100[,-1]

rownames(BM27_120522_WNN100) <- BM27_120522_WNN100$uid
BM27 <-BM27_120522_WNN100[,-1]

rownames(WF26_031423_WNN100) <- WF26_031423_WNN100$uid
WF26 <-WF26_031423_WNN100[,-1]

rownames(WM34_120522_WNN100) <- WM34_120522_WNN100$uid
WM34 <-WM34_120522_WNN100[,-1]

BF21_tibble <-as_tibble(t(BF21))
BM27_tibble <-as_tibble(t(BM27))
WF26_tibble <-as_tibble(t(WF26))
WM34_tibble <-as_tibble(t(WM34))

#read in Seurat objects created during single-cell processing
#ADT counts have been filtered
BF21_032123_cite <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/old/Seurat_objects/BF21_032123TNC-filtered.rds")
BM27_120522_cite <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/old/Seurat_objects/BM27_120522TNC-filtered.rds")
WF26_031423_cite <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/old/Seurat_objects/WF26_031423TNC-filtered.rds")
WM34_120522_cite <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/old/Seurat_objects/WM34_120522TNC-filtered.rds")

#export ADT names to compare with flow, generate overlap list
write.csv(rownames(BF21_032123_cite@assays$ADT), file = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/to_merge_fcs/BF21_row_names.csv")

overlap_markers_df <- read.csv("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/overlap_markers.csv", header = FALSE)
overlap_markers<- overlap_markers_df$V1

#normalize ADTs
BF21_032123_cite <- NormalizeData(BF21_032123_cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)
BM27_120522_cite <- NormalizeData(BM27_120522_cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)
WF26_031423_cite <- NormalizeData(WF26_031423_cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)
WM34_120522_cite <- NormalizeData(WM34_120522_cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)

#cd34
BF21_CD34_Cite <- NormalizeData(BF21_CD34_Cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)

BM27_CD34_Cite <- NormalizeData(BM27_CD34_Cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)

WF26_CD34_Cite <- NormalizeData(WF26_CD34_Cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)

WM34_CD34_Cite <- NormalizeData(WM34_CD34_Cite, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = F) %>%
  ScaleData(verbose = F)


# Extract CITE-seq data (normalized) and add cell barcodes
BF21_032123_cite_data <- BF21_032123_cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
# Add the "id" column with cell barcodes
BF21_032123_cite_data <- BF21_032123_cite_data %>%
  mutate(id = colnames(BF21_032123_cite))


BM27_120522_cite_data <- BM27_120522_cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
# Add the "id" column with cell barcodes
BM27_120522_cite_data <- BM27_120522_cite_data %>%
  mutate(id = colnames(BM27_120522_cite))

WF26_031423_cite_data <- WF26_031423_cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
# Add the "id" column with cell barcodes
WF26_031423_cite_data <- WF26_031423_cite_data %>%
  mutate(id = colnames(WF26_031423_cite))

WM34_120522_cite_data <- WM34_120522_cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
# Add the "id" column with cell barcodes
WM34_120522_cite_data <- WM34_120522_cite_data %>%
  mutate(id = colnames(WM34_120522_cite))

#CD34 data
BF21_CD34_Cite <- BF21_CD34_Cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
BM27_CD34_Cite <- BM27_CD34_Cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
WF26_CD34_Cite <- WF26_CD34_Cite@assays$ADT@data %>%
  t() %>%
  as_tibble()
WM34_CD34_Cite <- WM34_CD34_Cite@assays$ADT@data %>%
  t() %>%
  as_tibble()


# Add columns to tibble for cyCombine processing
BF21_032123_cite_data <- BF21_032123_cite_data %>%
  mutate('batch' = 'BF21_032123_cite_data',
         'sample' = 'BF21_032123_cite',
         'label' = NULL)
BM27_120522_cite_data <- BM27_120522_cite_data %>%
  mutate('batch' = 'BM27_120522_cite_data',
         'sample' = 'BM27_120522_cite',
         'label' = NULL)
WF26_031423_cite_data <- WF26_031423_cite_data %>%
  mutate('batch' = 'WF26_031423_cite_data',
         'sample' = 'WF26_031423_cite',
         'label' = NULL)
WM34_120522_cite_data <- WM34_120522_cite_data %>%
  mutate('batch' = 'WM34_120522_cite_data',
         'sample' = 'WM34_120522_cite',
         'label' = NULL)

#CD34 data
BF21_CD34_Cite <- BF21_CD34_Cite %>%
  mutate('batch' = 'BF21_032123_CD34_cite',
         'sample' = 'BF21_032123_CD34_cite',
         'label' = NULL)

BM27_CD34_Cite <- BM27_CD34_Cite %>%
  mutate('batch' = 'BM27_120522_CD34_cite',
         'sample' = 'BM27_120522_CD34_cite',
         'label' = NULL)

WF26_CD34_Cite <- WF26_CD34_Cite %>%
  mutate('batch' = 'WF26_031423_CD34_cite',
         'sample' = 'WF26_031423_CD34_cite',
         'label' = NULL)

WM34_CD34_Cite <- WM34_CD34_Cite %>%
  mutate('batch' = 'WM34_120522_CD34_cite',
         'sample' = 'WM34_120522_CD34_cite',
         'label' = NULL)

#alternative if read in as txt above
BF21_tibble <- BF21_tibble %>%
  mutate('batch' = 'BF21_032123_cite',
         'sample' = 'BF21_032123_cite',
         'label' = NULL)

BM27_tibble <- BM27_tibble %>%
  mutate('batch' = 'BM27_120522_cite',
         'sample' = 'BM27_120522_cite',
         'label' = NULL)

WF26_tibble <- WF26_tibble %>%
  mutate('batch' = 'WF26_0314233_cite',
         'sample' = 'WF26_031423_cite',
         'label' = NULL)

WM34_tibble <- WM34_tibble %>%
  mutate('batch' = 'WM34_120522_cite',
         'sample' = 'WM34_120522_cite',
         'label' = NULL)


#read in the 4 donor's flow data

BF21_flowset <- compile_fcs("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BF21")

BM27_flowset <- compile_fcs("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BM27")

WF26_flowset <- compile_fcs("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WF26")

WM34_flowset <- compile_fcs("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WM34")

raw_data_BF21 <-convert_flowset(
  BF21_flowset,
  metadata = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BF21/metadata_BF21.csv", 
  filename_col = "Filename",
  batch_ids = "Batch",
  sample_ids = "Patient_id",
  clean_colnames = FALSE,
  down_sample = FALSE,
)

raw_data_BM27 <-convert_flowset(
  BM27_flowset,
  metadata = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BM27/metadata_BM27.csv", 
  filename_col = "Filename",
  batch_ids = "Batch",
  sample_ids = "Patient_id",
  clean_colnames = FALSE,
  down_sample = FALSE,
)

raw_data_WF26 <-convert_flowset(
  WF26_flowset,
  metadata = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WF26/metadata_WF26.csv", 
  filename_col = "Filename",
  batch_ids = "Batch",
  sample_ids = "Patient_id",
  clean_colnames = FALSE,
  down_sample = FALSE,
)

raw_data_WM34 <-convert_flowset(
  WM34_flowset,
  metadata = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WM34/metadata_WM34.csv", 
  filename_col = "Filename",
  batch_ids = "Batch",
  sample_ids = "Patient_id",
  clean_colnames = FALSE,
  down_sample = FALSE,
)

#read in 4 donors flow data as one flow set
donors4_flowset <- compile_fcs("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/merged")


raw_data <- convert_flowset(
  donors4_flowset,
  metadata = paste0(data_dir, "/metadata_4donors_flow.csv"), 
  filename_col = "Filename",
  batch_ids = "Batch",
  sample_ids = "Patient_id",
  clean_colnames = FALSE,
  down_sample = FALSE,
)

fcs_markers <- get_markers(raw_data)
uncorrected_4donors <- transform_asinh(raw_data, markers = fcs_markers,cofactor=6000, derand = FALSE)
write.csv(colnames(uncorrected_4donors), file = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/col_names.csv")
modified_colnames <- read.csv("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/col_names_new.csv", header = FALSE)
colnames(uncorrected_4donors) <- modified_colnames$V1

saveRDS(uncorrected_4donors, file = file.path(merge_with_cite_fcs_dir, "uncorrected_4donors_flow.RDS"))

#7/17/23
uncorrected_BF21 <- transform_asinh(raw_data_BF21, markers = fcs_markers,cofactor=6000, derand = FALSE)
uncorrected_BM27 <- transform_asinh(raw_data_BM27, markers = fcs_markers,cofactor=6000, derand = FALSE)
uncorrected_WF26 <- transform_asinh(raw_data_WF26, markers = fcs_markers,cofactor=6000, derand = FALSE)
uncorrected_WM34 <- transform_asinh(raw_data_WM34, markers = fcs_markers,cofactor=6000, derand = FALSE)

colnames(uncorrected_BF21) <- modified_colnames$V1
colnames(uncorrected_BM27) <- modified_colnames$V1
colnames(uncorrected_WF26) <- modified_colnames$V1
colnames(uncorrected_WM34) <- modified_colnames$V1


#make all CITE and flow to merge as one big tibble (original code)
uncorrected_total <- bind_rows(uncorrected_4donors[,overlap_markers], BF21_CD34_cite[,overlap_markers],BM27_CD34_cite[,overlap_markers],WF26_CD34_cite[,overlap_markers],WM34_CD34_cite[,overlap_markers])%>%
  mutate(id = 1:n())

#updated 08/23/23 to combine datasets while preserving cell barcodes
combined_data <- bind_rows(
  uncorrected_4donors[, overlap_markers],
  BF21_032123_cite_data[, overlap_markers],
  BM27_120522_cite_data[, overlap_markers],
  WF26_031423_cite_data[, overlap_markers],
  WM34_120522_cite_data[, overlap_markers]
) %>%
  mutate(dataset = rep(c("uncorrected_4donors", "BF21", "BM27", "WF26", "WM34"), 
                       times = c(nrow(uncorrected_4donors), 
                                 nrow(BF21_032123_cite_data),
                                 nrow(BM27_120522_cite_data),
                                 nrow(WF26_031423_cite_data),
                                 nrow(WM34_120522_cite_data))),
         cell_barcode = c(rownames(uncorrected_4donors),
                          BF21_032123_cite_data$id,
                          BM27_120522_cite_data$id,
                          WF26_031423_cite_data$id,
                          WM34_120522_cite_data$id)) %>%
  mutate(id = 1:n())  # Assign a unique ID


#make on tibble per donor
uncorrected_total_BF21 <- bind_rows(uncorrected_BF21[,overlap_markers], BF21_CD34_Cite[,overlap_markers])%>%
  mutate(id = 1:n())
uncorrected_total_BM27 <- bind_rows(uncorrected_BM27[,overlap_markers], BM27_CD34_Cite[,overlap_markers])%>%
  mutate(id = 1:n())
uncorrected_total_WF26 <- bind_rows(uncorrected_WF26[,overlap_markers], WF26_CD34_Cite[,overlap_markers])%>%
  mutate(id = 1:n())
uncorrected_total_WM34 <- bind_rows(uncorrected_WM34[,overlap_markers], WM34_CD34_Cite[,overlap_markers])%>%
  mutate(id = 1:n())

saveRDS(combined_data,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_updated_082323/uncorrected_total.rds")

# Run batch correction
corrected_total_BF21 <- uncorrected_total_BF21 %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

corrected_total_BM27 <- uncorrected_total_BM27 %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

corrected_total_WF26 <- uncorrected_total_WF26 %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

corrected_total_WM34 <- uncorrected_total_WM34 %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

corrected_total <- uncorrected_total %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

corrected_total <- combined_data %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

corrected_total$cell_barcode <- combined_data$cell_barcode

saveRDS(corrected_total,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_updated_082323/corrected_total.rds")


#evaluate batch correction
# Generate labels for EMD evaluation
labels <- corrected_total %>%
  create_som(seed = seed,
             xdim = 8,
             ydim = 8,
             markers = get_markers(corrected_total)[1:132])
test_corrected <- corrected_total
test_corrected$som <-"som"

labels_BF21 <- corrected_total_BF21 %>%
  create_som(seed = seed,
             xdim = 8,
             ydim = 8,
             markers = get_markers(corrected_total_BF21))
test_corrected_BF21 <- corrected_total_BF21
test_corrected_BF21$som <-"som"

labels_BM27 <- corrected_total_BM27 %>%
  create_som(seed = seed,
             xdim = 8,
             ydim = 8,
             markers = get_markers(corrected_total_BM27))
test_corrected_BM27 <- corrected_total_BM27
test_corrected_BM27$som <-"som"

labels_WF26 <- corrected_total_WF26 %>%
  create_som(seed = seed,
             xdim = 8,
             ydim = 8,
             markers = get_markers(corrected_total_WF26))
test_corrected_WF26 <- corrected_total_WF26
test_corrected_WF26$som <-"som"

labels_WM34 <- corrected_total_WM34 %>%
  create_som(seed = seed,
             xdim = 8,
             ydim = 8,
             markers = get_markers(corrected_total_WM34))
test_corrected_WM34 <- corrected_total_WM34
test_corrected_WM34$som <-"som"


# Add labels
test_corrected <- test_corrected %>%
  dplyr::mutate(som = labels)
celltype_col <- "som"

test_corrected_BF21 <- test_corrected_BF21 %>%
  dplyr::mutate(som = labels_BF21)
celltype_col <- "som"

test_corrected_BM27<- test_corrected_BM27 %>%
  dplyr::mutate(som = labels_BM27)
celltype_col <- "som"

test_corrected_WF26<- test_corrected_WF26 %>%
  dplyr::mutate(som = labels_WF26)

test_corrected_WM34<- test_corrected_WM34 %>%
  dplyr::mutate(som = labels_WM34)

test_uncorrected <-uncorrected_total
test_uncorrected$id <- seq(1, nrow(test_uncorrected))

test_uncorrected <-combined_data
test_uncorrected$id <- seq(1, nrow(test_uncorrected))

test_uncorrected_BF21 <-uncorrected_total_BF21
test_uncorrected_BF21$id <- seq(1, nrow(test_uncorrected_BF21))

test_uncorrected_BM27 <-uncorrected_total_BM27
test_uncorrected_BM27$id <- seq(1, nrow(test_uncorrected_BM27))

test_uncorrected_WF26 <-uncorrected_total_WF26
test_uncorrected_WF26$id <- seq(1, nrow(test_uncorrected_WF26))

test_uncorrected_WM34 <-uncorrected_total_WM34
test_uncorrected_WM34$id <- seq(1, nrow(test_uncorrected_WM34))

test_uncorrected <- test_corrected %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(test_uncorrected, by="id")

test_uncorrected_BF21 <- test_corrected_BF21 %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(test_uncorrected_BF21, by="id")

test_uncorrected_BM27 <- test_corrected_BM27 %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(test_uncorrected_BM27, by="id")

test_uncorrected_WF26 <- test_corrected_WF26 %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(test_uncorrected_WF26, by="id")

test_uncorrected_WM34 <- test_corrected_WM34 %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(test_uncorrected_WM34, by="id")

# Evaluate Earth Movers Distance and show plots
emd <- test_uncorrected %>%
  evaluate_emd(test_corrected, markers = get_markers(test_corrected)[1:132], cell_col = celltype_col)
cowplot::plot_grid(emd$violin, emd$scatterplot)
ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_updated_082323/emd_correction.pdf")


emd_BF21 <- test_uncorrected_BF21 %>%
  evaluate_emd(test_corrected_BF21, markers = get_markers(test_corrected_BF21), cell_col = celltype_col)
cowplot::plot_grid(emd_BF21$violin, emd_BF21$scatterplot)
ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BF21/emd_correction.pdf")

emd_BM27 <- test_uncorrected_BM27 %>%
  evaluate_emd(test_corrected_BM27, markers = get_markers(test_corrected_BM27), cell_col = celltype_col)
cowplot::plot_grid(emd_BM27$violin, emd_BM27$scatterplot)
ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BM27/emd_correction.pdf")

emd_WF26 <- test_uncorrected_WF26 %>%
  evaluate_emd(test_corrected_WF26, markers = get_markers(test_corrected_WF26), cell_col = celltype_col)
cowplot::plot_grid(emd_WF26$violin, emd_WF26$scatterplot)
ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WF26/emd_correction.pdf")

emd_WM34 <- test_uncorrected_WM34 %>%
  evaluate_emd(test_corrected_WM34, markers = get_markers(test_corrected_WM34), cell_col = celltype_col)
cowplot::plot_grid(emd_WM34$violin, emd_WM34$scatterplot)
ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WM34/emd_correction.pdf")

saveRDS(test_corrected_BF21,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BF21/corrected_total_BF21.rds")
saveRDS(test_uncorrected_BF21,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BF21/uncorrected_total_BF21.rds")

saveRDS(test_corrected_BM27,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BM27/corrected_total_BM27.rds")
saveRDS(test_uncorrected_BM27,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BM27/uncorrected_total_BM27.rds")

saveRDS(test_corrected_WM34,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WM34/corrected_total_WM34.rds")
saveRDS(test_uncorrected_WM34,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WM34/uncorrected_total_WM34.rds")

saveRDS(test_corrected_WF26,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WF26/corrected_total_WF26.rds")
saveRDS(test_uncorrected_WF26,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WF26/uncorrected_total_WF26.rds")

saveRDS(test_corrected,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/merged/test_corrected.rds")
saveRDS(test_uncorrected,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/merged/test_uncorrected.rds")


# MAD
mad <- test_uncorrected %>%
  cyCombine::evaluate_mad(test_corrected,
                          markers = get_markers(test_corrected),
                          cell_col = celltype_col,
                          filter_limit = NULL)

cat('The MAD score is: ', mad$score, '\n')

density_grid <- plot_density(test_uncorrected, test_corrected, ncol = 4, xlims = c(-2, 9))

ggsave(
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/density_plot_CITE_flow.pdf",
  plot = density_grid,
  width = 280,  # Set the desired width
  height = 280, # Set the desired height
  limitsize=FALSE
)

#modify function to save as a multi-page pdf
plot_density <- function(uncorrected,
                         corrected,
                         markers = NULL,
                         directory=NULL,
                         filename = NULL,
                         y = "batch",
                         dataset_names = NULL,
                         ncol = 1,
                         format = 1,
                         markers_per_page=1) {
  
  # Check for packages
  cyCombine:::missing_package("ggridges", "CRAN")
  cyCombine:::missing_package("ggplot2", "CRAN")
  cyCombine:::missing_package("cowplot", "CRAN")
  
  # Get markers
  if (is.null(markers)) {
    markers <- overlap_markers
  }
  
  # Select only the overlapping markers and necessary columns from uncorrected and corrected data
  uncorrected <- uncorrected %>%
    dplyr::select(dplyr::all_of(c("batch", markers)))  # Select only the necessary columns
  
  corrected <- corrected %>%
    dplyr::select(dplyr::all_of(c("batch", markers)))  # Select only the necessary columns

  
  # Rename datasets
  if (is.null(dataset_names)) {
    dataset_names <- c('Uncorrected', 'Corrected')
  } else if (length(dataset_names) != 2) {
    dataset_names <- c('Dataset 1', 'Dataset 2')
  }
  
  # Define how densities are stacked
  if (y == 'Type') {
    batch1 <- dataset_names[1]
    batch2 <- dataset_names[2]
  } else {
    batch1 <- as.factor(uncorrected[[y]])
    batch2 <- as.factor(corrected[[y]])
  }
  
  # Combine uncorrected and corrected data
  df <- rbind(
    mutate(uncorrected, Type = dataset_names[1]),
    mutate(corrected, Type = dataset_names[2])
  )
  
  # Create a list to store the plots
  plot_list <- list()
  
  # Loop through markers and create plots
  for (c in 1:length(markers)) {
    marker_df <- df %>%
      dplyr::select(all_of(markers[c]), Type, batch)
    
    p <- ggplot2::ggplot(marker_df, ggplot2::aes_string(x = markers[c], y = y)) +
      ggridges::geom_density_ridges(ggplot2::aes(color = Type, fill = Type), alpha = 0.4) +
      ggplot2::scale_fill_manual(values = c("#FF6666", "#3333CC"), guide = guide_legend(title = "Batch")) +
      ggplot2::scale_color_manual(values = c("#FF6666", "#3333CC"), guide = guide_legend(title = "Batch")) +
      ggplot2::labs(y = y) +
      ggplot2::theme_bw()
    
    plot_list[[c]] <- p
  }
    
  # Save plots to multiple pages...
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / markers_per_page)
  
  for (page in 1:num_pages) {
    start_marker <- (page - 1) * markers_per_page + 1
    end_marker <- min(page * markers_per_page, num_plots)
    current_plots <- plot_list[start_marker:end_marker]
    
    pdf_file <- if (!is.null(filename) && !is.null(directory)) {
      current_markers <- markers[start_marker:end_marker] # Get the current markers for this page
      marker_names <- paste(current_markers, collapse = "_") # Concatenate marker names with underscore
      file.path(directory, sprintf("%s_page%d_%s.pdf", filename, page,marker_names))
    } else {
      NULL
    }
    
    # Save plots to PDF file (one page at a time)
    if (!is.null(pdf_file)) {
      pdf(pdf_file, width = 14, height = 10)
      gridExtra::grid.arrange(grobs = current_plots, ncol = ncol)
      dev.off()
    }
  }

}


plots_total4 <- plot_density(test_uncorrected, test_corrected, filename = "donors_CD34",
                                        directory = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/merged", markers_per_page = 1, ncol = 1)


plots_BF21 <- plot_density(test_uncorrected_BF21, test_corrected_BF21, filename = "BF21",
                          directory = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BF21", markers_per_page = 1, ncol = 1)

plots_BM27 <- plot_density(test_uncorrected_BM27, test_corrected_BM27, filename = "BM27",
                      directory = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/BM27", markers_per_page = 1, ncol = 1)

plots_WF26 <- plot_density(test_uncorrected_WF26, test_corrected_WF26, filename = "WF26",
                      directory = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WF26", markers_per_page = 1, ncol = 1)

plots_WM34 <- plot_density(test_uncorrected_WM34, test_corrected_WM34, filename = "WM34",
                      directory = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/WM34", markers_per_page = 1, ncol = 1)

test_uncorrected_BF21 <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/BF21/uncorrected_total_BF21.rds")
test_corrected_BF21 <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/BF21/corrected_total_BF21.rds")
test_uncorrected_BM27 <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/BM27/uncorrected_total_BM27.rds")
test_corrected_BM27 <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/BM27/corrected_total_BM27.rds")
test_uncorrected_WF26 <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/WF26/uncorrected_total_WF26.rds")
test_corrected_WF26 <- readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/WF26/corrected_total_WF26.rds")
test_uncorrected_WM34 <-readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/WM34/uncorrected_total_WM34.rds")
test_corrected_WM34 <-readRDS("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_new_WNN/WM34/corrected_total_WM34.rds")


uncorrected_sliced <- test_uncorrected %>%
  dplyr::group_by(batch) %>%
  dplyr::slice_sample(n = 40000) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(id)

corrected_sliced <- test_corrected %>%
  dplyr::filter(id %in% uncorrected_sliced$id)

uncorrected_sliced_BF21 <- test_uncorrected_BF21 %>%
  dplyr::group_by(batch) %>%
  dplyr::slice_sample(n = 40000) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(id)

corrected_sliced_BF21 <- test_corrected_BF21 %>%
  dplyr::filter(id %in% uncorrected_sliced_BF21$id)

uncorrected_sliced_BM27 <- test_uncorrected_BM27 %>%
  dplyr::group_by(batch) %>%
  dplyr::slice_sample(n = 40000) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(id)

corrected_sliced_BM27 <- test_corrected_BM27 %>%
  dplyr::filter(id %in% uncorrected_sliced_BM27$id)

uncorrected_sliced_WF26 <- test_uncorrected_WF26 %>%
  dplyr::group_by(batch) %>%
  dplyr::slice_sample(n = 40000) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(id)

corrected_sliced_WF26 <- test_corrected_WF26 %>%
  dplyr::filter(id %in% uncorrected_sliced_WF26$id)

uncorrected_sliced_WM34 <- test_uncorrected_WM34 %>%
  dplyr::group_by(batch) %>%
  dplyr::slice_sample(n = 40000) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(id)

corrected_sliced_WM34 <- test_corrected_WM34 %>%
  dplyr::filter(id %in% uncorrected_sliced_WM34$id)

# UMAP plot uncorrected
umap1 <- uncorrected_sliced %>%
  plot_dimred(name = "uncorrected", type = "umap")

# UMAP plot corrected
umap2 <- corrected_sliced %>%
  plot_dimred(name = "corrected", type = "umap")

# Show plots
cowplot::plot_grid(umap1, umap2)


ggsave(
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_CD34_only/merged/UMAP_4_donors_Cite_flow.pdf",
  width = 16,   # Set the desired width
  height = 8,  # Set the desired height
  units = "in" # Specify the units (inches in this case)
)

umap1_BF21 <- uncorrected_sliced_BF21 %>%
  plot_dimred(name = "uncorrected_BF21", type = "umap")
umap2_BF21 <- corrected_sliced_BF21 %>%
  plot_dimred(name = "corrected_BF21", type = "umap")
UMAP_BF21 <- cowplot::plot_grid(umap1_BF21, umap2_BF21)
ggsave(plot=UMAP_BF21,
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_TNC_only/BF21/UMAP_overlay_Cite_flow.pdf",
  width = 16,   # Set the desired width
  height = 8,  # Set the desired height
  units = "in" # Specify the units (inches in this case)
)

umap1_BM27 <- uncorrected_sliced_BM27 %>%
  plot_dimred(name = "uncorrected_BM27", type = "umap")
umap2_BM27 <- corrected_sliced_BM27 %>%
  plot_dimred(name = "corrected_BM27", type = "umap")
UMAP_BM27 <- cowplot::plot_grid(umap1_BM27, umap2_BM27)
ggsave(plot=UMAP_BM27,
       filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_TNC_only/BM27/UMAP_overlay_Cite_flow.pdf",
       width = 16,   # Set the desired width
       height = 8,  # Set the desired height
       units = "in" # Specify the units (inches in this case)
)

umap1_WF26 <- uncorrected_sliced_WF26 %>%
  plot_dimred(name = "uncorrected_WF26", type = "umap")
umap2_WF26 <- corrected_sliced_WF26 %>%
  plot_dimred(name = "corrected_WF26", type = "umap")
UMAP_WF26 <- cowplot::plot_grid(umap1_WF26, umap2_WF26)
ggsave(plot=UMAP_WF26,
       filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_TNC_only/WF26/UMAP_overlay_Cite_flow.pdf",
       width = 16,   # Set the desired width
       height = 8,  # Set the desired height
       units = "in" # Specify the units (inches in this case)
)

umap1_WM34 <- uncorrected_sliced_WM34 %>%
  plot_dimred(name = "uncorrected_WM34", type = "umap")
umap2_WM34 <- corrected_sliced_WM34 %>%
  plot_dimred(name = "corrected_WM34", type = "umap")
UMAP_WM34 <- cowplot::plot_grid(umap1_WM34, umap2_WM34)
ggsave(plot=UMAP_WM34,
       filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_TNC_only/WM34/UMAP_overlay_Cite_flow.pdf",
       width = 16,   # Set the desired width
       height = 8,  # Set the desired height
       units = "in" # Specify the units (inches in this case)
)


saveRDS(test_corrected,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/cite_flow_corrected.rds")
saveRDS(test_uncorrected,"/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/cite_flow_uncorrected.rds")

uncor_df <- cbind.data.frame(umap1$data[,1:2], uncorrected_sliced)
cor_df <- cbind.data.frame(umap2$data[,1:2], corrected_sliced)

# Loop over markers to make the plots for both the uncorrected and corrected data
expr_plots <- list()
for (m in get_markers(corrected_sliced)) {
  
  # Find range for marker to make plots comparable
  range <- summary(c(uncor_df[,m], cor_df[,m]))[c(1,6)]
  
  expr_plots[[paste0(m, '1')]] <- ggplot(uncor_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes_string(color = m), alpha = 0.3, size = 0.4) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle('UMAP - uncorrected') + 
    scale_color_gradientn(m, limits  = range, colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(50)) +
    facet_wrap(~batch)
  
  expr_plots[[paste0(m, '2')]] <- ggplot(cor_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes_string(color = m), alpha = 0.3, size = 0.4) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle('UMAP - corrected') + 
    scale_color_gradientn(m, limits = range, colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(50)) +
    facet_wrap(~batch)
  
}

# Show plots
UMAP_plots_CITE_Flow <- cowplot::plot_grid(plotlist = expr_plots, ncol = 2)

pdf("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/UMAP_plots_CITE_Flow_new.pdf")

for (marker_name in names(expr_plots)) {
  # Start a new page for each plot
  if (marker_name != names(expr_plots)[1]) {
    # Add a blank page
    plot.new()
    # Set the layout to a single plot occupying the entire page
    layout(matrix(1))
  }
  
  # Save each plot to a PDF file
  ggsave(expr_plots[[marker_name]], filename = paste0("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/UMAP_", marker_name, ".pdf"))
}


# Close the PDF file
dev.off()


#
#
#Pre-processing the spectral flow cytometry data for integration and correction of just flow

#due to issues with marker names, ran CyCombine functions separately
#key is to preserve the marker names in the original FCS files, otherwise it will trigger warning and issues with asinh transform

flowset <- compile_fcs(data_dir)

raw_data <- convert_flowset(
  flowset,
  metadata = paste0(data_dir, "/metadata_9donors.csv"), 
  filename_col = "Filename",
  batch_ids = "Batch",
  sample_ids = "Patient_id",
  clean_colnames = FALSE,
  down_sample = FALSE,
)

#sfc_markers %in% colnames(fcs_data)

fcs_markers <- get_markers(raw_data)

uncorrected <- transform_asinh(raw_data, markers = fcs_markers,cofactor=6000, derand = FALSE)

# Clustering with kohonen 6x6, since we are only interested in relatively high-level labels
#optional
som_ <- uncorrected %>%
  dplyr::select(dplyr::all_of(fcs_markers)) %>%
  as.matrix() %>%
  kohonen::som(grid = kohonen::somgrid(xdim = 6, ydim = 6),
               rlen = 10,
               dist.fcts = "euclidean")


cell_clustering_som <- som_$unit.classif
codes <- som_$codes[[1]]

# Meta-clustering with ConsensusClusterPlus
mc <- ConsensusClusterPlus::ConsensusClusterPlus(t(codes), maxK = 25, reps = 100,
                                                 pItem = 0.9, pFeature = 1, plot = F,
                                                 clusterAlg = "hc", innerLinkage = "average", 
                                                 finalLinkage = "average",
                                                 distance = "euclidean", seed = seed)

saveRDS(uncorrected, file = file.path(data_dir, "uncorrected.RDS"))



# Get cluster ids for each cell - here we choose to look at 20 meta-clusters but we will merge some in the next step
code_clustering1 <- mc[[20]]$consensusClass
uncorrected$label <- as.factor(code_clustering1[cell_clustering_som])

# Down-sampling for plot
set.seed(seed)
spectral_sliced <- uncorrected%>%
  dplyr::slice_sample(n = 50000)

# Let us visualize these clusters on a UMAP
umap <- spectral_sliced %>% plot_dimred(name = 'SFC clusters')
spectral_sliced <- cbind.data.frame(umap$data, spectral_sliced)

color_clusters <- c(RColorBrewer::brewer.pal(12, 'Paired'), RColorBrewer::brewer.pal(8, 'Dark2'))

ggplot(spectral_sliced, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = label), alpha = 0.3, size = 0.4) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = color_clusters) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))

ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/flow_uncorrected_7_donors.pdf")

#We see significant batch effect between the donors, so now we are going to do batch correction on just the flow data to see how we can fix it

#ignore for now, for correcting with CITE-seq
corrected <- batch_correct(
  df = uncorrected,
  covar = "Set",
  markers = fcs_markers,
  norm_method = "rank",
  rlen = 100, #higher values are recommended if 10 does not appear to perform well
  seed = 840 #Recommended to use your own random seed
)

#batch correct for flow data
seed = 840 #same seed asset in the beginning
corrected <- uncorrected %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

#saveRDS(uncorrected, file = file.path(data_dir, "corrected_progress.RDS")

labels <- corrected %>%
  create_som(seed = seed,
            xdim = 8,
            ydim = 8,
            markers = get_markers(corrected))
        
corrected <- corrected %>%
  dplyr::mutate(som = labels)

celltype_col <- "som"

uncorrected <- corrected %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(uncorrected, by = "id")

saveRDS(corrected, file = file.path(data_dir, "corrected_9donors.RDS"))
saveRDS(uncorrected, file = file.path(data_dir, "uncorrected_9donors.RDS"))

# Evaluate Earth Movers Distance

emd <- uncorrected %>%
  evaluate_emd(corrected, markers = get_markers(corrected), cell_col = celltype_col)

# Show plots
cowplot::plot_grid(emd$violin, emd$scatterplot)
ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/flow_corrected_9_donors.pdf")


# MAD
mad <- uncorrected %>%
  cyCombine::evaluate_mad(corrected,
                          markers = get_markers(corrected),
                          cell_col = celltype_col,
                          filter_limit = NULL)

#However, as one should always do we will also visualize the correction with plots. First, the marker distributions before and after:

# Capture the output of plot_density() in a variable
density_grid <- (plot_density(uncorrected, corrected, ncol = 4, xlims = c(-2, 9)))


ggsave(
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/density_plot_un_vs_corrected_9donors.pdf",
  plot = density_grid,
  width = 200,  # Set the desired width
  height = 200, # Set the desired height
  limitsize=FALSE
  )


plot_density_modified <- function(...) {
  density_plots <- plot_density(...)
  
  # Open PDF device
  pdf("/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/density_plot_un_vs_corrected_9donors_new.pdf")
  
  # Loop through the density plots and save each one to the PDF device
  for (i in seq_along(density_plots)) {
    # Start a new page
    if (i > 1) {
      # Add a new page to the PDF device
      grDevices::newpage()
    }
    
    # Save the current plot to the PDF device
    print(density_plots[[i]])
  }
  
  # Close the PDF device
  dev.off()
  
  # Return the density plots
  density_plots
}

# Use the modified plot_density function
density_grid <- plot_density_modified(uncorrected, corrected, ncol = 4, xlims = c(-2, 9))

# Access the individual density plots in density_grid if needed


#convert corrected df to FCS file again; first make sure you read in the panel, without the panel, we cannot create the fcs file later
panel_info <- read_csv(file.path(data_dir, "panel.csv"))

corrected_sce <- df2SCE(
  corrected,
  non_markers = NULL,
  sample_col = "sample",
  panel = panel_info,
  panel_channel = "Channel",
  panel_antigen = "Marker",
  panel_type = "Type",
  transform_cofactor = 5
)

saveRDS(corrected_sce, file = file.path(data_dir, "corrected_sce_9donors.RDS"))

sce2FCS(
  corrected_sce,
  outdir = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/Corrected_FCS",
  assay = "counts",
  keep_dr = TRUE,
  keep_cd = TRUE,
  randomize = TRUE
)

# UMAP representations
#Down-sampling for plot (10,000 cells from each platform)
#set.seed(seed)
uncorrected_sliced <- uncorrected %>%
  dplyr::group_by(batch) %>%
  dplyr::slice_sample(n = 10000) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(id)

corrected_sliced <- corrected %>%
  dplyr::filter(id %in% uncorrected_sliced$id)


# UMAP plot uncorrected
umap1 <- uncorrected_sliced %>%
  plot_dimred(name = "uncorrected", type = "umap")
#ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/UMAP_uncorrected.pdf")

# UMAP plot corrected
umap2 <- corrected_sliced %>%
  plot_dimred(name = "corrected", type = "umap")
#ggsave("/data/salomonis2/LabFiles/Annie/CyCombine/UMAP_corrected.pdf")

# Show overlayed plots

cowplot::plot_grid(umap1, umap2)

# Adjust the size of the output device
ggsave(
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/UMAP_overlay_9donors.pdf",
  width = 12,   # Set the desired width
  height = 6,  # Set the desired height
  units = "in" # Specify the units (inches in this case)
)

#if assign names manually
spectral$label <- spectral_cl_names[spectral$label]

#optional: pull up violin plots to annotate some clusters
spectral_sliced <- flow_obj_9donors %>%
  dplyr::slice_sample(n = 20000)
umap <- spectral_sliced %>% plot_dimred(name = 'SFC clusters')
spectral_sliced <- cbind.data.frame(umap$data, spectral_sliced)

ggplot(spectral_sliced, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = factor(label)), alpha = 0.3, size = 0.4) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = expanded_colors) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))

ggsave(
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/umap_9donors_clusters.pdf",
  width = 18,   # Set the desired width
  height = 12,  # Set the desired height
  units = "in" # Specify the units (inches in this case)
)

# And now let us use violin plots for the assignment of labels for each population
p <- list()
for (m in new_markers) {
  p[[m]] <- ggplot(flow_obj_9donors, aes_string(x = 'label', y = m)) +
    geom_violin(aes(fill = factor(label))) +
    scale_fill_manual(values = expanded_colors) +
    theme_bw()
}

pdf("/data/salomonis2/LabFiles/Annie/CyCombine/test_violin_plots.pdf", onefile = TRUE)

for (i in seq(length(p))) {
  grob <- ggplotGrob(p[[i]])
  grid.arrange(grob)
}

dev.off()

p2 <- ggpubr::ggarrange(plotlist = p, common.legend = T)

ggsave(
  p2,
  filename = "/data/salomonis2/LabFiles/Annie/CyCombine/violinPlots_9donors_clusters.pdf",
  width = 320,   # Set the desired width
  height = 320,  # Set the desired height
  limitsize=FALSE
)



colnames(flow_obj_9donors)[colnames(flow_obj_9donors) == "HLA-DR"] <- "HLA_DR"
colnames(flow_obj_9donors)[colnames(flow_obj_9donors) == "InfinityMarker_TIM-4"] <- "InfinityMarker_TIM_4"
colnames(flow_obj_9donors)[colnames(flow_obj_9donors) == "InfinityMarker_HLA-DR"] <- "InfinityMarker_HLA_DR"

new_markers <- get_markers(flow_obj_9donors)


# Initialize a list to store the UMAP plots
umap_plots <- list()

# Iterate over the markers
for (marker in new_markers) {
  
  # Generate the UMAP plot colored by the marker
  umap_plot <- ggplot(spectral_sliced, aes(x = UMAP1, y = UMAP2, color = .data[[marker]])) +
    geom_point(alpha = 0.3, size = 0.4) +
    theme_bw() +
    ggtitle(paste("UMAP -", marker)) +
    scale_color_gradientn(limits = range(spectral_sliced[[marker]]),
                          colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 150, name = "Spectral")))(133)) +
    facet_wrap(~batch)
  
  # Store the UMAP plot in the list
  umap_plots[[marker]] <- umap_plot
}

# Display the UMAP plots
umap_markers <- cowplot::plot_grid(plotlist = umap_plots, ncol = 4)


pdf("/data/salomonis2/LabFiles/Annie/CyCombine/test_plots.pdf", onefile = TRUE)

for (i in seq(length(umap_plots))) {
  grob <- ggplotGrob(umap_plots[[i]])
  grid.arrange(grob)
}

dev.off()

colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD133_S16016B"] <- "Hu.CD133"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD14_M5E2"] <- "Hu.CD14"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD15_W6D3"] <- "Hu.CD15"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD226_TX25"] <- "Hu.CD226"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD305_LAIR1"] <- "Hu.CD305"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD38_HIT2"] <- "Hu.CD38"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD45_2D1"] <- "Hu.CD45"
colnames(BF21_tibble)[colnames(BF21_tibble) == "Hu.CD4_RPA.T4"] <- "Hu.CD4"

colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD105_43A3"] <- "Hu.CD105"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD138_DL.101"] <- "Hu.CD138"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD133_S16016B"] <- "Hu.CD133"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD14_M5E2"] <- "Hu.CD14"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD15_W6D3"] <- "Hu.CD15"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD226_TX25"] <- "Hu.CD226"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD305_LAIR1"] <- "Hu.CD305"

colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD38_HIT2"] <- "Hu.CD38"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD45_2D1"] <- "Hu.CD45"
colnames(BM27_tibble)[colnames(BM27_tibble) == "Hu.CD4_RPA.T4"] <- "Hu.CD4"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD105_43A3"] <- "Hu.CD105"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD138_DL.101"] <- "Hu.CD138"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD133_S16016B"] <- "Hu.CD133"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD14_M5E2"] <- "Hu.CD14"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD15_W6D3"] <- "Hu.CD15"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD226_TX25"] <- "Hu.CD226"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD305_LAIR1"] <- "Hu.CD305"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD38_HIT2"] <- "Hu.CD38"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD45_2D1"] <- "Hu.CD45"
colnames(WF26_tibble)[colnames(WF26_tibble) == "Hu.CD4_RPA.T4"] <- "Hu.CD4"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD105_43A3"] <- "Hu.CD105"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD138_DL.101"] <- "Hu.CD138"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD133_S16016B"] <- "Hu.CD133"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD14_M5E2"] <- "Hu.CD14"

colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD15_W6D3"] <- "Hu.CD15"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD226_TX25"] <- "Hu.CD226"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD305_LAIR1"] <- "Hu.CD305"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD38_HIT2"] <- "Hu.CD38"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD45_2D1"] <- "Hu.CD45"
colnames(WM34_tibble)[colnames(WM34_tibble) == "Hu.CD4_RPA.T4"] <- "Hu.CD4"

colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(BF21_032123_cite)[colnames(BF21_032123_cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"

colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(BM27_120522_cite)[colnames(BM27_120522_cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(WF26_031423_cite)[colnames(WF26_031423_cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(WM34_120522_cite)[colnames(WM34_120522_cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(test_uncorrected_4donors)[colnames(test_uncorrected_4donors) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(test_corrected_4donors)[colnames(test_corrected_4donors) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"


colnames(uncorrected_9donors)[colnames(uncorrected_9donors) == "HLA-DR"] <- "HLA_DR"
colnames(uncorrected_9donors)[colnames(uncorrected_9donors) == "InfinityMarker_FR-b"] <- "InfinityMarker_FR_b"
colnames(uncorrected_9donors)[colnames(uncorrected_9donors) == "InfinityMarker_HLA-A-B-C"] <- "InfinityMarker_HLA_A_B_C"
colnames(uncorrected_9donors)[colnames(uncorrected_9donors) == "InfinityMarker_HLA-DR-DP-DQ"] <- "InfinityMarker_HLA_DR_DP_DQ"
colnames(uncorrected_9donors)[colnames(uncorrected_9donors) == "InfinityMarker_IgG-FC"] <- "InfinityMarker_IgG_FC"
colnames(uncorrected_9donors)[colnames(uncorrected_9donors) == "InfinityMarker_TIM-4"] <- "InfinityMarker_TIM_4"

colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(BF21_CD34_Cite)[colnames(BF21_CD34_Cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(BM27_CD34_Cite)[colnames(BM27_CD34_Cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(WF26_CD34_Cite)[colnames(WF26_CD34_Cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"

colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD105-43A3"] <- "Hu.CD105"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD138-DL.101"] <- "Hu.CD138"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD133-S16016B"] <- "Hu.CD133"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD14-M5E2"] <- "Hu.CD14"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD15-W6D3"] <- "Hu.CD15"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD226-TX25"] <- "Hu.CD226"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD305-LAIR1"] <- "Hu.CD305"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD38-HIT2"] <- "Hu.CD38"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD45-2D1"] <- "Hu.CD45"
colnames(WM34_CD34_Cite)[colnames(WM34_CD34_Cite) == "Hu.CD4-RPA.T4"] <- "Hu.CD4"



corrected_total <- combined_data %>%
  select(-dataset) %>%
  select(-cell_barcode) %>%
  batch_correct(seed = seed,
                xdim = 8,
                ydim = 8,
                norm_method = 'rank',
                ties.method = 'average')

#08-23-23 output UMAP coordinates
umap_corrected <- corrected_sliced %>%
  plot_dimred(name = "corrected", type = "umap")

# Extract UMAP coordinates from the umap_corrected object
umap_coords <- umap_corrected$data

# Add UMAP coordinates to the corrected_sliced object
corrected_sliced$umap_x <- umap_coords[, 1]$UMAP1
corrected_sliced$umap_y <- umap_coords[, 2]$UMAP2

#take batches we want
batches_to_export <- c("BF21_032123_cite_data","BM27_120522_cite_data","WF26_031423_cite_data","WM34_120522_cite_data")
selected_batches_data <- corrected_sliced[corrected_sliced$batch %in% batches_to_export, ]
selected_attributes <-c("label","sample","cell_barcode","umap_x","umap_y")
# Create a matrix containing the selected attributes
attribute_matrix <- selected_batches_data[selected_attributes]

# Define the file path for the output text matrix
output_file <- "/data/salomonis2/LabFiles/Annie/CyCombine/updated_0723/4_donors_updated_082323/output_matrix.txt"

# Write the matrix to a tab-delimited text file
write.table(attribute_matrix, file = output_file, sep = "\t", row.names = FALSE)
