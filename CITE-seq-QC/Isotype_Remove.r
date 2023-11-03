suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
  make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  folder containing ADT mtx..."),
  make_option(c("-q", "--quantile"), type="double", default=0.9, help="percentile cutoff of a isotype being positive ", metavar="number"),
  make_option(c("-n", "--numPos"), type="integer", default=1, help="cells with more than this number of positive isotype are regarded as aggregation ", metavar="number"),
  make_option(c("-o", "--outname"), action = "store", default = "output", type = "character", help="Output file(s) base name .")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p, arguments$q , arguments$n , arguments$o)
message(" Loading the required packages .... ")
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DropletUtils))

message("creating function  in the  environment to load mtx file ...")
loadMex <- function(matrix_dir, barcodetsv="barcodes.tsv.gz", 
                    featuretsv="features.tsv.gz", matmtx = "matrix.mtx.gz", removeLane = FALSE) {
  library(Matrix)
  barcode.path <- paste0(matrix_dir,'/', barcodetsv)
  
  if(file.access(paste0(matrix_dir,'/features.tsv.gz'), mode=4) > -1  ){
    features.path <- paste0(matrix_dir,'/features.tsv.gz')
  }else if(file.access(paste0(matrix_dir,'/genes.tsv.gz'), mode=4) > -1){
    features.path <- paste0(matrix_dir,'/genes.tsv.gz')
  }else{
    features.path <- paste0(matrix_dir,'/', featuretsv)
  }
  matrix.path <- paste0(matrix_dir,'/', matmtx)
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  if(removeLane){
    colnames(mat) = substr(barcode.names$V1, 1, 16)
  }else{
    colnames(mat) = barcode.names$V1
  }
  
  rownames(mat) = feature.names$V1
  
  mat
  
}
message("reading mtx file ...")
mtx <- loadMex(arguments$p)
print(paste("sample has" ,nrow(mtx), "of ADTs and" ,ncol(mtx), "of cellbarcodes"))
message("loading and creating ADT-Gene Rank  variable ...")
genADTRank <- function(df_raw_FeatByCell){
  mydf = df_raw_FeatByCell
  
  res = data.frame(t(apply(mydf,1, quantile, c(0.25,0.5,0.75,1))))
  
  res$Sum = rowSums(mydf)
  sumTotal = sum(res$Sum)
  res$PercentCount = formatC(res$Sum*100/sumTotal, format = 'f', digits = 2)
  
  
  res = res[order(res$Sum, decreasing = T),]
  colnames(res)[1:4] = c("25pct","Median","75pct","Max")
  
  res
}
rank=  genADTRank(mtx); UID = rownames(rank)
rank = cbind(UID, rank)
write.table(rank, file = paste(arguments$o, "_geneADT-Rank.txt", sep = ""),  row.names = F, col.names = T, sep ="\t", quote = F)
#isotype aggregation removal
## detect protein aggregates based on isotype counts 
## parameters
##   df_isotype_by_cell: a data frame or matrix with ONLY isotype info, antibody in row and cells in column
##   cutoff_quantile:    percentile cutoff of a isotype being positive
##   cutoff_num_pos:     cells with more than this number of positive isotype are regarded as aggregation  

is_isotype_agg = function(df_isotype_by_cell, cutoff_quantile = arguments$q, cutoff_num_pos=arguments$n){
  table_positive  = apply(df_isotype_by_cell, 1, function(x){
    x > quantile(x, cutoff_quantile)
  })
  rowSums(table_positive) > cutoff_num_pos
}
#### few steps for removing isotypes .. --- ###
iso.mtx <- mtx[grep(pattern = "isotype|Isotype", x = rownames(mtx)),]
print(paste("the iso.mtx consist of", nrow(iso.mtx), "isotypes"))
is_isotype_agg_data = is_isotype_agg(iso.mtx) ##### feed only isotype by cell matrix

message("removing  isotype control from mtx ....")
mex_clean <- mtx[,which(!is_isotype_agg_data )] 
print(paste("sample has" ,nrow(mex_clean), "of ADTs and" ,ncol(mex_clean), "of cellbarcodes"))
print(paste(ncol(mtx) - ncol(mex_clean) , " cell barcodes got removed ..."))

message("Remove cells by aggregates considering all ADTs, greater than 1% total reads ...")
is_agg_mtx_ADT<- function(df){
  table_positive  = apply(df, 1, function(x){
    x > 0.1 * rowSums(df)
  })
  table_positive
  #now ADT = column, cell = row
} 

mex_clean2 <- is_agg_mtx_ADT(mex_clean)

all_cells = rowSums(mex_clean2)
#bad_cells = all_cells[which(all_cells > 0)]
good_cells = data.frame(all_cells[which(all_cells < 1)])
barcodes = rownames(good_cells)
output = mex_clean[, barcodes]  ### this is final // run sc 
message(paste("final mtx consist of", ncol(output), "of cell barcodes .."))

genADTRank = function(df_raw_FeatByCell){
  mydf = df_raw_FeatByCell
  
  res = data.frame(t(apply(mydf,1, quantile, c(0.25,0.5,0.75,1))))
  
  res$Sum = rowSums(mydf)
  sumTotal = sum(res$Sum)
  res$PercentCount = formatC(res$Sum*100/sumTotal, format = 'f', digits = 2)
  
  
  res = res[order(res$Sum, decreasing = T),]
  colnames(res)[1:4] = c("25pct","Median","75pct","Max")
  
  res
}

output_rank = genADTRank(output); UID= rownames(output_rank)
output_rank  =  cbind(UID, output_rank)
fwrite(output_rank, file = paste(arguments$o, "_FinalgeneADT-Rank.txt", sep = ""),  row.names = F, col.names = T, sep ="\t", quote = F)
message("writting new mtx file  in Isotypes-Removed-ADT folder")
write10xCounts(x = output, path = "Isotypes-Removed-ADT", version='3')
message("writting the text file ....")
output  =  as.data.frame(as.matrix(output))
#colnames(output) = paste(arguments$o, colnames(output) , sep = "_" )
output = rownames_to_column(output, var = "UID")
fwrite(output, file = paste(arguments$o, "ADT_clean.txt", sep = "_"),  row.names = F, col.names = T, sep ="\t", quote = F)
sessionInfo()
quit()


