options(future.globals.maxSize =90e+09)
library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)
library(future)
library(reshape2)
set.seed(as.numeric(snakemake@params[["seed"]])) #seed arg

output_file <- paste(getwd(),snakemake@output[["marker_file"]],sep="/")
#input_rds <- paste(getwd(),as.character(snakemake@input[["seuratObj"]]),sep="/")
n_cores <- as.numeric(snakemake@params[["n_cores"]])
optimal_res_file <- paste(getwd(),snakemake@input[["optimal_res_file"]],sep="/")
contrast <- snakemake@params[["contrast_DE"]]

type <- snakemake@params[["type"]]

baseline <- contrast[[2]]
contrast <- contrast[[1]]

optimal_res <- as.double(read.table(optimal_res_file)[1,1])

input_dir <- paste(strsplit(optimal_res_file,"/")[[1]][-length(strsplit(optimal_res_file,"/")[[1]])],collapse="/")

input_rds <- grep(paste(sprintf("^%s",optimal_res),"SeuratObj.rds",sep="_"),list.files(input_dir),value = T)[1]

seuratObj = readRDS(paste(input_dir,input_rds,sep="/"))


out_dir <- paste(strsplit(output_file,"/")[[1]][-length(strsplit(output_file,"/")[[1]])],collapse="/")
if(!(file.exists(out_dir))) {
  dir.create(out_dir,showWarnings = FALSE, recursive = TRUE)
}

meta_list <- sapply(colnames(seuratObj@meta.data),function(x){all(contrast %in% unique(seuratObj@meta.data[[x]]))})

if(sum(meta_list) > 1){
	warning("contrast meta data is ambiguous")
}
if(sum(meta_list) < 1){
	stop("contrasts are not in any metadata columns, or there is more than the two contrasts in the column")
}
meta_name <- names(which(meta_list == TRUE))[1]
Idents(seuratObj) <- meta_name
plan("multiprocess", workers = n_cores)
cluster_markers <- FindMarkers(seuratObj, assay = type, ident.1 = contrast, ident.2 = baseline,only.pos = F, logfc.threshold=0.0, min.pct = 0)
write.table(cluster_markers, file=output_file, sep = "\t", row.names=T, quote=F)
