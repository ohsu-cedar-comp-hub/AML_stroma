options(future.globals.maxSize =90e+09)

set.seed(as.numeric(snakemake@params[["seed"]]) #seed arg

output_file <- snakemake@output[["marker_file"]]
input_rds <- as.character(snakemake@input)
n_cores <- as.numeric(snakemake@params[["n_cores"]])
optimal_res_file <- snakemake@input[["optimal_res_file"]]
contrast <- snakemake@params[["contrast_DE"]]
type <- snakemake@params[["type"]]
meta_name <- names(contrast)
baseline <- contrast[[2]]
contrast <- contrast[[1]]

library(Seurat)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)
library(future)
library(reshape2)

optimal_res <- as.double(read.table(optimal_res_file)[1,1])

filename <- grep(sprintf("%s_SeuratObj.rds",optimal_res),list.files(path = input_dir),value = T)
seuratObj = readRDS(paste(input_dir,filename,sep = "/"))


directory_out <- params$destDir
if(!(file.exists(directory_out))) {
  dir.create(directory_out,showWarnings = FALSE, recursive = TRUE)
}
setwd(output_dir)
file_path <- output_dir

Idents(seuratObj) <- meta_name
plan("multiprocess", workers = n_cores)
cluster_markers <- FindMarkers(seuratObj, assay = type, ident.1 = contrast, ident.2 = baseline,only.pos = F, logfc.threshold=0.0, min.pct = 0)
write.table(cluster_markers, file = output_file, sep = "\t", row.names=T, quote=F)
