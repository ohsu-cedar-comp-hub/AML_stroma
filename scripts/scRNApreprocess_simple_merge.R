
library(Seurat)
library(data.table)
library(knitr)
library(ggplot2)
library(DT)
library(grid)
library(gridExtra)
library(dplyr)
library(future)
options(future.globals.maxSize = 60000 * 1024^2)



# Load library to allow parallelization
# Note: options are multiprocess, sequential
run2process <- as.character(snakemake@input)
seurat_out <- as.character(snakemake@output[["seurat_obj"]])
out_dir <- as.character(snakemake@params[["directory"]])

n_Dims <- as.numeric(snakemake@params[["n_Dims"]])
n_cores <- as.numeric(snakemake@params[["n_cores"]])
run2process <- read.delim(run2process, header = T)
integrate_samples <- ifelse(dim(run2process)[1] > 1,TRUE,FALSE)

set.seed(as.numeric(snakemake@params[['seed']]))
seuratObj_master <- NULL
min.cells <- 3
min.features <- 200
input.data <- NULL


run2process <- run2process[order(run2process[,2]),] #order by sampleID

if(!(file.exists(out_dir))) {
  dir.create(out_dir,showWarnings = FALSE, recursive = TRUE)
}

for (i in 1:dim(run2process)[1]){
	
	input.data <- Read10X(data.dir = as.character(run2process[i,1]))
	# Initialize the Seurat object with the raw (non-normalized data).
	new_seurat_Obj <- CreateSeuratObject(counts = input.data, project = as.character(run2process[i,2]), min.cells = min.cells, min.features = min.features)
	
	#merge expression count data with the same sample ID
	if (i > 1){
		seuratObj_master <- merge(seuratObj_master,new_seurat_Obj, project = as.character(run2process[i,2]))
		
	}else{
		seuratObj_master <- new_seurat_Obj
	}

}

#get mito
seuratObj_master[["percent.mt"]] <- PercentageFeatureSet(object = seuratObj_master, pattern = "^MT-")


### Calculate SCTransform
#this is used to calculate the variable features to be used for intgration

seuratObj_master <- SCTransform(seuratObj_master, verbose = FALSE, return.only.var.genes = FALSE)


plan("multiprocess", workers = n_cores)
seuratObj_master <- NormalizeData(seuratObj_master, verbose = FALSE, assay = "RNA")
seuratObj_master <- FindVariableFeatures(seuratObj_master, selection.method = "vst", nfeatures = 2000, verbose = FALSE, assay = "RNA")
seuratObj_master <- ScaleData(seuratObj_master, verbose = FALSE, assay = "RNA")
seuratObj_master <- RunPCA(seuratObj_master, npcs = n_Dims, verbose = FALSE, assay = "RNA")
seuratObj_master <- RunUMAP(seuratObj_master, reduction = "pca", dims = 1:n_Dims, assay = "RNA")



saveRDS(object = seuratObj_master, file = seurat_out)






