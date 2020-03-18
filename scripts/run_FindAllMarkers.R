# Integrated analysis
options(future.globals.maxSize =90e+09)


resolution = as.double(snakemake@params[["resolution"]])
########################################## finish it
set.seed(as.numeric(snakemake@params[["seed"]])) #seed arg
out_dir <- paste(getwd(),as.character(snakemake@params[["directory"]]),sep="/")
input_rds <- paste(getwd(),snakemake@input[["seuratObj"]],sep="/")
max_nFeature <- as.numeric(snakemake@params[["max_nFeature"]])
min_nFeature <- as.numeric(snakemake@params[["min_nFeature"]])
max_mito <- as.numeric(snakemake@params[["max_mito"]])
n_Dims <- as.numeric(snakemake@params[["n_Dims"]])
n_cores <- as.numeric(snakemake@params[["n_cores"]])
seurat_out <- paste(getwd(),snakemake@output[["seuratObj"]],sep="/")
RNA_out <- paste(getwd(),snakemake@output[["RNA"]],sep="/")
integrated_out <- paste(getwd(),snakemake@output[["Integrated"]],sep="/")

library(Seurat)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)
library(future)
library(reshape2)

#patient_path <- strsplit(getwd(),"/")
#patient_path <- paste(patient_path[[1]][1:length(patient_path[[1]])-1],collapse = "/")
#file_path <- paste(patient_path,'integrated_analysis',sep = "/")
#dir.create(file_path, showWarnings = FALSE)

print(out_dir)
print(input_rds)
print(seurat_out)
print(RNA_out)
print(integrated_out)

if(!(file.exists(out_dir))) {
  dir.create(out_dir,showWarnings = FALSE, recursive = TRUE)
}
setwd(out_dir)
file_path <- out_dir
################################# functions ################################



simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


################################################## analyis start ############################



integrated = readRDS(file = input_rds)
#


DefaultAssay(object = integrated) <- "RNA"


integrated <- subset(integrated, subset = nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature & percent.mt < max_mito )
DefaultAssay(object = integrated) <- "integrated"

#integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)
options(future.globals.maxSize =10e+09)


# UMAP and Clustering
Ndims <- n_Dims

integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:Ndims)
integrated <- FindClusters(integrated, resolution = resolution)
integrated <- RunUMAP(integrated, reduction.name = "umap3d", n.components = 3,reduction = "pca" ,dims = 1:Ndims)

integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:Ndims, reduction.name = "umap")



plots <- DimPlot(integrated, group.by = "orig.ident")
#plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, byrow = TRUE, override.aes = list(size = 2.5))))
ggsave(sprintf("%s_dimPlot.pdf",integrated@project.name), plot = plots, path = file_path, device = "pdf")

######################################  Cell Cycle ###################################################
s.genes <- cc.genes$s.genes
#s.genes = unname(sapply(sapply(s.genes, tolower),simpleCap))#for mouse (or non-human)  
g2m.genes <- cc.genes$g2m.genes
#g2m.genes = unname(sapply(sapply(g2m.genes, tolower),simpleCap))


integrated <- CellCycleScoring(integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrated <- RunPCA(integrated, features = c(s.genes, g2m.genes))

pdf(paste(resolution,sep='_','integrated_cell_phase_pca.pdf'))
PCAPlot(integrated)
dev.off()

pdf(paste(resolution,sep='_','integrated_cell_phase_umap_test.pdf'))
plots <- DimPlot(integrated, group.by ='Phase',split.by= 'condition') + facet_wrap(c("condition","Phase"), nrow =2)
plots
dev.off()


integrated$CC.Difference <- integrated$S.Score - integrated$G2M.Score
integrated <- ScaleData(integrated, vars.to.regress = "CC.Difference", features = rownames(integrated))
integrated <- RunPCA(integrated, features = VariableFeatures(integrated), nfeatures.print = 10)
integrated <- RunPCA(integrated, features = c(s.genes, g2m.genes))
integrated <- FindClusters(integrated, resolution = resolution)

saveRDS(integrated, paste(resolution, sep='_','integrated__cc_normalized__seurat_obj.rds'))

pdf(paste(resolution,sep='_','integrated_cell_phase_post_norm_pca.pdf'))
PCAPlot(integrated)
dev.off()

pdf(paste(resolution,sep='_','integrated_cell_phase_post_norm_umap.pdf'))
DimPlot(integrated, group.by = 'Phase')
dev.off()
########################################################  Find All Markers Integrated ###############################



Idents(integrated) <- paste("integrated_snn_res",resolution,sep=".")
plan("multiprocess", workers = n_cores)
cluster_markers <- FindAllMarkers(integrated, only.pos = T, logfc.threshold=0.2, assay = "integrated") #using integrated
saveRDS(cluster_markers, paste(file_path,paste(resolution,sep='_','integrated_cluster_markers.rds'),sep="/"))
write.table(cluster_markers, file = integrated_out, sep = "\t", row.names=T, quote=F)


##############################################  Log Normalize ############################################################
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)
integrated <- FindVariableFeatures(integrated, selection = "vst",nfeatures = length(integrated@assays$integrated@var.features))
integrated <- ScaleData(integrated, features = rownames(integrated))


############################################## Find All Markers RNA ###############################
Idents(integrated) <- paste("integrated_snn_res",resolution,sep=".")
DefaultAssay(integrated) <- "RNA"
plan("multiprocess", workers = n_cores)
cluster_markers <- FindAllMarkers(integrated, only.pos = T, logfc.threshold=0.2) #using RNA
saveRDS(cluster_markers, paste(file_path,paste(resolution,sep='_','RNA_cluster_markers.rds'),sep="/"))
write.table(cluster_markers, file = RNA_out, sep = "\t", row.names=T, quote=F)



saveRDS(integrated, seurat_out)

print("COMPLETED")