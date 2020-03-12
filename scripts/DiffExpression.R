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

seuratObj@misc[[sprintf("optimal_res_%s",type)]] <- optimal_res

directory_out <- params$destDir

dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)
file_path <- output_dir

Idents(seuratObj) <- meta_name
meta_markers <- FindAllMarkers(seuratObj, assay = type, ident.1 = contrast, ident.2 = baseline,only.pos = F, logfc.threshold=0.0, min.pct = 0)
write.table(meta_markers, file = output_file, sep = "\t", row.names=T, quote=F)

cluster_markers <- FindAllMarkers(seuratObj, assay = type, ident.1 = contrast, ident.2 = baseline,only.pos = F, logfc.threshold=0.0, min.pct = 0)


for (meta in integrated@misc$meta_tab){
	
	for (uni_m in unique(integrated[[meta]])){
		
		
		
		meta_idx <- rownames(integrated@meta.data[integrated@meta.data[[meta]]==as.character(uni_m),])
		meta_UMAP <- DimPlot(integrated[,flow_idx], reduction = "umap", cells = flow_idx, group.by = meta, split.by = "Phase")
		pdf(paste(resolution,sep='_',sprintf('integrated_umap_%s_Phase_%s.pdf',meta,uni_m)),width=15,height=7)
		meta_UMAP
		dev.off()
		
		current_meta <- table(Idents(integrated[,meta_idx]), integrated[,meta_idx]$Phase)
		meta_percent <- current_meta/rowSums(current_meta) *100.0
		melt_meta <- melt(meta_percent,value.name = 'Percentage',varnames = c('Cluster','CC'))


		pdf(paste(resolution,sep='_',sprintf('integrated_%s_stack.pdf',meta)))
		p4 <- ggplot(aes(y = Percentage, x = Cluster, fill = CC), data = melt_meta, cumulative = TRUE) + ggtitle(paste(meta,uni_m,sep = " - ")) +theme_minimal()+geom_col() +
		  geom_text(aes(label = sprintf("%0.2f", round(Percentage, digits = 2))), color = "white", size = 3, position = position_stack(vjust = 0.5))
		p4
		dev.off()

	}
	plan("multiprocess", workers = n_cores)
	
	
	cluster_markers <- FindMarkers(integrated, assay = "RNA", ident.1 = "flow", ident.2 = "cultured",only.pos = F, logfc.threshold=0.0, min.pct = 0) #using RNA assay
	write.table(cluster_markers, file = paste(resolution, sep="_", "RNA_seurat_bulk_markers.tsv"), sep = "\t", row.names=T, quote=F)
}


for (meta in integrated@misc$meta_tab){


	p1 <- DimPlot(integrated, reduction = "umap", group.by = meta, split.by = "Phase") 
	pdf(paste(resolution,sep='_',sprintf('integrated_umap_%s_split_by_phase.pdf',meta)),width=15,height=7)
	p1
	dev.off()

	p4 <- DimPlot(integrated, reduction = "umap", group.by = "Phase", split.by=meta) + facet_wrap(c(meta,"Phase"), nrow =2)
	pdf(paste(resolution,sep='_',sprintf('integrated_umap_%s_Phase_facet_wrap.pdf',meta)),width=15,height=7)
	p4
	dev.off()


	#stromal = c("CD90","CXCL12","FGF2")
	#pdf(paste(resolution,sep='_','integrated_stromal_feature_plot.pdf'), width=15, height=10)
	#p6 = FeaturePlot(integrated,reduction = "umap", features = stromal, ncol=3)
	#p6
	#dev.off()

	p2 <- DimPlot(integrated, reduction = "umap", group.by = meta, label = F)
	pdf(paste(resolution,sep='_',sprintf('integrated_umap_group_by_%s.pdf',meta)))
	p2
	dev.off()


	for (uni_m in unique(integrated[[meta]])){
		meta_idx <- rownames(integrated@meta.data[integrated@meta.data[[meta]]==as.character(uni_m),])
		meta_UMAP <- DimPlot(integrated[,flow_idx], reduction = "umap", cells = flow_idx, group.by = meta, split.by = "Phase")
		pdf(paste(resolution,sep='_',sprintf('integrated_umap_%s_Phase_%s.pdf',meta,uni_m)),width=15,height=7)
		meta_UMAP
		dev.off()
		
		current_meta <- table(Idents(integrated[,meta_idx]), integrated[,meta_idx]$Phase)
		meta_percent <- current_meta/rowSums(current_meta) *100.0
		melt_meta <- melt(meta_percent,value.name = 'Percentage',varnames = c('Cluster','CC'))


		pdf(paste(resolution,sep='_',sprintf('integrated_%s_stack.pdf',meta)))
		p4 <- ggplot(aes(y = Percentage, x = Cluster, fill = CC), data = melt_meta, cumulative = TRUE) + ggtitle(paste(meta,uni_m,sep = " - ")) +theme_minimal()+geom_col() +
		  geom_text(aes(label = sprintf("%0.2f", round(Percentage, digits = 2))), color = "white", size = 3, position = position_stack(vjust = 0.5))
		p4
		dev.off()

	}
}