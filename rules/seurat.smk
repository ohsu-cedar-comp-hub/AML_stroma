

rule postprocess:
	input:
		naive_seurat = "results/{wave}/preprocessing/simple_seurat_{wave}.rds",
		optimal_res_file="results/{wave}/intermediate/{type}_optimal_resolution.txt", 
		down="results/{wave}/GO_global/{type}/{contrast_DE}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt",
		up="results/{wave}/GO_global/{type}/{contrast_DE}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt",
		bulk_markers="results/{wave}/global_diffexp/{type}_{contrast_DE}_optimal_marker.tsv"
	output:
		html = "results/{wave}/post/{adjp}_{FC}_{contrast_DE}/optimal_{type}_summary.html",
		seuratObj = "results/{wave}/post_{type}/{adjp}_{FC}_{contrast_DE}/optimal_seuratObj.rds"
	params:
		directory = lambda wildcards:"results/{wave}/post_{type}/{adjp}_{FC}_{contrast_DE}".format(wave=wildcards.wave,type=wildcards.type,adjp=wildcards.adjp,FC=wildcards.FC,contrast_DE=wildcards.contrast_DE),
		seed=config["seed"],
		input_directory="results/{wave}/post_processing",
		type = lambda wildcards: "{}".format(wildcards.type)
	conda:
		"../envs/seurat.yaml"
	shell:
		"Rscript rmarkdown::render('../scripts/resolution_post_processing.Rmd', output_file = {output.html}, params=list(optimal_res = {input.optimal_res_file}, destDir = {params.directory}, seed={params.seed}, input_dir={params.input_directory}, GOdown={input.down}, GOup = {input.up}, type={params.type}, markers={input.bulk_markers}))"	
		


rule global_diffexp:
	input:
		optimal_res_file="results/{wave}/intermediate/{type}_optimal_resolution.txt" # extract optimal_res output and get rds from intermediate files
	output:
		marker_file = "results/{wave}/global_diffexp/{type}_{contrast_DE}_optimal_marker.tsv"
	params:
		seed=config["seed"],
		input_directory = "results/{wave}/intermediate",
		out_dir = lambda wildcards: "results/{}/global_diffexp".format(wildcards.wave),
		contrast_DE=get_contrast_global, #named list
		type = lambda wildcards: "{}".format(wildcards.type)
	conda:
		"../envs/seurat.yaml"
	script:
		"../scripts/Global_DE.R"	

		
rule GO_global:
	input:
		degFile="results/{wave}/global_diffexp/{type}_{contrast_DE}_optimal_marker.tsv"
	output:
		down="results/{wave}/GO_global/{type}/{contrast_DE}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt",
		up="results/{wave}/GO_global/{type}/{contrast_DE}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"
	params:
		seed = config["seed"],
		assembly = config["assembly"],
		FC = config["FC"],
		adjp = config["adjp"],
		out_dir = lambda w: "results/GO_global/{}".format(w.contrast_DE),
		up_barplot_out = lambda w: "results/GO_global/{}/{}.upFC.{}.adjp.{}.BP_GO_barplot.pdf".format(w.type, w.contrast_DE, w.FC, w.adjp),
		up_dag_out = lambda w: "results/GO_global/{}/{}.upFC.{}.adjp.{}.BP_GO_dag".format(w.type, w.contrast_DE, w.FC, w.adjp),
		down_barplot_out = lambda w: "results/GO_global/{}/{}.downFC.{}.adjp.{}.BP_GO_barplot.pdf".format(w.type, w.contrast_DE, w.FC, w.adjp),
		down_dag_out = lambda w: "results/GO_global/{}/{}.downFC.{}.adjp.{}.BP_GO_barplot.pdf".format(w.type, w.contrast_DE, w.FC, w.adjp),
		up_consolidated_out = lambda w: "results/GO_global/{}/{}.upFC.{}.adjp.{}.BP_GO_consolidated.tsv".format(w.type, w.contrast_DE, w.FC, w.adjp),
		down_consolidated_out = lambda w: "results/GO_global/{}/{}.downFC.{}.adjp.{}.BP_GO_consolidated.tsv".format(w.type, w.contrast_DE, w.FC, w.adjp)
	conda:
		"../envs/runGO.yaml"
	script:
		"../scripts/runGO_singlecell.R"


		
rule optimal_res:
	input:
		RNA=expand("results/{{wave}}/intermediate/{resolution}_RNA_cluster_markers.tsv",resolution = resolutions),
		Integrated=expand("results/{{wave}}/intermediate/{resolution}_integrated_cluster_markers.tsv",resolution = resolutions)
	output:
		expand("results/{{wave}}/intermediate/{type}_optimal_resolution.txt",type=["RNA","integrated"]) #double brackets escape the expanse
	params:
		out_dir = lambda wildcards: "results/{}/intermediate".format(wildcards.wave)
	conda:
		"../envs/seurat.yaml"
	script:
		"../scripts/optimal_res.py"
		
rule intermediate_processing:
	input:
		"results/{wave}/preprocessing/SeuratObj_{wave}.rds"
	output:	
		RNA="results/{wave}/intermediate/{resolution}_RNA_cluster_markers.tsv",
		Integrated="results/{wave}/intermediate/{resolution}_integrated_cluster_markers.tsv",
		seuratObj="results/{wave}/intermediate/{resolution}_SeuratObj.rds"
	params:
		directory= lambda wildcards: "results/{wave}/intermediate".format(wave = wildcards.wave),
		resolution = lambda wildcards: "{}".format(wildcards.resolution),
		seed=config["seed"],
		max_nFeature=config["max_nFeature"],
		min_nFeature=config["min_nFeature"],
		max_mito=config["max_mito"],
		n_Dims=config["n_Dims"],
		n_cores=config["n_cores"]
	conda:
		"../envs/seurat.yaml"
	script:
		"../scripts/run_FindAllMarkers.R"
		
rule simple_merge:
	input:
		"input/{wave}_metadata.txt" #input is a metafile.tsv with tab delim, first column is 10x files locations and 2nd column is underscore separated metadata
	output:
		seurat_obj="results/{wave}/preprocessing/simple_seurat_{wave}.rds",
	params:
		directory= lambda wildcards: "results/{wave}/preprocessing".format(wave = wildcards.wave),
		seed=config["seed"],
		n_Dims=config["n_Dims"],
		n_cores=config["n_cores"]
	conda:
		"../envs/seurat.yaml"
	script:
		"../scripts/scRNApreprocess_simple_merge.R"
		
rule preprocessing:
	input:
		"input/{wave}_metadata.txt" #input is a metafile.tsv with tab delim, first column is 10x files locations and 2nd column is underscore separated metadata
	output:
		seurat_obj="results/{wave}/preprocessing/SeuratObj_{wave}.rds",
		html="results/{wave}/preprocessing/preprocessing.html"
	params:
		directory= lambda wildcards: "results/{wave}/preprocessing".format(wave = wildcards.wave),
		seed=config["seed"],
		reference=config["reference"]
	conda:
		"../envs/seurat.yaml"
	shell:
		"""Rscript rmarkdown::render("../scripts/scRNApreprocessing.Rmd", output_file = {output.html}, params=list(run2process = {input}, destDir = {params.directory}, seed={params.seed}, reference_set={params.reference}))"""