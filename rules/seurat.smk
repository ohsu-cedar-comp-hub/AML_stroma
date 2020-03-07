

rule postprocess:
	input:
		optimal_res_file="results/{contrast}/intermediate/{type}_optimal_resolution.txt", 
		down="results/{contrast}/GO_global/{type}/{contrast_DE}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt",
        up="results/{contrast}/GO_global/{type}/{contrast_DE}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt",
		bulk_markers="results/{contrast}/global_diffexp/{type}_{contrast_DE}_optimal_marker.tsv"
	output:
		html = "results/{contrast}/post/optimal_{type}_summary.html",
		directory = "results/{contrast}/post_{type}",
		seuratObj = "results/{constrast}/post_{type}/optimal_seuratObj.rds"
	params:
		seed=config["seed"]
		input_directory="results/{contrast}/post_processing"
		type = lambda wildcards: "{}".format(wildcards.type)
	conda:
		"../envs/monocle3.yml"
	shell:
		"Rscript rmarkdown::render("../scripts/post_processing.Rmd", output_file = {output.html}, params=list(optimal_res = {input.optimal_res_file}, destDir = {output.directory}, seed={params.seed}, input_dir={params.input_directory}, GOdown={input.down}, GOup = {input.up}, type={params.type}, markers={input.bulk_markers}))"	
		


rule global_diffexp:
	input:
		optimal_res_file="results/{contrast}/intermediate/{type}_optimal_resolution.txt" # extract optimal_res output and get rds from intermediate files
	output:
		marker_file = "results/{contrast}/global_diffexp/{type}_{contrast_DE}_optimal_marker.tsv"
	params:
		seed=config["seed"],
		input_directory = "results/{contrast}/intermediate",
		out_dir = lambda wildcards: "results/{}/global_diffexp".format(wildcards.contrast),
		contrast_DE=get_contrast_global, #named list
		type = lambda wildcards: "{}".format(wildcards.type)
	conda:
		"../envs/monocle3.yml"
	script:
		"../scripts/Global_DE.R"	

		
rule GO_global:
    input:
        degFile="results/{contrast}/global_diffexp/{type}_{contrast_DE}_optimal_marker.tsv"
    output:
        down="results/{contrast}/GO_global/{type}/{contrast_DE}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt",
        up="results/{contrast}/GO_global/{type}/{contrast_DE}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"
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
        up_consolidated_out = lambda w: "results/{GO_global/}/{}.upFC.{}.adjp.{}.BP_GO_consolidated.tsv".format(w.type, w.contrast_DE, w.FC, w.adjp),
        down_consolidated_out = lambda w: "results/GO_global/{}/{}.downFC.{}.adjp.{}.BP_GO_consolidated.tsv".format(w.type, w.contrast_DE, w.FC, w.adjp)
    conda:
        "../envs/runGO.yaml"
    script:
        "../scripts/runGO_singlecell.R"


		
rule optimal_res:
	input:
		RNA=expand("results/{{contrast}}/intermediate/{resolution}_RNA_cluster_markers.tsv",resolution = resolutions),
		Integrated=expand("results/{{contrast}}/intermediate/{resolution}_integrated_cluster_markers.tsv",resolution = resolutions)
	output:
		expand("results/{{contrast}}/{type}_optimal_resolution.txt",type=["RNA","integrated"]) #double brackets escape the expanse
	params:
		out_dir = lambda wildcards: "results/{}".format(wildcards.contrast)
	conda:
		"../envs/monocle3.yml"
	script:
		"../scripts/optimal_res.py"
		
rule intermediate_processing:
	input:
		"results/{contrast}/preprocessing/SeuratObj_{contrast}.rds"
	output:
		directory="results/{contrast}/intermediate",	
		RNA="results/{contrast}/intermediate/{resolution}_RNA_cluster_markers.tsv",
		Integrated="results/{contrast}/intermediate/{resolution}_integrated_cluster_markers.tsv"
		seuratObj="results/{contrast}/intermediate/{resolution}_SeuratObj.rds"
	params:
		resolution = lambda wildcards: "{}".format(wildcards.resolution),
		seed=config["seed"],
		max_nFeature=config["max_nFeature"],
		min_nFeature=config["min_nFeature"],
		max_mito=config["max_mito"],
		n_Dims=config["n_Dims"],
		n_cores=config["n_cores"]
	conda:
		"../envs/monocle3.yml"
	script:
		"../scripts/run_FindAllMarkers.R"
		
		
rule preprocessing:
	input:
		"input/{contrast}.tsv" #input is a metafile.tsv with tab delim, first column is 10x files locations and 2nd column is underscore separated metadata
	output:
		directory="results/{contrast}/preprocessing",
		seurat_obj="results/{contrast}/preprocessing/SeuratObj_{contrast}.rds",
		html="results/{contrast}/preprocessing/preprocessing.html"
	params:
		seed=config["seed"]
	conda:
		"../envs/monocle3.yml"
	shell:
		"Rscript rmarkdown::render("../scripts/scRNApreprocessing.Rmd", output_file = {output.html}, params=list(run2process = {input}, destDir = {output.directory}, seed={params.seed}))"