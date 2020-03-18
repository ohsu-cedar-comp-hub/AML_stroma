import os
import pandas as pd
import numpy as np



####################
# GLOBAL VARIABLES #
####################

RNA_markers = snakemake.input.RNA
integrated_markers = snakemake.input.Integrated
directory = snakemake.params.out_dir 

############
# FUNCTION #
############

def optimize_resolution(marker_files, marker_type):
	
	marker_frames = [pd.read_csv(x,sep='\t',index_col=0) for x in marker_files]
	resolution_dict = {}
	for resolution_file, frame in zip(marker_files,marker_frames):
		
		resolution = resolution_file.split('_')[0].split("/")[-1]
		
		clusters = frame.groupby('cluster')
		
		resolution_frame = pd.DataFrame(columns =frame.cluster.unique(), index=frame.cluster.unique())
		for row_cluster,row_cluster_frame in clusters:
			for col_cluster,col_cluster_frame in clusters:
				percent = (len(np.intersect1d(row_cluster_frame.gene, col_cluster_frame.gene))/len(row_cluster_frame.gene))*100.0
				resolution_frame.loc[row_cluster,col_cluster] = percent
		resolution_dict[resolution] = resolution_frame
	nomination = None
	optimal_resolution = None
	for resolution,resolution_frame in resolution_dict.items():
		eval_metric = (resolution_frame.median().median())
		eval_average_metric = (resolution_frame.sum() - 100).sum()/resolution_frame.sum().shape[0]

		if nomination == None:
			nomination = eval_metric
			optimal_resolution = resolution
		if eval_metric < nomination:
			nomination = eval_metric
			
			optimal_resolution = resolution
			
	
	resolution_top_dict = {}
	for resolution_file, frame in zip(marker_files,marker_frames):
		resolution = resolution_file.split('_')[0].split("/")[-1]
		clusters = frame.groupby('cluster')
		
		resolution_frame = pd.DataFrame(columns =frame.cluster.unique(), index=frame.cluster.unique())
		for row_cluster,row_cluster_frame in clusters:
			for col_cluster,col_cluster_frame in clusters:
				row_cluster_frame = row_cluster_frame.sort_values('avg_logFC', ascending = False).head(n = 100)
				col_cluster_frame = col_cluster_frame.sort_values('avg_logFC', ascending = False).head(n = 100)
				percent = (len(np.intersect1d(row_cluster_frame.gene, col_cluster_frame.gene))/len(row_cluster_frame.gene))*100.0
				resolution_frame.loc[row_cluster,col_cluster] = percent
		resolution_top_dict[resolution] = resolution_frame
	nomination = None
	optimal_resolution = None
	for resolution,resolution_frame in resolution_top_dict.items():
		eval_average_metric = (resolution_frame.sum() - 100).sum()/resolution_frame.sum().shape[0]
		eval_metric = (resolution_frame.median().median())
		#
		if nomination == None:
			nomination = eval_metric
			optimal_resolution = resolution
		if eval_metric < nomination:
			nomination = eval_metric
			
			optimal_resolution = resolution
			
	#print('Selected {optimal_resolution} as optimal resolution'.format(optimal_resolution=optimal_resolution))
	
	with open("{}/{}_optimal_resolution.txt".format(directory,marker_type), "x") as fh:
		fh.write("{}\n".format(optimal_resolution))
		
		
###############
#	 MAIN	#
###############

def main():
	optimize_resolution(RNA_markers,"RNA")
	optimize_resolution(integrated_markers,"integrated")
main()