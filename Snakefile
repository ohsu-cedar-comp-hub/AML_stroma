__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import numpy as np
import datetime
import sys
import os
import pandas as pd
import json

resolutions = np.arange(0.3,1.3,.1)
resolutions = np.round(resolutions,2)

WAVES,=glob_wildcards('input/{wave}_metadata.txt')

def get_contrast_global(wildcards):
	"""Return each contrast provided in the configuration file"""
	return config["diffexp"]["global_contrasts"][wildcards.contrast_DE]

def get_contrast_local(wildcards):
	"""Return each contrast provided in the configuration file"""
	return config["diffexp"]["local_contrasts"][wildcards.contrast_DE]
	
timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"

with open('cluster.json') as json_file:
	json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
	if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
		log_out = os.path.join(os.getcwd(), 'logs', rule)
		os.makedirs(log_out)
		print(log_out)


def message(mes):
	sys.stderr.write("|--- " + mes + "\n")

for wave in WAVES:
	message("10x files in " + wave + " will be processed")

rule all:
	input:
		expand(["results/{wave}/intermediate/{resolution}_RNA_cluster_markers.tsv","results/{wave}/intermediate/{resolution}_integrated_cluster_markers.tsv"],wave=WAVES,resolution=resolutions),
		expand(["results/{wave}/GO_global/{type}/{contrast_DE}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/{wave}/GO_global/{type}/{contrast_DE}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], wave=WAVES,FC=config["FC"], adjp=config["adjp"], contrast_DE=config["diffexp"]["global_contrasts"], type = ["RNA","integrated"]),
		expand("results/{wave}/post/{adjp}_{FC}_{contrast_DE}/optimal_{type}_summary.html",wave = WAVES,FC=config["FC"], adjp=config["adjp"], contrast_DE=config["diffexp"]["global_contrasts"], type = ["RNA","integrated"])

include: "rules/seurat.smk"
