import os
import pegasus as pg



#RNA_markers = snakemake.input.RNA
#integrated_markers = snakemake.input.integrated

####################
# GLOBAL VARIABLES #
####################
args = get_args()

adata = pg.read_input("MantonBM_nonmix_subset.h5sc")
#directory = snakemake.params.out_dir 

############
# FUNCTION #
############



###############
#     MAIN    #
###############

def main():
    
main()