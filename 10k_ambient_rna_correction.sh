#!/usr/bin/env bash
###################################################################
#Script Name	  : 10k_ambient_rna_correction.sh
#Description	  : perform ambient RNA correction on the gene expression counts
#Args           : 
#Author       	: Roy Oelen
#               
###################################################################

# locations of files and programs
OUTPUT_FILE='10k_PBMC_Multiome_nextgem_Chromium_X_cellbent_feature_bc_matrix.h5'
INPUT_FILE='10k_PBMC_Multiome_nextgem_Chromium_X_raw_feature_bc_matrix.h5'
export TMPDIR='/groups/umcg-franke-scrna/tmp02/external_datasets/10x_multiome_10k_pbmcs/tmp/'
# number of cells
NCELLS=10000

# load CUDA
ml CUDA/11.7.0

# load the environment
~/miniconda3/bin/activate cellbender_env

# perform correction
~/miniconda3/envs/cellbender_env/bin/cellbender \
    remove-background \
        --input=${INPUT_FILE} \
        --output=${OUTPUT_FILE} \
        --expected-cells=${NCELLS} \
        --cuda \
        --checkpoint ./ckpt.tar.gz

# done
echo 'done'