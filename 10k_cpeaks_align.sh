#!/usr/bin/env bash
###################################################################
#Script Name	  : 10k_cpeaks_align.sh
#Description	  : align fragments to cpeaks reference
#Args           : 
#Author       	: Roy Oelen
#               
###################################################################

# locations of files and programs
OUTPUT_DIR='/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/10k_PBMC_Multiome_nextgem_Chromium_X_cpeaks_atac_fragments.tsv.gz'
INPUT_FILE='/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz'
BARCODES_FILE=''

# where the cpeaks path is
CPEAKS_DIR='/groups/umcg-franke-scrna/tmp02/software/cPeaks/'

# number of cores
CORES=4

# activate the environment
~/miniconda3/bin/activate cpeaks_env

# go to the cpeaks directory
cd ${CPEAKS_DIR}

~/miniconda3/envs/cpeaks_env/bin/python ${CPEAKS_DIR}'main.py' \
    --fragment_path ${INPUT_FILE} \
    --barcode_path ${BARCODES_FILE} \
    --output ${OUTPUT_DIR} \
    --num_cores ${CORES}