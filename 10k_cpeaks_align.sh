#!/usr/bin/env bash
###################################################################
#Script Name	  : 10k_cpeaks_align.sh
#Description	  : align fragments to cpeaks reference
#Args           : 
#Author       	: Roy Oelen
#               
###################################################################

# locations of files and programs
OUTPUT_DIR='/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/cpeaks_aligned_atac_fragments/'
INPUT_FILE='/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz'
BARCODES_FILE='/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/filtered_feature_bc_matrix/barcodes.tsv'

# where the cpeaks path is
CPEAKS_DIR='/groups/umcg-franke-scrna/tmp04/software/cPeaks/'

# number of cores
CORES=4

# activate the environment
~/miniconda3/bin/activate cpeaks_env

# go to the cpeaks directory
cd ${CPEAKS_DIR}

# do the alignment to cpeaks
~/miniconda3/envs/cpeaks_env/bin/python ${CPEAKS_DIR}'main.py' \
    --fragment_path ${INPUT_FILE} \
    --barcode_path ${BARCODES_FILE} \
    --output ${OUTPUT_DIR} \
    --num_cores ${CORES}

# zipping in the cpeaks pipeline does not use bgzip, as such we have to rezip before we can index
gunzip 10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz
bgzip 10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz
tabix -p bed 10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz

# let them know we finished
echo 'done'