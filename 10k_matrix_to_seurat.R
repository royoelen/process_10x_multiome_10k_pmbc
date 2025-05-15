#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 10k_matrix_to_seurat.R
# Function: read the matrix/features/barcodes into a Seurat object
############################################################################################################################


####################
# libraries        #
####################

# required to create object
library(Seurat)


####################
# Functions        #
####################

#' get a Seurat object from a combination of matrix, barcodes and features files
#' 
#' @param counts_dir directory containing matrix/barcodes/features
#' @param min.cells minimal cells for a gene to be included
#' @param min.features minimal unique transcript for a cell to be included
#' @returns a Seurat object
#' 
add_data <- function(counts_dir='./', min.cells = 0, min.features = 0) {
  # read the actual counts
  counts <- ReadMtx(mtx = paste(counts_dir, 'matrix.mtx.gz', sep = ''),
                    cells = paste(counts_dir, 'barcodes.tsv.gz', sep = ''),
                    features = paste(counts_dir, 'features.tsv.gz', sep = ''),
                    feature.column = 1)
  
  # create some metadata, for now, we'll first just store the lane here
  barcodes_short <- gsub('(-\\d+)', '', colnames(counts))
  
  metadata <- data.frame(barcode=barcodes_short,
                         barcode_1=colnames(counts))
  
  # create the actual object
  seurat_new <- Seurat::CreateSeuratObject(counts = counts,
                                           min.cells = min.cells,
                                           min.features = min.features,
                                           project = "10x_10k_mo",
                                           meta.data = metadata)
  return(seurat_new)
}


####################
# Main Code        #
####################

# where the inputs are
input_directory <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10x_formatted/'
# where to save the object
seurat_object_raw_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10k_seurat_raw.rds'

# get the object
seurat_object <- add_data(input_directory)
# save the result
saveRDS(seurat_object, seurat_object_raw_loc)
