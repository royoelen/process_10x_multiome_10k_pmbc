#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 10k_multimodal_seurat.R
# Function: merge the ATAC and RNA data
############################################################################################################################


####################
# libraries        #
####################

# required to create object
library(Seurat)
library(Signac)
# and to output the matrices
library(Matrix)
# to zip
library(R.utils)


####################
# Functions        #
####################

#' Make a POSIX-safe String
#'
#' This function takes an input string and replaces any character that is not a letter, digit, hyphen, or underscore with an underscore.
#'
#' @param input_string A character string that needs to be converted to a POSIX-safe format.
#' @return A character string where all non-POSIX-safe characters are replaced with underscores.
#' @examples
#' make_posix_safe("example string!") # Returns "example_string_"
#' make_posix_safe("another@string#") # Returns "another_string_"
#' @export
make_posix_safe <- function(input_string) {
  # replace any character that is not a letter, digit, hyphen, or underscore with an underscore
  posix_safe_string <- gsub("[^a-zA-Z0-9_-]", "_", input_string)
  return(posix_safe_string)
}


#' Write Data Sample
#'
#' This function writes expression data, chromatin data, and metadata to specified folders, ensuring that the sample name is POSIX-safe. It also creates MD5 checksums for the written files.
#'
#' @param expression_data A matrix containing the expression data.
#' @param chromatin_data A matrix containing the chromatin data.
#' @param metadata A data frame containing the metadata.
#' @param sample_name A character string representing the sample name.
#' @param output_folder A character string specifying the output folder path.
#' @param expression_name A character string for the expression data folder name. Default is 'RNA'.
#' @param chromatin_name A character string for the chromatin data folder name. Default is 'peaks'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' write_data_sample(expression_data, chromatin_data, metadata, "sample1", "/path/to/output")
#' }
#' @export
write_data_sample <- function(expression_data, chromatin_data, metadata, sample_name, output_folder, expression_name='RNA', chromatin_name='peaks', binarize_chromatin_matrix=T) {
  # make the sample name posix safe
  sample_name_posix <- make_posix_safe(sample_name)
  # warn if this make the name different
  if (sample_name != sample_name_posix) {
    warning(paste('sample name was not POSIX safe,', sample_name, 'was renamed to', sample_name_posix))
  }
  # create folder for the expression data
  expression_data_folder <- paste0(output_folder, '/', sample_name_posix, '/', expression_name, '/')
  # the chromatin folder as well
  chromatin_data_folder <- paste0(output_folder, '/', sample_name_posix, '/', chromatin_name, '/')
  # create the directories
  dir.create(expression_data_folder, recursive = T)
  dir.create(chromatin_data_folder, recursive = T)
  
  # add barcode as explicit column
  metadata <- cbind(data.frame('barcode' = rownames(metadata)), metadata)
  # write the metadata
  metadata_loc <- paste0(output_folder, '/', sample_name_posix, '/metadata.tsv.gz')
  # write it
  write.table(metadata, gzfile(metadata_loc), row.names = F, col.names = T, quote = F, sep = '\t')
  # and make an md5
  mdfiver::create_md5_for_file(metadata_loc)
  
  # write the expression data
  expression_barcodes_loc <- paste0(expression_data_folder, 'barcodes.tsv.gz')
  write.table(data.frame(x = colnames(expression_data)), gzfile(expression_barcodes_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(expression_barcodes_loc)
  expression_features_loc <- paste0(expression_data_folder, 'features.tsv.gz')
  write.table(data.frame(x = rownames(expression_data)), gzfile(expression_features_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(expression_features_loc)
  # write matrix
  expression_matrix_loc <- paste0(expression_data_folder, 'matrix.mtx')
  writeMM(expression_data, expression_matrix_loc)
  # Compress the MTX file using gzip
  expression_matrix_gz_loc <- paste0(expression_data_folder, 'matrix.mtx.gz')
  gzip(expression_matrix_loc, expression_matrix_gz_loc, overwrite = TRUE)
  # and make md5 checksum
  mdfiver::create_md5_for_file(expression_matrix_gz_loc)
  
  # binarize the chromatin data if requested
  if (binarize_chromatin_matrix) {
    chromatin_data@x[chromatin_data@x > 1] <- 1
  }
  
  # and write the chromatin data
  chromatin_barcodes_loc <- paste0(chromatin_data_folder, 'barcodes.tsv.gz')
  write.table(data.frame(x = colnames(chromatin_data)), gzfile(chromatin_barcodes_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(chromatin_barcodes_loc)
  chromatin_features_loc <- paste0(chromatin_data_folder, 'features.tsv.gz')
  write.table(data.frame(x = rownames(chromatin_data)), gzfile(chromatin_features_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(chromatin_features_loc)
  # write matrix
  chromatin_matrix_loc <- paste0(chromatin_data_folder, 'matrix.mtx')
  writeMM(chromatin_data, chromatin_matrix_loc)
  # Compress the MTX file using gzip
  chromatin_matrix_gz_loc <- paste0(chromatin_data_folder, 'matrix.mtx.gz')
  gzip(chromatin_matrix_loc, chromatin_matrix_gz_loc, overwrite = TRUE)
  # and make md5 checksum
  mdfiver::create_md5_for_file(chromatin_matrix_gz_loc)
  return(0)
}


#' Deconstruct Seurat Object
#'
#' This function deconstructs a Seurat object by extracting RNA and chromatin data, subsetting to specified regions and genes, and writing the data for each sample to the specified output directory.
#'
#' @param seurat_object A Seurat object containing the data.
#' @param regions A character vector of regions to subset from the chromatin data.
#' @param genes A character vector of genes to subset from the RNA data.
#' @param output_dir A character string specifying the output directory path.
#' @param rna_assay A character string for the RNA assay name. Default is 'RNA'.
#' @param chromatin_assay A character string for the chromatin assay name. Default is 'peaks'.
#' @param rna_layer A character string for the RNA data layer. Default is 'counts'.
#' @param chromatin_layer A character string for the chromatin data layer. Default is 'counts'.
#' @param seurat_assignment_column A character string for the column in metadata that contains sample assignments. Default is 'sample_final'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' deconstruct_seurat_object(seurat_object, regions, genes, "/path/to/output")
#' }
#' @export
deconstruct_seurat_object <- function(seurat_object, regions, genes, output_dir, rna_assay='RNA', chromatin_assay='peaks', rna_layer='counts', chromatin_layer='counts', seurat_assignment_column='sample_final', binarize_chromatin_matrix=T, remove_empty_entries=T, verbose=T) {
  # get the metadata
  metadata <- seurat_object@meta.data
  # get the RNA data
  expression_data <- GetAssayData(seurat_object, assay = rna_assay, layer = rna_layer)
  # get the chromatin data
  chromatin_data <- GetAssayData(seurat_object, assay = chromatin_assay, layer = chromatin_layer)
  # remove empty entries
  if (remove_empty_entries) {
    expression_data <- expression_data[rowSums(expression_data) > 0, ]
    chromatin_data <- chromatin_data[rowSums(chromatin_data) > 0, ]
  }
  # check which regions we have
  regions_have <- intersect(rownames(chromatin_data), regions)
  # report on what is missing if present
  if (verbose & length(regions_have) < length(regions)) {
    message(paste('only', as.character(length(regions_have)), 'out of ', length(regions), 'peaks present in data'))
  }
  # check which genes we have
  genes_have <- intersect(rownames(expression_data), genes)
  # report on what is missing if present
  if (verbose & length(genes_have) < length(genes)) {
    message(paste('only', as.character(length(genes_have)), 'out of ', length(genes), 'genes present in data'))
  }
  # subset to those regions and genes
  expression_data <- expression_data[genes_have, ]
  chromatin_data <- chromatin_data[regions_have, ]
  # check each sample
  samples_present <- unique(metadata[[seurat_assignment_column]])
  # and remove entries where we have no data
  samples_present <- samples_present[!is.na(samples_present)]
  # check each sample
  for (sample_name in samples_present) {
    if (verbose) {
      message(paste('starting sample', sample_name))
    }
    # get indices for sample
    indices_sample <- !is.na(metadata[[seurat_assignment_column]]) & metadata[[seurat_assignment_column]] == sample_name
    # subset the data to that sample
    expression_data_sample <- expression_data[, indices_sample]
    chromatin_data_sample <- chromatin_data[, indices_sample]
    metadata_sample <- metadata[indices_sample, ]
    # check if we have enough data for this sample
    if (is.null(nrow(metadata_sample)) | nrow(metadata_sample) < 2) {
      warning(paste('skipping', sample_name, 'due to too few cells'))
    }
    else if(is.null(nrow(expression_data_sample)) | nrow(expression_data_sample) < 2) {
      warning(paste('skipping', sample_name, 'due to too few expressing genes'))
    }
    else if(is.null(nrow(chromatin_data_sample)) | nrow(chromatin_data_sample) < 2) {
      warning(paste('skipping', sample_name, 'due to too few regions with observations'))
    }
    else {
      # write for this sample
      write_data_sample(expression_data = expression_data_sample, 
                        chromatin_data = chromatin_data_sample, 
                        metadata = metadata_sample, 
                        sample_name = sample_name, 
                        output_folder = output_dir, 
                        expression_name = rna_assay, 
                        chromatin_name = chromatin_assay, 
                        binarize_chromatin_matrix = binarize_chromatin_matrix)
    }
    if(verbose) {
      message(paste('finished sample', sample_name))
    }
  }
  return(0)
}


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)


####################
# Main code        #
####################

# where the expression data is (that was processed)
rna_object_annotated_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10k_seurat_annotated.rds'
# where the chromatin data is
atac_processed_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/signac/10x_mo_10k_pbmc_mo_signac_processed.rds'
# and the fragments
fragments_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/rounding/10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz'
# location of the multiome object
multiome_object_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/multimodal/10k_multimodal_object_raw.rds'
multiome_processed_object_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/multimodal/10k_multimodal_object_processed.rds'
# and the one processed and split
multiome_disassembled_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/deconstructed/'
# cre pairs to test
cre_pairs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/confinements/region_to_peak_variant_overlaps.tsv.gz'

# read the objects
rna_object <- readRDS(rna_object_annotated_loc)
atac_object <- readRDS(atac_processed_loc)

# create the chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = GetAssayData(atac_object, layer = 'counts', assay = 'peaks'),
  sep = c(":", "-"),
  fragments = fragments_loc,
  min.cells = 10,
  min.features = 200
)
# and temporary object
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = atac_object@meta.data[colnames(chrom_assay), ]
)
# now get the barcodes present in both
barcodes_both <- intersect(rownames(atac@meta.data), rownames(rna_object@meta.data))
# subset the atac one
atac <- atac[, barcodes_both]
# subset the RNA one as well
rna <- rna_object[, barcodes_both]
# now add the atac metadata to the RNA object
rna@meta.data <- cbind(rna@meta.data, atac@meta.data[match(rownames(rna@meta.data), rownames(atac@meta.data)), c(-1)])
# finally add the chromatin data
rna[['peaks']] <- CreateChromatinAssay(counts = GetAssayData(atac, layer = 'counts', assay = 'peaks'))
# finally remove all extra data
multiome <- DietSeurat(rna)
# set default assay
DefaultAssay(multiome) <- 'RNA'
# remove the assays we don't want
multiome[['SCT']] <- NULL
multiome[['prediction.score.celltype.azi.l1']] <- NULL
multiome[['prediction.score.celltype.azi.l3']] <- NULL
multiome[['celltype.azi.ADT']] <- NULL
# save the raw object
saveRDS(multiome, multiome_object_loc)

# try to do processing of the peaks data
DefaultAssay(multiome) <- "peaks"
# We exclude the first dimension as this is typically correlated with sequencing depth
multiome <- RunTFIDF(multiome)
multiome <- FindVariableFeatures(multiome, assay = 'peaks')
multiome <- RunSVD(multiome)
multiome <- RunUMAP(multiome, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", return.model = T)
multiome <- FindNeighbors(multiome, dims = 2:30, graph.name = 'peaks_snn', reduction = 'lsi')
multiome <- FindClusters(multiome, resolution = 1.2, verbose = FALSE, graph.name = 'peaks_snn')

# next do RNA
DefaultAssay(multiome) <- "RNA"
# perform visualization and clustering steps
multiome <- NormalizeData(multiome)
multiome <- FindVariableFeatures(multiome)
multiome <- ScaleData(multiome)
multiome <- RunPCA(multiome, verbose = FALSE)
multiome <- FindNeighbors(multiome, dims = 1:30)
multiome <- FindClusters(multiome, resolution = 1.2, verbose = FALSE)
multiome <- RunUMAP(multiome, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = T)

# finally do WNN
multiome <- FindMultiModalNeighbors(
  multiome, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30), modality.weight.name = c("RNA.weight", "ATAC.weigth")
)
multiome <- RunUMAP(multiome, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = T)
multiome <- FindClusters(multiome, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

# save this processed object
saveRDS(multiome, multiome_processed_object_loc)


################################
# Export for naive CRE mapping #
################################

# we only have one donor, so let's add that data
multiome@meta.data[['sample']] <- 'S1'

# read the CREs to test
cre_pairs <- read.table(cre_pairs_loc, header = F, sep = '\t')

# start the procedure
deconstruct_seurat_object(
  seurat_object = multiome, 
  regions = unique(cre_pairs[['V1']]), 
  genes = unique(cre_pairs[['V2']]), 
  output_dir = multiome_disassembled_loc, 
  rna_assay = 'RNA', 
  chromatin_assay = 'peaks', 
  rna_layer = 'counts', 
  chromatin_layer = 'counts', 
  seurat_assignment_column = 'sample', 
  binarize_chromatin_matrix = T, 
  verbose = T
)


###############################
# make course-specific object #
###############################

# remove any cell type annotation
multiome@meta.data[['predicted.celltype.azi.l1.score']] <- NULL
multiome@meta.data[['predicted.celltype.azi.l1']] <- NULL
multiome@meta.data[['predicted.celltype.azi.l2.score']] <- NULL
multiome@meta.data[['predicted.celltype.azi.l2']] <- NULL
multiome@meta.data[['predicted.celltype.azi.l3.score']] <- NULL
multiome@meta.data[['predicted.celltype.azi.l3']] <- NULL

# and save the object
multiome_processed_object_noct_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/multimodal/10k_multimodal_object_processed_nocelltype.rds'
saveRDS(multiome, multiome_processed_object_noct_loc)
