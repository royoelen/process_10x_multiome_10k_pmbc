#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 10k_azimuth_reference_map.R
# Function: reference mapping of multiome data using Azimuth reference PBMC dataset
############################################################################################################################

####################
# libraries        #
####################

library(SeuratData)
library(SeuratDisk)
library(Seurat)


####################
# Functions        #
####################

#' do Azimuth reference mapping of query onto a reference
#' 
#' @param reference the reference Seurat object
#' @param query the query Seurat object
#' @returns a Seurat object
#' 
do_reference_mapping <- function(reference, query){
  # find transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    recompute.residuals = FALSE
  )

  # do the reference mapping
  query <- MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = list(
      celltype.azi.l1 = "celltype.l1",
      celltype.azi.l2 = "celltype.l2",
      celltype.azi.l3 = "celltype.l3",
      celltype.azi.ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  return(query)
}


#' do default Seurat preprocessing on object
#' 
#' @param cbmc Seurat object
#' @returns a Seurat object
#' 
process_lane <- function(cbmc) {
  # filter by UMI
  cbmc <- cbmc[, cbmc@meta.data[['nCount_RNA']] >= 200]
  cbmc <- cbmc[, cbmc@meta.data[['nFeature_RNA']] >= 200]
  # first do RNA
  DefaultAssay(cbmc) <- "RNA"
  # perform visualization and clustering steps
  cbmc <- NormalizeData(cbmc)
  cbmc <- FindVariableFeatures(cbmc)
  cbmc <- ScaleData(cbmc)
  cbmc <- RunPCA(cbmc, verbose = FALSE)
  cbmc <- FindNeighbors(cbmc, dims = 1:30, reduction = 'pca')
  cbmc <- FindClusters(cbmc, resolution = 1.2, verbose = FALSE)
  cbmc <- RunUMAP(cbmc, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = T)
    
  # now do SCT
  cbmc <- SCTransform(cbmc)
  cbmc <- RunPCA(cbmc, verbose = FALSE)
  cbmc <- FindNeighbors(cbmc, dims = 1:30, reduction = 'pca')
  cbmc <- FindClusters(cbmc, resolution = 1.2, verbose = FALSE)
  cbmc <- RunUMAP(cbmc, dims = 1:30, reduction.name = "umap.sct", reduction.key = "sctUMAP_", return.model = T)

  return(cbmc)
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

# location of the reference
reference_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/azimuth_pbmc_reference/multi.h5seurat'
# read the reference
reference <- LoadH5Seurat(reference_loc)
# update to latest version, as our query is also in this format
reference <- UpdateSeuratObject(reference)

# location of the query
seurat_object_raw_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10k_seurat_raw.rds'
# and where to store the result
seurat_object_annotated_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10k_seurat_annotated.rds'

# get the object
seurat_object <- readRDS(seurat_object_raw_loc)
# process
seurat_object <- process_lane(seurat_object)

# perform reference mapping
seurat_object <- do_reference_mapping(reference, seurat_object)

# save object
saveRDS(seurat_object, seurat_object_annotated_loc)

# save the metadata as well, this will contain the celltype annotations which we can load into Pycistopic or Scanpy
celltype_annotation_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10k_seurat_annotated_metadata.tsv.gz'
write.table(seurat_object@meta.data, gzfile(celltype_annotation_loc), row.names = F, col.names = T, sep = '\t', quote = F)
