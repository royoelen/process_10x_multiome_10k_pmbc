#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 10k_cpeaks_to_signac.R
# Function: load cpeaks aligned data and do ATAC QC
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Signac)
# this database is needed for annotations
library(EnsDb.Hsapiens.v86)
# for writing the matrix
library(Matrix)

####################
# Functions        #
####################


read_cpeaks <- function(cpeaks_loc, cellranger_loc, annotations, atac_fragments_append='/outs/atac_fragments.tsv.gz') {
  # read the features
  features <- read.table(paste(cpeaks_loc, '/features.tsv.gz', sep = ''), sep = '\t', header = F)
  features[['peak']] <- paste(features[[1]], ':',  features[[2]], '-', features[[3]], sep = '')
  write.table(features, gzfile(paste(cpeaks_loc, '/named_features.tsv.gz', sep = '')), sep = '\t', row.names = F, col.names = F, quote = F)
  # read the actual counts
  peaks <- ReadMtx(mtx = paste(cpeaks_loc, '/atac_matrix.mtx.gz', sep = ''),
                   cells = paste(cpeaks_loc, '/barcodes.tsv.gz', sep = ''),
                   features = paste(cpeaks_loc, '/named_features.tsv.gz', sep = ''),
                   feature.column = 4)
  
  # create some metadata, for now, we'll first just store the lane here
  barcodes_short <- gsub('(-\\d+)', '', colnames(peaks))
  # create metadata
  metadata <- data.frame(barcode=barcodes_short,
                         barcode_1=colnames(peaks))
  # get the fragments
  fragments_loc <- paste(cellranger_loc, '/', atac_fragments_append, sep = '')
  # create the chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = peaks,
    sep = c(":", "-"),
    fragments = fragments_loc,
    min.cells = 10,
    min.features = 200
  )
  # create object
  seurat_object <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata,
    project = '10x_10k_mo'
  )
  # set the annotations to the object now
  Annotation(seurat_object) <- annotations
  # return result
  return(seurat_object)
}


export_atac_counts_object <- function(signac_object, out_folder) {
  # features
  features_gz <- gzfile(paste(out_folder, 'features.tsv.gz', sep = ''))
  write.table(data.frame(x = rownames(signac_object@assays$peaks@counts)), features_gz, row.names = F, col.names = F, quote = F)
  # barcodes
  barcodes_gz <- gzfile(paste(out_folder, 'barcodes.tsv.gz', sep = ''))
  write.table(data.frame(x = colnames(signac_object@assays$peaks@counts)), barcodes_gz, row.names = F, col.names = F, quote = F)
  # metadata
  metadata_gz <- gzfile(paste(out_folder, 'metadata.tsv.gz', sep = ''))
  write.table(cbind(data.frame(bc = rownames(signac_object@meta.data)), signac_object@meta.data), metadata_gz, row.names = F, col.names = T, quote = F, sep = '\t')
  # and finally the count matrix
  counts_gz <- paste(out_folder, 'matrix.mtx', sep = '')
  writeMM(signac_object@assays$peaks@counts, counts_gz)
  # upon success
  return(0)
}


###################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)

# get annotations from ensemble database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# we use UCSC gencode
seqlevelsStyle(annotations) <- "UCSC"
# and this was aligned on build 38
genome(annotations) <- "hg38"


####################
# Main Code        #
####################


# location of the cpeaks objects
cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/cpeaks_aligned_atac_fragments/'
# get cellranger loc
cellranger_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/'
# location of the gene annotations
gtf_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
# fragment loc
frag_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/'
# location of the arc metadata
arc_metadata_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/10k_PBMC_Multiome_nextgem_Chromium_X_per_barcode_metrics.csv'
# location where to save objects
objects_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/signac/'
object_raw_loc <- paste0(objects_loc, '10x_mo_10k_pbmc_mo_signac_raw.rds')
object_processed_loc <- paste0(objects_loc, '10x_mo_10k_pbmc_mo_signac_processed.rds')
# location of deconstructed data for pycistopic
deconstructed_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/pycistopic/deconstructed/'

# read the object
object_lane <- read_cpeaks(cpeaks_loc, cellranger_loc = frag_loc, annotations = annotations, atac_fragments_append = '10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz')

# read the barcode metadata
arc_metadata <- read.table(arc_metadata_loc, header = T, sep = ',')
rownames(arc_metadata) <- arc_metadata[['barcode']]
# add info to object
object_lane <- AddMetaData(object_lane, arc_metadata['atac_peak_region_fragments'], 'peak_region_fragments')
object_lane <- AddMetaData(object_lane, arc_metadata['atac_fragments'], 'atac_fragments')

# compute nucleosome signal score per cell
object_lane <- NucleosomeSignal(object = object_lane)
# compute TSS enrichment score per cell
object_lane <- TSSEnrichment(object = object_lane, fast = T)

# now the blacklist region
object_lane$blacklist_fraction <- FractionCountsInRegion(
  object = object_lane, 
  assay = 'peaks',
  regions = blacklist_hg38
)

# more QC
object_lane$pct_reads_in_peaks <- object_lane$peak_region_fragments / object_lane$atac_fragments * 100
object_lane$blacklist_ratio <- object_lane$blacklist_fraction / object_lane$peak_region_fragments
object_lane$nucleosome_group <- ifelse(object_lane$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# save the raw object
saveRDS(object_lane, object_raw_loc)

# remove outliers
ncol(object_lane)
# [1] 10968
object_lane <- subset(x = object_lane, subset = nCount_peaks > 3000)
ncol(object_lane)
# [1] 10307
object_lane <- subset(x = object_lane, subset = nCount_peaks < 30000)
ncol(object_lane)
# [1] 9286
object_lane <- subset(x = object_lane, subset = pct_reads_in_peaks > 15)
ncol(object_lane)
# [1] 9285
object_lane <- subset(x = object_lane, subset = blacklist_ratio < 0.05)
ncol(object_lane)
# [1] 9285
object_lane <- subset(x = object_lane, subset = nucleosome_signal < 4)
ncol(object_lane)
# [1] 9285
object_lane <- subset(x = object_lane, subset = TSS.enrichment > 3)
ncol(object_lane)
# [1] 9282

# do normalization
object_lane <- RunTFIDF(object_lane)
object_lane <- FindTopFeatures(object_lane, min.cutoff = 'q0')
object_lane <- RunSVD(object_lane)
# clustering and UMAP
object_lane <- RunUMAP(object = object_lane, reduction = 'lsi', dims = 1:30)
object_lane <- FindNeighbors(object = object_lane, reduction = 'lsi', dims = 1:30)
object_lane <- FindClusters(object = object_lane, verbose = FALSE, algorithm = 3)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene_activities <- GeneActivity(object_lane)
object_lane[['activity']] <- CreateAssayObject(counts = gene_activities)
object_lane <- NormalizeData(
  object = object_lane,
  assay = 'activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(object_lane$nCount_activity)
)

# save result
saveRDS(object_lane, object_processed_loc)

# and export
export_atac_counts_object(object_lane, deconstructed_loc)
