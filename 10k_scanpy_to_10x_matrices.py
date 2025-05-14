"""
10k_scanpy_to_10x_matrices.py

This script is for converting the Scanpy object to 10x format to load it into Seurat for celltype annotation

authors: Roy Oelen
"""

###########
# imports #
###########

# import the libraries required
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sparse
import scipy.io as sio
import scipy
import gzip
import os


#############
# functions #
#############

def scanpy_object_to_10x_format(scanpy_object, output_loc):
    """
        
        Parameters
        ----------
        scanpy_object_loc : AnnData
            the AnnData Scanpy object to deconstruct
        output_loc : String
            the location of where to put the deconstructed parts
        
        Returns
        -------
        result
           0 if succesfull
    """
    # make the genes unique
    scanpy_object.var_names_make_unique()
    # location of outputs
    matrix_loc = ''.join([output_loc, 'matrix.mtx.gz'])
    features_loc = ''.join([output_loc, 'features.tsv.gz'])
    barcodes_loc = ''.join([output_loc, 'barcodes.tsv.gz'])
    metadata_loc = ''.join([output_loc, 'metadata.tsv.gz'])
    # write barcodes and features
    with gzip.open(barcodes_loc, 'wb') as f:
        pd.DataFrame(data = {'barcodes' : scanpy_object.obs_names.tolist()}).to_csv(f, sep = '\t', header = False, index = False)
    with gzip.open(features_loc, 'wb') as f:
        pd.DataFrame(data = {'features' : scanpy_object.var_names.tolist()}).to_csv(f, sep = '\t', header = False, index = False)
    # write the metadata
    with gzip.open(metadata_loc, 'wb') as f:
        scanpy_object.obs.to_csv(f, sep = '\t', header = True, index = True, index_label = 'barcode')
    # write the matrix
    with gzip.open(matrix_loc, 'wb') as f:
        scipy.io.mmwrite(f, scanpy_object.X.transpose())
    # let us know we were succesfull
    return 0


##################
# files and dirs #
##################

# set up paths
h5_object_loc = '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/cellbender/10k_PBMC_Multiome_nextgem_Chromium_X_cellbent_feature_bc_matrix_filtered.h5'
output_directory = '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/seurat/10x_formatted/'


#############
# main code #
#############

# load the files
scanpy_object = sc.read_10x_h5(h5_object_loc)
# make variable names unique
scanpy_object.var_names_make_unique()
# convert
scanpy_object_to_10x_format(scanpy_object, output_directory)
