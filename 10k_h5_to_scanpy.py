"""
10k_h5_to_scanpy.py

This script is for converting the output of cellbender into Scanpy objects

authors: Roy Oelen
"""

###########
# imports #
###########

# libraries
import scipy.io
import numpy as np
import pandas as pd
from glob import glob
import pathlib
import os
import scanpy as sc
import re

##################
# files and dirs #
##################

# set up paths
matrix_loc='10k_PBMC_Multiome_nextgem_Chromium_X_cellbent_feature_bc_matrix.h5'
scanpy_object_loc='10k_PBMC_Multiome_nextgem_Chromium_X_gex.h5ad'


#############
# load data #
#############

# load the files
scanpy_object = sc.read_10x_h5(matrix_loc)

# make variable names unique
scanpy_object.var_names_make_unique()

# backup the raw expression in the way that scenic+ wants
scanpy_object.raw = scanpy_object


###################
# quality control #
###################

# get QC metrics
sc.pp.calculate_qc_metrics(
    scanpy_object, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# removing cells with few reads
scanpy_object = scanpy_object[scanpy_object.obs.n_genes_by_counts < 2500, :].copy()
# and high MT content
scanpy_object = scanpy_object[scanpy_object.obs.pct_counts_mt < 5, :].copy()


#################
# normalization #
#################

# normalization
sc.pp.normalize_total(scanpy_object, target_sum=1e4)
sc.pp.log1p(scanpy_object)
sc.pp.highly_variable_genes(scanpy_object, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(scanpy_object, max_value=10)


#####################################
# dimensional recution and clusters #
#####################################

# dimensional reduction
sc.tl.pca(scanpy_object, svd_solver="arpack")
sc.pp.neighbors(scanpy_object, n_neighbors=10, n_pcs=40)
sc.tl.umap(scanpy_object)
# clustering
sc.tl.leiden(
    scanpy_object,
    resolution=0.9,
    random_state=0,
    n_iterations=2,
    directed=False,
)
# marker genes
sc.tl.rank_genes_groups(scanpy_object, "leiden", method="t-test")


#################
# export result #
#################

# write to file
scanpy_object.write(scanpy_object_loc)
