"""

10k_calculate_atac_rna_betas.py

This script is used to combine calculate betas and p values between chromatin and epxression data

authors: Roy Oelen

example usage:

python 10k_calculate_atac_rna_betas.py \
    --expression_folder /groups/umcg-franke-scrna/tmp02/external_datasets/10x_multiome_10k_pbmcs/deconstructed/S1/RNA/ \
    --chromatin_folder /groups/umcg-franke-scrna/tmp02/external_datasets/10x_multiome_10k_pbmcs/deconstructed/S1/peaks/ \
    --output_folder /groups/umcg-franke-scrna/tmp02/external_datasets/10x_multiome_10k_pbmcs/naive_cres/S1/ \
    --use_gpu \
    --check_order_intersect \
    --cre_loc /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/confinements/region_to_peak_variant_overlaps.tsv.gz \
    --n_perm 10

"""

#############
# libraries #
#############

# standard format for CPU
import numpy as np
# standard format for GPU
import cupy as cp
# we need these matrices, coo to get the matrix, and csr to actually use
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
# we read the mtx.gz files
from scipy.io import mmread
# to test if we are dealing with a scipy or cupy implementation
from scipy.sparse import issparse
# the cupy implementation of the matrix
from cupyx.scipy.sparse import csr_matrix as cp_csr_matrix
from cupyx.scipy.sparse import coo_matrix as cp_coo_matrix
# to test if we are dealing with a scipy or cupy implementation
from cupyx.scipy.sparse import issparse as cp_issparse
# we use the scipy versions of these tools in CPU-only scenarios
from sklearn.preprocessing import PowerTransformer
from sklearn.preprocessing import StandardScaler
from scipy.stats import yeojohnson
# for binomial regression using GPU
import torch
import torch.nn as nn
import torch.optim as optim
# to get betas and p
import statsmodels.api as sm
# multithreading on CPU instead
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
# the CRE confinement can be loaded using pandas
import pandas as pd
# files might be zipped
import gzip
# we need to suppress some warnings for not overloading with text
import os
import sys
# finally, also parse the arguments
import argparse
# for permutations
import random
# for warning
import warnings


#############
# functions #
#############

def load_matrix_cpu(matrix_loc):
    """
    Load a sparse matrix from a .mtx.gz file using SciPy.

    Parameters:
    matrix_loc (str): The file path to the .mtx.gz file.

    Returns:
    scipy.sparse.coo_matrix: The loaded sparse matrix in COO format.
    """
    # load the matrix as coo
    with gzip.open(matrix_loc, 'rt') as f:
        scipy_sparse_matrix = mmread(f)
    return scipy_sparse_matrix


def load_matrix_gpu(matrix_loc):
    """
    Load a sparse matrix from a .mtx.gz file using SciPy and convert it to a CuPy sparse matrix.

    Parameters:
    matrix_loc (str): The file path to the .mtx.gz file.

    Returns:
    cupyx.scipy.sparse.csr_matrix: The loaded sparse matrix in CSR format.
    """
    # load using scipy cpu function
    scipy_sparse_matrix = load_matrix_cpu(matrix_loc)
    # convert the data type to be compatible with CuPy
    scipy_sparse_matrix.data = scipy_sparse_matrix.data.astype(cp.float64)
    # convert the SciPy coo sparse matrix to a CuPy sparse matrix
    cupy_sparse_matrix = cp_csr_matrix(scipy_sparse_matrix)
    # clear memory
    del scipy_sparse_matrix
    return cupy_sparse_matrix


def load_matrix(matrix_loc):
    """
    Load a sparse matrix from a .mtx.gz file using either CPU or GPU implementation.

    Parameters:
    matrix_loc (str): The file path to the .mtx.gz file.

    Returns:
    scipy.sparse.csr_matrix or cupyx.scipy.sparse.csr_matrix: The loaded sparse matrix in CSR format.
    """
    # use GPU or CPU implementation
    if use_gpu:
        return load_matrix_gpu(matrix_loc)
    else:
        return load_matrix_cpu(matrix_loc).tocsr()

    
def yeo_johnson_normalization_and_scale_row_cpu(row):
    """
    Apply Yeo-Johnson normalization and scaling to a single row using CPU.

    Parameters:
    row (scipy.sparse.csr_matrix): A single row of the sparse matrix.

    Returns:
    numpy.ndarray: The normalized and scaled row.
    """
    # init transformer object
    transformer = PowerTransformer(method='yeo-johnson', standardize=True)
    # init scaler object
    scaler = StandardScaler()
    # convert sparse row to dense format
    row_dense = row.toarray().flatten()
    # use the transformer
    row_transformed = transformer.fit_transform(row_dense.reshape(-1, 1)).flatten()
    # use the scaler
    row_scaled = scaler.fit_transform(row_transformed.reshape(-1, 1)).flatten()
    return row_scaled


def yeo_johnson_normalization_and_scale_cpu(matrix):
    """
    Apply Yeo-Johnson normalization and scaling to each row of the matrix using CPU.

    Parameters:
    matrix (scipy.sparse.csr_matrix): The sparse matrix.

    Returns:
    numpy.ndarray: The normalized and scaled matrix.
    """
    # setup threadpool
    with ProcessPoolExecutor() as executor:
        # do row per thread
        normalized_rows = list(executor.map(yeo_johnson_normalization_and_scale_row_cpu, matrix))
    # put restult in numpy array
    return np.array(normalized_rows)


def yeo_johnson_transform(x, lmbda):
    """
    Apply the Yeo-Johnson transformation to a CuPy array.

    Parameters:
    x (cupy.ndarray): The input array.
    lmbda (float): The lambda parameter for the Yeo-Johnson transformation.

    Returns:
    cupy.ndarray: The transformed array.
    """
    """Apply the Yeo-Johnson transformation to a CuPy array."""
    
    # create output array that is the same shape as x, but with zeroes
    out = cp.zeros_like(x)
    # create booleans of the positive (or zero) and negative values
    pos = x >= 0
    neg = ~pos

    # if the supplied lambda is 0, the transformation for non-negative values is the natural logarithm of x + 1
    if lmbda == 0:
        out[pos] = cp.log1p(x[pos])
    # otherwise it is a bit more complicated
    else:
        out[pos] = (cp.power(x[pos] + 1, lmbda) - 1) / lmbda
    # if the lambda is 2, the transformation for negative values is the negative natural logarithm of 1 - x
    if lmbda == 2:
        out[neg] = -cp.log1p(-x[neg])
    # otherwise it is again a bit more complicated
    else:
        out[neg] = -((cp.power(-x[neg] + 1, 2 - lmbda) - 1) / (2 - lmbda))
    
    return out


def handle_nan_inf(matrix):
    """
    Replace NaN and infinite values in the matrix with zeros.

    Parameters:
    matrix (cupyx.scipy.sparse.csr_matrix): The sparse matrix.

    Returns:
    cupyx.scipy.sparse.csr_matrix: The matrix with NaN and infinite values replaced by zeros.
    """
    # make all NAN zero
    matrix.data[cp.isnan(matrix.data)] = 0
    # make all infinite zero
    matrix.data[cp.isinf(matrix.data)] = 0
    return matrix


def yeo_johnson_normalization_and_scale_row_gpu(row):
    """
    Apply Yeo-Johnson normalization and scaling to a single row using GPU.

    Parameters:
    row (cupyx.scipy.sparse.csr_matrix): A single row of the sparse matrix.

    Returns:
    cupy.ndarray: The normalized and scaled row.
    """
    # make the row into dense format
    row_dense = row.toarray().flatten()
    
    # find the optimal lambda for the Yeo-Johnson transformation
    _, lmbda = yeojohnson(row_dense.get())
    
    # apply the Yeo-Johnson transformation using the optimal lambda we determined
    row_transformed = yeo_johnson_transform(cp.array(row_dense), lmbda)
    
    # Scale the transformed row using CuPy, same implementation as StandardScaler of scipy
    row_mean = cp.mean(row_transformed)
    row_std = cp.std(row_transformed) + 1e-8  # Add epsilon to avoid division by zero
    row_scaled = (row_transformed - row_mean) / row_std
    
    return row_scaled


def yeo_johnson_normalization_and_scale_gpu(cp_matrix):
    """
    Apply Yeo-Johnson normalization and scaling to each row of the matrix using GPU.

    Parameters:
    cp_matrix (cupyx.scipy.sparse.csr_matrix): The sparse matrix.

    Returns:
    cupy.ndarray: The normalized and scaled matrix.
    """
    # do each row, cupy should leverage paralellism automatically
    normalized_rows = cp.array([yeo_johnson_normalization_and_scale_row_gpu(row) for row in cp_matrix])
    return normalized_rows


def yeo_johnson_normalization_and_scale(matrix):
    """
    Apply Yeo-Johnson normalization and scaling to each row of the matrix using either CPU or GPU.

    Parameters:
    matrix (scipy.sparse.csr_matrix or cupyx.scipy.sparse.csr_matrix): The sparse matrix.

    Returns:
    numpy.ndarray or cupy.ndarray: The normalized and scaled matrix.
    """
    # depending on whether we have a GPU, do the CPU or GPU implementation
    if use_gpu:
        return yeo_johnson_normalization_and_scale_gpu(matrix)
    else:
        return yeo_johnson_normalization_and_scale_cpu(matrix)


def binomial_regression_gpu(binary_array_2d, normal_array_2d, epochs=100):
    """
    Perform binomial regression on 2D arrays using GPU acceleration with PyTorch.

    Parameters:
    binary_array_2d (cupy.ndarray): A 2D CuPy array of binary response variables.
    normal_array_2d (cupy.ndarray): A 2D CuPy array of predictor variables.
    epochs (int): The number of epochs for training the logistic regression model. Default is 100.

    Returns:
    dict: A dictionary containing three keys:
          - 'p': A 1D CuPy array of p-values for each row.
          - 'beta': A 1D CuPy array of beta coefficients for each row.
          - 'std_err': A 1D CuPy array of standard errors for each row.

    Example:
    >>> binary_array_2d = cp.array([[0, 1, 0], [1, 0, 1]])
    >>> normal_array_2d = cp.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    >>> results = binomial_regression_gpu(binary_array_2d, normal_array_2d)
    >>> print(results['beta'], results['p'], results['std_err'])
    [0.5, 0.7] [0.045, 0.032] [0.1, 0.15]

    Notes:
    - This function uses PyTorch to perform logistic regression on each row of the input arrays.
    - The logistic regression model is trained using binary cross-entropy loss and the Adam optimizer.
    - The beta coefficients are extracted from the trained model, and p-values and standard errors are computed using statsmodels.
    - If a singular matrix error occurs during the fitting process, the p-value and standard error are set to NaN.
    """
    # create custom class of regression model that uses GPU
    class LogisticRegressionModel(nn.Module):
        # constructor
        def __init__(self):
            super(LogisticRegressionModel, self).__init__()
            self.linear = nn.Linear(1, 1)
        # move through array
        def forward(self, x):
            return torch.sigmoid(self.linear(x))
    
    # Initialize the arrays of p values, betas, and standard errors
    betas = cp.zeros(normal_array_2d.shape[0])
    p_values = cp.ones(normal_array_2d.shape[0])
    std_errs = cp.zeros(normal_array_2d.shape[0])
    # do each row
    for i in range(normal_array_2d.shape[0]):
        # create GPU logistic regression model
        model = LogisticRegressionModel().cuda()
        # cross entropy loss
        criterion = nn.BCELoss()
        # optimizer
        optimizer = optim.Adam(model.parameters(), lr=0.01)
        # setup the training data, so the original data
        X_train = torch.tensor(normal_array_2d[i].get(), dtype=torch.float32).reshape(-1, 1).cuda()
        y_train = torch.tensor(binary_array_2d[i].get(), dtype=torch.float32).reshape(-1, 1).cuda()
        # do multiple iterations to train the model
        for epoch in range(epochs):
            # do a model
            model.train()
            # set gradient to zero
            optimizer.zero_grad()
            # do the predictions
            outputs = model(X_train)
            # get the loss
            loss = criterion(outputs, y_train)
            # and update the model
            loss.backward()
            optimizer.step()
        
        # get beta
        beta = model.linear.weight.data.cpu().numpy()[0][0]
        # put into list
        betas[i] = beta
        
        # Compute p-values and standard errors using statsmodels
        try:
            # constant term for intercept
            X_sm = sm.add_constant(normal_array_2d[i].get())
            # make model using data and intercept
            logit_model = sm.Logit(binary_array_2d[i].get(), X_sm)
            # get the result of the model
            result = logit_model.fit(disp=0)
            # extract p-value for the predictor
            p_value = result.pvalues[1]
            # extract standard error for the predictor
            std_err = result.bse[1]
        except np.linalg.LinAlgError:
            # set p-value and standard error to NaN if Singular matrix error occurs and we don't have a good fit
            p_value = np.nan
            std_err = np.nan

        # add p value and standard error in list
        p_values[i] = p_value
        std_errs[i] = std_err
    
    # put results in a dictionary and return those
    results = {'p': p_values, 'beta': betas, 'std_err': std_errs}
    return results


def binomial_regression_single_row_cpu(X_row, y_row):
    """
    Perform binomial regression for a single row of data.

    Parameters:
    X_row (array-like): A 1D array of predictor variables for a single observation.
    y_row (array-like): A 1D array of binary response variables for a single observation.

    Returns:
    tuple: A tuple containing the beta coefficient, the p-value, and the standard error for the predictor variable.
           If a singular matrix error occurs, all values are set to NaN.

    Example:
    >>> X_row = [1.0, 2.0, 3.0]
    >>> y_row = [0, 1, 0]
    >>> beta, p_value, std_err = binomial_regression_single_row_cpu(X_row, y_row)
    >>> print(beta, p_value, std_err)
    0.5 0.045 0.1

    Notes:
    - This function uses statsmodels to perform logistic regression.
    - A constant term is added to the predictor variables to account for the intercept.
    - If a singular matrix error occurs during the fitting process, all values are set to NaN.
    """
    # do binomial regression for a single row
    try:
        # constant term for intercept
        X_sm = sm.add_constant(X_row)
        # make model use data and use intercept
        logit_model = sm.Logit(y_row, X_sm)
        # get result of the model
        result = logit_model.fit(disp=0)
        # extract the coefficient/beta
        beta = result.params[1]
        # get the p-value
        p_value = result.pvalues[1]
        # get the standard error
        std_err = result.bse[1]
    except np.linalg.LinAlgError:
        # set both beta, p, and std_err to NAN if we encounter an error
        beta = np.nan
        p_value = np.nan
        std_err = np.nan
    return beta, p_value, std_err


def binomial_regression_cpu(binary_array_2d, normal_array_2d):
    """
    Perform binomial regression on 2D arrays using multithreading on the CPU.

    Parameters:
    binary_array_2d (array-like): A 2D array of binary response variables.
    normal_array_2d (array-like): A 2D array of predictor variables.

    Returns:
    dict: A dictionary containing three keys:
          - 'p': A 1D array of p-values for each row.
          - 'beta': A 1D array of beta coefficients for each row.
          - 'std_err': A 1D array of standard errors for each row.

    Example:
    >>> binary_array_2d = np.array([[0, 1, 0], [1, 0, 1]])
    >>> normal_array_2d = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    >>> results = binomial_regression_cpu(binary_array_2d, normal_array_2d)
    >>> print(results['beta'], results['p'], results['std_err'])
    [0.5, 0.7] [0.045, 0.032] [0.1, 0.15]

    Notes:
    - This function uses ThreadPoolExecutor for multithreading to perform regression on each row concurrently.
    - The helper function `binomial_regression_single_row_cpu` is used to perform regression on a single row.
    - The results are stored in arrays and returned as a dictionary.
    """
    # initialize the arrays of p values, betas, and standard errors, given the shape of the input
    betas = np.zeros(normal_array_2d.shape[0])
    p_values = np.ones(normal_array_2d.shape[0])
    std_errs = np.zeros(normal_array_2d.shape[0])

    # use ThreadPoolExecutor for multithreading
    with ThreadPoolExecutor() as executor:
        # store results in list
        futures = []
        # check each row in the first array (both arrays are the same size)
        for i in range(normal_array_2d.shape[0]):
            # do the row
            futures.append(executor.submit(binomial_regression_single_row_cpu, normal_array_2d[i].flatten(), binary_array_2d[i].flatten()))
        # go through all results
        for i, future in enumerate(futures):
            # grab the beta, p-value, and standard error from the result of the thread
            beta, p_value, std_err = future.result()
            # put into the list
            betas[i] = beta
            p_values[i] = p_value
            std_errs[i] = std_err
    
    # put all results in a dictionary and return those
    results = {'p': p_values, 'beta': betas, 'std_err': std_errs}
    return results


def binomial_regression(binary_array_2d, normal_array_2d):
    """
    Perform binomial regression on 2D arrays, using GPU acceleration if available.

    Parameters:
    binary_array_2d (array-like): A 2D array of binary response variables.
    normal_array_2d (array-like): A 2D array of predictor variables.

    Returns:
    dict: A dictionary containing two keys:
          - 'p': A 1D array of p-values for each row.
          - 'beta': A 1D array of beta coefficients for each row.

    Example:
    >>> binary_array_2d = np.array([[0, 1, 0], [1, 0, 1]])
    >>> normal_array_2d = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    >>> results = binomial_regression(binary_array_2d, normal_array_2d)
    >>> print(results['beta'], results['p'])
    [0.5, 0.7] [0.045, 0.032]

    Notes:
    - This function selects between GPU and CPU implementations of binomial regression based on the availability of GPU.
    - The GPU implementation uses PyTorch for logistic regression, while the CPU implementation uses statsmodels with multithreading.
    """
    # depending on GPU, use one or the other method
    if use_gpu:
        return binomial_regression_gpu(binary_array_2d, normal_array_2d)
    else:
        return binomial_regression_cpu(np.array(binary_array_2d), normal_array_2d)


def binomial_regression_chunked_sparse(binary_sparse_matrix, normal_array_2d, chunk_size=1000):
    """
    Perform binomial regression on sparse 2D arrays in chunks.

    Parameters:
    binary_sparse_matrix (csr_matrix): A sparse matrix of binary response variables.
    normal_array_2d (array-like): A 2D array of predictor variables.
    chunk_size (int): The number of rows to process in each chunk.

    Returns:
    dict: A dictionary containing three keys:
          - 'p': A 1D array of p-values for each row.
          - 'beta': A 1D array of beta coefficients for each row.
          - 'std_err': A 1D array of standard errors for each row.
    """
    # grab the number of rows
    num_rows = binary_sparse_matrix.shape[0]
    betas = np.zeros(num_rows)
    p_values = np.ones(num_rows)
    std_errs = np.zeros(num_rows)
    
    for start in range(0, num_rows, chunk_size):
        end = min(start + chunk_size, num_rows)
        chunk_binary_sparse = binary_sparse_matrix[start:end]
        chunk_normal = normal_array_2d[start:end]
        
        # Convert the sparse chunk to dense format
        chunk_binary_dense = chunk_binary_sparse.todense()
        
        # Perform binomial regression on the chunk
        chunk_results = binomial_regression(chunk_binary_dense, chunk_normal)
        
        # Store the results
        if use_gpu:
            # we need to use .get for cupy arrays
            betas[start:end] = chunk_results['beta'].get()
            p_values[start:end] = chunk_results['p'].get()
            std_errs[start:end] = chunk_results['std_err'].get()
        else:
            betas[start:end] = chunk_results['beta']
            p_values[start:end] = chunk_results['p']
            std_errs[start:end] = chunk_results['std_err']
    
    results = {'p': p_values, 'beta': betas, 'std_err': std_errs}
    return results


def shuffle_csr_columns(matrix, seed=None):
    """
    Shuffle the columns of a CSR (Compressed Sparse Row) matrix.

    Parameters:
    matrix (csr_matrix): The input CSR matrix whose columns need to be shuffled.
    seed (int, optional): The seed for the random number generator. Default is None.

    Returns:
    csr_matrix: A new CSR matrix with shuffled columns.

    Example:
    >>> from scipy.sparse import csr_matrix
    >>> import numpy as np
    >>> matrix = csr_matrix(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
    >>> shuffled_matrix = shuffle_csr_columns(matrix, seed=42)
    >>> print(shuffled_matrix.toarray())
    [[3 1 2]
     [6 4 5]
     [9 7 8]]
    """
    # get columns
    num_cols = matrix.shape[1]
    # get shuffled indices
    if use_gpu:
        # with a seed
        if seed is not None:
            cp.random.seed(seed)
        shuffled_indices = cp.random.permutation(num_cols)
    else:
        if seed is not None:
            np.random.seed(seed)
        shuffled_indices = np.random.permutation(num_cols)
    # shuffle based on the indices
    shuffled_matrix = matrix[:, shuffled_indices]
    return shuffled_matrix

###################
# parse arguments #
###################

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--expression_folder', type = str, help = 'location of the expression files (string)')
parser.add_argument('-a', '--chromatin_folder', type = str, help = 'location of the chromatin files (string)')
parser.add_argument('-o', '--output_folder', type = str, help = 'location of the  (string)')
parser.add_argument('-g', '--use_gpu', action = 'store_true', help = 'use GPU acceleration')
parser.add_argument('-i', '--check_order_intersect', action = 'store_true', help = 'check the barcodes to intersect and order barcodes')
parser.add_argument('-c', '--cre_loc', type = str, help = 'location of the CRE confinement file (string)')
parser.add_argument('-n', '--nrow_chunk', type = int, help = 'number of rows per chunk, zero means no chunking (integer)', default = 0)
parser.add_argument('-s', '--sample_name', type = str, help = 'sample name to add to all output tables, leave parameter out for no header (string)', default = None)
parser.add_argument('-p', '--n_perm', type = int, help = 'number of permutations to perform (integer)', default = 0)
parser.add_argument('-d', '--seeds_file_loc', type = str, help = 'list of permutation seeds to use (string)', default = None)
args = parser.parse_args()

# whether we use GPU or not
use_gpu = args.use_gpu
# whether we need to check for ordering and intersects of the two modalities
check_order_and_intersect = args.check_order_intersect
# location of the rna data
rna_data_loc = args.expression_folder
# location of the atac data
atac_data_loc = args.chromatin_folder
# location of the CREs to test
cre_loc = args.cre_loc
# location where to place the outputs
output_folder = args.output_folder
# how many rows per chunk
row_chunk_size = args.nrow_chunk
# the name of the sample to put as header
sample_name = args.sample_name
# get the number of permutations
n_perm = args.n_perm
# get the seeds file
seeds_file_loc = args.seeds_file_loc

# location of the RNA matrix
rna_matrix_loc = ''.join([rna_data_loc, '/matrix.mtx.gz'])
# location of the ATAC matrix
atac_matrix_loc = ''.join([atac_data_loc, '/matrix.mtx.gz'])

# location of features and barcodes
rna_barcodes_loc = ''.join([rna_data_loc, '/barcodes.tsv.gz'])
rna_features_loc = ''.join([rna_data_loc, '/features.tsv.gz'])
atac_barcodes_loc = ''.join([atac_data_loc, '/barcodes.tsv.gz'])
atac_features_loc = ''.join([atac_data_loc, '/features.tsv.gz'])

# specifically the betas and p values
output_beta_loc = ''.join([output_folder, '/beta.txt.gz'])
output_p_loc = ''.join([output_folder, '/p.txt.gz'])
output_se_loc = ''.join([output_folder, '/se.txt.gz'])


#################
# load all data #
#################

# load the two matrices
rna_matrix = load_matrix(rna_matrix_loc)
atac_matrix = load_matrix(atac_matrix_loc)

# load the barcodes, we know how many to expect
rna_barcodes = values = [None] * rna_matrix.shape[1]
# open the file and put each barcode in the list
with gzip.open(rna_barcodes_loc, 'r') as file:
    for i, line in enumerate(file):
        rna_barcodes[i] = line.strip()
# repeat for the ATAC barcodes
atac_barcodes = values = [None] * atac_matrix.shape[1]
# open the file and put each barcode in the list
with gzip.open(atac_barcodes_loc, 'r') as file:
    for i, line in enumerate(file):
        atac_barcodes[i] = line.strip()

# load the features as well, again we know how many to expect
rna_features = values = [None] * rna_matrix.shape[0]
# open the file and put each feature in the list
with gzip.open(rna_features_loc, 'r') as file:
    for i, line in enumerate(file):
        rna_features[i] = line.strip().decode('utf-8')
# repeat for the ATAC features
atac_features = values = [None] * atac_matrix.shape[0]
# open the file and put each feature in the list
with gzip.open(atac_features_loc, 'r') as file:
    for i, line in enumerate(file):
        atac_features[i] = line.strip().decode('utf-8')


#######################
# intersect and order #
#######################

# now we will make sure we have the same barcodes in the same order, if that was not already the case
if check_order_and_intersect:
    # check which are present in both
    barcodes_both = list(set(rna_barcodes) & set(atac_barcodes))
    # get the indices in the original lists
    indices_barcodes_rna = [i for i, value in enumerate(rna_barcodes) if value in barcodes_both]
    indices_barcodes_atac = [i for i, value in enumerate(atac_barcodes) if value in barcodes_both]
    # subset those barcodes now
    rna_barcodes = [rna_barcodes[i] for i in indices_barcodes_rna]
    atac_barcodes = [atac_barcodes[i] for i in indices_barcodes_atac]
    # and the two matrices
    rna_matrix = rna_matrix[:, indices_barcodes_rna]
    atac_matrix = atac_matrix[:, indices_barcodes_atac]


#####################
# overlap with CREs #
#####################

# load those CREs
cres = pd.read_csv(cre_loc, header = None, delimiter = '\t')
# set the column names for clarity sake
cres.columns = ['region', 'gene']

# extract those regions and genes
regions = cres['region'].values
genes = cres['gene'].values
# map the names of the features to their index
region_to_num = {char: i for i, char in enumerate(atac_features)}
gene_to_num = {char: i for i, char in enumerate(rna_features)}
# get indices where we have the region and the gene
valid_indices = [(i, region_to_num[region], gene_to_num[gene]) 
                 for i, (region, gene) in enumerate(zip(regions, genes)) 
                 if region in region_to_num and gene in gene_to_num]
# Separate the indices and their original positions
if valid_indices:
    original_positions, valid_indices_regions, valid_indices_genes = zip(*valid_indices)
else:
    original_positions, valid_indices_regions, valid_indices_genes = [], [], []

# Subset the region matrix
subset_regions = atac_matrix[list(valid_indices_regions), :] if valid_indices_regions else atac_matrix[:0, :]
# Subset the gene matrix
subset_genes = rna_matrix[list(valid_indices_genes), :] if valid_indices_genes else rna_matrix[:0, :]


#########################
# normalize and regress #
#########################

# do yeo johnson normalization on the rna data
rna_normal = yeo_johnson_normalization_and_scale(subset_genes)

# do the regression
regression_results = None
# if the chunk size is bigger than zero, we do chunking
if row_chunk_size > 0:
    regression_results = binomial_regression_chunked_sparse(subset_regions, rna_normal, chunk_size = row_chunk_size)
# otherwise we do everything all in once
else:
    regression_results = binomial_regression(subset_regions.todense(), rna_normal)

# add the results back in the order of the confinement file
ordered_ps = [None] * len(regions)
ordered_betas = [None] * len(regions)
ordered_ses = [None] * len(regions)

# Place the results back in the original order
for i in range(0, len(original_positions), 1):
    # get the resulting p
    p = regression_results['p'][i]
    # and beta
    beta = regression_results['beta'][i]
    # and se
    std_err = regression_results['std_err'][i]
    # get position in original list
    i_original = original_positions[i]
    # put in the lists at those positions
    ordered_ps[i_original] = p
    ordered_betas[i_original] = beta
    ordered_ses[i_original] = std_err

# write the lists to gzipped text files
with gzip.open(output_p_loc, 'wt') as f:
    if sample_name is not None:
        f.write(f"{sample_name}\n")
    for item in ordered_ps:
        f.write(f"{item}\n")

with gzip.open(output_beta_loc, 'wt') as f:
    if sample_name is not None:
        f.write(f"{sample_name}\n")
    for item in ordered_betas:
        f.write(f"{item}\n")

with gzip.open(output_se_loc, 'wt') as f:
    if sample_name is not None:
        f.write(f"{sample_name}\n")
    for item in ordered_ses:
        f.write(f"{item}\n")

####################
# permutation runs #
####################

# now do permutations
if n_perm > 0:
    # set up the seeds
    seeds = []
    # get from file if a file was supplied
    if seeds_file_loc is not None:
        # open seeds file
        with gzip.open(seeds_file_loc, 'r') as file:
            for i, line in enumerate(file):
                seeds.append(int(line.strip()))
    # check how many seeds we have
    if len(seeds) < n_perm:
        # warn we are adding seeds if some were supplied
        if seeds_file_loc is not None:
            warnings.warn(' '.join([str(n_perm), 'requested, but only', str(len(seeds)), 'seeds supplied, will add more seeds\n'])) 
        # get as many seeds as there are permutations
        for i in range(len(seeds), n_perm):
            # create a seed from zero to max
            seeds.append(random.randint(0, 2**32 - 1))
    elif len(seeds) > n_perm:
        # warn there are too many seeds
        warnings.warn(' '.join([str(n_perm), 'permutations requested, but', str(len(seeds)), 'seeds supplied, will only use first', str(n_perm), 'seeds\n'])) 
        seeds = seeds[0 : (n_perm -1)]

    # put results in a list
    permuted_results = []
    # do each seed
    for seed in seeds:
        # get a permuted ATAC matrix
        permuted_subset_regions = shuffle_csr_columns(subset_regions, seed)
        # and do the actual analysis now
        regression_results_permuted = None
        # if the chunk size is bigger than zero, we do chunking
        if row_chunk_size > 0:
            regression_results_permuted = binomial_regression_chunked_sparse(permuted_subset_regions, rna_normal, chunk_size = row_chunk_size)
        # otherwise we do everything all in once
        else:
            regression_results_permuted = binomial_regression(permuted_subset_regions.todense(), rna_normal)
        # put in list
        permuted_results.append(regression_results_permuted)

    # now get the p, se and beta for each permutation
    perm_stats = []
    for perm_result in permuted_results:
        # add the results back in the order of the confinement file
        perm_ordered_ps = [None] * len(regions)
        perm_ordered_betas = [None] * len(regions)
        perm_ordered_ses = [None] * len(regions)
    
        # Place the results back in the original order
        for i in range(0, len(original_positions), 1):
            # get the resulting p
            perm_p = perm_result['p'][i]
            # and beta
            perm_beta = perm_result['beta'][i]
            # and se
            perm_std_err = perm_result['std_err'][i]
            # get position in original list
            i_original = original_positions[i]
            # put in the lists at those positions
            perm_ordered_ps[i_original] = perm_p
            perm_ordered_betas[i_original] = perm_beta
            perm_ordered_ses[i_original] = perm_std_err
        # add these to the perm stats
        perm_stats.append({'p' : perm_ordered_ps, 'beta' : perm_ordered_betas, 'std_err' : perm_ordered_ses})

    # open the file handle to the seeds file
    output_seeds_loc = ''.join([output_folder, '/used_seeds.txt.gz'])
    with gzip.open(output_seeds_loc, 'wt') as sf:
        # write each seed and its results
        for i in range(0, len(seeds)):
            # write the seed
            sf.write(f"{str(seeds[i])}\n")
            # now get the path to the permuted outputs
            perm_beta_loc = ''.join([output_folder, '/perm_', str(i), '_beta.txt.gz'])
            perm_p_loc = ''.join([output_folder, '/perm_', str(i), '_p.txt.gz'])
            perm_se_loc = ''.join([output_folder, '/perm_', str(i), '_se.txt.gz'])
            # write those permuted results
            with gzip.open(perm_p_loc, 'wt') as f:
                if sample_name is not None:
                    f.write(f"{sample_name}\n")
                for item in perm_stats[i]['p']:
                    f.write(f"{item}\n")
    
            with gzip.open(perm_beta_loc, 'wt') as f:
                if sample_name is not None:
                    f.write(f"{sample_name}\n")
                for item in perm_stats[i]['beta']:
                    f.write(f"{item}\n")
    
            with gzip.open(perm_se_loc, 'wt') as f:
                if sample_name is not None:
                    f.write(f"{sample_name}\n")
                for item in perm_stats[i]['std_err']:
                    f.write(f"{item}\n")
