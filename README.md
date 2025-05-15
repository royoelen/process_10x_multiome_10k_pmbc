# processing of public 10x multiome dataset of 10k PBMCs

This repo describes processing a public PBMC dataset, generated with 10x multiome.

## sofware used

This is the sofware used to process the data

Seurat (R)
Signac (R)
Azimuth (R)
Scanpy (python)
Pycistopic (python)
SCENIC+ (python)
cPeaks (python)
Cellbender (python)

## steps

These are the steps to process the data


### download the data

'*10k_download_data.sh*'    download the data from the 10x website


### expression data

'*10k_ambient_rna_correction.sh*    use Cellbender for ambient RNA correction\
'*10k_h5_to_scanpy.py*    convert the Cellbender corrected matrix to a Scanpy object to do QC\
'*10k_scanpy_to_10x_matrices.py*    export the raw scanpy object into Seurat compatible matrices\
'*10k_matrix_to_seurat.R*'  convert the matrices into a Seurat object\
'*10k_azimuth_reference_map.R*' do reference mapping of the expression data to a reference dataset for cell type annotation


### chromatin data

'*10k_round_fragments.py*'  round fragments coming from cellranger ARC, as there are paired reads, which should be counted as one\
'*10k_cpeaks_align.sh*' align the counted fragments on the cpeaks reference\
'*10k_cpeaks_to_signac.R*'  load the cpeaks counted fragments into Signac to do QC


### SCENIC+ pipeline

'*10k_pycistopic_modelling.ipynb*'  use Pycistopic to model topics based on the chromatin data, and find differentially open regions\
'*10k_scenicplus_config.yaml*'  run SCENIC+ using the pycistopic object, the scanpy object and the DARs/topics\
'*10k_scenic_add_region_info_to_cres.R*'    add info from the topics/dars/celltypes to the SCENIC+ output and merge direct and extended


### result
'*eRegulons_both.tsv.gz*'   gzipped version of eRegulon outputs for direct and extended
