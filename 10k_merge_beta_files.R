#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 10k_merge_beta_files.R
# Function: merge the beta calculation files
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(qvalue)
library(mdfiver)

####################
# Functions        #
####################


merge_single_sample_files <- function(base_loc, betas_file='betas.tsv.gz', ses_file='ses.tsv.gz', ps_file='ps.tsv.gz', out_file='10k_merged.tsv.gz') {
  # location of the betas
  betas_loc <- paste(base_loc, betas_file, sep = '/')
  ps_loc <- paste(base_loc, ps_file, sep = '/')
  ses_loc <- paste(base_loc, ses_file, sep = '/')
  
  # load the data
  betas <- fread(betas_loc, header = T, sep = '\t', na.strings = c('NA', 'NaN', 'None'))
  ps <- fread(ps_loc, header = T, sep = '\t', na.strings = c('NA', 'NaN', 'None'))
  ses <- fread(ses_loc, header = T, sep = '\t', na.strings = c('NA', 'NaN', 'None'))
  
  # we only have one sample, so we can make the sample column the name of the table instead
  colnames(betas) <- c('region', 'gene', 'beta')
  colnames(ps) <- c('region', 'gene', 'p')
  colnames(ses) <- c('region', 'gene', 'se')
  
  # now we'll merge it all
  all <- merge(betas, ses, by.x = c('region', 'gene'), by.y = c('region', 'gene'))
  all <- merge(all, ps, by.x = c('region', 'gene'), by.y = c('region', 'gene'))
  
  # we'll add multiple testing correction
  all[['qvalue']] <- qvalue(all[['p']])$qvalue
  all[['bonferroni']] <- p.adjust(all[['p']], method = 'bonferroni')
  all[['BH']] <- p.adjust(all[['p']], method = 'BH')
  
  # where to save the results
  merged_loc <- paste(base_loc, out_file, sep = '/')
  # and save the result
  write.table(all, gzfile(merged_loc), row.names = F, col.names = T, sep = '\t', quote = F)
  # make checksum
  mdfiver::create_md5_for_file(merged_loc)
  return(0)
}


merge_single_sample_files_per_folder <- function(folders_loc, folders=NULL, betas_file='betas.tsv.gz', ses_file='ses.tsv.gz', ps_file='ps.tsv.gz', out_file='10k_merged.tsv.gz') {
  # list all folders
  folders_to_process <- list.dirs(folders_loc, full.names = F, recursive = F)
  # subset to specific folders if requested
  if (!is.null(folders)) {
    folders_to_process <- intersect(folders_to_process, folders)
  }
  # do each folder
  for (in_folder in folders_to_process) {
    # do each folder
    merge_single_sample_files(
      base_loc = paste0(folders_loc, '/', in_folder, '/'), 
      betas_file = betas_file, 
      ses_file = ses_file, 
      ps_file = ps_file, 
      out_file = out_file
    )
  }
  return(0)
}


####################
# Settings         #
####################


####################
# Main code        #
####################

# location of the betas
betas_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/naive_cres/'

# do each folder
merge_single_sample_files_per_folder(betas_loc, folders=c('B', 'CD4_T', 'CD8_T', 'DC', 'Mono', 'NK'))
