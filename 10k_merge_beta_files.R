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


####################
# Settings         #
####################


####################
# Main code        #
####################

# location of the betas
betas_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/naive_cres/betas.tsv.gz'
ps_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/naive_cres/ps.tsv.gz'
ses_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/naive_cres/ses.tsv.gz'

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
merged_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/10x_multiome_10k_pbmcs/naive_cres/10k_merged.tsv.gz'
# and save the result
write.table(all, gzfile(merged_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# make checksum
mdfiver::create_md5_for_file(merged_loc)
