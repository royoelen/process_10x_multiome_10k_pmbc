#!/usr/bin/env bash
###################################################################
#Script Name	  : 10k_aggregated_creqtl_inputs.sh
#Description	  : aggregate all the CRE QTL inputs
#Args           : cre location, folder containing per-sample p/beta/se, folder to output aggregated matrices, optional number of permutations to check for
#Author       	: Roy Oelen
#example        : 10k_aggregated_creqtl_inputs.sh \
#                 /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/confinements/region_to_peak_variant_overlaps.tsv.gz \
#                 /groups/umcg-franke-scrna/tmp02/external_datasets/10x_multiome_10k_pbmcs/naive_cres/ \
#                 /groups/umcg-franke-scrna/tmp02/external_datasets/10x_multiome_10k_pbmcs/naive_cres/ \
#                 10
###################################################################

# the parameters supplied
CRE_LOC=$1
SUMMARIES_LOC=$2
OUTPUT_LOC=$3
N_PERM=$4
# default permutations is 0
if [ -z "${4}" ];
    then
    N_PERM='0'
fi

# generate random UUID for the temporary file
uuid=$(uuidgen)
# generate temporary beta/p/se files
tmp_ps=${OUTPUT_LOC}'/'${uuid}'.ps.tsv'
tmp_betas=${OUTPUT_LOC}'/'${uuid}'.betas.tsv'
tmp_ses=${OUTPUT_LOC}'/'${uuid}'.ses.tsv'
# and full files
ps_full=${OUTPUT_LOC}'/ps.tsv'
betas_full=${OUTPUT_LOC}'/betas.tsv'
ses_full=${OUTPUT_LOC}'/ses.tsv'

# keep track of the columns
aggregated_columns_loc=${OUTPUT_LOC}'/starting_columns.tsv'

# create the base file
echo -en "region\tgene" > ${aggregated_columns_loc}

# get all directories
dirlist=(${SUMMARIES_LOC}*)

# init the tables to put the outputs is
zcat ${CRE_LOC} > ${ps_full}
zcat ${CRE_LOC} > ${betas_full}
zcat ${CRE_LOC} > ${ses_full}

# make a list for the permutations
perm_numbers=[]
# to be filled if there were permutations
if [ "$N_PERM" != "0" ]; then
    # permutations were zero-based, so we need to stop one earlier
    perm_stop="$((N_PERM-1))"
    # use seq to put them in the array
    perm_numbers=($(seq 0 1 ${perm_stop}))
fi
# initialize all the permutation files
for perm in ${perm_numbers[@]}; do
    zcat ${CRE_LOC} > ${OUTPUT_LOC}'/perm_'${perm}'_betas.tsv'
    zcat ${CRE_LOC} > ${OUTPUT_LOC}'/perm_'${perm}'_ses.tsv'
    zcat ${CRE_LOC} > ${OUTPUT_LOC}'/perm_'${perm}'_ps.tsv'
    echo -en "region\tgene" > ${OUTPUT_LOC}'perm_'${perm}'_starting_columns.tsv'
done

# check each directory
for sum_dir in ${dirlist[@]}; do
    # get base name of directory, which is sample
    sample_name=$(basename "$sum_dir")
    # check if it is a directory
    sample_dir=${SUMMARIES_LOC}'/'${sample_name}
    if [ -d "$sample_dir" ];
        then
        # paste the full path to the betas file
        betas_file=${sum_dir}'/beta.txt.gz'
        # next, check if it has output, by searching for the betas file
        if [ -e "$betas_file" ];
            then
            # add to the header if we have this file
            echo -en '\t'${sample_name} >> ${aggregated_columns_loc}
            # do the same for the betas and ps
            ps_file=${sum_dir}'/p.txt.gz'
            ps_file_extracted=${sum_dir}'/p.txt'
            # unzip file
            gunzip -c ${ps_file} > ${ps_file_extracted}
            # paste into temporary file
            paste ${ps_full} ${ps_file_extracted} > ${tmp_ps}
            # remove old full file
            rm ${ps_full}
            # make tmp file new file
            mv ${tmp_ps} ${ps_full}
            # remove unzipped file
            rm ${ps_file_extracted}
           
            # do the same for ses
            ses_file=${sum_dir}'/se.txt.gz'
            ses_file_extracted=${sum_dir}'/se.txt'
            gunzip -c ${ses_file} > ${ses_file_extracted}
            paste ${ses_full} ${ses_file_extracted} > ${tmp_ses}
            rm ${ses_full}
            mv ${tmp_ses} ${ses_full}
            rm ${ses_file_extracted}

            # and finally the betas
            betas_file_extracted=${sum_dir}'/beta.txt'
            gunzip -c ${betas_file} > ${betas_file_extracted}
            paste ${betas_full} ${betas_file_extracted} > ${tmp_betas}
            rm ${betas_full}
            mv ${tmp_betas} ${betas_full}
            rm ${betas_file_extracted}

            # now try to do the same for the permutations
            for perm in ${perm_numbers[@]}; do
                # paste the full path to the betas file
                perm_betas_file=${sum_dir}'/perm_'${perm}'_beta.txt.gz'
                # next, check if it has output, by searching for the betas file
                if [ -e "$perm_betas_file" ]; then
                    # add to the header if we have this file
                    echo -en '\t'${sample_name} >> ${OUTPUT_LOC}'/perm_'${perm}'_starting_columns.tsv'
                    
                    # do the same for the betas and ps
                    perm_ps_file=${sum_dir}'/perm_'${perm}'_p.txt.gz'
                    perm_ps_file_extracted=${sum_dir}'/perm_'${perm}'_p.txt'
                    # unzip file
                    gunzip -c ${perm_ps_file} > ${perm_ps_file_extracted}
                    # generate temporary beta/p/se files
                    perm_tmp_ps=${OUTPUT_LOC}'/'${uuid}'.perm_'${perm}'.ps.tsv'
                    perm_tmp_betas=${OUTPUT_LOC}'/'${uuid}'.perm_'${perm}'.betas.tsv'
                    perm_tmp_ses=${OUTPUT_LOC}'/'${uuid}'.perm_'${perm}'.ses.tsv'
                    # and the full files
                    perm_ps_full=${OUTPUT_LOC}'/perm_'${perm}'_ps.tsv'
                    perm_betas_full=${OUTPUT_LOC}'/perm_'${perm}'_betas.tsv'
                    perm_ses_full=${OUTPUT_LOC}'/perm_'${perm}'_ses.tsv'
                    
                    # paste into temporary file
                    paste ${perm_ps_full} ${perm_ps_file_extracted} > ${perm_tmp_ps}
                    # remove old full file
                    rm ${perm_ps_full}
                    # make tmp file new file
                    mv ${perm_tmp_ps} ${perm_ps_full}
                    # remove unzipped file
                    rm ${perm_ps_file_extracted}

                    # do the same for ses
                    perm_ses_file=${sum_dir}'/perm_'${perm}'_se.txt.gz'
                    perm_ses_file_extracted=${sum_dir}'/perm_'${perm}'_se.txt'
                    gunzip -c ${perm_ses_file} > ${perm_ses_file_extracted}
                    paste ${perm_ses_full} ${perm_ses_file_extracted} > ${perm_tmp_ses}
                    rm ${perm_ses_full}
                    mv ${perm_tmp_ses} ${perm_ses_full}
                    rm ${perm_ses_file_extracted}

                    # and beta finally
                    perm_betas_file_extracted=${sum_dir}'/perm_'${perm}'_beta.txt'
                    gunzip -c ${perm_betas_file} > ${perm_betas_file_extracted}
                    paste ${perm_betas_full} ${perm_betas_file_extracted} > ${perm_tmp_betas}
                    rm ${perm_betas_full}
                    mv ${perm_tmp_betas} ${perm_betas_full}
                    rm ${perm_betas_file_extracted}
                fi
            done
        fi
    fi
done

# add newline to the headers
echo -en '\n' >> ${aggregated_columns_loc}

# also add the headers now
cat ${aggregated_columns_loc} > ${tmp_ps}
# and the content
cat ${ps_full} >> ${tmp_ps}
# remove old full file
rm ${ps_full}
# and make the temporary file the new one
mv ${tmp_ps} ${ps_full}

# same for the betas
cat ${aggregated_columns_loc} > ${tmp_betas}
cat ${betas_full} >> ${tmp_betas}
rm ${betas_full}
mv ${tmp_betas} ${betas_full}

# and the ses
cat ${aggregated_columns_loc} > ${tmp_ses}
cat ${ses_full} >> ${tmp_ses}
rm ${ses_full}
mv ${tmp_ses} ${ses_full}

# finally zip them all
bgzip ${ps_full}
bgzip ${betas_full}
bgzip ${ses_full}

# and clean up
rm ${aggregated_columns_loc}

# permutations as well
for perm in ${perm_numbers[@]}; do
# generate temporary beta/p/se files    
    perm_tmp_ps=${OUTPUT_LOC}'/'${uuid}'.perm_'${perm}'.ps.tsv'
    perm_tmp_betas=${OUTPUT_LOC}'/'${uuid}'.perm_'${perm}'.betas.tsv'
    perm_tmp_ses=${OUTPUT_LOC}'/'${uuid}'.perm_'${perm}'.ses.tsv'
    # and the full files
    perm_ps_full=${OUTPUT_LOC}'/perm_'${perm}'_ps.tsv'
    perm_betas_full=${OUTPUT_LOC}'/perm_'${perm}'_betas.tsv'
    perm_ses_full=${OUTPUT_LOC}'/perm_'${perm}'_ses.tsv'
    # add newline to the headers
    echo -en '\n' >> ${OUTPUT_LOC}'perm_'${perm}'_starting_columns.tsv'
    # also add the headers now
    cat ${OUTPUT_LOC}'perm_'${perm}'_starting_columns.tsv' > ${perm_tmp_ps}
    # and the content
    cat ${perm_ps_full} >> ${perm_tmp_ps}
    # remove old full file
    rm ${perm_ps_full}
    # and make the temporary file the new one
    mv ${perm_tmp_ps} ${perm_ps_full}

    # same for the betas
    cat ${OUTPUT_LOC}'perm_'${perm}'_starting_columns.tsv' > ${perm_tmp_betas}
    cat ${perm_betas_full} >> ${perm_tmp_betas}
    rm ${perm_betas_full}
    mv ${perm_tmp_betas} ${perm_betas_full}

    # and the ses
    cat ${OUTPUT_LOC}'perm_'${perm}'_starting_columns.tsv' > ${perm_tmp_ses}
    cat ${perm_ses_full} >> ${perm_tmp_ses}
    rm ${perm_ses_full}
    mv ${perm_tmp_ses} ${perm_ses_full}

    # finally zip them all
    bgzip ${perm_ps_full}
    bgzip ${perm_betas_full}
    bgzip ${perm_ses_full}

    # and clean
    rm ${OUTPUT_LOC}'perm_'${perm}'_starting_columns.tsv'

done