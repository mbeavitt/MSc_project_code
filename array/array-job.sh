#!/bin/bash

################################################################################
#                                                                              #
# Array Job - process all files in a directory                                 #
#                                                                              #
################################################################################

# Grid Engine options
#
#$ -cwd
#
# Task range. Tasks need to go from 1 to the number of files in the target 
# directory.
#
#$ -l h_vmem=2G
#$ -t 1-48
#$ -e error
#$ -o output

# Target directory
TARGET_DIR=/exports/eddie/scratch/s1653324/transcriptome_data/raw_data/20220726

# Get list of files in target directory
files=$(ls -1 ${TARGET_DIR}/)

# Get file to be processed by *this* task 
# extract the Nth file in the list of files, $files, where N == $SGE_TASK_ID
this_file=$(echo "${files}" | sed -n ${SGE_TASK_ID}p)
echo Processing file: ${this_file} on $HOSTNAME

# Process file
$1 "${TARGET_DIR}/${this_file}"
