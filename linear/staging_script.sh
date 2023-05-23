#!/bin/bash
#
# Grid Engine options
#
#$ -cwd
#
# Task range. Tasks need to go from 1 to the number of files in the target 
# directory.
#
#$ -l h_vmem=8G
#$ -q staging
#$ -e staging_error
#$ -o staging_output

rsync "/exports/csce/eddie/biology/groups/mooney/Michael Beavitt/RNAseq Data - 19 Jan/archive.zip" "/exports/eddie/scratch/s1653324/"



