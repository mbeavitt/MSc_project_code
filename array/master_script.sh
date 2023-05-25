#!/bin/bash

echo "Loading modules..."
source ~/.bashrc
source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a
module load igmm/apps/cutadapt/1.16

echo "Clearing scratch space..."
rm -r /exports/eddie/scratch/s1653324/*

echo "Staging data"
qsub -sync y /home/s1653324/Code/RNAseq/linear/staging_script.sh

echo "Unpacking sequence files"
qsub -sync y /home/s1653324/Code/RNAseq/linear/unpack_script.sh

echo "Removing archive"
rm archive.zip

echo "Running fastqc... (first run)"
mkdir -p ${SCRATCH_SPACE}/fastqc_output/pre-trim/
mkdir -p ${SCRATCH_SPACE}/fastqc_output/post-trim/
qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/fastqc-run.sh "4G" 4

echo "Running cutadapt..."
qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/cutadapt-run.sh "4G" 1

echo "Running fastqc... (second run)"
qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/fastqc-run.sh "4G" 4

echo "Unzipping files..."
qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/un-gz.sh "200M" 1

echo "Concatenating lanes..."
qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/cat_files.sh "200M" 1

echo "Building genome..."
qsub -sync y /home/s1653324/Code/RNAseq/linear/generate.sh

echo "Aligning sequence files to genome..."
mkdir -p ${SCRATCH_SPACE}/analysis_genomes
mv /home/s1653324/rna_analysis_genomes/* ${SCRATCH_SPACE}/analysis_genomes
mkdir -p ${SCRATCH_SPACE}/alignment_output
qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/align_sequences.sh "40G" 1

echo "Done!"
