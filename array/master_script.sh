#!/bin/bash

echo "Loading modules..."
source ~/.bashrc
source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a
module load igmm/apps/cutadapt/1.16

#echo "Clearing scratch space..."
#rm -r /exports/eddie/scratch/s1653324/*

#echo "Staging data"
#qsub -sync y /home/s1653324/Code/RNAseq/linear/staging_script.sh

#echo "Unpacking sequence files"
#qsub -sync y /home/s1653324/Code/RNAseq/linear/unpack_script.sh

#echo "Removing archive"
#rm ${SCRATCH_SPACE}/archive.zip

#echo "Running fastqc... (first run)"
#mkdir -p ${SCRATCH_SPACE}/fastqc_output/pre-trim/
#mkdir -p ${SCRATCH_SPACE}/fastqc_output/post-trim/
#qsub -cwd -l h_vmem=4G -pe sharedmem 4 -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/fastqc-run.sh

#echo "Running cutadapt..."
#qsub -cwd -l h_vmem=1G -pe sharedmem 1 -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/cutadapt-run.sh

#echo "Running fastqc... (second run)"
#qsub -cwd -l h_vmem=4G -pe sharedmem 4 -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/fastqc-run_trimmed.sh

#echo "Unzipping files..."
#qsub -cwd -l h_vmem=2G -pe sharedmem 1 -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/un-gz.sh

#find /exports/eddie/scratch/s1653324/raw_data/20220726 -maxdepth 2 -type f -name "*.gz" -delete
#find /exports/eddie/scratch/s1653324/raw_data/20220726 -maxdepth 2 -type f -name "*.count" -delete

#echo "Concatenating lanes..."
#qsub -cwd -l h_vmem=2G -pe sharedmem 1 -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/cat_files.sh

#echo "Building genome..."
#mkdir -p ${SCRATCH_SPACE}/analysis_genomes
#cp /home/s1653324/rna_analysis_genomes/* ${SCRATCH_SPACE}/analysis_genomes
#qsub -sync y /home/s1653324/Code/RNAseq/linear/generate.sh

echo "Aligning sequence files to genome..."
mkdir -p ${SCRATCH_SPACE}/alignment_output
qsub -wd ${SCRATCH_SPACE} -l h_vmem=40G -pe sharedmem 1 -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/align_sequences.sh

echo "Done!"

