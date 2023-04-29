#!/bin/bash

# Change the following path to the root directory containing the directories with fastq.gz files
ROOT_DIR="/exports/eddie/scratch/s1653324/transcriptome_data/unzipped_files/raw_data/20220726/"

# Change the following path to the FastQC executable if it's not in the PATH
FASTQC_EXEC="fastqc"

# Change the following path to the directory where you want to save the FastQC output
OUTPUT_DIR="/exports/eddie/scratch/s1653324/transcriptome_data/unzipped_files/raw_data/20220726/fastqc_output"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Find directories containing fastq.gz files
DIRS_WITH_FASTQ=$(find "${ROOT_DIR}" -type f -name "*.fastq.gz" -exec dirname {} \; | sort -u)

echo DIRS_WITH_FASTQ

# Loop through the directories containing fastq.gz files
for current_dir in ${DIRS_WITH_FASTQ}; do
    # Find all fastq.gz files in the current directory
    FASTQ_FILES=$(find "${current_dir}" -type f -name "*.fastq.gz")
    
    # Run FastQC on each fastq.gz file
    for fastq_file in ${FASTQ_FILES}; do
        echo "Running FastQC on ${fastq_file}"
        "${FASTQC_EXEC}" -o "${OUTPUT_DIR}" "${fastq_file}"
    done
done

echo "FastQC analysis completed."

