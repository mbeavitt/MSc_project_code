#!/bin/bash

echo "Loading modules..."
source ~/.bashrc
source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a
module load igmm/apps/cutadapt/1.16

# Clear the scratch space
rm -r /exports/eddie/scratch/s1653324/*

mkdir /exports/eddie/scratch/s1653324/transcriptome_data

qsub -sync y /home/s1653324/Code/RNAseq/linear/staging_script.sh

if [ $# -ne 1 ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

directory="$1"

unique_part=$(for file in "$directory"/*.fastq.gz; do echo $file | cut -f1 -d "." | cut -f6 -d "_"; done | uniq)

mkdir -p /exports/eddie/scratch/s1653324/fastqc_output/pre-trim/
mkdir -p /exports/eddie/scratch/s1653324/fastqc_output/post-trim/

qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/fastqc-run.sh 4 4

echo "running cutadapt..."

qsub -sync y /home/s1653324/Code/RNAseq/array/array-job.sh /home/s1653324/Code/RNAseq/array/cutadapt-run.sh 4 1

echo "unzipping output..."
gunzip "$directory"/trimmed*.gz

# "Cleaning up compressed files..."
rm "$directory"/*.gz

# Concatenate the R1 and R2 files from different lanes
echo "Concatenating files..."
cat "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_1.fastq > "$directory"/trimmed_merged_220721_A00291_0453_AH7JGYDRX2_"$unique_part"_1.fastq
cat "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_2.fastq > "$directory"/trimmed_merged_220721_A00291_0453_AH7JGYDRX2_"$unique_part"_2.fastq
rm "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_1.fastq
rm "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_2.fastq

echo "Aligning sequence files to genome..."
for seqfile in "$directory"/*.fastq; do
    filename=$(basename "$seqfile")
    identifier=$(echo "$filename" | cut -f1 -d "." | cut -f8 -d "_")

    echo "Aligning sequence ${unique_part}_${identifier}"

    STAR --genomeDir "$SEQUENCE_PATH/genomes/genom_out" \
        --runThreadN 4 \
        --readFilesIn "$seqfile" \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outSAMunmapped Within \
        --outFileNamePrefix "$directory/${unique_part}_${identifier}" \
        --outSAMtype BAM SortedByCoordinate
done

mkdir -p "$SEQUENCE_PATH"/aligned_reads

echo "Moving bam files.."
mv "$directory"/*.bam "$SEQUENCE_PATH"/aligned_reads

echo "done!"

