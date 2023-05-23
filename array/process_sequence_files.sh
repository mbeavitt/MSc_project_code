#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

directory="$1"

echo "Loading modules..."
module load igmm/apps/STAR/2.7.8a
module load igmm/apps/cutadapt/1.16

echo "running cutadapt_single.py from code folder..."
python3 ~/Code/RNAseq/linear/cutadapt_single.py

echo "unzipping output..."
gunzip "$directory"/trimmed*.gz

# "Cleaning up compressed files..."
rm "$directory"/*.gz

unique_part=$(for file in "$directory"/*.fastq; do echo $file | cut -f1 -d "." | cut -f7 -d "_"; done | uniq)

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

echo "done!"

