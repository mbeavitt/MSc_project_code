#!/bin/bash

 echo "Loading modules..."
 module load igmm/apps/STAR/2.7.8a
 module load igmm/apps/cutadapt/1.16

target_dir="$1"

# echo "running cutadapt_single.py from code folder..."
# python3 ~/Code/RNAseq/linear/cutadapt_single.py

# echo "unzipping output..."
# gunzip trimmed*.gz

# "Cleaning up compressed files..."
# rm *.gz

# unique_part=$(for file in *.fastq; do echo $file | cut -f1 -d "." | cut -f7 -d "_";  done | uniq)

# Concatenate the R1 and R2 files from different lanes
#echo "Concatenating files..."
#cat trimmed_220721_A00291_0453_AH7JGYDRX2_*_${unique_part}_1.fastq > trimmed_merged_220721_A00291_0453_AH7JGYDRX2_${unique_part}_1.fastq
#cat trimmed_220721_A00291_0453_AH7JGYDRX2_*_${unique_part}_2.fastq > trimmed_merged_220721_A00291_0453_AH7JGYDRX2_${unique_part}_2.fastq
#rm trimmed_220721_A00291_0453_AH7JGYDRX2_*_${unique_part}_1.fastq
#rm trimmed_220721_A00291_0453_AH7JGYDRX2_*_${unique_part}_2.fastq

echo "Aligning sequence files to genome..."
for seqfile in *.fastq; do

    identifier=$(echo $seqfile | cut -f1 -d "." | cut -f8 -d "_")

    echo "Aligning sequence ${uniquepart}_${identifier}"

    STAR --genomeDir "$SEQUENCE_PATH/genomes/genom_out" \
	--runThreadN 4 \
	--readFilesIn $seqfile \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--outSAMunmapped Within \
	--outFileNamePrefix ${uniquepart}_${identifier} \
	--outSAMtype BAM SortedByCoordinate

done

echo "done!"


