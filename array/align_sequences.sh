#!/bin/bash

source ~/.bashrc
source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a

directory=$1

for seqfile in "$directory"/*.fastq; do
    filename=$(basename "$seqfile")
    identifier=$(echo "$filename" | cut -f1 -d "." | cut -f8 -d "_")
    unique_part=$(for file in ${directory}/*.fastq; do echo $file | cut -f1 -d "." | cut -f7 -d "_"; done | uniq)
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

mv ${directory}/*.bam ${SCRATCH_SPACE}/alignment_output
