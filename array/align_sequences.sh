#!/bin/bash

source ~/.bashrc
source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a

directory=$1

unique_part=$(for file in ${directory}/*.fastq; do echo $file | cut -f1 -d "." | cut -f8 -d "_"; done | uniq)

echo "Aligning sequence ${unique_part}"


STAR --genomeDir "$SCRATCH_SPACE/genom_out" \
     --runThreadN 4 \
     --readFilesIn ${directory}/*${unique_part}*_1.fastq ${directory}/*${unique_part}*_2.fastq \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --outSAMunmapped Within \
     --outFileNamePrefix "$directory/${unique_part}" \
     --outSAMtype BAM SortedByCoordinate

mv ${directory}/*.bam ${SCRATCH_SPACE}/alignment_output
