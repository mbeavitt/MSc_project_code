#!/bin/bash

directory=$1
unique_part=$(for file in ${directory}/*.fastq; do echo $file | cut -f1 -d "." | cut -f6 -d "_"; done | uniq)

cat "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_1.fastq > "$directory"/trimmed_merged_220721_A00291_0453_AH7JGYDRX2_"$unique_part"_1.fastq
cat "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_2.fastq > "$directory"/trimmed_merged_220721_A00291_0453_AH7JGYDRX2_"$unique_part"_2.fastq
rm "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_1.fastq
rm "$directory"/trimmed_220721_A00291_0453_AH7JGYDRX2_*_"$unique_part"_2.fastq

