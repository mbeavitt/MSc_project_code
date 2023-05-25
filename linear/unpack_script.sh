#!/bin/bash
#
# Grid Engine options
#
#$ -cwd
#
# Task range. Tasks need to go from 1 to the number of files in the target
# directory.
#
#$ -l h_vmem=16G
#$ -e unpack_error
#$ -o unpack_output

source ~/.bashrc

unzip -l ${SCRATCH_SPACE}/archive.zip > ${SCRATCH_SPACE}/files.txt

# Counting total lines in the file
total_lines=$(wc -l < ${SCRATCH_SPACE}/files.txt)

# Calculating the lines to be read
lines_to_read=$(($total_lines - 3))

# Reading file content except the first 4 lines and the last 2 lines
content=$(tail -n $lines_to_read ${SCRATCH_SPACE}/files.txt | head -n -3)

paths=()

# Looping over the content and appending 4th column to paths array
while read -r line; do
    path=$(echo $line | awk '{print $4}')
    if [[ $path == raw_data/20220726/B* ]]; then
        paths+=($path)
    fi
done <<< "$content"

# Writing paths array to paths.txt
printf "%s\n" "${paths[@]}" > ${SCRATCH_SPACE}/paths.txt

while IFS= read -r path
do
    unzip $SCRATCH_SPACE/archive.zip $path -d $SCRATCH_SPACE
done < "${SCRATCH_SPACE}/paths.txt"


# Correcting sample label mistake
mkdir ${SCRATCH_SPACE}/raw_data/20220726/B195

mv ${SCRATCH_SPACE}/raw_data/20220726/B165/*18410RE0009L01* ${SCRATCH_SPACE}/raw_data/20220726/B195

# Cleaning up
rm ${SCRATCH_SPACE}/paths.txt
rm ${SCRATCH_SPACE}/files.txt
