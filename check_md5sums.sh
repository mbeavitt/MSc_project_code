#!/bin/bash
# Written by Michael Beavitt, aided by GPT-4

# Prompt: how do I use a md5sums.txt file to check the integrity of downloaded and unzipped data?

# Usage: ./check_md5sums.sh md5sums.txt

md5sums_file="$1"

# Check if the md5sums file is provided
if [ -z "$md5sums_file" ]; then
    echo "Usage: $0 <md5sums.txt>"
    exit 1
fi

# Check if the md5sums file exists
if [ ! -f "$md5sums_file" ]; then
    echo "Error: md5sums file not found."
    exit 1
fi

# Loop through the md5sums file and check the MD5 hashes of the files
all_files_valid=true
while IFS= read -r line; do
    expected_md5=$(echo "$line" | awk '{print $1}')
    file_path=$(echo "$line" | awk '{print $2}')

    if [ ! -f "$file_path" ]; then
        echo "Error: File not found: $file_path"
        all_files_valid=false
        continue
    fi

    actual_md5=$(md5sum "$file_path" | awk '{print $1}')

    if [ "$expected_md5" == "$actual_md5" ]; then
        echo "Valid: $file_path"
    else
        echo "Invalid: $file_path"
        echo "Expected MD5: $expected_md5"
        echo "Actual MD5: $actual_md5"
        all_files_valid=false
    fi
done < "$md5sums_file"

if [ "$all_files_valid" = true ]; then
    echo "All files are valid."
else
    echo "Some files are invalid."
    exit 1
fi
