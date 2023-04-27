#!/bin/bash


# Written by Michael Beavitt aided by GPT-4

# Prompt:
# Can you write me a looping bash script that will take file paths from a text file containing one path per line,
# and extract each of those files from a zip archive? (the first command might be e.g. 'unzip archive.zip file1.csv')

# Usage: ./extract_files.sh archive.zip file_paths.txt

archive="$1"
file_paths="$2"

# Check if both arguments are provided
if [ -z "$archive" ] || [ -z "$file_paths" ]; then
    echo "Usage: $0 <archive.zip> <file_paths.txt>"
    exit 1
fi

# Check if the archive and file paths list exist
if [ ! -f "$archive" ] || [ ! -f "$file_paths" ]; then
    echo "Error: Archive or file paths list not found."
    exit 1
fi

# Loop through the file paths and extract each file from the archive
while IFS= read -r file_path; do
    echo "Extracting $file_path from $archive..."
    unzip "$archive" "$file_path"
done < "$file_paths"

