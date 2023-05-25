# Set the target directory
target_dir="$1"

# Iterate through all .fastq.gz files in the target directory
for file in "${target_dir}"/trimmed*.fastq.gz
do
    # Extract the base filename without extension
    base_filename=$(basename "$file" .fastq.gz)
    
    # Call gunzip on the file and output to the target directory with the base filename
    gunzip -c "$file" > "${target_dir}/${base_filename}.fastq"
done

rm ${target_dir}/*.gz
