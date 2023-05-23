# Set the target directory
target_dir="$1"

source /etc/profile.d/modules.sh

# Iterate through all .fastq.gz files in the target directory
for file in "${target_dir}"/trimmed*.fastq.gz
do
    # Extract the base filename without extension
    base_filename=$(basename "$file" .fastq.gz)
    
    # Call gunzip on the file and output to specified directory with the base filename
    gunzip -c "$file" > "/exports/eddie/scratch/s1653324/transcriptome_data/raw_data/20220726/samples/${base_filename}.fastq"
done

