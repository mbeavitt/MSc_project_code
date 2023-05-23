# Set the target directory
target_dir="$1"

# Load the fastqc module
source /etc/profile.d/modules.sh
module load igmm/apps/FastQC/0.11.9

# Iterate through all .fastq.gz files in the target directory
for file in "${target_dir}"/trimmed*.fastq.gz
do
    # Call fastqc on the file
    fastqc "$file" --threads 4 -o "/exports/eddie/scratch/s1653324/transcriptome_data/raw_data/trimmed_fastqc_output"
done

