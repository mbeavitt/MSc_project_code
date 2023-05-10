# Set the target directory
target_dir="$1"

# Load the fastqc module
source /etc/profile.d/modules.sh
module load igmm/apps/cutadapt/1.9.1

# Iterate through all .fastq.gz files in the target directory
for file in "${target_dir}"/*.fastq.gz
do
    # Call fastqc on the file
    fastqc "$file" --threads 4 -o "/exports/eddie/scratch/s1653324/transcriptome_data/unzipped_files/raw_data/20220726/fastqc_output"
done

