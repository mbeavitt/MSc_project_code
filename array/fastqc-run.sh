# Set the target directory
target_dir="$1"

echo $target_dir

# Load the fastqc module
source /etc/profile.d/modules.sh
module load igmm/apps/FastQC/0.11.9

echo $(ls ${target_dir})

# Iterate through all .fastq.gz files in the target directory
for file in "${target_dir}"/2*.fastq.gz
do
    # Call fastqc on the file
    fastqc "$file" --threads 4 -o "/exports/eddie/scratch/s1653324/fastqc_output/pre-trim/"
done

