# Set the target directory
target_dir="$1"

# Load the fastqc module
source /etc/profile.d/modules.sh
module load igmm/apps/cutadapt/1.16

python3 ~/Code/RNAseq/array/cutadapt_script.py ${target_dir}
