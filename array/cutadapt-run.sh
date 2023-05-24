# Set the target directory
path="$1"

# Load the cutadapt module
source /etc/profile.d/modules.sh
module load igmm/apps/cutadapt/1.16

# Grab identifier
unique_part=$(for file in "$path"/*.fastq.gz; do echo $file | cut -f1 -d "." | cut -f6 -d "_"; done | uniq)

# First command
cutadapt \
--minimum-length 20 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ${path}/trimmed_220721_A00291_0453_AH7JGYDRX2_1_${unique_part}_1.fastq.gz \
-p ${path}/trimmed_220721_A00291_0453_AH7JGYDRX2_1_${unique_part}_2.fastq.gz \
${path}/220721_A00291_0453_AH7JGYDRX2_1_${unique_part}_1.fastq.gz \
${path}/220721_A00291_0453_AH7JGYDRX2_1_${unique_part}_2.fastq.gz

# Second command
cutadapt \
--minimum-length 20 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ${path}/trimmed_220721_A00291_0453_AH7JGYDRX2_2_${unique_part}_1.fastq.gz \
-p ${path}/trimmed_220721_A00291_0453_AH7JGYDRX2_2_${unique_part}_2.fastq.gz \
${path}/220721_A00291_0453_AH7JGYDRX2_2_${unique_part}_1.fastq.gz \
${path}/220721_A00291_0453_AH7JGYDRX2_2_${unique_part}_2.fastq.gz
