# Grid Engine options
#
#$ -wd /exports/eddie/scratch/s1653324/
#
# Task range. Tasks need to go from 1 to the number of files in the target 
# directory.
#
#$ -l h_vmem=10G
#$ -pe sharedmem 8
#$ -e error
#$ -o output

source ~/.bashrc
source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a

mkdir ${SCRATCH_SPACE}/genom_out

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ${SCRATCH_SPACE}/genom_out/ --genomeFastaFiles ${SCRATCH_SPACE}/analysis_genomes/combined.fa --sjdbGTFfile ${SCRATCH_SPACE}/analysis_genomes/combined.gtf
