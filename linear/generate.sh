# Grid Engine options
#
#$ -cwd
#
# Task range. Tasks need to go from 1 to the number of files in the target 
# directory.
#
#$ -l h_vmem=10G
#$ -pe sharedmem 8
#$ -e error
#$ -o output

source /etc/profile.d/modules.sh
module load igmm/apps/STAR/2.7.8a

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir genom_out/ --genomeFastaFiles combined.fa --sjdbGTFfile combined.gtf
