#!/exports/applications/apps/SL7/anaconda/5.3.1/bin/python3

import sys
import subprocess


path = "/exports/eddie/scratch/s1653324/transcriptome_data/raw_data/20220726/test_folder"

command = [
    "cutadapt",
    "--minimum-length", "20",
    "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "-o", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_1_18410RE0020L01_1.fastq.gz",
    "-p", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_1_18410RE0020L01_2.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_1_18410RE0020L01_1.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_1_18410RE0020L01_2.fastq.gz"
]

subprocess.run(command)


command2 = [ 
    "cutadapt",
    "--minimum-length", "20",
    "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "-o", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_2_18410RE0020L01_1.fastq.gz",
    "-p", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_2_18410RE0020L01_2.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_2_18410RE0020L01_1.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_2_18410RE0020L01_2.fastq.gz"
]

subprocess.run(command2)

