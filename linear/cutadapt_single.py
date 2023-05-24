#!/exports/applications/apps/SL7/anaconda/5.3.1/bin/python3

import sys
import subprocess

path = sys.argv[1]
identifier = sys.argv[2]

command = [
    "cutadapt",
    "--minimum-length", "20",
    "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "-o", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_1_" + identifier + "_1.fastq.gz",
    "-p", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_1_" + identifier + "_2.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_1_" + identifier + "_1.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_1_" + identifier + "_2.fastq.gz"
]

subprocess.run(command)


command2 = [
    "cutadapt",
    "--minimum-length", "20",
    "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "-o", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_2_" + identifier + "_1.fastq.gz",
    "-p", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_2_" + identifier + "_2.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_2_" + identifier + "_1.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_2_" + identifier + "_2.fastq.gz"
]

subprocess.run(command2)

