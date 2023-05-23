#!/exports/applications/apps/SL7/anaconda/5.3.1/bin/python3
'''Usage: this_script.py <directory containing fastq.gz files>'''

import sys
import subprocess

path = sys.argv[1]

# Run the OS command and capture its output
output = subprocess.check_output(["ls", "-l", path])

# Decode the output from bytes to a string
output = output.decode("utf-8").strip().split('\n')[1:]

dirs = []

for i in output:
    filename = i.split()[-1]
    if filename.startswith('220721') and filename.endswith('.gz'):
        unique_part = filename.split('_')[-2]

command = [
    "cutadapt",
    "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "-o", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_1_" + unique_part + "_1.fastq.gz",
    "-p", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_1_" + unique_part + "_2.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_1_" + unique_part + "_1.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_1_" + unique_part + "_2.fastq.gz"
]

subprocess.run(command)


command2 = [ 
    "cutadapt",
    "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "-o", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_2_" + unique_part + "_1.fastq.gz",
    "-p", path + "/trimmed_220721_A00291_0453_AH7JGYDRX2_2_" + unique_part + "_2.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_2_" + unique_part + "_1.fastq.gz",
    path + "/220721_A00291_0453_AH7JGYDRX2_2_" + unique_part + "_2.fastq.gz"
]

subprocess.run(command2)

