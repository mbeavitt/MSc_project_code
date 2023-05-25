#!/bin/bash

qsub -sync y /home/s1653324/Code/RNAseq/linear/staging_script.sh
qsub -sync y /home/s1653324/Code/RNAseq/linear/unpack_script.sh
