#!/bin/bash

unzip -l $SCRATCH_SPACE/archive.zip > $SCRATCH_SPACE/files.txt

path_parse.py $SCRATCH_SPACE/files.txt

while IFS= read -r path
do
    unzip $SCRATCH_SPACE/archive.zip $path -d $SCRATCH_SPACE
done < "paths.txt"

