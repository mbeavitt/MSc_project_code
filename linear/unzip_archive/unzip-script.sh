#!/bin/bash

unzip -l archive.zip > files.txt

path_parse.py files.txt

while IFS= read -r path
do
    unzip archive.zip $path
done < "paths.txt"

