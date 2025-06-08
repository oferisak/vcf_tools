#!/usr/bin/env bash

# Ask the user for the folder and number of files
read -p "Enter directory to search for .vcf.gz files: " DIR
read -p "How many files would you like to select? " X

# Find full paths, sort, take the first X, and write to selected_vcfs.txt
find "$DIR" -maxdepth 1 -type f -name '*.vcf.gz' \
  | sort \
  | head -n "$X" \
  > selected_vcfs.txt

echo "Wrote the first $X .vcf.gz file paths from $DIR to selected_vcfs.txt"
