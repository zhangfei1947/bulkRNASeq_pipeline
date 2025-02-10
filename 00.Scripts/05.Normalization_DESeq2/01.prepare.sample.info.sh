#!/bin/bash

quanpath=/scratch/group/lilab/Phil/20250107_isogrp/04.Quantification

output_file="/scratch/group/lilab/Phil/20250107_isogrp/05.Normalization_DESeq2/sampleinfo.csv"

if [ -e "$output_file" ]; then
    echo "Error: $output_file already exists. Please remove it or choose a different output file."
    exit 1
fi

# Create the directory if it doesn't exist
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"

# Print the header
echo "sample	condition" > "$output_file"

# Find the directories, process them, and sort the output
find $quanpath -mindepth 1  -type d | while read -r dir; do
  # Extract the directory name without the path
  dir_name=$(basename "$dir")

  # Extract the base name (without postfix 1, 2, 3)
  base_name=$(echo "$dir_name" | sed 's/[123]$//')

  # Store the lines in an array for sorting
  printf "%s\t%s\n" "$dir_name" "$base_name"
done | sort > temp_output.tsv # Sort the output and write to temp file

cat temp_output.tsv >> "$output_file" #Append the sorted output to the output file
rm temp_output.tsv #Remove the temp file

echo "Output written to $output_file"
