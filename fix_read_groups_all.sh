#!/bin/bash

# Define input and output directories
input_dir="linvilla"
output_dir="linvilla_fixed"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through BAM files and add unique read groups
for bamfile in "$input_dir"/*.bam; do
    filename=$(basename "$bamfile")
    
	# use all bam files
   sample="${filename%.bam}"
    output_bam="$output_dir/${filename%.bam}_RGfixed.bam"

    echo "Processing $filename -> $output_bam with sample name $sample"

    # Add read groups using samtools
    samtools addreplacerg -r "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:lib1" -o "$output_bam" "$bamfile"

    # Index the new BAM file
    samtools index "$output_bam"
done

echo "✅ Read group fixing complete!"
