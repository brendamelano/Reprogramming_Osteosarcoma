#! /bin/env bash

#$ -S /bin/bash                     #-- the shell for the job
#$ -l h_rt=65:00:00
#$ -l scratch=65G                   #-- SGE resources (home and scratch disks)
#$ -q long.q
#$ -pe smp 4
#$ -l mem_free=70G                  #-- submits on nodes with enough free memory (required)    
#$ -o ~/perturb_seq/OS384
#$ -e ~/perturb_seq/OS384


cd ~/perturb_seq/OS384


# File paths
guide_file="guide_sequences.txt"
sequence_file="OS384_combined_guide_reads_filtered_2X.txt"
barcode_file="OS384_cell_barcodes.txt"
output_file="OS384_guide_cell_barcode_matches_counts.txt"


# Clear the output file
> "$output_file"


# Read cell barcodes and guides into arrays
mapfile -t barcodes < "$barcode_file"
mapfile -t guides < "$guide_file"


# Declare an associative array for counting unique guide-barcode pairs
declare -A guide_barcode_count


# Process each sequence to find matching guide and barcode pairs
while read -r sequence; do
    for guide in "${guides[@]}"; do
        if [[ "$sequence" == *"$guide"* ]]; then
            for barcode in "${barcodes[@]}"; do
                if [[ "$sequence" == *"$barcode"* ]]; then
                    # Increment the count for the guide-barcode pair
                    guide_barcode_count["$guide,$barcode"]=$((guide_barcode_count["$guide,$barcode"]+1))
                    break # Assuming each sequence contains at most one barcode
                fi
            done
            break # Assuming each sequence contains at most one guide
        fi
    done
done < "$sequence_file"

# Output the counts to the file
for key in "${!guide_barcode_count[@]}"; do
    IFS=',' read -r guide barcode <<< "$key"
    echo -e "$barcode\t$guide\t${guide_barcode_count[$key]}" >> "$output_file"
done

# Sort the output file for better readability (optional)
sort "$output_file" -o "$output_file"
