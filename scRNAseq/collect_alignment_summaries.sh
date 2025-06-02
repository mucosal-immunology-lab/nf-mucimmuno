#!/bin/bash

# Define the output file
output_file="AlignmentSummary.csv"

# Clear the output file if it exists
> "$output_file"

# Print header to the output file
echo "Sample,Number of Reads,Reads With Valid Barcodes,Estimated Number of Cells,Reads Mapped to Genome: Unique" >> "$output_file"

# Loop through each sample folder in the STARsolo directory
for sample_folder in results/STARsolo/*; do
    if [ -d "$sample_folder" ]; then
        sample_name=$(basename "$sample_folder")
        summary_file="$sample_folder/$sample_name.Solo.out/GeneFull_Ex50pAS/Summary.csv"

        if [ -f "$summary_file" ]; then
            # Extract relevant values using awk
            number_of_reads=$(awk -F, '$1=="Number of Reads"{print $2}' "$summary_file")
            valid_barcodes=$(awk -F, '$1=="Reads With Valid Barcodes"{print $2}' "$summary_file")
            estimated_cells=$(awk -F, '$1=="Estimated Number of Cells"{print $2}' "$summary_file")
            unique_reads=$(awk -F, '$1=="Reads Mapped to Genome: Unique"{print $2}' "$summary_file")

            # Append results in CSV format
            echo "$sample_name,$number_of_reads,$valid_barcodes,$estimated_cells,$unique_reads" >> "$output_file"

        else
            echo "Warning: Summary file not found for sample $sample_name"
        fi
    fi
done

echo "Alignment summary has been written to $output_file."
