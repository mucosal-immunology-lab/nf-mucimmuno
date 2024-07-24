#!/bin/bash

# Define the output file
output_file="AlignmentSummary.txt"

# Clear the output file if it exists
> "$output_file"

# Loop through each sample folder in the STARsolo directory
for sample_folder in results/STARsolo/*; do
    if [ -d "$sample_folder" ]; then
        sample_name=$(basename "$sample_folder")
        summary_file="$sample_folder/$sample_name.Solo.out/GeneFull_Ex50pAS/Summary.csv"

        if [ -f "$summary_file" ]; then
            # Add the sample name to the output file
            echo "Sample: $sample_name" >> "$output_file"

            # Add the first two lines of the Summary.csv file to the output file
            head -n 2 "$summary_file" >> "$output_file"

            # Extract and add the estimated number of cells line
            grep "Estimated Number of Cells" "$summary_file" >> "$output_file"
 
            # Add a blank line for readability between samples
            echo "" >> "$output_file"

        else
            echo "Warning: Summary file not found for sample $sample_name"
        fi
    fi
done

echo "Alignment summary has been written to $output_file."
