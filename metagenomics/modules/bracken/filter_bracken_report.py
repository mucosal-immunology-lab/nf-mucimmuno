#!/usr/bin/env python3
import csv
import argparse
import os
import sys

def main():
    p = argparse.ArgumentParser(
        description="Slim down a combined Bracken report to just taxid, name, and counts per sample."
    )
    p.add_argument("-i", "--input", required=True,
                   help="Combined Bracken TSV input file")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV file")
    args = p.parse_args()

    # Open input and prepare output
    with open(args.input, newline='') as inf:
        reader = csv.reader(inf, delimiter='\t')
        header = next(reader)

        # Find indices
        try:
            idx_taxid = header.index("taxonomy_id")
            idx_name  = header.index("name")
        except ValueError as e:
            sys.exit(f"Error: missing column in header: {e}")

        # Identify all _num columns, keep their indices and derive sample names
        num_cols = []
        sample_names = []
        for i, col in enumerate(header):
            if col.endswith("_num"):
                num_cols.append(i)
                sample_names.append(col[:-4])  # strip "_num"

        if not num_cols:
            sys.exit("Error: no columns ending in '_num' found in header")

        # Build new header: taxid, name, sample1, sample2, ...
        new_header = ["taxonomy_id", "name"] + sample_names

        # Write output
        with open(args.output, "w", newline='') as outf:
            writer = csv.writer(outf, delimiter='\t')
            writer.writerow(new_header)

            for row in reader:
                outrow = [row[idx_taxid], row[idx_name]]
                for ci in num_cols:
                    outrow.append(row[ci])
                writer.writerow(outrow)

if __name__ == "__main__":
    main()
