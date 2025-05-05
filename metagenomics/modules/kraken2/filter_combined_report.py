#!/usr/bin/env python3
"""
filter_combined_report.py

Reads a combined Kraken2 report (tab-delimited) and outputs a reduced table containing only:
  - taxid
  - name (leading whitespace removed)
  - one column per sample (the `_lvl` counts, with `_lvl` suffix removed)

Extraneous header lines beginning with '#' are skipped entirely.
"""
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter combined Kraken2 report to taxid, name, and per-sample level counts"
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Path to the combined Kraken2 report TSV file'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Path for the filtered output TSV file'
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Read all lines from input
    try:
        lines = open(args.input, 'r').read().splitlines()
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

    header_cols = None
    sample_cols = []  # list of (index, sample_name)
    taxid_idx = None
    name_idx = None

    # Identify the header row (starts with '#perc') and parse column names
    for line in lines:
        if line.startswith('#perc'):
            cols = line.lstrip('#').split('\t')
            header_cols = [c.strip() for c in cols]
            for i, col in enumerate(header_cols):
                if col == 'taxid':
                    taxid_idx = i
                elif col == 'name':
                    name_idx = i
                elif col.endswith('_lvl') and col != 'tot_lvl':
                    sample_cols.append((i, col[:-4]))  # remove '_lvl'
            break

    if header_cols is None:
        sys.exit('Error: header line not found in input')
    if taxid_idx is None or name_idx is None or not sample_cols:
        sys.exit('Error: required columns (taxid, name, *_lvl) not found in header')

    # Write filtered output
    try:
        with open(args.output, 'w') as out:
            # Construct and write output header
            out_header = ['taxid', 'name'] + [name for _, name in sample_cols]
            out.write('\t'.join(out_header) + '\n')

            # Process each non-header, non-comment line
            for line in lines:
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) <= max(taxid_idx, name_idx):
                    continue
                # Strip whitespace around taxid and name
                taxid = parts[taxid_idx].strip()
                name = parts[name_idx].strip()

                # Extract level counts
                lvl_values = []
                for idx, _ in sample_cols:
                    val = parts[idx].strip() if idx < len(parts) else ''
                    lvl_values.append(val)

                out.write('\t'.join([taxid, name] + lvl_values) + '\n')
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

if __name__ == '__main__':
    main()
