#!/usr/bin/env python3
import csv
import argparse
import os
import sys

def load_taxonomy(db_path):
    """Load taxonomy nodes and names from a Kraken2/Bracken database."""
    nodes_path = os.path.join(db_path, "taxonomy", "nodes.dmp")
    names_path = os.path.join(db_path, "taxonomy", "names.dmp")

    if not (os.path.exists(nodes_path) and os.path.exists(names_path)):
        sys.exit(
            f"Error: could not find nodes.dmp or names.dmp in {os.path.join(db_path, 'taxonomy')}"
        )

    parents = {}
    ranks = {}
    with open(nodes_path) as nf:
        for line in nf:
            parts = [p.strip() for p in line.split("|")]
            taxid, parent, rank = parts[0], parts[1], parts[2]
            parents[taxid] = parent
            ranks[taxid] = rank

    names = {}
    with open(names_path) as nf:
        for line in nf:
            parts = [p.strip() for p in line.split("|")]
            taxid, name_txt, _, name_class = parts[:4]
            if name_class == "scientific name":
                names[taxid] = name_txt

    return parents, ranks, names

def lineage_for_taxid(taxid, parents, ranks, names, desired_ranks):
    lineage = {r: "" for r in desired_ranks}
    current = taxid
    visited = set()
    while current in parents and current not in visited:
        rank = ranks.get(current)
        rank_key = "kingdom" if rank == "superkingdom" else rank
        if rank_key in desired_ranks and not lineage[rank_key]:
            lineage[rank_key] = names.get(current, "")
            if all(lineage[r] for r in desired_ranks):
                break
        visited.add(current)
        parent = parents.get(current)
        if parent == current:
            break
        current = parent
    return lineage

def main():
    p = argparse.ArgumentParser(
        description="Slim down a combined Bracken report and append higher-level taxonomy ranks (genus through kingdom).",
    )
    p.add_argument("-i", "--input", required=True,
                   help="Combined Bracken TSV input file")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV file")
    p.add_argument("--taxonomy-db", required=True,
                   help="Path to Kraken2/Bracken database (with taxonomy/nodes.dmp and taxonomy/names.dmp)")
    args = p.parse_args()

    parents, ranks, names = load_taxonomy(args.taxonomy_db)
    desired_ranks = ["genus", "family", "order", "class", "phylum", "kingdom", "domain"]

    # Open input and read all rows
    with open(args.input, newline="") as inf:
        reader = csv.reader(inf, delimiter="\t")
        header = next(reader)
        rows = list(reader)

    # Find indices
    try:
        idx_taxid = header.index("taxonomy_id")
        idx_name = header.index("name")
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

    taxids = {row[idx_taxid] for row in rows}
    lineages = {
        tid: lineage_for_taxid(tid, parents, ranks, names, desired_ranks)
        for tid in taxids
    }

    # Build new header: taxid, name, ranks, sample1, sample2, ...
    new_header = ["taxonomy_id", "name"] + desired_ranks + sample_names

    # Write output
    with open(args.output, "w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(new_header)

        for row in rows:
            lineage = lineages.get(row[idx_taxid], {r: "" for r in desired_ranks})
            outrow = [row[idx_taxid], row[idx_name]] + [lineage[r] for r in desired_ranks]
            for ci in num_cols:
                outrow.append(row[ci])
            writer.writerow(outrow)

if __name__ == "__main__":
    main()
