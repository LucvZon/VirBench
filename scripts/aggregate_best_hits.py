#!/usr/bin/env python3

import sys
import csv
import argparse
from pathlib import Path

def main(output_file, pident_threshold, input_files):
    # The header list (3 new columns + 22 existing columns)
    headers = [
        "assembly_type", "assembler", "sample",
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", 
        "taxon_name", "species", "genus", "family", "order", "class", 
        "phylum", "kingdom", "superkingdom"
    ]

    with open(output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(headers) # Write the header line

        for file_path_str in input_files:
            p = Path(file_path_str)
            
            # Extract variables based on directory structure:
            # results/7_annotation/{assembly_type}/{sample}/post_processed/{assembler}_annotated_contigs.tsv
            assembler = p.name.replace("_annotated_contigs.tsv", "")
            sample = p.parent.parent.name
            assembly_type = p.parent.parent.parent.name
            
            try:
                with open(p, 'r') as in_f:
                    reader = csv.reader(in_f, delimiter='\t')
                    for row in reader:
                        # Skip empty rows or malformed rows
                        if not row or len(row) < 3:
                            continue
                        
                        try:
                            # pident is the 3rd column (index 2)
                            pident = float(row[2])
                        except ValueError:
                            continue # Skip row if pident is somehow not a number
                        
                        # Apply filter
                        if pident >= pident_threshold:
                            # Prepend the 3 new columns and write the row
                            writer.writerow([assembly_type, assembler, sample] + row)
                            
            except Exception as e:
                print(f"Warning: Could not process file {p}. Error: {e}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate filtered best hits from DIAMOND annotation TSVs.")
    parser.add_argument('-o', '--output', required=True, help="Output aggregated TSV file")
    parser.add_argument('-t', '--threshold', type=float, default=90.0, help="Minimum pident threshold (default: 90.0)")
    parser.add_argument('input_files', nargs='+', help="Input {assembler}_annotated_contigs.tsv files")
    
    args = parser.parse_args()
    main(args.output, args.threshold, args.input_files)
