#!/usr/bin/env python3
"""
Merge circRNA counts from CIRIquant/CIRI2 outputs into a count matrix

This script processes the output files from circRNA detection tools and
creates a unified count matrix suitable for differential expression analysis.

Usage:
    python merge_circrna_counts.py

Input: CIRIquant GTF files or CIRI2 output files
Output: Tab-separated count matrix
"""

import pandas as pd
import glob
import os
import sys
import re
from pathlib import Path

def parse_ciriquant_gtf(gtf_file):
    """Parse CIRIquant GTF output file"""
    circrna_data = {}

    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                if feature_type != 'circRNA':
                    continue

                attributes = fields[8]

                # Extract circRNA ID and BSJ count
                circ_id_match = re.search(r'circ_id "([^"]+)"', attributes)
                bsj_match = re.search(r'bsj "([^"]+)"', attributes)

                if circ_id_match and bsj_match:
                    circ_id = circ_id_match.group(1)
                    try:
                        bsj_count = int(bsj_match.group(1))
                        circrna_data[circ_id] = bsj_count
                    except ValueError:
                        continue

    except Exception as e:
        print(f"Error parsing {gtf_file}: {e}")

    return circrna_data

def parse_ciri2_output(ciri_file):
    """Parse CIRI2 output file"""
    circrna_data = {}

    try:
        df = pd.read_csv(ciri_file, sep='\t')

        # Check for expected columns
        if 'circRNA_ID' in df.columns and '#junction_reads' in df.columns:
            for _, row in df.iterrows():
                circ_id = row['circRNA_ID']
                junction_reads = row['#junction_reads']

                try:
                    circrna_data[circ_id] = int(junction_reads)
                except (ValueError, TypeError):
                    continue

    except Exception as e:
        print(f"Error parsing {ciri_file}: {e}")

    return circrna_data

def main():
    # Get output directory from snakemake or command line
    if 'snakemake' in globals():
        input_files = snakemake.input
        output_file = snakemake.output.count_matrix
        log_file = snakemake.log[0]
    else:
        # Fallback for direct execution
        input_files = glob.glob("../results/circrna/05_ciri_quant/*/*.gtf")
        if not input_files:
            input_files = glob.glob("../results/circrna/05_ciri_quant/*/*_ciri.out")
        output_file = "../results/circrna/05_ciri_quant/merged_count_matrix.tsv"
        log_file = "../results/circrna/logs/merge_counts.log"

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Initialize logging
    log_messages = []
    log_messages.append(f"Starting circRNA count matrix merging")
    log_messages.append(f"Found {len(input_files)} input files")

    count_dict = {}
    processed_samples = 0

    for file_path in input_files:
        # Extract sample name from file path
        sample_name = os.path.basename(os.path.dirname(file_path))
        log_messages.append(f"Processing sample: {sample_name} from {file_path}")

        # Determine file type and parse accordingly
        if file_path.endswith('.gtf'):
            sample_data = parse_ciriquant_gtf(file_path)
        elif file_path.endswith('_ciri.out'):
            sample_data = parse_ciri2_output(file_path)
        else:
            log_messages.append(f"Unknown file format: {file_path}")
            continue

        if sample_data:
            processed_samples += 1
            log_messages.append(f"Found {len(sample_data)} circRNAs in {sample_name}")

            # Add to master dictionary
            for circ_id, count in sample_data.items():
                if circ_id not in count_dict:
                    count_dict[circ_id] = {}
                count_dict[circ_id][sample_name] = count
        else:
            log_messages.append(f"No circRNA data found in {file_path}")

    # Create count matrix
    if count_dict:
        count_df = pd.DataFrame(count_dict).T.fillna(0).astype(int)
        count_df.index.name = 'circRNA_ID'

        # Sort columns (sample names) for consistency
        count_df = count_df.sort_index(axis=1)

        # Save count matrix
        count_df.to_csv(output_file, sep='\t')

        log_messages.append(f"Successfully created count matrix with {len(count_df)} circRNAs and {len(count_df.columns)} samples")
        log_messages.append(f"Count matrix saved to: {output_file}")

        # Print summary statistics
        total_counts = count_df.sum().sum()
        mean_circs_per_sample = count_df.sum(axis=0).mean()
        mean_counts_per_circ = count_df.sum(axis=1).mean()

        log_messages.append(f"Summary statistics:")
        log_messages.append(f"  - Total junction reads: {total_counts}")
        log_messages.append(f"  - Mean circRNAs per sample: {mean_circs_per_sample:.1f}")
        log_messages.append(f"  - Mean counts per circRNA: {mean_counts_per_circ:.1f}")

    else:
        log_messages.append("No circRNA data found in any input files!")
        log_messages.append("Creating empty count matrix")

        # Create empty matrix
        empty_df = pd.DataFrame()
        empty_df.index.name = 'circRNA_ID'
        empty_df.to_csv(output_file, sep='\t')

    # Write log file
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_messages))

    # Print to stdout for immediate feedback
    for msg in log_messages:
        print(msg)

if __name__ == "__main__":
    main()