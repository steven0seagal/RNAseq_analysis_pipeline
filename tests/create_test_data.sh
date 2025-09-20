#!/bin/bash

# ==============================================================================
# Test Data Creation Script
#
# This script creates minimal test FASTQ files for pipeline validation
# and testing purposes
#
# Usage: ./tests/create_test_data.sh
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TEST_DATA_DIR="$PROJECT_ROOT/data/test_data"

# Create test data directory
mkdir -p "$TEST_DATA_DIR"

echo "Creating test FASTQ files..."

# Sample names based on config/samples.tsv
SAMPLES=("Cancer_1" "Cancer_2" "Cancer_3" "Healthy_1" "Healthy_2" "Healthy_3")

# Create minimal FASTQ files for testing
for sample in "${SAMPLES[@]}"; do
    echo "Creating test data for sample: $sample"

    # Create R1 file
    cat > "$TEST_DATA_DIR/${sample}_R1.fastq" << EOF
@read1_${sample}
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_${sample}
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3_${sample}
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

    # Create R2 file
    cat > "$TEST_DATA_DIR/${sample}_R2.fastq" << EOF
@read1_${sample}
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_${sample}
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3_${sample}
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

    # Compress files
    gzip "$TEST_DATA_DIR/${sample}_R1.fastq"
    gzip "$TEST_DATA_DIR/${sample}_R2.fastq"
done

# Create a minimal reference genome for testing
echo "Creating minimal reference genome..."
cat > "$TEST_DATA_DIR/test_genome.fa" << EOF
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
>chr2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
EOF

# Create a minimal GTF annotation
echo "Creating minimal GTF annotation..."
cat > "$TEST_DATA_DIR/test_annotation.gtf" << EOF
chr1	test	gene	1	100	.	+	.	gene_id "GENE001"; gene_name "TESTGENE1"; gene_type "protein_coding";
chr1	test	transcript	1	100	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; gene_name "TESTGENE1"; gene_type "protein_coding";
chr1	test	exon	1	100	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; gene_name "TESTGENE1"; gene_type "protein_coding";
chr2	test	gene	1	100	.	+	.	gene_id "GENE002"; gene_name "TESTGENE2"; gene_type "protein_coding";
chr2	test	transcript	1	100	.	+	.	gene_id "GENE002"; transcript_id "TRANS002"; gene_name "TESTGENE2"; gene_type "protein_coding";
chr2	test	exon	1	100	.	+	.	gene_id "GENE002"; transcript_id "TRANS002"; gene_name "TESTGENE2"; gene_type "protein_coding";
EOF

# Create test configuration files
echo "Creating test configuration files..."

# Copy and modify samples.tsv for test data
sed 's/data\/raw_fastq/data\/test_data/g' "$PROJECT_ROOT/config/samples.tsv" > "$TEST_DATA_DIR/test_samples.tsv"

# Copy metadata.tsv
cp "$PROJECT_ROOT/config/metadata.tsv" "$TEST_DATA_DIR/test_metadata.tsv"

# Create test config.yaml
cat > "$TEST_DATA_DIR/test_config.yaml" << EOF
# Test configuration for minimal dataset

reference:
  genome_fasta: "data/test_data/test_genome.fa"
  annotation_gtf: "data/test_data/test_annotation.gtf"
  star_index: "data/test_data/star_index"

resources:
  threads: 2
  memory_gb: 4

params:
  fastp:
    min_len: 10
  star:
    sjdbOverhang: 49  # Shorter for test reads
  featureCounts:
    strand: 0  # Unstranded for simplicity

analysis:
  thresholds:
    padj: 0.05
    log2fc: 0.5  # Lower threshold for test data
    baseMean: 1
EOF

# Create a test runner script
cat > "$TEST_DATA_DIR/run_test.sh" << 'EOF'
#!/bin/bash

# Test runner script
echo "Running pipeline test with minimal data..."

# Set up test environment
export THREADS=2
export MEMORY=4
export SKIP_DOWNLOAD=true

# Create test results directory
mkdir -p ../results_test

# Run pipeline with test data
cd ..
./scripts/run_rnaseq_pipeline.sh \
    --threads 2 \
    --memory 4 \
    --skip-download \
    --no-interpretation

echo "Test run completed. Check results_test/ for output."
EOF

chmod +x "$TEST_DATA_DIR/run_test.sh"

# Create validation script
cat > "$TEST_DATA_DIR/validate_test.py" << 'EOF'
#!/usr/bin/env python3

"""
Validation script for test data pipeline run
"""

import os
import sys
import pandas as pd

def validate_test_results():
    """Validate that test pipeline produced expected outputs"""

    results_dir = "../results_test"

    # Check if results directory exists
    if not os.path.exists(results_dir):
        print("ERROR: Results directory not found")
        return False

    # Expected output files
    expected_files = [
        "01_fastqc_raw",
        "02_trimmed_fastq",
        "03_fastqc_trimmed",
        "05_featurecounts/raw_counts.tsv",
        "06_multiqc/multiqc_report.html"
    ]

    # Check for expected files
    missing_files = []
    for file_path in expected_files:
        full_path = os.path.join(results_dir, file_path)
        if not os.path.exists(full_path):
            missing_files.append(file_path)

    if missing_files:
        print(f"ERROR: Missing output files: {missing_files}")
        return False

    # Validate count matrix
    counts_file = os.path.join(results_dir, "05_featurecounts/raw_counts.tsv")
    try:
        counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
        print(f"Count matrix shape: {counts_df.shape}")

        if counts_df.shape[0] < 1:
            print("ERROR: No genes in count matrix")
            return False

        if counts_df.shape[1] < 6:  # 6 samples expected
            print("ERROR: Incorrect number of samples in count matrix")
            return False

    except Exception as e:
        print(f"ERROR: Could not read count matrix: {e}")
        return False

    print("✓ Test validation PASSED")
    return True

if __name__ == "__main__":
    success = validate_test_results()
    sys.exit(0 if success else 1)
EOF

chmod +x "$TEST_DATA_DIR/validate_test.py"

echo "✓ Test data creation completed!"
echo ""
echo "Created files:"
echo "- FASTQ files for 6 samples (3 Cancer, 3 Healthy)"
echo "- Minimal reference genome and annotation"
echo "- Test configuration files"
echo "- Test runner script: $TEST_DATA_DIR/run_test.sh"
echo "- Validation script: $TEST_DATA_DIR/validate_test.py"
echo ""
echo "To run tests:"
echo "1. cd $TEST_DATA_DIR"
echo "2. ./run_test.sh"
echo "3. ./validate_test.py"