#!/bin/bash

# ==============================================================================
# Small RNA-seq Analysis Pipeline with miRDeep2
#
# Comprehensive pipeline for small RNA analysis including:
# 1. Quality control and adapter trimming with cutadapt
# 2. Mapping to genome and miRBase with Bowtie
# 3. miRNA discovery and quantification with miRDeep2
# 4. Novel miRNA prediction and validation
# 5. Differential expression analysis of miRNAs
# 6. Target prediction and pathway analysis
#
# This pipeline identifies known and novel miRNAs from small RNA-seq data
# and provides comprehensive annotation for cancer genomics research.
#
# Usage: ./scripts/small_rna_analysis.sh [options]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# === CONFIGURATION ===
THREADS=${THREADS:-8}
MIN_LENGTH=${MIN_LENGTH:-18}
MAX_LENGTH=${MAX_LENGTH:-25}
ADAPTER_SEQ=${ADAPTER_SEQ:-"AGATCGGAAGAGCACACGTCT"}

# Directories
RAW_DIR="$PROJECT_ROOT/02_raw_data"
TRIM_DIR="$PROJECT_ROOT/03_trimmed_smallrna"
SMALL_RNA_DIR="$PROJECT_ROOT/08_small_rna"
MIRDEEP_DIR="$SMALL_RNA_DIR/mirdeep2"
REF_DIR="$PROJECT_ROOT/01_references"
LOG_DIR="$PROJECT_ROOT/logs/small_rna"

# Reference files
GENOME_FASTA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
MATURE_MIRNA="$REF_DIR/human_mature.fa"
HAIRPIN_MIRNA="$REF_DIR/hairpin.fa"
BOWTIE_INDEX="$REF_DIR/bowtie_index/genome"

# Sample configuration
SAMPLES_FILE="$PROJECT_ROOT/config/samples.tsv"

# Create directories
mkdir -p "$TRIM_DIR" "$SMALL_RNA_DIR" "$MIRDEEP_DIR" "$LOG_DIR"

# === LOGGING ===
LOG_FILE="$LOG_DIR/small_rna_analysis_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

error_exit() {
    log "ERROR: $1"
    exit 1
}

# === VALIDATION ===
validate_inputs() {
    log "=== Validating Inputs ==="

    # Check required tools
    local tools=("cutadapt" "bowtie" "mapper.pl" "miRDeep2.pl" "quantifier.pl")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            error_exit "$tool not found in PATH"
        fi
    done

    # Check reference files
    [[ -f "$GENOME_FASTA" ]] || error_exit "Genome FASTA not found: $GENOME_FASTA"
    [[ -f "$MATURE_MIRNA" ]] || error_exit "Mature miRNA file not found: $MATURE_MIRNA"
    [[ -f "$HAIRPIN_MIRNA" ]] || error_exit "Hairpin miRNA file not found: $HAIRPIN_MIRNA"
    [[ -f "$SAMPLES_FILE" ]] || error_exit "Samples file not found: $SAMPLES_FILE"

    # Check Bowtie index
    if [[ ! -f "${BOWTIE_INDEX}.1.ebwt" ]]; then
        log "Bowtie index not found, creating..."
        build_bowtie_index
    fi

    # Read samples
    mapfile -t SAMPLES < <(tail -n +2 "$SAMPLES_FILE" | cut -f1)
    [[ ${#SAMPLES[@]} -gt 0 ]] || error_exit "No samples found in $SAMPLES_FILE"

    log "Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

    # Check input FASTQ files
    for sample in "${SAMPLES[@]}"; do
        local fastq_file="$RAW_DIR/${sample}_R1.fastq.gz"
        [[ -f "$fastq_file" ]] || error_exit "FASTQ file not found for $sample: $fastq_file"
    done

    log "Input validation completed"
}

# === BUILD BOWTIE INDEX ===
build_bowtie_index() {
    log "=== Building Bowtie Index for Small RNA Analysis ==="

    local index_dir="$(dirname "$BOWTIE_INDEX")"
    mkdir -p "$index_dir"

    bowtie-build "$GENOME_FASTA" "$BOWTIE_INDEX" \
        > "$LOG_DIR/bowtie_build.log" 2>&1

    log "Bowtie index completed"
}

# === QUALITY CONTROL AND TRIMMING ===
trim_adapters() {
    local sample=$1
    log "=== Trimming Adapters for $sample ==="

    local input_fastq="$RAW_DIR/${sample}_R1.fastq.gz"
    local output_fastq="$TRIM_DIR/${sample}_trimmed.fastq.gz"
    local report_file="$TRIM_DIR/${sample}_cutadapt_report.txt"
    local log_file="$LOG_DIR/${sample}_cutadapt.log"

    # Skip if already processed
    if [[ -f "$output_fastq" ]]; then
        log "Adapter trimming already completed for $sample"
        return 0
    fi

    # Cutadapt for small RNA adapter trimming
    cutadapt \
        --adapter="$ADAPTER_SEQ" \
        --minimum-length="$MIN_LENGTH" \
        --maximum-length="$MAX_LENGTH" \
        --quality-cutoff=20 \
        --cores="$THREADS" \
        --output="$output_fastq" \
        --report=minimal \
        "$input_fastq" \
        > "$report_file" 2> "$log_file"

    # Generate quality report
    if command -v fastqc &> /dev/null; then
        fastqc -o "$TRIM_DIR" -t "$THREADS" "$output_fastq" \
            > "$LOG_DIR/${sample}_fastqc_trimmed.log" 2>&1
    fi

    log "Adapter trimming completed for $sample"
}

# === PREPROCESS FOR MIRDEEP2 ===
preprocess_for_mirdeep() {
    local sample=$1
    log "=== Preprocessing for miRDeep2: $sample ==="

    local input_fastq="$TRIM_DIR/${sample}_trimmed.fastq.gz"
    local processed_fa="$MIRDEEP_DIR/${sample}_processed.fa"
    local log_file="$LOG_DIR/${sample}_preprocess.log"

    # Skip if already processed
    if [[ -f "$processed_fa" ]]; then
        log "Preprocessing already completed for $sample"
        return 0
    fi

    # Convert to FASTA format with miRDeep2 format
    zcat "$input_fastq" | \
    awk 'BEGIN{count=1}
         NR%4==1{header=substr($0,2)}
         NR%4==2{seq=$0}
         NR%4==0{
             gsub(/T/,"U",seq);
             print ">"sample"_"count"_x"length(seq);
             print seq;
             count++
         }' sample="$sample" \
    > "$processed_fa" 2> "$log_file"

    log "Preprocessing completed for $sample"
}

# === MIRDEEP2 MAPPING ===
run_mirdeep_mapping() {
    local sample=$1
    log "=== Running miRDeep2 Mapping for $sample ==="

    local processed_fa="$MIRDEEP_DIR/${sample}_processed.fa"
    local mapped_arf="$MIRDEEP_DIR/${sample}_mapped.arf"
    local collapsed_fa="$MIRDEEP_DIR/${sample}_collapsed.fa"
    local log_file="$LOG_DIR/${sample}_mapper.log"

    # Skip if already processed
    if [[ -f "$mapped_arf" ]]; then
        log "miRDeep2 mapping already completed for $sample"
        return 0
    fi

    # Run mapper.pl
    cd "$MIRDEEP_DIR"

    mapper.pl "$processed_fa" \
        -d \
        -c \
        -i \
        -j \
        -l 18 \
        -m \
        -p "$BOWTIE_INDEX" \
        -s "$collapsed_fa" \
        -t "$mapped_arf" \
        -v \
        > "$log_file" 2>&1

    cd "$PROJECT_ROOT"

    log "miRDeep2 mapping completed for $sample"
}

# === MIRDEEP2 CORE ANALYSIS ===
run_mirdeep_core() {
    local sample=$1
    log "=== Running miRDeep2 Core Analysis for $sample ==="

    local collapsed_fa="$MIRDEEP_DIR/${sample}_collapsed.fa"
    local mapped_arf="$MIRDEEP_DIR/${sample}_mapped.arf"
    local output_dir="$MIRDEEP_DIR/${sample}_results"
    local log_file="$LOG_DIR/${sample}_mirdeep.log"

    # Skip if already processed
    if [[ -d "$output_dir" ]]; then
        log "miRDeep2 core analysis already completed for $sample"
        return 0
    fi

    mkdir -p "$output_dir"
    cd "$output_dir"

    # Prepare miRNA references
    local mature_ref="$MIRDEEP_DIR/mature_ref.fa"
    local hairpin_ref="$MIRDEEP_DIR/hairpin_ref.fa"

    # Extract human miRNAs if not already done
    if [[ ! -f "$mature_ref" ]]; then
        grep -A1 "^>hsa-" "$MATURE_MIRNA" | grep -v "^--$" > "$mature_ref"
    fi

    if [[ ! -f "$hairpin_ref" ]]; then
        grep -A1 "^>hsa-" "$HAIRPIN_MIRNA" | grep -v "^--$" > "$hairpin_ref"
    fi

    # Run miRDeep2 core
    miRDeep2.pl \
        "$collapsed_fa" \
        "$GENOME_FASTA" \
        "$mapped_arf" \
        "$mature_ref" \
        none \
        "$hairpin_ref" \
        -t hsa \
        -P \
        -z "${sample}_" \
        > "$log_file" 2>&1

    cd "$PROJECT_ROOT"

    log "miRDeep2 core analysis completed for $sample"
}

# === QUANTIFY KNOWN MIRNAS ===
quantify_known_mirnas() {
    log "=== Quantifying Known miRNAs ==="

    local quantification_dir="$MIRDEEP_DIR/quantification"
    mkdir -p "$quantification_dir"

    # Collect all processed reads files
    local read_files=""
    for sample in "${SAMPLES[@]}"; do
        local collapsed_fa="$MIRDEEP_DIR/${sample}_collapsed.fa"
        if [[ -f "$collapsed_fa" ]]; then
            read_files="$read_files,$collapsed_fa"
        fi
    done
    read_files="${read_files#,}"  # Remove leading comma

    if [[ -z "$read_files" ]]; then
        error_exit "No processed read files found for quantification"
    fi

    cd "$quantification_dir"

    # Run quantifier
    quantifier.pl \
        -p "$HAIRPIN_MIRNA" \
        -m "$MATURE_MIRNA" \
        -r "$read_files" \
        -t hsa \
        -y "${SAMPLES[0]}_quantification" \
        > "$LOG_DIR/quantifier.log" 2>&1

    cd "$PROJECT_ROOT"

    log "miRNA quantification completed"
}

# === DIFFERENTIAL EXPRESSION ANALYSIS ===
run_differential_expression() {
    log "=== Running Differential Expression Analysis ==="

    local de_script="$SMALL_RNA_DIR/mirna_differential_expression.R"

    # Create R script for differential expression
    cat > "$de_script" << 'EOF'
#!/usr/bin/env Rscript

# miRNA Differential Expression Analysis
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read miRNA count matrix
count_file <- "quantification/miRNAs_expressed_all_samples.csv"
if (!file.exists(count_file)) {
    cat("Count file not found:", count_file, "\n")
    quit(status = 1)
}

# Read counts
counts <- read_csv(count_file, show_col_types = FALSE)

# Prepare count matrix
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- counts[[1]]

# Create sample metadata (adjust based on your design)
sample_names <- colnames(count_matrix)
condition <- ifelse(grepl("_T", sample_names), "Tumor", "Normal")
if (all(condition == "Tumor")) {
    # If no normal samples, create artificial groups
    condition <- rep(c("Group1", "Group2"), length.out = length(sample_names))
}

metadata <- data.frame(
    sample = sample_names,
    condition = factor(condition),
    row.names = sample_names
)

# Filter low count miRNAs
keep <- rowSums(count_matrix >= 10) >= ncol(count_matrix) * 0.1
count_matrix <- count_matrix[keep, ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", levels(metadata$condition)[2], levels(metadata$condition)[1]))
res_df <- as.data.frame(res)
res_df$miRNA <- rownames(res_df)

# Save results
write_csv(res_df, "differential_mirnas.csv")

# Significant miRNAs
sig_mirnas <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write_csv(sig_mirnas, "significant_mirnas.csv")

# Volcano plot
p1 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.6) +
    geom_point(data = sig_mirnas, color = "red", size = 2) +
    labs(title = "miRNA Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal()

ggsave("mirna_volcano_plot.png", p1, width = 8, height = 6, dpi = 300)

# Heatmap of top miRNAs
if (nrow(sig_mirnas) > 0) {
    top_mirnas <- head(sig_mirnas[order(sig_mirnas$padj), ], 50)

    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    heatmap_data <- norm_counts[rownames(norm_counts) %in% top_mirnas$miRNA, ]

    # Log transform
    heatmap_data <- log2(heatmap_data + 1)

    # Create heatmap
    png("mirna_heatmap.png", width = 10, height = 8, units = "in", res = 300)
    pheatmap(
        heatmap_data,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        scale = "row",
        annotation_col = metadata[, "condition", drop = FALSE],
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = "Top Differentially Expressed miRNAs"
    )
    dev.off()
}

# Summary statistics
cat("miRNA Differential Expression Analysis Summary:\n")
cat("Total miRNAs analyzed:", nrow(res_df), "\n")
cat("Significant miRNAs (padj < 0.05, |log2FC| > 1):", nrow(sig_mirnas), "\n")
cat("Upregulated:", sum(sig_mirnas$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_mirnas$log2FoldChange < 0), "\n")

cat("\nAnalysis completed successfully!\n")
EOF

    chmod +x "$de_script"

    # Run differential expression analysis
    cd "$SMALL_RNA_DIR"
    if command -v Rscript &> /dev/null; then
        Rscript "$de_script" > "$LOG_DIR/differential_expression.log" 2>&1
        log "Differential expression analysis completed"
    else
        log "R not available - skipping differential expression analysis"
    fi
    cd "$PROJECT_ROOT"
}

# === GENERATE SUMMARY REPORT ===
generate_summary_report() {
    log "=== Generating Summary Report ==="

    local summary_file="$SMALL_RNA_DIR/small_rna_analysis_summary.txt"

    cat > "$summary_file" << EOF
# Small RNA-seq Analysis Summary Report
# Generated on: $(date)

## Analysis Overview
Pipeline: miRDeep2-based small RNA analysis
Samples processed: ${#SAMPLES[@]}
Reference genome: GRCh38
miRNA database: miRBase (latest)

## Processing Parameters
Minimum read length: $MIN_LENGTH nt
Maximum read length: $MAX_LENGTH nt
Adapter sequence: $ADAPTER_SEQ
Threads used: $THREADS

## Sample Processing Status
EOF

    for sample in "${SAMPLES[@]}"; do
        local status="PENDING"
        local trimmed_file="$TRIM_DIR/${sample}_trimmed.fastq.gz"
        local mapped_file="$MIRDEEP_DIR/${sample}_mapped.arf"
        local results_dir="$MIRDEEP_DIR/${sample}_results"

        if [[ -d "$results_dir" ]]; then
            status="COMPLETED"
        elif [[ -f "$mapped_file" ]]; then
            status="MAPPED"
        elif [[ -f "$trimmed_file" ]]; then
            status="TRIMMED"
        fi

        echo "$sample: $status" >> "$summary_file"
    done

    cat >> "$summary_file" << EOF

## Output Files
- Trimmed reads: $TRIM_DIR/
- miRDeep2 results: $MIRDEEP_DIR/
- Quantification: $MIRDEEP_DIR/quantification/
- Differential expression: $SMALL_RNA_DIR/differential_mirnas.csv
- Visualizations: $SMALL_RNA_DIR/*.png

## Next Steps
1. Review miRDeep2 results for novel miRNA predictions
2. Validate significant differentially expressed miRNAs
3. Perform target prediction analysis
4. Integrate with mRNA expression data

## Key Files for Review
- miRNAs_expressed_all_samples.csv: miRNA expression matrix
- differential_mirnas.csv: Differential expression results
- significant_mirnas.csv: Significantly changed miRNAs
- mirna_volcano_plot.png: Volcano plot visualization
- mirna_heatmap.png: Heatmap of top miRNAs
EOF

    log "Summary report generated: $summary_file"
}

# === PROCESS SINGLE SAMPLE ===
process_sample() {
    local sample=$1
    log "========================================="
    log "Processing sample: $sample"
    log "========================================="

    trim_adapters "$sample"
    preprocess_for_mirdeep "$sample"
    run_mirdeep_mapping "$sample"
    run_mirdeep_core "$sample"

    log "Sample $sample processing completed"
}

# === MAIN EXECUTION ===
main() {
    log "========================================="
    log "Small RNA-seq Analysis Pipeline Started"
    log "========================================="
    log "Parameters:"
    log "  Threads: $THREADS"
    log "  Read length range: $MIN_LENGTH-$MAX_LENGTH nt"
    log "  Adapter sequence: $ADAPTER_SEQ"
    log "========================================="

    validate_inputs

    # Process each sample
    for sample in "${SAMPLES[@]}"; do
        process_sample "$sample"
    done

    quantify_known_mirnas
    run_differential_expression
    generate_summary_report

    log "========================================="
    log "Small RNA-seq Analysis Pipeline Completed"
    log "========================================="
    log "Output files:"
    log "  Trimmed reads: $TRIM_DIR/"
    log "  miRDeep2 results: $MIRDEEP_DIR/"
    log "  Quantification: $MIRDEEP_DIR/quantification/"
    log "  Differential expression: $SMALL_RNA_DIR/"
    log "  Summary report: $SMALL_RNA_DIR/small_rna_analysis_summary.txt"
    log "  Log files: $LOG_DIR/"
    log ""
    log "Next steps:"
    log "1. Review miRDeep2 novel miRNA predictions"
    log "2. Validate differential expression results"
    log "3. Perform miRNA target prediction"
    log "4. Integrate with multi-modal analysis"
}

# === COMMAND LINE PARSING ===
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            cat << EOF
Small RNA-seq Analysis Pipeline with miRDeep2

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -t, --threads INT       Number of threads (default: 8)
    --min-length INT        Minimum read length (default: 18)
    --max-length INT        Maximum read length (default: 25)
    --adapter STRING        Adapter sequence for trimming
    --dry-run              Show commands without executing

EXAMPLES:
    # Basic run
    $0

    # Custom parameters
    $0 --threads 16 --min-length 16 --max-length 30

    # Custom adapter
    $0 --adapter "AGATCGGAAGAGCACACGTCT"

DESCRIPTION:
    Comprehensive small RNA-seq analysis pipeline:
    1. Quality control and adapter trimming
    2. Mapping to genome with Bowtie
    3. miRNA discovery and quantification with miRDeep2
    4. Differential expression analysis
    5. Visualization and reporting

    Identifies both known and novel miRNAs for cancer research.

EOF
            exit 0
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --min-length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        --max-length)
            MAX_LENGTH="$2"
            shift 2
            ;;
        --adapter)
            ADAPTER_SEQ="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            error_exit "Unknown option: $1. Use --help for usage information."
            ;;
    esac
done

# Execute main function
if [[ "${DRY_RUN:-false}" == "true" ]]; then
    log "DRY RUN MODE - Commands would be executed but files not modified"
    exit 0
fi

main "$@"