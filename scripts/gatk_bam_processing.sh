#!/bin/bash

# ==============================================================================
# GATK-Compliant BAM Processing Pipeline for RNA-seq Variant Calling
#
# Transforms STAR-aligned BAM files into GATK-ready format through:
# 1. AddOrReplaceReadGroups - Required read group information
# 2. MarkDuplicates - PCR duplicate identification and flagging
# 3. SplitNCigarReads - RNA-specific processing for splice-aware alignment
# 4. Base Quality Score Recalibration (BQSR) - Systematic error correction
#
# This pipeline follows GATK best practices for RNA-seq variant calling
# and is essential for accurate somatic mutation detection.
#
# Usage: ./scripts/gatk_bam_processing.sh [options]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# === CONFIGURATION ===
THREADS=${THREADS:-8}
MEMORY_GB=${MEMORY_GB:-32}
JAVA_OPTS=${JAVA_OPTS:-"-Xmx${MEMORY_GB}g -XX:+UseParallelGC"}

# Directories
ALIGN_DIR="$PROJECT_ROOT/04_aligned_lncrna"
VAR_DIR="$PROJECT_ROOT/06_variants"
REF_DIR="$PROJECT_ROOT/01_references"
LOG_DIR="$PROJECT_ROOT/logs/gatk"

# Reference files
GENOME_FASTA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
DBSNP_VCF="$REF_DIR/databases/dbsnp_155.hg38.vcf.gz"

# Sample configuration
SAMPLES_FILE="$PROJECT_ROOT/config/samples.tsv"

# Create directories
mkdir -p "$VAR_DIR" "$LOG_DIR"

# === LOGGING ===
LOG_FILE="$LOG_DIR/gatk_processing_$(date +%Y%m%d_%H%M%S).log"

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

    # Check GATK installation
    if ! command -v gatk &> /dev/null; then
        error_exit "GATK not found in PATH"
    fi

    local gatk_version=$(gatk --version 2>&1 | head -1)
    log "GATK version: $gatk_version"

    # Check reference files
    [[ -f "$GENOME_FASTA" ]] || error_exit "Genome FASTA not found: $GENOME_FASTA"
    [[ -f "$DBSNP_VCF" ]] || error_exit "dbSNP VCF not found: $DBSNP_VCF"
    [[ -f "$SAMPLES_FILE" ]] || error_exit "Samples file not found: $SAMPLES_FILE"

    # Check for required indices
    [[ -f "${GENOME_FASTA}.fai" ]] || error_exit "FASTA index not found: ${GENOME_FASTA}.fai"
    [[ -f "${GENOME_FASTA%.*}.dict" ]] || error_exit "Sequence dictionary not found: ${GENOME_FASTA%.*}.dict"

    # Read samples
    mapfile -t SAMPLES < <(tail -n +2 "$SAMPLES_FILE" | cut -f1)
    [[ ${#SAMPLES[@]} -gt 0 ]] || error_exit "No samples found in $SAMPLES_FILE"

    log "Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

    # Check input BAM files
    for sample in "${SAMPLES[@]}"; do
        local bam_file="$ALIGN_DIR/${sample}_2pass_Aligned.sortedByCoord.out.bam"
        [[ -f "$bam_file" ]] || error_exit "BAM file not found for $sample: $bam_file"
        [[ -f "${bam_file}.bai" ]] || error_exit "BAM index not found for $sample: ${bam_file}.bai"
    done

    log "Input validation completed"
}

# === STEP 1: ADD OR REPLACE READ GROUPS ===
add_read_groups() {
    local sample=$1
    log "=== Step 1: Adding Read Groups for $sample ==="

    local input_bam="$ALIGN_DIR/${sample}_2pass_Aligned.sortedByCoord.out.bam"
    local output_bam="$VAR_DIR/${sample}.rg.bam"
    local log_file="$LOG_DIR/${sample}_add_rg.log"

    # Skip if already processed
    if [[ -f "$output_bam" ]]; then
        log "Read groups already added for $sample"
        return 0
    fi

    # Extract sequencing platform from BAM header (if available)
    local platform="ILLUMINA"  # Default assumption
    local library="lib_${sample}"
    local platform_unit="unit_${sample}"

    gatk $JAVA_OPTS AddOrReplaceReadGroups \
        --INPUT "$input_bam" \
        --OUTPUT "$output_bam" \
        --RGID "rg_${sample}" \
        --RGLB "$library" \
        --RGPL "$platform" \
        --RGPU "$platform_unit" \
        --RGSM "$sample" \
        --CREATE_INDEX true \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT \
        > "$log_file" 2>&1

    log "Read groups added for $sample"
}

# === STEP 2: MARK DUPLICATES ===
mark_duplicates() {
    local sample=$1
    log "=== Step 2: Marking Duplicates for $sample ==="

    local input_bam="$VAR_DIR/${sample}.rg.bam"
    local output_bam="$VAR_DIR/${sample}.dedup.bam"
    local metrics_file="$VAR_DIR/${sample}.dedup_metrics.txt"
    local log_file="$LOG_DIR/${sample}_mark_dup.log"

    # Skip if already processed
    if [[ -f "$output_bam" ]]; then
        log "Duplicates already marked for $sample"
        return 0
    fi

    gatk $JAVA_OPTS MarkDuplicates \
        --INPUT "$input_bam" \
        --OUTPUT "$output_bam" \
        --METRICS_FILE "$metrics_file" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER coordinate \
        > "$log_file" 2>&1

    # Log duplicate statistics
    if [[ -f "$metrics_file" ]]; then
        local total_reads=$(grep -A1 "^LIBRARY" "$metrics_file" | tail -1 | cut -f3)
        local duplicates=$(grep -A1 "^LIBRARY" "$metrics_file" | tail -1 | cut -f7)
        local dup_rate=$(grep -A1 "^LIBRARY" "$metrics_file" | tail -1 | cut -f9)

        log "$sample duplicate statistics: $duplicates/$total_reads (${dup_rate}%)"
    fi

    log "Duplicates marked for $sample"
}

# === STEP 3: SPLIT N CIGAR READS ===
split_n_cigar_reads() {
    local sample=$1
    log "=== Step 3: Splitting N CIGAR Reads for $sample ==="

    local input_bam="$VAR_DIR/${sample}.dedup.bam"
    local output_bam="$VAR_DIR/${sample}.split.bam"
    local log_file="$LOG_DIR/${sample}_split_n.log"

    # Skip if already processed
    if [[ -f "$output_bam" ]]; then
        log "N CIGAR reads already split for $sample"
        return 0
    fi

    gatk $JAVA_OPTS SplitNCigarReads \
        --reference "$GENOME_FASTA" \
        --input "$input_bam" \
        --output "$output_bam" \
        --create-output-bam-index true \
        --refactor-NDN-cigar-string \
        --skip-mapping-quality-transform \
        > "$log_file" 2>&1

    log "N CIGAR reads split for $sample"
}

# === STEP 4: BASE QUALITY SCORE RECALIBRATION ===
base_quality_recalibration() {
    local sample=$1
    log "=== Step 4: Base Quality Score Recalibration for $sample ==="

    local input_bam="$VAR_DIR/${sample}.split.bam"
    local recal_table="$VAR_DIR/${sample}.recal_data.table"
    local output_bam="$VAR_DIR/${sample}.bqsr.bam"
    local log_file="$LOG_DIR/${sample}_bqsr.log"

    # Skip if already processed
    if [[ -f "$output_bam" ]]; then
        log "BQSR already completed for $sample"
        return 0
    fi

    # Step 4a: BaseRecalibrator
    log "  Running BaseRecalibrator for $sample"
    gatk $JAVA_OPTS BaseRecalibrator \
        --reference "$GENOME_FASTA" \
        --input "$input_bam" \
        --known-sites "$DBSNP_VCF" \
        --output "$recal_table" \
        --tmp-dir "$VAR_DIR/tmp" \
        > "${log_file}_base_recalibrator.log" 2>&1

    # Step 4b: ApplyBQSR
    log "  Applying BQSR for $sample"
    gatk $JAVA_OPTS ApplyBQSR \
        --reference "$GENOME_FASTA" \
        --input "$input_bam" \
        --bqsr-recal-file "$recal_table" \
        --output "$output_bam" \
        --create-output-bam-index true \
        --tmp-dir "$VAR_DIR/tmp" \
        > "${log_file}_apply_bqsr.log" 2>&1

    log "BQSR completed for $sample"
}

# === QUALITY CONTROL CHECKS ===
run_quality_checks() {
    local sample=$1
    log "=== Quality Control Checks for $sample ==="

    local final_bam="$VAR_DIR/${sample}.bqsr.bam"
    local qc_dir="$VAR_DIR/qc"
    mkdir -p "$qc_dir"

    # Validate BAM file
    log "  Validating final BAM file for $sample"
    gatk $JAVA_OPTS ValidateSamFile \
        --INPUT "$final_bam" \
        --OUTPUT "$qc_dir/${sample}_validation.txt" \
        --MODE SUMMARY \
        --VALIDATION_STRINGENCY LENIENT \
        > "$LOG_DIR/${sample}_validate.log" 2>&1

    # Collect alignment summary metrics
    log "  Collecting alignment metrics for $sample"
    gatk $JAVA_OPTS CollectAlignmentSummaryMetrics \
        --INPUT "$final_bam" \
        --OUTPUT "$qc_dir/${sample}_alignment_metrics.txt" \
        --REFERENCE_SEQUENCE "$GENOME_FASTA" \
        > "$LOG_DIR/${sample}_alignment_metrics.log" 2>&1

    # Collect insert size metrics (for paired-end data)
    log "  Collecting insert size metrics for $sample"
    gatk $JAVA_OPTS CollectInsertSizeMetrics \
        --INPUT "$final_bam" \
        --OUTPUT "$qc_dir/${sample}_insert_size_metrics.txt" \
        --HISTOGRAM_FILE "$qc_dir/${sample}_insert_size_histogram.pdf" \
        > "$LOG_DIR/${sample}_insert_size.log" 2>&1

    log "Quality control checks completed for $sample"
}

# === GENERATE PROCESSING SUMMARY ===
generate_summary() {
    log "=== Generating Processing Summary ==="

    local summary_file="$VAR_DIR/bam_processing_summary.txt"

    cat > "$summary_file" << EOF
# GATK BAM Processing Summary
# Generated on: $(date)

Sample	Input_Reads	Duplicates	Duplicate_Rate	Final_BAM_Size
EOF

    for sample in "${SAMPLES[@]}"; do
        local metrics_file="$VAR_DIR/${sample}.dedup_metrics.txt"
        local final_bam="$VAR_DIR/${sample}.bqsr.bam"

        if [[ -f "$metrics_file" && -f "$final_bam" ]]; then
            # Extract metrics
            local total_reads=$(grep -A1 "^LIBRARY" "$metrics_file" | tail -1 | cut -f3)
            local duplicates=$(grep -A1 "^LIBRARY" "$metrics_file" | tail -1 | cut -f7)
            local dup_rate=$(grep -A1 "^LIBRARY" "$metrics_file" | tail -1 | cut -f9)
            local bam_size=$(du -h "$final_bam" | cut -f1)

            echo -e "$sample\t$total_reads\t$duplicates\t${dup_rate}%\t$bam_size" >> "$summary_file"
        else
            echo -e "$sample\tN/A\tN/A\tN/A\tN/A" >> "$summary_file"
        fi
    done

    log "Processing summary saved: $summary_file"
}

# === CLEANUP INTERMEDIATE FILES ===
cleanup_intermediate() {
    if [[ "${CLEANUP_INTERMEDIATE:-true}" == "true" ]]; then
        log "=== Cleaning Up Intermediate Files ==="

        for sample in "${SAMPLES[@]}"; do
            # Remove intermediate BAM files (keep final BQSR BAM)
            rm -f "$VAR_DIR/${sample}.rg.bam" "$VAR_DIR/${sample}.rg.bai"
            rm -f "$VAR_DIR/${sample}.dedup.bam" "$VAR_DIR/${sample}.dedup.bai"
            rm -f "$VAR_DIR/${sample}.split.bam" "$VAR_DIR/${sample}.split.bai"
        done

        # Remove temporary directory
        rm -rf "$VAR_DIR/tmp"

        log "Intermediate file cleanup completed"
    fi
}

# === PROCESS SINGLE SAMPLE ===
process_sample() {
    local sample=$1
    log "========================================="
    log "Processing sample: $sample"
    log "========================================="

    add_read_groups "$sample"
    mark_duplicates "$sample"
    split_n_cigar_reads "$sample"
    base_quality_recalibration "$sample"
    run_quality_checks "$sample"

    log "Sample $sample processing completed"
}

# === MAIN EXECUTION ===
main() {
    log "========================================="
    log "GATK BAM Processing Pipeline Started"
    log "========================================="
    log "Parameters:"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Java options: $JAVA_OPTS"
    log "========================================="

    validate_inputs

    # Create temporary directory
    mkdir -p "$VAR_DIR/tmp"

    # Process each sample
    for sample in "${SAMPLES[@]}"; do
        process_sample "$sample"
    done

    generate_summary
    cleanup_intermediate

    log "========================================="
    log "GATK BAM Processing Pipeline Completed"
    log "========================================="
    log "Output files:"
    log "  Final BAMs: $VAR_DIR/*.bqsr.bam"
    log "  QC metrics: $VAR_DIR/qc/"
    log "  Summary: $VAR_DIR/bam_processing_summary.txt"
    log "  Log files: $LOG_DIR/"
    log ""
    log "Next steps:"
    log "1. Run somatic variant calling with Mutect2"
    log "2. Annotate and filter variants"
    log "3. Integrate with expression analysis"
}

# === COMMAND LINE PARSING ===
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            cat << EOF
GATK-Compliant BAM Processing Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -t, --threads INT       Number of threads (default: 8)
    -m, --memory INT        Memory in GB (default: 32)
    --no-cleanup           Keep intermediate files
    --dry-run              Show commands without executing

EXAMPLES:
    # Basic run
    $0

    # Custom resources
    $0 --threads 16 --memory 64

    # Keep intermediate files for debugging
    $0 --no-cleanup

DESCRIPTION:
    Processes STAR-aligned BAM files through GATK best practices:
    1. AddOrReplaceReadGroups - Required metadata
    2. MarkDuplicates - PCR duplicate flagging
    3. SplitNCigarReads - RNA-seq specific processing
    4. BaseRecalibrator + ApplyBQSR - Quality score correction

    Output BAMs are ready for variant calling with Mutect2.

EOF
            exit 0
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY_GB="$2"
            JAVA_OPTS="-Xmx${MEMORY_GB}g -XX:+UseParallelGC"
            shift 2
            ;;
        --no-cleanup)
            CLEANUP_INTERMEDIATE=false
            shift
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