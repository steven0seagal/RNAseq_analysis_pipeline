#!/bin/bash

# ==============================================================================
# STAR 2-Pass Alignment Pipeline for Enhanced Splice Junction Detection
#
# Implements GATK best practices for RNA-seq alignment:
# 1. First pass: Discover splice junctions from all samples
# 2. Build augmented genome index with novel junctions
# 3. Second pass: Re-align with improved junction sensitivity
#
# This approach maximizes alignment accuracy around splice sites and
# is critical for accurate variant calling in RNA-seq data.
#
# Usage: ./scripts/run_star_2pass.sh [options]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# === CONFIGURATION ===
THREADS=${THREADS:-8}
MEMORY_GB=${MEMORY_GB:-32}
READ_LENGTH=${READ_LENGTH:-100}

# Directories
TRIM_DIR="$PROJECT_ROOT/03_trimmed"
ALIGN_DIR="$PROJECT_ROOT/04_aligned_lncrna"
REF_DIR="$PROJECT_ROOT/01_references"
LOG_DIR="$PROJECT_ROOT/logs/star_2pass"

# Reference files
GENOME_FASTA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="$REF_DIR/gencode.v41.primary_assembly.annotation.gtf"
STAR_INDEX_1PASS="$REF_DIR/star_index_1pass"
STAR_INDEX_2PASS="$REF_DIR/star_index_2pass"

# Sample configuration
SAMPLES_FILE="$PROJECT_ROOT/config/samples.tsv"

# Create directories
mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$REF_DIR"

# === LOGGING ===
LOG_FILE="$LOG_DIR/star_2pass_$(date +%Y%m%d_%H%M%S).log"

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

    # Check reference files
    [[ -f "$GENOME_FASTA" ]] || error_exit "Genome FASTA not found: $GENOME_FASTA"
    [[ -f "$GENOME_GTF" ]] || error_exit "Genome GTF not found: $GENOME_GTF"
    [[ -f "$SAMPLES_FILE" ]] || error_exit "Samples file not found: $SAMPLES_FILE"

    # Read samples
    mapfile -t SAMPLES < <(tail -n +2 "$SAMPLES_FILE" | cut -f1)
    [[ ${#SAMPLES[@]} -gt 0 ]] || error_exit "No samples found in $SAMPLES_FILE"

    log "Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

    # Check trimmed FASTQ files
    for sample in "${SAMPLES[@]}"; do
        local r1_file="$TRIM_DIR/${sample}_R1.paired.fastq.gz"
        local r2_file="$TRIM_DIR/${sample}_R2.paired.fastq.gz"

        [[ -f "$r1_file" ]] || error_exit "R1 file not found for $sample: $r1_file"
        [[ -f "$r2_file" ]] || error_exit "R2 file not found for $sample: $r2_file"
    done

    log "Input validation completed"
}

# === STAR INDEX GENERATION (1ST PASS) ===
build_star_index_1pass() {
    if [[ -d "$STAR_INDEX_1PASS" ]] && [[ -n "$(ls -A "$STAR_INDEX_1PASS" 2>/dev/null)" ]]; then
        log "STAR 1st pass index already exists: $STAR_INDEX_1PASS"
        return 0
    fi

    log "=== Building STAR 1st Pass Index ==="

    mkdir -p "$STAR_INDEX_1PASS"

    local sjdb_overhang=$((READ_LENGTH - 1))

    STAR \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_1PASS" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GENOME_GTF" \
        --sjdbOverhang "$sjdb_overhang" \
        --limitGenomeGenerateRAM $((MEMORY_GB * 1000000000)) \
        --outTmpDir "$STAR_INDEX_1PASS/tmp" \
        > "$LOG_DIR/star_index_1pass.log" 2>&1

    log "STAR 1st pass index completed"
}

# === STAR 1ST PASS ALIGNMENT ===
run_star_1st_pass() {
    log "=== Running STAR 1st Pass Alignment ==="

    for sample in "${SAMPLES[@]}"; do
        log "1st pass alignment for sample: $sample"

        local r1_file="$TRIM_DIR/${sample}_R1.paired.fastq.gz"
        local r2_file="$TRIM_DIR/${sample}_R2.paired.fastq.gz"
        local output_prefix="$ALIGN_DIR/${sample}_1pass_"

        # Check if already completed
        if [[ -f "${output_prefix}SJ.out.tab" ]]; then
            log "1st pass already completed for $sample"
            continue
        fi

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_1PASS" \
            --readFilesIn "$r1_file" "$r2_file" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$output_prefix" \
            --outSAMtype None \
            --outSAMattributes Standard \
            --limitBAMsortRAM $((MEMORY_GB * 1000000000)) \
            > "$LOG_DIR/${sample}_1pass.log" 2>&1

        log "1st pass completed for $sample"
    done

    log "All 1st pass alignments completed"
}

# === COLLECT SPLICE JUNCTIONS ===
collect_splice_junctions() {
    log "=== Collecting Splice Junctions from 1st Pass ==="

    local sj_files=()
    for sample in "${SAMPLES[@]}"; do
        local sj_file="$ALIGN_DIR/${sample}_1pass_SJ.out.tab"
        if [[ -f "$sj_file" ]]; then
            sj_files+=("$sj_file")
        else
            error_exit "Splice junction file not found: $sj_file"
        fi
    done

    log "Found ${#sj_files[@]} splice junction files"

    # Filter splice junctions (GATK recommendations)
    # Keep junctions with at least 3 uniquely mapped reads
    local filtered_sj_dir="$ALIGN_DIR/filtered_sj"
    mkdir -p "$filtered_sj_dir"

    for sj_file in "${sj_files[@]}"; do
        local basename=$(basename "$sj_file" .tab)
        local filtered_file="$filtered_sj_dir/${basename}_filtered.tab"

        # Filter: column 7 (unique reads) >= 3, column 6 (motif) != 0
        awk '$7 >= 3 && $6 != 0' "$sj_file" > "$filtered_file"

        local original_count=$(wc -l < "$sj_file")
        local filtered_count=$(wc -l < "$filtered_file")
        log "$(basename "$sj_file"): $original_count -> $filtered_count junctions"
    done

    log "Splice junction filtering completed"
}

# === STAR INDEX GENERATION (2ND PASS) ===
build_star_index_2pass() {
    if [[ -d "$STAR_INDEX_2PASS" ]] && [[ -n "$(ls -A "$STAR_INDEX_2PASS" 2>/dev/null)" ]]; then
        log "STAR 2nd pass index already exists: $STAR_INDEX_2PASS"
        return 0
    fi

    log "=== Building STAR 2nd Pass Index with Novel Junctions ==="

    mkdir -p "$STAR_INDEX_2PASS"

    # Collect all filtered splice junction files
    local filtered_sj_files=()
    for sample in "${SAMPLES[@]}"; do
        local filtered_file="$ALIGN_DIR/filtered_sj/${sample}_1pass_SJ.out_filtered.tab"
        if [[ -f "$filtered_file" ]]; then
            filtered_sj_files+=("$filtered_file")
        fi
    done

    local sjdb_overhang=$((READ_LENGTH - 1))

    STAR \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_2PASS" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbFileChrStartEnd "${filtered_sj_files[@]}" \
        --sjdbOverhang "$sjdb_overhang" \
        --limitGenomeGenerateRAM $((MEMORY_GB * 1000000000)) \
        --limitSjdbInsertNsj 2000000 \
        --outTmpDir "$STAR_INDEX_2PASS/tmp" \
        > "$LOG_DIR/star_index_2pass.log" 2>&1

    log "STAR 2nd pass index completed"

    # Report junction statistics
    local total_junctions=0
    for sj_file in "${filtered_sj_files[@]}"; do
        local count=$(wc -l < "$sj_file")
        total_junctions=$((total_junctions + count))
    done

    log "Total novel junctions incorporated: $total_junctions"
}

# === STAR 2ND PASS ALIGNMENT ===
run_star_2nd_pass() {
    log "=== Running STAR 2nd Pass Alignment ==="

    for sample in "${SAMPLES[@]}"; do
        log "2nd pass alignment for sample: $sample"

        local r1_file="$TRIM_DIR/${sample}_R1.paired.fastq.gz"
        local r2_file="$TRIM_DIR/${sample}_R2.paired.fastq.gz"
        local output_prefix="$ALIGN_DIR/${sample}_2pass_"

        # Check if already completed
        if [[ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]]; then
            log "2nd pass already completed for $sample"
            continue
        fi

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_2PASS" \
            --readFilesIn "$r1_file" "$r2_file" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$output_prefix" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --limitBAMsortRAM $((MEMORY_GB * 1000000000)) \
            --quantMode TranscriptomeSAM GeneCounts \
            > "$LOG_DIR/${sample}_2pass.log" 2>&1

        # Index the BAM file
        if [[ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]]; then
            samtools index "${output_prefix}Aligned.sortedByCoord.out.bam"
            log "BAM file indexed for $sample"
        fi

        log "2nd pass completed for $sample"
    done

    log "All 2nd pass alignments completed"
}

# === ALIGNMENT STATISTICS ===
generate_alignment_stats() {
    log "=== Generating Alignment Statistics ==="

    local stats_file="$ALIGN_DIR/alignment_statistics.txt"

    cat > "$stats_file" << EOF
# STAR 2-Pass Alignment Statistics
# Generated on: $(date)

Sample	Total_Reads	Uniquely_Mapped	Multi_Mapped	Unmapped	Mapping_Rate
EOF

    for sample in "${SAMPLES[@]}"; do
        local log_file="$ALIGN_DIR/${sample}_2pass_Log.final.out"

        if [[ -f "$log_file" ]]; then
            # Extract key statistics from STAR log
            local total_reads=$(grep "Number of input reads" "$log_file" | awk '{print $NF}')
            local unique_mapped=$(grep "Uniquely mapped reads number" "$log_file" | awk '{print $NF}')
            local multi_mapped=$(grep "Number of reads mapped to multiple loci" "$log_file" | awk '{print $NF}')
            local unmapped=$(grep "Number of reads unmapped: too short" "$log_file" | awk '{print $NF}')

            # Calculate mapping rate
            local mapping_rate="N/A"
            if [[ -n "$total_reads" && -n "$unique_mapped" && "$total_reads" -gt 0 ]]; then
                mapping_rate=$(awk "BEGIN {printf \"%.2f%%\", ($unique_mapped/$total_reads)*100}")
            fi

            echo -e "$sample\t$total_reads\t$unique_mapped\t$multi_mapped\t$unmapped\t$mapping_rate" >> "$stats_file"

            log "$sample: $mapping_rate mapping rate"
        else
            log "Warning: Log file not found for $sample"
        fi
    done

    log "Alignment statistics saved: $stats_file"
}

# === CLEANUP ===
cleanup_intermediate_files() {
    if [[ "${CLEANUP_INTERMEDIATE:-true}" == "true" ]]; then
        log "=== Cleaning Up Intermediate Files ==="

        # Remove 1st pass alignment files (keep splice junctions)
        for sample in "${SAMPLES[@]}"; do
            rm -f "$ALIGN_DIR/${sample}_1pass_Aligned.out.sam"
            rm -f "$ALIGN_DIR/${sample}_1pass_Log.out"
            rm -f "$ALIGN_DIR/${sample}_1pass_Log.progress.out"
        done

        # Remove temporary directories
        rm -rf "$STAR_INDEX_1PASS/tmp" "$STAR_INDEX_2PASS/tmp"

        log "Intermediate file cleanup completed"
    fi
}

# === MAIN EXECUTION ===
main() {
    log "========================================"
    log "STAR 2-Pass Alignment Pipeline Started"
    log "========================================"
    log "Parameters:"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Read length: $READ_LENGTH"
    log "========================================"

    validate_inputs
    build_star_index_1pass
    run_star_1st_pass
    collect_splice_junctions
    build_star_index_2pass
    run_star_2nd_pass
    generate_alignment_stats
    cleanup_intermediate_files

    log "========================================"
    log "STAR 2-Pass Alignment Pipeline Completed"
    log "========================================"
    log "Output files:"
    log "  Final alignments: $ALIGN_DIR/*_2pass_Aligned.sortedByCoord.out.bam"
    log "  Alignment stats: $ALIGN_DIR/alignment_statistics.txt"
    log "  Log files: $LOG_DIR/"
    log ""
    log "Next steps:"
    log "1. Run GATK preprocessing pipeline"
    log "2. Perform gene quantification"
    log "3. Run variant calling"
}

# === COMMAND LINE PARSING ===
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            cat << EOF
STAR 2-Pass Alignment Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -t, --threads INT       Number of threads (default: 8)
    -m, --memory INT        Memory in GB (default: 32)
    -r, --read-length INT   Read length (default: 100)
    --no-cleanup           Keep intermediate files
    --dry-run              Show commands without executing

EXAMPLES:
    # Basic run
    $0

    # Custom resources
    $0 --threads 16 --memory 64

    # Keep intermediate files
    $0 --no-cleanup

DESCRIPTION:
    Implements STAR 2-pass alignment following GATK best practices:
    1. First pass discovers novel splice junctions
    2. Second pass uses augmented index for improved accuracy
    3. Results in high-quality alignments suitable for variant calling

EOF
            exit 0
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY_GB="$2"
            shift 2
            ;;
        -r|--read-length)
            READ_LENGTH="$2"
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