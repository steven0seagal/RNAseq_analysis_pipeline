#!/bin/bash

# ==============================================================================
# Comprehensive RNA-Seq Analysis Pipeline (Bash Script)
#
# Description: This script performs a complete RNA-seq analysis workflow,
#              from raw FASTQ files to a gene count matrix ready for analysis.
#
# Steps:
#   1. Environment validation and setup
#   2. Reference genome download and preparation
#   3. Raw FASTQ Quality Control (FastQC)
#   4. Adapter and Quality Trimming (fastp)
#   5. Post-Trimming Quality Control (FastQC)
#   6. Genome Indexing (STAR) - if not already done
#   7. Read Alignment (STAR)
#   8. Gene-Level Quantification (featureCounts)
#   9. Aggregate QC Reporting (MultiQC)
#   10. Automated interpretation (optional)
#
# Requirements: Conda environment with required bioinformatics tools
# Usage: ./scripts/run_rnaseq_pipeline.sh [options]
# ==============================================================================

set -euo pipefail  # Exit on any error, undefined variable, or pipe failure

# === SCRIPT CONFIGURATION ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"

# === DEFAULT PARAMETERS ===
THREADS=${THREADS:-8}
READ_LENGTH=${READ_LENGTH:-100}
MEMORY=${MEMORY:-32}
STRAND=${STRAND:-2}  # 0=unstranded, 1=stranded, 2=reverse-stranded
MIN_READ_LENGTH=${MIN_READ_LENGTH:-25}
SKIP_DOWNLOAD=${SKIP_DOWNLOAD:-false}
RUN_INTERPRETATION=${RUN_INTERPRETATION:-true}

# === LOGGING SETUP ===
LOG_DIR="$RESULTS_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# === HELP FUNCTION ===
show_help() {
    cat << EOF
RNA-Seq Analysis Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -t, --threads INT       Number of threads to use (default: 8)
    -r, --read-length INT   Read length for STAR index (default: 100)
    -m, --memory INT        Memory in GB (default: 32)
    -s, --strand INT        Strandness: 0=unstranded, 1=stranded, 2=reverse (default: 2)
    --min-length INT        Minimum read length after trimming (default: 25)
    --skip-download         Skip reference genome download
    --no-interpretation     Skip automated result interpretation
    --dry-run              Show commands without executing

EXAMPLES:
    # Basic run with default parameters
    $0

    # Run with 16 threads and custom parameters
    $0 --threads 16 --memory 64 --strand 0

    # Quick run without downloading references
    $0 --skip-download --threads 4

REQUIREMENTS:
    - Conda environment with: fastqc, fastp, star, subread, multiqc, samtools
    - Raw FASTQ files in data/raw_fastq/ directory
    - Sample configuration in config/samples.tsv
    - At least 8GB RAM and 4 CPU cores recommended

OUTPUT:
    Results will be saved in the 'results/' directory with timestamped logs.

EOF
}

# === COMMAND LINE PARSING ===
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -r|--read-length)
            READ_LENGTH="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -s|--strand)
            STRAND="$2"
            shift 2
            ;;
        --min-length)
            MIN_READ_LENGTH="$2"
            shift 2
            ;;
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --no-interpretation)
            RUN_INTERPRETATION=false
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

# === ENVIRONMENT VALIDATION ===
validate_environment() {
    log "=== Validating Environment ==="

    # Check if we're in the correct directory
    if [[ ! -f "$PROJECT_ROOT/README.md" ]]; then
        error_exit "Please run this script from the project root directory"
    fi

    # Check required tools
    local required_tools=("fastqc" "fastp" "STAR" "featureCounts" "multiqc" "samtools")
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            error_exit "Required tool '$tool' not found. Please install conda environment."
        fi
    done

    # Check for sample configuration
    if [[ ! -f "$CONFIG_DIR/samples.tsv" ]]; then
        error_exit "Sample configuration file not found: $CONFIG_DIR/samples.tsv"
    fi

    # Check for input data
    if [[ ! -d "$DATA_DIR/raw_fastq" ]] || [[ -z "$(ls -A "$DATA_DIR/raw_fastq" 2>/dev/null)" ]]; then
        error_exit "No FASTQ files found in $DATA_DIR/raw_fastq/"
    fi

    log "Environment validation completed successfully"
}

# === READ SAMPLE CONFIGURATION ===
read_samples() {
    log "=== Reading Sample Configuration ==="

    # Read samples from config file, skipping header
    mapfile -t SAMPLES < <(tail -n +2 "$CONFIG_DIR/samples.tsv" | cut -f1)

    if [[ ${#SAMPLES[@]} -eq 0 ]]; then
        error_exit "No samples found in $CONFIG_DIR/samples.tsv"
    fi

    log "Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

    # Validate that FASTQ files exist for all samples
    for sample in "${SAMPLES[@]}"; do
        local r1_file="$DATA_DIR/raw_fastq/${sample}_R1.fastq.gz"
        local r2_file="$DATA_DIR/raw_fastq/${sample}_R2.fastq.gz"

        if [[ ! -f "$r1_file" ]] || [[ ! -f "$r2_file" ]]; then
            error_exit "Missing FASTQ files for sample: $sample"
        fi
    done

    log "All FASTQ files validated"
}

# === REFERENCE GENOME SETUP ===
setup_reference() {
    log "=== Setting Up Reference Genome ==="

    local ref_dir="$DATA_DIR/reference"
    mkdir -p "$ref_dir"

    local genome_fa="$ref_dir/GRCh38.primary_assembly.genome.fa"
    local annotation_gtf="$ref_dir/gencode.v38.primary_assembly.annotation.gtf"

    if [[ "$SKIP_DOWNLOAD" == "false" ]]; then
        # Download genome if not present
        if [[ ! -f "$genome_fa" ]]; then
            log "Downloading human reference genome (GRCh38)..."
            wget -q --show-progress -O "${genome_fa}.gz" \
                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"
            gunzip "${genome_fa}.gz"
        fi

        # Download annotation if not present
        if [[ ! -f "$annotation_gtf" ]]; then
            log "Downloading gene annotation (GENCODE v38)..."
            wget -q --show-progress -O "${annotation_gtf}.gz" \
                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
            gunzip "${annotation_gtf}.gz"
        fi
    fi

    # Validate reference files exist
    if [[ ! -f "$genome_fa" ]] || [[ ! -f "$annotation_gtf" ]]; then
        error_exit "Reference files not found. Run without --skip-download to download them."
    fi

    log "Reference genome setup completed"
}

# === CREATE OUTPUT DIRECTORIES ===
create_directories() {
    log "=== Creating Output Directories ==="

    local dirs=(
        "$RESULTS_DIR/01_fastqc_raw"
        "$RESULTS_DIR/02_trimmed_fastq"
        "$RESULTS_DIR/03_fastqc_trimmed"
        "$RESULTS_DIR/04_star_alignment"
        "$RESULTS_DIR/05_featurecounts"
        "$RESULTS_DIR/06_multiqc"
        "$LOG_DIR"
    )

    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
    done

    log "Output directories created"
}

# === STEP 1: RAW FASTQC ===
run_fastqc_raw() {
    log "=== Step 1: Running FastQC on Raw FASTQ Files ==="

    for sample in "${SAMPLES[@]}"; do
        log "Processing sample: $sample"

        local r1_file="$DATA_DIR/raw_fastq/${sample}_R1.fastq.gz"
        local r2_file="$DATA_DIR/raw_fastq/${sample}_R2.fastq.gz"
        local output_dir="$RESULTS_DIR/01_fastqc_raw"

        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            echo "DRY RUN: fastqc -t $THREADS -o $output_dir $r1_file $r2_file"
        else
            fastqc -t "$THREADS" -o "$output_dir" "$r1_file" "$r2_file" 2>&1 | tee -a "$LOG_FILE"
        fi
    done

    log "Raw FastQC analysis completed"
}

# === STEP 2: ADAPTER TRIMMING ===
run_trimming() {
    log "=== Step 2: Adapter and Quality Trimming with fastp ==="

    for sample in "${SAMPLES[@]}"; do
        log "Trimming sample: $sample"

        local input_r1="$DATA_DIR/raw_fastq/${sample}_R1.fastq.gz"
        local input_r2="$DATA_DIR/raw_fastq/${sample}_R2.fastq.gz"
        local output_r1="$RESULTS_DIR/02_trimmed_fastq/${sample}_R1.trimmed.fastq.gz"
        local output_r2="$RESULTS_DIR/02_trimmed_fastq/${sample}_R2.trimmed.fastq.gz"
        local html_report="$RESULTS_DIR/02_trimmed_fastq/${sample}.fastp.html"
        local json_report="$RESULTS_DIR/02_trimmed_fastq/${sample}.fastp.json"

        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            echo "DRY RUN: fastp -i $input_r1 -I $input_r2 -o $output_r1 -O $output_r2 --html $html_report --json $json_report -l $MIN_READ_LENGTH -w $THREADS"
        else
            fastp \
                -i "$input_r1" \
                -I "$input_r2" \
                -o "$output_r1" \
                -O "$output_r2" \
                --html "$html_report" \
                --json "$json_report" \
                -l "$MIN_READ_LENGTH" \
                -w "$THREADS" \
                2>&1 | tee -a "$LOG_FILE"
        fi
    done

    log "Adapter trimming completed"
}

# === STEP 3: POST-TRIMMING FASTQC ===
run_fastqc_trimmed() {
    log "=== Step 3: Running FastQC on Trimmed FASTQ Files ==="

    for sample in "${SAMPLES[@]}"; do
        log "QC check for trimmed sample: $sample"

        local r1_file="$RESULTS_DIR/02_trimmed_fastq/${sample}_R1.trimmed.fastq.gz"
        local r2_file="$RESULTS_DIR/02_trimmed_fastq/${sample}_R2.trimmed.fastq.gz"
        local output_dir="$RESULTS_DIR/03_fastqc_trimmed"

        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            echo "DRY RUN: fastqc -t $THREADS -o $output_dir $r1_file $r2_file"
        else
            fastqc -t "$THREADS" -o "$output_dir" "$r1_file" "$r2_file" 2>&1 | tee -a "$LOG_FILE"
        fi
    done

    log "Post-trimming FastQC analysis completed"
}

# === STEP 4: STAR GENOME INDEXING ===
build_star_index() {
    log "=== Step 4: Building STAR Genome Index ==="

    local star_index_dir="$DATA_DIR/reference/star_index"
    local genome_fa="$DATA_DIR/reference/GRCh38.primary_assembly.genome.fa"
    local annotation_gtf="$DATA_DIR/reference/gencode.v38.primary_assembly.annotation.gtf"

    if [[ -d "$star_index_dir" ]] && [[ -n "$(ls -A "$star_index_dir" 2>/dev/null)" ]]; then
        log "STAR index already exists, skipping indexing"
        return 0
    fi

    mkdir -p "$star_index_dir"

    local sjdb_overhang=$((READ_LENGTH - 1))

    if [[ "${DRY_RUN:-false}" == "true" ]]; then
        echo "DRY RUN: STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $star_index_dir --genomeFastaFiles $genome_fa --sjdbGTFfile $annotation_gtf --sjdbOverhang $sjdb_overhang"
    else
        log "Building STAR index (this may take 30-60 minutes)..."
        STAR \
            --runThreadN "$THREADS" \
            --runMode genomeGenerate \
            --genomeDir "$star_index_dir" \
            --genomeFastaFiles "$genome_fa" \
            --sjdbGTFfile "$annotation_gtf" \
            --sjdbOverhang "$sjdb_overhang" \
            --limitGenomeGenerateRAM $((MEMORY * 1000000000)) \
            2>&1 | tee -a "$LOG_FILE"
    fi

    log "STAR genome indexing completed"
}

# === STEP 5: READ ALIGNMENT ===
run_star_alignment() {
    log "=== Step 5: Aligning Reads with STAR ==="

    local star_index_dir="$DATA_DIR/reference/star_index"

    for sample in "${SAMPLES[@]}"; do
        log "Aligning sample: $sample"

        local input_r1="$RESULTS_DIR/02_trimmed_fastq/${sample}_R1.trimmed.fastq.gz"
        local input_r2="$RESULTS_DIR/02_trimmed_fastq/${sample}_R2.trimmed.fastq.gz"
        local output_prefix="$RESULTS_DIR/04_star_alignment/${sample}_"

        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            echo "DRY RUN: STAR --runThreadN $THREADS --genomeDir $star_index_dir --readFilesIn $input_r1 $input_r2 --readFilesCommand zcat --outFileNamePrefix $output_prefix --outSAMtype BAM SortedByCoordinate"
        else
            STAR \
                --runThreadN "$THREADS" \
                --genomeDir "$star_index_dir" \
                --readFilesIn "$input_r1" "$input_r2" \
                --readFilesCommand zcat \
                --outFileNamePrefix "$output_prefix" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outFilterType BySJout \
                --outFilterMultimapNmax 20 \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                2>&1 | tee -a "$LOG_FILE"
        fi
    done

    log "STAR alignment completed"
}

# === STEP 6: GENE QUANTIFICATION ===
run_feature_counts() {
    log "=== Step 6: Gene Quantification with featureCounts ==="

    local annotation_gtf="$DATA_DIR/reference/gencode.v38.primary_assembly.annotation.gtf"
    local output_file="$RESULTS_DIR/05_featurecounts/raw_counts.tsv"

    # Collect all BAM files
    local bam_files=()
    for sample in "${SAMPLES[@]}"; do
        bam_files+=("$RESULTS_DIR/04_star_alignment/${sample}_Aligned.sortedByCoord.out.bam")
    done

    if [[ "${DRY_RUN:-false}" == "true" ]]; then
        echo "DRY RUN: featureCounts -T $THREADS -p -s $STRAND -a $annotation_gtf -o $output_file ${bam_files[*]}"
    else
        featureCounts \
            -T "$THREADS" \
            -p \
            -s "$STRAND" \
            -a "$annotation_gtf" \
            -o "$output_file" \
            "${bam_files[@]}" \
            2>&1 | tee -a "$LOG_FILE"
    fi

    log "Gene quantification completed"
}

# === STEP 7: MULTIQC REPORT ===
run_multiqc() {
    log "=== Step 7: Generating MultiQC Report ==="

    if [[ "${DRY_RUN:-false}" == "true" ]]; then
        echo "DRY RUN: multiqc $RESULTS_DIR -o $RESULTS_DIR/06_multiqc --force"
    else
        multiqc "$RESULTS_DIR" -o "$RESULTS_DIR/06_multiqc" --force 2>&1 | tee -a "$LOG_FILE"
    fi

    log "MultiQC report generated"
}

# === STEP 8: AUTOMATED INTERPRETATION ===
run_interpretation() {
    if [[ "$RUN_INTERPRETATION" == "false" ]]; then
        log "Skipping automated interpretation"
        return 0
    fi

    log "=== Step 8: Running Automated Result Interpretation ==="

    local counts_file="$RESULTS_DIR/05_featurecounts/raw_counts.tsv"
    local metadata_file="$CONFIG_DIR/metadata.tsv"

    # Check if interpretation scripts exist
    local r_script="$SCRIPT_DIR/interpret_results.R"
    local python_script="$SCRIPT_DIR/interpret_results.py"

    if [[ -f "$r_script" ]] && [[ -f "$metadata_file" ]]; then
        log "Running R-based analysis..."
        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            echo "DRY RUN: Rscript $r_script $counts_file $metadata_file $RESULTS_DIR/analysis_R"
        else
            Rscript "$r_script" "$counts_file" "$metadata_file" "$RESULTS_DIR/analysis_R" 2>&1 | tee -a "$LOG_FILE" || log "R analysis failed or not available"
        fi
    fi

    if [[ -f "$python_script" ]] && [[ -f "$metadata_file" ]]; then
        log "Running Python-based analysis..."
        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            echo "DRY RUN: python $python_script $counts_file $metadata_file $RESULTS_DIR/analysis_python"
        else
            python "$python_script" "$counts_file" "$metadata_file" "$RESULTS_DIR/analysis_python" 2>&1 | tee -a "$LOG_FILE" || log "Python analysis failed or not available"
        fi
    fi

    log "Automated interpretation completed"
}

# === PIPELINE SUMMARY ===
print_summary() {
    log "=== Pipeline Execution Summary ==="
    log "Start time: $START_TIME"
    log "End time: $(date)"
    log "Total runtime: $(($(date +%s) - $(date -d "$START_TIME" +%s))) seconds"
    log "Samples processed: ${#SAMPLES[@]}"
    log "Results directory: $RESULTS_DIR"
    log "Log file: $LOG_FILE"

    # Count output files
    local count_files=$(find "$RESULTS_DIR" -type f | wc -l)
    log "Total output files generated: $count_files"

    # Check for key output files
    local key_files=(
        "$RESULTS_DIR/05_featurecounts/raw_counts.tsv"
        "$RESULTS_DIR/06_multiqc/multiqc_report.html"
    )

    log "Key output files:"
    for file in "${key_files[@]}"; do
        if [[ -f "$file" ]]; then
            log "  ✓ $(basename "$file")"
        else
            log "  ✗ $(basename "$file") - MISSING"
        fi
    done

    log ""
    log "=== Next Steps ==="
    log "1. Review MultiQC report: $RESULTS_DIR/06_multiqc/multiqc_report.html"
    log "2. Examine count matrix: $RESULTS_DIR/05_featurecounts/raw_counts.tsv"
    log "3. Run differential expression analysis if not already done"
    log "4. Check analysis results in: $RESULTS_DIR/analysis_*/"
    log ""
    log "Pipeline completed successfully!"
}

# === MAIN EXECUTION ===
main() {
    START_TIME=$(date)

    log "========================================"
    log "RNA-Seq Analysis Pipeline Started"
    log "========================================"
    log "Version: 1.0.0"
    log "Start time: $START_TIME"
    log "Parameters:"
    log "  Threads: $THREADS"
    log "  Read length: $READ_LENGTH"
    log "  Memory: ${MEMORY}GB"
    log "  Strand: $STRAND"
    log "  Min read length: $MIN_READ_LENGTH"
    log "  Skip download: $SKIP_DOWNLOAD"
    log "  Run interpretation: $RUN_INTERPRETATION"
    log "========================================"

    # Execute pipeline steps
    validate_environment
    read_samples
    setup_reference
    create_directories
    run_fastqc_raw
    run_trimming
    run_fastqc_trimmed
    build_star_index
    run_star_alignment
    run_feature_counts
    run_multiqc
    run_interpretation

    print_summary
}

# === SCRIPT EXECUTION ===
# Check if script is being sourced or executed
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi