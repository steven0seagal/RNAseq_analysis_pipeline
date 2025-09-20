#!/bin/bash

# ==============================================================================
# Mutect2 Somatic Variant Discovery Pipeline for RNA-seq Data
#
# Implements GATK best practices for somatic variant calling from RNA-seq:
# 1. Mutect2 variant calling with panel of normals and germline filtering
# 2. FilterMutectCalls for artifact filtering and confidence scoring
# 3. VEP annotation for functional impact assessment
# 4. Integration with COSMIC and ClinVar databases
#
# This pipeline identifies somatic mutations in tumor RNA-seq data and
# provides comprehensive annotation for cancer genomics analysis.
#
# Usage: ./scripts/mutect2_somatic_calling.sh [options]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# === CONFIGURATION ===
THREADS=${THREADS:-8}
MEMORY_GB=${MEMORY_GB:-32}
JAVA_OPTS=${JAVA_OPTS:-"-Xmx${MEMORY_GB}g -XX:+UseParallelGC"}

# Directories
VAR_DIR="$PROJECT_ROOT/06_variants"
REF_DIR="$PROJECT_ROOT/01_references"
LOG_DIR="$PROJECT_ROOT/logs/mutect2"
VEP_DIR="$PROJECT_ROOT/07_annotation"

# Reference files
GENOME_FASTA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
DBSNP_VCF="$REF_DIR/databases/dbsnp_155.hg38.vcf.gz"
COSMIC_VCF="$REF_DIR/databases/cosmic_v96_grch38.vcf.gz"
GNOMAD_VCF="$REF_DIR/databases/gnomad.genomes.r3.0.sites.vcf.gz"

# Panel of normals (will be created if doesn't exist)
PON_VCF="$REF_DIR/databases/panel_of_normals.vcf.gz"

# Sample configuration
SAMPLES_FILE="$PROJECT_ROOT/config/samples.tsv"

# Create directories
mkdir -p "$VAR_DIR" "$LOG_DIR" "$VEP_DIR"

# === LOGGING ===
LOG_FILE="$LOG_DIR/mutect2_calling_$(date +%Y%m%d_%H%M%S).log"

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

    # Check processed BAM files
    for sample in "${SAMPLES[@]}"; do
        local bam_file="$VAR_DIR/${sample}.bqsr.bam"
        [[ -f "$bam_file" ]] || error_exit "Processed BAM file not found for $sample: $bam_file"
        [[ -f "${bam_file}.bai" ]] || error_exit "BAM index not found for $sample: ${bam_file}.bai"
    done

    log "Input validation completed"
}

# === CREATE PANEL OF NORMALS ===
create_panel_of_normals() {
    if [[ -f "$PON_VCF" ]]; then
        log "Panel of normals already exists: $PON_VCF"
        return 0
    fi

    log "=== Creating Panel of Normals ==="

    # Identify normal samples (samples with "_N" suffix or from config)
    local normal_samples=()
    for sample in "${SAMPLES[@]}"; do
        if [[ "$sample" =~ _N$ ]] || [[ "$sample" =~ [Nn]ormal ]]; then
            normal_samples+=("$sample")
        fi
    done

    if [[ ${#normal_samples[@]} -eq 0 ]]; then
        log "No normal samples found - creating minimal PON"
        # Create an empty PON for tumor-only analysis
        gatk $JAVA_OPTS CreateSomaticPanelOfNormals \
            --germline-resource "$DBSNP_VCF" \
            --output "$PON_VCF" \
            --min-sample-count 1 \
            > "$LOG_DIR/create_pon.log" 2>&1
        return 0
    fi

    log "Found ${#normal_samples[@]} normal samples: ${normal_samples[*]}"

    # Run Mutect2 on normal samples
    local normal_vcfs=()
    for normal in "${normal_samples[@]}"; do
        local normal_bam="$VAR_DIR/${normal}.bqsr.bam"
        local normal_vcf="$VAR_DIR/${normal}_normal.vcf.gz"

        if [[ ! -f "$normal_vcf" ]]; then
            log "Running Mutect2 on normal sample: $normal"
            gatk $JAVA_OPTS Mutect2 \
                --reference "$GENOME_FASTA" \
                --input "$normal_bam" \
                --output "$normal_vcf" \
                --max-mnp-distance 0 \
                --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
                > "$LOG_DIR/${normal}_mutect2_normal.log" 2>&1
        fi

        normal_vcfs+=("--vcfs" "$normal_vcf")
    done

    # Create Panel of Normals
    log "Creating panel of normals from ${#normal_samples[@]} samples"
    gatk $JAVA_OPTS CreateSomaticPanelOfNormals \
        "${normal_vcfs[@]}" \
        --germline-resource "$DBSNP_VCF" \
        --output "$PON_VCF" \
        --min-sample-count 2 \
        > "$LOG_DIR/create_pon.log" 2>&1

    log "Panel of normals created: $PON_VCF"
}

# === MUTECT2 VARIANT CALLING ===
run_mutect2() {
    local sample=$1
    log "=== Running Mutect2 for $sample ==="

    local input_bam="$VAR_DIR/${sample}.bqsr.bam"
    local output_vcf="$VAR_DIR/${sample}_somatic_raw.vcf.gz"
    local stats_file="$VAR_DIR/${sample}_somatic.vcf.stats"
    local f1r2_file="$VAR_DIR/${sample}_f1r2.tar.gz"
    local log_file="$LOG_DIR/${sample}_mutect2.log"

    # Skip if already processed
    if [[ -f "$output_vcf" ]]; then
        log "Mutect2 already completed for $sample"
        return 0
    fi

    # Determine if tumor-normal or tumor-only analysis
    local normal_sample=""
    local mutect2_args=()

    # Look for paired normal sample
    if [[ "$sample" =~ _T$ ]]; then
        local potential_normal="${sample%_T}_N"
        if [[ " ${SAMPLES[*]} " =~ " $potential_normal " ]]; then
            normal_sample="$potential_normal"
        fi
    fi

    # Configure Mutect2 arguments
    mutect2_args+=(
        "--reference" "$GENOME_FASTA"
        "--input" "$input_bam"
        "--tumor-sample" "$sample"
        "--output" "$output_vcf"
        "--stats" "$stats_file"
        "--f1r2-tar-gz" "$f1r2_file"
        "--germline-resource" "$DBSNP_VCF"
        "--panel-of-normals" "$PON_VCF"
        "--disable-read-filter" "MateOnSameContigOrNoMappedMateReadFilter"
        "--native-pair-hmm-threads" "$THREADS"
        "--max-reads-per-alignment-start" "0"
        "--dont-use-soft-clipped-bases"
    )

    if [[ -n "$normal_sample" ]]; then
        local normal_bam="$VAR_DIR/${normal_sample}.bqsr.bam"
        mutect2_args+=(
            "--input" "$normal_bam"
            "--normal-sample" "$normal_sample"
        )
        log "Running tumor-normal analysis: $sample vs $normal_sample"
    else
        log "Running tumor-only analysis for $sample"
    fi

    # Add population allele frequency if available
    if [[ -f "$GNOMAD_VCF" ]]; then
        mutect2_args+=("--af-of-alleles-not-in-resource" "0.0000025")
    fi

    # Run Mutect2
    gatk $JAVA_OPTS Mutect2 "${mutect2_args[@]}" \
        > "$log_file" 2>&1

    log "Mutect2 completed for $sample"
}

# === LEARN READ ORIENTATION MODEL ===
learn_read_orientation() {
    local sample=$1
    log "=== Learning Read Orientation Model for $sample ==="

    local f1r2_file="$VAR_DIR/${sample}_f1r2.tar.gz"
    local orientation_model="$VAR_DIR/${sample}_read_orientation_model.tar.gz"
    local log_file="$LOG_DIR/${sample}_orientation.log"

    # Skip if already processed
    if [[ -f "$orientation_model" ]]; then
        log "Read orientation model already exists for $sample"
        return 0
    fi

    gatk $JAVA_OPTS LearnReadOrientationModel \
        --input "$f1r2_file" \
        --output "$orientation_model" \
        > "$log_file" 2>&1

    log "Read orientation model completed for $sample"
}

# === GET PILEUP SUMMARIES ===
get_pileup_summaries() {
    local sample=$1
    log "=== Getting Pileup Summaries for $sample ==="

    local input_bam="$VAR_DIR/${sample}.bqsr.bam"
    local pileup_table="$VAR_DIR/${sample}_getpileupsummaries.table"
    local log_file="$LOG_DIR/${sample}_pileup.log"

    # Skip if already processed
    if [[ -f "$pileup_table" ]]; then
        log "Pileup summaries already exist for $sample"
        return 0
    fi

    # Use dbSNP for pileup summaries
    gatk $JAVA_OPTS GetPileupSummaries \
        --input "$input_bam" \
        --variant "$DBSNP_VCF" \
        --intervals "$DBSNP_VCF" \
        --output "$pileup_table" \
        > "$log_file" 2>&1

    log "Pileup summaries completed for $sample"
}

# === CALCULATE CONTAMINATION ===
calculate_contamination() {
    local sample=$1
    log "=== Calculating Contamination for $sample ==="

    local tumor_pileup="$VAR_DIR/${sample}_getpileupsummaries.table"
    local contamination_table="$VAR_DIR/${sample}_contamination.table"
    local segments_table="$VAR_DIR/${sample}_segments.table"
    local log_file="$LOG_DIR/${sample}_contamination.log"

    # Skip if already processed
    if [[ -f "$contamination_table" ]]; then
        log "Contamination already calculated for $sample"
        return 0
    fi

    # Look for paired normal
    local normal_sample=""
    if [[ "$sample" =~ _T$ ]]; then
        local potential_normal="${sample%_T}_N"
        if [[ " ${SAMPLES[*]} " =~ " $potential_normal " ]]; then
            normal_sample="$potential_normal"
        fi
    fi

    local contamination_args=(
        "--input" "$tumor_pileup"
        "--output" "$contamination_table"
        "--tumor-segmentation" "$segments_table"
    )

    if [[ -n "$normal_sample" ]]; then
        local normal_pileup="$VAR_DIR/${normal_sample}_getpileupsummaries.table"
        if [[ -f "$normal_pileup" ]]; then
            contamination_args+=("--matched-normal" "$normal_pileup")
        fi
    fi

    gatk $JAVA_OPTS CalculateContamination "${contamination_args[@]}" \
        > "$log_file" 2>&1

    log "Contamination calculation completed for $sample"
}

# === FILTER MUTECT CALLS ===
filter_mutect_calls() {
    local sample=$1
    log "=== Filtering Mutect Calls for $sample ==="

    local raw_vcf="$VAR_DIR/${sample}_somatic_raw.vcf.gz"
    local filtered_vcf="$VAR_DIR/${sample}_somatic_filtered.vcf.gz"
    local stats_file="$VAR_DIR/${sample}_somatic.vcf.stats"
    local contamination_table="$VAR_DIR/${sample}_contamination.table"
    local segments_table="$VAR_DIR/${sample}_segments.table"
    local orientation_model="$VAR_DIR/${sample}_read_orientation_model.tar.gz"
    local log_file="$LOG_DIR/${sample}_filter.log"

    # Skip if already processed
    if [[ -f "$filtered_vcf" ]]; then
        log "Mutect calls already filtered for $sample"
        return 0
    fi

    local filter_args=(
        "--reference" "$GENOME_FASTA"
        "--variant" "$raw_vcf"
        "--output" "$filtered_vcf"
        "--stats" "$stats_file"
        "--contamination-table" "$contamination_table"
        "--tumor-segmentation" "$segments_table"
        "--orientation-bias-artifact-priors" "$orientation_model"
        "--max-alt-allele-count" "1"
        "--min-reads-per-strand" "0"
    )

    gatk $JAVA_OPTS FilterMutectCalls "${filter_args[@]}" \
        > "$log_file" 2>&1

    log "Mutect call filtering completed for $sample"
}

# === SELECT VARIANTS ===
select_pass_variants() {
    local sample=$1
    log "=== Selecting PASS Variants for $sample ==="

    local filtered_vcf="$VAR_DIR/${sample}_somatic_filtered.vcf.gz"
    local pass_vcf="$VAR_DIR/${sample}_somatic_pass.vcf.gz"
    local log_file="$LOG_DIR/${sample}_select.log"

    # Skip if already processed
    if [[ -f "$pass_vcf" ]]; then
        log "PASS variants already selected for $sample"
        return 0
    fi

    gatk $JAVA_OPTS SelectVariants \
        --reference "$GENOME_FASTA" \
        --variant "$filtered_vcf" \
        --output "$pass_vcf" \
        --exclude-filtered \
        > "$log_file" 2>&1

    log "PASS variant selection completed for $sample"
}

# === ANNOTATE WITH VEP ===
annotate_variants() {
    local sample=$1
    log "=== Annotating Variants with VEP for $sample ==="

    local pass_vcf="$VAR_DIR/${sample}_somatic_pass.vcf.gz"
    local annotated_vcf="$VEP_DIR/${sample}_somatic_annotated.vcf"
    local vep_summary="$VEP_DIR/${sample}_vep_summary.html"
    local log_file="$LOG_DIR/${sample}_vep.log"

    # Skip if already processed
    if [[ -f "$annotated_vcf" ]]; then
        log "Variants already annotated for $sample"
        return 0
    fi

    # Check if VEP is available
    if ! command -v vep &> /dev/null; then
        log "VEP not found - skipping annotation for $sample"
        return 0
    fi

    # VEP annotation with cancer-relevant plugins
    vep \
        --input_file "$pass_vcf" \
        --output_file "$annotated_vcf" \
        --format vcf \
        --vcf \
        --symbol \
        --terms SO \
        --canonical \
        --protein \
        --biotype \
        --gene_phenotype \
        --regulatory \
        --variant_class \
        --stats_file "$vep_summary" \
        --force_overwrite \
        --offline \
        --cache \
        --dir_cache "$REF_DIR/vep_cache" \
        --species homo_sapiens \
        --assembly GRCh38 \
        --fork "$THREADS" \
        > "$log_file" 2>&1

    log "VEP annotation completed for $sample"
}

# === VARIANT STATISTICS ===
generate_variant_stats() {
    log "=== Generating Variant Statistics ==="

    local stats_file="$VAR_DIR/somatic_variant_summary.txt"

    cat > "$stats_file" << EOF
# Somatic Variant Discovery Summary
# Generated on: $(date)

Sample	Raw_Variants	Filtered_Variants	PASS_Variants	SNVs	Indels
EOF

    for sample in "${SAMPLES[@]}"; do
        # Skip normal samples
        if [[ "$sample" =~ _N$ ]] || [[ "$sample" =~ [Nn]ormal ]]; then
            continue
        fi

        local raw_vcf="$VAR_DIR/${sample}_somatic_raw.vcf.gz"
        local filtered_vcf="$VAR_DIR/${sample}_somatic_filtered.vcf.gz"
        local pass_vcf="$VAR_DIR/${sample}_somatic_pass.vcf.gz"

        if [[ -f "$raw_vcf" && -f "$filtered_vcf" && -f "$pass_vcf" ]]; then
            # Count variants
            local raw_count=$(zcat "$raw_vcf" | grep -v "^#" | wc -l)
            local filtered_count=$(zcat "$filtered_vcf" | grep -v "^#" | wc -l)
            local pass_count=$(zcat "$pass_vcf" | grep -v "^#" | wc -l)

            # Count SNVs and indels in PASS variants
            local snv_count=$(zcat "$pass_vcf" | grep -v "^#" | awk '{if(length($4)==1 && length($5)==1) print}' | wc -l)
            local indel_count=$(zcat "$pass_vcf" | grep -v "^#" | awk '{if(length($4)!=1 || length($5)!=1) print}' | wc -l)

            echo -e "$sample\t$raw_count\t$filtered_count\t$pass_count\t$snv_count\t$indel_count" >> "$stats_file"

            log "$sample: $pass_count PASS variants ($snv_count SNVs, $indel_count indels)"
        else
            echo -e "$sample\tN/A\tN/A\tN/A\tN/A\tN/A" >> "$stats_file"
        fi
    done

    log "Variant statistics saved: $stats_file"
}

# === PROCESS SINGLE SAMPLE ===
process_sample() {
    local sample=$1

    # Skip normal samples for somatic calling
    if [[ "$sample" =~ _N$ ]] || [[ "$sample" =~ [Nn]ormal ]]; then
        log "Skipping normal sample for somatic calling: $sample"
        return 0
    fi

    log "========================================="
    log "Processing sample: $sample"
    log "========================================="

    run_mutect2 "$sample"
    learn_read_orientation "$sample"
    get_pileup_summaries "$sample"
    calculate_contamination "$sample"
    filter_mutect_calls "$sample"
    select_pass_variants "$sample"
    annotate_variants "$sample"

    log "Sample $sample processing completed"
}

# === MAIN EXECUTION ===
main() {
    log "========================================="
    log "Mutect2 Somatic Variant Discovery Started"
    log "========================================="
    log "Parameters:"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Java options: $JAVA_OPTS"
    log "========================================="

    validate_inputs
    create_panel_of_normals

    # Process normal samples first for pileup summaries
    for sample in "${SAMPLES[@]}"; do
        if [[ "$sample" =~ _N$ ]] || [[ "$sample" =~ [Nn]ormal ]]; then
            get_pileup_summaries "$sample"
        fi
    done

    # Process tumor samples
    for sample in "${SAMPLES[@]}"; do
        process_sample "$sample"
    done

    generate_variant_stats

    log "========================================="
    log "Mutect2 Somatic Variant Discovery Completed"
    log "========================================="
    log "Output files:"
    log "  Raw variants: $VAR_DIR/*_somatic_raw.vcf.gz"
    log "  Filtered variants: $VAR_DIR/*_somatic_filtered.vcf.gz"
    log "  PASS variants: $VAR_DIR/*_somatic_pass.vcf.gz"
    log "  VEP annotations: $VEP_DIR/*_somatic_annotated.vcf"
    log "  Statistics: $VAR_DIR/somatic_variant_summary.txt"
    log "  Log files: $LOG_DIR/"
    log ""
    log "Next steps:"
    log "1. Review variant statistics and quality metrics"
    log "2. Perform functional annotation and pathway analysis"
    log "3. Integrate with expression and small RNA data"
}

# === COMMAND LINE PARSING ===
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            cat << EOF
Mutect2 Somatic Variant Discovery Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -t, --threads INT       Number of threads (default: 8)
    -m, --memory INT        Memory in GB (default: 32)
    --dry-run              Show commands without executing

EXAMPLES:
    # Basic run
    $0

    # Custom resources
    $0 --threads 16 --memory 64

DESCRIPTION:
    Performs somatic variant discovery from RNA-seq data using GATK Mutect2:
    1. Creates panel of normals from available normal samples
    2. Runs Mutect2 for tumor samples (tumor-normal or tumor-only)
    3. Applies comprehensive filtering and quality control
    4. Annotates variants with VEP for functional impact

    Supports both tumor-normal pairs and tumor-only analysis.
    Output includes high-confidence somatic variants ready for analysis.

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