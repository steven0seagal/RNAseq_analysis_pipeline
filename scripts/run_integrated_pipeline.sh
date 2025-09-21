#!/bin/bash
set -e

#================================================================================
# INTEGRATED RNA-SEQ AND circRNA ANALYSIS PIPELINE
#================================================================================
# This script integrates both standard RNA-seq and circRNA analysis
# combining the existing pipeline with the new circRNA functionality
#
# USAGE:
# ./scripts/run_integrated_pipeline.sh [OPTIONS]
#
# OPTIONS:
#   --mode         Analysis mode: rna, circrna, multimodal, or all (default: all)
#   --samples      Sample names (comma-separated)
#   --conditions   Condition names (comma-separated)
#   --threads      Number of threads to use (default: 8)
#   --help         Show this help message
#================================================================================

# Default parameters
MODE="all"
THREADS=8
SAMPLES=""
CONDITIONS=""
HELP=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)
            MODE="$2"
            shift 2
            ;;
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        --conditions)
            CONDITIONS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Show help if requested
if [[ "$HELP" == "true" ]]; then
    cat << EOF
INTEGRATED RNA-SEQ AND circRNA ANALYSIS PIPELINE

USAGE:
    ./scripts/run_integrated_pipeline.sh [OPTIONS]

OPTIONS:
    --mode         Analysis mode:
                   - rna: Standard RNA-seq analysis only
                   - circrna: circRNA analysis only
                   - multimodal: Both RNA-seq and circRNA
                   - all: Complete multimodal analysis (default)
    --samples      Sample names (comma-separated)
                   Example: "Cancer_1,Cancer_2,Healthy_1,Healthy_2"
    --conditions   Condition names (comma-separated)
                   Example: "Cancer,Cancer,Healthy,Healthy"
    --threads      Number of threads to use (default: 8)
    --help         Show this help message

EXAMPLES:
    # Run complete analysis
    ./scripts/run_integrated_pipeline.sh --mode all

    # Run only circRNA analysis
    ./scripts/run_integrated_pipeline.sh --mode circrna

    # Run with custom samples
    ./scripts/run_integrated_pipeline.sh --samples "S1,S2,S3" --conditions "A,A,B"

REQUIREMENTS:
    - Conda environment with required tools
    - Input FASTQ files in data/raw_fastq/
    - Reference files in data/reference/

OUTPUT:
    - Standard RNA-seq results: results/
    - circRNA analysis results: results/circrna/
    - Integrated analysis: results/multimodal/
EOF
    exit 0
fi

echo "==================================================================="
echo "INTEGRATED RNA-SEQ AND circRNA ANALYSIS PIPELINE"
echo "==================================================================="
echo "Mode: $MODE"
echo "Threads: $THREADS"
echo "Samples: ${SAMPLES:-Auto-detect}"
echo "Conditions: ${CONDITIONS:-Auto-detect}"
echo ""

# Function to check if conda environment is active
check_environment() {
    if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
        echo "Warning: No conda environment detected"
        echo "Please activate the rnaseq-pipeline environment:"
        echo "  conda activate rnaseq-pipeline"
        echo ""
    else
        echo "Active conda environment: $CONDA_DEFAULT_ENV"
    fi
}

# Function to check required tools
check_tools() {
    local tools=("fastqc" "cutadapt" "bwa" "samtools" "snakemake")
    local missing_tools=()

    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        echo "Error: Missing required tools: ${missing_tools[*]}"
        echo "Please install missing tools or activate the correct conda environment"
        exit 1
    fi

    echo "✓ All required tools are available"
}

# Function to check input files
check_inputs() {
    if [[ ! -d "data/raw_fastq" ]]; then
        echo "Error: data/raw_fastq directory not found"
        echo "Please create the directory and place your FASTQ files there"
        exit 1
    fi

    local fastq_files=(data/raw_fastq/*.fastq.gz)
    if [[ ${#fastq_files[@]} -eq 0 ]] || [[ ! -f "${fastq_files[0]}" ]]; then
        echo "Error: No FASTQ files found in data/raw_fastq/"
        echo "Please place your paired-end FASTQ files in data/raw_fastq/"
        echo "Expected format: {sample}_R1.fastq.gz, {sample}_R2.fastq.gz"
        exit 1
    fi

    echo "✓ Found FASTQ files in data/raw_fastq/"
}

# Function to auto-detect samples if not provided
auto_detect_samples() {
    if [[ -z "$SAMPLES" ]]; then
        echo "Auto-detecting samples from FASTQ files..."
        SAMPLES=$(ls data/raw_fastq/*_R1.fastq.gz | sed 's|data/raw_fastq/||g' | sed 's|_R1.fastq.gz||g' | tr '\n' ',' | sed 's/,$//')
        echo "Detected samples: $SAMPLES"
    fi

    if [[ -z "$CONDITIONS" ]]; then
        echo "Auto-detecting conditions from sample names..."
        # Simple heuristic: if sample contains "cancer", "tumor", "case" -> Cancer
        # if sample contains "healthy", "normal", "control" -> Healthy
        local sample_array
        IFS=',' read -r -a sample_array <<< "$SAMPLES"
        local condition_array=()

        for sample in "${sample_array[@]}"; do
            if [[ $sample =~ [Cc]ancer|[Tt]umor|[Cc]ase ]]; then
                condition_array+=("Cancer")
            elif [[ $sample =~ [Hh]ealthy|[Nn]ormal|[Cc]ontrol ]]; then
                condition_array+=("Healthy")
            else
                # Default assignment based on position
                if [[ ${#condition_array[@]} -lt $((${#sample_array[@]} / 2)) ]]; then
                    condition_array+=("Condition_A")
                else
                    condition_array+=("Condition_B")
                fi
            fi
        done

        CONDITIONS=$(IFS=',' ; echo "${condition_array[*]}")
        echo "Detected conditions: $CONDITIONS"
    fi
}

# Function to run standard RNA-seq analysis
run_rnaseq_analysis() {
    echo ""
    echo "--- Running Standard RNA-seq Analysis ---"

    # Check if results already exist
    if [[ -f "results/06_multiqc/multiqc_report.html" ]]; then
        echo "Standard RNA-seq analysis results already exist. Skipping..."
        return
    fi

    # Run Snakemake workflow
    cd workflow
    echo "Running Snakemake RNA-seq workflow..."
    snakemake --use-conda --cores "$THREADS" --printshellcmds

    cd ..
    echo "✓ Standard RNA-seq analysis completed"
}

# Function to run circRNA analysis
run_circrna_analysis() {
    echo ""
    echo "--- Running circRNA Analysis ---"

    # Check if results already exist
    if [[ -f "results/circrna/analysis_summary.txt" ]]; then
        echo "circRNA analysis results already exist. Skipping..."
        return
    fi

    # Option 1: Run bash script
    echo "Running circRNA analysis (bash script)..."
    export THREADS="$THREADS"
    export REF_GENOME="data/reference/GRCh38.primary_assembly.genome.fa"
    export GTF_ANNOTATION="data/reference/gencode.v38.primary_assembly.annotation.gtf"
    export MIRNA_FASTA="data/reference/mature_mirna.fa"

    ./scripts/circrna_pipeline.sh "data/raw_fastq" "$SAMPLES" "$CONDITIONS"

    # Option 2: Run Snakemake workflow (alternative)
    # cd workflow
    # echo "Running Snakemake circRNA workflow..."
    # snakemake --use-conda --cores "$THREADS" --printshellcmds -s Snakefile_circrna
    # cd ..

    echo "✓ circRNA analysis completed"
}

# Function to run multimodal integration
run_multimodal_integration() {
    echo ""
    echo "--- Running Multimodal Integration ---"

    # Check if we have both RNA-seq and circRNA results
    if [[ ! -f "results/05_featurecounts/raw_counts.tsv" ]]; then
        echo "Warning: Standard RNA-seq results not found. Running RNA-seq analysis first..."
        run_rnaseq_analysis
    fi

    if [[ ! -f "results/circrna/05_ciri_quant/merged_count_matrix.tsv" ]]; then
        echo "Warning: circRNA results not found. Running circRNA analysis first..."
        run_circrna_analysis
    fi

    # Run multimodal integration script
    mkdir -p results/multimodal
    if [[ -f "scripts/multimodal_integration.R" ]]; then
        echo "Running multimodal integration..."
        Rscript scripts/multimodal_integration.R
        echo "✓ Multimodal integration completed"
    else
        echo "Warning: Multimodal integration script not found"
    fi
}

# Function to generate final report
generate_report() {
    echo ""
    echo "--- Generating Final Report ---"

    REPORT_FILE="results/integrated_analysis_report.html"

    cat > "$REPORT_FILE" << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Integrated RNA-seq and circRNA Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }
        .header { background-color: #f0f8ff; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
        .results { background-color: #f9f9f9; }
        .warning { background-color: #fff3cd; border-color: #ffeaa7; }
        h1, h2, h3 { color: #2c3e50; }
        .file-link { color: #3498db; text-decoration: none; }
        .file-link:hover { text-decoration: underline; }
        ul { padding-left: 20px; }
        .status { font-weight: bold; }
        .status.success { color: #27ae60; }
        .status.warning { color: #f39c12; }
        .status.error { color: #e74c3c; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Integrated RNA-seq and circRNA Analysis Report</h1>
        <p><strong>Generated:</strong> {TIMESTAMP}</p>
        <p><strong>Analysis Mode:</strong> {MODE}</p>
        <p><strong>Samples:</strong> {SAMPLES}</p>
    </div>

    <div class="section">
        <h2>📊 Analysis Summary</h2>
        <p>This report summarizes the results of an integrated RNA-seq and circRNA analysis pipeline,
        combining differential gene expression analysis with circular RNA discovery and functional annotation.</p>
    </div>

    <div class="section results">
        <h2>📁 Key Results</h2>

        <h3>Standard RNA-seq Analysis</h3>
        <ul>
            <li><a href="06_multiqc/multiqc_report.html" class="file-link">Quality Control Report (MultiQC)</a></li>
            <li><a href="05_featurecounts/raw_counts.tsv" class="file-link">Gene Expression Count Matrix</a></li>
            <li><a href="analysis_R/deseq2_full_results.csv" class="file-link">Differential Expression Results</a></li>
            <li><a href="analysis_R/volcano_plot.png" class="file-link">Volcano Plot</a></li>
            <li><a href="analysis_R/pca_plot.png" class="file-link">PCA Plot</a></li>
        </ul>

        <h3>circRNA Analysis</h3>
        <ul>
            <li><a href="circrna/05_ciri_quant/merged_count_matrix.tsv" class="file-link">circRNA Count Matrix</a></li>
            <li><a href="circrna/06_de_analysis/DE_results.csv" class="file-link">Differential circRNA Expression</a></li>
            <li><a href="circrna/06_de_analysis/volcano_plot.pdf" class="file-link">circRNA Volcano Plot</a></li>
            <li><a href="circrna/07_functional_annotation/miranda_predictions.txt" class="file-link">miRNA Binding Predictions</a></li>
            <li><a href="circrna/analysis_summary.txt" class="file-link">circRNA Analysis Summary</a></li>
        </ul>

        <h3>Multimodal Integration</h3>
        <ul>
            <li><a href="multimodal/" class="file-link">Integrated Analysis Results</a></li>
        </ul>
    </div>

    <div class="section">
        <h2>🔬 Analysis Pipeline</h2>
        <p>The integrated pipeline consists of the following stages:</p>
        <ol>
            <li><strong>Quality Control:</strong> FastQC analysis of raw and trimmed reads</li>
            <li><strong>Read Processing:</strong> Adapter trimming and quality filtering</li>
            <li><strong>Alignment:</strong>
                <ul>
                    <li>STAR alignment for standard RNA-seq</li>
                    <li>BWA-MEM alignment for circRNA detection</li>
                </ul>
            </li>
            <li><strong>Quantification:</strong>
                <ul>
                    <li>featureCounts for gene expression</li>
                    <li>CIRIquant for circRNA identification and quantification</li>
                </ul>
            </li>
            <li><strong>Differential Analysis:</strong> DESeq2 for both mRNAs and circRNAs</li>
            <li><strong>Functional Annotation:</strong> miRNA binding site prediction for circRNAs</li>
            <li><strong>Integration:</strong> Multi-modal data integration and visualization</li>
        </ol>
    </div>

    <div class="section">
        <h2>📋 Next Steps</h2>
        <ul>
            <li>Review quality control reports to ensure data quality</li>
            <li>Examine differential expression results for both mRNAs and circRNAs</li>
            <li>Investigate functional annotations for significant circRNAs</li>
            <li>Perform pathway enrichment analysis on differentially expressed genes</li>
            <li>Validate interesting candidates experimentally (qPCR, Western blot, etc.)</li>
            <li>Consider integration with other omics data if available</li>
        </ul>
    </div>

    <div class="section warning">
        <h2>⚠️ Important Notes</h2>
        <ul>
            <li>All computational predictions require experimental validation</li>
            <li>Statistical significance does not guarantee biological relevance</li>
            <li>Consider batch effects if samples were processed at different times</li>
            <li>Pathway analysis may be limited by current knowledge in databases</li>
        </ul>
    </div>

    <div class="section">
        <h2>📚 References</h2>
        <ul>
            <li>Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.</li>
            <li>Gao, Y., Zhang, J., Zhao, F. (2018). Circular RNA identification based on multiple seed matching. Briefings in Bioinformatics, 19(5), 803-810.</li>
            <li>Zhang, X.O., Dong, R., Zhang, Y., et al. (2016). Diverse alternative back-splicing and alternative splicing landscape of circular RNAs. Genome Research, 26(9), 1277-1287.</li>
        </ul>
    </div>
</body>
</html>
EOF

    # Replace placeholders
    sed -i "s/{TIMESTAMP}/$(date)/" "$REPORT_FILE"
    sed -i "s/{MODE}/$MODE/" "$REPORT_FILE"
    sed -i "s/{SAMPLES}/$SAMPLES/" "$REPORT_FILE"

    echo "✓ Final report generated: $REPORT_FILE"
}

# Main execution
main() {
    echo "Starting integrated analysis pipeline..."

    # Pre-flight checks
    check_environment
    check_tools
    check_inputs
    auto_detect_samples

    # Create results directory
    mkdir -p results

    # Run analysis based on mode
    case "$MODE" in
        "rna")
            run_rnaseq_analysis
            ;;
        "circrna")
            run_circrna_analysis
            ;;
        "multimodal")
            run_rnaseq_analysis
            run_circrna_analysis
            run_multimodal_integration
            ;;
        "all")
            run_rnaseq_analysis
            run_circrna_analysis
            run_multimodal_integration
            ;;
        *)
            echo "Error: Unknown mode '$MODE'"
            echo "Valid modes: rna, circrna, multimodal, all"
            exit 1
            ;;
    esac

    # Generate final report
    generate_report

    echo ""
    echo "==================================================================="
    echo "✅ INTEGRATED ANALYSIS COMPLETED SUCCESSFULLY"
    echo "==================================================================="
    echo "Results are available in: results/"
    echo "View the report: results/integrated_analysis_report.html"
    echo ""
    echo "Summary of outputs:"
    if [[ "$MODE" =~ ^(rna|multimodal|all)$ ]]; then
        echo "  - RNA-seq results: results/"
    fi
    if [[ "$MODE" =~ ^(circrna|multimodal|all)$ ]]; then
        echo "  - circRNA results: results/circrna/"
    fi
    if [[ "$MODE" =~ ^(multimodal|all)$ ]]; then
        echo "  - Integrated results: results/multimodal/"
    fi
    echo ""
}

# Run main function
main