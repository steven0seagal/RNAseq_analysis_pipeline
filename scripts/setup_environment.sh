#!/bin/bash

# ==============================================================================
# RNA-Seq Pipeline Environment Setup Script
#
# This script automates the installation and setup of the conda environment
# for the RNA-seq analysis pipeline
#
# Usage: ./scripts/setup_environment.sh [options]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
ENV_NAME="rnaseq-pipeline"

# === COLORS FOR OUTPUT ===
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# === LOGGING FUNCTIONS ===
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# === HELP FUNCTION ===
show_help() {
    cat << EOF
RNA-Seq Pipeline Environment Setup

DESCRIPTION:
    This script sets up the conda environment for the RNA-seq analysis pipeline.
    It will install all required bioinformatics tools and dependencies.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help          Show this help message
    -n, --name NAME     Environment name (default: rnaseq-pipeline)
    -f, --force         Force recreation of existing environment
    -m, --mamba         Use mamba instead of conda (faster)
    --no-test           Skip environment testing
    --minimal           Install minimal environment (core tools only)

EXAMPLES:
    # Basic setup
    $0

    # Force reinstall with custom name
    $0 --force --name my-rnaseq-env

    # Fast setup with mamba
    $0 --mamba

REQUIREMENTS:
    - Conda or Miniconda installed
    - Internet connection for package downloads
    - At least 8GB free disk space

EOF
}

# === COMMAND LINE PARSING ===
FORCE_INSTALL=false
USE_MAMBA=false
SKIP_TESTING=false
MINIMAL_INSTALL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -n|--name)
            ENV_NAME="$2"
            shift 2
            ;;
        -f|--force)
            FORCE_INSTALL=true
            shift
            ;;
        -m|--mamba)
            USE_MAMBA=true
            shift
            ;;
        --no-test)
            SKIP_TESTING=true
            shift
            ;;
        --minimal)
            MINIMAL_INSTALL=true
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# === ENVIRONMENT DETECTION ===
detect_package_manager() {
    if command -v mamba &> /dev/null && [[ "$USE_MAMBA" == "true" ]]; then
        PACKAGE_MANAGER="mamba"
        log_info "Using mamba for package management"
    elif command -v conda &> /dev/null; then
        PACKAGE_MANAGER="conda"
        log_info "Using conda for package management"
    else
        log_error "Neither conda nor mamba found. Please install Miniconda or Anaconda first."
        echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
}

# === ENVIRONMENT MANAGEMENT ===
check_existing_environment() {
    if conda env list | grep -q "^${ENV_NAME} "; then
        if [[ "$FORCE_INSTALL" == "true" ]]; then
            log_warning "Removing existing environment: $ENV_NAME"
            conda env remove -n "$ENV_NAME" -y
            return 1
        else
            log_warning "Environment '$ENV_NAME' already exists."
            echo "Use --force to recreate it, or activate it with:"
            echo "  conda activate $ENV_NAME"
            return 0
        fi
    fi
    return 1
}

# === ENVIRONMENT CREATION ===
create_environment() {
    log_info "Creating conda environment: $ENV_NAME"

    if [[ "$MINIMAL_INSTALL" == "true" ]]; then
        log_info "Installing minimal environment..."
        create_minimal_environment
    else
        log_info "Installing full environment..."
        create_full_environment
    fi
}

create_minimal_environment() {
    $PACKAGE_MANAGER create -n "$ENV_NAME" -y \
        -c conda-forge -c bioconda \
        python=3.9 \
        fastqc=0.11.9 \
        fastp=0.23.2 \
        star=2.7.10a \
        subread=2.0.3 \
        samtools=1.17 \
        multiqc=1.14 \
        snakemake=7.18.2
}

create_full_environment() {
    # Use the environment.yml file
    if [[ -f "$PROJECT_ROOT/environment.yml" ]]; then
        log_info "Using environment.yml file"
        $PACKAGE_MANAGER env create -f "$PROJECT_ROOT/environment.yml" -n "$ENV_NAME"
    else
        log_warning "environment.yml not found, creating basic environment"
        create_minimal_environment
    fi
}

# === ENVIRONMENT TESTING ===
test_environment() {
    if [[ "$SKIP_TESTING" == "true" ]]; then
        log_info "Skipping environment testing"
        return 0
    fi

    log_info "Testing environment installation..."

    # Activate environment for testing
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"

    # Test core bioinformatics tools
    local tools=("fastqc" "fastp" "STAR" "featureCounts" "multiqc" "samtools" "snakemake")
    local failed_tools=()

    for tool in "${tools[@]}"; do
        if command -v "$tool" &> /dev/null; then
            local version
            case $tool in
                "fastqc")
                    version=$($tool --version 2>&1 | head -1)
                    ;;
                "STAR")
                    version=$($tool --version 2>&1)
                    ;;
                "featureCounts")
                    version=$($tool -v 2>&1 | head -1)
                    ;;
                "multiqc")
                    version=$($tool --version 2>&1)
                    ;;
                "samtools")
                    version=$($tool --version 2>&1 | head -1)
                    ;;
                "snakemake")
                    version=$($tool --version 2>&1)
                    ;;
                *)
                    version=$($tool --version 2>&1 | head -1)
                    ;;
            esac
            log_success "$tool: $version"
        else
            log_error "$tool: NOT FOUND"
            failed_tools+=("$tool")
        fi
    done

    # Test Python packages
    log_info "Testing Python packages..."
    python -c "
import sys
packages = ['pandas', 'numpy', 'matplotlib', 'seaborn', 'sklearn']
failed = []

for package in packages:
    try:
        __import__(package)
        print(f'✓ {package}')
    except ImportError:
        print(f'✗ {package}')
        failed.append(package)

if failed:
    print(f'Failed packages: {failed}')
    sys.exit(1)
"

    # Test R packages (if R is installed)
    if command -v R &> /dev/null; then
        log_info "Testing R packages..."
        R --slave -e "
        packages <- c('DESeq2', 'ggplot2', 'pheatmap')
        missing <- packages[!packages %in% rownames(installed.packages())]
        if (length(missing) > 0) {
            cat('Missing R packages:', missing, '\n')
        } else {
            cat('All R packages available\n')
        }
        "
    fi

    if [[ ${#failed_tools[@]} -eq 0 ]]; then
        log_success "All tools installed successfully!"
        return 0
    else
        log_error "Some tools failed to install: ${failed_tools[*]}"
        return 1
    fi
}

# === REFERENCE DATA SETUP ===
setup_reference_data() {
    log_info "Setting up reference data directories..."

    local ref_dir="$PROJECT_ROOT/data/reference"
    mkdir -p "$ref_dir"

    # Create download script
    cat > "$ref_dir/download_references.sh" << 'EOF'
#!/bin/bash
# Script to download reference genome and annotation

echo "Downloading human reference genome (GRCh38)..."
wget -q --show-progress -O GRCh38.primary_assembly.genome.fa.gz \
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"

echo "Downloading gene annotation (GENCODE v38)..."
wget -q --show-progress -O gencode.v38.primary_assembly.annotation.gtf.gz \
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"

echo "Extracting files..."
gunzip *.gz

echo "Reference files downloaded successfully!"
EOF

    chmod +x "$ref_dir/download_references.sh"
    log_success "Reference data setup complete"
    log_info "To download reference files, run: $ref_dir/download_references.sh"
}

# === POST-INSTALLATION SETUP ===
post_installation_setup() {
    log_info "Performing post-installation setup..."

    # Create activation script
    cat > "$PROJECT_ROOT/activate_env.sh" << EOF
#!/bin/bash
# Activation script for RNA-seq pipeline environment

echo "Activating RNA-seq pipeline environment..."
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo "Environment activated! Available tools:"
echo "  - fastqc --version"
echo "  - STAR --version"
echo "  - snakemake --version"
echo ""
echo "To run the pipeline:"
echo "  ./scripts/run_rnaseq_pipeline.sh"
echo ""
echo "To run Snakemake workflow:"
echo "  cd workflow && snakemake --use-conda --cores 8"
EOF

    chmod +x "$PROJECT_ROOT/activate_env.sh"

    # Create example data directory structure
    mkdir -p "$PROJECT_ROOT/data/raw_fastq"
    mkdir -p "$PROJECT_ROOT/data/test_data"

    # Create a simple validation script
    cat > "$PROJECT_ROOT/validate_installation.sh" << EOF
#!/bin/bash
# Validation script for RNA-seq pipeline

echo "Validating RNA-seq pipeline installation..."

# Check environment
if conda env list | grep -q "$ENV_NAME"; then
    echo "✓ Conda environment '$ENV_NAME' exists"
else
    echo "✗ Conda environment '$ENV_NAME' not found"
    exit 1
fi

# Check scripts
scripts=("scripts/run_rnaseq_pipeline.sh" "scripts/interpret_results.R" "scripts/interpret_results.py")
for script in "\${scripts[@]}"; do
    if [[ -x "\$script" ]]; then
        echo "✓ \$script is executable"
    else
        echo "✗ \$script not found or not executable"
    fi
done

# Check configuration files
configs=("config/config.yaml" "config/samples.tsv" "config/metadata.tsv")
for config in "\${configs[@]}"; do
    if [[ -f "\$config" ]]; then
        echo "✓ \$config exists"
    else
        echo "✗ \$config not found"
    fi
done

echo "Validation complete!"
EOF

    chmod +x "$PROJECT_ROOT/validate_installation.sh"

    log_success "Post-installation setup complete"
}

# === MAIN EXECUTION ===
main() {
    log_info "=========================================="
    log_info "RNA-Seq Pipeline Environment Setup"
    log_info "=========================================="
    log_info "Environment name: $ENV_NAME"
    log_info "Force install: $FORCE_INSTALL"
    log_info "Package manager: detecting..."

    # Detect package manager
    detect_package_manager

    # Check for existing environment
    if check_existing_environment; then
        log_success "Environment '$ENV_NAME' is ready to use!"
        log_info "Activate with: conda activate $ENV_NAME"
        exit 0
    fi

    # Create environment
    log_info "Creating new environment..."
    create_environment

    # Test installation
    test_environment

    # Setup reference data
    setup_reference_data

    # Post-installation setup
    post_installation_setup

    # Final instructions
    log_success "=========================================="
    log_success "Installation completed successfully!"
    log_success "=========================================="
    echo ""
    log_info "Next steps:"
    echo "1. Activate the environment:"
    echo "   conda activate $ENV_NAME"
    echo ""
    echo "2. Download reference data (optional):"
    echo "   cd data/reference && ./download_references.sh"
    echo ""
    echo "3. Place your FASTQ files in data/raw_fastq/"
    echo ""
    echo "4. Edit configuration files in config/"
    echo ""
    echo "5. Run the pipeline:"
    echo "   ./scripts/run_rnaseq_pipeline.sh"
    echo ""
    echo "6. Or use Snakemake:"
    echo "   cd workflow && snakemake --use-conda --cores 8"
    echo ""
    log_info "For help, see README.md or run scripts with --help"
}

# === SCRIPT EXECUTION ===
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi