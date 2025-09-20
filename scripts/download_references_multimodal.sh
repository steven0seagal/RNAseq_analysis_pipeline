#!/bin/bash

# ==============================================================================
# Enhanced Reference Data Download Script for Multi-Modal RNA-Seq Analysis
#
# Downloads all necessary reference files for cancer genomics analysis:
# - Human genome and comprehensive annotation (GENCODE v41)
# - dbSNP database for GATK BQSR
# - COSMIC cancer mutations database
# - miRBase mature miRNA sequences
# - Additional databases for variant annotation
#
# Usage: ./scripts/download_references_multimodal.sh
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
REF_DIR="$PROJECT_ROOT/01_references"
DB_DIR="$REF_DIR/databases"

# Create reference directories
mkdir -p "$REF_DIR" "$DB_DIR"

# Logging
LOG_FILE="$PROJECT_ROOT/logs/download_references.log"
mkdir -p "$(dirname "$LOG_FILE")"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

log "Starting enhanced reference data download..."

# === PRIMARY REFERENCE GENOME AND ANNOTATION ===
log "=== Downloading Primary Reference Files ==="

# Human genome GRCh38 (primary assembly)
GENOME_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz"
GENOME_FILE="$REF_DIR/GRCh38.primary_assembly.genome.fa"

if [[ ! -f "$GENOME_FILE" ]]; then
    log "Downloading human reference genome (GRCh38)..."
    wget -q --show-progress -O "${GENOME_FILE}.gz" "$GENOME_URL"
    gunzip "${GENOME_FILE}.gz"
    log "Genome download completed"
else
    log "Genome file already exists: $GENOME_FILE"
fi

# GENCODE comprehensive gene annotation v41
ANNOTATION_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz"
ANNOTATION_FILE="$REF_DIR/gencode.v41.primary_assembly.annotation.gtf"

if [[ ! -f "$ANNOTATION_FILE" ]]; then
    log "Downloading GENCODE v41 comprehensive annotation..."
    wget -q --show-progress -O "${ANNOTATION_FILE}.gz" "$ANNOTATION_URL"
    gunzip "${ANNOTATION_FILE}.gz"
    log "Annotation download completed"
else
    log "Annotation file already exists: $ANNOTATION_FILE"
fi

# === VARIANT CALLING REFERENCES ===
log "=== Downloading Variant Calling References ==="

# dbSNP database for GATK BQSR
DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b155_GRCh38p7/VCF/GATK/dbsnp_155.hg38.vcf.gz"
DBSNP_FILE="$DB_DIR/dbsnp_155.hg38.vcf.gz"

if [[ ! -f "$DBSNP_FILE" ]]; then
    log "Downloading dbSNP v155 for GATK BQSR..."
    wget -q --show-progress -O "$DBSNP_FILE" "$DBSNP_URL"

    # Download index
    wget -q --show-progress -O "${DBSNP_FILE}.tbi" "${DBSNP_URL}.tbi"
    log "dbSNP download completed"
else
    log "dbSNP file already exists: $DBSNP_FILE"
fi

# COSMIC cancer mutations database (requires registration - provide instructions)
COSMIC_FILE="$DB_DIR/cosmic_v96_grch38.vcf.gz"
if [[ ! -f "$COSMIC_FILE" ]]; then
    log "NOTE: COSMIC database requires manual download due to licensing:"
    log "1. Register at https://cancer.sanger.ac.uk/cosmic/register"
    log "2. Download CosmicCodingMuts.vcf.gz"
    log "3. Place as: $COSMIC_FILE"
    log "Creating placeholder file..."
    touch "$COSMIC_FILE"
fi

# === SMALL RNA REFERENCES ===
log "=== Downloading Small RNA References ==="

# miRBase mature miRNA sequences
MIRBASE_URL="ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz"
MIRBASE_FILE="$REF_DIR/mature.fa"

if [[ ! -f "$MIRBASE_FILE" ]]; then
    log "Downloading miRBase mature miRNA sequences..."
    wget -q --show-progress -O "${MIRBASE_FILE}.gz" "$MIRBASE_URL"
    gunzip "${MIRBASE_FILE}.gz"
    log "miRBase download completed"
else
    log "miRBase file already exists: $MIRBASE_FILE"
fi

# Human-specific miRNAs only
HUMAN_MIRBASE_FILE="$REF_DIR/human_mature.fa"
if [[ ! -f "$HUMAN_MIRBASE_FILE" && -f "$MIRBASE_FILE" ]]; then
    log "Extracting human miRNAs from miRBase..."
    grep -A1 "^>hsa-" "$MIRBASE_FILE" | grep -v "^--$" > "$HUMAN_MIRBASE_FILE"
    log "Human miRNA extraction completed"
fi

# miRNA hairpin sequences
HAIRPIN_URL="ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz"
HAIRPIN_FILE="$REF_DIR/hairpin.fa"

if [[ ! -f "$HAIRPIN_FILE" ]]; then
    log "Downloading miRBase hairpin sequences..."
    wget -q --show-progress -O "${HAIRPIN_FILE}.gz" "$HAIRPIN_URL"
    gunzip "${HAIRPIN_FILE}.gz"
    log "Hairpin download completed"
else
    log "Hairpin file already exists: $HAIRPIN_FILE"
fi

# === FUNCTIONAL ANNOTATION DATABASES ===
log "=== Setting up Functional Annotation References ==="

# GO and KEGG databases are typically handled by R packages
# Create database information file
cat > "$DB_DIR/database_info.txt" << EOF
# Database Information for Multi-Modal RNA-Seq Analysis
# Generated on: $(date)

## Primary References
Genome: GRCh38.primary_assembly.genome.fa
Annotation: gencode.v41.primary_assembly.annotation.gtf
Source: GENCODE Release 41

## Variant Calling References
dbSNP: dbsnp_155.hg38.vcf.gz (Build 155)
COSMIC: cosmic_v96_grch38.vcf.gz (Version 96 - manual download required)

## Small RNA References
miRBase Mature: mature.fa (Latest release)
miRBase Hairpin: hairpin.fa (Latest release)
Human miRNAs: human_mature.fa (Extracted from mature.fa)

## Functional Databases (R/Bioconductor)
- org.Hs.eg.db: Human genome annotation
- GO.db: Gene Ontology
- KEGG.db: KEGG pathways
- TxDb.Hsapiens.UCSC.hg38.knownGene: Transcript database

## Usage Notes
1. COSMIC database requires registration and manual download
2. Functional databases are installed via R/Bioconductor packages
3. All files use GRCh38/hg38 coordinate system
4. STAR index will be built automatically during pipeline execution
EOF

# === CREATE INDEX PREPARATION SCRIPT ===
log "=== Creating Index Preparation Script ==="

cat > "$REF_DIR/prepare_indices.sh" << 'EOF'
#!/bin/bash

# Index preparation script for multi-modal RNA-seq analysis
# Run this after downloading reference files

set -euo pipefail

REF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Preparing reference indices..."

# FASTA index for GATK
if [[ ! -f "$REF_DIR/GRCh38.primary_assembly.genome.fa.fai" ]]; then
    echo "Creating FASTA index..."
    samtools faidx "$REF_DIR/GRCh38.primary_assembly.genome.fa"
fi

# Sequence dictionary for GATK
if [[ ! -f "$REF_DIR/GRCh38.primary_assembly.genome.dict" ]]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary \
        -R "$REF_DIR/GRCh38.primary_assembly.genome.fa" \
        -O "$REF_DIR/GRCh38.primary_assembly.genome.dict"
fi

# Bowtie index for small RNA analysis
if [[ ! -f "$REF_DIR/bowtie_index/genome.1.ebwt" ]]; then
    echo "Creating Bowtie index for small RNA analysis..."
    mkdir -p "$REF_DIR/bowtie_index"
    bowtie-build "$REF_DIR/GRCh38.primary_assembly.genome.fa" "$REF_DIR/bowtie_index/genome"
fi

echo "Index preparation completed!"
EOF

chmod +x "$REF_DIR/prepare_indices.sh"

# === CREATE VALIDATION SCRIPT ===
cat > "$REF_DIR/validate_references.sh" << 'EOF'
#!/bin/bash

# Validation script for reference files

REF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Validating reference files..."

files=(
    "GRCh38.primary_assembly.genome.fa"
    "gencode.v41.primary_assembly.annotation.gtf"
    "databases/dbsnp_155.hg38.vcf.gz"
    "mature.fa"
    "hairpin.fa"
    "human_mature.fa"
)

all_present=true

for file in "${files[@]}"; do
    if [[ -f "$REF_DIR/$file" ]]; then
        size=$(du -h "$REF_DIR/$file" | cut -f1)
        echo "✓ $file ($size)"
    else
        echo "✗ $file - MISSING"
        all_present=false
    fi
done

if [[ "$all_present" == "true" ]]; then
    echo ""
    echo "✓ All reference files are present!"
    echo "Next steps:"
    echo "1. Run ./prepare_indices.sh to create index files"
    echo "2. Manually download COSMIC database if needed"
    echo "3. Start pipeline execution"
else
    echo ""
    echo "✗ Some reference files are missing. Please re-run download script."
fi
EOF

chmod +x "$REF_DIR/validate_references.sh"

log "=== Reference Download Summary ==="
log "Reference files downloaded to: $REF_DIR"
log "Database files located in: $DB_DIR"
log ""
log "Next steps:"
log "1. Run: $REF_DIR/validate_references.sh"
log "2. Run: $REF_DIR/prepare_indices.sh"
log "3. Manually download COSMIC database if needed"
log ""
log "Reference download script completed!"

# Run validation
"$REF_DIR/validate_references.sh"