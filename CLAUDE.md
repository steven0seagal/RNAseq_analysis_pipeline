# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a comprehensive RNA-seq differential expression analysis pipeline supporting multiple analysis modes: standard RNA-seq, circular RNA (circRNA) discovery, and multimodal analysis (somatic variants, small RNA). The pipeline is implemented with both bash scripts and Snakemake workflows for flexibility and scalability.

## Key Commands

### Environment Setup
```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate rnaseq-pipeline

# For multimodal analysis features
conda env create -f environment_multimodal.yml
```

### Running Analysis Pipelines

**Standard RNA-seq Analysis:**
```bash
# Quick bash pipeline
./scripts/run_rnaseq_pipeline.sh [--threads 8] [--strand 2]

# Snakemake workflow (from project root)
cd workflow && snakemake --use-conda --cores 8

# Run from workflow directory
cd workflow && snakemake --use-conda --cores 8 -s Snakefile
```

**circRNA Analysis:**
```bash
# Standalone circRNA pipeline
./scripts/circrna_pipeline.sh [fastq_dir] [samples] [conditions]

# Snakemake circRNA workflow
cd workflow && snakemake --use-conda --cores 8 -s Snakefile_circrna
```

**Multimodal Analysis:**
```bash
# 2-pass STAR alignment for variant calling
./scripts/run_star_2pass.sh

# BAM processing for GATK
./scripts/gatk_bam_processing.sh

# Somatic variant calling
./scripts/mutect2_somatic_calling.sh

# Small RNA analysis
./scripts/small_rna_analysis.sh

# Integrated multimodal pipeline
./scripts/run_integrated_pipeline.sh --mode all

# Multimodal data integration
Rscript scripts/multimodal_integration.R
```

### Analysis and Visualization

**R-based analysis:**
```bash
Rscript scripts/interpret_results.R \
  results/05_featurecounts/raw_counts.tsv \
  config/metadata.tsv \
  results/analysis_R
```

**Python-based analysis:**
```bash
python scripts/interpret_results.py \
  results/05_featurecounts/raw_counts.tsv \
  config/metadata.tsv \
  results/analysis_python
```

**circRNA differential expression:**
```bash
Rscript scripts/circrna_deseq2.R \
  results/circrna_analysis/06_de_analysis/circrna_count_matrix.csv \
  config/metadata.tsv \
  results/circrna_analysis/06_de_analysis
```

### Testing and Validation
```bash
# Create test dataset
./tests/create_test_data.sh

# Validate Snakemake workflow
cd workflow && snakemake --lint
```

## Architecture and Key Components

### Pipeline Structure

The pipeline follows a three-tier architecture:

1. **Bash Scripts Layer** (`scripts/`): Direct execution scripts for quick runs and individual steps
2. **Snakemake Workflow Layer** (`workflow/`): Scalable, reproducible workflows with automatic parallelization
3. **Analysis Scripts Layer** (`scripts/*.py`, `scripts/*.R`): Statistical analysis and visualization in both R and Python

### Core Workflow Steps

**Standard RNA-seq (both bash and Snakemake):**
1. Quality control (FastQC)
2. Adapter trimming (fastp)
3. Post-trimming QC (FastQC)
4. Genome alignment (STAR)
5. Gene quantification (featureCounts)
6. QC aggregation (MultiQC)
7. Differential expression (DESeq2/PyDESeq2)
8. Functional enrichment (clusterProfiler/GSEApy)

**circRNA Analysis Pipeline:**
1. Quality control (FastQC)
2. Adapter trimming (Cutadapt)
3. BWA-MEM alignment for split-read detection
4. circRNA identification (CIRIquant/CIRI2)
5. Back-splice junction quantification
6. Differential expression (DESeq2)
7. Functional annotation (miRNA binding sites, ceRNA networks)

**Multimodal Analysis:**
1. 2-pass STAR alignment
2. GATK BAM processing and variant calling (MuTect2)
3. Small RNA quantification
4. Integrated multi-omics correlation analysis

### Directory Structure

```
RNAseq_analysis_pipeline/
├── scripts/               # Executable bash and analysis scripts
│   ├── run_rnaseq_pipeline.sh      # Main bash pipeline
│   ├── circrna_pipeline.sh         # circRNA analysis
│   ├── run_integrated_pipeline.sh  # Multimodal integration
│   ├── interpret_results.py        # Python DE analysis
│   ├── interpret_results.R         # R DE analysis
│   └── multimodal_integration.R    # Multi-omics integration
├── workflow/              # Snakemake workflows
│   ├── Snakefile                   # Standard RNA-seq
│   ├── Snakefile_circrna          # circRNA workflow
│   ├── Snakefile_multimodal       # Multimodal workflow
│   └── envs/                      # Conda environment specs per tool
├── config/                # Configuration files
│   ├── config.yaml                # Main pipeline parameters
│   ├── config_circrna.yaml        # circRNA-specific settings
│   ├── config_multimodal.yaml     # Multimodal settings
│   ├── samples.tsv                # Sample definitions
│   └── metadata.tsv               # Sample metadata for DE analysis
├── data/                  # Input data (not tracked in git)
│   ├── raw_fastq/                 # Raw FASTQ files (sample_R1.fastq.gz, sample_R2.fastq.gz)
│   └── reference/                 # Reference genome and annotations
└── results/               # Output directory
    ├── 01_fastqc_raw/
    ├── 02_trimmed_fastq/
    ├── 03_fastqc_trimmed/
    ├── 04_star_alignment/
    ├── 05_featurecounts/          # Contains raw_counts.tsv
    ├── 06_multiqc/
    ├── analysis_R/                # R analysis outputs
    ├── analysis_python/           # Python analysis outputs
    ├── circrna_analysis/          # circRNA results
    └── multimodal/                # Multimodal integration results
```

### Configuration Files

**`config/config.yaml`**: Central configuration for standard RNA-seq pipeline
- Reference genome paths
- Tool parameters (fastp, STAR, featureCounts)
- Analysis thresholds (padj, log2FC)
- Resource allocation (threads, memory)

**`config/samples.tsv`**: Sample sheet with format:
```
sample       condition
Cancer_1     Cancer
Healthy_1    Healthy
```

**`config/metadata.tsv`**: Metadata for differential expression analysis with at minimum:
```
sample       condition
```

### Workflow-Specific Notes

**Snakemake Workflows:**
- Must be run from the `workflow/` directory
- Use `--use-conda` flag to activate per-rule conda environments
- Config file paths are relative to workflow directory (use `../config/`)
- Input/output paths are relative to workflow directory (use `../data/`, `../results/`)

**Bash Pipelines:**
- Run from project root directory
- Use environment variables for parameter customization (THREADS, MEMORY, STRAND)
- Automatic logging to `results/logs/`

### Python Analysis Scripts

**`scripts/interpret_results.py`:**
- Performs differential expression with PyDESeq2
- Generates PCA plots, volcano plots, heatmaps
- Optional functional enrichment with GSEApy
- Requires: `pydeseq2`, `gseapy`, `pandas`, `matplotlib`, `seaborn`

**`scripts/generate_visualizations.py`:**
- Creates custom publication-quality figures
- Supports multiple plot types and formats

### R Analysis Scripts

**`scripts/interpret_results.R`:**
- DESeq2-based differential expression analysis
- clusterProfiler functional enrichment (GO, KEGG)
- High-quality visualizations

**`scripts/circrna_deseq2.R`:**
- Specialized differential expression for circRNA count data
- Handles back-splice junction quantifications

**`scripts/multimodal_integration.R`:**
- Integrates multiple data types (RNA-seq, variants, circRNA)
- Network analysis and correlation studies

### Key Input/Output Files

**Inputs:**
- Raw FASTQ: `data/raw_fastq/{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz`
- Reference: `data/reference/GRCh38.primary_assembly.genome.fa`
- Annotation: `data/reference/gencode.v38.primary_assembly.annotation.gtf`

**Key Outputs:**
- Count matrix: `results/05_featurecounts/raw_counts.tsv`
- QC report: `results/06_multiqc/multiqc_report.html`
- DE results: `results/analysis_R/deseq2_results.csv` or `results/analysis_python/deseq2_results.csv`
- Plots: PCA, volcano plots, heatmaps in respective analysis directories

### Important Implementation Details

**Tool-Specific Settings:**
- STAR `sjdbOverhang` should be ReadLength - 1 (default 99 for 100bp reads)
- featureCounts `strand` parameter: 0=unstranded, 1=stranded, 2=reverse-stranded (default 2 for Illumina)
- fastp performs automatic adapter detection by default

**Analysis Thresholds (modifiable in config.yaml):**
- Adjusted p-value: 0.05
- Log2 fold change: 1.0
- Minimum base mean expression: 10

**circRNA-Specific:**
- Uses BWA-MEM for alignment (better split-read detection than STAR)
- CIRIquant for back-splice junction identification
- Requires longer minimum read length (35bp) due to junction spanning

**Multimodal Analysis:**
- Requires 2-pass STAR alignment for improved variant calling
- GATK best practices for BAM processing (MarkDuplicates, AddOrReplaceReadGroups, BaseRecalibrator)
- MuTect2 for somatic variant calling in tumor-normal pairs

### Common Patterns

**Adding New Samples:**
1. Place FASTQ files in `data/raw_fastq/` with naming: `{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz`
2. Add sample entries to `config/samples.tsv`
3. Add metadata to `config/metadata.tsv`

**Custom Comparisons:**
- Edit `config/config.yaml` under `analysis.contrasts` section
- Format: `"Condition1_vs_Condition2"`

**Cluster Execution:**
```bash
# SLURM
cd workflow && snakemake --profile slurm --jobs 100

# SGE
cd workflow && snakemake --cluster "qsub -pe smp {threads}" --jobs 50
```

## Tool Versions

Core tools (see `environment.yml` for complete list):
- STAR: 2.7.10a
- featureCounts (Subread): 2.0.3
- FastQC: 0.11.9
- fastp: 0.23.2
- MultiQC: 1.14
- DESeq2: 1.38.0 (Bioconductor)
- PyDESeq2: 0.3.2
- Snakemake: 7.18.2
- GATK: 4.4.0
- CIRIquant: latest
- Python: 3.9
- R: 4.2.0
