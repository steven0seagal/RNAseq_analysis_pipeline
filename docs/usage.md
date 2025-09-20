# Usage Guide

This guide provides comprehensive instructions for using the RNA-Seq Analysis Pipeline for differential gene expression analysis.

## Table of Contents

- [Quick Start](#quick-start)
- [Data Preparation](#data-preparation)
- [Configuration](#configuration)
- [Pipeline Execution](#pipeline-execution)
- [Output Analysis](#output-analysis)
- [Advanced Usage](#advanced-usage)
- [Examples](#examples)

## Quick Start

### 1. Activate Environment

```bash
conda activate rnaseq-pipeline
# or
source activate_env.sh
```

### 2. Prepare Your Data

```bash
# Place FASTQ files in data/raw_fastq/
# Format: {sample}_R1.fastq.gz, {sample}_R2.fastq.gz
cp /path/to/your/fastq/* data/raw_fastq/

# Edit sample configuration
vim config/samples.tsv
vim config/metadata.tsv
```

### 3. Run Pipeline

```bash
# Option A: Bash pipeline (simple)
./scripts/run_rnaseq_pipeline.sh

# Option B: Snakemake workflow (advanced)
cd workflow
snakemake --use-conda --cores 8

# Option C: Multimodal analysis
./scripts/run_star_2pass.sh  # 2-pass STAR alignment
./scripts/multimodal_integration.R  # Integrated analysis

# Option D: Docker container
docker-compose up
```

### 4. Analyze Results

```bash
# View quality control report
open results/06_multiqc/multiqc_report.html

# Check differential expression results
ls results/analysis_R/
ls results/analysis_python/
```

## Data Preparation

### Input Data Requirements

#### FASTQ Files
- **Format**: Paired-end, gzipped FASTQ files
- **Naming**: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- **Quality**: Illumina quality scores (Phred+33)
- **Location**: `data/raw_fastq/` directory

Example file structure:
```
data/raw_fastq/
├── Cancer_1_R1.fastq.gz
├── Cancer_1_R2.fastq.gz
├── Cancer_2_R1.fastq.gz
├── Cancer_2_R2.fastq.gz
├── Healthy_1_R1.fastq.gz
├── Healthy_1_R2.fastq.gz
└── ...
```

#### Reference Data
- **Genome**: Human GRCh38 (auto-downloaded)
- **Annotation**: GENCODE v38 (auto-downloaded)
- **Custom**: Edit `config/config.yaml` for other organisms

### Sample Metadata

#### samples.tsv Format
```tsv
sample	condition	batch	replicate	notes
Cancer_1	Cancer	1	1	Tumor tissue sample
Cancer_2	Cancer	1	2	Tumor tissue sample
Cancer_3	Cancer	2	3	Tumor tissue sample
Healthy_1	Healthy	1	1	Normal tissue control
Healthy_2	Healthy	1	2	Normal tissue control
Healthy_3	Healthy	2	3	Normal tissue control
```

#### metadata.tsv Format
```tsv
sample_id	condition	batch	replicate	seq_date	library_prep	notes
Cancer_1	Cancer	1	1	2024-01-15	TruSeq	Patient C001
Cancer_2	Cancer	1	2	2024-01-15	TruSeq	Patient C002
Healthy_1	Healthy	1	1	2024-01-15	TruSeq	Control H001
```

### Data Quality Considerations

#### Minimum Requirements
- **Read length**: ≥50 bp
- **Read depth**: ≥20M reads per sample
- **Replicates**: ≥3 per condition
- **Quality**: >80% bases with Q30+

#### Recommended Specifications
- **Read length**: 100-150 bp
- **Read depth**: 30-50M reads per sample
- **Replicates**: 3-6 per condition
- **Strand specificity**: Known and documented

## Configuration

### Global Configuration (config.yaml)

Key parameters to review:

```yaml
# Computational resources
resources:
  threads: 8          # Adjust based on your system
  memory_gb: 32       # Adjust based on available RAM

# Analysis parameters
params:
  fastp:
    min_len: 25       # Minimum read length after trimming

  star:
    sjdbOverhang: 99  # Read length - 1

  featureCounts:
    strand: 2         # 0=unstranded, 1=stranded, 2=reverse
```

### Sample-Specific Configuration

#### Batch Effects
If your samples were processed in batches:

```yaml
# In config.yaml
advanced:
  batch_correction:
    enabled: true
    batch_column: "batch"
```

#### Strand Specificity
Determine library preparation protocol:

| Protocol | Strand Setting |
|----------|----------------|
| Unstranded | 0 |
| Forward stranded | 1 |
| Reverse stranded (TruSeq) | 2 |

### Custom Reference Genomes

For non-human organisms, edit `config.yaml`:

```yaml
reference:
  genome_fasta: "data/reference/custom_genome.fa"
  annotation_gtf: "data/reference/custom_annotation.gtf"
```

## Pipeline Execution

### Method 1: Bash Script

#### Basic Usage
```bash
./scripts/run_rnaseq_pipeline.sh
```

#### Advanced Options
```bash
# Show help
./scripts/run_rnaseq_pipeline.sh --help

# Custom parameters
./scripts/run_rnaseq_pipeline.sh \
    --threads 16 \
    --memory 64 \
    --strand 0 \
    --min-length 30

# Skip reference download (if already present)
./scripts/run_rnaseq_pipeline.sh --skip-download

# Dry run (show commands without executing)
./scripts/run_rnaseq_pipeline.sh --dry-run
```

### Method 2: Snakemake Workflow

#### Basic Usage
```bash
cd workflow
snakemake --use-conda --cores 8
```

#### Advanced Options
```bash
# Dry run to check workflow
snakemake --dry-run

# Generate workflow visualization
snakemake --dag | dot -Tpng > workflow_dag.png

# Run specific rule
snakemake results/05_featurecounts/raw_counts.tsv --use-conda

# Cluster execution (SLURM)
snakemake --profile slurm --jobs 100

# Clean up intermediate files
snakemake clean --use-conda
```

### Method 3: Docker

#### Interactive Mode
```bash
# Start container
docker run -it --rm \
    -v $(pwd)/data:/app/data \
    -v $(pwd)/results:/app/results \
    rnaseq-pipeline bash

# Inside container
./scripts/run_rnaseq_pipeline.sh
```

#### Docker Compose
```bash
# Start all services
docker-compose up -d

# Execute pipeline
docker-compose exec pipeline ./scripts/run_rnaseq_pipeline.sh

# Execute multimodal pipeline
docker-compose exec pipeline ./scripts/run_star_2pass.sh

# Access Jupyter notebook
open http://localhost:8888

# Stop services
docker-compose down
```

## Monitoring Progress

### Log Files
```bash
# Check main pipeline log
tail -f results/logs/pipeline_*.log

# Check individual tool logs
tail -f results/logs/star_align/Cancer_1.log
tail -f results/logs/fastqc_raw/Cancer_1.log
```

### Resource Usage
```bash
# Monitor CPU and memory
htop

# Check disk usage
df -h
du -sh results/

# Monitor specific processes
ps aux | grep STAR
ps aux | grep snakemake
```

### Snakemake Progress
```bash
# Show rule execution status
snakemake --summary

# Show failed jobs
snakemake --detailed-summary

# Rerun failed jobs
snakemake --rerun-incomplete
```

## Output Analysis

### Directory Structure
```
results/
├── 01_fastqc_raw/          # Raw data QC
├── 02_trimmed_fastq/       # Cleaned reads
├── 03_fastqc_trimmed/      # Post-trim QC
├── 04_star_alignment/      # Aligned reads
├── 05_featurecounts/       # Gene counts
├── 06_multiqc/            # Summary QC report
├── analysis_R/            # R analysis results
├── analysis_python/       # Python analysis results
├── multimodal/            # Multimodal analysis results
│   ├── variants/          # Somatic variants
│   ├── small_rna/         # Small RNA analysis
│   └── integration/       # Integrated results
└── logs/                  # Pipeline logs
```

### Key Output Files

#### Quality Control
- `06_multiqc/multiqc_report.html`: Comprehensive QC summary
- `01_fastqc_raw/`: Individual FastQC reports for raw data
- `03_fastqc_trimmed/`: Individual FastQC reports after trimming

#### Expression Data
- `05_featurecounts/raw_counts.tsv`: Gene count matrix
- `05_featurecounts/raw_counts.tsv.summary`: Counting statistics

#### Differential Expression
- `analysis_R/deseq2_full_results.csv`: Complete DESeq2 results
- `analysis_R/up_regulated_genes.csv`: Upregulated genes
- `analysis_R/down_regulated_genes.csv`: Downregulated genes

#### Visualizations
- `analysis_R/pca_plot.png`: Principal component analysis
- `analysis_R/volcano_plot.png`: Volcano plot
- `analysis_R/heatmap_top50_degs.png`: Heatmap of top DEGs

#### Functional Analysis
- `analysis_R/GO_enrichment_*.csv`: Gene Ontology results
- `analysis_R/KEGG_enrichment_*.csv`: KEGG pathway results

#### Multimodal Analysis
- `multimodal/variants/somatic_variants.vcf`: Somatic variants
- `multimodal/small_rna/small_rna_counts.tsv`: Small RNA expression
- `multimodal/integration/multimodal_report.html`: Integrated analysis report

### Interpreting Results

#### Quality Control Metrics

**Good Quality Indicators:**
- >90% reads with Q30+
- <5% adapter contamination
- >80% alignment rate
- Even distribution across gene features

**Warning Signs:**
- <80% reads with Q30+
- >10% adapter contamination
- <70% alignment rate
- Batch effects in PCA plot

#### Differential Expression

**Statistical Thresholds:**
- Adjusted p-value < 0.05
- |log2 fold change| > 1.0
- Base mean > 10

**Result Interpretation:**
- **Positive log2FC**: Upregulated in test vs control
- **Negative log2FC**: Downregulated in test vs control
- **padj**: Multiple testing corrected p-value

## Advanced Usage

### Multimodal Analysis

#### Somatic Variant Calling
```bash
# Prepare paired tumor-normal samples
# Edit config/pairs.tsv with sample pairs

# Run somatic calling
./scripts/mutect2_somatic_calling.sh

# Process BAM files with GATK
./scripts/gatk_bam_processing.sh
```

#### Small RNA Analysis
```bash
# Run small RNA pipeline
./scripts/small_rna_analysis.sh

# Generate custom visualizations
python scripts/generate_visualizations.py
```

#### Multimodal Integration
```bash
# Integrate RNA-seq with other omics data
Rscript scripts/multimodal_integration.R

# Use multimodal configuration
# Edit config/config_multimodal.yaml
```

### Custom Analysis Scripts

#### Running R Analysis Separately
```bash
# Activate environment
conda activate rnaseq-pipeline

# Run R analysis
Rscript scripts/interpret_results.R \
    results/05_featurecounts/raw_counts.tsv \
    config/metadata.tsv \
    results/custom_analysis_R
```

#### Running Python Analysis Separately
```bash
# Run Python analysis
python scripts/interpret_results.py \
    results/05_featurecounts/raw_counts.tsv \
    config/metadata.tsv \
    results/custom_analysis_python
```

### Cluster Computing

#### SLURM Configuration
```bash
# Create SLURM profile
mkdir -p ~/.config/snakemake/slurm

# Edit cluster configuration
# See Snakemake documentation for details
```

#### PBS/Torque
```bash
# Submit to PBS queue
qsub -l nodes=1:ppn=8,walltime=24:00:00 run_pipeline.sh
```

### Cloud Computing

#### AWS Batch
```bash
# Use AWS Batch for large-scale analysis
# Configure AWS credentials and submit jobs
```

#### Google Cloud
```bash
# Use Google Cloud Life Sciences API
# Submit pipeline to managed compute environments
```

### Multiple Comparisons

For complex experimental designs with multiple conditions:

1. **Edit samples.tsv** to include all conditions
2. **Modify config.yaml** to specify contrasts:
```yaml
analysis:
  contrasts:
    - "Treatment_vs_Control"
    - "TimePoint2_vs_TimePoint1"
    - "Condition1_vs_Condition2"
```

### Time Series Analysis

For time-course experiments:

1. **Organize samples** by time points
2. **Use specialized R packages** (DESeq2 LRT test)
3. **Modify analysis scripts** for temporal patterns

## Troubleshooting

### Common Issues

#### 1. Out of Memory Errors
```bash
# Reduce memory usage
export STAR_limitGenomeGenerateRAM=16000000000

# Or edit config.yaml
resources:
  memory_gb: 16
```

#### 2. Disk Space Issues
```bash
# Clean intermediate files
rm -rf results/01_fastqc_raw/
rm -rf results/02_trimmed_fastq/

# Or use Snakemake cleanup
snakemake clean
```

#### 3. Tool Not Found
```bash
# Ensure environment is activated
conda activate rnaseq-pipeline

# Reinstall missing tools
conda install -c bioconda fastqc
```

### Performance Optimization

#### I/O Optimization
```bash
# Use fast local storage for temporary files
export TMPDIR=/fast/local/storage

# Parallelize I/O operations
snakemake --cores 16 --max-jobs-per-second 1
```

#### Memory Optimization
```bash
# Monitor memory usage
free -h
htop

# Adjust parameters for available memory
# Edit config.yaml accordingly
```

## Best Practices

### Reproducibility
1. **Version control**: Use git for tracking changes
2. **Environment management**: Pin software versions
3. **Documentation**: Record all parameters and modifications
4. **Data backup**: Maintain copies of raw data

### Quality Control
1. **Always check MultiQC report** before proceeding
2. **Validate sample metadata** for accuracy
3. **Examine PCA plots** for batch effects
4. **Review alignment statistics** for quality

### Data Management
1. **Organize directories** consistently
2. **Use meaningful sample names**
3. **Document experimental design**
4. **Archive intermediate results**

## Getting Help

For additional support:

1. Check the [FAQ](faq.md)
2. Browse [GitHub Issues](https://github.com/username/RNAseq_analysis_pipeline/issues)
3. Join [Community Discussions](https://github.com/username/RNAseq_analysis_pipeline/discussions)
4. Read the [Troubleshooting Guide](troubleshooting.md)

## Next Steps

After completing your analysis:

1. **Validate results** with additional experiments
2. **Perform functional analysis** of significant genes
3. **Generate publication figures** from output plots
4. **Share your workflow** for reproducible research