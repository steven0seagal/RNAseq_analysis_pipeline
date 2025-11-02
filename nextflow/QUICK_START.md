# Quick Start Guide - Nextflow RNA-seq Pipeline

## 🚀 Ready to Run!

The standard RNA-seq workflow is **fully functional** and ready for testing.

## Installation

### 1. Install Nextflow

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Move to your PATH (optional)
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

**Required**: Nextflow ≥23.04.0, Java ≥11

### 2. Choose Execution Mode

Pick one profile based on your environment:

- **`conda`** - Uses Conda for dependencies (recommended for development)
- **`docker`** - Uses Docker containers (recommended for production)
- **`singularity`** - Uses Singularity (recommended for HPC)

## Basic Usage

### Step 1: Prepare Input Samplesheet

Create a CSV file with your samples:

```csv
sample,fastq_1,fastq_2,condition
Cancer_1,/path/to/Cancer_1_R1.fastq.gz,/path/to/Cancer_1_R2.fastq.gz,Cancer
Cancer_2,/path/to/Cancer_2_R1.fastq.gz,/path/to/Cancer_2_R2.fastq.gz,Cancer
Healthy_1,/path/to/Healthy_1_R1.fastq.gz,/path/to/Healthy_1_R2.fastq.gz,Healthy
Healthy_2,/path/to/Healthy_2_R1.fastq.gz,/path/to/Healthy_2_R2.fastq.gz,Healthy
```

See `assets/samplesheet_example.csv` for reference.

### Step 2: Prepare Metadata (for DESeq2)

Create a TSV file for differential expression:

```tsv
sample	condition
Cancer_1	Cancer
Cancer_2	Cancer
Healthy_1	Healthy
Healthy_2	Healthy
```

### Step 3: Run the Pipeline

```bash
cd nextflow/

nextflow run main.nf \
    --input samples.csv \
    --mode rnaseq \
    --outdir results \
    --genome_fasta ../data/reference/GRCh38.primary_assembly.genome.fa \
    --annotation_gtf ../data/reference/gencode.v38.primary_assembly.annotation.gtf \
    --metadata metadata.tsv \
    -profile conda
```

### Step 4: Resume After Interruption

If the pipeline stops, resume with:

```bash
nextflow run main.nf \
    --input samples.csv \
    --mode rnaseq \
    -profile conda \
    -resume
```

## Common Use Cases

### Test Run (Help Message)

```bash
nextflow run main.nf --help
```

### Skip Specific Steps

```bash
# Skip FastQC
nextflow run main.nf --input samples.csv --mode rnaseq --skip_fastqc -profile conda

# Skip trimming
nextflow run main.nf --input samples.csv --mode rnaseq --skip_trimming -profile conda

# Skip DESeq2
nextflow run main.nf --input samples.csv --mode rnaseq --skip_deseq2 -profile conda
```

### Custom Parameters

```bash
nextflow run main.nf \
    --input samples.csv \
    --mode rnaseq \
    --outdir my_results \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf \
    --metadata metadata.tsv \
    --strand_specificity 2 \
    --padj_threshold 0.01 \
    --log2fc_threshold 2.0 \
    --control_condition Healthy \
    --treatment_condition Cancer \
    --max_cpus 16 \
    --max_memory 64.GB \
    -profile conda \
    -resume
```

### Using Docker

```bash
nextflow run main.nf \
    --input samples.csv \
    --mode rnaseq \
    -profile docker
```

### Using Singularity (HPC)

```bash
nextflow run main.nf \
    --input samples.csv \
    --mode rnaseq \
    -profile singularity
```

## Output Structure

```
results/
├── fastqc/                       # Quality control reports
│   ├── *_fastqc.html
│   └── *_fastqc.zip
├── fastp/                        # Trimming reports
│   ├── *.fastp.html
│   └── *.fastp.json
├── star/                         # Alignment files
│   ├── *.bam
│   ├── *.bam.bai
│   └── *Log.final.out
├── featurecounts/                # Gene counts
│   ├── *.featureCounts.txt
│   └── *.featureCounts.txt.summary
├── multiqc/                      # Aggregated QC
│   └── multiqc_report.html
├── deseq2/                       # Differential expression
│   ├── deseq2_results.csv
│   ├── deseq2_results_significant.csv
│   ├── up_regulated_genes.csv
│   ├── down_regulated_genes.csv
│   ├── pca_plot.png
│   ├── volcano_plot.png
│   ├── heatmap_top50.png
│   └── ma_plot.png
└── pipeline_info/                # Execution info
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.html
```

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Samplesheet CSV | Required |
| `--mode` | Analysis mode | `rnaseq` |
| `--outdir` | Output directory | `results` |
| `--genome_fasta` | Genome FASTA file | Required |
| `--annotation_gtf` | Gene annotation GTF | Required |
| `--metadata` | Metadata TSV for DE | `config/metadata.tsv` |
| `--strand_specificity` | Strandness (0/1/2) | `2` |
| `--padj_threshold` | Adjusted p-value cutoff | `0.05` |
| `--log2fc_threshold` | Log2FC cutoff | `1.0` |
| `--control_condition` | Control group name | `Healthy` |
| `--treatment_condition` | Treatment group name | `Cancer` |

## Troubleshooting

### Issue: "Command not found: nextflow"
**Solution**: Install Nextflow or add to PATH

### Issue: "No such file or directory"
**Solution**: Check that paths in samplesheet are absolute or relative to run directory

### Issue: Conda environment conflicts
**Solution**: Use `mamba` profile instead:
```bash
nextflow run main.nf -profile mamba
```

### Issue: Out of memory
**Solution**: Increase memory allocation:
```bash
nextflow run main.nf --max_memory 128.GB
```

### Issue: Need to restart failed step
**Solution**: Use `-resume` flag:
```bash
nextflow run main.nf -resume
```

## Advanced Features

### Check Pipeline Syntax
```bash
nextflow run main.nf --help
```

### Dry Run (Show What Would Execute)
```bash
nextflow run main.nf -preview
```

### View Execution DAG
After running, check: `results/pipeline_info/pipeline_dag_*.html`

### Clean Work Directory
```bash
# Remove temporary work files
nextflow clean -f

# Keep only last run
nextflow clean -f -k
```

## Next Steps

- Review MultiQC report: `results/multiqc/multiqc_report.html`
- Check DESeq2 results: `results/deseq2/deseq2_results.csv`
- View visualizations: `results/deseq2/*.png`
- Examine execution reports: `results/pipeline_info/`

## Need Help?

- **Documentation**: See `README.md` for detailed information
- **Implementation Details**: See `IMPLEMENTATION_STATUS.md`
- **Full Roadmap**: See `../nextflow.plan`
- **Original Pipeline**: Check `../workflow/Snakefile` for comparison

## What's Working Now

✅ **Fully Functional**:
- Standard RNA-seq workflow (FASTQ → DESeq2)
- Quality control (FastQC, fastp, MultiQC)
- Alignment (STAR)
- Quantification (featureCounts)
- Differential expression (DESeq2 with plots)

🚧 **Coming Soon**:
- circRNA analysis workflow
- Multimodal analysis workflow
- Integrated workflow

---

**Ready to test?** Start with a small dataset and use `-resume` for quick iterations!
