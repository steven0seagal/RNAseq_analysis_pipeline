# RNA-seq Analysis Pipeline - Nextflow Implementation

This directory contains the Nextflow (DSL2) implementation of the RNA-seq analysis pipeline, migrated from the original Snakemake workflows.

## Project Status

### ✅ Completed (Phase 1 - Initial Setup)

1. **Directory Structure** - Complete Nextflow project structure created
2. **Configuration Files** - Core configuration files implemented:
   - `nextflow.config` - Main configuration with parameters and profiles
   - `conf/base.config` - Process resource labels and defaults
   - `conf/modules.config` - Module-specific configurations
3. **Analysis Scripts** - Copied existing Python and R analysis scripts to `bin/`
4. **Core Modules** - Implemented essential nf-core modules:
   - ✅ `FASTQC` - Quality control for FASTQ files
   - ✅ `FASTP` - Adapter trimming and quality filtering
   - ✅ `STAR_GENOMEGENERATE` - STAR genome index generation
   - ✅ `STAR_ALIGN` - Read alignment with STAR
   - ✅ `SAMTOOLS_INDEX` - BAM indexing
   - ✅ `SUBREAD_FEATURECOUNTS` - Gene-level quantification
   - ✅ `MULTIQC` - QC report aggregation
5. **Main Entry Point** - `main.nf` with workflow routing logic

### 🚧 In Progress / To Do

Based on the detailed `nextflow.plan` file, the following major tasks remain:

#### Phase 2: Standard RNA-seq Workflow
- [ ] Create `subworkflows/local/qc_trimming.nf`
- [ ] Create `subworkflows/local/alignment_quantification.nf`
- [ ] Create `subworkflows/local/differential_expression.nf`
- [ ] Create `workflows/rnaseq.nf` (main standard RNA-seq workflow)
- [ ] Create DESeq2 analysis module
- [ ] Test end-to-end standard RNA-seq pipeline

#### Phase 3: circRNA Analysis Workflow
- [ ] Create BWA modules (index, mem)
- [ ] Create circRNA-specific modules:
  - CIRIquant/CIRI2 detection
  - Count matrix merging
  - Differential expression
  - Functional annotation (miRanda)
- [ ] Create `workflows/circrna.nf`
- [ ] Test circRNA workflow

#### Phase 4: Multimodal Analysis Workflow
- [ ] Create GATK4 modules for BAM processing
- [ ] Create Mutect2 variant calling modules
- [ ] Create small RNA analysis modules
- [ ] Create multimodal integration module
- [ ] Create `workflows/multimodal.nf`
- [ ] Test multimodal workflow

#### Phase 5: Testing & Validation
- [ ] Create test profiles with small datasets
- [ ] Implement unit tests for modules
- [ ] Implement integration tests for workflows
- [ ] Compare outputs with Snakemake version
- [ ] Performance benchmarking

#### Phase 6: Documentation
- [ ] Complete usage documentation
- [ ] Document output files
- [ ] Add troubleshooting guide
- [ ] Create tutorial examples

## Directory Structure

```
nextflow/
├── main.nf                      # Main entry point
├── nextflow.config              # Main configuration
├── README.md                    # This file
├── modules/
│   ├── local/                   # Pipeline-specific modules
│   │   ├── circrna/            # circRNA-specific processes
│   │   ├── multimodal/         # Multimodal-specific processes
│   │   ├── deseq2/             # DESeq2 analysis
│   │   └── utils/              # Utility processes
│   └── nf-core/                # Standard nf-core modules
│       ├── fastqc/             # ✅ Implemented
│       ├── fastp/              # ✅ Implemented
│       ├── star/               # ✅ Implemented (genomegenerate, align)
│       ├── samtools/           # ✅ Implemented (index)
│       ├── subread/            # ✅ Implemented (featurecounts)
│       ├── multiqc/            # ✅ Implemented
│       ├── bwa/                # To be implemented
│       ├── gatk4/              # To be implemented
│       ├── salmon/             # To be implemented
│       └── ...
├── subworkflows/
│   └── local/                  # Composed subworkflows
│       ├── qc_trimming.nf      # To be implemented
│       ├── alignment_*.nf      # To be implemented
│       └── ...
├── workflows/
│   ├── rnaseq.nf               # Standard RNA-seq workflow (to be implemented)
│   ├── circrna.nf              # circRNA workflow (to be implemented)
│   ├── multimodal.nf           # Multimodal workflow (to be implemented)
│   └── integrated.nf           # Combined workflow (to be implemented)
├── bin/                         # Executable scripts
│   ├── interpret_results.py    # ✅ Copied from scripts/
│   ├── interpret_results.R     # ✅ Copied from scripts/
│   └── ...
├── conf/                        # Configuration files
│   ├── base.config             # ✅ Process defaults
│   ├── modules.config          # ✅ Module configs
│   ├── test.config             # To be created
│   └── ...
├── assets/                      # Assets and schemas
├── lib/                         # Groovy helper functions
├── docs/                        # Documentation
└── tests/                       # Test workflows
```

## Quick Start (Once Complete)

### Installation

1. Install Nextflow (version ≥23.04.0):
```bash
curl -s https://get.nextflow.io | bash
```

2. Choose your dependency management method:
   - **Conda/Mamba** (recommended for development)
   - **Docker** (recommended for production)
   - **Singularity** (recommended for HPC)

### Running the Pipeline

#### Standard RNA-seq Analysis
```bash
nextflow run main.nf \\
    --input samplesheet.csv \\
    --mode rnaseq \\
    --outdir results \\
    --genome_fasta data/reference/GRCh38.primary_assembly.genome.fa \\
    --annotation_gtf data/reference/gencode.v38.primary_assembly.annotation.gtf \\
    --metadata config/metadata.tsv \\
    -profile conda
```

#### circRNA Analysis
```bash
nextflow run main.nf \\
    --input samplesheet.csv \\
    --mode circrna \\
    --outdir results_circrna \\
    --genome_fasta data/reference/GRCh38.primary_assembly.genome.fa \\
    --annotation_gtf data/reference/gencode.v38.primary_assembly.annotation.gtf \\
    -profile conda
```

#### Multimodal Analysis
```bash
nextflow run main.nf \\
    --input samplesheet.csv \\
    --mode multimodal \\
    --outdir results_multimodal \\
    --genome_fasta data/reference/GRCh38.primary_assembly.genome.fa \\
    --annotation_gtf data/reference/gencode.v38.primary_assembly.annotation.gtf \\
    --dbsnp data/reference/databases/dbsnp_155.hg38.vcf.gz \\
    -profile conda
```

### Input Samplesheet Format

Create a comma-separated samplesheet with the following columns:

```csv
sample,fastq_1,fastq_2,condition
Cancer_1,data/raw_fastq/Cancer_1_R1.fastq.gz,data/raw_fastq/Cancer_1_R2.fastq.gz,Cancer
Cancer_2,data/raw_fastq/Cancer_2_R1.fastq.gz,data/raw_fastq/Cancer_2_R2.fastq.gz,Cancer
Healthy_1,data/raw_fastq/Healthy_1_R1.fastq.gz,data/raw_fastq/Healthy_1_R2.fastq.gz,Healthy
Healthy_2,data/raw_fastq/Healthy_2_R1.fastq.gz,data/raw_fastq/Healthy_2_R2.fastq.gz,Healthy
```

## Configuration

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to input samplesheet | null (required) |
| `--mode` | Analysis mode | 'rnaseq' |
| `--outdir` | Output directory | 'results' |
| `--genome_fasta` | Genome FASTA file | null (required) |
| `--annotation_gtf` | Gene annotation GTF | null (required) |
| `--metadata` | Metadata for DE analysis | 'config/metadata.tsv' |
| `--strand_specificity` | Strandness (0/1/2) | 2 |
| `--padj_threshold` | Adjusted p-value cutoff | 0.05 |
| `--log2fc_threshold` | Log2 fold-change cutoff | 1.0 |

See `nextflow.config` for all available parameters.

### Profiles

- `conda` - Use Conda for dependencies
- `mamba` - Use Mamba (faster than Conda)
- `docker` - Use Docker containers
- `singularity` - Use Singularity containers
- `test` - Run with test data (small dataset)

### Resource Configuration

Resources are allocated based on process labels:
- `process_single`: 1 CPU, 6 GB RAM
- `process_low`: 2 CPUs, 12 GB RAM
- `process_medium`: 6 CPUs, 36 GB RAM
- `process_high`: 12 CPUs, 72 GB RAM

Modify `conf/base.config` to adjust resource allocations.

## Development

### Adding New Modules

1. Create module directory: `modules/nf-core/tool_name/`
2. Create `main.nf` with process definition
3. Add module configuration to `conf/modules.config`
4. Test module independently

Example module structure:
```nextflow
process TOOL_NAME {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::tool=version"
    container "..."

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.output"), emit: output
    path "versions.yml"              , emit: versions

    script:
    """
    tool_command \\
        --input $input_file \\
        --output output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version)
    END_VERSIONS
    """
}
```

### Creating Subworkflows

Subworkflows compose multiple processes:

```nextflow
include { PROCESS_A } from '../../modules/nf-core/tool_a/main'
include { PROCESS_B } from '../../modules/nf-core/tool_b/main'

workflow SUBWORKFLOW_NAME {
    take:
    ch_input

    main:
    PROCESS_A ( ch_input )
    PROCESS_B ( PROCESS_A.out.result )

    emit:
    output = PROCESS_B.out.result
}
```

## Troubleshooting

### Common Issues

1. **Java version issues**: Nextflow requires Java 11+
   ```bash
   java -version  # Check version
   ```

2. **Conda environment conflicts**: Use isolated environments
   ```bash
   nextflow run main.nf -profile conda -with-conda
   ```

3. **Resume not working**: Ensure work directory is preserved
   ```bash
   nextflow run main.nf -resume
   ```

## Comparison with Snakemake Version

### Advantages of Nextflow

- **Better parallelization**: Channel-based data flow
- **Cloud-native**: Built-in support for AWS, Google Cloud, Azure
- **Resumability**: Efficient caching and resume mechanism
- **Portability**: First-class container support
- **Modularity**: Easy to share and reuse modules

### Migration Notes

Key differences from Snakemake:
- Rules → Processes
- Wildcards → Channel metadata
- `expand()` → Channel operators
- File-based → Channel-based data flow
- `config.yaml` → `nextflow.config`

## References

- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/
- **nf-core Guidelines**: https://nf-co.re/developers
- **Original Snakemake Workflows**: `../workflow/`
- **Migration Plan**: `../nextflow.plan`

## Contributing

When contributing to this Nextflow implementation:

1. Follow nf-core style guidelines
2. Add tests for new modules
3. Update documentation
4. Use semantic versioning
5. Follow the migration plan in `../nextflow.plan`

## License

MIT License - see parent directory LICENSE file

## Contact

For questions or issues, please open an issue on GitHub or contact the maintainers.

---

**Status**: 🚧 Early Development - Phase 1 Complete (Setup & Core Modules)

**Last Updated**: 2025-11-02
