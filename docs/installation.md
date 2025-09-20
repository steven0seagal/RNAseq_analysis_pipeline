# Installation Guide

This guide provides detailed instructions for installing and setting up the RNA-Seq Analysis Pipeline on different systems.

## Table of Contents

- [System Requirements](#system-requirements)
- [Quick Installation](#quick-installation)
- [Manual Installation](#manual-installation)
- [Docker Installation](#docker-installation)
- [Troubleshooting](#troubleshooting)
- [Verification](#verification)

## System Requirements

### Minimum Requirements
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+) or macOS (10.14+)
- **RAM**: 8 GB minimum, 32 GB recommended
- **CPU**: 4 cores minimum, 8+ cores recommended
- **Storage**: 50 GB free space minimum, 200 GB recommended
- **Internet**: Required for downloading reference genomes and software

### Software Dependencies
- **Conda/Miniconda**: Package manager for bioinformatics tools
- **Python**: 3.8+ (installed via conda)
- **R**: 4.0+ (installed via conda)
- **Docker**: Optional, for containerized deployment

## Quick Installation

### Option 1: Automated Setup (Recommended)

```bash
# Clone the repository
git clone https://github.com/username/RNAseq_analysis_pipeline.git
cd RNAseq_analysis_pipeline

# Run automated setup
./scripts/setup_environment.sh

# Activate environment
conda activate rnaseq-pipeline
```

### Option 2: Using Conda Environment File

```bash
# Clone repository
git clone https://github.com/username/RNAseq_analysis_pipeline.git
cd RNAseq_analysis_pipeline

# Create environment from file
conda env create -f environment.yml

# Activate environment
conda activate rnaseq-pipeline
```

## Manual Installation

### Step 1: Install Miniconda

If you don't have conda installed:

```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Initialize conda
source ~/.bashrc
```

For macOS:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

### Step 2: Add Bioconda Channels

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Step 3: Create Environment

```bash
# Create new environment
conda create -n rnaseq-pipeline python=3.9

# Activate environment
conda activate rnaseq-pipeline
```

### Step 4: Install Core Tools

```bash
# Quality control tools
conda install -c bioconda fastqc=0.11.9 fastp=0.23.2 multiqc=1.14

# Alignment tools
conda install -c bioconda star=2.7.10a samtools=1.17

# Quantification tools
conda install -c bioconda subread=2.0.3 salmon=1.9.0

# Workflow management
conda install -c bioconda snakemake=7.18.2

# Data science packages
conda install pandas numpy matplotlib seaborn scikit-learn
```

### Step 5: Install R Packages

```bash
# Install R and Bioconductor packages
conda install r-base=4.2.0 r-essentials
conda install -c bioconda bioconductor-deseq2 bioconductor-clusterprofiler
conda install r-ggplot2 r-pheatmap r-ggrepel
```

### Step 6: Install Python Bioinformatics Packages

```bash
# Install via pip (within conda environment)
pip install pydeseq2 gseapy biopython pysam
```

## Docker Installation

### Prerequisites
- Docker installed and running
- Docker Compose (optional, for multi-container setup)

### Option 1: Pull Pre-built Image

```bash
# Pull from GitHub Container Registry (when available)
docker pull ghcr.io/username/rnaseq_analysis_pipeline:latest

# Run container
docker run -it --rm -v $(pwd)/data:/app/data ghcr.io/username/rnaseq_analysis_pipeline:latest
```

### Option 2: Build from Source

```bash
# Clone repository
git clone https://github.com/username/RNAseq_analysis_pipeline.git
cd RNAseq_analysis_pipeline

# Build Docker image
docker build -t rnaseq-pipeline -f docker/Dockerfile .

# Run container
docker run -it --rm -v $(pwd)/data:/app/data rnaseq-pipeline
```

### Option 3: Docker Compose

```bash
# Start all services
docker-compose -f docker/docker-compose.yml up -d

# Access pipeline container
docker-compose exec pipeline bash

# Access Jupyter service
# Open http://localhost:8888 in browser

# Stop services
docker-compose down
```

## Cloud Installation

### AWS EC2

```bash
# Launch EC2 instance (recommended: m5.4xlarge or larger)
# Connect to instance
ssh -i your-key.pem ec2-user@your-instance-ip

# Install git and download pipeline
sudo yum update -y
sudo yum install -y git
git clone https://github.com/username/RNAseq_analysis_pipeline.git

# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Setup pipeline
cd RNAseq_analysis_pipeline
./scripts/setup_environment.sh
```

### Google Cloud Platform

```bash
# Create VM instance
gcloud compute instances create rnaseq-vm \
    --image-family=ubuntu-2004-lts \
    --image-project=ubuntu-os-cloud \
    --machine-type=n1-standard-8 \
    --boot-disk-size=200GB

# Connect and install
gcloud compute ssh rnaseq-vm
# Follow standard Linux installation steps
```

## Reference Data Setup

### Automatic Download

```bash
# Download human reference genome and annotation
cd data/reference
./download_references.sh
```

### Manual Download

```bash
# Create reference directory
mkdir -p data/reference
cd data/reference

# Download human genome (GRCh38)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Download gene annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz
```

### Alternative Reference Genomes

For mouse (GRCm39):
```bash
# Mouse genome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz

# Mouse annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz
```

## Troubleshooting

### Common Issues

#### 1. Conda Installation Fails

**Error**: `conda: command not found`

**Solution**:
```bash
# Ensure conda is in PATH
export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc

# Or reinstall conda
rm -rf ~/miniconda3
# Download and install again
```

#### 2. Bioconda Packages Not Found

**Error**: `PackagesNotFoundError: The following packages are not available`

**Solution**:
```bash
# Update conda
conda update conda

# Add channels in correct order
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### 3. STAR Index Building Fails

**Error**: `EXITING because of FATAL ERROR: not enough memory for BAM sorting`

**Solution**:
```bash
# Reduce memory usage or increase available RAM
# Edit config.yaml:
resources:
  memory_gb: 16  # Reduce from 32
```

#### 4. R Package Installation Issues

**Error**: `there is no package called 'DESeq2'`

**Solution**:
```bash
# Install missing Bioconductor packages
conda install -c bioconda bioconductor-deseq2
conda install -c bioconda bioconductor-org.hs.eg.db
```

#### 5. Permission Denied Errors

**Error**: `Permission denied: '/opt/conda/...'`

**Solution**:
```bash
# Fix conda permissions
sudo chown -R $USER:$USER ~/miniconda3

# Or use user installation
conda create --prefix ./envs/rnaseq python=3.9
```

### Platform-Specific Issues

#### macOS

- **Issue**: Apple Silicon (M1/M2) compatibility
- **Solution**: Use conda with osx-arm64 channel or use Docker

```bash
# For Apple Silicon
conda config --env --set subdir osx-arm64
```

#### CentOS/RHEL

- **Issue**: glibc version compatibility
- **Solution**: Use Docker or update system packages

```bash
# Update system
sudo yum update -y
sudo yum install -y gcc gcc-c++ make
```

## Verification

### Test Installation

```bash
# Activate environment
conda activate rnaseq-pipeline

# Run validation script
./validate_installation.sh

# Test core tools
fastqc --version
STAR --version
snakemake --version

# Test Python packages
python -c "import pandas, numpy, matplotlib; print('Python packages OK')"

# Test R packages
R --slave -e "library(DESeq2); cat('R packages OK\n')"
```

### Test with Example Data

```bash
# Create test FASTQ files
mkdir -p data/test_data
# (Test data creation script will be provided)

# Run pipeline dry-run
./scripts/run_rnaseq_pipeline.sh --dry-run

# Run Snakemake dry-run
cd workflow
snakemake --dry-run
```

## Performance Optimization

### System Tuning

```bash
# Increase file descriptor limits
echo "* soft nofile 65535" | sudo tee -a /etc/security/limits.conf
echo "* hard nofile 65535" | sudo tee -a /etc/security/limits.conf

# Optimize I/O scheduler for SSDs
echo noop | sudo tee /sys/block/sda/queue/scheduler
```

### Memory Management

```bash
# Monitor memory usage
htop
free -h

# Set environment variables for memory-intensive tools
export TMPDIR=/tmp
export STAR_limitGenomeGenerateRAM=32000000000  # 32GB
```

## Getting Help

If you encounter issues not covered here:

1. Check the [FAQ](faq.md)
2. Search existing [GitHub Issues](https://github.com/username/RNAseq_analysis_pipeline/issues)
3. Create a new issue with:
   - System information (`uname -a`)
   - Conda environment (`conda list`)
   - Error messages (full stack trace)
   - Steps to reproduce

## Next Steps

After successful installation:

1. Read the [Usage Guide](usage.md)
2. Review [Configuration Options](configuration.md)
3. Try the [Tutorial Examples](examples/)
4. Join our [Community Discussions](https://github.com/username/RNAseq_analysis_pipeline/discussions)