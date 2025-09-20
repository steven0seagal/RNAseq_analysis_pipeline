# ==============================================================================
# RNA-Seq Analysis Pipeline Docker Container
#
# This Dockerfile creates a containerized environment for RNA-seq analysis
# with all required bioinformatics tools and dependencies
#
# Build: docker build -t rnaseq-pipeline .
# Run: docker run -v /path/to/data:/data rnaseq-pipeline
# ==============================================================================

# Use conda-forge base image with miniconda
FROM condaforge/mambaforge:latest

# Metadata
LABEL maintainer="RNA-Seq Pipeline Team"
LABEL description="Complete RNA-seq differential expression analysis pipeline"
LABEL version="1.0.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH="/opt/conda/envs/rnaseq-pipeline/bin:$PATH"

# Update system packages
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    tree \
    htop \
    vim \
    less \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /app

# Copy environment file
COPY environment.yml /app/

# Create conda environment
RUN mamba env create -f environment.yml && \
    mamba clean -afy

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "rnaseq-pipeline", "/bin/bash", "-c"]

# Copy pipeline files
COPY scripts/ /app/scripts/
COPY workflow/ /app/workflow/
COPY config/ /app/config/
COPY README.md /app/

# Make scripts executable
RUN chmod +x /app/scripts/*.sh && \
    chmod +x /app/scripts/*.R && \
    chmod +x /app/scripts/*.py

# Create data directories
RUN mkdir -p /app/data/raw_fastq \
    /app/data/reference \
    /app/data/test_data \
    /app/results

# Create entrypoint script
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /app/entrypoint.sh && \
    echo 'conda activate rnaseq-pipeline' >> /app/entrypoint.sh && \
    echo 'exec "$@"' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/app/entrypoint.sh"]

# Default command
CMD ["bash"]

# Expose ports for Jupyter (if needed)
EXPOSE 8888

# Set volumes for data
VOLUME ["/app/data", "/app/results"]

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD conda list -n rnaseq-pipeline | grep -q "fastqc" || exit 1