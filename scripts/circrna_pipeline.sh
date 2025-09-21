#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error.
set -o pipefail # The return value of a pipeline is the status of the last command to exit with a non-zero status.

#================================================================================
# COMPREHENSIVE circRNA-SEQ ANALYSIS PIPELINE (BASH SCRIPT)
#================================================================================
# This script performs a full circRNA analysis workflow, from raw FASTQ files
# to functional annotation of differentially expressed circRNAs.
#
# Based on the comprehensive circRNA-seq analysis framework detailed in research3
#
# USAGE:
# 1. Activate the conda environment: conda activate rnaseq-pipeline
# 2. Edit the CONFIGURATION section below.
# 3. Make the script executable: chmod +x circrna_pipeline.sh
# 4. Run the script: ./circrna_pipeline.sh
#================================================================================

#===========================#
# --- CONFIGURATION --- #
#===========================#

# --- Paths to Input Data ---
# Directory containing paired-end FASTQ files (e.g., sample1_R1.fastq.gz, sample1_R2.fastq.gz)
FASTQ_DIR="${1:-data/raw_fastq}"
# A comma-separated list of sample names (prefixes of FASTQ files)
SAMPLES="${2:-Cancer_1,Cancer_2,Cancer_3,Healthy_1,Healthy_2,Healthy_3}"
# A comma-separated list of conditions corresponding to the samples
CONDITIONS="${3:-Cancer,Cancer,Cancer,Healthy,Healthy,Healthy}"

# --- Paths to Reference Files ---
REF_GENOME="${REF_GENOME:-data/reference/GRCh38.primary_assembly.genome.fa}"
BWA_INDEX="${BWA_INDEX:-data/reference/GRCh38.primary_assembly.genome.fa}"
GTF_ANNOTATION="${GTF_ANNOTATION:-data/reference/gencode.v38.primary_assembly.annotation.gtf}"
MIRNA_FASTA="${MIRNA_FASTA:-data/reference/mature_mirna.fa}"

# --- Analysis Parameters ---
THREADS="${THREADS:-8}"
ADAPTER_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
MIN_READ_LEN=35
P_VALUE_CUTOFF=0.05
LOG2FC_CUTOFF=1.0

# --- Output Directories ---
OUTPUT_DIR="results/circrna_analysis"

#================================================================================
# --- SCRIPT EXECUTION ---
#================================================================================

echo "--- Starting Comprehensive circRNA-seq Analysis Pipeline ---"
echo "Input directory: ${FASTQ_DIR}"
echo "Samples: ${SAMPLES}"
echo "Conditions: ${CONDITIONS}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# --- Create Directory Structure ---
mkdir -p ${OUTPUT_DIR}/{01_fastqc_raw,02_trimmed_reads,03_fastqc_trimmed,04_alignment,05_ciri_quant,06_de_analysis,07_functional_annotation}

# --- Convert comma-separated strings to arrays ---
IFS=',' read -r -a SAMPLE_ARRAY <<< "$SAMPLES"
IFS=',' read -r -a CONDITION_ARRAY <<< "$CONDITIONS"

echo "Processing ${#SAMPLE_ARRAY[@]} samples..."

#================================================#
# STAGE 1 & 2: QC AND TRIMMING                   #
#================================================#
echo "--- STAGE 1 & 2: Running FastQC and Cutadapt ---"
for sample in "${SAMPLE_ARRAY[@]}"; do
    R1_IN="${FASTQ_DIR}/${sample}_R1.fastq.gz"
    R2_IN="${FASTQ_DIR}/${sample}_R2.fastq.gz"
    R1_OUT="${OUTPUT_DIR}/02_trimmed_reads/${sample}_R1.trimmed.fastq.gz"
    R2_OUT="${OUTPUT_DIR}/02_trimmed_reads/${sample}_R2.trimmed.fastq.gz"

    echo "Processing sample: ${sample}"

    # Check if input files exist
    if [[ ! -f "${R1_IN}" ]] || [[ ! -f "${R2_IN}" ]]; then
        echo "Warning: Input files for ${sample} not found. Skipping..."
        continue
    fi

    # Raw FastQC
    if [[ ! -f "${OUTPUT_DIR}/01_fastqc_raw/${sample}_R1_fastqc.html" ]]; then
        echo "  Running FastQC on raw reads..."
        fastqc -t ${THREADS} -o ${OUTPUT_DIR}/01_fastqc_raw ${R1_IN} ${R2_IN}
    fi

    # Adapter and Quality Trimming
    if [[ ! -f "${R1_OUT}" ]]; then
        echo "  Running Cutadapt for trimming..."
        cutadapt -j ${THREADS} \
            -a ${ADAPTER_FWD} -A ${ADAPTER_REV} \
            -q 20,20 -m ${MIN_READ_LEN} \
            -o ${R1_OUT} -p ${R2_OUT} \
            ${R1_IN} ${R2_IN} > ${OUTPUT_DIR}/02_trimmed_reads/${sample}_trimming_report.txt
    fi

    # Trimmed FastQC
    if [[ ! -f "${OUTPUT_DIR}/03_fastqc_trimmed/${sample}_R1.trimmed_fastqc.html" ]]; then
        echo "  Running FastQC on trimmed reads..."
        fastqc -t ${THREADS} -o ${OUTPUT_DIR}/03_fastqc_trimmed ${R1_OUT} ${R2_OUT}
    fi
done

#================================================#
# STAGE 3: ALIGNMENT WITH BWA-MEM                #
#================================================#
echo "--- STAGE 3: Aligning reads with BWA-MEM ---"

# First, check if BWA index exists. If not, create it.
if [[ ! -f "${BWA_INDEX}.bwt" ]]; then
    echo "BWA index not found. Creating index for ${REF_GENOME}..."
    bwa index ${REF_GENOME}
fi

for sample in "${SAMPLE_ARRAY[@]}"; do
    R1_TRIMMED="${OUTPUT_DIR}/02_trimmed_reads/${sample}_R1.trimmed.fastq.gz"
    R2_TRIMMED="${OUTPUT_DIR}/02_trimmed_reads/${sample}_R2.trimmed.fastq.gz"
    SAM_OUT="${OUTPUT_DIR}/04_alignment/${sample}.sam"
    BAM_OUT="${OUTPUT_DIR}/04_alignment/${sample}.bam"

    # Skip if input files don't exist
    if [[ ! -f "${R1_TRIMMED}" ]] || [[ ! -f "${R2_TRIMMED}" ]]; then
        echo "Warning: Trimmed files for ${sample} not found. Skipping alignment..."
        continue
    fi

    if [[ ! -f "${SAM_OUT}" ]]; then
        echo "Aligning sample: ${sample}"
        # BWA-MEM alignment. CIRI2 requires the SAM file.
        bwa mem -t ${THREADS} -T 19 ${BWA_INDEX} ${R1_TRIMMED} ${R2_TRIMMED} > ${SAM_OUT}
    fi

    # Convert to BAM for other potential uses (optional but good practice)
    if [[ ! -f "${BAM_OUT}" ]]; then
        samtools view -@ ${THREADS} -Sb ${SAM_OUT} > ${BAM_OUT}
    fi
done

#================================================#
# STAGE 4: circRNA IDENTIFICATION & QUANTIFICATION #
#================================================#
echo "--- STAGE 4: Running CIRIquant for circRNA identification ---"

# Create a config file for CIRIquant
CIRI_CONFIG="${OUTPUT_DIR}/05_ciri_quant/ciri_config.yml"
mkdir -p ${OUTPUT_DIR}/05_ciri_quant

echo "tools:" > ${CIRI_CONFIG}
echo "  bwa: $(which bwa)" >> ${CIRI_CONFIG}
echo "  hisat2: $(which hisat2 || echo 'hisat2')" >> ${CIRI_CONFIG}
echo "  stringtie: $(which stringtie || echo 'stringtie')" >> ${CIRI_CONFIG}
echo "  samtools: $(which samtools)" >> ${CIRI_CONFIG}
echo "reference:" >> ${CIRI_CONFIG}
echo "  fasta: ${REF_GENOME}" >> ${CIRI_CONFIG}
echo "  gtf: ${GTF_ANNOTATION}" >> ${CIRI_CONFIG}
echo "  bwa_index: ${BWA_INDEX}" >> ${CIRI_CONFIG}

for sample in "${SAMPLE_ARRAY[@]}"; do
    R1_TRIMMED="${OUTPUT_DIR}/02_trimmed_reads/${sample}_R1.trimmed.fastq.gz"
    R2_TRIMMED="${OUTPUT_DIR}/02_trimmed_reads/${sample}_R2.trimmed.fastq.gz"

    # Skip if input files don't exist
    if [[ ! -f "${R1_TRIMMED}" ]] || [[ ! -f "${R2_TRIMMED}" ]]; then
        echo "Warning: Trimmed files for ${sample} not found. Skipping CIRIquant..."
        continue
    fi

    if [[ ! -d "${OUTPUT_DIR}/05_ciri_quant/${sample}" ]]; then
        echo "Running CIRIquant for sample: ${sample}"

        # Check if CIRIquant is available, if not use CIRI2 + CIRIquant workflow
        if command -v CIRIquant &> /dev/null; then
            CIRIquant -t ${THREADS} \
                -1 ${R1_TRIMMED} -2 ${R2_TRIMMED} \
                --config ${CIRI_CONFIG} \
                -p ${sample} \
                -o ${OUTPUT_DIR}/05_ciri_quant/${sample}
        else
            echo "CIRIquant not found. Using CIRI2 for circRNA detection..."
            # Alternative: Use CIRI2 directly
            mkdir -p ${OUTPUT_DIR}/05_ciri_quant/${sample}
            SAM_FILE="${OUTPUT_DIR}/04_alignment/${sample}.sam"
            if [[ -f "${SAM_FILE}" ]]; then
                CIRI2.pl -I ${SAM_FILE} -O ${OUTPUT_DIR}/05_ciri_quant/${sample}/${sample}_ciri.out \
                    -F ${REF_GENOME} -A ${GTF_ANNOTATION} -T ${THREADS} -G ${OUTPUT_DIR}/05_ciri_quant/${sample}/${sample}.log
            fi
        fi
    fi
done

#================================================#
# STAGE 5: MERGE COUNTS & DIFFERENTIAL EXPRESSION #
#================================================#
echo "--- STAGE 5: Merging counts and running DESeq2 ---"

# Create a sample information file for the R script
SAMPLE_INFO_FILE="${OUTPUT_DIR}/06_de_analysis/sample_info.tsv"
echo -e "sample\tcondition" > ${SAMPLE_INFO_FILE}
for i in "${!SAMPLE_ARRAY[@]}"; do
    echo -e "${SAMPLE_ARRAY[$i]}\t${CONDITION_ARRAY[$i]}" >> ${SAMPLE_INFO_FILE}
done

# Create a Python script to merge CIRIquant outputs into a count matrix
MERGE_SCRIPT="${OUTPUT_DIR}/06_de_analysis/merge_counts.py"
COUNT_MATRIX="${OUTPUT_DIR}/06_de_analysis/circRNA_count_matrix.tsv"

cat << 'EOF' > ${MERGE_SCRIPT}
import pandas as pd
import glob
import os
import sys

# Get the output directory from command line argument
output_dir = sys.argv[1] if len(sys.argv) > 1 else "results/circrna_analysis"

# Look for CIRIquant GTF files or CIRI2 output files
ciri_files = glob.glob(f"{output_dir}/05_ciri_quant/*/*.gtf")
if not ciri_files:
    ciri_files = glob.glob(f"{output_dir}/05_ciri_quant/*/*_ciri.out")

print(f"Found {len(ciri_files)} circRNA result files")

count_dict = {}

for f in ciri_files:
    sample_name = os.path.basename(os.path.dirname(f))
    print(f"Processing {sample_name}: {f}")

    try:
        if f.endswith('.gtf'):
            # Process CIRIquant GTF output
            df = pd.read_csv(f, sep='\t', header=None, comment='#')
            if len(df.columns) >= 9:
                df.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
                df_circ = df[df['type'] == 'circRNA'].copy()

                if not df_circ.empty:
                    df_circ['circ_id'] = df_circ['attributes'].str.extract(r'circ_id "([^"]+)"')
                    df_circ['bsj_reads'] = df_circ['attributes'].str.extract(r'bsj "([^"]+)"').astype(int)

                    for _, row in df_circ.iterrows():
                        if pd.notna(row['circ_id']):
                            if row['circ_id'] not in count_dict:
                                count_dict[row['circ_id']] = {}
                            count_dict[row['circ_id']][sample_name] = row['bsj_reads']
        else:
            # Process CIRI2 output
            df = pd.read_csv(f, sep='\t')
            if 'circRNA_ID' in df.columns and '#junction_reads' in df.columns:
                for _, row in df.iterrows():
                    circ_id = row['circRNA_ID']
                    junction_reads = row['#junction_reads']
                    if circ_id not in count_dict:
                        count_dict[circ_id] = {}
                    count_dict[circ_id][sample_name] = junction_reads
    except Exception as e:
        print(f"Error processing {f}: {e}")
        continue

if count_dict:
    count_df = pd.DataFrame(count_dict).T.fillna(0).astype(int)
    count_df.index.name = 'circRNA_ID'
    output_file = f"{output_dir}/06_de_analysis/circRNA_count_matrix.tsv"
    count_df.to_csv(output_file, sep='\t')
    print(f"Count matrix saved to {output_file}")
    print(f"Found {len(count_df)} circRNAs across {len(count_df.columns)} samples")
else:
    print("No circRNA data found!")
    # Create empty matrix
    count_df = pd.DataFrame()
    count_df.to_csv(f"{output_dir}/06_de_analysis/circRNA_count_matrix.tsv", sep='\t')
EOF

python ${MERGE_SCRIPT} ${OUTPUT_DIR}

# Create an R script for DESeq2 analysis
DESEQ_SCRIPT="${OUTPUT_DIR}/06_de_analysis/run_deseq2.R"
DE_RESULTS_FILE="${OUTPUT_DIR}/06_de_analysis/DE_results.csv"
VOLCANO_PLOT="${OUTPUT_DIR}/06_de_analysis/volcano_plot.pdf"

cat << EOF > ${DESEQ_SCRIPT}
suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
})

# Set working directory and load data
output_dir <- "${OUTPUT_DIR}"
count_file <- file.path(output_dir, "06_de_analysis", "circRNA_count_matrix.tsv")
sample_file <- file.path(output_dir, "06_de_analysis", "sample_info.tsv")

if (!file.exists(count_file)) {
    cat("Count matrix file not found:", count_file, "\n")
    quit(status=1)
}

if (!file.exists(sample_file)) {
    cat("Sample info file not found:", sample_file, "\n")
    quit(status=1)
}

# Load data
count_data <- read.table(count_file, header=TRUE, row.names=1, sep="\t")
sample_info <- read.table(sample_file, header=TRUE, row.names=1, sep="\t")

cat("Loaded count matrix with", nrow(count_data), "circRNAs and", ncol(count_data), "samples\n")
cat("Sample conditions:", unique(sample_info\$condition), "\n")

if (nrow(count_data) == 0) {
    cat("No circRNAs found in count matrix. Creating empty results.\n")
    empty_results <- data.frame(
        circRNA_ID = character(0),
        log2FoldChange = numeric(0),
        pvalue = numeric(0),
        padj = numeric(0)
    )
    write.csv(empty_results, file="${DE_RESULTS_FILE}", row.names=FALSE)
    quit(status=0)
}

# Ensure column names match row names
count_data <- count_data[, rownames(sample_info), drop=FALSE]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat("After filtering:", nrow(dds), "circRNAs retained\n")

if (nrow(dds) < 10) {
    cat("Too few circRNAs for differential expression analysis\n")
    quit(status=0)
}

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res) %>%
          rownames_to_column("circRNA_ID") %>%
          arrange(padj)

# Write results to file
write.csv(res_df, file="${DE_RESULTS_FILE}", row.names=FALSE)

# Create a volcano plot
res_df <- res_df %>%
  mutate(significant = ifelse(padj < ${P_VALUE_CUTOFF} & abs(log2FoldChange) > ${LOG2FC_CUTOFF}, "yes", "no"))

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_bw(base_size = 14) +
  labs(title = "Volcano Plot of Differential circRNA Expression",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  geom_vline(xintercept = c(-${LOG2FC_CUTOFF}, ${LOG2FC_CUTOFF}), linetype = "dashed") +
  geom_hline(yintercept = -log10(${P_VALUE_CUTOFF}), linetype = "dashed")

ggsave("${VOLCANO_PLOT}", plot = p, width = 8, height = 6)

cat("DESeq2 analysis completed successfully\n")
cat("Results saved to:", "${DE_RESULTS_FILE}", "\n")
cat("Volcano plot saved to:", "${VOLCANO_PLOT}", "\n")
EOF

# Run R script if R is available
if command -v Rscript &> /dev/null; then
    echo "Running DESeq2 analysis..."
    Rscript ${DESEQ_SCRIPT}
else
    echo "R not found. Skipping differential expression analysis."
fi

#================================================#
# STAGE 6: FUNCTIONAL ANNOTATION                 #
#================================================#
echo "--- STAGE 6: Running functional annotation on significant circRNAs ---"

# Check if DE results exist
if [[ -f "${DE_RESULTS_FILE}" ]]; then
    # Extract significant circRNAs
    SIG_CIRCS_BED="${OUTPUT_DIR}/07_functional_annotation/significant_circs.bed"

    # Parse circRNA IDs and convert to BED format
    awk -F',' -v pval=${P_VALUE_CUTOFF} -v fc=${LOG2FC_CUTOFF} '
    NR > 1 && $5 < pval && ($2 > fc || $2 < -fc) {
        # Parse circRNA ID (format: chr:start-end)
        split($1, parts, ":")
        chr = parts[1]
        split(parts[2], coords, "-")
        start = coords[1]
        end = coords[2]
        print chr "\t" start "\t" end "\t" $1
    }' ${DE_RESULTS_FILE} > ${SIG_CIRCS_BED}

    if [[ -s "${SIG_CIRCS_BED}" ]]; then
        echo "Found $(wc -l < ${SIG_CIRCS_BED}) significant circRNAs"

        # Extract sequences if reference genome exists
        if [[ -f "${REF_GENOME}" ]]; then
            SIG_CIRCS_FASTA="${OUTPUT_DIR}/07_functional_annotation/significant_circs.fasta"
            bedtools getfasta -fi ${REF_GENOME} -bed ${SIG_CIRCS_BED} -fo ${SIG_CIRCS_FASTA}

            # Predict miRNA binding sites with miRanda if available
            if [[ -f "${MIRNA_FASTA}" ]] && command -v miranda &> /dev/null; then
                MIRANDA_OUT="${OUTPUT_DIR}/07_functional_annotation/miranda_predictions.txt"
                miranda ${MIRNA_FASTA} ${SIG_CIRCS_FASTA} -out ${MIRANDA_OUT}
                echo "miRNA binding site predictions completed"
            else
                echo "miRanda or miRNA reference not available. Skipping miRNA prediction."
            fi
        else
            echo "Reference genome not found. Skipping sequence extraction."
        fi
    else
        echo "No significant circRNAs found for functional annotation"
    fi
else
    echo "Differential expression results not found. Skipping functional annotation."
fi

#================================================#
# STAGE 7: SUMMARY REPORT                       #
#================================================#
echo "--- Generating Analysis Summary ---"

SUMMARY_FILE="${OUTPUT_DIR}/analysis_summary.txt"

cat > ${SUMMARY_FILE} << EOF
circRNA-seq Analysis Pipeline Summary
=====================================

Analysis completed on: $(date)
Samples processed: ${SAMPLES}
Output directory: ${OUTPUT_DIR}

Pipeline Stages Completed:
1. Quality Control (FastQC)
2. Adapter Trimming (Cutadapt)
3. Read Alignment (BWA-MEM)
4. circRNA Identification (CIRIquant/CIRI2)
5. Differential Expression (DESeq2)
6. Functional Annotation (miRanda)

Key Output Files:
- Raw count matrix: ${COUNT_MATRIX}
- Differential expression results: ${DE_RESULTS_FILE}
- Volcano plot: ${VOLCANO_PLOT}
- Functional annotations: ${OUTPUT_DIR}/07_functional_annotation/

Sample Information:
EOF

for i in "${!SAMPLE_ARRAY[@]}"; do
    echo "- ${SAMPLE_ARRAY[$i]}: ${CONDITION_ARRAY[$i]}" >> ${SUMMARY_FILE}
done

if [[ -f "${COUNT_MATRIX}" ]]; then
    echo "" >> ${SUMMARY_FILE}
    echo "Statistics:" >> ${SUMMARY_FILE}
    circrnas=$(tail -n +2 ${COUNT_MATRIX} | wc -l)
    samples=$(head -1 ${COUNT_MATRIX} | awk '{print NF-1}')
    echo "- Total circRNAs detected: $circrnas" >> ${SUMMARY_FILE}
    echo "- Total samples: $samples" >> ${SUMMARY_FILE}
fi

echo "" >> ${SUMMARY_FILE}
echo "Next Steps:" >> ${SUMMARY_FILE}
echo "1. Review quality control reports in 01_fastqc_raw/ and 03_fastqc_trimmed/" >> ${SUMMARY_FILE}
echo "2. Examine differential expression results in ${DE_RESULTS_FILE}" >> ${SUMMARY_FILE}
echo "3. Investigate functional annotations for significant circRNAs" >> ${SUMMARY_FILE}
echo "4. Validate interesting candidates experimentally" >> ${SUMMARY_FILE}

echo "--- Pipeline Finished Successfully ---"
echo "Results are located in: ${OUTPUT_DIR}"
echo "Summary report: ${SUMMARY_FILE}"