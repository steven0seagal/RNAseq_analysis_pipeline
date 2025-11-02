#!/usr/bin/env Rscript

# ==============================================================================
# Automated RNA-Seq Downstream Analysis and Visualization in R
#
# Description: This script performs differential expression analysis using DESeq2,
#              generates publication-quality visualizations (PCA, Volcano, Heatmap),
#              and conducts functional enrichment analysis (GO & KEGG) using
#              clusterProfiler.
#
# Usage: Rscript interpret_results.R <counts_file> <metadata_file> <output_dir>
#
# Example: Rscript interpret_results.R raw_counts.tsv metadata.tsv deseq2_results
#
# Requirements: DESeq2, pheatmap, RColorBrewer, ggplot2, ggrepel,
#               clusterProfiler, org.Hs.eg.db, AnnotationDbi
# ==============================================================================

# === LOAD REQUIRED LIBRARIES ===
suppressPackageStartupMessages({
    library(DESeq2)
    library(pheatmap)
    library(RColorBrewer)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tibble)
    library(readr)

    # Optional libraries for enrichment analysis
    tryCatch({
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(AnnotationDbi)
        ENRICHMENT_AVAILABLE <- TRUE
    }, error = function(e) {
        cat("Note: Enrichment analysis packages not available. Skipping enrichment analysis.\n")
        ENRICHMENT_AVAILABLE <<- FALSE
    })
})

# === COMMAND LINE ARGUMENT PARSING ===
args <- commandArgs(trailingOnly = TRUE)

# Default parameters if no arguments provided (for testing)
if (length(args) == 0) {
    cat("No arguments provided. Using default test parameters.\n")
    args <- c("results/05_featurecounts/raw_counts.tsv",
              "config/metadata.tsv",
              "results/analysis_R")
}

if (length(args) != 3) {
    cat("ERROR: Incorrect number of arguments provided.\n")
    cat("Usage: Rscript interpret_results.R <counts_file> <metadata_file> <output_dir>\n")
    cat("Example: Rscript interpret_results.R raw_counts.tsv metadata.tsv deseq2_results\n")
    quit(status = 1)
}

counts_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]

# === ANALYSIS PARAMETERS ===
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1.0
MIN_BASEMEAN <- 10

# === CREATE OUTPUT DIRECTORY ===
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === LOGGING FUNCTION ===
log_message <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", timestamp, message))
}

log_message("Starting RNA-Seq differential expression analysis")
log_message(sprintf("Counts file: %s", counts_file))
log_message(sprintf("Metadata file: %s", metadata_file))
log_message(sprintf("Output directory: %s", output_dir))

# ==============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ==============================================================================
log_message("Loading and preprocessing data")

# Validate input files exist
if (!file.exists(counts_file)) {
    stop(sprintf("Counts file not found: %s", counts_file))
}

if (!file.exists(metadata_file)) {
    stop(sprintf("Metadata file not found: %s", metadata_file))
}

# Load count data from featureCounts
log_message("Loading count data")
count_data <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Extract only the count columns (assuming featureCounts format)
# featureCounts output has: Geneid, Chr, Start, End, Strand, Length, then sample columns
if (ncol(count_data) < 7) {
    stop("Invalid count file format. Expected featureCounts output format.")
}

# Extract sample count columns (skip first 6 metadata columns)
count_matrix <- count_data[, 7:ncol(count_data)]

# Clean up column names (remove path and .bam extension)
colnames(count_matrix) <- gsub(".*\\/", "", colnames(count_matrix))  # Remove path
colnames(count_matrix) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(count_matrix))  # Remove STAR suffix
colnames(count_matrix) <- gsub("\\.bam$", "", colnames(count_matrix))  # Remove .bam extension

log_message(sprintf("Loaded count matrix with %d genes and %d samples", nrow(count_matrix), ncol(count_matrix)))
log_message(sprintf("Sample names: %s", paste(colnames(count_matrix), collapse = ", ")))

# Load metadata
log_message("Loading sample metadata")
metadata <- read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

log_message(sprintf("Loaded metadata for %d samples", nrow(metadata)))

# Ensure sample order matches between counts and metadata
common_samples <- intersect(colnames(count_matrix), rownames(metadata))
if (length(common_samples) == 0) {
    stop("No matching samples found between count matrix and metadata")
}

if (length(common_samples) < ncol(count_matrix)) {
    missing_samples <- setdiff(colnames(count_matrix), rownames(metadata))
    log_message(sprintf("Warning: %d samples in count matrix not found in metadata: %s",
                       length(missing_samples), paste(missing_samples, collapse = ", ")))
}

# Subset to common samples
count_matrix <- count_matrix[, common_samples]
metadata <- metadata[common_samples, ]

log_message(sprintf("Proceeding with %d common samples", length(common_samples)))

# Ensure count matrix contains only integers
count_matrix <- round(count_matrix)

# ==============================================================================
# 2. DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================
log_message("Running DESeq2 analysis")

# Check if condition column exists
if (!"condition" %in% colnames(metadata)) {
    stop("'condition' column not found in metadata file")
}

# Set condition as factor with appropriate reference level
unique_conditions <- unique(metadata$condition)
log_message(sprintf("Found conditions: %s", paste(unique_conditions, collapse = ", ")))

if (length(unique_conditions) != 2) {
    stop(sprintf("Expected exactly 2 conditions, found %d: %s",
                length(unique_conditions), paste(unique_conditions, collapse = ", ")))
}

# Set reference level (Healthy if available, otherwise first alphabetically)
if ("Healthy" %in% unique_conditions) {
    reference_level <- "Healthy"
    test_level <- setdiff(unique_conditions, "Healthy")
} else {
    reference_level <- sort(unique_conditions)[1]
    test_level <- sort(unique_conditions)[2]
}

metadata$condition <- factor(metadata$condition, levels = c(reference_level, test_level))
log_message(sprintf("Comparison: %s vs %s (reference)", test_level, reference_level))

# Create DESeqDataSet object
design_formula <- ~ condition

# Add batch effect if batch column exists
if ("batch" %in% colnames(metadata) && length(unique(metadata$batch)) > 1) {
    log_message("Batch column detected. Including batch effect in design.")
    design_formula <- ~ batch + condition
}

log_message(sprintf("Design formula: %s", deparse(design_formula)))

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = design_formula)

# Pre-filter genes with very low counts
keep <- rowSums(counts(dds)) >= MIN_BASEMEAN
dds <- dds[keep,]

log_message(sprintf("Filtered to %d genes with sufficient expression", nrow(dds)))

# Run DESeq2
log_message("Running DESeq2 differential expression analysis")
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", test_level, reference_level))
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)  # Remove rows with NA values

log_message(sprintf("DESeq2 analysis complete. %d genes tested.", nrow(res_df)))

# Save full results table
full_results_file <- file.path(output_dir, "deseq2_full_results.csv")
write.csv(res_df, full_results_file)
log_message(sprintf("Full results saved: %s", full_results_file))

# ==============================================================================
# 3. IDENTIFY AND SAVE SIGNIFICANT GENE LISTS
# ==============================================================================
log_message("Identifying significant genes")

# Apply significance thresholds
sig_genes <- subset(res_df, padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD)

# Separate into up- and down-regulated
up_regulated <- subset(sig_genes, log2FoldChange > 0)
down_regulated <- subset(sig_genes, log2FoldChange < 0)

log_message(sprintf("Found %d significantly upregulated genes", nrow(up_regulated)))
log_message(sprintf("Found %d significantly downregulated genes", nrow(down_regulated)))

# Add gene symbols if possible
if (ENRICHMENT_AVAILABLE) {
    tryCatch({
        # Extract Ensembl IDs (remove version numbers)
        get_ensembl_id <- function(ids) {
            sapply(strsplit(ids, "\\."), `[`, 1)
        }

        if (nrow(up_regulated) > 0) {
            up_ensembl_ids <- get_ensembl_id(rownames(up_regulated))
            up_regulated$symbol <- mapIds(org.Hs.eg.db,
                                        keys = up_ensembl_ids,
                                        column = "SYMBOL",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")
        }

        if (nrow(down_regulated) > 0) {
            down_ensembl_ids <- get_ensembl_id(rownames(down_regulated))
            down_regulated$symbol <- mapIds(org.Hs.eg.db,
                                          keys = down_ensembl_ids,
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")
        }
    }, error = function(e) {
        log_message("Warning: Could not add gene symbols")
    })
}

# Save gene lists
if (nrow(up_regulated) > 0) {
    up_file <- file.path(output_dir, "up_regulated_genes.csv")
    write.csv(up_regulated, up_file)
    log_message(sprintf("Upregulated genes saved: %s", up_file))
}

if (nrow(down_regulated) > 0) {
    down_file <- file.path(output_dir, "down_regulated_genes.csv")
    write.csv(down_regulated, down_file)
    log_message(sprintf("Downregulated genes saved: %s", down_file))
}

# ==============================================================================
# 4. GENERATE VISUALIZATIONS
# ==============================================================================
log_message("Generating visualizations")

# === PCA PLOT ===
log_message("Creating PCA plot")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = name), size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    ggtitle("Principal Component Analysis") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "bottom"
    ) +
    scale_color_brewer(type = "qual", palette = "Set1")

pca_file <- file.path(output_dir, "pca_plot.png")
ggsave(pca_file, plot = pca_plot, width = 8, height = 6, dpi = 300)
log_message(sprintf("PCA plot saved: %s", pca_file))

# === VOLCANO PLOT ===
log_message("Creating volcano plot")

# Prepare data for volcano plot
volcano_data <- res_df
volcano_data$diffexpressed <- "NO"
volcano_data$diffexpressed[volcano_data$log2FoldChange > LOG2FC_THRESHOLD & volcano_data$padj < PADJ_THRESHOLD] <- "UP"
volcano_data$diffexpressed[volcano_data$log2FoldChange < -LOG2FC_THRESHOLD & volcano_data$padj < PADJ_THRESHOLD] <- "DOWN"

# Add gene labels for top significant genes
volcano_data$delabel <- NA
top_genes <- head(volcano_data[volcano_data$diffexpressed != "NO", ], 20)
if (nrow(top_genes) > 0 && ENRICHMENT_AVAILABLE) {
    tryCatch({
        top_ensembl_ids <- get_ensembl_id(rownames(top_genes))
        volcano_data[rownames(top_genes), "delabel"] <- mapIds(org.Hs.eg.db,
                                                              keys = top_ensembl_ids,
                                                              column = "SYMBOL",
                                                              keytype = "ENSEMBL",
                                                              multiVals = "first")
    }, error = function(e) {
        volcano_data[rownames(top_genes), "delabel"] <- rownames(top_genes)
    })
}

volcano_plot <- ggplot(data = volcano_data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    geom_text_repel(max.overlaps = 10, size = 3) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
                      name = "Expression",
                      labels = c("DOWN" = "Downregulated", "NO" = "Not significant", "UP" = "Upregulated")) +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), col = "black", linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD), col = "black", linetype = "dashed", alpha = 0.7) +
    ggtitle(sprintf("Volcano Plot: %s vs %s", test_level, reference_level)) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "bottom"
    )

volcano_file <- file.path(output_dir, "volcano_plot.png")
ggsave(volcano_file, plot = volcano_plot, width = 10, height = 8, dpi = 300)
log_message(sprintf("Volcano plot saved: %s", volcano_file))

# === HEATMAP OF TOP DEGS ===
if (nrow(sig_genes) > 1) {
    log_message("Creating heatmap of top differentially expressed genes")

    # Select top 50 most significant genes
    n_genes <- min(50, nrow(sig_genes))
    top_sig_genes <- head(sig_genes[order(sig_genes$padj), ], n_genes)

    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    heatmap_matrix <- norm_counts[rownames(top_sig_genes), ]

    # Scale by rows (genes)
    heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

    # Prepare annotation
    col_annotation <- data.frame(
        Condition = metadata[colnames(heatmap_matrix), "condition"]
    )
    rownames(col_annotation) <- colnames(heatmap_matrix)

    # Color palette
    annotation_colors <- list(
        Condition = setNames(brewer.pal(length(unique_conditions), "Set1")[1:length(unique_conditions)],
                           unique_conditions)
    )

    heatmap_file <- file.path(output_dir, "heatmap_top50_degs.png")
    png(heatmap_file, width = 10, height = 12, units = "in", res = 300)

    pheatmap(heatmap_matrix_scaled,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             annotation_col = col_annotation,
             annotation_colors = annotation_colors,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = sprintf("Top %d Differentially Expressed Genes", n_genes),
             fontsize = 10)

    dev.off()
    log_message(sprintf("Heatmap saved: %s", heatmap_file))
} else {
    log_message("No significant genes found for heatmap")
}

# ==============================================================================
# 5. FUNCTIONAL ENRICHMENT ANALYSIS
# ==============================================================================
if (ENRICHMENT_AVAILABLE && (nrow(up_regulated) > 0 || nrow(down_regulated) > 0)) {
    log_message("Performing functional enrichment analysis")

    # Function to get clean Ensembl IDs
    get_ensembl_id <- function(id_vector) {
        sapply(strsplit(id_vector, "\\."), `[`, 1)
    }

    # GO enrichment for upregulated genes
    if (nrow(up_regulated) > 0) {
        log_message("GO enrichment analysis for upregulated genes")
        tryCatch({
            up_ensembl <- get_ensembl_id(rownames(up_regulated))

            ego_up <- enrichGO(gene = up_ensembl,
                              OrgDb = org.Hs.eg.db,
                              keyType = "ENSEMBL",
                              ont = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

            if (nrow(ego_up) > 0) {
                go_up_file <- file.path(output_dir, "GO_enrichment_upregulated.csv")
                write.csv(as.data.frame(ego_up), go_up_file)
                log_message(sprintf("GO enrichment (up) saved: %s", go_up_file))

                # Create dotplot
                go_up_plot_file <- file.path(output_dir, "GO_enrichment_upregulated_plot.png")
                png(go_up_plot_file, width = 12, height = 8, units = "in", res = 300)
                print(dotplot(ego_up, showCategory = 20) + ggtitle("GO Enrichment - Upregulated Genes"))
                dev.off()
            }
        }, error = function(e) {
            log_message(sprintf("GO enrichment for upregulated genes failed: %s", e$message))
        })
    }

    # GO enrichment for downregulated genes
    if (nrow(down_regulated) > 0) {
        log_message("GO enrichment analysis for downregulated genes")
        tryCatch({
            down_ensembl <- get_ensembl_id(rownames(down_regulated))

            ego_down <- enrichGO(gene = down_ensembl,
                                OrgDb = org.Hs.eg.db,
                                keyType = "ENSEMBL",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)

            if (nrow(ego_down) > 0) {
                go_down_file <- file.path(output_dir, "GO_enrichment_downregulated.csv")
                write.csv(as.data.frame(ego_down), go_down_file)
                log_message(sprintf("GO enrichment (down) saved: %s", go_down_file))

                # Create dotplot
                go_down_plot_file <- file.path(output_dir, "GO_enrichment_downregulated_plot.png")
                png(go_down_plot_file, width = 12, height = 8, units = "in", res = 300)
                print(dotplot(ego_down, showCategory = 20) + ggtitle("GO Enrichment - Downregulated Genes"))
                dev.off()
            }
        }, error = function(e) {
            log_message(sprintf("GO enrichment for downregulated genes failed: %s", e$message))
        })
    }

    # KEGG pathway enrichment
    if (nrow(up_regulated) > 0) {
        log_message("KEGG pathway analysis for upregulated genes")
        tryCatch({
            up_entrez <- mapIds(org.Hs.eg.db,
                               keys = get_ensembl_id(rownames(up_regulated)),
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")
            up_entrez <- up_entrez[!is.na(up_entrez)]

            if (length(up_entrez) > 0) {
                kegg_up <- enrichKEGG(gene = up_entrez,
                                     organism = "hsa",
                                     pvalueCutoff = 0.05)

                if (nrow(kegg_up) > 0) {
                    kegg_up_file <- file.path(output_dir, "KEGG_enrichment_upregulated.csv")
                    write.csv(as.data.frame(kegg_up), kegg_up_file)
                    log_message(sprintf("KEGG enrichment (up) saved: %s", kegg_up_file))
                }
            }
        }, error = function(e) {
            log_message(sprintf("KEGG enrichment for upregulated genes failed: %s", e$message))
        })
    }
} else {
    log_message("Skipping enrichment analysis (packages not available or no significant genes)")
}

# ==============================================================================
# 6. GENERATE SUMMARY REPORT
# ==============================================================================
log_message("Generating analysis summary report")

summary_file <- file.path(output_dir, "analysis_summary.txt")
cat("RNA-Seq Differential Expression Analysis Summary\n", file = summary_file)
cat("================================================\n\n", file = summary_file, append = TRUE)

cat(sprintf("Analysis Date: %s\n", Sys.Date()), file = summary_file, append = TRUE)
cat(sprintf("Comparison: %s vs %s (reference)\n\n", test_level, reference_level), file = summary_file, append = TRUE)

cat("Input Data:\n", file = summary_file, append = TRUE)
cat(sprintf("- Count file: %s\n", counts_file), file = summary_file, append = TRUE)
cat(sprintf("- Metadata file: %s\n", metadata_file), file = summary_file, append = TRUE)
cat(sprintf("- Total genes: %d\n", nrow(count_matrix)), file = summary_file, append = TRUE)
cat(sprintf("- Total samples: %d\n", ncol(count_matrix)), file = summary_file, append = TRUE)
cat(sprintf("- Genes after filtering: %d\n\n", nrow(dds)), file = summary_file, append = TRUE)

cat("Analysis Parameters:\n", file = summary_file, append = TRUE)
cat(sprintf("- Adjusted p-value threshold: %g\n", PADJ_THRESHOLD), file = summary_file, append = TRUE)
cat(sprintf("- Log2 fold change threshold: %g\n", LOG2FC_THRESHOLD), file = summary_file, append = TRUE)
cat(sprintf("- Minimum base mean: %g\n\n", MIN_BASEMEAN), file = summary_file, append = TRUE)

cat("Results:\n", file = summary_file, append = TRUE)
cat(sprintf("- Total significant genes: %d\n", nrow(sig_genes)), file = summary_file, append = TRUE)
cat(sprintf("- Upregulated genes: %d\n", nrow(up_regulated)), file = summary_file, append = TRUE)
cat(sprintf("- Downregulated genes: %d\n\n", nrow(down_regulated)), file = summary_file, append = TRUE)

cat("Output Files:\n", file = summary_file, append = TRUE)
cat("- deseq2_full_results.csv: Complete DESeq2 results\n", file = summary_file, append = TRUE)
cat("- up_regulated_genes.csv: Significantly upregulated genes\n", file = summary_file, append = TRUE)
cat("- down_regulated_genes.csv: Significantly downregulated genes\n", file = summary_file, append = TRUE)
cat("- pca_plot.png: Principal component analysis\n", file = summary_file, append = TRUE)
cat("- volcano_plot.png: Volcano plot of differential expression\n", file = summary_file, append = TRUE)
cat("- heatmap_top50_degs.png: Heatmap of top differentially expressed genes\n", file = summary_file, append = TRUE)

if (ENRICHMENT_AVAILABLE) {
    cat("- GO_enrichment_*.csv: Gene Ontology enrichment results\n", file = summary_file, append = TRUE)
    cat("- KEGG_enrichment_*.csv: KEGG pathway enrichment results\n", file = summary_file, append = TRUE)
}

log_message(sprintf("Analysis complete! Results saved in: %s", output_dir))
log_message(sprintf("Summary report: %s", summary_file))

# Print session info for reproducibility
session_info_file <- file.path(output_dir, "session_info.txt")
capture.output(sessionInfo(), file = session_info_file)

log_message("R script execution completed successfully")