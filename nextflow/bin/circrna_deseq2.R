#!/usr/bin/env Rscript
# Differential Expression Analysis for circRNAs using DESeq2
# Based on the comprehensive circRNA-seq analysis framework

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# Function to log messages
log_message <- function(msg) {
    cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

main <- function() {
    log_message("Starting circRNA differential expression analysis")

    # Get parameters from Snakemake or set defaults
    if (exists("snakemake")) {
        count_file <- snakemake@input$counts
        sample_file <- snakemake@input$sample_info
        output_dir <- snakemake@params$output_dir
        results_file <- snakemake@output$results
        volcano_file <- snakemake@output$volcano
        pca_file <- snakemake@output$pca
        pval_cutoff <- as.numeric(snakemake@params$pval_cutoff)
        log2fc_cutoff <- as.numeric(snakemake@params$log2fc_cutoff)
        log_file <- snakemake@log[[1]]
    } else {
        # Default values for standalone execution
        count_file <- "../results/circrna/05_ciri_quant/merged_count_matrix.tsv"
        sample_file <- "../results/circrna/06_de_analysis/sample_info.tsv"
        output_dir <- "../results/circrna/06_de_analysis"
        results_file <- file.path(output_dir, "DE_results.csv")
        volcano_file <- file.path(output_dir, "volcano_plot.pdf")
        pca_file <- file.path(output_dir, "pca_plot.pdf")
        pval_cutoff <- 0.05
        log2fc_cutoff <- 1.0
        log_file <- file.path(output_dir, "deseq2.log")
    }

    # Ensure output directory exists
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Redirect output to log file
    if (!is.null(log_file)) {
        sink(log_file, append = FALSE)
        on.exit(sink())
    }

    log_message(paste("Count file:", count_file))
    log_message(paste("Sample file:", sample_file))
    log_message(paste("Output directory:", output_dir))

    # Check if input files exist
    if (!file.exists(count_file)) {
        stop("Count matrix file not found: ", count_file)
    }

    if (!file.exists(sample_file)) {
        stop("Sample info file not found: ", sample_file)
    }

    # Load count data
    log_message("Loading count matrix...")
    count_data <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t")
    log_message(paste("Loaded count matrix with", nrow(count_data), "circRNAs and", ncol(count_data), "samples"))

    # Load sample information
    log_message("Loading sample information...")
    sample_info <- read.table(sample_file, header = TRUE, row.names = 1, sep = "\t")
    log_message(paste("Sample conditions:", paste(unique(sample_info$condition), collapse = ", ")))

    # Check if we have data
    if (nrow(count_data) == 0) {
        log_message("No circRNAs found in count matrix. Creating empty results.")
        empty_results <- data.frame(
            circRNA_ID = character(0),
            baseMean = numeric(0),
            log2FoldChange = numeric(0),
            lfcSE = numeric(0),
            stat = numeric(0),
            pvalue = numeric(0),
            padj = numeric(0)
        )
        write.csv(empty_results, file = results_file, row.names = FALSE)
        return()
    }

    # Ensure sample order matches between count matrix and sample info
    common_samples <- intersect(colnames(count_data), rownames(sample_info))
    log_message(paste("Common samples:", length(common_samples)))

    if (length(common_samples) < 2) {
        stop("Need at least 2 samples with matching names in count matrix and sample info")
    }

    count_data <- count_data[, common_samples, drop = FALSE]
    sample_info <- sample_info[common_samples, , drop = FALSE]

    # Check if we have at least 2 conditions
    if (length(unique(sample_info$condition)) < 2) {
        stop("Need at least 2 different conditions for differential expression analysis")
    }

    log_message("Creating DESeqDataSet...")

    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info,
        design = ~ condition
    )

    log_message("Pre-filtering low count circRNAs...")

    # Pre-filter low count circRNAs
    # Keep circRNAs with at least 10 total counts
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]

    log_message(paste("After filtering:", nrow(dds), "circRNAs retained"))

    if (nrow(dds) < 10) {
        log_message("Warning: Very few circRNAs for analysis")
        if (nrow(dds) == 0) {
            stop("No circRNAs passed filtering")
        }
    }

    # Set reference level (first condition alphabetically)
    dds$condition <- relevel(dds$condition, ref = sort(unique(dds$condition))[1])
    log_message(paste("Reference condition:", levels(dds$condition)[1]))

    log_message("Running DESeq2 analysis...")

    # Run DESeq
    dds <- DESeq(dds)

    # Get results
    res <- results(dds, alpha = pval_cutoff)
    log_message("DESeq2 analysis completed")

    # Convert to data frame and add circRNA IDs
    res_df <- as.data.frame(res) %>%
        rownames_to_column("circRNA_ID") %>%
        arrange(padj, pvalue)

    log_message(paste("Total circRNAs tested:", nrow(res_df)))

    # Count significant results
    sig_results <- res_df %>%
        filter(padj < pval_cutoff & abs(log2FoldChange) > log2fc_cutoff)

    log_message(paste("Significant circRNAs (padj <", pval_cutoff, ", |log2FC| >", log2fc_cutoff, "):", nrow(sig_results)))

    if (nrow(sig_results) > 0) {
        up_regulated <- sum(sig_results$log2FoldChange > 0)
        down_regulated <- sum(sig_results$log2FoldChange < 0)
        log_message(paste("  - Up-regulated:", up_regulated))
        log_message(paste("  - Down-regulated:", down_regulated))
    }

    # Save results
    log_message(paste("Saving results to:", results_file))
    write.csv(res_df, file = results_file, row.names = FALSE)

    # Create visualizations
    log_message("Creating visualizations...")

    # PCA Plot
    log_message("Creating PCA plot...")
    if (nrow(dds) >= 2) {
        tryCatch({
            # Transform data for PCA
            vsd <- vst(dds, blind = FALSE)

            # PCA plot
            pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
            percentVar <- round(100 * attr(pca_data, "percentVar"))

            pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
                geom_point(size = 3) +
                geom_text(aes(label = name), vjust = -1) +
                xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                ggtitle("PCA of circRNA Expression") +
                theme_bw() +
                theme(
                    plot.title = element_text(hjust = 0.5),
                    legend.position = "bottom"
                )

            ggsave(pca_file, pca_plot, width = 8, height = 6)
            log_message(paste("PCA plot saved to:", pca_file))
        }, error = function(e) {
            log_message(paste("Error creating PCA plot:", e$message))
        })
    }

    # Volcano Plot
    log_message("Creating volcano plot...")
    tryCatch({
        # Add significance column
        res_df <- res_df %>%
            mutate(
                significant = case_when(
                    is.na(padj) ~ "Not tested",
                    padj < pval_cutoff & log2FoldChange > log2fc_cutoff ~ "Up-regulated",
                    padj < pval_cutoff & log2FoldChange < -log2fc_cutoff ~ "Down-regulated",
                    TRUE ~ "Not significant"
                )
            )

        # Create volcano plot
        volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
            geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
            scale_color_manual(
                values = c(
                    "Up-regulated" = "#d62728",
                    "Down-regulated" = "#2ca02c",
                    "Not significant" = "#7f7f7f",
                    "Not tested" = "#bcbd22"
                )
            ) +
            geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", alpha = 0.5) +
            geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", alpha = 0.5) +
            labs(
                title = "Volcano Plot of Differential circRNA Expression",
                x = "Log2 Fold Change",
                y = "-Log10 P-value",
                color = "Significance"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.position = "bottom"
            )

        # Add text annotations for most significant circRNAs
        if (nrow(sig_results) > 0 && nrow(sig_results) <= 10) {
            top_circs <- head(sig_results, 5)
            volcano_plot <- volcano_plot +
                geom_text(
                    data = top_circs,
                    aes(label = circRNA_ID),
                    size = 3,
                    vjust = -0.5,
                    hjust = 0.5,
                    check_overlap = TRUE
                )
        }

        ggsave(volcano_file, volcano_plot, width = 10, height = 8)
        log_message(paste("Volcano plot saved to:", volcano_file))
    }, error = function(e) {
        log_message(paste("Error creating volcano plot:", e$message))
    })

    # Create heatmap for top significant circRNAs if we have enough
    if (nrow(sig_results) >= 5 && nrow(sig_results) <= 50) {
        log_message("Creating heatmap for top significant circRNAs...")
        tryCatch({
            # Get top 20 most significant circRNAs
            top_circs <- head(sig_results, min(20, nrow(sig_results)))
            top_circ_ids <- top_circs$circRNA_ID

            # Get normalized counts
            norm_counts <- counts(dds, normalized = TRUE)
            top_counts <- norm_counts[top_circ_ids, , drop = FALSE]

            # Log2 transform (add pseudocount)
            log_counts <- log2(top_counts + 1)

            # Create annotation for samples
            annotation_col <- data.frame(
                Condition = sample_info$condition
            )
            rownames(annotation_col) <- rownames(sample_info)

            # Create heatmap
            heatmap_file <- file.path(output_dir, "heatmap_top_circrnas.pdf")
            pdf(heatmap_file, width = 10, height = 8)
            pheatmap(
                log_counts,
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                annotation_col = annotation_col,
                show_rownames = TRUE,
                show_colnames = TRUE,
                color = colorRampPalette(c("blue", "white", "red"))(100),
                main = "Top Differentially Expressed circRNAs"
            )
            dev.off()
            log_message(paste("Heatmap saved to:", heatmap_file))
        }, error = function(e) {
            log_message(paste("Error creating heatmap:", e$message))
        })
    }

    log_message("Differential expression analysis completed successfully")
}

# Run the analysis
if (!interactive()) {
    main()
}