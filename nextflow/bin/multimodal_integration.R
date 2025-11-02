#!/usr/bin/env Rscript

# ==============================================================================
# Multi-Modal Cancer Genomics Integration Analysis
#
# Comprehensive R script for integrating multiple data types:
# 1. mRNA expression data (DESeq2 results)
# 2. Somatic variant data (Mutect2 VCF files)
# 3. miRNA expression data (miRDeep2 results)
# 4. Functional annotation and pathway analysis
# 5. Multi-omics visualization and reporting
#
# This script creates an integrated view of cancer genomics data
# for comprehensive molecular characterization.
#
# Usage: Rscript multimodal_integration.R [--config config.yaml]
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(VariantAnnotation)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(enrichplot)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(VennDiagram)
    library(RColorBrewer)
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
    library(vcfR)
    library(igraph)
    library(networkD3)
    library(plotly)
    library(DT)
    library(knitr)
    library(rmarkdown)
})

# === CONFIGURATION ===
cat("Multi-Modal Cancer Genomics Integration Analysis\n")
cat("===============================================\n")

# Default parameters
PVALUE_CUTOFF <- 0.05
FC_CUTOFF <- 1.0
MIN_VARIANT_SUPPORT <- 5
OUTPUT_FORMAT <- "png"
DPI <- 300

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_file <- "config/analysis_config.yaml"

if (length(args) > 0) {
    if (args[1] == "--config" && length(args) > 1) {
        config_file <- args[2]
    }
}

# Set working directory to project root
if (file.exists("scripts")) {
    setwd(".")
} else if (file.exists("../scripts")) {
    setwd("..")
} else {
    stop("Cannot find project root directory")
}

# Create output directories
dir.create("09_multimodal", showWarnings = FALSE, recursive = TRUE)
dir.create("09_multimodal/plots", showWarnings = FALSE, recursive = TRUE)
dir.create("09_multimodal/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("09_multimodal/reports", showWarnings = FALSE, recursive = TRUE)

# === UTILITY FUNCTIONS ===

# Logging function
log_message <- function(message) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(timestamp, message, "\n")

    # Write to log file
    log_file <- "09_multimodal/analysis.log"
    write(paste(timestamp, message), file = log_file, append = TRUE)
}

# Save plot function
save_plot <- function(plot_obj, filename, width = 10, height = 8) {
    filepath <- file.path("09_multimodal/plots", paste0(filename, ".", OUTPUT_FORMAT))

    if (OUTPUT_FORMAT == "png") {
        png(filepath, width = width, height = height, units = "in", res = DPI)
        print(plot_obj)
        dev.off()
    } else if (OUTPUT_FORMAT == "pdf") {
        pdf(filepath, width = width, height = height)
        print(plot_obj)
        dev.off()
    }

    log_message(paste("Plot saved:", filepath))
}

# === DATA LOADING FUNCTIONS ===

load_mrna_data <- function() {
    log_message("Loading mRNA expression data...")

    # Load DESeq2 results
    dds_file <- "05_quantification/dds_object.rds"
    results_file <- "05_quantification/deseq2_results.csv"

    if (file.exists(dds_file)) {
        dds <- readRDS(dds_file)
        log_message("Loaded DESeq2 object")
    } else {
        log_message("DESeq2 object not found, attempting to load from CSV")
        dds <- NULL
    }

    if (file.exists(results_file)) {
        mrna_results <- read_csv(results_file, show_col_types = FALSE)
        log_message(paste("Loaded", nrow(mrna_results), "mRNA expression results"))
    } else {
        log_message("mRNA results file not found")
        mrna_results <- NULL
    }

    return(list(dds = dds, results = mrna_results))
}

load_variant_data <- function() {
    log_message("Loading somatic variant data...")

    # Find VCF files
    vcf_files <- list.files("06_variants", pattern = "*_somatic_pass.vcf.gz$", full.names = TRUE)

    if (length(vcf_files) == 0) {
        log_message("No somatic variant VCF files found")
        return(NULL)
    }

    variant_data <- list()

    for (vcf_file in vcf_files) {
        sample_name <- basename(vcf_file) %>% str_remove("_somatic_pass.vcf.gz")

        tryCatch({
            vcf <- read.vcfR(vcf_file)

            # Extract variant information
            variants <- data.frame(
                CHROM = getCHROM(vcf),
                POS = getPOS(vcf),
                REF = getREF(vcf),
                ALT = getALT(vcf),
                QUAL = getQUAL(vcf),
                sample = sample_name,
                stringsAsFactors = FALSE
            )

            # Add variant type
            variants$variant_type <- ifelse(
                nchar(variants$REF) == 1 & nchar(variants$ALT) == 1,
                "SNV", "INDEL"
            )

            variant_data[[sample_name]] <- variants
            log_message(paste("Loaded", nrow(variants), "variants for", sample_name))

        }, error = function(e) {
            log_message(paste("Error loading VCF for", sample_name, ":", e$message))
        })
    }

    if (length(variant_data) > 0) {
        all_variants <- do.call(rbind, variant_data)
        log_message(paste("Total variants loaded:", nrow(all_variants)))
        return(all_variants)
    } else {
        return(NULL)
    }
}

load_mirna_data <- function() {
    log_message("Loading miRNA expression data...")

    # Load miRNA differential expression results
    mirna_file <- "08_small_rna/differential_mirnas.csv"

    if (file.exists(mirna_file)) {
        mirna_results <- read_csv(mirna_file, show_col_types = FALSE)
        log_message(paste("Loaded", nrow(mirna_results), "miRNA expression results"))
        return(mirna_results)
    } else {
        log_message("miRNA results file not found")
        return(NULL)
    }
}

# === INTEGRATION ANALYSIS FUNCTIONS ===

identify_significant_features <- function(mrna_data, variant_data, mirna_data) {
    log_message("Identifying significant molecular features...")

    significant_features <- list()

    # Significant mRNAs
    if (!is.null(mrna_data$results)) {
        sig_mrnas <- mrna_data$results %>%
            filter(!is.na(padj), padj < PVALUE_CUTOFF, abs(log2FoldChange) > FC_CUTOFF) %>%
            mutate(feature_type = "mRNA")

        significant_features$mrnas <- sig_mrnas
        log_message(paste("Significant mRNAs:", nrow(sig_mrnas)))
    }

    # High-impact variants
    if (!is.null(variant_data)) {
        sig_variants <- variant_data %>%
            filter(QUAL > MIN_VARIANT_SUPPORT) %>%
            mutate(feature_type = "variant")

        significant_features$variants <- sig_variants
        log_message(paste("High-quality variants:", nrow(sig_variants)))
    }

    # Significant miRNAs
    if (!is.null(mirna_data)) {
        sig_mirnas <- mirna_data %>%
            filter(!is.na(padj), padj < PVALUE_CUTOFF, abs(log2FoldChange) > FC_CUTOFF) %>%
            mutate(feature_type = "miRNA")

        significant_features$mirnas <- sig_mirnas
        log_message(paste("Significant miRNAs:", nrow(sig_mirnas)))
    }

    return(significant_features)
}

perform_pathway_analysis <- function(significant_features) {
    log_message("Performing pathway enrichment analysis...")

    pathway_results <- list()

    if (!is.null(significant_features$mrnas)) {
        # Get gene symbols
        gene_symbols <- significant_features$mrnas$gene_name
        gene_symbols <- gene_symbols[!is.na(gene_symbols)]

        if (length(gene_symbols) > 10) {
            tryCatch({
                # GO enrichment
                go_results <- enrichGO(
                    gene = gene_symbols,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2
                )

                if (nrow(go_results) > 0) {
                    pathway_results$go <- go_results
                    log_message(paste("GO enrichment terms:", nrow(go_results)))
                }

                # KEGG enrichment
                kegg_results <- enrichKEGG(
                    gene = gene_symbols,
                    organism = "hsa",
                    keyType = "SYMBOL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2
                )

                if (nrow(kegg_results) > 0) {
                    pathway_results$kegg <- kegg_results
                    log_message(paste("KEGG pathways:", nrow(kegg_results)))
                }

            }, error = function(e) {
                log_message(paste("Error in pathway analysis:", e$message))
            })
        }
    }

    return(pathway_results)
}

# === VISUALIZATION FUNCTIONS ===

create_overview_plot <- function(significant_features) {
    log_message("Creating multi-modal overview plot...")

    # Count features by type
    feature_counts <- data.frame(
        type = character(),
        count = numeric(),
        direction = character(),
        stringsAsFactors = FALSE
    )

    if (!is.null(significant_features$mrnas)) {
        mrna_up <- sum(significant_features$mrnas$log2FoldChange > 0)
        mrna_down <- sum(significant_features$mrnas$log2FoldChange < 0)

        feature_counts <- rbind(feature_counts,
                               data.frame(type = "mRNA", count = mrna_up, direction = "Up"),
                               data.frame(type = "mRNA", count = mrna_down, direction = "Down"))
    }

    if (!is.null(significant_features$variants)) {
        var_count <- nrow(significant_features$variants)
        feature_counts <- rbind(feature_counts,
                               data.frame(type = "Variants", count = var_count, direction = "Total"))
    }

    if (!is.null(significant_features$mirnas)) {
        mirna_up <- sum(significant_features$mirnas$log2FoldChange > 0)
        mirna_down <- sum(significant_features$mirnas$log2FoldChange < 0)

        feature_counts <- rbind(feature_counts,
                               data.frame(type = "miRNA", count = mirna_up, direction = "Up"),
                               data.frame(type = "miRNA", count = mirna_down, direction = "Down"))
    }

    if (nrow(feature_counts) > 0) {
        p <- ggplot(feature_counts, aes(x = type, y = count, fill = direction)) +
            geom_bar(stat = "identity", position = "dodge") +
            scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "Total" = "gray")) +
            labs(title = "Multi-Modal Cancer Genomics Overview",
                 subtitle = "Significant molecular features across data types",
                 x = "Data Type", y = "Count", fill = "Direction") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        save_plot(p, "multimodal_overview", width = 10, height = 6)
        return(p)
    }

    return(NULL)
}

create_volcano_plots <- function(significant_features) {
    log_message("Creating volcano plots...")

    plots <- list()

    # mRNA volcano plot
    if (!is.null(significant_features$mrnas)) {
        mrna_volcano <- ggplot(significant_features$mrnas,
                              aes(x = log2FoldChange, y = -log10(padj))) +
            geom_point(alpha = 0.6, size = 1) +
            geom_point(data = significant_features$mrnas %>%
                          filter(padj < PVALUE_CUTOFF, abs(log2FoldChange) > FC_CUTOFF),
                      color = "red", size = 2) +
            geom_hline(yintercept = -log10(PVALUE_CUTOFF), linetype = "dashed", alpha = 0.5) +
            geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed", alpha = 0.5) +
            labs(title = "mRNA Expression Volcano Plot",
                 x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
            theme_minimal()

        save_plot(mrna_volcano, "mrna_volcano_plot", width = 8, height = 6)
        plots$mrna <- mrna_volcano
    }

    # miRNA volcano plot
    if (!is.null(significant_features$mirnas)) {
        mirna_volcano <- ggplot(significant_features$mirnas,
                               aes(x = log2FoldChange, y = -log10(padj))) +
            geom_point(alpha = 0.6, size = 1) +
            geom_point(data = significant_features$mirnas %>%
                          filter(padj < PVALUE_CUTOFF, abs(log2FoldChange) > FC_CUTOFF),
                      color = "red", size = 2) +
            geom_hline(yintercept = -log10(PVALUE_CUTOFF), linetype = "dashed", alpha = 0.5) +
            geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed", alpha = 0.5) +
            labs(title = "miRNA Expression Volcano Plot",
                 x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
            theme_minimal()

        save_plot(mirna_volcano, "mirna_volcano_plot", width = 8, height = 6)
        plots$mirna <- mirna_volcano
    }

    return(plots)
}

create_pathway_plots <- function(pathway_results) {
    log_message("Creating pathway enrichment plots...")

    plots <- list()

    if (!is.null(pathway_results$go)) {
        # GO enrichment plot
        go_plot <- dotplot(pathway_results$go, showCategory = 20) +
            ggtitle("GO Biological Process Enrichment")

        save_plot(go_plot, "go_enrichment", width = 12, height = 8)
        plots$go <- go_plot

        # GO network plot
        if (nrow(pathway_results$go) > 5) {
            tryCatch({
                go_network <- emapplot(pairwise_termsim(pathway_results$go),
                                      showCategory = 15) +
                    ggtitle("GO Term Network")

                save_plot(go_network, "go_network", width = 12, height = 10)
                plots$go_network <- go_network
            }, error = function(e) {
                log_message(paste("Error creating GO network plot:", e$message))
            })
        }
    }

    if (!is.null(pathway_results$kegg)) {
        # KEGG enrichment plot
        kegg_plot <- dotplot(pathway_results$kegg, showCategory = 20) +
            ggtitle("KEGG Pathway Enrichment")

        save_plot(kegg_plot, "kegg_enrichment", width = 12, height = 8)
        plots$kegg <- kegg_plot
    }

    return(plots)
}

create_integration_heatmap <- function(significant_features, mrna_data) {
    log_message("Creating integration heatmap...")

    if (!is.null(mrna_data$dds) && !is.null(significant_features$mrnas)) {
        # Get top significant genes
        top_genes <- significant_features$mrnas %>%
            arrange(padj) %>%
            head(50) %>%
            pull(gene_name)

        top_genes <- top_genes[!is.na(top_genes)]

        if (length(top_genes) > 5) {
            tryCatch({
                # Get normalized counts
                norm_counts <- counts(mrna_data$dds, normalized = TRUE)

                # Filter to top genes
                heatmap_data <- norm_counts[rownames(norm_counts) %in% top_genes, ]

                # Log transform
                heatmap_data <- log2(heatmap_data + 1)

                # Create annotation
                sample_info <- colData(mrna_data$dds)
                annotation_col <- data.frame(
                    condition = sample_info$condition,
                    row.names = rownames(sample_info)
                )

                # Create heatmap
                png("09_multimodal/plots/integration_heatmap.png",
                    width = 12, height = 10, units = "in", res = DPI)

                pheatmap(
                    heatmap_data,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    scale = "row",
                    annotation_col = annotation_col,
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    main = "Top Differentially Expressed Genes",
                    fontsize_row = 8
                )

                dev.off()
                log_message("Integration heatmap saved")

            }, error = function(e) {
                log_message(paste("Error creating integration heatmap:", e$message))
            })
        }
    }
}

# === REPORT GENERATION ===

generate_html_report <- function(significant_features, pathway_results) {
    log_message("Generating HTML report...")

    # Create summary statistics
    summary_stats <- list(
        mrna_total = ifelse(!is.null(significant_features$mrnas), nrow(significant_features$mrnas), 0),
        mrna_up = ifelse(!is.null(significant_features$mrnas),
                        sum(significant_features$mrnas$log2FoldChange > 0), 0),
        mrna_down = ifelse(!is.null(significant_features$mrnas),
                          sum(significant_features$mrnas$log2FoldChange < 0), 0),
        variants_total = ifelse(!is.null(significant_features$variants),
                               nrow(significant_features$variants), 0),
        mirna_total = ifelse(!is.null(significant_features$mirnas),
                            nrow(significant_features$mirnas), 0),
        mirna_up = ifelse(!is.null(significant_features$mirnas),
                         sum(significant_features$mirnas$log2FoldChange > 0), 0),
        mirna_down = ifelse(!is.null(significant_features$mirnas),
                           sum(significant_features$mirnas$log2FoldChange < 0), 0),
        go_terms = ifelse(!is.null(pathway_results$go), nrow(pathway_results$go), 0),
        kegg_pathways = ifelse(!is.null(pathway_results$kegg), nrow(pathway_results$kegg), 0)
    )

    # Save summary statistics
    write_csv(as.data.frame(summary_stats), "09_multimodal/tables/summary_statistics.csv")

    # Save detailed results tables
    if (!is.null(significant_features$mrnas)) {
        write_csv(significant_features$mrnas, "09_multimodal/tables/significant_mrnas.csv")
    }

    if (!is.null(significant_features$variants)) {
        write_csv(significant_features$variants, "09_multimodal/tables/significant_variants.csv")
    }

    if (!is.null(significant_features$mirnas)) {
        write_csv(significant_features$mirnas, "09_multimodal/tables/significant_mirnas.csv")
    }

    if (!is.null(pathway_results$go)) {
        write_csv(as.data.frame(pathway_results$go), "09_multimodal/tables/go_enrichment.csv")
    }

    if (!is.null(pathway_results$kegg)) {
        write_csv(as.data.frame(pathway_results$kegg), "09_multimodal/tables/kegg_enrichment.csv")
    }

    log_message("Results tables saved to 09_multimodal/tables/")
}

# === MAIN ANALYSIS WORKFLOW ===

main <- function() {
    log_message("Starting multi-modal cancer genomics analysis...")

    # Load all data types
    mrna_data <- load_mrna_data()
    variant_data <- load_variant_data()
    mirna_data <- load_mirna_data()

    # Identify significant features across all data types
    significant_features <- identify_significant_features(mrna_data, variant_data, mirna_data)

    # Perform pathway analysis
    pathway_results <- perform_pathway_analysis(significant_features)

    # Create visualizations
    create_overview_plot(significant_features)
    create_volcano_plots(significant_features)
    create_pathway_plots(pathway_results)
    create_integration_heatmap(significant_features, mrna_data)

    # Generate reports
    generate_html_report(significant_features, pathway_results)

    log_message("Multi-modal analysis completed successfully!")
    log_message("Results saved in 09_multimodal/ directory")

    # Print summary
    cat("\n=== ANALYSIS SUMMARY ===\n")
    if (!is.null(significant_features$mrnas)) {
        cat("Significant mRNAs:", nrow(significant_features$mrnas), "\n")
    }
    if (!is.null(significant_features$variants)) {
        cat("High-quality variants:", nrow(significant_features$variants), "\n")
    }
    if (!is.null(significant_features$mirnas)) {
        cat("Significant miRNAs:", nrow(significant_features$mirnas), "\n")
    }
    if (!is.null(pathway_results$go)) {
        cat("GO enrichment terms:", nrow(pathway_results$go), "\n")
    }
    if (!is.null(pathway_results$kegg)) {
        cat("KEGG pathways:", nrow(pathway_results$kegg), "\n")
    }
    cat("\nPlots and tables saved in 09_multimodal/\n")
}

# Execute main analysis
if (!interactive()) {
    main()
}