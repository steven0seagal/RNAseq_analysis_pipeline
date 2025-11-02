process DESEQ2_ANALYSIS {
    tag "DESeq2"
    label 'process_medium'

    conda "conda-forge::r-base=4.2.0 bioconda::bioconductor-deseq2=1.38.0 bioconda::bioconductor-clusterprofiler=4.6.0 conda-forge::r-ggplot2 conda-forge::r-pheatmap conda-forge::r-ggrepel"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path counts
    path metadata
    val control_condition
    val treatment_condition
    val padj_threshold
    val log2fc_threshold

    output:
    path "deseq2_results.csv"              , emit: results
    path "deseq2_results_significant.csv"  , emit: significant
    path "up_regulated_genes.csv"          , emit: up_genes
    path "down_regulated_genes.csv"        , emit: down_genes
    path "pca_plot.png"                    , emit: pca_plot
    path "volcano_plot.png"                , emit: volcano_plot
    path "heatmap_top50.png"               , emit: heatmap
    path "ma_plot.png"                     , emit: ma_plot, optional: true
    path "*.pdf"                           , emit: plots_pdf, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(DESeq2)
        library(ggplot2)
        library(pheatmap)
        library(RColorBrewer)
        library(ggrepel)
    })

    # Set parameters
    PADJ_THRESHOLD <- ${padj_threshold}
    LOG2FC_THRESHOLD <- ${log2fc_threshold}
    CONTROL <- "${control_condition}"
    TREATMENT <- "${treatment_condition}"

    # Load count data
    cat("Loading count data...\\n")
    counts_data <- read.table("${counts}", header = TRUE, row.names = 1, sep = "\\t", comment.char = "#")

    # Remove annotation columns (Chr, Start, End, Strand, Length)
    if (ncol(counts_data) > 6) {
        counts_matrix <- as.matrix(counts_data[, 7:ncol(counts_data)])
    } else {
        counts_matrix <- as.matrix(counts_data)
    }

    # Clean column names (remove paths and .bam extension)
    colnames(counts_matrix) <- gsub(".*\\\\/", "", colnames(counts_matrix))
    colnames(counts_matrix) <- gsub("\\\\..*\\\\.bam\$", "", colnames(counts_matrix))
    colnames(counts_matrix) <- gsub("\\\\.bam\$", "", colnames(counts_matrix))

    # Load metadata
    cat("Loading metadata...\\n")
    metadata <- read.table("${metadata}", header = TRUE, sep = "\\t", row.names = 1)

    # Ensure samples match between counts and metadata
    common_samples <- intersect(colnames(counts_matrix), rownames(metadata))
    counts_matrix <- counts_matrix[, common_samples]
    metadata <- metadata[common_samples, , drop = FALSE]

    cat(paste("Analyzing", ncol(counts_matrix), "samples\\n"))
    cat(paste("Control condition:", CONTROL, "\\n"))
    cat(paste("Treatment condition:", TREATMENT, "\\n"))

    # Create DESeq2 dataset
    cat("Creating DESeq2 dataset...\\n")
    dds <- DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = metadata,
        design = ~ condition
    )

    # Filter low count genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat(paste("Retained", sum(keep), "genes after filtering\\n"))

    # Run DESeq2 analysis
    cat("Running differential expression analysis...\\n")
    dds <- DESeq(dds)

    # Set reference level
    dds\$condition <- relevel(dds\$condition, ref = CONTROL)

    # Get results
    res <- results(dds,
                   contrast = c("condition", TREATMENT, CONTROL),
                   alpha = PADJ_THRESHOLD)

    # Convert to dataframe and add gene names
    res_df <- as.data.frame(res)
    res_df\$gene <- rownames(res_df)
    res_df <- res_df[order(res_df\$padj), ]

    # Save all results
    write.csv(res_df, "deseq2_results.csv", row.names = FALSE)

    # Filter significant genes
    sig_genes <- subset(res_df, !is.na(padj) & padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD)
    write.csv(sig_genes, "deseq2_results_significant.csv", row.names = FALSE)

    # Separate up and down regulated
    up_genes <- subset(sig_genes, log2FoldChange > 0)
    down_genes <- subset(sig_genes, log2FoldChange < 0)
    write.csv(up_genes, "up_regulated_genes.csv", row.names = FALSE)
    write.csv(down_genes, "down_regulated_genes.csv", row.names = FALSE)

    cat(paste("Found", nrow(sig_genes), "significant genes\\n"))
    cat(paste("  Up-regulated:", nrow(up_genes), "\\n"))
    cat(paste("  Down-regulated:", nrow(down_genes), "\\n"))

    # ===== VISUALIZATIONS =====

    # 1. PCA Plot
    cat("Generating PCA plot...\\n")
    vsd <- vst(dds, blind = FALSE)
    png("pca_plot.png", width = 1200, height = 1000, res = 150)
    pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text_repel(size = 3, max.overlaps = 20) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("PCA Plot - Sample Clustering") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
              legend.position = "right")
    dev.off()

    # 2. Volcano Plot
    cat("Generating volcano plot...\\n")
    png("volcano_plot.png", width = 1400, height = 1200, res = 150)

    # Prepare data for volcano plot
    volcano_data <- res_df
    volcano_data\$diffexpressed <- "NO"
    volcano_data\$diffexpressed[volcano_data\$log2FoldChange > LOG2FC_THRESHOLD & volcano_data\$padj < PADJ_THRESHOLD] <- "UP"
    volcano_data\$diffexpressed[volcano_data\$log2FoldChange < -LOG2FC_THRESHOLD & volcano_data\$padj < PADJ_THRESHOLD] <- "DOWN"
    volcano_data\$delabel <- NA

    # Label top genes
    top_genes <- head(volcano_data[order(volcano_data\$padj), ], 15)
    volcano_data\$delabel[rownames(volcano_data) %in% rownames(top_genes)] <- top_genes\$gene

    ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_text_repel(size = 3, max.overlaps = 15) +
        scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
                          labels = c(paste("Down (", nrow(down_genes), ")", sep=""),
                                   "Not Significant",
                                   paste("Up (", nrow(up_genes), ")", sep=""))) +
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), col = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(PADJ_THRESHOLD), col = "black", linetype = "dashed") +
        labs(title = "Volcano Plot - Differential Expression",
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value",
             color = "Expression") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
              legend.position = "right")
    dev.off()

    # 3. Heatmap of top 50 DEGs
    cat("Generating heatmap...\\n")
    if (nrow(sig_genes) > 0) {
        top_genes_heatmap <- head(sig_genes\$gene, 50)
        mat <- assay(vsd)[top_genes_heatmap, ]
        mat <- t(scale(t(mat)))

        png("heatmap_top50.png", width = 1400, height = 1600, res = 150)
        pheatmap(mat,
                 annotation_col = as.data.frame(colData(dds)["condition"]),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 main = "Top 50 Differentially Expressed Genes",
                 fontsize_row = 8,
                 fontsize_col = 10)
        dev.off()
    }

    # 4. MA Plot
    cat("Generating MA plot...\\n")
    png("ma_plot.png", width = 1200, height = 1000, res = 150)
    plotMA(res, ylim = c(-5, 5),
           main = "MA Plot - Mean Expression vs Log Fold Change",
           colNonSig = "gray60",
           colSig = "red",
           colLine = "blue")
    dev.off()

    # Version information
    cat("\\nGenerating version information...\\n")
    version_file <- file("versions.yml", "w")
    cat('"${task.process}":\\n', file = version_file)
    cat(paste0('    deseq2: "', as.character(packageVersion("DESeq2")), '"\\n'), file = version_file)
    cat(paste0('    r-base: "', R.version\$major, ".", R.version\$minor, '"\\n'), file = version_file)
    close(version_file)

    cat("\\nDESeq2 analysis complete!\\n")
    """

    stub:
    """
    touch deseq2_results.csv
    touch deseq2_results_significant.csv
    touch up_regulated_genes.csv
    touch down_regulated_genes.csv
    touch pca_plot.png
    touch volcano_plot.png
    touch heatmap_top50.png
    touch ma_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: "1.38.0"
        r-base: "4.2.0"
    END_VERSIONS
    """
}
