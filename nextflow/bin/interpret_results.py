#!/usr/bin/env python3

"""
Automated RNA-Seq Downstream Analysis and Visualization in Python

Description: This script performs differential expression analysis using PyDESeq2,
             generates publication-quality visualizations (PCA, Volcano, Heatmap),
             and conducts functional enrichment analysis using GSEApy.

Usage: python interpret_results.py <counts_file> <metadata_file> <output_dir>

Example: python interpret_results.py raw_counts.tsv metadata.tsv python_results

Requirements: pydeseq2, pandas, numpy, matplotlib, seaborn, sklearn, gseapy
"""

import sys
import os
import logging
import warnings
from datetime import datetime
from pathlib import Path

# Data manipulation and analysis
import pandas as pd
import numpy as np

# Statistical analysis
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

# Bioinformatics
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    print("Warning: PyDESeq2 not available. Using alternative analysis methods.")
    PYDESEQ2_AVAILABLE = False

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    print("Warning: GSEApy not available. Skipping enrichment analysis.")
    GSEAPY_AVAILABLE = False

# Suppress warnings
warnings.filterwarnings('ignore')

# === ANALYSIS PARAMETERS ===
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0
MIN_BASEMEAN = 10

# === SETUP LOGGING ===
def setup_logging(output_dir):
    """Setup logging for the analysis"""
    log_file = os.path.join(output_dir, 'analysis.log')
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# === UTILITY FUNCTIONS ===
def validate_inputs(counts_file, metadata_file):
    """Validate that input files exist and are readable"""
    if not os.path.exists(counts_file):
        raise FileNotFoundError(f"Counts file not found: {counts_file}")

    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    return True

def clean_ensembl_ids(gene_ids):
    """Remove version numbers from Ensembl gene IDs"""
    return [gene_id.split('.')[0] for gene_id in gene_ids]

# === DATA LOADING AND PREPROCESSING ===
def load_count_data(counts_file):
    """Load and preprocess count data from featureCounts output"""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading count data from: {counts_file}")

    # Load the count data
    count_df = pd.read_csv(counts_file, sep='\\t', index_col=0)

    # Extract only the sample count columns (skip first 5 metadata columns in featureCounts output)
    # featureCounts columns: Geneid, Chr, Start, End, Strand, Length, then sample columns
    if count_df.shape[1] < 6:
        raise ValueError("Invalid count file format. Expected featureCounts output format.")

    # Extract sample count columns
    count_matrix = count_df.iloc[:, 5:]  # Skip first 5 metadata columns

    # Clean up column names
    count_matrix.columns = [col.replace('_Aligned.sortedByCoord.out.bam', '')
                           for col in count_matrix.columns]
    count_matrix.columns = [col.split('/')[-1].replace('.bam', '')
                           for col in count_matrix.columns]

    # Ensure counts are integers
    count_matrix = count_matrix.round().astype(int)

    logger.info(f"Loaded count matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} samples")
    logger.info(f"Sample names: {', '.join(count_matrix.columns)}")

    return count_matrix

def load_metadata(metadata_file):
    """Load sample metadata"""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading metadata from: {metadata_file}")

    metadata_df = pd.read_csv(metadata_file, sep='\\t', index_col=0)
    logger.info(f"Loaded metadata for {metadata_df.shape[0]} samples")

    return metadata_df

def align_samples(count_matrix, metadata_df):
    """Align samples between count matrix and metadata"""
    logger = logging.getLogger(__name__)

    # Find common samples
    common_samples = set(count_matrix.columns) & set(metadata_df.index)

    if len(common_samples) == 0:
        raise ValueError("No matching samples found between count matrix and metadata")

    if len(common_samples) < len(count_matrix.columns):
        missing_samples = set(count_matrix.columns) - set(metadata_df.index)
        logger.warning(f"Samples in count matrix not found in metadata: {missing_samples}")

    # Subset to common samples
    common_samples = list(common_samples)
    count_matrix = count_matrix[common_samples]
    metadata_df = metadata_df.loc[common_samples]

    logger.info(f"Proceeding with {len(common_samples)} common samples")

    return count_matrix, metadata_df

# === DIFFERENTIAL EXPRESSION ANALYSIS ===
def run_pydeseq2_analysis(count_matrix, metadata_df):
    """Run differential expression analysis using PyDESeq2"""
    logger = logging.getLogger(__name__)
    logger.info("Running PyDESeq2 differential expression analysis")

    # Check for condition column
    if 'condition' not in metadata_df.columns:
        raise ValueError("'condition' column not found in metadata")

    # Get unique conditions
    unique_conditions = metadata_df['condition'].unique()
    logger.info(f"Found conditions: {', '.join(unique_conditions)}")

    if len(unique_conditions) != 2:
        raise ValueError(f"Expected exactly 2 conditions, found {len(unique_conditions)}")

    # Set reference level (Healthy if available, otherwise first alphabetically)
    if 'Healthy' in unique_conditions:
        reference_level = 'Healthy'
        test_level = [c for c in unique_conditions if c != 'Healthy'][0]
    else:
        reference_level = sorted(unique_conditions)[0]
        test_level = sorted(unique_conditions)[1]

    logger.info(f"Comparison: {test_level} vs {reference_level} (reference)")

    # Filter low count genes
    keep_genes = count_matrix.sum(axis=1) >= MIN_BASEMEAN
    count_matrix_filtered = count_matrix[keep_genes]
    logger.info(f"Filtered to {count_matrix_filtered.shape[0]} genes with sufficient expression")

    # Create DeseqDataSet
    dds = DeseqDataSet(
        counts=count_matrix_filtered.T,  # PyDESeq2 expects samples as rows
        metadata=metadata_df,
        design_factors="condition",
        ref_level=["condition", reference_level]
    )

    # Run DESeq2
    dds.deseq2()

    # Get results
    stat_res = DeseqStats(dds, contrast=["condition", test_level, reference_level])
    stat_res.summary()

    # Extract results as DataFrame
    results_df = stat_res.results_df.copy()

    logger.info(f"DESeq2 analysis complete. {results_df.shape[0]} genes tested.")

    return results_df, dds

def run_alternative_analysis(count_matrix, metadata_df):
    """Run alternative differential expression analysis if PyDESeq2 is not available"""
    logger = logging.getLogger(__name__)
    logger.info("Running alternative differential expression analysis")

    # This is a simplified analysis - in practice, you'd want to use proper statistical methods
    # Get conditions
    unique_conditions = metadata_df['condition'].unique()
    if len(unique_conditions) != 2:
        raise ValueError(f"Expected exactly 2 conditions, found {len(unique_conditions)}")

    condition1_samples = metadata_df[metadata_df['condition'] == unique_conditions[0]].index
    condition2_samples = metadata_df[metadata_df['condition'] == unique_conditions[1]].index

    # Calculate log2 fold change and p-values using t-test
    results = []

    for gene in count_matrix.index:
        group1 = count_matrix.loc[gene, condition1_samples] + 1  # Add pseudocount
        group2 = count_matrix.loc[gene, condition2_samples] + 1

        # Calculate means
        mean1 = group1.mean()
        mean2 = group2.mean()

        # Calculate log2 fold change
        log2fc = np.log2(mean2 / mean1)

        # Calculate p-value using t-test
        try:
            _, pval = stats.ttest_ind(group1, group2)
        except:
            pval = 1.0

        results.append({
            'gene': gene,
            'baseMean': (mean1 + mean2) / 2,
            'log2FoldChange': log2fc,
            'pvalue': pval,
            'padj': pval  # Simplified - should use proper multiple testing correction
        })

    results_df = pd.DataFrame(results).set_index('gene')

    # Apply multiple testing correction
    from scipy.stats import false_discovery_control
    results_df['padj'] = false_discovery_control(results_df['pvalue'])

    logger.info(f"Alternative analysis complete. {results_df.shape[0]} genes tested.")

    return results_df, None

# === VISUALIZATION FUNCTIONS ===
def create_pca_plot(count_matrix, metadata_df, output_dir):
    """Create PCA plot"""
    logger = logging.getLogger(__name__)
    logger.info("Creating PCA plot")

    # Normalize and transform data
    # Add pseudocount and log transform
    log_counts = np.log2(count_matrix + 1)

    # Standardize features
    scaler = StandardScaler()
    scaled_counts = scaler.fit_transform(log_counts.T)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_counts)

    # Create PCA DataFrame
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=count_matrix.columns)
    pca_df = pca_df.join(metadata_df['condition'])

    # Create plot
    plt.figure(figsize=(10, 8))

    # Plot points
    for condition in pca_df['condition'].unique():
        condition_data = pca_df[pca_df['condition'] == condition]
        plt.scatter(condition_data['PC1'], condition_data['PC2'],
                   label=condition, s=100, alpha=0.7)

    # Add sample labels
    for i, (idx, row) in enumerate(pca_df.iterrows()):
        plt.annotate(idx, (row['PC1'], row['PC2']),
                    xytext=(5, 5), textcoords='offset points', fontsize=9)

    plt.xlabel(f'PC1: {pca.explained_variance_ratio_[0]:.1%} variance')
    plt.ylabel(f'PC2: {pca.explained_variance_ratio_[1]:.1%} variance')
    plt.title('Principal Component Analysis')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    pca_file = os.path.join(output_dir, 'pca_plot.png')
    plt.savefig(pca_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"PCA plot saved: {pca_file}")

    return pca_file

def create_volcano_plot(results_df, output_dir, test_condition="Test", reference_condition="Reference"):
    """Create volcano plot"""
    logger = logging.getLogger(__name__)
    logger.info("Creating volcano plot")

    # Prepare data
    plot_data = results_df.copy()
    plot_data = plot_data.dropna()

    # Add significance categories
    plot_data['significant'] = 'Not Significant'
    plot_data.loc[(plot_data['log2FoldChange'] > LOG2FC_THRESHOLD) &
                  (plot_data['padj'] < PADJ_THRESHOLD), 'significant'] = 'Upregulated'
    plot_data.loc[(plot_data['log2FoldChange'] < -LOG2FC_THRESHOLD) &
                  (plot_data['padj'] < PADJ_THRESHOLD), 'significant'] = 'Downregulated'

    # Create plot
    plt.figure(figsize=(12, 10))

    # Color mapping
    colors = {'Not Significant': 'gray', 'Upregulated': 'red', 'Downregulated': 'blue'}

    for category in ['Not Significant', 'Downregulated', 'Upregulated']:
        data = plot_data[plot_data['significant'] == category]
        plt.scatter(data['log2FoldChange'], -np.log10(data['padj']),
                   c=colors[category], label=category, alpha=0.6, s=20)

    # Add threshold lines
    plt.axvline(x=LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.7)
    plt.axvline(x=-LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.7)
    plt.axhline(y=-np.log10(PADJ_THRESHOLD), color='black', linestyle='--', alpha=0.7)

    # Label top genes
    top_genes = plot_data[plot_data['significant'] != 'Not Significant'].nlargest(10, 'padj')
    for idx, row in top_genes.iterrows():
        plt.annotate(idx[:15], (row['log2FoldChange'], -np.log10(row['padj'])),
                    xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8)

    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title(f'Volcano Plot: {test_condition} vs {reference_condition}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    volcano_file = os.path.join(output_dir, 'volcano_plot.png')
    plt.savefig(volcano_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Volcano plot saved: {volcano_file}")

    return volcano_file

def create_heatmap(count_matrix, metadata_df, significant_genes, output_dir, n_genes=50):
    """Create heatmap of top differentially expressed genes"""
    logger = logging.getLogger(__name__)
    logger.info("Creating heatmap of top differentially expressed genes")

    if len(significant_genes) == 0:
        logger.warning("No significant genes found for heatmap")
        return None

    # Select top genes
    n_genes = min(n_genes, len(significant_genes))
    top_genes = significant_genes.head(n_genes).index

    # Get normalized counts (log2 + 1 transformation)
    heatmap_data = np.log2(count_matrix.loc[top_genes] + 1)

    # Z-score normalization across genes (rows)
    heatmap_data_scaled = (heatmap_data.T - heatmap_data.mean(axis=1)) / heatmap_data.std(axis=1)
    heatmap_data_scaled = heatmap_data_scaled.T

    # Create annotation for samples
    sample_colors = []
    unique_conditions = metadata_df['condition'].unique()
    color_palette = sns.color_palette("Set1", len(unique_conditions))
    condition_colors = dict(zip(unique_conditions, color_palette))

    for sample in heatmap_data_scaled.columns:
        condition = metadata_df.loc[sample, 'condition']
        sample_colors.append(condition_colors[condition])

    # Create the plot
    plt.figure(figsize=(12, 16))

    # Create heatmap
    sns.heatmap(heatmap_data_scaled,
                cmap='RdBu_r',
                center=0,
                xticklabels=True,
                yticklabels=False,
                cbar_kws={'label': 'Z-score'})

    # Add color bar for conditions
    ax = plt.gca()
    for i, (condition, color) in enumerate(condition_colors.items()):
        ax.add_patch(Rectangle((0, -0.5 - i*0.3), len(heatmap_data_scaled.columns), 0.2,
                              facecolor=color, clip_on=False))
        ax.text(-1, -0.4 - i*0.3, condition, fontsize=10, ha='right', va='center')

    plt.title(f'Top {n_genes} Differentially Expressed Genes')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.tight_layout()

    heatmap_file = os.path.join(output_dir, 'heatmap_top_degs.png')
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Heatmap saved: {heatmap_file}")

    return heatmap_file

# === FUNCTIONAL ENRICHMENT ANALYSIS ===
def run_enrichment_analysis(gene_lists, output_dir):
    """Run functional enrichment analysis using GSEApy"""
    logger = logging.getLogger(__name__)

    if not GSEAPY_AVAILABLE:
        logger.warning("GSEApy not available. Skipping enrichment analysis.")
        return

    logger.info("Running functional enrichment analysis")

    gene_sets = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
                 'KEGG_2021_Human', 'WikiPathway_2023_Human']

    for direction, genes in gene_lists.items():
        if len(genes) == 0:
            continue

        logger.info(f"Enrichment analysis for {direction} genes ({len(genes)} genes)")

        # Clean gene IDs (remove version numbers)
        clean_genes = clean_ensembl_ids(genes)

        try:
            # Run enrichment analysis
            enr = gp.enrichr(
                gene_list=clean_genes,
                gene_sets=gene_sets,
                organism='human',
                outdir=os.path.join(output_dir, f'enrichment_{direction}'),
                cutoff=0.05
            )

            logger.info(f"Enrichment analysis completed for {direction} genes")

        except Exception as e:
            logger.error(f"Enrichment analysis failed for {direction} genes: {e}")

# === MAIN ANALYSIS FUNCTION ===
def main(counts_file, metadata_file, output_dir):
    """Main analysis function"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Setup logging
    logger = setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("RNA-Seq Differential Expression Analysis - Python")
    logger.info("=" * 60)
    logger.info(f"Analysis started: {datetime.now()}")
    logger.info(f"Counts file: {counts_file}")
    logger.info(f"Metadata file: {metadata_file}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"PyDESeq2 available: {PYDESEQ2_AVAILABLE}")
    logger.info(f"GSEApy available: {GSEAPY_AVAILABLE}")

    # Validate inputs
    validate_inputs(counts_file, metadata_file)

    # Load data
    count_matrix = load_count_data(counts_file)
    metadata_df = load_metadata(metadata_file)
    count_matrix, metadata_df = align_samples(count_matrix, metadata_df)

    # Run differential expression analysis
    if PYDESEQ2_AVAILABLE:
        results_df, dds = run_pydeseq2_analysis(count_matrix, metadata_df)
    else:
        results_df, dds = run_alternative_analysis(count_matrix, metadata_df)

    # Save full results
    results_file = os.path.join(output_dir, 'full_results.csv')
    results_df.to_csv(results_file)
    logger.info(f"Full results saved: {results_file}")

    # Identify significant genes
    significant_genes = results_df[
        (results_df['padj'] < PADJ_THRESHOLD) &
        (abs(results_df['log2FoldChange']) > LOG2FC_THRESHOLD)
    ].sort_values('padj')

    up_regulated = significant_genes[significant_genes['log2FoldChange'] > 0]
    down_regulated = significant_genes[significant_genes['log2FoldChange'] < 0]

    logger.info(f"Significant genes found: {len(significant_genes)}")
    logger.info(f"Upregulated: {len(up_regulated)}")
    logger.info(f"Downregulated: {len(down_regulated)}")

    # Save gene lists
    if len(up_regulated) > 0:
        up_file = os.path.join(output_dir, 'up_regulated_genes.csv')
        up_regulated.to_csv(up_file)
        logger.info(f"Upregulated genes saved: {up_file}")

    if len(down_regulated) > 0:
        down_file = os.path.join(output_dir, 'down_regulated_genes.csv')
        down_regulated.to_csv(down_file)
        logger.info(f"Downregulated genes saved: {down_file}")

    # Generate visualizations
    logger.info("Generating visualizations...")

    # PCA plot
    create_pca_plot(count_matrix, metadata_df, output_dir)

    # Volcano plot
    conditions = metadata_df['condition'].unique()
    test_condition = conditions[0] if 'Healthy' not in conditions else [c for c in conditions if c != 'Healthy'][0]
    reference_condition = 'Healthy' if 'Healthy' in conditions else conditions[1]
    create_volcano_plot(results_df, output_dir, test_condition, reference_condition)

    # Heatmap
    if len(significant_genes) > 0:
        create_heatmap(count_matrix, metadata_df, significant_genes, output_dir)

    # Enrichment analysis
    gene_lists = {
        'upregulated': up_regulated.index.tolist(),
        'downregulated': down_regulated.index.tolist()
    }
    run_enrichment_analysis(gene_lists, output_dir)

    # Generate summary report
    summary_file = os.path.join(output_dir, 'analysis_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("RNA-Seq Differential Expression Analysis Summary (Python)\\n")
        f.write("=" * 60 + "\\n\\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        f.write(f"Comparison: {test_condition} vs {reference_condition} (reference)\\n\\n")

        f.write("Input Data:\\n")
        f.write(f"- Count file: {counts_file}\\n")
        f.write(f"- Metadata file: {metadata_file}\\n")
        f.write(f"- Total genes: {count_matrix.shape[0]}\\n")
        f.write(f"- Total samples: {count_matrix.shape[1]}\\n\\n")

        f.write("Analysis Parameters:\\n")
        f.write(f"- Adjusted p-value threshold: {PADJ_THRESHOLD}\\n")
        f.write(f"- Log2 fold change threshold: {LOG2FC_THRESHOLD}\\n")
        f.write(f"- Minimum base mean: {MIN_BASEMEAN}\\n\\n")

        f.write("Results:\\n")
        f.write(f"- Total significant genes: {len(significant_genes)}\\n")
        f.write(f"- Upregulated genes: {len(up_regulated)}\\n")
        f.write(f"- Downregulated genes: {len(down_regulated)}\\n\\n")

        f.write("Output Files:\\n")
        f.write("- full_results.csv: Complete analysis results\\n")
        f.write("- up_regulated_genes.csv: Significantly upregulated genes\\n")
        f.write("- down_regulated_genes.csv: Significantly downregulated genes\\n")
        f.write("- pca_plot.png: Principal component analysis\\n")
        f.write("- volcano_plot.png: Volcano plot of differential expression\\n")
        f.write("- heatmap_top_degs.png: Heatmap of top differentially expressed genes\\n")

        if GSEAPY_AVAILABLE:
            f.write("- enrichment_*/: Functional enrichment analysis results\\n")

    logger.info(f"Analysis summary saved: {summary_file}")
    logger.info("=" * 60)
    logger.info("Analysis completed successfully!")
    logger.info(f"Results saved in: {output_dir}")

if __name__ == '__main__':
    # Parse command line arguments
    if len(sys.argv) == 1:
        # Default parameters for testing
        print("No arguments provided. Using default test parameters.")
        sys.argv = ['interpret_results.py',
                   'results/05_featurecounts/raw_counts.tsv',
                   'config/metadata.tsv',
                   'results/analysis_python']

    if len(sys.argv) != 4:
        print("ERROR: Incorrect number of arguments provided.")
        print("Usage: python interpret_results.py <counts_file> <metadata_file> <output_dir>")
        print("Example: python interpret_results.py raw_counts.tsv metadata.tsv python_results")
        sys.exit(1)

    counts_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_dir = sys.argv[3]

    try:
        main(counts_file, metadata_file, output_dir)
    except Exception as e:
        print(f"ERROR: Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)