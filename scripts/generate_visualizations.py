#!/usr/bin/env python3

"""
Comprehensive Visualization Suite for Multi-Modal Cancer Genomics Analysis

This script generates publication-ready visualizations for:
- mRNA expression analysis results
- Somatic variant landscapes
- miRNA expression patterns
- Multi-modal data integration
- Interactive plots and dashboards

Usage: python generate_visualizations.py [--config config.yaml]
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import yaml
import os
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-ready plots
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans']
})

class MultiModalVisualizer:
    """
    Comprehensive visualization suite for multi-modal cancer genomics data
    """

    def __init__(self, config_file=None):
        """Initialize the visualizer with configuration"""
        self.config = self.load_config(config_file)
        self.setup_directories()
        self.data = {}

    def load_config(self, config_file):
        """Load configuration from YAML file"""
        default_config = {
            'output_dir': '10_visualizations',
            'dpi': 300,
            'figure_format': 'png',
            'interactive_format': 'html',
            'color_palette': 'husl',
            'font_size': 12
        }

        if config_file and os.path.exists(config_file):
            with open(config_file, 'r') as f:
                user_config = yaml.safe_load(f)
            default_config.update(user_config.get('visualization', {}))

        return default_config

    def setup_directories(self):
        """Create output directories"""
        self.output_dir = Path(self.config['output_dir'])
        self.plots_dir = self.output_dir / 'plots'
        self.interactive_dir = self.output_dir / 'interactive'
        self.reports_dir = self.output_dir / 'reports'

        for directory in [self.plots_dir, self.interactive_dir, self.reports_dir]:
            directory.mkdir(parents=True, exist_ok=True)

    def log(self, message):
        """Simple logging function"""
        print(f"[INFO] {message}")

    def load_data(self):
        """Load all available data files"""
        self.log("Loading multi-modal data...")

        # mRNA expression data
        mrna_file = "05_quantification/deseq2_results.csv"
        if os.path.exists(mrna_file):
            self.data['mrna'] = pd.read_csv(mrna_file)
            self.log(f"Loaded {len(self.data['mrna'])} mRNA expression results")

        # Variant data summary
        variant_file = "06_variants/somatic_variant_summary.txt"
        if os.path.exists(variant_file):
            self.data['variants'] = pd.read_csv(variant_file, sep='\t')
            self.log(f"Loaded variant data for {len(self.data['variants'])} samples")

        # miRNA expression data
        mirna_file = "08_small_rna/differential_mirnas.csv"
        if os.path.exists(mirna_file):
            self.data['mirna'] = pd.read_csv(mirna_file)
            self.log(f"Loaded {len(self.data['mirna'])} miRNA expression results")

        # Multi-modal summary
        summary_file = "09_multimodal/tables/summary_statistics.csv"
        if os.path.exists(summary_file):
            self.data['summary'] = pd.read_csv(summary_file)
            self.log("Loaded multi-modal summary statistics")

        return len(self.data) > 0

    def create_expression_plots(self):
        """Create mRNA and miRNA expression visualizations"""
        self.log("Creating expression analysis plots...")

        # mRNA Volcano Plot
        if 'mrna' in self.data:
            self.create_volcano_plot(self.data['mrna'], 'mRNA', 'mrna_volcano')

        # miRNA Volcano Plot
        if 'mirna' in self.data:
            self.create_volcano_plot(self.data['mirna'], 'miRNA', 'mirna_volcano')

        # Combined expression heatmap
        if 'mrna' in self.data and 'mirna' in self.data:
            self.create_combined_heatmap()

        # MA plots
        if 'mrna' in self.data:
            self.create_ma_plot(self.data['mrna'], 'mRNA', 'mrna_ma')

    def create_volcano_plot(self, data, data_type, filename):
        """Create volcano plot for expression data"""
        if 'log2FoldChange' not in data.columns or 'padj' not in data.columns:
            self.log(f"Skipping volcano plot for {data_type} - missing required columns")
            return

        # Filter valid data
        plot_data = data.dropna(subset=['log2FoldChange', 'padj'])
        plot_data = plot_data[plot_data['padj'] > 0]  # Remove zero p-values for log transformation

        if len(plot_data) == 0:
            self.log(f"No valid data for {data_type} volcano plot")
            return

        # Create static plot
        fig, ax = plt.subplots(figsize=(10, 8))

        # Calculate significance
        significant = (plot_data['padj'] < 0.05) & (abs(plot_data['log2FoldChange']) > 1)
        upregulated = significant & (plot_data['log2FoldChange'] > 1)
        downregulated = significant & (plot_data['log2FoldChange'] < -1)

        # Plot points
        ax.scatter(plot_data.loc[~significant, 'log2FoldChange'],
                  -np.log10(plot_data.loc[~significant, 'padj']),
                  c='lightgray', alpha=0.6, s=20, label='Non-significant')

        ax.scatter(plot_data.loc[upregulated, 'log2FoldChange'],
                  -np.log10(plot_data.loc[upregulated, 'padj']),
                  c='red', alpha=0.7, s=30, label='Upregulated')

        ax.scatter(plot_data.loc[downregulated, 'log2FoldChange'],
                  -np.log10(plot_data.loc[downregulated, 'padj']),
                  c='blue', alpha=0.7, s=30, label='Downregulated')

        # Add threshold lines
        ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        # Formatting
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-Log10 Adjusted P-value')
        ax.set_title(f'{data_type} Expression Volcano Plot')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Add statistics
        n_up = sum(upregulated)
        n_down = sum(downregulated)
        ax.text(0.02, 0.98, f'Upregulated: {n_up}\nDownregulated: {n_down}',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()
        plt.savefig(self.plots_dir / f'{filename}.{self.config["figure_format"]}',
                   dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()

        # Create interactive version
        self.create_interactive_volcano(plot_data, data_type, filename)

    def create_interactive_volcano(self, data, data_type, filename):
        """Create interactive volcano plot using Plotly"""
        # Calculate significance
        significant = (data['padj'] < 0.05) & (abs(data['log2FoldChange']) > 1)
        upregulated = significant & (data['log2FoldChange'] > 1)
        downregulated = significant & (data['log2FoldChange'] < -1)

        # Create figure
        fig = go.Figure()

        # Add traces for different categories
        categories = [
            (data.loc[~significant], 'Non-significant', 'lightgray'),
            (data.loc[upregulated], 'Upregulated', 'red'),
            (data.loc[downregulated], 'Downregulated', 'blue')
        ]

        for subset, label, color in categories:
            if len(subset) > 0:
                hover_text = []
                for idx, row in subset.iterrows():
                    gene_name = row.get('gene_name', row.get('miRNA', str(idx)))
                    hover_text.append(
                        f"Gene: {gene_name}<br>"
                        f"Log2FC: {row['log2FoldChange']:.2f}<br>"
                        f"P-adj: {row['padj']:.2e}<br>"
                        f"Base Mean: {row.get('baseMean', 'N/A')}"
                    )

                fig.add_trace(go.Scatter(
                    x=subset['log2FoldChange'],
                    y=-np.log10(subset['padj']),
                    mode='markers',
                    name=label,
                    marker=dict(color=color, size=6, opacity=0.7),
                    text=hover_text,
                    hoverinfo='text'
                ))

        # Add threshold lines
        max_x = max(abs(data['log2FoldChange'].min()), data['log2FoldChange'].max())
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", opacity=0.5)
        fig.add_vline(x=1, line_dash="dash", line_color="black", opacity=0.5)
        fig.add_vline(x=-1, line_dash="dash", line_color="black", opacity=0.5)

        # Update layout
        fig.update_layout(
            title=f'{data_type} Expression Volcano Plot (Interactive)',
            xaxis_title='Log2 Fold Change',
            yaxis_title='-Log10 Adjusted P-value',
            width=800,
            height=600,
            template='plotly_white'
        )

        # Save interactive plot
        pyo.plot(fig, filename=str(self.interactive_dir / f'{filename}_interactive.html'),
                auto_open=False)

    def create_ma_plot(self, data, data_type, filename):
        """Create MA plot (M vs A plot)"""
        if 'log2FoldChange' not in data.columns or 'baseMean' not in data.columns:
            self.log(f"Skipping MA plot for {data_type} - missing required columns")
            return

        # Filter valid data
        plot_data = data.dropna(subset=['log2FoldChange', 'baseMean', 'padj'])
        plot_data = plot_data[plot_data['baseMean'] > 0]

        if len(plot_data) == 0:
            return

        fig, ax = plt.subplots(figsize=(10, 8))

        # Calculate significance
        significant = (plot_data['padj'] < 0.05) & (abs(plot_data['log2FoldChange']) > 1)

        # Plot points
        ax.scatter(np.log10(plot_data.loc[~significant, 'baseMean']),
                  plot_data.loc[~significant, 'log2FoldChange'],
                  c='lightgray', alpha=0.6, s=15, label='Non-significant')

        ax.scatter(np.log10(plot_data.loc[significant, 'baseMean']),
                  plot_data.loc[significant, 'log2FoldChange'],
                  c='red', alpha=0.7, s=20, label='Significant')

        # Add fold change threshold lines
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax.axhline(y=-1, color='black', linestyle='--', alpha=0.5)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)

        # Formatting
        ax.set_xlabel('Log10 Base Mean Expression')
        ax.set_ylabel('Log2 Fold Change')
        ax.set_title(f'{data_type} MA Plot')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(self.plots_dir / f'{filename}.{self.config["figure_format"]}',
                   dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()

    def create_combined_heatmap(self):
        """Create combined heatmap of top mRNA and miRNA changes"""
        self.log("Creating combined expression heatmap...")

        # Get top significant genes from both datasets
        top_mrnas = self.data['mrna'][
            (self.data['mrna']['padj'] < 0.05) &
            (abs(self.data['mrna']['log2FoldChange']) > 1)
        ].nlargest(25, 'log2FoldChange')

        top_mirnas = self.data['mirna'][
            (self.data['mirna']['padj'] < 0.05) &
            (abs(self.data['mirna']['log2FoldChange']) > 1)
        ].nlargest(15, 'log2FoldChange')

        if len(top_mrnas) == 0 and len(top_mirnas) == 0:
            self.log("No significant features found for combined heatmap")
            return

        # Prepare data for heatmap
        heatmap_data = []
        labels = []

        for _, row in top_mrnas.iterrows():
            heatmap_data.append([row['log2FoldChange']])
            labels.append(f"{row.get('gene_name', 'Unknown')} (mRNA)")

        for _, row in top_mirnas.iterrows():
            heatmap_data.append([row['log2FoldChange']])
            labels.append(f"{row.get('miRNA', 'Unknown')} (miRNA)")

        if len(heatmap_data) == 0:
            return

        heatmap_array = np.array(heatmap_data)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(6, max(8, len(labels) * 0.3)))

        sns.heatmap(heatmap_array, yticklabels=labels, xticklabels=['Log2 Fold Change'],
                   cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                   cbar_kws={'label': 'Log2 Fold Change'}, ax=ax)

        ax.set_title('Top Differentially Expressed Features')
        plt.tight_layout()
        plt.savefig(self.plots_dir / f'combined_expression_heatmap.{self.config["figure_format"]}',
                   dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()

    def create_variant_plots(self):
        """Create variant analysis visualizations"""
        if 'variants' not in self.data:
            self.log("No variant data available for plotting")
            return

        self.log("Creating variant analysis plots...")

        variants = self.data['variants']

        # Variant counts per sample
        self.create_variant_barplot(variants)

        # Variant type distribution
        self.create_variant_type_plot(variants)

    def create_variant_barplot(self, variants):
        """Create bar plot of variant counts per sample"""
        # Handle missing values
        variants_clean = variants.fillna(0)

        fig, axes = plt.subplots(2, 2, figsize=(15, 10))

        # PASS variants
        if 'PASS_Variants' in variants_clean.columns:
            axes[0, 0].bar(variants_clean['Sample'], variants_clean['PASS_Variants'])
            axes[0, 0].set_title('PASS Variants per Sample')
            axes[0, 0].set_ylabel('Number of Variants')
            axes[0, 0].tick_params(axis='x', rotation=45)

        # SNVs vs Indels
        if 'SNVs' in variants_clean.columns and 'Indels' in variants_clean.columns:
            x = np.arange(len(variants_clean))
            width = 0.35

            axes[0, 1].bar(x - width/2, variants_clean['SNVs'], width, label='SNVs')
            axes[0, 1].bar(x + width/2, variants_clean['Indels'], width, label='Indels')
            axes[0, 1].set_title('SNVs vs Indels per Sample')
            axes[0, 1].set_ylabel('Number of Variants')
            axes[0, 1].set_xticks(x)
            axes[0, 1].set_xticklabels(variants_clean['Sample'], rotation=45)
            axes[0, 1].legend()

        # Variant proportion pie chart (if only one sample or aggregated)
        if len(variants_clean) == 1:
            snvs = variants_clean['SNVs'].iloc[0] if 'SNVs' in variants_clean.columns else 0
            indels = variants_clean['Indels'].iloc[0] if 'Indels' in variants_clean.columns else 0

            if snvs + indels > 0:
                axes[1, 0].pie([snvs, indels], labels=['SNVs', 'Indels'], autopct='%1.1f%%')
                axes[1, 0].set_title('Variant Type Distribution')

        # Total variant burden
        if 'PASS_Variants' in variants_clean.columns:
            axes[1, 1].hist(variants_clean['PASS_Variants'], bins=min(10, len(variants_clean)),
                          alpha=0.7, edgecolor='black')
            axes[1, 1].set_title('Distribution of Variant Burden')
            axes[1, 1].set_xlabel('Number of PASS Variants')
            axes[1, 1].set_ylabel('Frequency')

        plt.tight_layout()
        plt.savefig(self.plots_dir / f'variant_analysis.{self.config["figure_format"]}',
                   dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()

    def create_variant_type_plot(self, variants):
        """Create detailed variant type analysis"""
        # Interactive variant plot
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Variants per Sample', 'SNV vs Indel Distribution',
                          'Variant Burden Distribution', 'Sample Comparison'),
            specs=[[{"type": "bar"}, {"type": "bar"}],
                   [{"type": "histogram"}, {"type": "scatter"}]]
        )

        # Bar plot of total variants
        if 'PASS_Variants' in variants.columns:
            fig.add_trace(
                go.Bar(x=variants['Sample'], y=variants['PASS_Variants'],
                      name='PASS Variants', marker_color='lightblue'),
                row=1, col=1
            )

        # SNV vs Indel comparison
        if 'SNVs' in variants.columns and 'Indels' in variants.columns:
            fig.add_trace(
                go.Bar(x=variants['Sample'], y=variants['SNVs'],
                      name='SNVs', marker_color='red', opacity=0.7),
                row=1, col=2
            )
            fig.add_trace(
                go.Bar(x=variants['Sample'], y=variants['Indels'],
                      name='Indels', marker_color='blue', opacity=0.7),
                row=1, col=2
            )

        # Update layout
        fig.update_layout(
            title_text="Variant Analysis Dashboard",
            showlegend=True,
            height=800
        )

        # Save interactive plot
        pyo.plot(fig, filename=str(self.interactive_dir / 'variant_analysis_interactive.html'),
                auto_open=False)

    def create_summary_dashboard(self):
        """Create comprehensive summary dashboard"""
        self.log("Creating summary dashboard...")

        if 'summary' in self.data:
            summary = self.data['summary']
        else:
            # Create summary from available data
            summary = self.create_summary_statistics()

        # Create multi-panel summary figure
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # Data type overview
        data_types = []
        counts = []

        if 'mrna' in self.data:
            significant_mrna = len(self.data['mrna'][
                (self.data['mrna']['padj'] < 0.05) &
                (abs(self.data['mrna']['log2FoldChange']) > 1)
            ])
            data_types.extend(['mRNA Up', 'mRNA Down'])
            mrna_up = len(self.data['mrna'][
                (self.data['mrna']['padj'] < 0.05) &
                (self.data['mrna']['log2FoldChange'] > 1)
            ])
            mrna_down = len(self.data['mrna'][
                (self.data['mrna']['padj'] < 0.05) &
                (self.data['mrna']['log2FoldChange'] < -1)
            ])
            counts.extend([mrna_up, mrna_down])

        if 'variants' in self.data and 'PASS_Variants' in self.data['variants'].columns:
            total_variants = self.data['variants']['PASS_Variants'].sum()
            data_types.append('Variants')
            counts.append(total_variants)

        if 'mirna' in self.data:
            mirna_up = len(self.data['mirna'][
                (self.data['mirna']['padj'] < 0.05) &
                (self.data['mirna']['log2FoldChange'] > 1)
            ])
            mirna_down = len(self.data['mirna'][
                (self.data['mirna']['padj'] < 0.05) &
                (self.data['mirna']['log2FoldChange'] < -1)
            ])
            data_types.extend(['miRNA Up', 'miRNA Down'])
            counts.extend([mirna_up, mirna_down])

        if data_types:
            colors = plt.cm.Set3(np.linspace(0, 1, len(data_types)))
            axes[0, 0].bar(data_types, counts, color=colors)
            axes[0, 0].set_title('Multi-Modal Analysis Overview')
            axes[0, 0].set_ylabel('Count')
            axes[0, 0].tick_params(axis='x', rotation=45)

        # Additional summary plots can be added here
        for i in range(2):
            for j in range(3):
                if (i, j) != (0, 0):
                    axes[i, j].text(0.5, 0.5, 'Additional\nAnalysis\nPlaceholder',
                                   ha='center', va='center', transform=axes[i, j].transAxes,
                                   fontsize=12, alpha=0.5)
                    axes[i, j].set_xticks([])
                    axes[i, j].set_yticks([])

        plt.tight_layout()
        plt.savefig(self.plots_dir / f'summary_dashboard.{self.config["figure_format"]}',
                   dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()

    def create_summary_statistics(self):
        """Create summary statistics from available data"""
        summary = {}

        if 'mrna' in self.data:
            mrna = self.data['mrna']
            sig_mrna = mrna[(mrna['padj'] < 0.05) & (abs(mrna['log2FoldChange']) > 1)]
            summary['mrna_total'] = len(mrna)
            summary['mrna_significant'] = len(sig_mrna)
            summary['mrna_up'] = len(sig_mrna[sig_mrna['log2FoldChange'] > 0])
            summary['mrna_down'] = len(sig_mrna[sig_mrna['log2FoldChange'] < 0])

        if 'variants' in self.data:
            variants = self.data['variants']
            if 'PASS_Variants' in variants.columns:
                summary['total_variants'] = variants['PASS_Variants'].sum()
            if 'SNVs' in variants.columns:
                summary['total_snvs'] = variants['SNVs'].sum()
            if 'Indels' in variants.columns:
                summary['total_indels'] = variants['Indels'].sum()

        if 'mirna' in self.data:
            mirna = self.data['mirna']
            sig_mirna = mirna[(mirna['padj'] < 0.05) & (abs(mirna['log2FoldChange']) > 1)]
            summary['mirna_total'] = len(mirna)
            summary['mirna_significant'] = len(sig_mirna)
            summary['mirna_up'] = len(sig_mirna[sig_mirna['log2FoldChange'] > 0])
            summary['mirna_down'] = len(sig_mirna[sig_mirna['log2FoldChange'] < 0])

        return pd.DataFrame([summary])

    def generate_report(self):
        """Generate comprehensive HTML report"""
        self.log("Generating comprehensive visualization report...")

        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Multi-Modal Cancer Genomics Visualization Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
                .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 30px; }}
                .section {{ margin: 30px 0; padding: 20px; border: 1px solid #dee2e6; border-radius: 5px; }}
                .plot-container {{ text-align: center; margin: 20px 0; }}
                .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
                h1 {{ color: #343a40; }}
                h2 {{ color: #495057; border-bottom: 2px solid #dee2e6; padding-bottom: 10px; }}
                h3 {{ color: #6c757d; }}
                .stats-table {{ width: 100%; border-collapse: collapse; }}
                .stats-table th, .stats-table td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
                .stats-table th {{ background-color: #f8f9fa; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Multi-Modal Cancer Genomics Visualization Report</h1>
                <p>Comprehensive visualization suite for integrated cancer genomics analysis</p>
                <p><strong>Generated:</strong> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>

            <div class="section">
                <h2>Analysis Overview</h2>
                <p>This report presents visualizations from multi-modal cancer genomics analysis including:</p>
                <ul>
                    <li>mRNA differential expression analysis</li>
                    <li>Somatic variant discovery and characterization</li>
                    <li>miRNA expression profiling</li>
                    <li>Integrated multi-modal analysis</li>
                </ul>
            </div>

            <div class="section">
                <h2>Expression Analysis</h2>
                <div class="plot-grid">
        """

        # Add plot references
        plot_files = list(self.plots_dir.glob(f'*.{self.config["figure_format"]}'))
        for plot_file in sorted(plot_files):
            if 'volcano' in plot_file.name or 'ma' in plot_file.name or 'heatmap' in plot_file.name:
                html_content += f"""
                    <div class="plot-container">
                        <h3>{plot_file.stem.replace('_', ' ').title()}</h3>
                        <img src="plots/{plot_file.name}" alt="{plot_file.stem}" style="max-width: 100%; height: auto;">
                    </div>
                """

        html_content += """
                </div>
            </div>

            <div class="section">
                <h2>Variant Analysis</h2>
                <div class="plot-grid">
        """

        # Add variant plots
        for plot_file in sorted(plot_files):
            if 'variant' in plot_file.name:
                html_content += f"""
                    <div class="plot-container">
                        <h3>{plot_file.stem.replace('_', ' ').title()}</h3>
                        <img src="plots/{plot_file.name}" alt="{plot_file.stem}" style="max-width: 100%; height: auto;">
                    </div>
                """

        html_content += """
                </div>
            </div>

            <div class="section">
                <h2>Interactive Visualizations</h2>
                <p>Interactive plots are available in the interactive/ directory:</p>
                <ul>
        """

        # Add interactive plot links
        interactive_files = list(self.interactive_dir.glob('*.html'))
        for interactive_file in sorted(interactive_files):
            html_content += f'<li><a href="interactive/{interactive_file.name}">{interactive_file.stem.replace("_", " ").title()}</a></li>\n'

        html_content += """
                </ul>
            </div>

            <div class="section">
                <h2>Summary Statistics</h2>
        """

        # Add summary statistics if available
        if 'summary' in self.data:
            summary_df = self.data['summary']
            html_content += summary_df.to_html(classes='stats-table', table_id='summary-table')

        html_content += """
            </div>

            <div class="section">
                <h2>Files Generated</h2>
                <ul>
        """

        # List all generated files
        all_files = list(self.plots_dir.glob('*')) + list(self.interactive_dir.glob('*'))
        for file_path in sorted(all_files):
            html_content += f'<li>{file_path.relative_to(self.output_dir)}</li>\n'

        html_content += """
                </ul>
            </div>
        </body>
        </html>
        """

        # Save report
        report_file = self.reports_dir / 'visualization_report.html'
        with open(report_file, 'w') as f:
            f.write(html_content)

        self.log(f"Visualization report saved: {report_file}")

    def run_complete_analysis(self):
        """Run the complete visualization analysis"""
        self.log("Starting comprehensive visualization analysis...")

        # Load data
        if not self.load_data():
            self.log("No data files found. Please run the analysis pipeline first.")
            return False

        # Create all visualizations
        self.create_expression_plots()
        self.create_variant_plots()
        self.create_summary_dashboard()

        # Generate report
        self.generate_report()

        self.log("Visualization analysis completed successfully!")
        self.log(f"Results saved in: {self.output_dir}")

        return True


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Generate comprehensive visualizations for multi-modal cancer genomics analysis')
    parser.add_argument('--config', help='Configuration file path')
    parser.add_argument('--output-dir', default='10_visualizations', help='Output directory')

    args = parser.parse_args()

    # Create visualizer
    visualizer = MultiModalVisualizer(args.config)

    if args.output_dir:
        visualizer.config['output_dir'] = args.output_dir
        visualizer.setup_directories()

    # Run analysis
    success = visualizer.run_complete_analysis()

    if success:
        print(f"\n✓ Visualization analysis completed successfully!")
        print(f"✓ Results available in: {visualizer.output_dir}")
        print(f"✓ Open {visualizer.reports_dir}/visualization_report.html to view the complete report")
    else:
        print("✗ Visualization analysis failed. Check the logs for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()