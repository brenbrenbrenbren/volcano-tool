"""
Export utilities for figures and data.
"""

import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from typing import Optional, Dict, Any
import io


class Exporter:
    """Handle export of figures and data."""

    def __init__(self, dpi: int = 300):
        self.dpi = dpi

    def export_figure(
        self,
        fig: go.Figure,
        filename: str,
        format: str = 'png',
        width: int = 1200,
        height: int = 800
    ) -> bytes:
        """
        Export Plotly figure to image bytes.

        Args:
            fig: Plotly figure
            filename: Base filename (without extension)
            format: 'png' or 'pdf'
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            Image bytes
        """
        scale = self.dpi / 96  # Plotly default is 96 DPI

        return fig.to_image(
            format=format,
            width=width,
            height=height,
            scale=scale
        )

    def export_dataframe(
        self,
        df: pd.DataFrame,
        filename: str,
        format: str = 'csv'
    ) -> bytes:
        """
        Export DataFrame to file bytes.

        Args:
            df: DataFrame to export
            filename: Base filename
            format: 'csv' or 'xlsx'

        Returns:
            File bytes
        """
        buffer = io.BytesIO()

        if format == 'csv':
            df.to_csv(buffer, index=False)
        elif format == 'xlsx':
            df.to_excel(buffer, index=False, engine='openpyxl')

        buffer.seek(0)
        return buffer.getvalue()

    def create_deg_report(
        self,
        datasets: Dict[str, pd.DataFrame],
        thresholds: Dict[str, float]
    ) -> pd.DataFrame:
        """
        Create a summary report of DEGs.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            thresholds: Dictionary with 'log2fc' and 'pvalue' thresholds

        Returns:
            Summary DataFrame
        """
        log2fc_thresh = thresholds.get('log2fc', 1.0)
        pval_thresh = thresholds.get('pvalue', 0.05)

        summary_rows = []

        for name, df in datasets.items():
            total_genes = len(df)

            sig_mask = (df['log2FC'].abs() >= log2fc_thresh) & (df['padj'] <= pval_thresh)
            up_mask = (df['log2FC'] >= log2fc_thresh) & (df['padj'] <= pval_thresh)
            down_mask = (df['log2FC'] <= -log2fc_thresh) & (df['padj'] <= pval_thresh)

            summary_rows.append({
                'Dataset': name,
                'Total Genes': total_genes,
                'Significant DEGs': sig_mask.sum(),
                'Upregulated': up_mask.sum(),
                'Downregulated': down_mask.sum(),
                'Top Upregulated': ', '.join(
                    df[up_mask].nlargest(5, 'log2FC')['Gene'].tolist()
                ),
                'Top Downregulated': ', '.join(
                    df[down_mask].nsmallest(5, 'log2FC')['Gene'].tolist()
                )
            })

        return pd.DataFrame(summary_rows)

    def create_full_export(
        self,
        datasets: Dict[str, pd.DataFrame],
        thresholds: Dict[str, float]
    ) -> pd.DataFrame:
        """
        Create comprehensive export with all genes across all datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            thresholds: Dictionary with thresholds

        Returns:
            Combined DataFrame
        """
        log2fc_thresh = thresholds.get('log2fc', 1.0)
        pval_thresh = thresholds.get('pvalue', 0.05)

        # Get all unique genes
        all_genes = set()
        for df in datasets.values():
            all_genes.update(df['Gene'].unique())

        # Build export table
        export_data = []

        for gene in sorted(all_genes):
            row = {'Gene': gene}
            sig_count = 0

            for name, df in datasets.items():
                gene_data = df[df['Gene'] == gene]

                if len(gene_data) > 0:
                    g = gene_data.iloc[0]
                    row[f'{name}_log2FC'] = g['log2FC']
                    row[f'{name}_padj'] = g.get('padj', g.get('pvalue', 1))

                    is_sig = abs(g['log2FC']) >= log2fc_thresh and row[f'{name}_padj'] <= pval_thresh
                    row[f'{name}_Significant'] = is_sig
                    if is_sig:
                        sig_count += 1
                else:
                    row[f'{name}_log2FC'] = None
                    row[f'{name}_padj'] = None
                    row[f'{name}_Significant'] = False

            row['Significant_In_N_Datasets'] = sig_count
            export_data.append(row)

        return pd.DataFrame(export_data)
