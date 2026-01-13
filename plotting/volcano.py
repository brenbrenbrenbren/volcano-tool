"""
Volcano plot visualization.
"""

import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Dict


def create_volcano_plot(
    df: pd.DataFrame,
    title: str = "Volcano Plot",
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    use_padj: bool = True,
    highlight_genes: Optional[List[str]] = None,
    dark_mode: bool = False,
    show_labels: bool = True,
    top_n_labels: int = 10,
    gene_index: Optional[Dict[str, int]] = None
) -> go.Figure:
    """
    Create an interactive volcano plot.

    Args:
        df: DataFrame with Gene, log2FC, pvalue, padj columns
        title: Plot title
        log2fc_threshold: Threshold for fold change significance
        pvalue_threshold: Threshold for p-value significance
        use_padj: Use adjusted p-value
        highlight_genes: List of genes to highlight
        dark_mode: Use dark theme
        show_labels: Show gene labels for top genes
        top_n_labels: Number of top genes to label
        gene_index: Optional dict of {gene_upper: row_index} for O(1) gene lookups

    Returns:
        Plotly Figure object
    """
    df = df.copy()
    pval_col = 'padj' if use_padj else 'pvalue'

    # Calculate -log10(pvalue)
    df['neg_log10_pval'] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Classify significance
    df['Significance'] = 'Not Significant'
    up_mask = (df['log2FC'] >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
    down_mask = (df['log2FC'] <= -log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
    df.loc[up_mask, 'Significance'] = 'Upregulated'
    df.loc[down_mask, 'Significance'] = 'Downregulated'

    # Colors and theme
    if dark_mode:
        colors = {
            'Not Significant': 'rgba(128, 128, 128, 0.4)',
            'Upregulated': '#FF6B6B',
            'Downregulated': '#4ECDC4'
        }
        template = 'plotly_dark'
        paper_bgcolor = '#1E1E1E'
        plot_bgcolor = '#2D2D2D'
        font_color = '#FFFFFF'
        gridcolor = '#404040'
    else:
        colors = {
            'Not Significant': 'rgba(180, 180, 180, 0.5)',
            'Upregulated': '#E74C3C',
            'Downregulated': '#3498DB'
        }
        template = 'plotly_white'
        paper_bgcolor = '#FFFFFF'
        plot_bgcolor = '#FAFAFA'
        font_color = '#2C3E50'
        gridcolor = '#E0E0E0'

    fig = go.Figure()

    # Add traces for each significance category
    for sig_type in ['Not Significant', 'Downregulated', 'Upregulated']:
        subset = df[df['Significance'] == sig_type]

        fig.add_trace(go.Scatter(
            x=subset['log2FC'],
            y=subset['neg_log10_pval'],
            mode='markers',
            name=sig_type,
            marker=dict(
                color=colors[sig_type],
                size=6 if sig_type == 'Not Significant' else 8,
                opacity=0.7 if sig_type == 'Not Significant' else 0.9,
                line=dict(width=0)
            ),
            text=subset['Gene'],
            hovertemplate=(
                '<b>%{text}</b><br>'
                'log2FC: %{x:.3f}<br>'
                f'-log10({pval_col}): %{{y:.3f}}<br>'
                '<extra></extra>'
            )
        ))

    # Highlight specific genes if provided
    if highlight_genes:
        # Use gene_index for O(1) lookups if available; otherwise fallback to O(n) scan
        if gene_index:
            # O(1) lookups using pre-built index
            highlight_indices = []
            for gene in highlight_genes:
                gene_upper = gene.upper().strip()
                if gene_upper in gene_index:
                    highlight_indices.append(gene_index[gene_upper])

            if highlight_indices:
                highlight_df = df.iloc[highlight_indices].copy()
            else:
                highlight_df = pd.DataFrame()
        else:
            # Fallback: O(n) scan if index not available
            highlight_df = df[df['Gene'].str.upper().isin([g.upper() for g in highlight_genes])]

        if len(highlight_df) > 0:
            fig.add_trace(go.Scatter(
                x=highlight_df['log2FC'],
                y=highlight_df['neg_log10_pval'],
                mode='markers+text',
                name='Highlighted',
                marker=dict(
                    color='#FFD700',
                    size=12,
                    symbol='star',
                    line=dict(width=1, color='black')
                ),
                text=highlight_df['Gene'],
                textposition='top center',
                textfont=dict(size=10, color=font_color),
                hovertemplate=(
                    '<b>%{text}</b><br>'
                    'log2FC: %{x:.3f}<br>'
                    f'-log10({pval_col}): %{{y:.3f}}<br>'
                    '<extra></extra>'
                )
            ))

    # Add labels for top significant genes
    if show_labels and top_n_labels > 0:
        sig_df = df[df['Significance'] != 'Not Significant'].copy()
        if len(sig_df) > 0:
            sig_df['abs_log2FC'] = sig_df['log2FC'].abs()
            top_genes = sig_df.nlargest(top_n_labels, 'neg_log10_pval')

            fig.add_trace(go.Scatter(
                x=top_genes['log2FC'],
                y=top_genes['neg_log10_pval'],
                mode='text',
                text=top_genes['Gene'],
                textposition='top center',
                textfont=dict(size=9, color=font_color),
                showlegend=False,
                hoverinfo='skip'
            ))

    # Add threshold lines
    pval_line = -np.log10(pvalue_threshold)

    fig.add_hline(
        y=pval_line,
        line_dash="dash",
        line_color="gray",
        opacity=0.7,
        annotation_text=f"p = {pvalue_threshold}",
        annotation_position="right"
    )

    fig.add_vline(
        x=log2fc_threshold,
        line_dash="dash",
        line_color="gray",
        opacity=0.7
    )

    fig.add_vline(
        x=-log2fc_threshold,
        line_dash="dash",
        line_color="gray",
        opacity=0.7
    )

    # Calculate counts for subtitle
    up_count = up_mask.sum()
    down_count = down_mask.sum()

    # Update layout
    fig.update_layout(
        title=dict(
            text=f"{title}<br><sup>{up_count} upregulated | {down_count} downregulated</sup>",
            font=dict(size=16, color=font_color)
        ),
        xaxis=dict(
            title=dict(text="log2(Fold Change)", font=dict(size=14, color=font_color)),
            tickfont=dict(color=font_color),
            gridcolor=gridcolor,
            zerolinecolor=gridcolor
        ),
        yaxis=dict(
            title=dict(text=f"-log10({'adjusted p-value' if use_padj else 'p-value'})", font=dict(size=14, color=font_color)),
            tickfont=dict(color=font_color),
            gridcolor=gridcolor,
            zerolinecolor=gridcolor
        ),
        template=template,
        paper_bgcolor=paper_bgcolor,
        plot_bgcolor=plot_bgcolor,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(color=font_color)
        ),
        hovermode='closest',
        margin=dict(l=60, r=40, t=80, b=60)
    )

    return fig


def create_multi_volcano(
    datasets: dict,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    use_padj: bool = True,
    dark_mode: bool = False,
    cols: int = 2
) -> go.Figure:
    """
    Create a grid of volcano plots for multiple datasets.

    Args:
        datasets: Dictionary of dataset name to DataFrame
        log2fc_threshold: Threshold for fold change
        pvalue_threshold: Threshold for p-value
        use_padj: Use adjusted p-value
        dark_mode: Use dark theme
        cols: Number of columns in grid

    Returns:
        Plotly Figure with subplots
    """
    from plotly.subplots import make_subplots
    import math

    n_datasets = len(datasets)
    rows = math.ceil(n_datasets / cols)

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=list(datasets.keys()),
        horizontal_spacing=0.08,
        vertical_spacing=0.12
    )

    # Theme settings
    if dark_mode:
        colors = {
            'NS': 'rgba(128, 128, 128, 0.3)',
            'Up': '#FF6B6B',
            'Down': '#4ECDC4'
        }
        paper_bgcolor = '#1E1E1E'
        plot_bgcolor = '#2D2D2D'
    else:
        colors = {
            'NS': 'rgba(180, 180, 180, 0.4)',
            'Up': '#E74C3C',
            'Down': '#3498DB'
        }
        paper_bgcolor = '#FFFFFF'
        plot_bgcolor = '#FAFAFA'

    pval_col = 'padj' if use_padj else 'pvalue'

    for idx, (name, df) in enumerate(datasets.items()):
        row = idx // cols + 1
        col = idx % cols + 1

        df = df.copy()
        df['neg_log10_pval'] = -np.log10(df[pval_col].clip(lower=1e-300))

        # Classify
        up_mask = (df['log2FC'] >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
        down_mask = (df['log2FC'] <= -log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
        ns_mask = ~up_mask & ~down_mask

        # Add traces
        for mask, color, sig_name in [
            (ns_mask, colors['NS'], 'NS'),
            (down_mask, colors['Down'], 'Down'),
            (up_mask, colors['Up'], 'Up')
        ]:
            subset = df[mask]
            fig.add_trace(
                go.Scatter(
                    x=subset['log2FC'],
                    y=subset['neg_log10_pval'],
                    mode='markers',
                    marker=dict(color=color, size=4, opacity=0.7),
                    text=subset['Gene'],
                    hovertemplate='<b>%{text}</b><br>log2FC: %{x:.2f}<extra></extra>',
                    showlegend=(idx == 0)
                ),
                row=row,
                col=col
            )

    fig.update_layout(
        paper_bgcolor=paper_bgcolor,
        plot_bgcolor=plot_bgcolor,
        height=300 * rows,
        showlegend=False,
        margin=dict(l=50, r=30, t=50, b=50)
    )

    return fig
