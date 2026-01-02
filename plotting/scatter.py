"""
Fold change vs fold change scatter plot for pairwise dataset comparison.
"""

import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Optional, List


def create_fc_scatter(
    datasets: dict,
    dataset_a: str,
    dataset_b: str,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    use_padj: bool = True,
    highlight_genes: Optional[List[str]] = None,
    show_diagonal: bool = True,
    dark_mode: bool = False
) -> go.Figure:
    """
    Create a scatter plot comparing log2FC between two datasets.

    Args:
        datasets: Dictionary of dataset name to DataFrame
        dataset_a: Name of first dataset (x-axis)
        dataset_b: Name of second dataset (y-axis)
        log2fc_threshold: Threshold for significance
        pvalue_threshold: P-value threshold
        use_padj: Use adjusted p-value
        highlight_genes: List of genes to highlight
        show_diagonal: Show y=x diagonal line
        dark_mode: Use dark theme

    Returns:
        Plotly Figure object
    """
    df_a = datasets[dataset_a]
    df_b = datasets[dataset_b]

    # Merge on gene
    merged = df_a[['Gene', 'log2FC', 'pvalue', 'padj']].merge(
        df_b[['Gene', 'log2FC', 'pvalue', 'padj']],
        on='Gene',
        suffixes=('_A', '_B')
    )

    if merged.empty:
        return _empty_scatter("No overlapping genes found", dark_mode)

    # Determine significance
    pval_col = 'padj' if use_padj else 'pvalue'

    merged['Sig_A'] = (
        (merged['log2FC_A'].abs() >= log2fc_threshold) &
        (merged[f'{pval_col}_A'] <= pvalue_threshold)
    )
    merged['Sig_B'] = (
        (merged['log2FC_B'].abs() >= log2fc_threshold) &
        (merged[f'{pval_col}_B'] <= pvalue_threshold)
    )

    # Classify concordance
    merged['Category'] = 'Not Significant'

    # Significant in both - concordant up
    mask_up_both = (
        merged['Sig_A'] & merged['Sig_B'] &
        (merged['log2FC_A'] > 0) & (merged['log2FC_B'] > 0)
    )
    merged.loc[mask_up_both, 'Category'] = 'Up in Both'

    # Significant in both - concordant down
    mask_down_both = (
        merged['Sig_A'] & merged['Sig_B'] &
        (merged['log2FC_A'] < 0) & (merged['log2FC_B'] < 0)
    )
    merged.loc[mask_down_both, 'Category'] = 'Down in Both'

    # Significant in both - discordant
    mask_discordant = (
        merged['Sig_A'] & merged['Sig_B'] &
        (merged['log2FC_A'] * merged['log2FC_B'] < 0)
    )
    merged.loc[mask_discordant, 'Category'] = 'Discordant'

    # Significant in one only
    merged.loc[merged['Sig_A'] & ~merged['Sig_B'], 'Category'] = f'Sig in {dataset_a} only'
    merged.loc[~merged['Sig_A'] & merged['Sig_B'], 'Category'] = f'Sig in {dataset_b} only'

    # Theme settings
    if dark_mode:
        colors = {
            'Not Significant': 'rgba(128, 128, 128, 0.3)',
            'Up in Both': '#FF6B6B',
            'Down in Both': '#4ECDC4',
            'Discordant': '#FFD93D',
            f'Sig in {dataset_a} only': '#9B59B6',
            f'Sig in {dataset_b} only': '#E67E22'
        }
        paper_bgcolor = '#1E1E1E'
        plot_bgcolor = '#2D2D2D'
        font_color = '#FFFFFF'
        gridcolor = '#404040'
        diagonal_color = 'rgba(255, 255, 255, 0.3)'
    else:
        colors = {
            'Not Significant': 'rgba(180, 180, 180, 0.3)',
            'Up in Both': '#E74C3C',
            'Down in Both': '#3498DB',
            'Discordant': '#F39C12',
            f'Sig in {dataset_a} only': '#9B59B6',
            f'Sig in {dataset_b} only': '#E67E22'
        }
        paper_bgcolor = '#FFFFFF'
        plot_bgcolor = '#FAFAFA'
        font_color = '#2C3E50'
        gridcolor = '#E0E0E0'
        diagonal_color = 'rgba(0, 0, 0, 0.2)'

    fig = go.Figure()

    # Add traces for each category
    categories_order = [
        'Not Significant',
        f'Sig in {dataset_a} only',
        f'Sig in {dataset_b} only',
        'Discordant',
        'Down in Both',
        'Up in Both'
    ]

    for category in categories_order:
        subset = merged[merged['Category'] == category]
        if len(subset) == 0:
            continue

        fig.add_trace(go.Scatter(
            x=subset['log2FC_A'],
            y=subset['log2FC_B'],
            mode='markers',
            name=f'{category} ({len(subset)})',
            marker=dict(
                color=colors.get(category, 'gray'),
                size=6 if category == 'Not Significant' else 8,
                opacity=0.5 if category == 'Not Significant' else 0.8
            ),
            text=subset['Gene'],
            hovertemplate=(
                '<b>%{text}</b><br>'
                f'{dataset_a} log2FC: %{{x:.3f}}<br>'
                f'{dataset_b} log2FC: %{{y:.3f}}<br>'
                '<extra></extra>'
            )
        ))

    # Add highlighted genes
    if highlight_genes:
        highlight_df = merged[merged['Gene'].str.upper().isin([g.upper() for g in highlight_genes])]
        if len(highlight_df) > 0:
            fig.add_trace(go.Scatter(
                x=highlight_df['log2FC_A'],
                y=highlight_df['log2FC_B'],
                mode='markers+text',
                name='Highlighted',
                marker=dict(
                    color='#FFD700',
                    size=14,
                    symbol='star',
                    line=dict(width=1, color='black')
                ),
                text=highlight_df['Gene'],
                textposition='top center',
                textfont=dict(size=10, color=font_color),
                hovertemplate=(
                    '<b>%{text}</b><br>'
                    f'{dataset_a}: %{{x:.3f}}<br>'
                    f'{dataset_b}: %{{y:.3f}}<br>'
                    '<extra></extra>'
                )
            ))

    # Add diagonal line
    if show_diagonal:
        axis_range = max(
            abs(merged['log2FC_A'].min()), abs(merged['log2FC_A'].max()),
            abs(merged['log2FC_B'].min()), abs(merged['log2FC_B'].max())
        ) * 1.1

        fig.add_trace(go.Scatter(
            x=[-axis_range, axis_range],
            y=[-axis_range, axis_range],
            mode='lines',
            name='y = x',
            line=dict(color=diagonal_color, dash='dash', width=1),
            hoverinfo='skip',
            showlegend=False
        ))

    # Add quadrant lines
    fig.add_hline(y=0, line_color=gridcolor, line_width=1)
    fig.add_vline(x=0, line_color=gridcolor, line_width=1)

    # Calculate correlation
    valid_data = merged.dropna(subset=['log2FC_A', 'log2FC_B'])
    if len(valid_data) > 2:
        correlation = valid_data['log2FC_A'].corr(valid_data['log2FC_B'])
        corr_text = f"r = {correlation:.3f}"
    else:
        corr_text = ""

    # Calculate counts for title
    up_both = (merged['Category'] == 'Up in Both').sum()
    down_both = (merged['Category'] == 'Down in Both').sum()
    discordant = (merged['Category'] == 'Discordant').sum()

    fig.update_layout(
        title=dict(
            text=(
                f"{dataset_a} vs {dataset_b}<br>"
                f"<sup>Up both: {up_both} | Down both: {down_both} | Discordant: {discordant} | {corr_text}</sup>"
            ),
            font=dict(size=16, color=font_color)
        ),
        xaxis=dict(
            title=dict(text=f"log2FC ({dataset_a})", font=dict(size=14, color=font_color)),
            tickfont=dict(color=font_color),
            gridcolor=gridcolor,
            zerolinecolor=gridcolor
        ),
        yaxis=dict(
            title=dict(text=f"log2FC ({dataset_b})", font=dict(size=14, color=font_color)),
            tickfont=dict(color=font_color),
            gridcolor=gridcolor,
            zerolinecolor=gridcolor,
            scaleanchor="x",
            scaleratio=1
        ),
        paper_bgcolor=paper_bgcolor,
        plot_bgcolor=plot_bgcolor,
        font=dict(color=font_color),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            font=dict(size=10, color=font_color)
        ),
        hovermode='closest',
        height=600,
        margin=dict(l=60, r=150, t=80, b=60)
    )

    return fig


def _empty_scatter(message: str, dark_mode: bool) -> go.Figure:
    """Create an empty scatter plot with a message."""
    bg_color = '#1E1E1E' if dark_mode else '#FFFFFF'
    font_color = '#FFFFFF' if dark_mode else '#2C3E50'

    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(size=16, color=font_color)
    )
    fig.update_layout(
        paper_bgcolor=bg_color,
        plot_bgcolor=bg_color,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        height=500
    )
    return fig


def get_concordance_stats(
    datasets: dict,
    dataset_a: str,
    dataset_b: str,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05
) -> dict:
    """
    Calculate concordance statistics between two datasets.

    Returns:
        Dictionary with concordance metrics
    """
    df_a = datasets[dataset_a]
    df_b = datasets[dataset_b]

    merged = df_a[['Gene', 'log2FC', 'padj']].merge(
        df_b[['Gene', 'log2FC', 'padj']],
        on='Gene',
        suffixes=('_A', '_B')
    )

    sig_a = (merged['log2FC_A'].abs() >= log2fc_threshold) & (merged['padj_A'] <= pvalue_threshold)
    sig_b = (merged['log2FC_B'].abs() >= log2fc_threshold) & (merged['padj_B'] <= pvalue_threshold)

    both_sig = sig_a & sig_b
    concordant = both_sig & (merged['log2FC_A'] * merged['log2FC_B'] > 0)
    discordant = both_sig & (merged['log2FC_A'] * merged['log2FC_B'] < 0)

    return {
        'total_genes': len(merged),
        'sig_in_a': sig_a.sum(),
        'sig_in_b': sig_b.sum(),
        'sig_in_both': both_sig.sum(),
        'concordant': concordant.sum(),
        'discordant': discordant.sum(),
        'concordance_rate': concordant.sum() / both_sig.sum() if both_sig.sum() > 0 else 0,
        'correlation': merged['log2FC_A'].corr(merged['log2FC_B'])
    }
