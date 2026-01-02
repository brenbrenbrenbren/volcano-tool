"""
Bar plot visualizations for gene and pathway analysis.
"""

import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Set


def create_gene_barplot(
    datasets: Dict[str, pd.DataFrame],
    gene_name: str,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    use_padj: bool = True,
    dark_mode: bool = False
) -> go.Figure:
    """
    Create a bar plot showing a single gene's log2FC across all datasets.

    Args:
        datasets: Dictionary of dataset name to DataFrame
        gene_name: Gene to plot
        log2fc_threshold: Threshold for significance coloring
        pvalue_threshold: P-value threshold for significance
        use_padj: Use adjusted p-value
        dark_mode: Use dark theme

    Returns:
        Plotly Figure object
    """
    gene_name_upper = gene_name.upper().strip()

    # Collect data
    data = []
    for name, df in datasets.items():
        gene_row = df[df['Gene'].str.upper() == gene_name_upper]

        if len(gene_row) > 0:
            row = gene_row.iloc[0]
            log2fc = row['log2FC']
            pval = row.get('padj' if use_padj else 'pvalue', row.get('pvalue', 1))

            is_sig = abs(log2fc) >= log2fc_threshold and pval <= pvalue_threshold

            data.append({
                'Dataset': name,
                'log2FC': log2fc,
                'pvalue': pval,
                'Significant': is_sig,
                'Direction': 'Up' if log2fc > 0 else 'Down' if log2fc < 0 else 'None'
            })
        else:
            data.append({
                'Dataset': name,
                'log2FC': 0,
                'pvalue': 1,
                'Significant': False,
                'Direction': 'None'
            })

    df_plot = pd.DataFrame(data)

    # Theme settings
    if dark_mode:
        color_up = '#FF6B6B'
        color_down = '#4ECDC4'
        color_ns = 'rgba(128, 128, 128, 0.5)'
        paper_bgcolor = '#1E1E1E'
        plot_bgcolor = '#2D2D2D'
        font_color = '#FFFFFF'
        gridcolor = '#404040'
    else:
        color_up = '#E74C3C'
        color_down = '#3498DB'
        color_ns = 'rgba(180, 180, 180, 0.7)'
        paper_bgcolor = '#FFFFFF'
        plot_bgcolor = '#FAFAFA'
        font_color = '#2C3E50'
        gridcolor = '#E0E0E0'

    # Assign colors
    colors = []
    for _, row in df_plot.iterrows():
        if not row['Significant']:
            colors.append(color_ns)
        elif row['log2FC'] > 0:
            colors.append(color_up)
        else:
            colors.append(color_down)

    fig = go.Figure(data=go.Bar(
        x=df_plot['Dataset'],
        y=df_plot['log2FC'],
        marker_color=colors,
        text=[f"{fc:.2f}" for fc in df_plot['log2FC']],
        textposition='outside',
        textfont=dict(size=10, color=font_color),
        hovertemplate=(
            '<b>%{x}</b><br>'
            'log2FC: %{y:.3f}<br>'
            '<extra></extra>'
        )
    ))

    # Add significance markers
    for i, row in df_plot.iterrows():
        if row['Significant']:
            marker_y = row['log2FC'] + (0.3 if row['log2FC'] >= 0 else -0.3)
            fig.add_annotation(
                x=row['Dataset'],
                y=marker_y,
                text="*",
                showarrow=False,
                font=dict(size=16, color=font_color)
            )

    # Add threshold lines
    fig.add_hline(y=log2fc_threshold, line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_hline(y=-log2fc_threshold, line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_hline(y=0, line_color=gridcolor, line_width=1)

    # Count significant datasets
    sig_count = df_plot['Significant'].sum()

    fig.update_layout(
        title=dict(
            text=f"<b>{gene_name}</b> across datasets<br><sup>Significant in {sig_count}/{len(datasets)} datasets (* p < {pvalue_threshold})</sup>",
            font=dict(size=16, color=font_color)
        ),
        xaxis=dict(
            title="Dataset",
            tickfont=dict(color=font_color),
            tickangle=45
        ),
        yaxis=dict(
            title=dict(text="log2(Fold Change)", font=dict(color=font_color)),
            tickfont=dict(color=font_color),
            gridcolor=gridcolor
        ),
        paper_bgcolor=paper_bgcolor,
        plot_bgcolor=plot_bgcolor,
        font=dict(color=font_color),
        height=450,
        margin=dict(l=60, r=40, t=80, b=100),
        showlegend=False
    )

    return fig


def create_pathway_barplot(
    datasets: Dict[str, pd.DataFrame],
    pathway_genes: Set[str],
    pathway_name: str,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    dark_mode: bool = False,
    show_all_genes: bool = False,
    max_genes: int = 30
) -> go.Figure:
    """
    Create a grouped bar plot for pathway genes across datasets.

    Args:
        datasets: Dictionary of dataset name to DataFrame
        pathway_genes: Set of gene symbols in the pathway
        pathway_name: Name of the pathway
        log2fc_threshold: Threshold for significance
        pvalue_threshold: P-value threshold
        dark_mode: Use dark theme
        show_all_genes: Show all genes or only significant ones
        max_genes: Maximum genes to display

    Returns:
        Plotly Figure object
    """
    pathway_genes_upper = {g.upper() for g in pathway_genes}

    # Collect data for all pathway genes across datasets
    all_data = []

    for gene in pathway_genes_upper:
        for name, df in datasets.items():
            gene_row = df[df['Gene'].str.upper() == gene]

            if len(gene_row) > 0:
                row = gene_row.iloc[0]
                log2fc = row['log2FC']
                pval = row.get('padj', row.get('pvalue', 1))
                is_sig = abs(log2fc) >= log2fc_threshold and pval <= pvalue_threshold

                all_data.append({
                    'Gene': row['Gene'],
                    'Dataset': name,
                    'log2FC': log2fc,
                    'pvalue': pval,
                    'Significant': is_sig
                })

    if not all_data:
        return _empty_bar(f"No pathway genes found in datasets", dark_mode)

    df_all = pd.DataFrame(all_data)

    # Filter genes
    if not show_all_genes:
        sig_genes = df_all[df_all['Significant']]['Gene'].unique()
        df_all = df_all[df_all['Gene'].isin(sig_genes)]

    if df_all.empty:
        return _empty_bar(f"No significant pathway genes found", dark_mode)

    # Limit genes
    unique_genes = df_all['Gene'].unique()
    if len(unique_genes) > max_genes:
        # Prioritize by max absolute fold change
        gene_max_fc = df_all.groupby('Gene')['log2FC'].apply(lambda x: x.abs().max())
        top_genes = gene_max_fc.nlargest(max_genes).index
        df_all = df_all[df_all['Gene'].isin(top_genes)]

    # Theme settings
    if dark_mode:
        paper_bgcolor = '#1E1E1E'
        plot_bgcolor = '#2D2D2D'
        font_color = '#FFFFFF'
        gridcolor = '#404040'
        base_colors = ['#FF6B6B', '#4ECDC4', '#FFD93D', '#9B59B6', '#E67E22',
                       '#1ABC9C', '#3498DB', '#E91E63', '#00BCD4', '#8BC34A']
    else:
        paper_bgcolor = '#FFFFFF'
        plot_bgcolor = '#FAFAFA'
        font_color = '#2C3E50'
        gridcolor = '#E0E0E0'
        base_colors = ['#E74C3C', '#3498DB', '#F39C12', '#9B59B6', '#E67E22',
                       '#1ABC9C', '#2980B9', '#E91E63', '#00BCD4', '#8BC34A']

    dataset_names = list(datasets.keys())
    color_map = {name: base_colors[i % len(base_colors)] for i, name in enumerate(dataset_names)}

    fig = go.Figure()

    # Sort genes by average fold change
    gene_order = df_all.groupby('Gene')['log2FC'].mean().sort_values(ascending=False).index.tolist()

    for dataset_name in dataset_names:
        dataset_data = df_all[df_all['Dataset'] == dataset_name].set_index('Gene')

        y_values = []
        for gene in gene_order:
            if gene in dataset_data.index:
                y_values.append(dataset_data.loc[gene, 'log2FC'])
            else:
                y_values.append(0)

        fig.add_trace(go.Bar(
            name=dataset_name,
            x=gene_order,
            y=y_values,
            marker_color=color_map[dataset_name],
            hovertemplate=(
                f'<b>{dataset_name}</b><br>'
                'Gene: %{x}<br>'
                'log2FC: %{y:.3f}<br>'
                '<extra></extra>'
            )
        ))

    # Add threshold lines
    fig.add_hline(y=log2fc_threshold, line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_hline(y=-log2fc_threshold, line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_hline(y=0, line_color=gridcolor, line_width=1)

    # Calculate pathway statistics
    total_pathway_genes = len(pathway_genes_upper)
    found_genes = df_all['Gene'].nunique()
    sig_genes = df_all[df_all['Significant']]['Gene'].nunique()

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{pathway_name}</b><br>"
                f"<sup>{found_genes}/{total_pathway_genes} genes found | {sig_genes} significant</sup>"
            ),
            font=dict(size=16, color=font_color)
        ),
        xaxis=dict(
            title="Gene",
            tickfont=dict(color=font_color, size=9),
            tickangle=45
        ),
        yaxis=dict(
            title=dict(text="log2(Fold Change)", font=dict(color=font_color)),
            tickfont=dict(color=font_color),
            gridcolor=gridcolor
        ),
        paper_bgcolor=paper_bgcolor,
        plot_bgcolor=plot_bgcolor,
        font=dict(color=font_color),
        barmode='group',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(size=10, color=font_color)
        ),
        height=500,
        margin=dict(l=60, r=40, t=100, b=120)
    )

    return fig


def _empty_bar(message: str, dark_mode: bool) -> go.Figure:
    """Create an empty bar plot with a message."""
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
        height=400
    )
    return fig
