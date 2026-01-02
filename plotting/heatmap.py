"""
Heatmap visualization for cross-dataset comparison.
"""

import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Dict, List, Optional
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


def create_heatmap(
    datasets: Dict[str, pd.DataFrame],
    genes: Optional[List[str]] = None,
    cluster_genes: bool = True,
    cluster_datasets: bool = False,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    filter_significant: bool = True,
    min_datasets: int = 1,
    dark_mode: bool = False,
    max_genes: int = 100,
    show_values: bool = False,
    normalize_rows: bool = False
) -> go.Figure:
    """
    Create a heatmap of log2FC values across datasets.

    Args:
        datasets: Dictionary of dataset name to DataFrame
        genes: Specific genes to include (if None, auto-select)
        cluster_genes: Hierarchically cluster genes
        cluster_datasets: Hierarchically cluster datasets
        log2fc_threshold: Threshold for significance
        pvalue_threshold: P-value threshold
        filter_significant: Only show significant genes
        min_datasets: Minimum datasets a gene must be significant in
        dark_mode: Use dark theme
        max_genes: Maximum number of genes to display
        show_values: Show FC values as text on cells
        normalize_rows: Z-score normalize each row (gene) for better comparison

    Returns:
        Plotly Figure object
    """
    # Build matrix of log2FC values
    matrix_data = _build_fc_matrix(datasets)

    if matrix_data.empty:
        return _empty_heatmap("No data available", dark_mode)

    # Filter genes if requested
    if genes is not None:
        genes_upper = [g.upper() for g in genes]
        matrix_data = matrix_data[matrix_data.index.str.upper().isin(genes_upper)]
    elif filter_significant:
        matrix_data = _filter_significant_genes(
            matrix_data,
            datasets,
            log2fc_threshold,
            pvalue_threshold,
            min_datasets
        )

    if matrix_data.empty:
        return _empty_heatmap("No significant genes found with current filters", dark_mode)

    # Limit genes if too many
    if len(matrix_data) > max_genes:
        # Sort by variance and take top genes
        variance = matrix_data.var(axis=1)
        top_genes = variance.nlargest(max_genes).index
        matrix_data = matrix_data.loc[top_genes]

    # Store original values for display before normalization
    original_matrix = matrix_data.copy()

    # Per-row normalization (z-score) if requested
    if normalize_rows:
        matrix_data = _normalize_rows(matrix_data)

    # Clustering
    if cluster_genes and len(matrix_data) > 2:
        order = _get_cluster_order(matrix_data, axis=0)
        matrix_data = matrix_data.iloc[order]
        original_matrix = original_matrix.iloc[order]

    if cluster_datasets and len(matrix_data.columns) > 2:
        order = _get_cluster_order(matrix_data, axis=1)
        matrix_data = matrix_data.iloc[:, order]
        original_matrix = original_matrix.iloc[:, order]

    # Theme settings
    if dark_mode:
        paper_bgcolor = '#1E1E1E'
        plot_bgcolor = '#2D2D2D'
        font_color = '#FFFFFF'
        colorscale = [
            [0, '#4ECDC4'],      # Down - cyan
            [0.5, '#2D2D2D'],    # Neutral - dark gray
            [1, '#FF6B6B']       # Up - red
        ]
        text_color = '#FFFFFF'
    else:
        paper_bgcolor = '#FFFFFF'
        plot_bgcolor = '#FAFAFA'
        font_color = '#2C3E50'
        colorscale = [
            [0, '#3498DB'],      # Down - blue
            [0.5, '#FFFFFF'],    # Neutral - white
            [1, '#E74C3C']       # Up - red
        ]
        text_color = '#2C3E50'

    # Determine color range (symmetric around 0)
    if normalize_rows:
        max_abs_val = 3  # Z-score typically ranges -3 to 3
        colorbar_title = "Z-score"
    else:
        max_abs_val = max(abs(matrix_data.values.min()), abs(matrix_data.values.max()))
        max_abs_val = min(max_abs_val, 5)  # Cap at 5 for visualization
        colorbar_title = "log2FC"

    # Build hover text with original values
    hover_text = []
    for i, gene in enumerate(matrix_data.index):
        row_text = []
        for j, dataset in enumerate(matrix_data.columns):
            orig_val = original_matrix.iloc[i, j]
            display_val = matrix_data.iloc[i, j]
            if pd.isna(orig_val):
                row_text.append(f"<b>{gene}</b><br>Dataset: {dataset}<br>log2FC: N/A")
            else:
                if normalize_rows:
                    row_text.append(
                        f"<b>{gene}</b><br>Dataset: {dataset}<br>"
                        f"log2FC: {orig_val:.3f}<br>Z-score: {display_val:.2f}"
                    )
                else:
                    row_text.append(
                        f"<b>{gene}</b><br>Dataset: {dataset}<br>log2FC: {orig_val:.3f}"
                    )
        hover_text.append(row_text)

    # Build text annotations for cells
    if show_values:
        text_values = []
        for i in range(len(original_matrix)):
            row_text = []
            for j in range(len(original_matrix.columns)):
                val = original_matrix.iloc[i, j]
                if pd.isna(val):
                    row_text.append("")
                else:
                    row_text.append(f"{val:.1f}")
            text_values.append(row_text)
    else:
        text_values = None

    # Create heatmap
    heatmap_trace = go.Heatmap(
        z=matrix_data.values,
        x=matrix_data.columns.tolist(),
        y=matrix_data.index.tolist(),
        colorscale=colorscale,
        zmin=-max_abs_val,
        zmax=max_abs_val,
        colorbar=dict(
            title=dict(text=colorbar_title, font=dict(color=font_color)),
            tickfont=dict(color=font_color)
        ),
        hovertemplate='%{customdata}<extra></extra>',
        customdata=hover_text,
    )

    # Add text if requested
    if show_values:
        heatmap_trace.text = text_values
        heatmap_trace.texttemplate = "%{text}"
        heatmap_trace.textfont = dict(size=8, color=text_color)

    fig = go.Figure(data=heatmap_trace)

    # Build title
    title_text = f"<b>Expression Heatmap</b> ({len(matrix_data)} genes × {len(matrix_data.columns)} datasets)"
    if normalize_rows:
        title_text += "<br><sup>Row-normalized (Z-score) for cross-assay comparison</sup>"

    # Update layout
    fig.update_layout(
        title=dict(
            text=title_text,
            font=dict(size=16, color=font_color),
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title=dict(text="Dataset", font=dict(color=font_color)),
            tickfont=dict(color=font_color, size=10),
            tickangle=45
        ),
        yaxis=dict(
            title=dict(text="Gene", font=dict(color=font_color)),
            tickfont=dict(color=font_color, size=9),
            autorange="reversed"
        ),
        paper_bgcolor=paper_bgcolor,
        plot_bgcolor=plot_bgcolor,
        font=dict(color=font_color),
        height=max(400, len(matrix_data) * 18 + 150),
        margin=dict(l=120, r=80, t=100, b=120)
    )

    return fig


def _build_fc_matrix(datasets: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build a matrix of log2FC values with genes as rows, datasets as columns."""
    all_data = {}

    for name, df in datasets.items():
        gene_fc = df.set_index('Gene')['log2FC'].to_dict()
        all_data[name] = gene_fc

    matrix = pd.DataFrame(all_data)
    matrix = matrix.dropna(how='all')  # Remove genes with no data

    return matrix


def _filter_significant_genes(
    matrix: pd.DataFrame,
    datasets: Dict[str, pd.DataFrame],
    log2fc_threshold: float,
    pvalue_threshold: float,
    min_datasets: int
) -> pd.DataFrame:
    """Filter to genes significant in at least min_datasets."""
    sig_counts = {}

    for gene in matrix.index:
        count = 0
        for name, df in datasets.items():
            gene_data = df[df['Gene'] == gene]
            if len(gene_data) > 0:
                row = gene_data.iloc[0]
                pval = row.get('padj', row.get('pvalue', 1))
                fc = abs(row.get('log2FC', 0))
                if fc >= log2fc_threshold and pval <= pvalue_threshold:
                    count += 1
        sig_counts[gene] = count

    sig_genes = [gene for gene, count in sig_counts.items() if count >= min_datasets]

    return matrix.loc[matrix.index.isin(sig_genes)]


def _normalize_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    """Z-score normalize each row (gene) for better cross-assay comparison."""
    matrix_norm = matrix.copy()

    for idx in matrix_norm.index:
        row = matrix_norm.loc[idx]
        valid_vals = row.dropna()
        if len(valid_vals) > 1:
            mean_val = valid_vals.mean()
            std_val = valid_vals.std()
            if std_val > 0:
                matrix_norm.loc[idx] = (row - mean_val) / std_val

    return matrix_norm


def _get_cluster_order(matrix: pd.DataFrame, axis: int = 0) -> List[int]:
    """Get hierarchical clustering order for rows (axis=0) or columns (axis=1)."""
    try:
        # Fill NaN with 0 for clustering
        matrix_filled = matrix.fillna(0)

        if axis == 0:  # Cluster rows
            if len(matrix_filled) < 2:
                return list(range(len(matrix_filled)))
            distances = pdist(matrix_filled.values, metric='euclidean')
        else:  # Cluster columns
            if len(matrix_filled.columns) < 2:
                return list(range(len(matrix_filled.columns)))
            distances = pdist(matrix_filled.T.values, metric='euclidean')

        linkage = hierarchy.linkage(distances, method='average')
        order = hierarchy.leaves_list(linkage)

        return list(order)
    except Exception:
        if axis == 0:
            return list(range(len(matrix)))
        return list(range(len(matrix.columns)))


def _empty_heatmap(message: str, dark_mode: bool) -> go.Figure:
    """Create an empty heatmap with a message."""
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
