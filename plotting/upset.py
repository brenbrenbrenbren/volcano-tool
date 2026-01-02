"""
UpSet plot visualization for multi-dataset overlap analysis.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from typing import Dict, Set, List, Tuple, Optional


def create_upset_plot(
    deg_sets: Dict[str, Set[str]],
    min_subset_size: int = 1,
    max_subsets: int = 40,
    dark_mode: bool = False,
    show_percentages: bool = True,
    show_gene_labels_threshold: int = 10
) -> Tuple[go.Figure, List[dict]]:
    """
    Create an UpSet plot for visualizing set intersections.

    Args:
        deg_sets: Dictionary mapping dataset names to sets of gene names
        min_subset_size: Minimum intersection size to display
        max_subsets: Maximum number of intersections to show
        dark_mode: Use dark theme
        show_percentages: Show percentage labels
        show_gene_labels_threshold: Show gene names for intersections <= this size

    Returns:
        Tuple of (Plotly Figure object, intersection data list)
    """
    if not deg_sets:
        return _empty_figure("No data provided", dark_mode), []

    # Calculate all intersections
    intersections = _calculate_intersections(deg_sets)

    # Filter by minimum size
    intersections = {k: v for k, v in intersections.items() if len(v) >= min_subset_size}

    if not intersections:
        return _empty_figure("No intersections meet minimum size threshold", dark_mode), []

    # Sort by size and limit
    sorted_intersections = sorted(intersections.items(), key=lambda x: len(x[1]), reverse=True)
    sorted_intersections = sorted_intersections[:max_subsets]

    # Build intersection data for return
    intersection_data = []
    for combo, genes in sorted_intersections:
        intersection_data.append({
            'datasets': list(combo),
            'genes': sorted(list(genes)),
            'count': len(genes)
        })

    # Prepare data
    dataset_names = list(deg_sets.keys())
    n_datasets = len(dataset_names)
    n_intersections = len(sorted_intersections)

    # Theme settings
    if dark_mode:
        bar_color = '#4ECDC4'
        dot_color_active = '#FFFFFF'
        dot_color_inactive = '#555555'
        line_color = '#FFFFFF'
        bg_color = '#1E1E1E'
        plot_bg = '#2D2D2D'
        font_color = '#FFFFFF'
        grid_color = '#404040'
        annotation_color = '#FFD93D'
    else:
        bar_color = '#3498DB'
        dot_color_active = '#2C3E50'
        dot_color_inactive = '#E0E0E0'
        line_color = '#2C3E50'
        bg_color = '#FFFFFF'
        plot_bg = '#FAFAFA'
        font_color = '#2C3E50'
        grid_color = '#E0E0E0'
        annotation_color = '#E74C3C'

    # Create subplots: bar chart on top, dot matrix below
    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.15, 0.85],
        row_heights=[0.65, 0.35],
        horizontal_spacing=0.02,
        vertical_spacing=0.05,
        specs=[
            [{"type": "bar"}, {"type": "bar"}],
            [{"type": "scatter"}, {"type": "scatter"}]
        ]
    )

    # Intersection sizes (top right - main bar chart)
    intersection_sizes = [len(genes) for _, genes in sorted_intersections]
    x_positions = list(range(n_intersections))

    # Build hover text with gene names for small intersections
    hover_texts = []
    for combo, genes in sorted_intersections:
        gene_list = sorted(list(genes))
        datasets_str = " ∩ ".join(combo)
        if len(gene_list) <= show_gene_labels_threshold:
            genes_str = ", ".join(gene_list)
            hover_texts.append(
                f"<b>{len(gene_list)} genes</b><br>"
                f"<b>Datasets:</b> {datasets_str}<br>"
                f"<b>Genes:</b> {genes_str}"
            )
        else:
            preview = ", ".join(gene_list[:5]) + f"... (+{len(gene_list)-5} more)"
            hover_texts.append(
                f"<b>{len(gene_list)} genes</b><br>"
                f"<b>Datasets:</b> {datasets_str}<br>"
                f"<b>Preview:</b> {preview}"
            )

    fig.add_trace(
        go.Bar(
            x=x_positions,
            y=intersection_sizes,
            marker_color=bar_color,
            text=intersection_sizes,
            textposition='outside',
            textfont=dict(size=10, color=font_color),
            hovertemplate='%{customdata}<extra></extra>',
            customdata=hover_texts,
            showlegend=False
        ),
        row=1, col=2
    )

    # Add gene name annotations for very small intersections (≤5 genes)
    for i, (combo, genes) in enumerate(sorted_intersections):
        if len(genes) <= 5 and len(genes) > 0:
            gene_list = sorted(list(genes))
            annotation_text = "<br>".join(gene_list)
            fig.add_annotation(
                x=i,
                y=len(genes),
                text=annotation_text,
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=1,
                arrowcolor=annotation_color,
                ax=0,
                ay=-40 - (len(genes) * 12),
                font=dict(size=8, color=annotation_color),
                align="center",
                xref="x2",
                yref="y2"
            )

    # Set sizes (top left - horizontal bar chart)
    set_sizes = [len(deg_sets[name]) for name in dataset_names]

    fig.add_trace(
        go.Bar(
            y=list(range(n_datasets)),
            x=set_sizes,
            orientation='h',
            marker_color=bar_color,
            text=set_sizes,
            textposition='outside',
            textfont=dict(size=10, color=font_color),
            hovertemplate='<b>%{y}</b><br>%{x} genes<extra></extra>',
            showlegend=False
        ),
        row=2, col=1
    )

    # Dot matrix (bottom right)
    for i, (combo, genes) in enumerate(sorted_intersections):
        combo_set = set(combo)

        for j, dataset in enumerate(dataset_names):
            is_active = dataset in combo_set

            fig.add_trace(
                go.Scatter(
                    x=[i],
                    y=[j],
                    mode='markers',
                    marker=dict(
                        size=14,
                        color=dot_color_active if is_active else dot_color_inactive,
                        line=dict(width=2, color=line_color) if is_active else dict(width=0)
                    ),
                    hoverinfo='skip',
                    showlegend=False
                ),
                row=2, col=2
            )

        # Draw connecting lines for active dots
        active_indices = [j for j, dataset in enumerate(dataset_names) if dataset in combo_set]
        if len(active_indices) > 1:
            fig.add_trace(
                go.Scatter(
                    x=[i] * len(active_indices),
                    y=active_indices,
                    mode='lines',
                    line=dict(color=line_color, width=2),
                    hoverinfo='skip',
                    showlegend=False
                ),
                row=2, col=2
            )

    # Update layout
    fig.update_layout(
        title=dict(
            text="<b>UpSet Plot: DEG Overlap Analysis</b>",
            font=dict(size=18, color=font_color),
            x=0.5,
            xanchor='center'
        ),
        paper_bgcolor=bg_color,
        plot_bgcolor=plot_bg,
        font=dict(color=font_color),
        height=500 + (n_datasets * 25),
        margin=dict(l=150, r=40, t=100, b=40),
        showlegend=False
    )

    # Update axes for main bar chart (top right)
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        row=1, col=2
    )
    fig.update_yaxes(
        title=dict(text="Intersection Size", font=dict(color=font_color)),
        tickfont=dict(color=font_color),
        gridcolor=grid_color,
        row=1, col=2
    )

    # Update axes for set size bar chart (top left)
    fig.update_xaxes(
        title=dict(text="Set Size", font=dict(color=font_color)),
        tickfont=dict(color=font_color),
        gridcolor=grid_color,
        row=2, col=1
    )
    fig.update_yaxes(
        tickmode='array',
        tickvals=list(range(n_datasets)),
        ticktext=dataset_names,
        tickfont=dict(color=font_color),
        row=2, col=1
    )

    # Update axes for dot matrix (bottom right)
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        range=[-0.5, n_intersections - 0.5],
        row=2, col=2
    )
    fig.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        range=[-0.5, n_datasets - 0.5],
        row=2, col=2
    )

    return fig, intersection_data


def _calculate_intersections(deg_sets: Dict[str, Set[str]]) -> Dict[Tuple[str, ...], Set[str]]:
    """
    Calculate all possible intersections between sets.

    Returns exclusive intersections (genes in exactly these sets and no others).
    """
    dataset_names = list(deg_sets.keys())
    all_genes = set.union(*deg_sets.values()) if deg_sets else set()

    # For each gene, find which sets it belongs to
    gene_membership = {}
    for gene in all_genes:
        members = tuple(sorted([name for name, genes in deg_sets.items() if gene in genes]))
        if members not in gene_membership:
            gene_membership[members] = set()
        gene_membership[members].add(gene)

    return gene_membership


def _empty_figure(message: str, dark_mode: bool) -> go.Figure:
    """Create an empty figure with a message."""
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
        yaxis=dict(visible=False)
    )
    return fig


def get_intersection_genes(
    deg_sets: Dict[str, Set[str]],
    selected_datasets: List[str],
    exclusive: bool = True
) -> Set[str]:
    """
    Get genes in the intersection of selected datasets.

    Args:
        deg_sets: Dictionary of dataset name to gene set
        selected_datasets: List of datasets to intersect
        exclusive: If True, return only genes in exactly these datasets

    Returns:
        Set of gene names
    """
    if not selected_datasets:
        return set()

    # Genes in all selected datasets
    intersection = set.intersection(*[deg_sets[name] for name in selected_datasets])

    if exclusive:
        # Remove genes that are also in other datasets
        other_datasets = [name for name in deg_sets if name not in selected_datasets]
        for name in other_datasets:
            intersection -= deg_sets[name]

    return intersection


def get_genes_in_n_datasets(
    deg_sets: Dict[str, Set[str]],
    min_datasets: int = 1
) -> Dict[str, List[str]]:
    """
    Get genes that appear in at least n datasets.

    Args:
        deg_sets: Dictionary of dataset name to gene set
        min_datasets: Minimum number of datasets

    Returns:
        Dictionary mapping gene name to list of datasets it appears in
    """
    gene_to_datasets = {}

    for dataset_name, genes in deg_sets.items():
        for gene in genes:
            if gene not in gene_to_datasets:
                gene_to_datasets[gene] = []
            gene_to_datasets[gene].append(dataset_name)

    # Filter to genes in at least min_datasets
    return {
        gene: datasets
        for gene, datasets in gene_to_datasets.items()
        if len(datasets) >= min_datasets
    }
